use flate2::read::MultiGzDecoder;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write, BufWriter};
use std::time::Instant;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

struct Record {
    atomic_ids: Vec<String>,
    path_info: f32,
    path_raf: f32,
}

#[derive(Clone)]
struct PhasedDist {
    p0: f32,
    p1: f32,
}

fn smart_open(filename: &str) -> Box<dyn BufRead> {
    let file = File::open(filename).expect("Cannot open file");
    if filename.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}

// Retained for once-per-variant metadata (AF, INFO) where allocation doesn't matter
fn format_float(val: f32) -> String {
    let s = format!("{:.3}", val);
    let trimmed = s.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

// Zero-allocation f32 formatter for the Dosage (DS) field
fn write_float_trim(w: &mut impl Write, val: f32) -> io::Result<()> {
    let mut buf = [0u8; 32];
    let mut cursor = std::io::Cursor::new(&mut buf[..]);
    write!(&mut cursor, "{:.3}", val)?;
    let len = cursor.position() as usize;
    let mut s = &buf[..len];
    while s.ends_with(b"0") { s = &s[..s.len()-1]; }
    if s.ends_with(b".") { s = &s[..s.len()-1]; }
    if s.is_empty() {
        w.write_all(b"0")
    } else {
        w.write_all(s)
    }
}

// Zero-allocation per-mille formatter for Genotype Probabilities (GP)
fn write_permille(w: &mut impl Write, v: i32) -> io::Result<()> {
    if v == 0 { w.write_all(b"0") }
    else if v == 1000 { w.write_all(b"1") }
    else if v % 100 == 0 { write!(w, "0.{}", v / 100) }
    else if v % 10 == 0 { write!(w, "0.{:02}", v / 10) }
    else { write!(w, "0.{:03}", v) }
}

fn process_group(
    group_lines: &[(String, String)],
    id_buffer: &HashMap<String, (u32, String, String, String, String)>,
    max_alleles: usize,
    out_handle: &mut impl Write,
) {
    if group_lines.is_empty() { return; }

    let first_line_fields: Vec<&str> = group_lines[0].0.trim_end().split('\t').collect();
    if first_line_fields.len() <= 9 {
        panic!("Error: VCF does not contain sample columns.");
    }
    let num_samples = first_line_fields.len() - 9;
    let chrom = first_line_fields[0].to_string();

    let num_alleles = group_lines.len();
    let mut records: Vec<Record> = Vec::with_capacity(num_alleles);
    let mut all_atomic_ids = HashSet::new();

    let mut hap_probs: Vec<Vec<(f32, f32)>> = vec![vec![(0.0, 0.0); num_alleles]; num_samples];

    for (a, (line, site_info)) in group_lines.iter().enumerate() {
        let fields: Vec<&str> = line.trim_end().split('\t').collect();

        // 1. Strictly enforce INFO and RAF metadata
        let mut path_info: Option<f32> = None;
        let mut path_raf: Option<f32> = None;
        for item in fields[7].split(';') {
            if let Some(v) = item.strip_prefix("INFO=") {
                path_info = Some(v.parse::<f32>().expect("Error: INFO is not a valid float"));
            } else if let Some(v) = item.strip_prefix("RAF=") {
                path_raf = Some(v.parse::<f32>().expect("Error: RAF is not a valid float"));
            }
        }
        let path_info = path_info.expect("Error: Missing INFO field in GLIMPSE2 output");
        let path_raf = path_raf.expect("Error: Missing RAF field in GLIMPSE2 output");

        let mut atomic_ids = Vec::new();
        for item in site_info.split(';') {
            if let Some(id_str) = item.strip_prefix("ID=") {
                let replaced_id = id_str.replace(',', ":");
                for j in replaced_id.split(':').map(|s| s.trim()) {
                    if id_buffer.contains_key(j) {
                        atomic_ids.push(j.to_string());
                        all_atomic_ids.insert(j.to_string());
                    } else {
                        panic!("Error: Variant ID '{}' not found in the ID buffer. Ensure your biallelic VCF contains this ID and the window size is large enough.", j);
                    }
                }
                break;
            }
        }

        records.push(Record { atomic_ids, path_info, path_raf });

        // 2. Strictly enforce FORMAT header
        let fmt: Vec<&str> = fields[8].split(':').collect();
        let gt_idx = fmt.iter().position(|&x| x == "GT").expect("Error: Missing GT in FORMAT string");
        let gp_idx = fmt.iter().position(|&x| x == "GP").expect("Error: Missing GP in FORMAT string");

        for s in 0..num_samples {
            let sample_data = fields[9 + s];
            if sample_data == "." {
                panic!("Error: GLIMPSE2 output missing sample data ('.')");
            }
            
            // Allocation-free indexing of the sample format tokens
            let mut gt_val = "";
            let mut gp_str = "";
            for (i, val) in sample_data.split(':').enumerate() {
                if i == gt_idx { gt_val = val; }
                else if i == gp_idx { gp_str = val; }
            }
            
            if gt_val == "" || gt_val == "." || gp_str == "" || gp_str == "." {
                panic!("Error: Missing GT or GP value for sample");
            }

            let mut gp_iter = gp_str.split(',');
            let _gp0 = gp_iter.next().expect("Error: Malformed GP (missing gp0)");
            let gp1 = gp_iter.next().expect("Error: Malformed GP (missing gp1)").parse::<f32>().expect("Error: GP1 is not a valid float");
            let gp2 = gp_iter.next().expect("Error: Malformed GP (missing gp2)").parse::<f32>().expect("Error: GP2 is not a valid float");

            let (p0, p1) = if gt_val.starts_with("1|0") {
                (gp2 + gp1, gp2)
            } else if gt_val.starts_with("0|1") {
                (gp2, gp2 + gp1)
            } else {
                (gp2 + (gp1 / 2.0), gp2 + (gp1 / 2.0))
            };

            hap_probs[s][a] = (p0.clamp(1e-5, 1.0 - 1e-5), p1.clamp(1e-5, 1.0 - 1e-5));
        }
    }

    let mut atomic_sample_dists: HashMap<String, Vec<PhasedDist>> = HashMap::new();
    for id in &all_atomic_ids {
        atomic_sample_dists.insert(id.clone(), vec![PhasedDist { p0: 0.0, p1: 0.0 }; num_samples]);
    }

    let mut allele_scores: Vec<(usize, f32)> = Vec::with_capacity(num_alleles);

    for s in 0..num_samples {
        allele_scores.clear();
        for a in 0..num_alleles {
            let score = hap_probs[s][a].0 + hap_probs[s][a].1;
            allele_scores.push((a, score));
        }
        allele_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        let m = allele_scores.len().min(max_alleles);
        let top_alleles: Vec<usize> = allele_scores.iter().take(m).map(|x| x.0).collect();

        let mut z0 = 1.0_f32;
        let mut z1 = 1.0_f32;
        let mut w0_map = HashMap::new();
        let mut w1_map = HashMap::new();

        for &a in &top_alleles {
            let (p0, p1) = hap_probs[s][a];
            let w0 = p0 / (1.0 - p0);
            let w1 = p1 / (1.0 - p1);
            w0_map.insert(a, w0);
            w1_map.insert(a, w1);
            z0 += w0;
            z1 += w1;
        }

        for &a in &top_alleles {
            let norm_p0 = w0_map[&a] / z0;
            let norm_p1 = w1_map[&a] / z1;
            for atomic_id in &records[a].atomic_ids {
                if let Some(dists) = atomic_sample_dists.get_mut(atomic_id) {
                    dists[s].p0 += norm_p0;
                    dists[s].p1 += norm_p1;
                }
            }
        }
    }

    let mut sorted_atomic_vars: Vec<_> = all_atomic_ids.into_iter().collect();
    sorted_atomic_vars.sort_by(|a, b| {
        let data_a = id_buffer.get(a);
        let data_b = id_buffer.get(b);
        match (data_a, data_b) {
            (Some(da), Some(db)) => {
                da.0.cmp(&db.0)
                    .then_with(|| da.2.cmp(&db.2))
                    .then_with(|| da.3.cmp(&db.3))
            }
            _ => std::cmp::Ordering::Equal,
        }
    });

    let empty_dists = vec![PhasedDist { p0: 0.0, p1: 0.0 }; num_samples];

    for assigned_id in sorted_atomic_vars {
        let var_data = &id_buffer[&assigned_id];
        let coord = var_data.0;

        let mut new_info = vec![format!("ID={}", assigned_id)];
        if !var_data.4.is_empty() {
            new_info.push(format!("RAF={}", var_data.4));
        }

        let dists = atomic_sample_dists.get(&assigned_id).unwrap_or(&empty_dists);
        let mut ds_sum = 0.0_f32;
        
        for s in 0..num_samples {
            let dist = &dists[s];
            let p0 = dist.p0.clamp(0.0, 1.0);
            let p1 = dist.p1.clamp(0.0, 1.0);

            let gp1_raw = p0 * (1.0 - p1) + (1.0 - p0) * p1;
            let gp2_raw = p0 * p1;

            let gp1_q = (gp1_raw * 1000.0).round() as i32 as f32 / 1000.0;
            let gp2_q = (gp2_raw * 1000.0).round() as i32 as f32 / 1000.0;
            ds_sum += gp1_q + 2.0 * gp2_q;
        }

        let n_tar_haps = 2.0 * num_samples as f32;
        let safe_n = n_tar_haps.max(1e-9);
        let af = ds_sum / safe_n;
        new_info.push(format!("AF={}", format_float(af)));
        
        // 3. Simplified INFO calculation with guaranteed data
        let mut info_num = 0.0_f32;
        let mut info_den = 0.0_f32;
        let mut sum_info = 0.0_f32;
        let mut count = 0.0_f32;
        
        for rec in &records {
            if rec.atomic_ids.contains(&assigned_id) {
                info_num += rec.path_info * rec.path_raf;
                info_den += rec.path_raf;
                sum_info += rec.path_info;
                count += 1.0;
            }
        }
        
        if info_den > 0.0 {
            new_info.push(format!("INFO={}", format_float(info_num / info_den)));
        } else {
            // Safe unweighted fallback if all RAFs are exactly 0.0
            new_info.push(format!("INFO={}", format_float(sum_info / count.max(1.0))));
        }

        write!(out_handle, "{}\t{}\t{}\t{}\t{}\t.\t.\t{}\tGT:DS:GP",
            chrom, coord, var_data.1, var_data.2, var_data.3, new_info.join(";")
        ).unwrap();

        for s in 0..num_samples {
            let dist = &dists[s];
            let p0 = dist.p0.clamp(0.0, 1.0);
            let p1 = dist.p1.clamp(0.0, 1.0);

            let hap0_gt = if p0 > 0.5 { '1' } else { '0' };
            let hap1_gt = if p1 > 0.5 { '1' } else { '0' };

            let gp0_raw = (1.0 - p0) * (1.0 - p1);
            let gp1_raw = p0 * (1.0 - p1) + (1.0 - p0) * p1;
            let gp2_raw = p0 * p1;

            let mut v0 = (gp0_raw * 1000.0).round() as i32;
            let mut v1 = (gp1_raw * 1000.0).round() as i32;
            let mut v2 = (gp2_raw * 1000.0).round() as i32;

            let diff = 1000 - (v0 + v1 + v2);
            if diff != 0 {
                if v0 >= v1 && v0 >= v2 { v0 += diff; }
                else if v1 >= v0 && v1 >= v2 { v1 += diff; }
                else { v2 += diff; }
            }

            write!(out_handle, "\t{}|{}:", hap0_gt, hap1_gt).unwrap();
            write_float_trim(out_handle, p0 + p1).unwrap();
            write!(out_handle, ":").unwrap();
            write_permille(out_handle, v0).unwrap();
            write!(out_handle, ",").unwrap();
            write_permille(out_handle, v1).unwrap();
            write!(out_handle, ",").unwrap();
            write_permille(out_handle, v2).unwrap();
        }
        writeln!(out_handle).unwrap();
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: cat <multiallelic VCF> | {} <biallelic ID VCF> <sites VCF> [max_alleles] [window_size]", args[0]);
        std::process::exit(1);
    }

    let vcf_path = &args[1];
    let sites_path = &args[2];
    let max_alleles: usize = if args.len() > 3 {
        args[3].parse().unwrap_or(10)
    } else {
        10
    };
    let window_size: u32 = if args.len() > 4 {
        args[4].parse().unwrap_or(500_000)
    } else {
        500_000
    };

    let mut id_iter = smart_open(vcf_path).lines().peekable();
    let mut sites_iter = smart_open(sites_path).lines();

    let mut next_site = sites_iter.next();
    while let Some(Ok(ref l)) = next_site {
        if l.starts_with('#') {
            next_site = sites_iter.next();
        } else {
            break;
        }
    }

    let stdout = io::stdout();
    let mut out_handle = BufWriter::new(stdout.lock());
    let stdin = io::stdin();

    let mut current_pos: Option<String> = None;
    let mut current_chrom: String = String::new();
    let mut group: Vec<(String, String)> = Vec::new();
    let mut records_processed: usize = 0;

    let mut id_buffer: HashMap<String, (u32, String, String, String, String)> = HashMap::new();
    let mut active_id_chrom = String::new();

    let start_time = Instant::now();
    eprintln!("Starting Phased Joint-Distribution projection (Max Alleles: {}, Window Size: {})...", max_alleles, window_size);

    for line_result in stdin.lock().lines() {
        let line = line_result.unwrap();

        if line.starts_with('#') {
            if line.starts_with("##INFO=<ID=RAF,") || 
               line.starts_with("##INFO=<ID=AF,") || 
               line.starts_with("##INFO=<ID=INFO,") {
                continue; 
            }

            if line.starts_with("#CHROM") {
                writeln!(out_handle, "##INFO=<ID=RAF,Number=A,Type=Float,Description=\"ALT allele frequency in the reference panel\">").unwrap();
                writeln!(out_handle, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"ALT allele frequency computed from rounded GLIMPSE2 output DS/GP field across target samples\">").unwrap();
                writeln!(out_handle, "##INFO=<ID=INFO,Number=A,Type=Float,Description=\"RAF-weighted average INFO score across the parent bubble for paths containing the variant\">").unwrap();
                writeln!(out_handle, "{}", line).unwrap();
            } else if !line.contains("INFO=<ID=AK") && !line.contains("FORMAT=<ID=GL") && !line.contains("FORMAT=<ID=KC") {
                writeln!(out_handle, "{}", line).unwrap();
            }
            continue;
        }

        let site_line = next_site.unwrap_or_else(|| panic!("Error: Sites VCF ran out of records before the main VCF stream")).unwrap();
        next_site = sites_iter.next();

        let site_fields: Vec<&str> = site_line.splitn(9, '\t').collect();
        if site_fields.len() < 5 {
            panic!("Error: Sites VCF line is malformed or missing required columns: {}", site_line);
        }
        let site_chrom = site_fields[0];
        let site_pos = site_fields[1];
        let site_ref = site_fields[3];
        let site_alt = site_fields[4];
        let site_info = if site_fields.len() > 7 { site_fields[7].to_string() } else { "".to_string() };

        let fields: Vec<&str> = line.splitn(6, '\t').collect();
        if fields.len() < 5 { continue; }
        let chrom = fields[0].to_string();
        let pos = fields[1].to_string();
        let ref_seq = fields[3];
        let alt_seq = fields[4];

        if chrom != site_chrom || pos != site_pos || ref_seq != site_ref || alt_seq != site_alt {
            panic!(
                "Error: Lockstep synchronization failed! Variants do not match.\nMain Stream: {}:{} {} -> {}\nSites Stream: {}:{} {} -> {}",
                chrom, pos, ref_seq, alt_seq, site_chrom, site_pos, site_ref, site_alt
            );
        }

        records_processed += 1;
        if records_processed % 10_000 == 0 {
            let elapsed = start_time.elapsed().as_secs();
            let hours = elapsed / 3600;
            let mins = (elapsed % 3600) / 60;
            let secs = elapsed % 60;
            eprintln!(
                "[{:02}:{:02}:{:02}] Processed {} input records... (Currently at {}:{})",
                hours, mins, secs, records_processed, chrom, pos
            );
        }

        if current_pos.is_none() {
            current_pos = Some(pos.clone());
            current_chrom = chrom.clone();
        }

        if pos != *current_pos.as_ref().unwrap() || chrom != current_chrom {
            process_group(&group, &id_buffer, max_alleles, &mut out_handle);
            group.clear();
            current_pos = Some(pos.clone());
            current_chrom = chrom.clone();
        }

        if group.is_empty() {
            let pos_u32 = pos.parse::<u32>().unwrap_or(0);

            if current_chrom != active_id_chrom {
                id_buffer.clear();
                active_id_chrom = current_chrom.clone();
            }

            id_buffer.retain(|_, v| v.0 + window_size >= pos_u32);

            while let Some(Ok(peek_line)) = id_iter.peek() {
                if peek_line.starts_with('#') {
                    id_iter.next();
                    continue;
                }
                let peek_fields: Vec<&str> = peek_line.splitn(3, '\t').collect();
                let peek_chrom = peek_fields[0];

                if peek_chrom != current_chrom {
                    if active_id_chrom != current_chrom {
                        id_iter.next();
                        continue;
                    } else {
                        break;
                    }
                }

                let peek_pos: u32 = peek_fields[1].parse().unwrap_or(0);
                if peek_pos > pos_u32 + window_size {
                    break;
                }
                if peek_pos + window_size < pos_u32 {
                    id_iter.next();
                    continue;
                }

                let pop_line = id_iter.next().unwrap().unwrap();
                let pop_fields: Vec<&str> = pop_line.split('\t').collect();
                let orig_id = pop_fields[2].to_string();
                let ref_seq = pop_fields[3].to_string();
                let alt_seq = pop_fields[4].to_string();
                
                let mut id_val = String::new();
                let mut af_val = String::new();
                let mut ac_val = String::new();
                let mut an_val = String::new();
                
                for item in pop_fields[7].split(';') {
                    if let Some(v) = item.strip_prefix("ID=") {
                        id_val = v.to_string();
                    } else if let Some(v) = item.strip_prefix("AF=") {
                        af_val = v.to_string();
                    } else if let Some(v) = item.strip_prefix("AC=") {
                        ac_val = v.to_string();
                    } else if let Some(v) = item.strip_prefix("AN=") {
                        an_val = v.to_string();
                    }
                }
                
                if af_val.is_empty() && !ac_val.is_empty() && !an_val.is_empty() {
                    if let (Ok(ac), Ok(an)) = (ac_val.parse::<f32>(), an_val.parse::<f32>()) {
                        if an > 0.0 {
                            af_val = format_float((ac / an) as f32);
                        }
                    }
                }
                
                if !id_val.is_empty() {
                    id_buffer.insert(id_val, (peek_pos, orig_id, ref_seq, alt_seq, af_val));
                }
            }
        }

        group.push((line, site_info));
    }

    if !group.is_empty() {
        process_group(&group, &id_buffer, max_alleles, &mut out_handle);
    }

    out_handle.flush().unwrap();

    let total_elapsed = start_time.elapsed().as_secs();
    let t_hours = total_elapsed / 3600;
    let t_mins = (total_elapsed % 3600) / 60;
    let t_secs = total_elapsed % 60;
    eprintln!(
        "Finished! Processed a total of {} input records in {:02}:{:02}:{:02}.",
        records_processed, t_hours, t_mins, t_secs
    );
}
