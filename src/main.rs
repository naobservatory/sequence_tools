use std::io;
use std::io::Read;
use std::io::Write;
use clap::Parser;
use std::collections::HashMap;

#[derive(Parser)]
struct Cli {
    /// The patterns to look for.  Each can be separated by commas, and each
    /// comma-separated group will be highlighted in a different color.
    /// Ex: ABCD,EFGH IJKL,MNOP
    #[clap(multiple=true)]
    needle_groups: Vec<String>,

    /// Maximum substitions
    #[clap(long, value_parser, default_value_t = 3)]
    max_subs: u8,
}

fn ignore(b: u8) -> bool {
    b == b'\n'
}

#[derive(Eq, Hash, PartialEq)]
struct PrintSpot {
    loc: i64,
    before: bool,
    detail: bool,
}

const COLORS: &'static [&'static str] = &[
    "\x1b[1;31m",  // red
    "\x1b[1;32m",  // green
    "\x1b[1;33m",  // yellow
    "\x1b[1;34m",  // blue
    "\x1b[1;35m",  // magenta
    "\x1b[1;36m",  // cyan
    "\x1b[1;37m",  // white
];

const BEGIN_UNDERLINE: &'static str = "\x1b[4m";
const END_UNDERLINE: &'static str = "\x1b[24m";
const END_COLORS: &'static str = "\x1b[0m";

fn emit(b: u8, file_byte_index: i64, to_print: &mut HashMap<PrintSpot, Vec<u8>>) {
    if file_byte_index < 0 {
        return;
    }

    if let Some(needs_printing) = to_print.remove(&PrintSpot{
        loc:file_byte_index, before:true, detail:false}) {
        io::stdout().write(&needs_printing).unwrap();
    }
    if let Some(needs_printing) = to_print.remove(&PrintSpot{
        loc:file_byte_index, before:true, detail:true}) {
        io::stdout().write(&needs_printing).unwrap();
    }
    io::stdout().write(&[b]).unwrap();
    if let Some(needs_printing) = to_print.remove(&PrintSpot{
        loc:file_byte_index, before:false, detail:true}) {
        io::stdout().write(&needs_printing).unwrap();
    }
    if let Some(needs_printing) = to_print.remove(&PrintSpot{
        loc:file_byte_index, before:false, detail:false}) {
        io::stdout().write(&needs_printing).unwrap();
    }
}

fn main() {
    let args = Cli::parse();

    let mut needle_groups: Vec<Vec<String>> = Vec::new();
    for needle_group_str in args.needle_groups {
        needle_groups.push(
            needle_group_str.split(",")
                .map(str::to_string)
                .collect::<Vec<String>>());
    }

    let mut max_needle_len = 0;
    for needle_group in &needle_groups {
        for needle in needle_group {
            if needle.len() > max_needle_len {
                max_needle_len = needle.len();
            }
        }
    }

    // Allow up to 1 ignored character within the haystack per needle instance.
    let max_ignored = 1;
    let max_lookback: i64 = (max_needle_len + max_ignored).try_into().unwrap();

    let mut hist: Vec<u8> = Vec::new();
    for _ in 0..(max_lookback) {
        hist.push(0);
    }
    let hist_len = hist.len();
    let mut hist_pos = hist.len();

    // The byte index in the input corresponding to where we're currently
    // writing.  We're always reading ahead of this by max_lookback bytes.
    let mut file_byte_index = -max_lookback - 1;

    let mut to_print: HashMap<PrintSpot, Vec<u8>> = HashMap::new();

    let mut mismatches = Vec::new();
    for b in io::stdin().bytes() {
        file_byte_index += 1;

        let b = b.unwrap();
        hist_pos = (hist_pos + 1) % hist_len;
        emit(hist[hist_pos], file_byte_index, &mut to_print);
        hist[hist_pos] = b;

        if b == b'\n' {
            continue;
        }

        for (needle_group_index,
             needle_group) in needle_groups.iter().enumerate() {
            for needle in needle_group {
            mismatches.clear();
            let mut matched = false;

            let mut needle_loc = needle.len() - 1;
            let mut hist_loc = hist_pos + hist_len;
            let mut ignored = 0;

            loop {
                while ignored < max_ignored && ignore(hist[hist_loc % hist_len]) {
                    hist_loc -= 1;
                    ignored += 1;
                }

                if hist[hist_loc % hist_len] != needle.as_bytes()[needle_loc] {
                    mismatches.push((needle.len() - needle_loc + ignored) as i64)
                }

                if mismatches.len() > args.max_subs.into() {
                    break;
                }

                if needle_loc == 0 {
                    matched = true;
                    break;
                }

                needle_loc -= 1;
                hist_loc -= 1;
            }
            if matched {
                let current_file_location = file_byte_index + max_lookback;
                to_print.insert(PrintSpot {
                    loc: current_file_location + 1 - (needle.len() + ignored)
                        as i64,
                    before: true,
                    detail: false,
                }, COLORS[needle_group_index].as_bytes().to_vec());
                to_print.insert(PrintSpot {
                    loc: current_file_location,
                    before: false,
                    detail: false,
                }, END_COLORS.as_bytes().to_vec());

                for mismatch in &mismatches {
                    to_print.insert(PrintSpot {
                        loc: current_file_location + 1 - mismatch,
                        before: true,
                        detail: true,
                    }, BEGIN_UNDERLINE.as_bytes().to_vec());
                    to_print.insert(PrintSpot {
                        loc: current_file_location - mismatch + 1,
                        before: false,
                        detail: true,
                    }, END_UNDERLINE.as_bytes().to_vec());
                }}
            }
        }
    }

    // Now we need to print the last of the history that hasn't already been
    // printed.
    for i in (hist_pos + 1)..(hist_pos + hist_len + 1) {
        file_byte_index += 1;
        emit(hist[i % hist_len], file_byte_index, &mut to_print);
    }
}
