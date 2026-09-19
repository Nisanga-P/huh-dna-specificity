import requests
from bs4 import BeautifulSoup
import re
import csv
import time
import os

# ---- config ----
BASE_URL = "https://bioinfo-mml.sjtu.edu.cn/oriTDB2/browse_oriT_result.php?acc={}"
OUTPUT_CSV = "oritdb_pairs.csv"

FALLBACK_START_ACC = 100001   # only used if the CSV doesn't exist yet / is empty
BATCH_SIZE = 200                # how many *valid* records to pull this run
SLEEP_SECONDS = 1.5             # conservative delay between requests
MAX_CONSECUTIVE_MISSES = 150    # stop early if we hit a long dead zone (end of valid IDs)

FIELDNAMES = [
    "oriTDB_accession", "oriTDB_id", "nic site",
    "relaxase_family", "relaxase_sequence", "host",
]


def get_resume_point():
    """Look at the existing CSV (if any) and figure out where to continue."""
    if not os.path.exists(OUTPUT_CSV):
        return FALLBACK_START_ACC

    max_acc = None
    with open(OUTPUT_CSV, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                acc = int(row["oriTDB_accession"])
                if max_acc is None or acc > max_acc:
                    max_acc = acc
            except (ValueError, KeyError):
                continue

    return (max_acc + 1) if max_acc is not None else FALLBACK_START_ACC


def grab(pattern, text, flags=re.MULTILINE):
    m = re.search(pattern, text, flags)
    return m.group(1).strip() if m else ""


import re

FIELDNAMES = [
    "oriTDB_accession", "oriT_sequence", "nic site",
    "relaxase_family", "relaxase_sequence", "host",
]


def parse_entry(text):
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    norm = "\n".join(lines)

    # oriTDB accession -- anchored on the unique "DB ID" label
    acc_match = re.search(r"DB ID\n(\d+)", norm)
    if not acc_match:
        return None
    oriTDB_accession = acc_match.group(1)

    # host / organism
    host_match = re.search(r"Organism\n(.*?)\nSequence Completeness", norm, re.DOTALL)
    host = " ".join(host_match.group(1).split("\n")) if host_match else ""

    # nic site location
    nic_match = re.search(r"Location of\nnic\nsite\n(\S+)", norm)
    nic_site = nic_match.group(1) if nic_match else ""

    # oriT DNA sequence -- header starts with ">oriT_", sequence is ACGTN only
    orit_seq_match = re.search(r">oriT\S*[^\n]*\n([ACGTNacgtn]+)", norm)
    oriT_sequence = orit_seq_match.group(1) if orit_seq_match else ""

    # relaxase protein sequence -- anchored specifically to the "Relaxase protein sequence" section
    relax_seq_match = re.search(
        r"Relaxase protein sequence\nDownload\nLength:.*?\n>[^\n]*\n((?:[A-Z]{10,}\n?)+)",
        norm, re.DOTALL
    )
    relaxase_sequence = relax_seq_match.group(1).replace("\n", "") if relax_seq_match else ""

    relaxase_family = ""
    fam_match = re.search(r">\S+.*?\b(\w+)\s+family relaxase", norm)
    if fam_match:
        relaxase_family = fam_match.group(1)
    else:
        name_match = re.search(r"Relaxase\n(\w+)_", norm)
        if name_match:
            relaxase_family = name_match.group(1)

    
    # relaxase family
    relaxase_family = ""
    fam_match = re.search(r">\S+.*?\b(\w+)\s+family relaxase", norm)
    if fam_match:
        relaxase_family = fam_match.group(1)
    else:
        name_match = re.search(r"Relaxase\n(\w+)_", norm)
        if name_match:
            relaxase_family = name_match.group(1)

    return {
        "oriTDB_accession": oriTDB_accession,
        "oriT_sequence": oriT_sequence,
        "nic site": nic_site,
        "relaxase_family": relaxase_family,
        "relaxase_sequence": relaxase_sequence,
        "host": host,
    }


def fetch_text(acc, session):
    resp = session.get(BASE_URL.format(acc), timeout=15)
    if resp.status_code != 200:
        return None
    soup = BeautifulSoup(resp.text, "html.parser")
    return soup.get_text(separator="\n")


def main():
    session = requests.Session()
    session.headers.update({
        "User-Agent": "Research script (contact: you@example.com) - pulling public oriT/relaxase records for academic use"
    })

    acc = get_resume_point()
    print(f"Resuming from accession {acc}")

    file_exists = os.path.exists(OUTPUT_CSV)
    written_this_run = 0
    consecutive_misses = 0

    with open(OUTPUT_CSV, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        if not file_exists:
            writer.writeheader()

        while written_this_run < BATCH_SIZE and consecutive_misses < MAX_CONSECUTIVE_MISSES:
            text = fetch_text(acc, session)
            entry = parse_entry(text) if text else None

            if entry:
                writer.writerow(entry)
                f.flush()
                written_this_run += 1
                consecutive_misses = 0
                print(f"  [{written_this_run}/{BATCH_SIZE}] acc={acc} ok")
            else:
                consecutive_misses += 1

            acc += 1
            time.sleep(SLEEP_SECONDS)

    print(f"\nDone this run: {written_this_run} rows written.")
    print(f"Stopped at accession {acc - 1}.")
    if consecutive_misses >= MAX_CONSECUTIVE_MISSES:
        print("Hit the consecutive-miss limit — may have reached the end of valid IDs, or the page format doesn't match the parser.")


if __name__ == "__main__":
    main()


# import requests
# from bs4 import BeautifulSoup
# import re

# url = "https://bioinfo-mml.sjtu.edu.cn/oriTDB2/browse_oriT_result.php?acc=116159"
# resp = requests.get(url, timeout=15, headers={
#     "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"
# })

# soup = BeautifulSoup(resp.text, "html.parser")
# text = soup.get_text(separator="\n")
# text = "\n".join(line.strip() for line in text.splitlines() if line.strip())

# print(text[:3000])


# sample_text = text[:3000]


# print(parse_entry(sample_text))