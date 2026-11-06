import os
import urllib.request

NCBI_ACCESSION = "NC_000932.1"
FASTA_FILE = "plant_genome.fasta"

IR_LENGTH = 12         
MIN_DISTANCE = 100     
MAX_DISTANCE = 5000    


def fetch_plant_genome(accession: str, output_file: str) -> str:
    if not os.path.exists(output_file):
        print(f"Se descarcă genomul vegetal ({accession}) de pe NCBI...")
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
        urllib.request.urlretrieve(url, output_file)
        print(f"Descărcare completă: salvat în '{output_file}'.")

    seq_chunks = []
    with open(output_file, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq_chunks.append(line.strip().upper())

    full_seq = "".join(seq_chunks)
    print(f"Genom vegetal încărcat cu succes (Lungime totală: {len(full_seq):,} baze).\n")
    return full_seq


COMPLEMENT = str.maketrans('ATCG', 'TAGC')

def reverse_complement(dna_str: str) -> str:
    return dna_str.translate(COMPLEMENT)[::-1]


def find_novel_transposons(genome: str, ir_len: int = 12, min_dist: int = 100, max_dist: int = 5000):
    genome_len = len(genome)
    candidates = []

    print(f"Se începe scanarea cu ferestre glisante (S1 len={ir_len}, distanță minimă={min_dist} bp)...")

    for s1_start in range(0, genome_len - ir_len - min_dist):
        s1_seq = genome[s1_start : s1_start + ir_len]

        if len(set(s1_seq)) < 3:
            continue

        s2_target = reverse_complement(s1_seq)

        s2_search_start = s1_start + ir_len + min_dist
        s2_search_end = min(s1_start + ir_len + max_dist, genome_len)
        
        search_region = genome[s2_search_start : s2_search_end]

        match_pos = search_region.find(s2_target)
        
        while match_pos != -1:
            s2_actual_start = s2_search_start + match_pos
            transposon_end = s2_actual_start + ir_len
            transposon_size = transposon_end - s1_start

            candidates.append({
                'start': s1_start,
                'end': transposon_end,
                'size': transposon_size,
                's1_ir': s1_seq,
                's2_ir': s2_target
            })

            match_pos = search_region.find(s2_target, match_pos + 1)

    unique_candidates = []
    seen = set()
    for c in candidates:
        key = (c['start'], c['end'])
        if key not in seen:
            seen.add(key)
            unique_candidates.append(c)

    print(f"Scanare finalizată! S-au găsit {len(unique_candidates):,} candidați de transpozoni.\n")
    return unique_candidates


if __name__ == "__main__":
    print("=== DUAL SLIDING WINDOW TRANSPOSON DETECTOR ===\n")

    plant_dna = fetch_plant_genome(NCBI_ACCESSION, FASTA_FILE)

    results = find_novel_transposons(plant_dna, ir_len=IR_LENGTH, min_dist=MIN_DISTANCE, max_dist=MAX_DISTANCE)

    print("=" * 95)
    print(f"{'#':<4} | {'Start (S1)':<12} | {'End (S2)':<12} | {'Lungime Transpozon':<20} | {'IR (S1) <---> IR (S2)'}")
    print("=" * 95)

    for idx, r in enumerate(results[:15], 1):
        ir_pair = f"{r['s1_ir']} <---> {r['s2_ir']}"
        print(f"{idx:<4} | {r['start']:<12} | {r['end']:<12} | {r['size']:<20} bp | {ir_pair}")

    print("=" * 95)