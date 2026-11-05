import urllib.request
import re

def detect_dna_repeats(sequence: str, min_unit: int = 3, max_unit: int = 5, min_total: int = 6, max_total: int = 10):
    """
    Detecteaza repetițiile (tandem repeats) dintr-o secventa ADN unde:
    - Unitățile au intre 3 si 5 baze.
    - Lungimea totala a regiunii repetitive este intre 6 si 10 baze.
    """
    results = []
    seq_len = len(sequence)

    for start in range(seq_len):
        for unit_len in range(min_unit, max_unit + 1):
            if start + unit_len > seq_len:
                continue
            
            unit = sequence[start : start + unit_len]
            
            count = 1
            pos = start + unit_len
            while pos + unit_len <= seq_len and sequence[pos : pos + unit_len] == unit:
                count += 1
                pos += unit_len
            
            total_len = count * unit_len
            
            if count >= 2 and min_total <= total_len <= max_total:
                matched_repeat = unit * count
                results.append({
                    'start_pos': start,
                    'end_pos': start + total_len,
                    'unit': unit,
                    'unit_len': unit_len,
                    'count': count,
                    'total_len': total_len,
                    'repeat_seq': matched_repeat
                })

    unique_results = []
    seen = set()
    for res in results:
        key = (res['start_pos'], res['total_len'], res['unit'])
        if key not in seen:
            seen.add(key)
            unique_results.append(res)

    return unique_results


if __name__ == "__main__":
    print("\n=== EXERCIȚIUL 2: DETECTARE REPETIȚII (STR) ===")
    
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text"
    with urllib.request.urlopen(url) as response:
        lines = response.read().decode('utf-8').splitlines()
    dna_sequence = "".join([l.strip().upper() for l in lines if not l.startswith(">")])[:2000]

    print(f"Secvență analizată: {len(dna_sequence)} nucleotide.")

    detected_repeats = detect_dna_repeats(dna_sequence, min_unit=3, max_unit=5, min_total=6, max_total=10)

    print(f"\nS-au găsit {len(detected_repeats)} repetiții valide (lungime totală 6-10 bp, unitate 3-5 bp):\n")
    print(f"{'Pozitie':<10} | {'Unitate (Motif)':<15} | {'Repetari':<10} | {'Lungime Totala':<15} | {'Secventa Repetata'}")
    print("-" * 80)

    for item in detected_repeats[:20]:  
        print(f"{item['start_pos']:<10} | {item['unit']:<15} | {item['count']:<10} | {item['total_len']:<15} | {item['repeat_seq']}")