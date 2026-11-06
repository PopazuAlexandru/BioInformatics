import random

NUCLEOTIDES = ['A', 'C', 'G', 'T']
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def reverse_complement(seq: str) -> str:
    return "".join(COMPLEMENT[b] for b in reversed(seq))

def random_dna(length: int) -> str:
    return "".join(random.choice(NUCLEOTIDES) for _ in range(length))


def generate_transposon_structure(dr_len: int = 5, ir_len: int = 6, core_len: int = 20):
    dr = random_dna(dr_len)                   
    ir_left = random_dna(ir_len)              
    ir_right = reverse_complement(ir_left)    
    core = random_dna(core_len)               
    
    transposon_seq = dr + ir_left + core + ir_right + dr
    return transposon_seq, dr, ir_left, ir_right


def build_synthetic_dna_with_nested_transposons():
    dna = random_dna(200)
    
    t1_seq, t1_dr, t1_irl, t1_irr = generate_transposon_structure(dr_len=5, ir_len=6, core_len=30)
    t2_seq, t2_dr, t2_irl, t2_irr = generate_transposon_structure(dr_len=4, ir_len=5, core_len=25)
    t3_seq, t3_dr, t3_irl, t3_irr = generate_transposon_structure(dr_len=5, ir_len=6, core_len=20)
    
    pos1 = 30
    dna = dna[:pos1] + t1_seq + dna[pos1:]
    
    pos2 = 180
    dna = dna[:pos2] + t2_seq + dna[pos2:]
    
    pos3 = pos1 + 15  
    dna = dna[:pos3] + t3_seq + dna[pos3:]
    
    print(f"Secvență de ADN sintetic generată! Lungime totală: {len(dna)} bp.\n")
    return dna

def detect_transposons(dna_sequence: str, min_dr: int = 4, max_dr: int = 6, min_ir: int = 5, max_ir: int = 7):
    seq_len = len(dna_sequence)
    detected = []

    for dr_len in range(min_dr, max_dr + 1):
        for i in range(seq_len - dr_len):
            dr_left = dna_sequence[i : i + dr_len]
            
            for j in range(i + dr_len + 20, min(i + dr_len + 160, seq_len - dr_len + 1)):
                dr_right = dna_sequence[j : j + dr_len]
                
                if dr_left == dr_right:
                    
                    for ir_len in range(min_ir, max_ir + 1):
                        ir_left_start = i + dr_len
                        ir_left = dna_sequence[ir_left_start : ir_left_start + ir_len]
                        
                        ir_right_end = j
                        ir_right = dna_sequence[ir_right_end - ir_len : ir_right_end]
                        
                        if ir_right == reverse_complement(ir_left):
                            start_pos = i
                            end_pos = j + dr_len
                            
                            detected.append({
                                'start': start_pos,
                                'end': end_pos,
                                'length': end_pos - start_pos,
                                'dr': dr_left,
                                'ir_left': ir_left,
                                'ir_right': ir_right
                            })

    unique_detected = []
    seen = set()
    for item in detected:
        key = (item['start'], item['end'])
        if key not in seen:
            seen.add(key)
            unique_detected.append(item)

    unique_detected.sort(key=lambda x: x['start'])
    return unique_detected


if __name__ == "__main__":
    print("=== SIMULARE ȘI DETECȚIE ELEMENTE TRANSPOZABILE (LAB 6) ===\n")
    
    synthetic_dna = build_synthetic_dna_with_nested_transposons()
    print(f"ADN Sintetic (fragment 1-100): {synthetic_dna[:100]}...\n")
    
    results = detect_transposons(synthetic_dna)
    
    print("=" * 95)
    print(f"{'Transpozon #':<14} | {'Start':<8} | {'End':<8} | {'Lungime':<10} | {'DR':<8} | {'IR Stânga / IR Dreapta'}")
    print("=" * 95)
    
    for idx, t in enumerate(results, 1):
        ir_str = f"{t['ir_left']} <---> {t['ir_right']}"
        print(f"Element {idx:<7} | {t['start']:<8} | {t['end']:<8} | {t['length']:<10} | {t['dr']:<8} | {ir_str}")
        
    print("=" * 95)