GENETIC_CODE_3LETTER = {
    # U
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',

    # C
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',

    # A
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',

    # G
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}


def translate_with_start_and_stop(sequence: str) -> str:

    sequence = sequence.strip().upper().replace(' ', '').replace('T', 'U')      
    
    start_index = sequence.find('AUG')
    
    if start_index == -1:
        return "Eroare: Nu s-a gasit niciun codon START (AUG) in secventa!"

    print(f"-> Codonul START (AUG) a fost gasit la pozitia {start_index}.")

    coding_sequence = sequence[start_index:]

    protein = []

    for i in range(0, len(coding_sequence) - len(coding_sequence) % 3, 3):
        codon = coding_sequence[i:i+3]
        amino_acid = GENETIC_CODE_3LETTER.get(codon, '???')

        if amino_acid == 'Stop':
            print(f"-> Codon STOP ({codon}) detectat. Traducerea s-a incheiat.")
            break
        
        protein.append(amino_acid)

    return "-".join(protein)

#test

if __name__ == "__main__":
    test_secventa = "AUA GCA AUG UUU UCU UAA GGG"
    rezultat = translate_with_start_and_stop(test_secventa)
    print("Rezultat proteina:", rezultat)