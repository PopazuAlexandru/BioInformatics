from collections import Counter


GENETIC_CODE = {
    'TTT': ('Phe', 'Fenilalanină', 'Lactate, carne, pește, fasole, nuci'),
    'TTC': ('Phe', 'Fenilalanină', 'Lactate, carne, pește, fasole, nuci'),
    'TTA': ('Leu', 'Leucină', 'Carne roșie, pui, pește, ouă, brânză, linte'),
    'TTG': ('Leu', 'Leucină', 'Carne roșie, pui, pește, ouă, brânză, linte'),
    'CTT': ('Leu', 'Leucină', 'Carne roșie, pui, pește, ouă, brânză, linte'),
    'CTC': ('Leu', 'Leucină', 'Carne roșie, pui, pește, ouă, brânză, linte'),
    'CTA': ('Leu', 'Leucină', 'Carne roșie, pui, pește, ouă, brânză, linte'),
    'CTG': ('Leu', 'Leucină', 'Carne roșie, pui, pește, ouă, brânză, linte'),
    'ATT': ('Ile', 'Izoleucină', 'Somon, ton, ouă, brânză de vaci, semințe dovleac'),
    'ATC': ('Ile', 'Izoleucină', 'Somon, ton, ouă, brânză de vaci, semințe dovleac'),
    'ATA': ('Ile', 'Izoleucină', 'Somon, ton, ouă, brânză de vaci, semințe dovleac'),
    'ATG': ('Met', 'Metionină (START)', 'Albuș de ou, pește, susan, nuci de Brazilia'),
    'GTT': ('Val', 'Valină', 'Carne, lactate, ciuperci, alune, soia'),
    'GTC': ('Val', 'Valină', 'Carne, lactate, ciuperci, alune, soia'),
    'GTA': ('Val', 'Valină', 'Carne, lactate, ciuperci, alune, soia'),
    'GTG': ('Val', 'Valină', 'Carne, lactate, ciuperci, alune, soia'),
    'TCT': ('Ser', 'Serină', 'Ouă, soia, supă de oase, alune'),
    'TCC': ('Ser', 'Serină', 'Ouă, soia, supă de oase, alune'),
    'TCA': ('Ser', 'Serină', 'Ouă, soia, supă de oase, alune'),
    'TCG': ('Ser', 'Serină', 'Ouă, soia, supă de oase, alune'),
    'CCT': ('Pro', 'Prolină', 'Gelatină, supă de oase, varză, albuș de ou'),
    'CCC': ('Pro', 'Prolină', 'Gelatină, supă de oase, varză, albuș de ou'),
    'CCA': ('Pro', 'Prolină', 'Gelatină, supă de oase, varză, albuș de ou'),
    'CCG': ('Pro', 'Prolină', 'Gelatină, supă de oase, varză, albuș de ou'),
    'ACT': ('Thr', 'Treonină', 'Brânză cottage, pui, pește, linte, susan'),
    'ACC': ('Thr', 'Treonină', 'Brânză cottage, pui, pește, linte, susan'),
    'ACA': ('Thr', 'Treonină', 'Brânză cottage, pui, pește, linte, susan'),
    'ACG': ('Thr', 'Treonină', 'Brânză cottage, pui, pește, linte, susan'),
    'GCT': ('Ala', 'Alanină', 'Carne, pește, lactate, ouă, soia'),
    'GCC': ('Ala', 'Alanină', 'Carne, pește, lactate, ouă, soia'),
    'GCA': ('Ala', 'Alanină', 'Carne, pește, lactate, ouă, soia'),
    'GCG': ('Ala', 'Alanină', 'Carne, pește, lactate, ouă, soia'),
    'TAT': ('Tyr', 'Tirozină', 'Brânză, carne de curcan, pește, semințe susan'),
    'TAC': ('Tyr', 'Tirozină', 'Brânză, carne de curcan, pește, semințe susan'),
    'TAA': ('STOP', 'Codon STOP', 'N/A'),
    'TAG': ('STOP', 'Codon STOP', 'N/A'),
    'CAT': ('His', 'Histidină', 'Carne, pește, nuci, semințe, cereale'),
    'CAC': ('His', 'Histidină', 'Carne, pește, nuci, semințe, cereale'),
    'CAA': ('Gln', 'Glutamină', 'Carne de vită, pui, pește, ouă, lapte'),
    'CAG': ('Gln', 'Glutamină', 'Carne de vită, pui, pește, ouă, lapte'),
    'AAT': ('Asn', 'Asparagină', 'Lactate, carne de vită, pui, ouă, nuci'),
    'AAC': ('Asn', 'Asparagină', 'Lactate, carne de vită, pui, ouă, nuci'),
    'AAA': ('Lys', 'Lizină', 'Carne de porc/pui, parmezan, tofu, fasole neagră'),
    'AAG': ('Lys', 'Lizină', 'Carne de porc/pui, parmezan, tofu, fasole neagră'),
    'GAT': ('Asp', 'Acid Aspartic', 'Năut, carne de vită, ouă, semințe, alune'),
    'GAC': ('Asp', 'Acid Aspartic', 'Năut, carne de vită, ouă, semințe, alune'),
    'GAA': ('Glu', 'Acid Glutamic', 'Roșii, parmezan, ciuperci, soia, carne'),
    'GAG': ('Glu', 'Acid Glutamic', 'Roșii, parmezan, ciuperci, soia, carne'),
    'TGT': ('Cys', 'Cisteină', 'Pui, curcan, usturoi, ceapă, fulgi de ovăz'),
    'TGC': ('Cys', 'Cisteină', 'Pui, curcan, usturoi, ceapă, fulgi de ovăz'),
    'TGA': ('STOP', 'Codon STOP', 'N/A'),
    'TGG': ('Trp', 'Triptofan', 'Curcan, ouă, brânză, ciocolată neagră, ovăz'),
    'CGT': ('Arg', 'Arginină', 'Curcan, semințe de floarea-soarelui, arahide'),
    'CGC': ('Arg', 'Arginină', 'Curcan, semințe de floarea-soarelui, arahide'),
    'CGA': ('Arg', 'Arginină', 'Curcan, semințe de floarea-soarelui, arahide'),
    'CGG': ('Arg', 'Arginină', 'Curcan, semințe de floarea-soarelui, arahide'),
    'AGT': ('Ser', 'Serină', 'Ouă, soia, supă de oase, alune'),
    'AGC': ('Ser', 'Serină', 'Ouă, soia, supă de oase, alune'),
    'AGA': ('Arg', 'Arginină', 'Curcan, semințe de floarea-soarelui, arahide'),
    'AGG': ('Arg', 'Arginină', 'Curcan, semințe de floarea-soarelui, arahide')
}


def analyze_genome_sliding_window(dna_sequence: str):
    """
    Analizeaza genomul folosind o fereastra glisanta de 3 pozitii (pas 1),
    determina Top 10 cei mai frecventi codoni si afiseaza aminoacizii 
    impreuna cu sursele lor alimentare.
    """
    dna_sequence = dna_sequence.strip().upper().replace(" ", "").replace("\n", "").replace("\r", "").replace("\t", "")
    
    if len(dna_sequence) < 3:
        print("Eroare: Secvența introduse este prea scurtă pentru o fereastră de 3 poziții.")
        return

    codons = [dna_sequence[i:i+3] for i in range(len(dna_sequence) - 2)]
    
    codon_counts = Counter(codons)
    total_windows = len(codons)
    
    top_10 = codon_counts.most_common(10)
    
    print("=" * 115)
    print(f"ANALIZĂ FRECVENȚĂ CODONI (SLIDING WINDOW - PAS DE 1)")
    print(f"Lungime genom: {len(dna_sequence)} nucleotide | Total ferestre analizate: {total_windows}")
    print("=" * 115)
    
    header = f"{'Rang':<5} | {'ADN':<5} | {'ARN':<5} | {'Frecvență':<10} | {'Procent':<8} | {'Aminoacid (Nume)':<22} | {'Surse Alimentare (Necesare)'}"
    print(header)
    print("-" * 115)
    
    for rank, (codon, freq) in enumerate(top_10, 1):
        rna_codon = codon.replace('T', 'U')
        
        aa_info = GENETIC_CODE.get(codon, ('???', 'Necunoscut', 'N/A'))
        aa_code, aa_name, foods = aa_info
        
        percentage = (freq / total_windows) * 100
        full_aa_str = f"{aa_code} ({aa_name})"
        
        print(f"{rank:<5} | {codon:<5} | {rna_codon:<5} | {freq:<10} | {percentage:6.2f}%  | {full_aa_str:<22} | {foods}")

    print("=" * 115)


#test
if __name__ == "__main__":
    sample_genome = """
    ATGCGTACTGCTAGCTAGCTAGCTGATCGATCGATCGATCGATCGATCGATCGATCGATC
    GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
    GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
    """
    
    analyze_genome_sliding_window(sample_genome)