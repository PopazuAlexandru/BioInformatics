import random
import urllib.request
import numpy as np
import matplotlib.pyplot as plt

def fetch_ncbi_sequence(accession: str = "NC_000913.3", target_len: int = 2000) -> str:
    """Descarca o secventa ADN de pe NCBI de lungime intre 1000 si 3000 nt."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    print(f"Descarcare secventa NCBI ({accession})...")
    
    with urllib.request.urlopen(url) as response:
        lines = response.read().decode('utf-8').splitlines()
        
    sequence = "".join([line.strip().upper() for line in lines if not line.startswith(">")])
    return sequence[:target_len]


def generate_dna_samples(sequence: str, num_samples: int = 10, min_len: int = 100, max_len: int = 2000) -> np.ndarray:
    """Preleveaza N probe aleatoare din secventa principala si le stocheaza intr-un array."""
    max_len = min(max_len, len(sequence))
    samples = []
    
    for _ in range(num_samples):
        sample_len = random.randint(min_len, max_len)
        max_start = len(sequence) - sample_len
        start_idx = random.randint(0, max_start)
        sample_fragment = sequence[start_idx : start_idx + sample_len]
        samples.append(sample_fragment)
        
    return np.array(samples, dtype=object)


def simulate_gel_electrophoresis(samples_array: np.ndarray, min_bp: int = 100, max_bp: int = 2000):
    """
    Simuleaza migrarea fragmentelor pe gelul de electroforeza.
    Formula logaritmica: Distanta de migrare este invers proportionala cu log10(lungime_bp).
    """
    lengths = np.array([len(fragment) for fragment in samples_array])
    
    log_min = np.log10(min_bp)
    log_max = np.log10(max_bp)
    
    migration_distances = 1.0 - (np.log10(lengths) - log_min) / (log_max - log_min)
    gel_positions = 0.15 + migration_distances * 0.75

    fig, ax = plt.subplots(figsize=(6, 8), facecolor='#1e1e1e')
    ax.set_facecolor('#121212')
    
    well = plt.Rectangle((0.35, 0.92), 0.3, 0.04, color='#333333', ec='white')
    ax.add_patch(well)
    ax.text(0.5, 0.94, "Godeu (Anod -)", color='white', ha='center', va='center', fontsize=10, fontweight='bold')
    
    for idx, (bp_size, pos) in enumerate(zip(lengths, gel_positions), 1):
        band = plt.Rectangle((0.3, pos), 0.4, 0.008, color='#00ff66', alpha=0.85)
        ax.add_patch(band)
        ax.text(0.72, pos + 0.002, f"Proba {idx}: {bp_size} bp", color='#00ff66', fontsize=9, va='center')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.invert_yaxis()  
    ax.set_xticks([])
    ax.set_yticks([0.15, 0.5, 0.9])
    ax.set_yticklabels(['Lungi (Încet)', 'Medii', 'Scurte (Rapid)'], color='white')
    ax.set_title("Simulare Gel Electroforeză ADN (Laborator 5)", color='white', fontsize=12, pad=15)
    
    plt.tight_layout()
    plt.show()

    print("\n--- REZULTATE SAMPLES & SIMULARE MIGRARE ---")
    print(f"{'Proba #':<8} | {'Lungime Fragment (bp)':<22} | {'Distanta Migrare Relativă'}")
    print("-" * 60)
    for idx, (bp, dist) in enumerate(zip(lengths, migration_distances), 1):
        print(f"Proba {idx:<3} | {bp:<22} bp | {dist * 100:6.2f}% din gel")


if __name__ == "__main__":
    dna_seq = fetch_ncbi_sequence("NC_000913.3", target_len=2000)
    samples = generate_dna_samples(dna_seq, num_samples=10, min_len=100, max_len=2000)
    simulate_gel_electrophoresis(samples)