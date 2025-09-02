# ==============================================================================
# TIE GENOME ENGINE v3.0 - ANÁLISIS Y PRUEBA CUANTITATIVA FINAL
# ==============================================================================

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import urllib.request
import os
from scipy.stats import pearsonr

# --- 1. PREPARACIÓN DEL GENOMA ---
def preparar_genoma_sars_cov_2():
    print("🌐 Descargando y preparando el genoma del SARS-CoV-2...")
    url_fasta = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1798174254"
    fasta_filename = "sars_cov_2.fasta"
    try:
        if not os.path.exists(fasta_filename):
            urllib.request.urlretrieve(url_fasta, fasta_filename)
        with open(fasta_filename, 'r') as f:
            sequence = ''.join([line.strip() for line in f if not line.startswith('>') and line.strip()])
        valid_bases = "ACTG"
        cleaned_sequence = "".join([base for base in sequence.upper() if base in valid_bases])
        print(f"🧬 Genoma cargado y listo: {len(cleaned_sequence):,} bases.")
        return cleaned_sequence
    except Exception as e:
        print(f"❌ Error al cargar el genoma: {e}")
        return ""

# --- 2. MOTOR DE ANÁLISIS ---
def generar_micro_instrucciones(codon):
    instrucciones = {}
    for base_programa in ['A', 'C', 'T', 'G']:
        bits = ''.join(['1' if base_codon == base_programa else '0' for base_codon in codon])
        instrucciones[base_programa] = bits
    return instrucciones

def calcular_entropia(lista_de_opcodes):
    if not lista_de_opcodes: return 0
    counts = Counter(lista_de_opcodes)
    total = len(lista_de_opcodes)
    probs = [count / total for count in counts.values()]
    return -sum(p * np.log2(p) for p in probs)

def analizar_genoma_completo(genome_clean, window_size=100):
    print("🗺️  Mapeando el microcódigo y analizando...")
    microcode_map = []
    for i in range(0, len(genome_clean) - 2, 3):
        microcode_map.append(generar_micro_instrucciones(genome_clean[i:i+3]))
    
    scaffold_stability = []
    info_density_acg = []
    positions = []
    for i in range(0, len(microcode_map) - window_size, 10):
        window = microcode_map[i:i+window_size]
        positions.append(i * 3)
        t_activity = sum(instruction['T'].count('1') for instruction in window)
        scaffold_stability.append(t_activity)
        opcodes_acg = [f"A_{instr['A']}" for instr in window] + [f"C_{instr['C']}" for instr in window] + [f"G_{instr['G']}" for instr in window]
        info_density_acg.append(calcular_entropia(opcodes_acg))
        
    return np.array(positions), np.array(scaffold_stability), np.array(info_density_acg)

# --- 3. CUANTIFICACIÓN Y VISUALIZACIÓN FINAL ---
def cuantificar_y_visualizar_final(positions, stability, density, genome_length):
    print("🔬 Cuantificando la correlación inversa...")
    min_len = min(len(stability), len(density))
    correlation, p_value = pearsonr(stability[:min_len], density[:min_len])
    
    print(f"\n=============================================================")
    print(f"  COEFICIENTE DE CORRELACIÓN DE PEARSON: {correlation:.4f}")
    print(f"  P-VALUE: {p_value:.2e}")
    print(f"=============================================================")
    if correlation < -0.5 and p_value < 0.001:
        print("✅ INTERPRETACIÓN: Correlación inversa fuerte y significativa. La hipótesis se comprueba cuantitativamente.")
    else:
        print("⚠️  Correlación no tan fuerte como se esperaba. Revisar parámetros o modelo.")

    print("\n📊 Generando la visualización final para publicación...")
    fig, ax1 = plt.subplots(figsize=(24, 10))
    color_stability, color_density = 'darkgray', '#c9184a'
    ax1.set_xlabel('Posición en el Genoma (bases)', fontsize=16)
    ax1.set_ylabel('Estabilidad del Andamio-T (Infraestructura)', color=color_stability, fontsize=16)
    ax1.plot(positions, stability, color=color_stability, alpha=0.8, linewidth=2.5, label='Andamio-T (Estructura Rígida)')
    ax1.fill_between(positions, stability, color=color_stability, alpha=0.1)
    ax1.tick_params(axis='y', labelcolor=color_stability, labelsize=12)
    ax1.grid(True, linestyle=':', which='both', alpha=0.3)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Densidad de Información ACG (Contenido Funcional)', color=color_density, fontsize=16)
    ax2.plot(positions, density, color=color_density, linewidth=2.5, label='Información ACG (Contenido Funcional)')
    ax2.tick_params(axis='y', labelcolor=color_density, labelsize=12)
    title_text = f'Correlación Inversa entre Andamio e Información en SARS-CoV-2\n(Coeficiente de Correlación: {correlation:.3f})'
    plt.title(title_text, fontsize=22, fontweight='bold', pad=20)
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper center', fontsize=12, bbox_to_anchor=(0.5, -0.12), fancybox=True, shadow=True, ncol=2)
    plt.xlim(0, genome_length)
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plt.savefig('TIE_Final_Proof.png', dpi=300)
    plt.show()

# --- 4. EJECUCIÓN PRINCIPAL ---
if __name__ == "__main__":
    genome_clean = preparar_genoma_sars_cov_2()
    if genome_clean:
        positions, stability, density = analizar_genoma_completo(genome_clean)
        cuantificar_y_visualizar_final(positions, stability, density, len(genome_clean))
