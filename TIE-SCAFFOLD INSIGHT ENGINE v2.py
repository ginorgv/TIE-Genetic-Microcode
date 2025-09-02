# ==============================================================================
# TIE-SCAFFOLD INSIGHT ENGINE v2.0 - NOTEBOOK AUTÓNOMO
# Descripción: Aísla y ejecuta el "Programa-T" del genoma del SARS-CoV-2,
#              simula su andamio estructural y genera un análisis visual
#              avanzado de su dinámica computacional.
# Autor: Raúl Vázquez Verde (Teoría) y Asistente AI (Implementación)
# ==============================================================================

# --- 1. IMPORTACIÓN DE LIBRERÍAS ---
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import urllib.request
import gzip
import os

print("✅ Librerías importadas.")

# --- 2. DESCARGA Y PREPARACIÓN DEL GENOMA ---

def preparar_genoma_sars_cov_2():
    """Descarga, descomprime y lee el genoma del SARS-CoV-2."""
    print("\n🌐 Descargando y preparando el genoma del SARS-CoV-2...")
    
    url_fasta = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1798174254"
    fasta_filename = "sars_cov_2.fasta"

    try:
        urllib.request.urlretrieve(url_fasta, fasta_filename)
        print("✅ Genoma FASTA descargado.")
        
        with open(fasta_filename, 'r') as f:
            sequence = ''.join([line.strip() for line in f if not line.startswith('>')])
        
        print(f"🧬 Genoma cargado y listo: {len(sequence):,} bases.")
        return sequence.upper()
        
    except Exception as e:
        print(f"❌ Error durante la descarga o lectura: {e}")
        print("⚠️  Usando una secuencia de ejemplo para la demostración.")
        return ("ATG" * 1000 + "GGG" * 500 + "TAA" * 200) * 10

# --- 3. MOTOR DE SIMULACIÓN DEL PROGRAMA-T ---

def generar_micro_instruccion_T(codon):
    """Extrae la micro-instrucción de 3 bits solo para el Programa-T."""
    return ''.join(['1' if base == 'T' else '0' for base in codon])

class T_Scaffold_Simulator:
    """Ejecuta el Programa-T de un genoma y simula el andamio estructural."""
    
    def __init__(self, genome_sequence):
        if not genome_sequence:
            raise ValueError("La secuencia del genoma no puede estar vacía.")
        self.genome = genome_sequence
        self.t_program_sequence = self._extract_t_program()
        
        self.T_ISA = {
            '000': {'action': 'wait', 'delta': 0},
            '001': {'action': 'extend_weak', 'delta': 1},
            '010': {'action': 'extend_weak', 'delta': 1},
            '011': {'action': 'compress_strong', 'delta': -2},
            '100': {'action': 'break', 'delta': 'reset'},
            '101': {'action': 'init', 'delta': 3},
            '110': {'action': 'extend_strong', 'delta': 2},
            '111': {'action': 'max_scaffold', 'delta': 5},
        }
    
    def _extract_t_program(self):
        """Extrae la secuencia completa de opcodes para el Programa-T."""
        print("\n extracting T-Program instruction sequence...")
        program_seq = []
        for i in range(0, len(self.genome) - 2, 3):
            codon = self.genome[i:i+3]
            program_seq.append(generar_micro_instruccion_T(codon))
        print(f"✅ Programa-T extraído: {len(program_seq):,} instrucciones.")
        return program_seq
        
    def simulate_scaffold_formation(self):
        """Ejecuta el programa y simula la formación del andamio 1D."""
        print("🏗️  Simulando formación del andamio estructural...")
        
        scaffold_y = [0]
        current_y = 0
        
        for opcode in self.t_program_sequence:
            instruction = self.T_ISA.get(opcode)
            
            if instruction['delta'] == 'reset':
                current_y = 0
            else:
                current_y += instruction['delta']
            
            scaffold_y.append(current_y)
            
        print("✅ Simulación completada.")
        return scaffold_y
        
    def visualize_advanced_analysis(self, scaffold_data):
        """
        Visualiza el andamio de forma legible, mostrando Energía y Volatilidad.
        """
        print("📊 Generando visualización AVANZADA y LEGIBLE del andamio...")
        
        positions = np.arange(len(scaffold_data)) * 3
        
        # --- CÁLCULO DE MÉTRICAS LEGIBLES ---
        window_size = 150  # Ventana de 150 codones (450 bases)
        
        # 1. Energía Estructural (Suavizado): Muestra la tendencia general de la actividad
        energia_estructural = np.convolve(np.abs(scaffold_data), np.ones(window_size)/window_size, mode='same')
        
        # 2. Volatilidad Estructural (Derivada): Muestra los cambios abruptos o "terremotos"
        volatilidad = np.abs(np.diff(scaffold_data, prepend=0))
        volatilidad_suavizada = np.convolve(volatilidad, np.ones(window_size)/window_size, mode='same')

        # --- CREACIÓN DEL GRÁFICO MEJORADO ---
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(22, 12), sharex=True, 
                                     gridspec_kw={'height_ratios': [3, 2]})
        
        fig.suptitle('Análisis Dinámico del Programa-T del SARS-CoV-2 (Andamiaje Estructural)', 
                     fontsize=20, fontweight='bold')

        # GRÁFICO 1: ENERGÍA ESTRUCTURAL
        ax1.plot(positions, energia_estructural, color='#0077b6', label='Energía Estructural (Actividad Acumulada)')
        ax1.fill_between(positions, energia_estructural, color='#0077b6', alpha=0.1)
        ax1.set_title('Fase 1: Energía Estructural - ¿Dónde construye el virus su andamio?', fontsize=16)
        ax1.set_ylabel('Energía Estructural (Computada)', color='#0077b6', fontsize=14)
        ax1.tick_params(axis='y', labelcolor='#0077b6')
        ax1.grid(True, linestyle=':', alpha=0.7)
        
        # Superponer los genes para contexto
        gene_boundaries = {
            'ORF1ab': (266, 21555, '#a2d2ff'), 
            'S (Spike)': (21563, 25384, '#ffb703'), 
            'E (Envelope)': (26245, 26472, '#ef476f'),
            'M (Membrane)': (26523, 27191, '#06d6a0'), 
            'N (Nucleocapsid)': (28274, 29533, '#7209b7')
        }
        for gene, (start, end, color) in gene_boundaries.items():
            ax1.axvspan(start, end, alpha=0.2, label=f'Gen {gene}', color=color)
        ax1.legend(loc='upper left')

        # GRÁFICO 2: VOLATILIDAD ESTRUCTURAL
        ax2.plot(positions, volatilidad_suavizada, color='#d00000', label='Volatilidad Estructural (Cambios Abruptos)')
        ax2.fill_between(positions, volatilidad_suavizada, color='#d00000', alpha=0.1)
        ax2.set_title('Fase 2: Volatilidad Estructural - ¿Dónde están los límites y comandos críticos?', fontsize=16)
        ax2.set_ylabel('Volatilidad (Computada)', color='#d00000', fontsize=14)
        ax2.tick_params(axis='y', labelcolor='#d00000')
        ax2.set_xlabel('Posición en el Genoma (bases)', fontsize=14)
        ax2.grid(True, linestyle=':', alpha=0.7)
        
        # Resaltar picos de volatilidad, que deberían ser los límites de los genes
        umbral_volatilidad = np.percentile(volatilidad_suavizada, 98) # Top 2% de cambios
        ax2.axhline(y=umbral_volatilidad, color='black', linestyle='--', label=f'Umbral de "Terremoto" Estructural (Top 2%)')
        picos_volatilidad = [pos for pos, vol in zip(positions, volatilidad_suavizada) if vol > umbral_volatilidad]
        ax2.scatter(picos_volatilidad, [umbral_volatilidad]*len(picos_volatilidad), color='black', zorder=5)
        
        ax2.legend(loc='upper left')
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig('TIE_scaffold_advanced_analysis.png', dpi=300)
        plt.show()

# --- 4. FUNCIÓN PRINCIPAL DE EJECUCIÓN ---

def main():
    """Función principal que ejecuta todo el pipeline de análisis."""
    
    # PASO 1: Obtener la secuencia del genoma
    genome_clean = preparar_genoma_sars_cov_2()
    
    if genome_clean:
        # PASO 2: Inicializar el simulador
        t_simulator = T_Scaffold_Simulator(genome_clean)
        
        # PASO 3: Simular la formación del andamio
        scaffold_structure = t_simulator.simulate_scaffold_formation()
        
        # PASO 4: Generar la visualización avanzada
        t_simulator.visualize_advanced_analysis(scaffold_structure)
        
        print("\n🎉 Experimento completado con éxito.")
        print("Se ha generado el archivo 'TIE_scaffold_advanced_analysis.png' con el análisis legible.")
        print("\n🎯 Próximos Pasos:")
        print("1. Analiza el gráfico de Volatilidad (rojo): ¿Los picos negros coinciden con los bordes de color de los genes?")
        print("2. Analiza el gráfico de Energía (azul): ¿La 'montaña' de actividad se correlaciona con el gen ORF1ab?")
        print("3. Si la respuesta a ambas es SÍ, tienes una validación visual directa de la TIE.")

# --- EJECUTAR EL SCRIPT COMPLETO ---
if __name__ == "__main__":
    main()
