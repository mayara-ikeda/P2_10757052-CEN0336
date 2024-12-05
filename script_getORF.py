import sys

# Dicionário do código genético: mapeia códons para seus respectivos aminoácidos
GENETIC_CODE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

# Códon de início e códons de parada
START_CODON = 'ATG'
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

# Função para gerar o complemento reverso de uma sequência de DNA
def reverse_complement(seq):
    complement = str.maketrans('ATCG', 'TAGC')  # Mapeia nucleotídeos para seus complementos
    return seq.translate(complement)[::-1]  # Inverte a sequência após complementar

# Função para encontrar ORFs (Open Reading Frames) em uma sequência
def find_orfs(sequence, frame_offset):
    orfs = []
    # Percorre a sequência na fase de leitura especificada (frame_offset)
    for start in range(frame_offset, len(sequence) - 2, 3):
        codon = sequence[start:start + 3]
        if codon == START_CODON:  # Verifica se o códon é um códon de início
            for end in range(start + 3, len(sequence) - 2, 3):
                stop_codon = sequence[end:end + 3]
                if stop_codon in STOP_CODONS:  # Verifica se encontrou um códon de parada
                    if (end + 3 - start) % 3 == 0:  # Verifica se o comprimento é múltiplo de 3
                        orfs.append((start, end + 3))  # Armazena o ORF (start e end)
                    break  # Encerra a busca no ORF atual
    return orfs

# Função para traduzir sequência de DNA para sequência de peptídeos
def translate_sequence(seq):
    return ''.join(GENETIC_CODE[seq[i:i + 3]] for i in range(0, len(seq), 3) if seq[i:i + 3] in GENETIC_CODE)

# Função para processar o arquivo .fasta e retornar um dicionário com os identificadores e sequências
def process_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = {}
        current_header = None
        current_seq = []
        for line in file:
            if line.startswith('>'):  # Identifica o início de um novo registro
                if current_header:  # Salva o registro anterior
                    sequences[current_header] = ''.join(current_seq)
                current_header = line.strip()  # Salva o cabeçalho
                current_seq = []  # Reinicia a sequência
            else:
                current_seq.append(line.strip().upper())  # Adiciona a sequência (em letras maiúsculas)
        if current_header:  # Salva o último registro
            sequences[current_header] = ''.join(current_seq)
    return sequences

# Função principal
def main():
    # Verifica se o arquivo de entrada foi especificado
    if len(sys.argv) != 2:
        print("Uso: python script_getORF.py <arquivo multifasta>")
        sys.exit(1)

    input_file = sys.argv[1]  # Lê o nome do arquivo de entrada
    sequences = process_fasta(input_file)  # Processa o arquivo FASTA

    # Abre os arquivos de saída
    with open('ORF.fna', 'w') as orf_fna, open('ORF.faa', 'w') as orf_faa:
        for header, sequence in sequences.items():
            best_orf = (None, None, 0, None)  # Armazena o melhor ORF encontrado (start, end, length, frame)
            # Verifica todas as 6 fases de leitura (3 diretas + 3 reversas)
            for frame in range(3):
                orfs = find_orfs(sequence, frame)  # ORFs na sequência direta
                orfs += find_orfs(reverse_complement(sequence), frame)  # ORFs na sequência reversa
                for start, end in orfs:
                    length = end - start
                    if length > best_orf[2]:  # Atualiza o melhor ORF se o novo for mais longo
                        best_orf = (start, end, length, frame + 1 if start < len(sequence) else frame + 4)

            start, end, length, frame = best_orf
            if start is not None:  # Garante que encontrou um ORF
                # Obtém a sequência correspondente ao melhor ORF
                orf_seq = sequence[start:end] if frame <= 3 else reverse_complement(sequence)[start:end]
                peptide = translate_sequence(orf_seq)  # Traduz o ORF para peptídeo
                coords = f"frame{frame}_{start + 1}_{end}"  # Formata as coordenadas para o identificador
                # Escreve nos arquivos de saída
                orf_fna.write(f"{header}_{coords}\n{orf_seq}\n")
                orf_faa.write(f"{header}_{coords}\n{peptide}\n")

# Executa a função principal apenas se o script for executado diretamente
if __name__ == "__main__":
    main()
