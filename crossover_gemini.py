from multiplealign.myalign import MyAlign

def count_residues(sequence):
    """
    Conta o número de resíduos (caracteres que não são '-') numa sequência.
    Exemplo: count_residues("--AC--") -> 2
    """
    return sum(1 for char in sequence if char != '-')

def index_at_residue(sequence, num_residues_left):
    """
    Encontra o índice da primeira posição do lado direito para ter um número específico
    de resíduos no lado esquerdo.
    Exemplo: index_at_residue("--AC--GT", 2) -> 4 (o 'G' está no índice 4)
    """
    current_residues = 0
    for i, char in enumerate(sequence):
        if char != '-':
            current_residues += 1
        if current_residues == num_residues_left:
            return i + 1  # O split ocorre APÓS este índice
    return len(sequence) # Se todos os resíduos estiverem no lado esquerdo

def split_at(alignment, split_indexes):
    """
    Divide um alinhamento em duas partes (esquerda e direita) com base nos índices fornecidos.
    Retorna a parte esquerda, a parte direita e uma lista com o número de resíduos
    na parte esquerda de cada sequência.
    """
    left_part_sequences = []
    right_part_sequences = []
    left_residues_counts = []

    for i in range(alignment.num_seqs()): # Usar num_seqs()
        seq = alignment[i] # Usar __getitem__ para obter a sequência
        split_idx = split_indexes[i]
        left_seq = seq[:split_idx]
        right_seq = seq[split_idx:]

        left_part_sequences.append(left_seq)
        right_part_sequences.append(right_seq)
        left_residues_counts.append(count_residues(left_seq))

    left_alignment = MyAlign(left_part_sequences, alignment.al_type)
    right_alignment = MyAlign(right_part_sequences, alignment.al_type)

    return left_alignment, right_alignment, left_residues_counts

def residues_to_indexes(alignment, target_residue_counts):
    """
    Dado um alinhamento e uma lista de contagens de resíduos desejadas para o lado esquerdo
    de cada sequência, retorna os índices de divisão correspondentes.
    """
    split_indexes = []
    for i in range(alignment.num_seqs()): # Usar num_seqs()
        seq = alignment[i] # Usar __getitem__ para obter a sequência
        split_indexes.append(index_at_residue(seq, target_residue_counts[i]))
    return split_indexes

def pad_alignment(alignment, pad_side):
    """
    Preenche um alinhamento com gaps ('-') para que todas as sequências tenham
    o mesmo comprimento (o comprimento da mais longa).
    'pad_side' pode ser 'l' para preencher à esquerda ou 'r' para preencher à direita.
    """
    if not alignment.listseqs: # Usar listseqs para verificar se está vazio
        return MyAlign([])

    max_len = max(len(seq) for seq in alignment.listseqs)
    padded_sequences = []

    for seq in alignment.listseqs: # Usar listseqs
        current_len = len(seq)
        padding_needed = max_len - current_len
        if pad_side == 'r':
            padded_sequences.append(seq + '-' * padding_needed)
        elif pad_side == 'l':
            padded_sequences.append('-' * padding_needed + seq)
        else:
            raise ValueError("pad_side deve ser 'l' (left) ou 'r' (right)")
    return MyAlign(padded_sequences, alignment.al_type)

def merge(left_alignment, right_alignment):
    """
    Une duas partes de alinhamentos (esquerda e direita) numa só.
    As sequências são unidas elemento a elemento.
    """
    merged_sequences = []
    for i in range(left_alignment.num_seqs()): # Usar num_seqs()
        merged_sequences.append(left_alignment[i] + right_alignment[i]) # Usar __getitem__
    return MyAlign(merged_sequences, left_alignment.al_type)

def is_all_gaps(alignment, column_index):
    """
    Verifica se uma coluna num alinhamento é composta apenas por gaps ('-').
    """
    if not alignment.listseqs: # Usar listseqs
        return True
    return all(seq[column_index] == '-' for seq in alignment.listseqs) # Usar listseqs

def clean_alignment(alignment):
    """
    Remove colunas de um alinhamento que são compostas apenas por gaps.
    """
    if not alignment.listseqs: # Usar listseqs
        return MyAlign([])
    
    columns_to_keep_indexes = []
    for col_idx in range(len(alignment)): # Usar len(alignment)
        if not is_all_gaps(alignment, col_idx):
            columns_to_keep_indexes.append(col_idx)

    cleaned_sequences = []
    for seq_idx in range(alignment.num_seqs()): # Usar num_seqs()
        new_seq = "".join(alignment[seq_idx][col_idx] for col_idx in columns_to_keep_indexes) # Usar __getitem__
        cleaned_sequences.append(new_seq)
    
    return MyAlign(cleaned_sequences, alignment.al_type)

# --- Função Principal de Crossover ---

def generate_offspring(parent1_alignment, parent2_alignment, split_column_index):
    """
    Executa o operador de crossover entre dois alinhamentos parentais para gerar dois descendentes.

    parent1_alignment: Objeto MyAlign do primeiro parental.
    parent2_alignment: Objeto MyAlign do segundo parental.
    split_column_index: O índice da coluna onde o primeiro parental será dividido.
    """
    # 1. Dividir o primeiro parental (align1) na coluna especificada
    p1_split_indexes = [split_column_index] * parent1_alignment.num_seqs() # Usar num_seqs()
    p1_left, p1_right, p1_left_residue_counts = split_at(parent1_alignment, p1_split_indexes)

    # 2. Determinar os índices de divisão para o segundo parental (align2)
    p2_split_indexes = residues_to_indexes(parent2_alignment, p1_left_residue_counts)
    p2_left, p2_right, _ = split_at(parent2_alignment, p2_split_indexes)

    # 3. Preencher as partes para que possam ser unidas
    
    # Determinar o comprimento máximo de cada parte para o preenchimento
    max_len_p1_left = max(len(s) for s in p1_left.listseqs) if p1_left.listseqs else 0
    max_len_p2_left = max(len(s) for s in p2_left.listseqs) if p2_left.listseqs else 0
    max_len_left_parts = max(max_len_p1_left, max_len_p2_left)

    max_len_p1_right = max(len(s) for s in p1_right.listseqs) if p1_right.listseqs else 0
    max_len_p2_right = max(len(s) for s in p2_right.listseqs) if p2_right.listseqs else 0
    max_len_right_parts = max(max_len_p1_right, max_len_p2_right)

    # Preencher a parte esquerda de P1 e P2 para o comprimento máximo das partes esquerdas
    padded_p1_left_seqs = [s + '-' * (max_len_left_parts - len(s)) for s in p1_left.listseqs]
    padded_p2_left_seqs = [s + '-' * (max_len_left_parts - len(s)) for s in p2_left.listseqs]
    padded_p1_left = MyAlign(padded_p1_left_seqs)
    padded_p2_left = MyAlign(padded_p2_left_seqs)

    # Preencher a parte direita de P1 e P2 para o comprimento máximo das partes direitas
    padded_p1_right_seqs = ['-' * (max_len_right_parts - len(s)) + s for s in p1_right.listseqs]
    padded_p2_right_seqs = ['-' * (max_len_right_parts - len(s)) + s for s in p2_right.listseqs]
    padded_p1_right = MyAlign(padded_p1_right_seqs)
    padded_p2_right = MyAlign(padded_p2_right_seqs)

    # 4. Criar os descendentes unindo as partes
    # Offspring 1: P1_left + P2_right
    offspring1_raw = merge(padded_p1_left, padded_p2_right)
    # Offspring 2: P2_left + P1_right
    offspring2_raw = merge(padded_p2_left, padded_p1_right)

    # 5. Limpar os alinhamentos resultantes (remover colunas apenas de gaps)
    offspring1 = clean_alignment(offspring1_raw)
    offspring2 = clean_alignment(offspring2_raw)

    return offspring1.listseqs, offspring2.listseqs

# --- Teste com os exemplos do guia ---
print("--- Exemplos Iniciais ---")
align1 = MyAlign(["ATC-G-G-TT","A-CCGC-ATC","A-CCG--ATC"],"PROTEIN") # Ajustado: lista de strings diretamente
align2 = MyAlign(["ATCGG---TT","AC-CGCA-TC","A---CCGATC"],"PROTEIN") # Ajustado: lista de strings diretamente

print("Parental 1:")
print(align1)
print("\nParental 2:")
print(align2)

of1, of2 = generate_offspring(align1, align2, 3)

print("\nOffspring 1 (esperado ATCGG---TT, A-CCGCA-TC, A-C--CGATC):")
print(of1)
print("\nOffspring 2 (esperado ATC--G-G-TT, AC--CGC-ATC, A--CCG--ATC):")
print(of2)

# --- Teste das Funções Auxiliares Individualmente ---
print("\n--- Testes de Funções Auxiliares ---")
print("count_residues('--AC--'):", count_residues("--AC--")) # Out: 2
print("index_at_residue('--AC--GT', 2):", index_at_residue("--AC--GT", 2)) # Out: 4

print("\nAlinhamento 1:")
print(align1)
l1, r1, residues = split_at(align1, [3]*align1.num_seqs()) # Usar num_seqs()
print("\nSplit_at (l1):")
print(l1) # Out: ATC, A-C, A-C
print("\nSplit_at (r1):")
print(r1) # Out: -G-G-TT, CGC-ATC, CG--ATC
print("\nResidues left (l1):", residues) # Out: [3, 2, 2]

print("\nAlinhamento 2:")
print(align2)
indexes = residues_to_indexes(align2, residues)
print("\nIndexes for align2 based on align1 residues:", indexes) # Out: [3, 2, 5]

l2, r2, _ = split_at(align2, indexes)
print("\nSplit_at (l2):")
print(l2) # Out: ATC, AC, A---C
print("\nSplit_at (r2):")
print(r2) # Out: GG---TT, -CGCA-TC, CGATC

print("\nPadding teste (l2 pad right):")
l2pad = pad_alignment(l2, "r")
print(l2pad) # Out: ATC--, AC---, A---C

print("\nPadding teste (r2 pad left):")
r2pad = pad_alignment(r2, "l")
print(r2pad) # Out: -GG---TT, -CGCA-TC, ---CGATC

print("\nMerge teste (l1, r2pad):")
l1r2pad_merged = merge(l1, r2pad)
print(l1r2pad_merged) # Out: ATC-GG---TT, A-C-CGCA-TC, A-C---CGATC

print("\nis_all_gaps(l1r2pad_merged, 2):", is_all_gaps(l1r2pad_merged, 2)) # False
print("is_all_gaps(l1r2pad_merged, 3):", is_all_gaps(l1r2pad_merged, 3)) # True

print("\nClean_alignment(l1r2pad_merged):")
cleaned_alignment = clean_alignment(l1r2pad_merged)
print(cleaned_alignment) # Out: ATCGG---TT, A-CCGCA-TC, A-C--CGATC