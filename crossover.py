# ==========================
# AUXILIARY FUNCTIONS
# ==========================

def count_residues(seq):
    """Conta o número de caracteres diferentes de '-'."""
    return sum(1 for c in seq if c != '-')


def index_at_residue(seq, n):
    """
    Retorna o índice (0-based) da primeira posição
    à direita do n-ésimo resíduo (ignorando gaps).
    """
    count = 0
    for i, c in enumerate(seq):
        if c != '-':
            count += 1
        if count == n:
            return i + 1  # posição depois do n-ésimo resíduo
    return len(seq)


def split_at(align, indexes):
    """
    Divide um MSA (lista de sequências) em duas partes,
    cortando cada linha no índice indicado em 'indexes'.
    Retorna: (left_align, right_align, residues_left)
    """
    left = []
    right = []
    residues = []
    for i, seq in enumerate(align):
        idx = indexes[i]
        left.append(seq[:idx])
        right.append(seq[idx:])
        residues.append(count_residues(seq[:idx]))
    return left, right, residues


def residues_to_indexes(align, residues):
    """
    Para cada sequência no MSA, encontra o índice onde
    há o mesmo número de resíduos (ignorando gaps)
    que na lista 'residues'.
    """
    return [index_at_residue(seq, r) for seq, r in zip(align, residues)]


def pad_alignment(part, side):
    """
    Adiciona gaps ('-') para tornar todas as linhas do MSA do mesmo comprimento.
    side='r' → gaps à direita
    side='l' → gaps à esquerda
    """
    max_len = max(len(seq) for seq in part)
    padded = []
    for seq in part:
        diff = max_len - len(seq)
        if side == 'r':
            padded.append(seq + '-' * diff)
        elif side == 'l':
            padded.append('-' * diff + seq)
        else:
            raise ValueError("side deve ser 'r' (right) ou 'l' (left)")
    return padded


def merge(left, right):
    """Concatena duas partes (listas de strings) coluna a coluna."""
    return [l + r for l, r in zip(left, right)]


def is_all_gaps(align, col_idx):
    """Verifica se a coluna col_idx é composta apenas por gaps ('-')."""
    return all(seq[col_idx] == '-' for seq in align)


def clean_alignment(align):
    """
    Remove colunas compostas apenas por gaps.
    Retorna um novo MSA (lista de strings).
    """
    if not align:
        return []
    cols_to_keep = [i for i in range(len(align[0])) if not is_all_gaps(align, i)]
    cleaned = [''.join(seq[i] for i in cols_to_keep) for seq in align]
    return cleaned

# ==========================
# MAIN CROSSOVER FUNCTION
# ==========================

def generate_offspring(align1, align2, split_residues):
    """
    Executa o crossover entre dois MSAs alinhados.
    split_residues: número de resíduos (não colunas) após o qual cortar.
    """
    nseqs = len(align1)

    # 1. Cortar o primeiro MSA
    l1, r1, residues_left = split_at(align1, [split_residues] * nseqs)

    # 2. Determinar índices equivalentes no segundo MSA
    indexes2 = residues_to_indexes(align2, residues_left)

    # 3. Cortar o segundo MSA
    l2, r2, _ = split_at(align2, indexes2)

    # 4. Fazer padding das partes
    l1pad = pad_alignment(l1, "r")
    r1pad = pad_alignment(r1, "l")
    l2pad = pad_alignment(l2, "r")
    r2pad = pad_alignment(r2, "l")

    # 5. Fazer crossover
    offspring1 = clean_alignment(merge(l1pad, r2pad))
    offspring2 = clean_alignment(merge(l2pad, r1pad))

    return offspring1, offspring2
