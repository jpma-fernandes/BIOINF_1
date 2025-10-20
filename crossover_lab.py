###cross2

def count_residues(seq):
    return sum(1 for char in seq if char != '-')


def index_at_residue(seq, num_residues):
    count = 0
    for i, char in enumerate(seq):
        if char != '-':
            count += 1
            if count == num_residues:
                return i + 1
    return len(seq)


def split_at(align, indexes):
    left_part = []
    right_part = []
    residues_left = []
    
    for i, seq in enumerate(align):
        index = indexes[i]
        left = seq[:index]
        right = seq[index:]
        
        left_part.append(left)
        right_part.append(right)
        residues_left.append(count_residues(left))
    
    return left_part, right_part, residues_left


def residues_to_indexes(align, residues):
    indexes = []
    for i, seq in enumerate(align):
        target_residues = residues[i]
        index = index_at_residue(seq, target_residues)
        indexes.append(index)
    
    return indexes


def pad_alignment(align, side):
    if not align:
        return align
    
    max_len = max(len(seq) for seq in align)
    padded = []
    
    for seq in align:
        gap_needed = max_len - len(seq)
        if side == "r":
            padded.append(seq + '-' * gap_needed)
        elif side == "l":
            padded.append('-' * gap_needed + seq)
        else:
            padded.append(seq)
    
    return padded


def merge(left, right):
    if len(left) != len(right):
        raise ValueError("As duas partes devem ter o mesmo número de sequências")
    
    merged = []
    for i in range(len(left)):
        merged.append(left[i] + right[i])
    
    return merged


def is_all_gaps(align, col):
    if not align or col >= len(align[0]):
        return False
    
    for seq in align:
        if col < len(seq) and seq[col] != '-':
            return False
    
    return True


def clean_alignment(align):
    if not align or len(align[0]) == 0:
        return align
    
    cols_to_keep = []
    for col in range(len(align[0])):
        if not is_all_gaps(align, col):
            cols_to_keep.append(col)
    
    cleaned = []
    for seq in align:
        new_seq = ''.join(seq[col] for col in cols_to_keep)
        cleaned.append(new_seq)
    
    return cleaned


def generate_offspring(align1, align2, crossover_point):
    import random
    
    if len(align1) != len(align2):
        return align1, align2
    
    if not align1:
        return align1, align2
    
    if crossover_point is None:
        max_len = len(align1[0])
        if max_len <= 1:
            return align1, align2
        crossover_point = random.randint(1, max_len - 1)
    
    indexes1 = [crossover_point] * len(align1)
    l1, r1, residues = split_at(align1, indexes1)
    
    indexes2 = residues_to_indexes(align2, residues)
    
    l2, r2, _ = split_at(align2, indexes2)
    
    l2_padded = pad_alignment(l2, "r")
    r2_padded = pad_alignment(r2, "l")
    
    offspring1 = merge(l1, r2_padded)
    offspring2 = merge(l2_padded, r1)
    
    offspring1 = clean_alignment(offspring1)
    offspring2 = clean_alignment(offspring2)
    
    return offspring1, offspring2