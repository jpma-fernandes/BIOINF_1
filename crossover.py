def get_residues_only(seq):
    return seq.replace('-', '')


def extract_offsets_and_residues(align):

    offsets = []
    residues_list = []
    
    for seq in align:
        #* Contar gaps iniciais
        offset = 0
        for c in seq:
            if c == '-':
                offset += 1
            else:
                break
        
        offsets.append(offset)
        residues_list.append(get_residues_only(seq))
    
    return offsets, residues_list


def reconstruct_from_offsets(offsets, residues_list):

    #* Criar sequências com offsets
    seqs = []
    for offset, residues in zip(offsets, residues_list):
        seq = '-' * offset + residues
        seqs.append(seq)
    
    #* Normalizar comprimentos
    max_len = max(len(s) for s in seqs)
    seqs = [s + '-' * (max_len - len(s)) for s in seqs]
    
    return remove_gap_only_columns(seqs)


def remove_gap_only_columns(align):
    if not align or len(align[0]) == 0:
        return align
    
    cols_to_keep = []
    for col_idx in range(len(align[0])):
        #* Manter coluna se pelo menos uma sequência tem resíduo
        if any(seq[col_idx] != '-' for seq in align):
            cols_to_keep.append(col_idx)
    
    #* Voltamos a contruir  alinhamento
    result = []
    for seq in align:
        new_seq = ''.join(seq[i] for i in cols_to_keep)
        result.append(new_seq)
    
    return result


def generate_offspring(align1, align2, crossover_point=None):

    import random
    
    offsets1, residues1 = extract_offsets_and_residues(align1)
    offsets2, residues2 = extract_offsets_and_residues(align2)
    
    # Verificar se resíduos são iguais
    if residues1 != residues2:
        # Se não são, não podemos fazer crossover
        return align1, align2
    
    #* Se não recebemos o ponto de crossover fazer random (só por segurança pq na main passamos sempre random tbm)
    if crossover_point is None:
        crossover_point = random.randint(1, len(offsets1) - 1)
    
    #* Criar novos offsets trocando no ponto de cross
    new_offsets1 = offsets1[:crossover_point] + offsets2[crossover_point:]
    new_offsets2 = offsets2[:crossover_point] + offsets1[crossover_point:]
    
    #* Realinhar
    offspring1 = reconstruct_from_offsets(new_offsets1, residues1)
    offspring2 = reconstruct_from_offsets(new_offsets2, residues2)
    
    return offspring1, offspring2


