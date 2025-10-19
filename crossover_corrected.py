# ==========================
# CROSSOVER CORRETO PARA MSA  
# ==========================
# ABORDAGEM: Trocar OFFSETS (gaps iniciais) entre pais
# Mantém resíduos SEMPRE iguais

def get_residues_only(seq):
    return seq.replace('-', '')


def extract_offsets_and_residues(align):

    offsets = []
    residues_list = []
    
    for seq in align:
        # Contar gaps iniciais
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

    # Criar sequências com offsets
    seqs = []
    for offset, residues in zip(offsets, residues_list):
        seq = '-' * offset + residues
        seqs.append(seq)
    
    # Normalizar comprimentos
    max_len = max(len(s) for s in seqs)
    seqs = [s + '-' * (max_len - len(s)) for s in seqs]
    
    # Remover colunas só de gaps
    return remove_gap_only_columns(seqs)


def remove_gap_only_columns(align):
    if not align or len(align[0]) == 0:
        return align
    
    cols_to_keep = []
    for col_idx in range(len(align[0])):
        # Manter coluna se pelo menos uma sequência tem resíduo
        if any(seq[col_idx] != '-' for seq in align):
            cols_to_keep.append(col_idx)
    
    # Reconstruir alinhamento
    result = []
    for seq in align:
        new_seq = ''.join(seq[i] for i in cols_to_keep)
        result.append(new_seq)
    
    return result


def crossover_offsets(align1, align2, crossover_point=None):

    import random
    
    # Extrair offsets e resíduos
    offsets1, residues1 = extract_offsets_and_residues(align1)
    offsets2, residues2 = extract_offsets_and_residues(align2)
    
    # Verificar se resíduos são iguais
    if residues1 != residues2:
        # Se não são, não podemos fazer crossover
        return align1, align2
    
    # Escolher ponto de crossover (qual sequência)
    if crossover_point is None:
        crossover_point = random.randint(1, len(offsets1) - 1)
    
    # Criar novos offsets trocando no ponto de corte
    new_offsets1 = offsets1[:crossover_point] + offsets2[crossover_point:]
    new_offsets2 = offsets2[:crossover_point] + offsets1[crossover_point:]
    
    # Reconstruir alinhamentos
    offspring1 = reconstruct_from_offsets(new_offsets1, residues1)
    offspring2 = reconstruct_from_offsets(new_offsets2, residues2)
    
    return offspring1, offspring2


def generate_offspring(align1, align2, split_point=None):

    # Usar sempre crossover de offsets (método correto)
    return crossover_offsets(align1, align2, split_point)


if __name__ == "__main__":
    print("="*70)
    print("TESTE DO CROSSOVER CORRIGIDO PARA MSA")
    print("="*70)
    
    # Teste 1: One-point crossover
    print("\n### TESTE 1: One-Point Crossover ###\n")
    
    pai1 = ['A--TGC', 'ATG--C']
    pai2 = ['ATGC--', 'A-TG-C']
    
    print("Pai 1:")
    for seq in pai1:
        print(f"  {seq}  (resíduos: {get_residues_only(seq)})")
    
    print("\nPai 2:")
    for seq in pai2:
        print(f"  {seq}  (resíduos: {get_residues_only(seq)})")
    
    filho1, filho2 = crossover_offsets(pai1, pai2, 1)
    
    print("\nFilho 1 (corte na sequência 1 - troca offsets):")
    for seq in filho1:
        print(f"  {seq}  (resíduos: {get_residues_only(seq)})")
    
    print("\nFilho 2:")
    for seq in filho2:
        print(f"  {seq}  (resíduos: {get_residues_only(seq)})")
    
    # Verificação
    print("\n" + "-"*70)
    residues_ok = True
    for i in range(len(pai1)):
        r_pai1 = get_residues_only(pai1[i])
        r_filho1 = get_residues_only(filho1[i])
        r_filho2 = get_residues_only(filho2[i])
        
        if r_pai1 != r_filho1 or r_pai1 != r_filho2:
            residues_ok = False
            print(f"✗ ERRO Seq {i}: Resíduos alterados!")
        else:
            print(f"✓ Seq {i}: Resíduos preservados ({r_pai1})")
    
    print("-"*70)
    
    print("\n" + "="*70)
    print("RESULTADO: Crossover CORRETO - Resíduos preservados!" if residues_ok else "ERRO!")
    print("="*70)
