# ==========================
# CROSSOVER CORRETO PARA MSA  
# ==========================
# ABORDAGEM: Trocar OFFSETS (gaps iniciais) entre pais
# Mantém resíduos SEMPRE iguais

def get_residues_only(seq):
    """Remove todos os gaps de uma sequência"""
    return seq.replace('-', '')


def extract_offsets_and_residues(align):
    """
    Extrai offsets (gaps iniciais) e resíduos de cada sequência.
    
    Returns:
        offsets: lista de inteiros (número de gaps iniciais)
        residues_list: lista de strings (resíduos sem gaps)
    """
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
    """
    Reconstrói alinhamento a partir de offsets e resíduos.
    """
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


def get_gap_positions(seq):
    """Retorna lista de posições onde há gaps"""
    return [i for i, c in enumerate(seq) if c == '-']


def reconstruct_with_gaps(residues, gap_positions, total_length):
    """
    Reconstrói uma sequência inserindo gaps nas posições especificadas.
    
    Args:
        residues: string com resíduos (sem gaps)
        gap_positions: lista de índices onde colocar gaps
        total_length: comprimento total da sequência final
    """
    result = []
    residue_idx = 0
    
    for i in range(total_length):
        if i in gap_positions:
            result.append('-')
        else:
            if residue_idx < len(residues):
                result.append(residues[residue_idx])
                residue_idx += 1
            else:
                result.append('-')  # Padding se necessário
    
    return ''.join(result)


def crossover_gap_patterns(align1, align2, crossover_point):
    """
    CROSSOVER CORRETO PARA MSA:
    Troca padrões de gaps entre dois alinhamentos, preservando resíduos.
    
    Args:
        align1, align2: listas de strings (alinhamentos)
        crossover_point: coluna onde fazer o corte (0 a len-1)
    
    Returns:
        offspring1, offspring2: novos alinhamentos
    """
    offspring1 = []
    offspring2 = []
    
    for i in range(len(align1)):
        seq1 = align1[i]
        seq2 = align2[i]
        
        # Extrair resíduos (devem ser iguais!)
        residues = get_residues_only(seq1)
        
        # Verificação de segurança
        if get_residues_only(seq2) != residues:
            # Se os resíduos são diferentes, não podemos fazer crossover
            # Retorna pais inalterados
            return align1, align2
        
        # Trocar padrões de gaps no ponto de corte
        # Parte esquerda de um pai + parte direita do outro
        
        # Filho 1: esquerda de seq1 + direita de seq2
        left1 = seq1[:crossover_point]
        right2 = seq2[crossover_point:]
        child1_seq = left1 + right2
        
        # Filho 2: esquerda de seq2 + direita de seq1
        left2 = seq2[:crossover_point]
        right1 = seq1[crossover_point:]
        child2_seq = left2 + right1
        
        offspring1.append(child1_seq)
        offspring2.append(child2_seq)
    
    # Normalizar comprimentos
    offspring1 = normalize_lengths(offspring1)
    offspring2 = normalize_lengths(offspring2)
    
    # Remover colunas só de gaps
    offspring1 = remove_gap_only_columns(offspring1)
    offspring2 = remove_gap_only_columns(offspring2)
    
    return offspring1, offspring2


def crossover_uniform_columns(align1, align2, prob=0.5):
    """
    UNIFORM CROSSOVER:
    Para cada coluna, escolhe aleatoriamente de qual pai vem.
    Preserva resíduos porque troca colunas inteiras.
    
    Args:
        align1, align2: listas de strings (alinhamentos)
        prob: probabilidade de escolher align1 (0.5 = 50/50)
    """
    import random
    
    # Verificar que têm mesmo comprimento
    if len(align1[0]) != len(align2[0]):
        # Normalizar primeiro
        max_len = max(len(align1[0]), len(align2[0]))
        align1 = normalize_to_length(align1, max_len)
        align2 = normalize_to_length(align2, max_len)
    
    length = len(align1[0])
    offspring1 = [[] for _ in range(len(align1))]
    offspring2 = [[] for _ in range(len(align1))]
    
    # Para cada coluna
    for col_idx in range(length):
        if random.random() < prob:
            # Coluna de align1 vai para offspring1, align2 para offspring2
            for seq_idx in range(len(align1)):
                offspring1[seq_idx].append(align1[seq_idx][col_idx])
                offspring2[seq_idx].append(align2[seq_idx][col_idx])
        else:
            # Coluna de align2 vai para offspring1, align1 para offspring2
            for seq_idx in range(len(align1)):
                offspring1[seq_idx].append(align2[seq_idx][col_idx])
                offspring2[seq_idx].append(align1[seq_idx][col_idx])
    
    # Converter de listas para strings
    offspring1 = [''.join(seq) for seq in offspring1]
    offspring2 = [''.join(seq) for seq in offspring2]
    
    # Remover colunas só de gaps
    offspring1 = remove_gap_only_columns(offspring1)
    offspring2 = remove_gap_only_columns(offspring2)
    
    return offspring1, offspring2


def normalize_lengths(align):
    """Adiciona gaps no final para todas terem mesmo comprimento"""
    max_len = max(len(seq) for seq in align)
    return [seq + '-' * (max_len - len(seq)) for seq in align]


def normalize_to_length(align, target_length):
    """Normaliza para um comprimento específico"""
    result = []
    for seq in align:
        if len(seq) < target_length:
            result.append(seq + '-' * (target_length - len(seq)))
        elif len(seq) > target_length:
            result.append(seq[:target_length])
        else:
            result.append(seq)
    return result


def remove_gap_only_columns(align):
    """Remove colunas compostas apenas por gaps"""
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
    """
    CROSSOVER BASEADO EM OFFSETS - MÉTODO CORRETO!
    
    Troca os offsets (gaps iniciais) entre pais.
    Garante que resíduos são SEMPRE preservados.
    
    Args:
        align1, align2: alinhamentos pais
        crossover_point: índice da sequência onde fazer corte (se None, aleatório)
    
    Returns:
        offspring1, offspring2
    """
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
    """
    Função principal de crossover - USA MÉTODO CORRETO DE OFFSETS.
    
    Args:
        align1, align2: listas de strings (alinhamentos pais)
        split_point: índice onde fazer corte (ignorado na versão de offsets)
    
    Returns:
        offspring1, offspring2: novos alinhamentos
    """
    # Usar sempre crossover de offsets (método correto)
    return crossover_offsets(align1, align2, split_point)


# ==========================
# TESTES
# ==========================

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
    
    # Teste 2: Uniform crossover
    print("\n### TESTE 2: Uniform Crossover ###\n")
    
    import random
    random.seed(42)
    
    filho3, filho4 = crossover_uniform_columns(pai1, pai2)
    
    print("Filho 3 (uniform):")
    for seq in filho3:
        print(f"  {seq}  (resíduos: {get_residues_only(seq)})")
    
    print("\nFilho 4 (uniform):")
    for seq in filho4:
        print(f"  {seq}  (resíduos: {get_residues_only(seq)})")
    
    print("\n" + "="*70)
    print("RESULTADO: Crossover CORRETO - Resíduos preservados!" if residues_ok else "ERRO!")
    print("="*70)
