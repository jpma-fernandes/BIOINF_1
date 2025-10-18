from multiplealign.myseq import MySeq
#from multiplealign.multiplealign import MultipleAlignment
from multiplealign.substmatrix import SubstMatrix
from multiplealign.pairwisealignment import PairwiseAlignment
import random
import numpy as np
import crossover_corrected as crossover

PROTEIN_TYPE = "PROTEIN"
pa = None

def read_fasta(filename):
    sequences = []
    previous_line = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            #* Quando encontramos uma sequencia nova guardamos a antiga
            if line.startswith('<') or line.startswith('&'):

                #* Guardar sequência caso exista sequencia (pula só no 1 caso)
                if previous_line is not None:
                    seq_str = ''.join(current_seq)
                    sequences.append(seq_str)
                
                #* Iniciar nova sequência
                previous_line = line
                current_seq = []

            #* linhas com a proteina guardamos a sequencia
            else:
                current_seq.append(line) 
    return sequences


def initialize_population(sequences, population_size=100, max_offset=50):
    population = []
    
    for i in range(population_size):
        alignment = create_random_alignment(sequences, max_offset)
        population.append(alignment)
    
    return population


def create_random_alignment(sequences, max_offset=50, insert_internal_gaps=True):
    aligned_seqs = []
    
    for seq in sequences:
        # 1. Offset inicial (gaps no início)
        offset = random.randint(0, max_offset)
        aligned_seq = '-' * offset + seq
        
        # 2. NOVO: Adicionar gaps aleatórios NO MEIO da sequência (mais diversidade!)
        if insert_internal_gaps and random.random() > 0.3:  # 70% de chance
            num_internal_gaps = random.randint(1, max_offset // 5)  # Poucos gaps internos
            for _ in range(num_internal_gaps):
                # Inserir gap numa posição aleatória (não no início nem fim)
                if len(aligned_seq) > 2:
                    pos = random.randint(1, len(aligned_seq) - 1)
                    aligned_seq = aligned_seq[:pos] + '-' + aligned_seq[pos:]
        
        aligned_seqs.append(aligned_seq)
    
    # Comprimento máximo
    max_length = max(len(seq) for seq in aligned_seqs)
    
    # Preenche o final com gaps para todas terem o mesmo comprimento
    for i in range(len(aligned_seqs)):
        gaps_needed = max_length - len(aligned_seqs[i])
        aligned_seqs[i] = aligned_seqs[i] + ('-' * gaps_needed)
    
    # Remove colunas que só têm gaps
    aligned_seqs = remove_gap_only_columns(aligned_seqs)
    
    return aligned_seqs


def remove_gap_only_columns(aligned_seqs):
    if not aligned_seqs:
        return aligned_seqs
    
    #* Neste ponto todas as sequencias tem o mesmo tamanho
    seq_length = len(aligned_seqs[0])
    
    #* Identifica colunas validas (que tem pelo menos um nao-gap)
    valid_columns = []
    for col_idx in range(seq_length):
        has_non_gap = False
        for seq in aligned_seqs:
            if seq[col_idx] != '-':
                has_non_gap = True
                break
        if has_non_gap:
            valid_columns.append(col_idx)
    
    #* Reconstroi as sequencias apenas com colunas validas
    result = []
    for seq in aligned_seqs:
        new_seq = ''.join(seq[col_idx] for col_idx in valid_columns)
        result.append(new_seq)
    
    return result


def score_MSA(msa):
    """
    Calcula Sum-of-Pairs score do MSA.
    CORRIGIDO: Usa o número correto de sequências e acessa BLOSUM62 corretamente.
    """
    score = 0
    n = len(msa)  # ✅ CORRIGIDO: sem o -1!
    
    #* Calcular score de todos os pares
    for i in range(n):
        for j in range(i + 1, n):  # j vai de i+1 até n-1 (correto agora!)
            pa.seq1 = msa[i]
            pa.seq2 = msa[j]
            score += pa.score_align()
    
    return score

def print_indv(indv):
    """Imprime um individuo (alinhamento)"""
    if isinstance(indv, list):
        for seq in indv:
            print(seq[:80] + "..." if len(seq) > 80 else seq)
    else:
        for seq in indv.seqs:
            print(seq[:80] + "..." if len(seq) > 80 else seq)


def select_parents(scored_population, num_parents):
    """
    Seleciona pais usando Roulette Wheel Selection (seleção proporcional ao fitness).
    
    Args:
        scored_population: lista de tuplas (individuo, score)
        num_parents: quantos pais selecionar
    
    Returns:
        lista de individuos selecionados
    """
    # Converter scores negativos para positivos (adicionar offset se necessário)
    scores = [score for _, score in scored_population]
    min_score = min(scores)
    
    # Se todos os scores são negativos, adiciona offset para torná-los positivos
    if min_score < 0:
        adjusted_scores = [score - min_score + 1 for score in scores]
    else:
        adjusted_scores = scores
    
    total_fitness = sum(adjusted_scores)
    
    # Probabilidades proporcionais ao fitness
    probabilities = [fitness / total_fitness for fitness in adjusted_scores]
    
    # Selecionar pais
    indices = np.random.choice(len(scored_population), size=num_parents, p=probabilities, replace=True)
    selected_parents = [scored_population[i][0] for i in indices]
    
    return selected_parents


def create_next_generation(current_population, population_size, elite_size=0.1, mutation_prob=0.5):
    """
    Cria a próxima geração através de elitismo + crossover + mutação.
    
    Args:
        current_population: lista de tuplas (individuo, score)
        population_size: tamanho da população
        elite_size: percentagem de elite a manter (0.2 = 20%)
        mutation_prob: probabilidade de mutação vs crossover (0.3 = 30%)
    
    Returns:
        nova população (lista de tuplas (individuo, score))
    """
    # 1. ELITISMO: Manter os melhores
    num_elite = max(1, int(population_size * elite_size))
    next_gen = current_population[:num_elite]  # Já está ordenado
    
    print(f"  Mantendo {num_elite} elite(s)")
    
    # 2. REPRODUÇÃO: Gerar filhos até completar população
    num_mutations = 0
    num_crossovers = 0
    
    while len(next_gen) < population_size:
        # Escolher entre mutação e crossover
        if random.random() < mutation_prob:
            # MUTAÇÃO
            parent = select_parents(current_population, 1)[0]
            offspring = mutate_split_gap_block(parent)
            num_mutations += 1
        else:
            # CROSSOVER
            parents = select_parents(current_population, 2)
            father, mother = parents[0], parents[1]
            
            max_residues = min(count_residues(father[0]), count_residues(mother[0]))
            if max_residues > 1:
                split_point = random.randint(1, max_residues - 1)
            else:
                split_point = 1
            
            offspring1, offspring2 = crossover.generate_offspring(father, mother, split_point)
            
            # Escolher o melhor dos dois
            score1 = score_MSA(offspring1)
            score2 = score_MSA(offspring2)
            offspring = offspring1 if score1 > score2 else offspring2
            num_crossovers += 1
        
        # Calcular score do filho
        offspring_score = score_MSA(offspring)
        next_gen.append((offspring, offspring_score))
    
    print(f"  Operações: {num_crossovers} crossovers, {num_mutations} mutações")
    
    return next_gen


def count_residues(seq):
    """Conta resíduos (não-gaps) numa sequência"""
    return sum(1 for c in seq if c != '-')


def mutate_split_gap_block(alignment):
    """
    MUTAÇÃO OBRIGATÓRIA: Split a randomly selected gap block.
    
    Escolhe uma sequência aleatória, encontra um bloco de gaps,
    e divide-o inserindo uma coluna nova.
    """
    import copy
    mutated = copy.deepcopy(alignment)
    
    # Escolher tipo de mutação aleatoriamente
    mutation_type = random.choice(['split_gap', 'remove_gap', 'insert_gap'])
    
    if mutation_type == 'split_gap':
        # MUTAÇÃO ORIGINAL: Split gap block
        seq_idx = random.randint(0, len(mutated) - 1)
        seq = mutated[seq_idx]
        
        gap_blocks = []
        i = 0
        while i < len(seq):
            if seq[i] == '-':
                start = i
                while i < len(seq) and seq[i] == '-':
                    i += 1
                end = i
                if end - start >= 2:
                    gap_blocks.append((start, end))
            else:
                i += 1
        
        if gap_blocks:
            block_start, block_end = random.choice(gap_blocks)
            split_pos = random.randint(block_start + 1, block_end - 1)
            mutated = [s[:split_pos] + '-' + s[split_pos:] for s in mutated]
        else:
            # Inserir gap aleatório
            pos = random.randint(0, len(seq))
            mutated = [s[:pos] + '-' + s[pos:] for s in mutated]
    
    elif mutation_type == 'remove_gap':
        # NOVA MUTAÇÃO: Remover gap aleatório de uma sequência
        seq_idx = random.randint(0, len(mutated) - 1)
        seq = mutated[seq_idx]
        
        gap_positions = [i for i, c in enumerate(seq) if c == '-']
        if gap_positions:
            remove_pos = random.choice(gap_positions)
            # Remover gap dessa coluna em todas as sequências
            mutated = [s[:remove_pos] + s[remove_pos+1:] if i == seq_idx else s 
                      for i, s in enumerate(mutated)]
            # Ajustar comprimentos
            max_len = max(len(s) for s in mutated)
            mutated = [s + '-' * (max_len - len(s)) for s in mutated]
    
    elif mutation_type == 'insert_gap':
        # NOVA MUTAÇÃO: Inserir gap em posição aleatória
        seq_idx = random.randint(0, len(mutated) - 1)
        pos = random.randint(0, len(mutated[0]))
        mutated[seq_idx] = mutated[seq_idx][:pos] + '-' + mutated[seq_idx][pos:]
        # Ajustar outras sequências
        for i in range(len(mutated)):
            if i != seq_idx and len(mutated[i]) < len(mutated[seq_idx]):
                mutated[i] += '-'
    
    # Remover colunas só de gaps
    mutated = remove_gap_only_columns(mutated)
    
    return mutated


def apply_mutation_or_crossover(current_population, mutation_prob=0.3):
    """
    Escolhe aleatoriamente entre mutação e crossover.
    
    Args:
        current_population: população atual (lista de tuplas)
        mutation_prob: probabilidade de mutação (0.3 = 30%)
    
    Returns:
        novo indivíduo (offspring)
    """
    if random.random() < mutation_prob:
        # MUTAÇÃO
        parent = select_parents(current_population, 1)[0]
        offspring = mutate_split_gap_block(parent)
        return offspring
    else:
        # CROSSOVER
        parents = select_parents(current_population, 2)
        father, mother = parents[0], parents[1]
        
        max_residues = min(count_residues(father[0]), count_residues(mother[0]))
        if max_residues > 1:
            split_point = random.randint(1, max_residues - 1)
        else:
            split_point = 1
        
        offspring1, offspring2 = crossover.generate_offspring(father, mother, split_point)
        
        # Retornar o melhor dos dois
        score1 = score_MSA(offspring1)
        score2 = score_MSA(offspring2)
        
        return offspring1 if score1 > score2 else offspring2


def run_genetic_algorithm(sequences, population_size=10, max_generations=100, 
                         no_improvement_limit=20, max_offset=10):
    """
    Executa o Algoritmo Genético completo.
    
    Args:
        sequences: sequências originais (sem gaps)
        population_size: tamanho da população
        max_generations: número máximo de gerações
        no_improvement_limit: parar se não melhorar por N gerações
        max_offset: máximo de gaps iniciais
    
    Returns:
        melhor alinhamento encontrado, seu score, histórico
    """
    print("=" * 80)
    print("ALGORITMO GENÉTICO - MULTIPLE SEQUENCE ALIGNMENT")
    print("=" * 80)
    print(f"População: {population_size} | Max Gerações: {max_generations} | Limite sem melhoria: {no_improvement_limit}")
    print()
    
    # 1. INICIALIZAR POPULAÇÃO (Geração 0)
    print("Inicializando população...")
    population = initialize_population(sequences, population_size, max_offset)
    
    # Avaliar e ordenar
    scored_population = sorted(
        [(ind, score_MSA(ind)) for ind in population],
        key=lambda x: x[1],
        reverse=True
    )
    
    best_ever = scored_population[0]
    best_score_history = [best_ever[1]]
    generations_without_improvement = 0
    
    print(f"Geração 0: Melhor score = {best_ever[1]:.2f}")
    print()
    
    # 2. LOOP DE GERAÇÕES
    for generation in range(1, max_generations + 1):
        print(f"--- Geração {generation} ---")
        
        # Criar próxima geração
        scored_population = create_next_generation(scored_population, population_size)
        
        # Ordenar por fitness
        scored_population.sort(key=lambda x: x[1], reverse=True)
        
        # Melhor da geração
        best_current = scored_population[0]
        best_score_history.append(best_current[1])
        
        # Verificar se houve melhoria
        worst_current = scored_population[-1]
        avg_score = sum(s for _, s in scored_population) / len(scored_population)
        
        if best_current[1] > best_ever[1]:
            improvement = best_current[1] - best_ever[1]
            print(f"  ✓ NOVO MELHOR! Score: {best_current[1]:.2f} (+{improvement:.2f}) | Média: {avg_score:.2f}")
            best_ever = best_current
            generations_without_improvement = 0
        else:
            generations_without_improvement += 1
            print(f"  Melhor: {best_current[1]:.2f} | Pior: {worst_current[1]:.2f} | Média: {avg_score:.2f} | Sem melhoria: {generations_without_improvement}")
        
        # Critério de paragem
        if generations_without_improvement >= no_improvement_limit:
            print()
            print(f">>> PARAGEM: Sem melhoria por {no_improvement_limit} gerações consecutivas")
            break
    
    print()
    print("=" * 80)
    print("RESULTADO FINAL")
    print("=" * 80)
    print(f"Melhor score encontrado: {best_ever[1]:.2f}")
    print(f"Gerações executadas: {generation}")
    print()
    print("Melhor alinhamento:")
    print_indv(best_ever[0])
    
    return best_ever[0], best_ever[1], best_score_history


if __name__ == "__main__":
    # 1. Carregar dados
    sequences = read_fasta('./cytochromes.fa')
    submat = SubstMatrix()
    submat.read_submat_file("blosum62.mat")
    pa = PairwiseAlignment(submat, -1)  # Gap penalty = -1 (menos penalidade = scores mais altos)
    
    print("Sequências originais:")
    for i, seq in enumerate(sequences[:3]):  # Mostrar só as 3 primeiras
        print(f"Seq {i}: {seq[:60]}..." if len(seq) > 60 else f"Seq {i}: {seq}")
    print(f"... (total de {len(sequences)} sequências)")
    print()
    
    # 2. Executar Algoritmo Genético
    best_alignment, best_score, history = run_genetic_algorithm(
        sequences,
        population_size=50,       # População maior = mais diversidade
        max_generations=100,      # Mais gerações
        no_improvement_limit=30,  # Esperar mais antes de parar
        max_offset=30
    )
    
    # 3. Mostrar evolução
    print()
    print("Evolução do melhor score:")
    for gen, score in enumerate(history[::5]):  # A cada 5 gerações
        print(f"Gen {gen*5:3d}: {score:.2f}")
    print(f"Gen {len(history)-1:3d}: {history[-1]:.2f} (final)")




    
    