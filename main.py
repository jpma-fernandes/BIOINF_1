from multiplealign.myseq import MySeq
#from multiplealign.multiplealign import MultipleAlignment
from multiplealign.myalign import MyAlign
from multiplealign.substmatrix import SubstMatrix
from multiplealign.pairwisealignment import PairwiseAlignment
import random
import numpy as np
import crossover as cross

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


def create_random_alignment(sequences, max_offset=50):
    aligned_seqs = []
    
    #* colocar offsets de comprimento random no inicio das sequencias
    for seq in sequences:
        offset = random.randint(0, max_offset)
        gaps = '-' * offset
        aligned_seq = gaps + seq
        aligned_seqs.append(aligned_seq)
    
    #* comprimento máximo
    max_length = max(len(seq) for seq in aligned_seqs)
    
    #* Preenche o final com gaps para todas terem o mesmo comprimento
    for i in range(len(aligned_seqs)):
        gaps_needed = max_length - len(aligned_seqs[i])
        aligned_seqs[i] = aligned_seqs[i] + ('-' * gaps_needed)
    
    #* Remove colunas que so tem gaps
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
    score = 0
    n = len(msa)
    
    #* Calcular score de todos os pares
    for i in range(n):
        for j in range(i + 1, n):
            pa.seq1 = msa[i]
            pa.seq2 = msa[j]
            score += pa.score_align()
    
    return score

def print_indv(indv):
    if isinstance(indv, list):
        for seq in indv:
            print(seq)
    else:
        for seq in indv.seqs:
            print(seq)


def select_parents(scored_population, num_parents):
    # Converter scores negativos para positivos (adicionar offset se necessário)
    scores = [score for _, score in scored_population]
    min_score = min(scores)
    
    # Se há scores negativos, adiciona offset para torná los positivos
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
    #* ELITISMO: Manter os melhores
    num_elite = max(1, int(population_size * elite_size))
    next_gen = current_population[:num_elite]  # Já está ordenado
    
    print(f"Keeping {num_elite} in next generation")
    
    num_mutations = 0
    num_crossovers = 0
    
    while len(next_gen) < population_size:
        #? Nao é para fazer crossover e depois mutação? paper diz que não
        # Escolhe entre mutação e crossover
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
            
            offspring1, offspring2 = cross.generate_offspring(father, mother, split_point)
            
            # Escolher o melhor dos dois
            score1 = score_MSA(offspring1)
            score2 = score_MSA(offspring2)
            offspring = offspring1 if score1 > score2 else offspring2
            num_crossovers += 1
        
        # Calcular score do filho
        offspring_score = score_MSA(offspring)
        next_gen.append((offspring, offspring_score))
    
    print(f"Operacoes: {num_crossovers} crossovers, {num_mutations} mutacoes")
    
    return next_gen


def count_residues(seq):
    return sum(1 for c in seq if c != '-')


def mutate_split_gap_block(alignment):
    import copy
    mutated = copy.deepcopy(alignment)
    
    #! Nossa operacao é  Split a randomly selected gap block
    
    seq_idx = random.randint(0, len(mutated) - 1)
    checked = 0

    while checked < len(mutated):
        seq = mutated[seq_idx]

        #* Encontrar blocos com pelo menos 2 gaps seguidos
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

        #* Se encontrámos blocos, paramos aqui
        if gap_blocks:
            break

        #* Caso contrário, avançamos circularmente para a próxima sequência
        seq_idx = (seq_idx + 1) % len(mutated)
        checked += 1

    #* Se nenhuma sequência tiver blocos de gaps, devolvemos mutated sendo igual ao pai
    if checked == len(mutated):
        return mutated
    
    #* escolher bloco que vamos partir
    block_start, block_end = random.choice(gap_blocks)

    #* dividir o tamanho do bloco em 2 outros blocos menores
    block_length = block_end - block_start
    left_length = random.randint(1, block_length - 1)
    right_length = block_length - left_length

    #* escolher a posição de inserção na esquerda
    left_shift = random.choice([-3, -2, -1])
    insert_left_pos = max(0, block_start + left_shift) #* evitar out of bounds
    
    #* escolher a posição de inserção na direita
    right_shift = random.choice([1, 2, 3])
    insert_right_pos = min(len(seq), block_end + right_shift) #* evitar out of bounds

    new_seq = (
        seq[:insert_left_pos] +              #* Tudo até ao inicio do insert da esquerda
        '-' * left_length +                  #* Inserimos o numero de gaps
        seq[insert_left_pos: block_start] +  #* repomos o que vem depois do gap até ao inicio do bloco deleted
        seq[block_end: insert_right_pos] +   #* Pulamos o bloco deleted e colocamos tudo desde o fim do deleted até à zona de inserção da direita
        '-' * right_length +                 #* inserimos o numero de gaps
        seq[insert_right_pos:]               #* repomos tudo o que vem depois do gap até ao fim
    )
    #print(f"\nSequence before mutation: ${seq}")
    #print(f"New seq: ${new_seq}\n")
    
    mutated[seq_idx] = ''.join(new_seq)
    
    # Remover colunas só de gaps
    mutated = remove_gap_only_columns(mutated)
    
    return mutated


def run_genetic_algorithm(sequences, population_size=10, max_generations=100, 
                         no_improvement_limit=20, max_offset=10):

    print("ALGORITMO GENÉTICO - MULTIPLE SEQUENCE ALIGNMENT")
    print(f"População: {population_size} | Max Gerações: {max_generations} | Limite sem melhoria: {no_improvement_limit}")
    print()

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
    
    print(f"Geracao 0: Melhor score = {best_ever[1]}")
    print()
    
    #* Loop gerações
    for generation in range(1, max_generations + 1):
        print(f"--- Geração {generation} ---")
        
        scored_population = create_next_generation(scored_population, population_size, 0.5)
        
        #* Sort por score
        scored_population.sort(key=lambda x: x[1], reverse=True)
        
        best_current = scored_population[0]
        best_score_history.append(best_current[1])
        
        #* Verificar se houve melhoria overall
        worst_current = scored_population[-1]
        avg_score = sum(s for _, s in scored_population) / len(scored_population)
        
        if best_current[1] > best_ever[1]:
            improvement = best_current[1] - best_ever[1]
            print(f"Novo best Score: {best_current[1]} (+{improvement}) | Média: {avg_score:.2f}")
            best_ever = best_current
            generations_without_improvement = 0
        else:
            generations_without_improvement += 1
            print(f"Melhor: {best_current[1]} | Pior: {worst_current[1]} | Média: {avg_score:.2f} | Sem melhoria: {generations_without_improvement}")
        
        #* Stop criteria
        if generations_without_improvement >= no_improvement_limit:
            print()
            print(f"Stopping Sem melhoria por {no_improvement_limit} gerações consecutivas")
            break
    
    print()
    print("RESULTADO FINAL")
    print(f"Melhor score encontrado: {best_ever[1]}")
    print(f"Gerações executadas: {generation}")
    print()
    print("Melhor alinhamento:")
    print_indv(best_ever[0])
    
    return best_ever[0], best_ever[1], best_score_history


if __name__ == "__main__":

    sequences = read_fasta('./cytochromes.fa')
    submat = SubstMatrix()
    submat.read_submat_file("blosum62.mat") 
    pa = PairwiseAlignment(submat, -8) 
    
    print("Sequências originais:")
    for i, seq in enumerate(sequences):
        print(f"Seq {i}: {seq}")
    print()
    
    #* Executar Algoritmo Genético
    best_alignment, best_score, history = run_genetic_algorithm(
        sequences,
        population_size=50, 
        max_generations=100,  
        no_improvement_limit=30, 
        max_offset=10
    )
    
    # 3. Mostrar evolução
    print()
    print("Evolução do melhor score:")
    for gen, score in enumerate(history[::5]):  # A cada 5 gerações
        print(f"Gen {gen*5:3d}: {score:.2f}")
    print(f"Gen {len(history)-1:3d}: {history[-1]:.2f} (final)")




    
    