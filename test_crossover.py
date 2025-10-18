import crossover

# Teste simples
pai1 = ['A--TGC', 'ATG--C']
pai2 = ['ATGC--', 'A-TGCX']

print("="*60)
print("TESTE DO CROSSOVER")
print("="*60)
print("\nPai 1:")
for seq in pai1:
    print(f"  {seq}")

print("\nPai 2:")
for seq in pai2:
    print(f"  {seq}")

print("\nSplit residues = 2 (cortar depois do 2º resíduo)")

filho1, filho2 = crossover.generate_offspring(pai1, pai2, 2)

print("\nFilho 1 (esquerda Pai1 + direita Pai2):")
for seq in filho1:
    print(f"  {seq}")

print("\nFilho 2 (esquerda Pai2 + direita Pai1):")
for seq in filho2:
    print(f"  {seq}")

# Verificar resíduos preservados
print("\n" + "="*60)
print("VERIFICAÇÃO: Resíduos preservados?")
print("="*60)

def get_residues(seq):
    return seq.replace('-', '')

print("\nPai 1 resíduos:")
for i, seq in enumerate(pai1):
    print(f"  Seq {i}: {get_residues(seq)}")

print("\nFilho 1 resíduos:")
for i, seq in enumerate(filho1):
    print(f"  Seq {i}: {get_residues(seq)}")

print("\nFilho 2 resíduos:")
for i, seq in enumerate(filho2):
    print(f"  Seq {i}: {get_residues(seq)}")

# Verificar se resíduos são iguais
print("\n✓ Resíduos preservados?" if all(
    get_residues(pai1[i]) == get_residues(filho1[i]) == get_residues(filho2[i])
    for i in range(len(pai1))
) else "\n✗ ERRO: Resíduos NÃO preservados!")
