# Como Iniciar um Repositório Git com Ficheiros Existentes

Este guia explica como criar um repositório Git e adicionar ficheiros existentes ao mesmo.

## Método 1: Iniciar um Repositório Local

### Passo 1: Navegar até à pasta do projeto
```bash
cd /caminho/para/seu/projeto
```

### Passo 2: Inicializar o repositório Git
```bash
git init
```

### Passo 3: Adicionar todos os ficheiros existentes
```bash
git add .
```

Ou para adicionar ficheiros específicos:
```bash
git add ficheiro1.txt ficheiro2.py
```

### Passo 4: Fazer o primeiro commit
```bash
git commit -m "Commit inicial com ficheiros existentes"
```

### Passo 5: (Opcional) Conectar a um repositório remoto
```bash
git remote add origin https://github.com/seu-usuario/seu-repositorio.git
git branch -M main
git push -u origin main
```

## Método 2: Clonar um Repositório Existente e Adicionar Ficheiros

### Passo 1: Clonar o repositório vazio
```bash
git clone https://github.com/seu-usuario/seu-repositorio.git
cd seu-repositorio
```

### Passo 2: Copiar os ficheiros existentes para a pasta clonada
```bash
cp -r /caminho/dos/ficheiros/existentes/* .
```

### Passo 3: Adicionar e fazer commit dos ficheiros
```bash
git add .
git commit -m "Adicionar ficheiros existentes ao repositório"
git push origin main
```

## Boas Práticas

### 1. Criar um ficheiro .gitignore
Antes de adicionar os ficheiros, crie um `.gitignore` para excluir ficheiros desnecessários:

```bash
# Criar .gitignore
cat > .gitignore << EOF
# Ficheiros temporários
*.tmp
*.log
*~

# Pastas de dependências
node_modules/
venv/
__pycache__/

# Ficheiros do sistema
.DS_Store
Thumbs.db

# Ficheiros de configuração local
.env
config.local.json
EOF
```

### 2. Verificar o que será adicionado
```bash
git status
```

### 3. Adicionar ficheiros seletivamente
Se tiver ficheiros que não deve incluir no repositório:
```bash
git add ficheiro1.txt
git add pasta/ficheiro2.py
```

### 4. Verificar mudanças antes do commit
```bash
git diff --staged
```

## Exemplo Completo

```bash
# 1. Navegar até à pasta com ficheiros existentes
cd /home/usuario/meu-projeto

# 2. Inicializar Git
git init

# 3. Criar .gitignore
echo "*.log" > .gitignore
echo "*.tmp" >> .gitignore

# 4. Adicionar .gitignore primeiro
git add .gitignore
git commit -m "Adicionar .gitignore"

# 5. Adicionar todos os outros ficheiros
git add .
git commit -m "Adicionar ficheiros existentes do projeto"

# 6. Conectar ao GitHub (se necessário)
git remote add origin https://github.com/seu-usuario/seu-repositorio.git
git branch -M main
git push -u origin main
```

## Comandos Úteis

### Ver histórico de commits
```bash
git log --oneline
```

### Ver estado atual do repositório
```bash
git status
```

### Ver diferenças nos ficheiros
```bash
git diff
```

### Remover ficheiros do staging (antes do commit)
```bash
git reset ficheiro.txt
```

### Desfazer último commit (mantendo as alterações)
```bash
git reset --soft HEAD~1
```

## Resolução de Problemas Comuns

### Problema: Adicionei ficheiros por engano
**Solução:**
```bash
# Antes do commit
git reset ficheiro-errado.txt

# Depois do commit (último commit)
git reset --soft HEAD~1
git reset ficheiro-errado.txt
git commit -m "Adicionar ficheiros corretos"
```

### Problema: O repositório remoto já tem conteúdo
**Solução:**
```bash
git pull origin main --allow-unrelated-histories
git push origin main
```

### Problema: Ficheiros grandes a causar problemas
**Solução:** Use Git LFS (Large File Storage)
```bash
git lfs install
git lfs track "*.zip"
git lfs track "*.mp4"
git add .gitattributes
```

## Recursos Adicionais

- [Documentação Oficial do Git](https://git-scm.com/doc)
- [GitHub Guides](https://guides.github.com/)
- [Tutorial Interativo do Git](https://learngitbranching.js.org/)
