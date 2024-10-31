# Motif Search

## The Brute-Force Motif Search
### Task 1
* In R, create a function `Score()`, that calculates the score for a consensus string.

* Input:
    * an array of starting indexes
    * `DNAStringSet` of sequences (for example file `seq_score.fasta`)
    * motif length
   
* Output:
    * the score for the consensus string
 
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
    BiocManager::install("Biostrings")
 
library(Biostrings)
 
Score <- function(indexes_array, sequences_file, motif_len) {
  #Načtení fasta souboru
  dna_set <- readDNAStringSet(sequences_file, format = "fasta")
  
  #Získání motivů
  motif_list <- c()
  for (i in 1:length(dna_set)) {
    seq = dna_set[i]
    start <- indexes_array[i]
    motif = seq[start:(start + motif_len - 1)]
    motif_list <- c(motif_list, motif)
  }
  motif_set <- DNAStringSet(motif_list)
  
  
}

sekvence <- DNAStringSet(c("ATCCGTA", "GTGCATA", "AAGCGTG", "ATGCGTG"))

# Vytvoření frekvenčního profilu
frekvencni_profil <- consensusMatrix(sekvence, baseOnly = TRUE)
print(frekvencni_profil)

 
### Task 2
* In R, create function `NextLeaf()` according to the following pseudocode.

* Input:
    * *l*-mer `a = (a1 a2 … al)`, where *l* is motif length; an array of starting indexes
    * `L`; number of DNA sequences
    * `k = n - l + 1`, where `n` is length of DNA sequences and `l` is motif length

* Output:
    * an array of starting indexes that corresponds to the next leaf in the tree

```
NextLeaf(a, L, k)
1   for i ← L to 1
2     if a[i] < k
3       a[i] ← a[i] + 1
4       return a
5     a[i] ← 1
6   return a
```

### Task 3
* In R, create a function `BFMotifSearch()` according to the following pseudocode.

* Input:
    * `DNA`; `DNAStringSet` of DNA sequences (for example file `seq_motif.fasta`)
    * `t`; number of DNA sequences
    * `n`; length of each DNA sequence
    * `l`; motif length

* Output:
    * an array of starting positions for each DNA sequence with the best score for the consensus string

```
BFMotifSearch(DNA, t, n, l)
1   s ← (1, 1, ... , 1)
2   bestScore ← Score(s, DNA, l)
3   while forever
4     s ← NextLeaf(s, t, n − l + 1)
5       if Score(s, DNA, l) > bestScore
6         bestScore ← Score(s, DNA, l)
7         bestMotif ← (s1, s2, . . . , st)
8       if s = (1, 1, . . . , 1)
9         return bestMotif
```

## The Branch-and-Bound Motif Search
### Task 4
* In R, create a function `NextVertex()` according to the following pseudocode.

* Input:
    * *l*-mer `a = (a1 a2 … al)`, where *l* is motif length; an array of starting indexes
    * `i`; level of vertex
    * `L`; number of DNA sequences
    * `k = n - l + 1`, where `n` is length of DNA sequences and `l` is motif length

* Output:
    * *l*-mer of the next vertex in the tree
    * current level of vertex

```
NextVertex(a, i, L, k)
1   if i < L
2     a[i + 1] ← 1
3     return (a, i + 1)
4   else
5     for j ← L to 1
6       if a[j] < k
7         a[j] ← a[j] + 1
8         return (a, j)
9   return (a, 0)
```

### Task 5
* In R, create a function `ByPass()` according to the following pseudocode.

* Input:
    * *l*-mer `a = (a1 a2 … al)`, where *l* is motif length; an array of starting indexes
    * `i`; level of vertex
    * `L`; number of DNA sequences
    * `k = n - l + 1`, where `n` is length of DNA sequences and `l` is motif length

* Output:
    * *l*-mer of the next leaf after a skip of the subtree
    * current level of vertex

```
ByPass(a, i, L, k)
1   for j ← i to 1
2     if a[j] < k
3       a[j] ← a[j] + 1
4       return (a, j)
5   return (a, 0)
```

### Task 6
* In R, create a function `BBMotifSearch()` according to the following pseudocode.

* Input:
    * `DNA`; `DNAStringSet` of DNA sequences (for example file `seq_motif.fasta`)
    * `t`; number of DNA sequences
    * `n`; length of each DNA sequence
    * `l`; motif length

* Output:
    * an array of starting positions for each DNA sequence with the best score for the consensus string

* Modify function `Score()` to calculate a partial consensus score for first `i` rows of `DNA`.

```
BBMotifSearch(DNA, t, n, l)
1   s ← (1, ... , 1)
2   bestScore ← 0
3   i ← 1
4   while i > 0
5     if i < t
6       optimisticScore ← Score(s, i, DNA, l) + (t - i) * l
7       if optimisticScore < bestScore
8         (s, i) ← ByPass(s, i, t, n - l + 1)
9       else
10        (s, i) ← NextVertex(s, i, t, n − l + 1)
11    else
12      if Score(s, t, DNA, l) > bestScore
13        bestScore ← Score(s, t, DNA, l)
14        bestMotif ← (s1, s2, ... , st)
15      (s, i) ← NextVertex(s, i, t, n − l + 1)
16  return bestMotif
```


<details>
<summary>Download files from GitHub</summary>
<details>
<summary>Basic Git settings</summary>

>* Configure the Git editor
>    ```bash
>    git config --global core.editor notepad
>    ```
>* Configure your name and email address
>    ```bash
>    git config --global user.name "Zuzana Nova"
>    git config --global user.email z.nova@vut.cz
>    ```
>* Check current settings
>    ```bash
>    git config --global --list
>    ```
>
</details>

* Create a fork on your GitHub account. 
  On the GitHub page of this repository find a <kbd>Fork</kbd> button in the upper right corner.
  
* Clone forked repository from your GitHub page to your computer:
```bash
git clone <fork repository address>
```
* In a local repository, set new remote for a project repository:
```bash
git remote add upstream https://github.com/mpa-prg/exercise_07.git
```

#### Send files to GitHub
Create a new commit and send new changes to your remote repository.
* Add file to a new commit.
```bash
git add <file_name>
```
* Create a new commit, enter commit message, save the file and close it.
```bash
git commit
```
* Send a new commit to your GitHub repository.
```bash
git push origin main
```

</details>
