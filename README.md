# GibbsSampler


# Description
This Python script implements the Gibbs Sampling algorithm to find motifs in a set of DNA sequences. It iteratively refines candidate motifs by probabilistically selecting k-mers based on profile matrices with pseudocounts, aiming to minimize motif score and identify conserved patterns.

# Usage
Example
```
import random

def GibbsSampler(Dna, k, t, N):
    def ProfileWithPseudocounts(motifs):
        profile = {'A': [], 'C': [], 'G': [], 'T': []}
        for i in range(len(motifs[0])):
            col = [motif[i] for motif in motifs]
            total = len(col)
            profile['A'].append((col.count('A') + 1) / (total + 4))
            profile['C'].append((col.count('C') + 1) / (total + 4))
            profile['G'].append((col.count('G') + 1) / (total + 4))
            profile['T'].append((col.count('T') + 1) / (total + 4))
        return profile

    def ProfileMostProbableKmer(text, k, profile):
        max_prob = -1
        most_probable = ""
        for i in range(len(text) - k + 1):
            kmer = text[i:i + k]
            prob = 1
            for j in range(k):
                prob *= profile[kmer[j]][j]
            if prob > max_prob:
                max_prob = prob
                most_probable = kmer
        return most_probable

    def Score(motifs):
        consensus = Consensus(motifs)
        score = 0
        for motif in motifs:
            score += sum(1 for i in range(len(motif)) if motif[i] != consensus[i])
        return score

    def Consensus(motifs):
        consensus = ""
        for i in range(len(motifs[0])):
            counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            for motif in motifs:
                counts[motif[i]] += 1
            consensus += max(counts, key=counts.get)
        return consensus

    def RandomMotifs(Dna, k):
        random_motifs = []
        for string in Dna:
            start = random.randint(0, len(string) - k)
            random_motifs.append(string[start:start + k])
        return random_motifs

    BestMotifs = RandomMotifs(Dna, k)
    BestMotifsScore = Score(BestMotifs)

    for _ in range(N):
        i = random.randint(0, t - 1)
        motifs_except_i = BestMotifs[:i] + BestMotifs[i+1:]
        profile = ProfileWithPseudocounts(motifs_except_i)
        motif_i = ProfileMostProbableKmer(Dna[i], k, profile)
        new_motifs = BestMotifs[:i] + [motif_i] + BestMotifs[i+1:]
        if Score(new_motifs) < BestMotifsScore:
            BestMotifs = new_motifs
            BestMotifsScore = Score(new_motifs)

    return BestMotifs

# Example usage:
k = 8
t = 5
N = 100
Dna = [
    "CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA",
    "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
    "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
    "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
    "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
]

BestMotifs = GibbsSampler(Dna, k, t, N)
for motif in BestMotifs:
    print(motif)

```
# Output
GTAACCGG
GTGCCAAG
GTACCGAG
GTTTCAGG
GCTCCACG

The script outputs a list of the best motifs found, one motif per DNA sequence.

# Function Descriptions
* GibbsSampler(Dna, k, t, N): Executes the Gibbs Sampling algorithm on a list of DNA sequences (Dna), searching for motifs of length k across t sequences with N iterations.
* ProfileWithPseudocounts(motifs): Computes a profile matrix with pseudocounts (+1) to avoid zero probabilities.
* ProfileMostProbableKmer(text, k, profile): Finds the most probable k-mer in a given text based on the profile matrix.
* Score(motifs): Calculates the score of a set of motifs, defined as the total number of mismatches from the consensus.
* Consensus(motifs): Determines the consensus string from a set of motifs.
* RandomMotifs(Dna, k): Selects random initial motifs from each DNA sequence.

# Applications
* Identification of conserved motifs in DNA sequences.
* Useful for discovering regulatory elements such as transcription factor binding sites.
* Applicable in bioinformatics pipelines for motif discovery and sequence analysis.
* Educational tool for understanding probabilistic motif finding algorithms.

# License
This project is licensed under the MIT License.

