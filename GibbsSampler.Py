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

k = 8
t = 5
N = 100
Dna = [
    "CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA","GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG","TAGTACCGAGACCGAAAGAAGTATACAGGCGT","TAGATCAAGTTTCAGGTGCACGTCGGTGAACC","AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
]

BestMotifs = GibbsSampler(Dna, k, t, N)
for motif in BestMotifs:
    print(motif)
