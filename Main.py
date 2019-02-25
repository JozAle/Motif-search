from Greedy import GreedyMotifSearch
from BruteForce import BruteForceMotifSearch
from Random import RandomMotifSearch

DNA = (
    'tagtggtcttttgagtgtagatctgaagggaaagtatttccaccagttcggggtcacccagcagggcagggtgacttaat',
    'cgcgactcggcgctcacagttatcgcacgtttagaccaaaacggagttggatccgaaactggagtttaatcggagtcctt',
    'gttacttgtgagcctggttagacccgaaatataattgttggctgcatagcggagctgacatacgagtaggggaaatgcgt',
    'aacatcaggctttgattaaacaatttaagcacgtaaatccgaattgacctgatgacaatacggaacatgccggctccggg',
    'accaccggataggctgcttattaggtccaaaaggtagtatcgtaataatggctcagccatgtcaatgtgcggcattccac',
    'tagattcgaatcgatcgtgtttctccctctgtgggttaacgaggggtccgaccttgctcgcatgtgccgaacttgtaccc',
    'gaaatggttcggtgcgatatcaggccgttctcttaacttggcggtgcagatccgaacgtctctgactgtgtagatccgta',
    'atgtatactagacattctaacgctcgcttattggcggagaccatttgctccactacaagaggctactgtgtagatccgta',
    'ttcttacacccttctttagatccaaacctgttggcgccatcttcttttcgagtccttgtacctccatttgctctgatgac',
    'ctacctatgtaaaacaacatctactaacgtagtccggtctttcctgatctgccctaacctacaggtcgatccgaaattcg'
    )

k = 77

with open('Greedy.txt', 'w') as plik:
    Score1, S1 = GreedyMotifSearch(DNA, k, len(DNA))
    print(str(k) + '-mery', 'Score wynosi: ' + str(Score1), 'Motyw:', file=plik)
with open('Greedy.txt', 'a') as plik:
    for i, Sekwencja in enumerate(DNA):
        print(Sekwencja[S1[i]:S1[i] + k], file=plik)

with open('Random.txt', 'w') as plik:
    Score2, S2 = RandomMotifSearch(DNA, k, len(DNA))
    print(str(k) + '-mery', 'Score wynosi: ' + str(Score2), 'Motyw:', file=plik)
with open('Random.txt', 'a') as plik:
    for i, Sekwencja in enumerate(DNA):
        print(Sekwencja[S2[i]:S2[i] + k], file=plik)

with open('Brute.txt', 'w') as plik:
    Score3, S3 = BruteForceMotifSearch(DNA, len(DNA), len(DNA[0]), k)
    print(str(k) + '-mery', 'Score wynosi: ' + str(Score3), 'Motyw:', file=plik)
with open('Brute.txt', 'a') as plik:
    for i, Sekwencja in enumerate(DNA):
        print(Sekwencja[S3[i]:S3[i]+k], file=plik)
