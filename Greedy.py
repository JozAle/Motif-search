from Score import Score, Profile, MostProbableFromProfile


def GreedyMotifSearch(DNA, k, t):
    BestMotifs = []
    Motifs = []
    for Sekwencja in DNA:  # Dla kazdej sekwencji w DNA
        BestMotifs.append(0)
        Motifs.append(0)  # Uzupelnij listy Motifs oraz BestMotifs zerami
    for j in range(0, len(DNA[0]) - k):  # Petla od 0 do konca sekwencji - dlugpsc meru
        Motifs[0] = j  # Każdy kolejny mer w pierwszej sekwencji
        for i in range(1, t):  # Pętla od drugiej sekwencji do ostatniej
            Profil = Profile(Motifs, DNA, k)  # Utworzenie profilu z motywow
            Motifs[i] = MostProbableFromProfile(DNA[i], Profil, k)  # Dodanie nowego meru na podstawie profilu
        if Score(Motifs, DNA, k) > Score(BestMotifs, DNA, k):  # Sprawdzenie czy nowy motyw ma wyzszy score niz
            #poprzedni najwyzszy
            BestMotifs = Motifs  # Zapisanie motywu jako najlepszego
    return Score(BestMotifs, DNA, k), BestMotifs  # Zwrocenie Score najlepszego motywu oraz najlepszego motywu

if __name__ == '__main__':
    DNA = (
        'tagtgg',
        'cgcgac',
        'aacatc',
        'tagatt',
        'gaaatg',
        'ttctta',
        'ctacct'
    )

    Score, S = GreedyMotifSearch(DNA, 3, len(DNA))
    print(Score)
    for i, Sekwencja in enumerate(DNA):
        print(Sekwencja[S[i]:S[i]+3])
