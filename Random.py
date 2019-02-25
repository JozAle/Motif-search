from Score import Score, Profile, MostProbableFromProfile
import random


def RandomMotifSearch(DNA, k, t):
    Motifs = []  # utworzenie listy Motifs
    for i in range(0, t):  # Petla iterujaca od 0 do ilosci sekwencji
        Indeks = random.randint(0, len(DNA[0]) - k)  # losowanie zmiennej indeks w przedziale od 0 do konca sekwencji -
        # dlugosc meru
        Motifs.append(Indeks)  # Dodanie wylosowanej liczby do listy Motifs
    BestMotifs = Motifs
    while True:  # Niekonczaca sie petla
        Profil = Profile(Motifs, DNA, k)  # Tworzenie profilu na podstawie listy motifs
        for n in range(0, t):  # Petla od 0 do ilosci sekwencji
            Motifs[n] = MostProbableFromProfile(DNA[n], Profil, k)  # Znalezienie najbardziej podobnych motywow na
            # podstawie profilu
        if Score(Motifs, DNA, k) > Score(BestMotifs, DNA, k):  # Sprawdzenie czy Score Motifs jest wyzszy od
            # wylosowanego
            BestMotifs = Motifs  # Jesli tak to BestMotifs = Motifs
        else:
            return Score(BestMotifs, DNA, k), BestMotifs  # Jesli nie to zwracany jest Score oraz BestMotifs


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

    Score, S = RandomMotifSearch(DNA, 3, len(DNA))
    print(Score)
    for i, Sekwencja in enumerate(DNA):
        print(Sekwencja[S[i]:S[i] + 3])
