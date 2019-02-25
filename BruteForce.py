from Score import Score

def BruteForceMotifSearch(DNA, t, n, l):
    BestScore = 0  # Zmienna przechowująca najwyzsza wartosc Score
    S = []  # Lista przechowujaca indeksy pozycji motywu
    for i in range(0, t):  # Zapelnienie S zerowymi indeksami na poczatek
        S.append(0)
    for i in range(0, pow((n-l), t) - 1):  # Petla dzialajaca tak dlugo az sprawdzi wszystkie mozliwosci
        for j in range(0, len(S)):  # Petla iterujaca po kazdym indeksie w S
            if S[j] == n-l:  # Jesli ktorys indeks osiagnie wartosc dlugosci sekwencji
                S[j] = 0  # Wyzeruj go
                S[j+1] += 1  # Zwieksz o 1 indeks o 1 wyzszy
        if Score(S, DNA, l) > BestScore:  # Jesli Score nowego motywu jest wyzszy niz bestscore
            BestScore = Score(S, DNA, l)  # Bestscore = Score nowego motywu
            BestMotif = S  # Zapisuje Motyw jako najlepszy
        S[0] += 1  # Zwieksza wartosc indeksu pierwszego o 1
    return BestScore, BestMotif  # Zwraca Najwyzszy score oraz Najlepszy motyw

if __name__ == '__main__':
    DNA = (
        'agtggtcttttgagtgtagat',
        'cgcgactcggcgctcacacat',
        'gaaatggttcggtgcgatatc',
        'aaatccgaattgacctgatga',
        'cggtgcgatatcaggccggtc',
        'agctgacatacgagtagggac'
    )

    Score, S = BruteForceMotifSearch(DNA, len(DNA), len(DNA[0]), 15)
    print(Score)
    for i, Sekwencja in enumerate(DNA):
        print(Sekwencja[S[i]:S[i]+15])
