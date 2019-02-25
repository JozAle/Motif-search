def Score(s, DNA, k):
    """
    compute the consensus SCORE of a given k-mer alignement given
    offsets into each DNA string.
    s = list of starting indices.
    DNA = list of nucleotide strings.
    k= Target Motif length
    """

    score = 0  # Zmienna score równa sie 0
    for i in range(k):  # Pętla iterująca po długości meru
        cnt = dict(zip('acgt', (0, 0, 0, 0)))  # Stworzenie słownika cnt składającego się z nukleotydów, ich liczby oraz
        # jego zresetowanie po każdym przejściu petli
        for j, s_val in enumerate(s):  # Pętla iterująca po każdej sekwencji
            base = DNA[j][s_val + i]  # Sczytanie nukleotydu na odpowiedniej pozycji w DNA
            cnt[base] += 1  # Zwiększenie liczby przy odpowiednim nukleotydzie w słowniku
        score += max(cnt.values())  # Dodanie najwyższej wartości ze słownika do wartości score
    return score  # Zwrócenie score


def MojScore(s, DNA, k):
    Alignement = []
    Profile = []
    Score = 0
    try:
        for i, Wiersz in enumerate(DNA):
            Alignement.append(Wiersz[s[i]:s[i] + k])
    except IndexError:
        print('Zmienna DNA jest pusta!')
    for j in range(0, k):
        Lista = []
        for i in range(0, 4):
            Lista.append(0)
        Profile.append(Lista)
    try:
        for j in range(0, k):
            for i in range(0, len(Alignement)):
                try:
                    if Alignement[i][j] == 'a':
                        Profile[j][0] += 1
                    elif Alignement[i][j] == 'c':
                        Profile[j][1] += 1
                    elif Alignement[i][j] == 'g':
                        Profile[j][2] += 1
                    else:
                        Profile[j][3] += 1
                    assert Alignement[i][j] == 'a' or Alignement[i][j] == 'c' or Alignement[i][j] == 'g' or \
                           Alignement[i][j] == 't', 'Zle znaki w DNA!'
                except TypeError:
                    print('Indeks i lub j powinien byc typu int!')
    except IndexError:
        print('Funkcja wyszla poza zakres listy Alignement lub lista Alignement jest pusta!')
    for Kolumna in Profile:
        Score += max(Kolumna)
    assert Score > 0, 'Funkcja nie zlicza Score!'
    return Score


def Profile(Motywy, DNA, k):
    Profile = []  # Stworzenie listy Profile
    for m in range(0, k):  # Petla iterujaca od 0 do dlugosci meru
        Lista = []  # Stworzenie listy Lista
        for n in range(0, 4):  # Petla iterujaca od 0 do 4
            Lista.append(0)  # Zapelniam Liste zerami
        Profile.append(Lista)  # Zapelniam liste Profile listami pelnymi zer
    for i in range(0, k):  # Petla iterujaca od 0 do dlugosci meru
        Ilosc_A = 0
        Ilosc_G = 0
        Ilosc_C = 0
        Ilosc_T = 0  # Zerowanie ilosci zliczanych znakow
        for j in range(0, len(Motywy)):  # Petla dzialajaca od 0 do konca meru
            if DNA[j][Motywy[j] + i] == 'a':
                Ilosc_A += 1
            elif DNA[j][Motywy[j] + i] == 'c':
                Ilosc_C += 1
            elif DNA[j][Motywy[j] + i] == 'g':
                Ilosc_G += 1
            else:  # Sprawdzannie czy znak jest a, c ,g lub t i nastepnie zwieksza odpowienia wartosc
                Ilosc_T += 1
        Profile[i][0] = Ilosc_A / len(DNA) + 1
        Profile[i][1] = Ilosc_C / len(DNA) + 1
        Profile[i][2] = Ilosc_G / len(DNA) + 1
        Profile[i][3] = Ilosc_T / len(DNA) + 1  # Wartosci sa dzielone przez ilosc sekwencji oraz jest do nich dodawane
        # 1 by umozliwic znalezienie merow innych niz identyczne
    return Profile  # Zwracany jest profil


def MostProbableFromProfile(Sekwencja, Profile, k):
    NajPrawdo = 0  # Stworzenie zmiennej przechowujacej wartosc najwyzszego prawdopodobienstwa
    for i in range(0, len(Sekwencja) - k):  # Petla iterujaca od 0 do dlugosci sekwencji - dlugosc meru
        Prawdo = 1  # Stworzenie zmiennej przechowujacej wartosc prawdopodbienstwa
        KMer = Sekwencja[i:i + k]  # Przypisanie aktualnie sprawdzanego meru
        for j, Nukleotyd in enumerate(KMer):  # Petla iterujaca od poczatku do konca zmiennej KMer
            if Nukleotyd == 'a':
                Prawdo = Prawdo * Profile[j][0]
            elif Nukleotyd == 'c':
                Prawdo = Prawdo * Profile[j][1]
            elif Nukleotyd == 'g':
                Prawdo = Prawdo * Profile[j][2]
            else:  # Sprawdzanie znaku oraz zwiekszanie prawdopodobienstwa o odpowiednia wartosc
                Prawdo = Prawdo * Profile[j][3]
        if Prawdo > NajPrawdo:  # Jesli prawdopodobienstwo jest wyzsze niz najwyzsze
            NajPrawdo = Prawdo  # Zmiennej najwyzszej przypisywana jest wartosc prawdopodobienstwa
            NajKMer = i  # Zmienna NajKMer zapisuje indeks meru najwyzszego prawdopodobienstwa
    return NajKMer  # Zwracany jest indeks o najwyzszym prawdopodobienstwie


if __name__ == '__main__':
    DNA = ('tagtggtcttttgagtgtagatctgaagggaaagtatttccaccagttcggggtcacccagcagggcagggtgacttaat',
           'cgcgactcggcgctcacagttatcgcacgtttagaccaaaacggagttggatccgaaactggagtttaatcggagtcctt',
           'gttacttgtgagcctggttagacccgaaatataattgttggctgcatagcggagctgacatacgagtaggggaaatgcgt',
           'aacatcaggctttgattaaacaatttaagcacgtaaatccgaattgacctgatgacaatacggaacatgccggctccggg'
           )

    DNA2 = ('cctgatctgccctaacctacaggtcgatccgaaattcg',
            'tcttttcgagtccttgtacctccatttgctctgatgac',
            'atttgctccactacaagaggctactgtgtagatccgta',
            'gtaataatggctcagccatgtcaatgtgcggcattcca'
            )

    S1 = [2, 1, 10, 5]
    k = 8
    print('Dzialanie Score: ' + str(Score(S1, DNA, k)))
    print('Dzialanie MojScore: ' + str(MojScore(S1, DNA, k)))
    print('Profil dla motywu: ')
    for i, Sekwencja in enumerate(DNA):
        print(Sekwencja[S1[i]:S1[i] + 8])
    print('to: ' + str(Profile(S1, DNA, k)))
    NajPra = []
    for Sekwencja in DNA2:
        NajPra.append(MostProbableFromProfile(Sekwencja, Profile(S1, DNA, k), k))
    print('Najbardziej podobny do tego profilu z DNA2 jest motyw: ')
    for i, Sekwencja in enumerate(DNA2):
        print(Sekwencja[NajPra[i]:NajPra[i] + 8])

    print('Jej profil to: ')
    print(Profile(NajPra, DNA2, k))
