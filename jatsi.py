from random import randint

def heita_noppia() -> list:
    lista = []
    for i in range(5):
        noppa = randint(1,6)
        lista.append(noppa)
    lista.sort()
    return lista

def heita_uudelleen(lista: list[int]) -> list[int,int]:
    lista3 = []
    for index in range(2):
        print(f"Noppien silmäluvut ensimmäisen heiton jälkeen: {lista}")
        print("Mitkä nopat haluat pitää? (0 heittää kaikki nopat uudelleen)")
        print("(Anna silmäluvut pilkulla erotettuna.)")
        luvut = input("silmäluvut: ")
        if luvut == "0":
            lista = heita_noppia()
            lista = lista[0:5-len(lista3)]
        else:
            luvut = luvut.split(",")
            for lukuja in luvut:
                lukuja = int(lukuja)
                lista3.append(lukuja)
            lista = heita_noppia()
            lista = lista[0:5-len(lista3)]
            if len(lista) == 0:
                break
            if index < 1:
                print("Heitetään uudelleen noppia, joita ei valittu.")
    lista4 = lista3 + lista
    lista4.sort()
    print("Viimeisen heiton tulos:", lista4)
    return lista4
    
def taulukko(pistelasku: dict):
    print("Taulukko:")
    print("ykköset", pistelasku.get("ykköset",""))
    print("kakkoset", pistelasku.get("kakkoset",""))
    print("kolmoset", pistelasku.get("kolmoset",""))
    print("neloset", pistelasku.get("neloset",""))
    print("viitoset", pistelasku.get("viitoset",""))
    print("kuutoset", pistelasku.get("kuutoset",""))
    print("yksi pari", pistelasku.get("yksi pari",""))
    print("kaksi paria", pistelasku.get("kaksi paria",""))
    print("kolmoisluku", pistelasku.get("kolmoisluku",""))
    print("neloisluku", pistelasku.get("neloisluku",""))
    print("pieni suora", pistelasku.get("pieni suora",""))
    print("suuri suora", pistelasku.get("suuri suora",""))
    print("täyskäsi", pistelasku.get("täyskäsi",""))
    print("sattuma", pistelasku.get("sattuma",""))
    print("yatzy", pistelasku.get("yatzy",""))

def annettavat_pisteet(yhdistelma: list, pistelasku: dict) -> int:
    index = 0
    while index < 1:
        try:
            while True:
                merkitaan = input("Mihin yhdistelmään tulos merkitään: ")
                if pistelasku.get(merkitaan,"") == "":
                    break
                else:
                    print("Yhdistelmä on valittu jo aikaisemmin!")
            if merkitaan == "ykköset":
                pisteet = yhdistelma.count(1)
            if merkitaan == "kakkoset":
                pisteet = yhdistelma.count(2)*2
            if merkitaan == "kolmoset":
                pisteet = yhdistelma.count(3)*3
            if merkitaan == "neloset":
                pisteet = yhdistelma.count(4)*4
            if merkitaan == "viitoset":
                pisteet = yhdistelma.count(5)*5
            if merkitaan == "kuutoset":
                pisteet = yhdistelma.count(6)*6
            if merkitaan == "yksi pari":
                pisteet = 0 
                for i in range(1,7):
                    if yhdistelma.count(i) >= 2:
                        pisteet = i *2
            if merkitaan == "kaksi paria":
                pisteet = 0
                kp_lista = []
                for i in range(1,7):
                    if yhdistelma.count(i) >= 2:
                        pari = i*2
                        kp_lista.append(pari)
                if len(kp_lista) >= 2:
                    pisteet = kp_lista[-1] + kp_lista[-2]
            if merkitaan == "kolmoisluku":
                pisteet = 0 
                for i in range(1,7):
                    if yhdistelma.count(i) >= 3:
                        pisteet = i *3 
            if merkitaan == "neloisluku":
                pisteet = 0
                for i in range(1,7):
                    if yhdistelma.count(i) >= 4:
                        pisteet = i *4 
            if merkitaan == "pieni suora":
                pisteet = 0 
                if [1,2,3,4,5] in yhdistelma:
                    pisteet = 15 
            if merkitaan == "suuri suora":
                pisteet = 0
                if [2,3,4,5,6] in yhdistelma:
                    pisteet = 20
            if merkitaan == "täyskäsi":
                pisteet = 0
                tk_lista = []
                for i in range(1,7):
                    if yhdistelma.count(i) == 3:
                        kolmonen = i * 3
                        tk_lista.append(kolmonen)
                    elif yhdistelma.count(i) == 2:
                        kakkonen = i * 2
                        tk_lista.append(kakkonen)
                if len(tk_lista) == 2:
                    pisteet = tk_lista[-1] + tk_lista[-2]
            if merkitaan == "sattuma":
                summa = 0
                for luku4 in yhdistelma:
                    summa += luku4
                pisteet = summa
            if merkitaan == "yatzy":
                pisteet = 0
                for i in range(1,7):
                    if yhdistelma.count(i) == 5:
                        pisteet = 50
            pistelasku[merkitaan] = pisteet
            index = 1
            return pisteet
        except:
            print("kirjoitit väärin yhdistelmän nimen! Yritä uudestaan:")

def loppupisteet(pistelasku: dict) -> int:
    summa1 = pistelasku["ykköset"]+pistelasku["kakkoset"]+pistelasku["kolmoset"]
    summa2 = pistelasku["neloset"] + pistelasku["viitoset"] + pistelasku["kuutoset"]
    summa3 = pistelasku["yksi pari"] + pistelasku["kaksi paria"] + pistelasku["kolmoisluku"]
    summa4 = pistelasku["neloisluku"] + pistelasku["pieni suora"] + pistelasku["suuri suora"]
    summa5 = pistelasku["täyskäsi"] + pistelasku["sattuma"] + pistelasku["yatzy"]
    if summa1 + summa2 >= 63:
        bonus = 50
    else:
        bonus = 0
    pisteet_lopuksi = summa1+ summa2 + summa3 + summa4 + summa5 + bonus
    return pisteet_lopuksi

def tiedoston_avaus_tai_luonti():
    try:
        with open("yatzy_aiemmat_tulokset.txt") as tiedosto:
            sisalto = tiedosto.read()
            print("Aiempien pelien tulokset.")
            print(sisalto)
    except FileNotFoundError:
        print("Sinulla ei ole aiempien pelien Tuloksia tallennettuna!")
        print("Luodaan tiedosto ""yatzy_aiemmat_tulokset.txt"", johon tallenetaan tulokset pelin päätyttyä.")
        with open("yatzy_aiemmat_tulokset.txt","w") as tiedosto:
            tiedosto.write("Tulokset vanhimmasta uusinpaan:" + "\n")

def tallenna_tulos_tiedostoon(loppu_pisteet: int):
    with open("yatzy_aiemmat_tulokset.txt","a") as tiedosto:
        tiedosto.write(str(loppu_pisteet) + "\n")



tiedoston_avaus_tai_luonti()
print("---------------------------------")
print("Aloitetaan peli!")
pistelasku = {}
while len(pistelasku) < 15:
    nopat = heita_noppia()
    uudelleen_heitot = heita_uudelleen(nopat)
    taulukko(pistelasku)
    x = annettavat_pisteet(uudelleen_heitot, pistelasku)
    print(f"Merkataan siihen kohtaan {x} pistettä.")
print("Peli päättyi, taulukko on täynnä")
pisteet_lopuksi = loppupisteet(pistelasku)
print(f"Bonuspisteet huomioon ottaen yhteenlaskettu pistemäärä on {pisteet_lopuksi}.")
print("Tallennetaan tulos tiedostoon")
tallenna_tulos_tiedostoon(pisteet_lopuksi)
