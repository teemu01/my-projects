#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#*******************************
# Reset memory
#*******************************
from IPython import get_ipython
get_ipython().run_line_magic('reset', '-sf') 

import numpy as np
import matplotlib.pyplot as plt
import math
from ase.visualize import view
from ase.io import read
from numpy.linalg import inv
from ase.build import make_supercell
from ase import Atoms
from numba import njit
import pandas as pd



def read_cell_file(datafile):
    """
    Funktio, joka lukee cell-tiedoston.
    
    Parametrit:
        datafile: String
            Tiedoston polku.
            
    palauttaa: np.ndarray
        Sisältää cell-tiedoston kaikki arvot.
    """
    
    with open(datafile, "r") as file:
        cell = file.read()
        cell = cell.split("\n")
        cell = np.array([c.split() for c in cell[1:-1]]).astype(float)
        return cell



@njit
def Sf_histogram(grid, q_val, super_pos, particles):
    """
    Funktio, joka luo rakennetekijä-histogrammin yksittäiselle q-arvolle, missä X-akselilla on etäisyys ja Y-akselilla etäisyysväleihin
    levittäytyneet rakennetekijän arvot.
    
    Parametrit:
        grid: np.ndarray
            Histogrammin X-akselin arvot, joka kuvaa atomiparien etäisyyksiä.
        q_val: int
            Rakennetekijässä esiintyvä yksittäinen q-arvo.
        super_pos: np.ndarray
            Supercellin kaikkien hiukkasten paikat avaruudessa vektoreina esitettynä.
        particles: int
            Supercellin hiukkasten lukumäärä.
    Palauttaa: np.ndarray
        Rakennetekijä-histogrammi yksittäiselle q-arvolle, mihin listätty arvot niille kuuluville etäisyysväleille.
    """
    


    grid_dr = grid[1]-grid[0]
    
    #supercellin hiukkasten lukumäärä
    N = super_pos.shape[0]
    
    q_D = super_pos*q_val
    histo = np.zeros_like(grid)
    
    #Iteroi läpi jokaisen alkuperäisen cellin hiukkaset
    for i in range(0, particles):
        #iteroi läpi jokaisen supercellin hiukkaset
        for j in range(0,N):
            if not i == j:
                norm = 0.0
                for ic in range(3):
                    # Kahden atomin välinen etäisyys Pythagoraan lauseella
                    norm += (q_D[i,ic]-q_D[j,ic])**2
                    
                dis = np.sqrt(norm)
                
                #laskee yhden indeksiaron rakennetekijälle
                increment = (1/particles)*np.sin(dis)/(dis)
                if int((dis/q_val)/grid_dr) > len(histo):
                    raise Exception(f"Gridin suurin arvo liian pieni: kahden hiukkasen välinen etäisyys iteraatiossa on {int(dis/q_val)}")
                #lisää arvon histogrammiin oikeaan etäisyysväliin                    
                histo[int((dis/q_val)/grid_dr)] += increment
                
    return histo
            




def structure_factor(traj, cell, supercell_size, q_values, snap):
    """
    Funktio, joka muodostaa supercellin annetusta trajektorista yhdelle snapshotille
    ja laskee rakennetekijä-histogrammin supercellistä kaikille q-arvoille.
    
    Parametrit:
        traj: np.ndarray
            trajektori-data
        cell: np.ndarray
            cell-data
        supercell_size: int
            supercellin koko, esim. parametrin arvo 3 tekee 3x3x3 supercellin.
            Luvun oltava pariton, jotta alkuperäinen cell on supercellin keskellä!
        snap: int
            kokonaisluku, joka kertoo trajektorin snapshotin numeron.
        q_values: list
            Lista q-arvoista, joista rakennetekijän arvo riippuu.
    palauttaa: np.ndarray, np.ndarray
        Histogramissa esiintyvä grid-muuttuja, jota käytetään painofunktiossa,
        sekä np.ndarray rakennetekijä-histogrammeista jokaiselle q-arvolle.
    """
    
    #tarkastaa onko supercellin koko pariton
    if supercell_size % 2 == 0:
        raise ValueError("Supercellin koko ei ole pariton")
        
    #valitsee yksittäisen snaphsotin trajektorista ja cell-datasta
    traj_frame = traj[snap]
    cell_frame = cell[snap]
    #cellvektorit
    a = cell_frame[2:5]
    b = cell_frame[5:8]
    c = cell_frame[8:11]

    #Vaihteluväli supercellin koolle
    s_max_val = int((supercell_size+1)/2)
    s_min_val = int((1-supercell_size)/2)

    particles = len(traj_frame)
    supercell = traj_frame.positions.copy()
    #Lisää alicellejä supercell-muuttujaan iteroiden.
    for i in range(s_min_val,s_max_val):
        for j in range(s_min_val,s_max_val):
            for k in range(s_min_val,s_max_val):
                if i == 0 and j == 0 and k == 0:
                    continue
                else:
                    supercell = np.vstack((supercell,traj_frame.positions + a*i+b*j+c*k))
    
    #histogrammin suurin mahdollinen arvo. On oltava suurempi kuin
    #atomien suurin mahdollinen etäisyys rakennetekijä yhtälössä.
    cell_diagonal = supercell_size*(12+12+12)
    grid = np.arange(0, cell_diagonal/2, 0.1)
    
    #Lisää rakennetekijä-histogrammit jokaiselle eri q-arvolle listaan.
    sf_sum = []
    for qi in q_values:
        print(f"q-val: {qi}")
        histo = Sf_histogram(grid, qi, supercell, particles)
        sf_sum.append(histo)
    
    sf_sum = np.array(sf_sum)
    return sf_sum,grid
    


#painofunktio. Eksponentin nimittäjässä olevaa arvoa tulee muuttaa riippuen supercellin koosta.
def weight_function(grid: list):
    """
    Funktio muodostaa painofunktion, jolla säädellään rakennetekijän arvoja.
    
    parametrit:
        grid: np.ndarray
            Tässä käytetään samaa grid-muuttuujaa kuin minkä structure_factor-funktio palauttaa
    palauttaa: np.ndarray
        Array painofunktion arvoista grid-muuttujassa esiintyvien etäisyyksien suhteen.
    """
    
    weight = []
    for r in grid:
        #Eksponentissa olevaa nimittäjää voi muuttaa siten, että painofunktion arvo ei mene nollaan liian nopeasti.
        f = math.exp(-0.5*(r+0.05)**2/100**2)
        weight.append(f)

    weight = np.array(weight)
    return weight



"""
Tässä skriptissä lasketaan rakennetekijä atoimeille usealle snapshotille
ja näistä otetaan keskiarvo, joka esitetään kuvaajassa.
Tässä rakennetekijä johdetaan muunnellusta Debyen sironta kaavaasta.
Sen sijaan, että Debyen kaavassa olisi mukana jokainen supercellissä esiintyvä atomipari,
niin siinä on vain atomiparit, joissa vähintään toinen atomi kuuluu alkuperäiseen celliin.
Alkuperäisestä celllistä joudutaan tekemään supercell, jotta diffraktiopiikit selvästi erottuisivat
kuvaajassa.

Tässä sovelletaan myös painofunktiota rakennetekijään. Jokaiselle q-arvolle lasketaan histogrammi, josta
nähdään rakennetekijän jakautumuminen eri etäisyyksillä yksittäiselle q-arvolle. Jokainen histogrammi
kerrotaan painofunktiolla, jolloin atomiparien etäisyyden merkitys muuttuu, niin, että kaukana olevien
atomiparien kontribuutio lähes mitätöityy rakennetekijässä ja lähellä olevien kasvaa. Näin supercellin reunoissa
olevat atomit eivät vaikuta tulokseen.


Lämpövärähtelyn seurauksena hiukkaset värähtelevät tasapainoasemansa ympärillä. Luultavasti tämä leventää
rakennetekijän diffraktiopiikkejä. Tämän lisäksi kuvaajassa diffraktiopiikit näyttävät liikkuvan jaksollisesti
ajan suhteen, mikä luultavasti johtuu cellin tilavuuden muutoksista. Tämä myös vaikuttaa siihen, että keskiarvossa
piikit ovat leveämpiä. Pienemmissä simulaatioissa rakennetekijä muuttuu enemmän ajan suhteen esimerkiksi tämä
näkyy 48 molekyylin 30GPa simulaatiossa. Lisäksi välillä piikkien määrä muuttuu. Tämä antaa keskiarvolle epätarkan tuloksen.
64 molekyylin 10GPa simulaatiossa sen sijaa rakennetekijän muutos on pienempää, mutta tähän saattaa myös vaikuttaa se, että
cell on kuutiomainen eli symmetrisempi kuin 48 molekyylin simulaation tetraedri.
"""





"""
Koodista tarvitsee muuttaa tarvittaessa:
    trajektori- ja cell-tiedostojen polut:
        traj- ja cell-muuttujat, jotka voi muuttaa pääohjelmassa
    paineen arvo:
        pressure-muuttuuja, jonka voi muuttaa pääohjelmassa
    supercellin koko:
        supercell_size, jonka voi muuttaa pääohjelmassa
    painofunktion eksponentin arvo:
        jonka voi muuttaa weight_function-funktiossa
"""


#simulaatiosta saadut datat

traj = read("/home/turbo/Desktop/N2/64_molecules/10Gpa/npt_coords.xyz", index=":") 
cell = read_cell_file("/home/turbo/Desktop/N2/64_molecules/10Gpa/npt_cell.txt")

# traj = read("/home/turbo/Desktop/N2/96_molecules/npt_coords.xyz", index=":") 
# cell = read_cell_file("/home/turbo/Desktop/N2/96_molecules/npt_cell.txt")


expt_e = pd.read_csv("/home/turbo/Desktop/visual_studio_code/epsilon_mittaus.csv",header = None)
expt_d = pd.read_csv("/home/turbo/Desktop/visual_studio_code/delta_data.csv",header = None)


#Muuttaa kokeellisen tuloksen, johon verrataan, ja myös kuvaajan otsikon
pressure = 10

#96molekyyliä 30GPa
traj = traj[9000:]
cell = cell[9000:]
traj = traj[0::500]
cell = cell[0::500]


#Datan rajaus eri simulaatiotuloksille
# #64molekyyliä 10GPa
# traj = traj[10000:20000]
# cell = cell[10000:20000]
# traj = traj[0::500]
# cell = cell[0::500]


#64molekyyliä 30Gpa
# traj = traj[15000:]
# cell = cell[15000:]
# traj = traj[0::500]
# cell = cell[0::500]

#48molekyyliä 30Gpa
# traj = traj[5000:]
# cell = cell[5000:]
# traj = traj[0::200]
# cell = cell[0::200]

#48molekyyliä 10Gpa
# traj = traj[5000:]
# cell = cell[5000:]
# traj = traj[0::500]
# cell = cell[0::500]



#antaa koon supercellille. Esim. 17 tekee 17x17x17 supercellin. Luvun on oltava pariton
#jotta alkuperäinen cell olisi supercellin keskellä.
supercell_size = 11
q_values = np.linspace(2.0,3.4,100)



#Laskee rakennetekijän jokaiselle snapshotille.
#Grid-muuttujan pitäisi olla sama jokaisella snapshotilla
aver_sf = []
for snap in range(0,len(cell)):
    print(f"snapshot: {snap}")
    sf_sum,grid = structure_factor(traj,cell, supercell_size, q_values, snap)
    aver_sf.append(sf_sum)
aver_sf = np.array(aver_sf)




#Kaikkien snapshottien histogrammeista otetaan keskiarvo erikseen jokaiselle q-arvolle ja nämä kerrotaan painofunktiolla.
#Tämän jälkeen jokaiselle histogrammille erikseen summataan sen arvot yhteen.
sf_plot = 1 + np.sum(weight_function(grid) * np.average(aver_sf,axis=0),axis = 1)
sf_plot_max = np.max(sf_plot)

#kuvaajan plottaus Rakennetekijän keskiarvosta
plt.figure()
plt.plot(q_values,sf_plot,label= "Sim.")
if pressure == 10:
    plt.plot(expt_d[0][0:61],expt_d[1][0:61]*(sf_plot_max/np.max(expt_d[1])),color="black",label="Expt.")  #delta
else:
    plt.plot(expt_e[0][0:61],expt_e[1][0:61]*(sf_plot_max/np.max(expt_e[1])),color="black",label="Expt.")  #epsilon
plt.legend()
plt.grid()
plt.xlabel(r"$q$ [Å$^{-1}]$")
plt.ylabel("$S(q)$ [arb. units]")
plt.title(f"{pressure} GPa")



plt.figure(figsize=(15,15))
# jokaisen snapshotin rakennetekijän plottaus yhteen kuvaajaan.
for j in range(0,len(aver_sf)):
    if j == 0:
        plt.plot(q_values, j*5+1+np.sum(aver_sf[j],axis=1), label=f"{cell[0,1]} fs")
    else:
        plt.plot(q_values, j*5+1+np.sum(aver_sf[j],axis=1))
plt.title(f"time interval {cell[1,1] - cell[0,1]} fs")
plt.grid()
plt.legend()



#np.save("Atom_center_10GPa_time_aver.npy",np.array([q_values,sf_plot]))