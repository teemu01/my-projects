# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 13:09:58 2024

@author: OMISTAJA
"""



import numpy as np
from scipy.stats import shapiro
from scipy.stats import ttest_ind
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


#Lukee tiedoston ja valitsee muuttuujat, joita pitää tutkia
df = pd.read_csv(r"C:\Users\OMISTAJA\Documents\visual_studio\python\dataAnalyysi/habits.data",sep=";",na_values=["?"])
variables = ["kohde", "jasen", "pvknro", "sp", "IKAL1", "ASALUE", "V33", "V5", "V68", "V69", "H1b_A", "H1h_A"]
df = df[variables]



categorical = ["kohde", "jasen", "pvknro", "sp", "IKAL1", "ASALUE","H1b_A", "H1h_A"]
quantitative = ["V33", "V5", "V68","V69"]


for col in categorical:
    print(col, end="\n\n")
    print(df[col].value_counts().sort_index(), end="\n\n")





df.loc[df.H1b_A > 2, "H1b_A"] = 2
df.loc[df.H1h_A > 2, "H1h_A"] = 2
df.loc[df.H1b_A  == 0,"H1b_A"] = 1
df.loc[df.H1h_A == 0, "H1h_A"] = 1



#Change all time to minutes in all household columns.

def time_change(hours):
    if isinstance(hours, str) and ":" in hours:  
        val=hours.split(":")
        min=int(val[1])
        hour=int(val[0])
       
        min_total=hour*60+min
        return int(min_total)
    if pd.isna(hours):
        return hours
    else:
        return int(hours)

#Here we apply the function “time_change” to the columns. This changes the actual time.

for col in ["V33", "V5", "V68", "V69"]:        
    df[col]=df[col].apply(time_change)
    df[col] = df[col].astype(np.float64)




# iteroi jokaisen 6-10 kolumnin alkion ja muuttaa luvun kokonaislukuksi, jonka yksikkö on minuutti
for i in range(len(df.index)):
    for j in range(6,10):
        try:
            df.iloc[i,j] = np.int64(int(df.iloc[i,j]))
        except:
            if type(df.iloc[i,j]) == float:
                continue
            else:
                a,b = df.iloc[i,j].split(":")
                df.iloc[i,j] = np.int64((int(a)*60+int(b)))


#koko dataframe int64-tyypiksi
df = df.astype(np.float64)

#muuttujien esiintymistiheyksiä
#(Kaikki Rivit vielä mukana)
freq = df.groupby(["sp", "ASALUE", "IKAL1"]).size()


testi = df.groupby(["H1b_A", "H1h_A"]).size()

df.iloc[:,2:5].value_counts().plot.bar()

# for i in categorical[1:-2]:
#     df[i] = pd.Categorical(df[i])


#plottausta
pd.plotting.scatter_matrix(df[quantitative])


df = df.dropna(subset=["V33","V5","V68","V69"])


#tiivistelmä datasta
summarize = df.describe()


#2, 3 ja 4 tehtävät
house = df.iloc[:,[0,6,7,8,9]].groupby(["kohde"],as_index=False,sort=False).sum().mean(axis=0).iloc[1:]
day_of_week = df.iloc[:,[2,6,7,8,9]].groupby(["pvknro"]).mean()
gender = df.iloc[:,[3,6,7,8,9]].groupby(["sp"]).mean()




#p-arvot alle 0.05, joten ei voida olettaa olevan normaaleja

pval_V33 = shapiro(df.V33).pvalue
pval_V5 = shapiro(df.V5).pvalue
pval_V68 = shapiro(df.V68).pvalue
pval_V69 = shapiro(df.V69).pvalue

pval_df = pd.DataFrame({"V33":[pval_V33],"V5":[pval_V5],"V68":[pval_V68],"V69":[pval_V69]})


counts = df.iloc[:,[-1,-2]].apply(pd.Series.value_counts)
from scipy.stats import chi2_contingency

chisquare_test = chi2_contingency(counts)[1]

pearson_CB= scipy.stats.pearsonr(df.V33,df.V5)[1]

c_mat1=df.iloc[:,[6,7,8,9]].corr("spearman")




correlation_spearman, pval_spearman = spearmanr(df.iloc[:,6:10])

correlation_spearman = pd.DataFrame(correlation_spearman, columns = ["Dining", "Cooking", "Reading", "Radio"])
pval_spearman = pd.DataFrame(pval_spearman, columns = ["Dining", "Cooking", "Reading", "Radio"])


from scipy.stats import wilcoxon

from scipy.stats import mannwhitneyu



weekday = df.loc[df.pvknro == 1] 
weekend = df.loc[df.pvknro == 2]

utests = mannwhitneyu(weekday.iloc[:,6:10],weekend.iloc[:,6:10])[1]

utest_df = pd.DataFrame( columns=["Dining", "Cooking", "Reading", "Radio"])
utest_df.loc[0] = utests

df[["V33","V5","V68","V69","H1b_A","H1h_A"]].corr("pearson")