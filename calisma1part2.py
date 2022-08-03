import pandas as pd
import numpy as np
from deneme import adict

#Needed names
continents = ["AFR","AMR","EAS","EUR","SAS"]
rslerlist = [ "rs1815739", "rs4644994", "rs4253778", "rs1042713", "rs8192678", "rs2070744", "rs11549465", "rs1800012", "rs12722", "rs1049434",
        "rs1800795", "rs6265", "rs4680"]

#Dataframes created for different ethnicities
AFR = pd.DataFrame()
AMR = pd.DataFrame()
EAS = pd.DataFrame()
EUR = pd.DataFrame()
SAS = pd.DataFrame()

#Dataframes filled with according data
for i in rslerlist:
    for j in continents:
        a = pd.read_csv("data/ethnicity/{}/{}.csv".format(i,j),names=rslerlist)
        a = a.iloc[:,1]
        #print(a)
        if j == "AFR":
            AFR = pd.concat([AFR,a],axis=1)
            #print(AFR)
        elif j == "AMR":
            AMR = pd.concat([AMR,a], axis=1)
        elif j == "EAS":
            EAS = pd.concat([EAS, a], axis=1)
        elif j == "EUR":
            EUR = pd.concat([EUR,a], axis=1)
        elif j == "SAS":
            SAS = pd.concat([SAS,a], axis=1)

#Column names changed to correct names
AFR.columns = rslerlist
AMR.columns = rslerlist
EAS.columns = rslerlist
EUR.columns = rslerlist
SAS.columns = rslerlist


#print(AFR.rs1815739.unique())

#print(AFR.isna().sum().sum())
#Frquency calculator for homo, hetero and genotype genes

total_number = []
for df in [AFR, AMR, EAS, EUR, SAS]:
    total_number.append(df.count(axis=1))

AFR_tot = len(total_number[0])
AMR_tot = len(total_number[1])
EAS_tot = len(total_number[2])
EUR_tot = len(total_number[3])
SAS_tot = len(total_number[4])

AFR_values = []
AMR_values = []
EAS_values = []
EUR_values = []
SAS_values = []

for i in rslerlist:

    AFR_values.append(AFR[i].value_counts().to_dict())
    AMR_values.append(AMR[i].value_counts().to_dict())
    EAS_values.append(EAS[i].value_counts().to_dict())
    EUR_values.append(EUR[i].value_counts().to_dict())
    SAS_values.append(SAS[i].value_counts().to_dict())

#print(SAS_values)
values_list = [AFR_values, AMR_values, EAS_values, EUR_values, SAS_values]

heterovals = [["A|C","C|A"],["A|T","T|A"],["A|G","G|A"],["C|T","T|C"],["G|C","C|G"],["T|G","G|T"]]

for i in values_list:
    for b in i:
        for j, k in heterovals:
            if j in b.keys():
                if k in b.keys():
                    new_val = b[j] + b[k]
                    del b[j], b[k]
                    new_row = {"{}".format(j): new_val}
                    b.update(new_row)

def freq_calc(df_values, tot):

    ''' dataframe, total number '''

    df_values_frequency = []
    df_values_allel_frequency = []

    for b in df_values:
        df_values_frequency_rs = {}
        df_values_allel_rs1 = {}
        df_values_allel_rs2 = {}
        for j in ["A", "C", "T", "G"]:
            for k in ["A", "C", "T", "G"]:

                if "{}|{}".format(j, k) in b.keys():

                    df_values_frequency_rs.update({"{}|{}".format(j, k): ((b["{}|{}".format(j, k)] / tot) * 100)})
                    if j == k:
                        df_values_allel_rs1.update({ j: ((b["{}|{}".format(j, k)] * 2) / tot) *50})
                    elif j != k:
                        df_values_allel_rs2.update({j: (b["{}|{}".format(j, k)] / tot) *50})
                        df_values_allel_rs2.update({k: (b["{}|{}".format(j, k)] / tot) *50})
        df_values_allel_rs = adict(df_values_allel_rs1, df_values_allel_rs2)
        df_values_allel_frequency.append(df_values_allel_rs)
        df_values_frequency.append(df_values_frequency_rs)

    return df_values_frequency, df_values_allel_frequency


AFR_frequency, AFR_allel_frequency = freq_calc(AFR_values, AFR_tot)
AMR_frequency, AMR_allel_frequency = freq_calc(AMR_values, AMR_tot)
EAS_frequency, EAS_allel_frequency = freq_calc(EAS_values, EAS_tot)
EUR_frequency, EUR_allel_frequency = freq_calc(EUR_values, EUR_tot)
SAS_frequency, SAS_allel_frequency = freq_calc(SAS_values, SAS_tot)
print(SAS_allel_frequency)







