
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
from sklearn import tree
import csv

# Input data files are available in the "../input/" directory.
# For example, running this (by clicking run or pressing Shift+Enter) will list all files under the input directory

import os
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn import metrics
import copy 
from sklearn.ensemble import RandomForestClassifier  

# for dirname, _, filenames in os.walk('/kaggle/input'):
#     for filename in filenames:
#         print(filename) 
#         print(os.path.join(dirname, filename))


train = pd.read_csv (r'train.csv')
test = pd.read_csv (r'test.csv')
sample = pd.read_csv (r'sample.csv')



sequences = train["Sequence"]
Lable = train["Lable"]
ID = train["ID"]

N = 250000000000
#Feature 1 - Composition (adding 20 columns)
f1_data = {}
f1f2_data = {}
f1f2f5_data = {}


composition = []
residues = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
for row in train.itertuples(): 
    sequence = row[3]
    length = len(sequence)
    if(length>2*N):
        sequence = sequence[:N] + sequence[-1*N:]
        length = 2*N
    temp = []
    for residue in residues:
        count = 0
        for i in sequence:
            if(i == residue):
                count+=1
        presence = count/length
        temp.append(presence*100)
    composition.append(temp)


for i in range(20):
    res = residues[i]
    temp = []
    for j in composition:
        temp.append(j[i])
    train[res] = temp
    f1_data[res] = temp
    f1f2_data[res] = temp
    f1f2f5_data[res] = temp
    
f1 = pd.DataFrame(data=f1_data)

# print("done1")

# Feature 2 - Dipeptide composition (adding 400 columns)
f2_data = {}

all_pairs = []
for i in range(20):
    res1 = residues[i]
    for j in range(20):
        res2 = residues[j]
        all_pairs.append(res1+"-"+res2)
        
dipep_data = []
for row in train.itertuples():
    sequence = row[3]
    length = len(sequence)
    if(length>2*N):
        sequence = sequence[:N] + sequence[-1*N:]
        length = 2*N
    pairs = {}
    for i in range(20):
        res1 = residues[i]
        for j in range(20):
            res2 = residues[j]
            pairs[res1+"-"+res2] = 0
    
    for i in range(length-1):
        pairs[sequence[i]+"-"+sequence[(i)+1]] += 1
    
    values = []
    for key,value in pairs.items():
        values.append((value/(length-2))*100)
    
    dipep_data.append(values)

for pair in range(400):
    temp = []
    for data in dipep_data:
        temp.append(data[pair])
    train[all_pairs[pair]] = temp
    f2_data[all_pairs[pair]] = temp
    f1f2_data[all_pairs[pair]] = temp
    f1f2f5_data[all_pairs[pair]] = temp

f1f2 = pd.DataFrame(data=f1f2_data)
f2 = pd.DataFrame(data=f2_data)

# print("done2")

        

# Feature 5 - mass, charge, pi (adding 3 columns)

f5_data = {}

mass = {"A" : 89.1,"R" : 174.2,"N" : 132.1,"D" : 133.1,"C" : 121.2,"E" :147.1, "Q" : 146.2, "G" : 146.2,"H" : 155.2,"I" : 131.2,"L" : 131.2,"K" : 146.2,"M" : 149.2, "F" : 165.2 , "P": 115.1, "S" : 105.1,"T" : 119.1,"W" : 204.2, "Y" : 181.2, "V" : 117.1}

masses = []
for row in train.itertuples(): 
    sequence = row[3]
    temp = 0
    for i in sequence:
        temp+=mass[i]
    masses.append(temp)

train["mass"] = masses
f5_data["mass"] = masses
f1f2f5_data["mass"] = masses

charges = []
for row in train.itertuples(): 
    sequence = row[3]
    protein = IP(sequence)
    temp = protein.charge_at_pH(7.0)
    charges.append(temp)

train["charge"] = charges
f5_data["charge"] = charges
f1f2f5_data["charge"] = charges
pIs = []
for row in train.itertuples(): 
    sequence = row[3]
    protein = IP(sequence)
    temp = protein.pi()
    pIs.append(temp)

train["pI"] = pIs
f5_data["pI"] = pIs
f1f2f5_data["pI"] = pIs

f5 = pd.DataFrame(data=f5_data)
f1f2f5 = pd.DataFrame(data=f1f2f5_data)

# print(type(f1f2f5))
# print(type(f1f2f5_data))
# print("done5")

# full_data = copy.deepcopy(f1f2f5_data)
# full_data["labels"] = Lable
# full_data_pd= pd.DataFrame(data=full_data)
# full_data_pd.to_csv("data_out.csv", index=False, header=True)





final_train = {}

count = 0
for i in train:
    if(i!="ID" and i!="Lable" and i!="Sequence" ):
        temp = train[i].tolist()
        final_train[i] = temp
        count+=1
           
train_df = pd.DataFrame(data=final_train)

# print(train_df)
# print(train_df.at[1481,"Lable"])


    
# accuracy_f1 = 0
# for i in range(3):
#     X_train, X_test, y_train, y_test = train_test_split(f1,Lable,test_size = 0.3)

#     # print(train_df.shape())
#     # print(Lable.shape())


#     cls = RandomForestClassifier(n_estimators = 250)
#     cls.fit(X_train,y_train)

#     pred = cls.predict(X_test)
#     accuracy_f1+=metrics.accuracy_score(y_test,y_pred=pred)

# accuracy_f2 = 0
# for i in range(3):
#     X_train, X_test, y_train, y_test = train_test_split(f2,Lable,test_size = 0.3)

#     # print(train_df.shape())
#     # print(Lable.shape())


#     cls = RandomForestClassifier(n_estimators = 250)
#     cls.fit(X_train,y_train)

#     pred = cls.predict(X_test)
#     accuracy_f2+=metrics.accuracy_score(y_test,y_pred=pred)
    

# accuracy_f1f2 = 0
# for i in range(3):
#     X_train, X_test, y_train, y_test = train_test_split(f1f2,Lable,test_size = 0.3)

#     # print(train_df.shape())
#     # print(Lable.shape())


#     cls = RandomForestClassifier(n_estimators = 250)
#     cls.fit(X_train,y_train)

#     pred = cls.predict(X_test)
#     accuracy_f1f2+=metrics.accuracy_score(y_test,y_pred=pred)

# accuracy_f1f2f5 = 0
# for i in range(3):
#     X_train, X_test, y_train, y_test = train_test_split(f1f2f5,Lable,test_size = 0.3)

#     # print(train_df.shape())
#     # print(Lable.shape())


#     cls = RandomForestClassifier(n_estimators = 250)
#     cls.fit(X_train,y_train)

#     pred = cls.predict(X_test)
#     accuracy_f1f2f5+=metrics.accuracy_score(y_test,y_pred=pred)
    
# print(accuracy_f1/3,accuracy_f2/3,accuracy_f1f2/3, accuracy_f1f2f5/3)


##########
test_data = {}
ID = test["ID"].tolist()
composition = []
for seq in test.itertuples(): 
    sequence = seq[2]
    length = len(sequence)
    temp = []
    for residue in residues:
        count = 0
        for i in sequence:
            if(i == residue):
                count+=1
        presence = count/length
        temp.append(presence*100)
    composition.append(temp)


for i in range(20):
    res = residues[i]
    temp = []
    for j in composition:
        temp.append(j[i])
    
    test_data[res] = temp

#####################


all_pairs = []
for i in range(20):
    res1 = residues[i]
    for j in range(20):
        res2 = residues[j]
        all_pairs.append(res1+"-"+res2)
        
dipep_data = []
for row in test.itertuples():
    sequence = row[2]
    length = len(sequence)
    pairs = {}
    for i in range(20):
        res1 = residues[i]
        for j in range(20):
            res2 = residues[j]
            pairs[res1+"-"+res2] = 0
    
    for i in range(length-1):
        pairs[sequence[i]+"-"+sequence[i+1]] += 1
    
    values = []
    for key,value in pairs.items():
        values.append((value/(length-1))*100)
    
    dipep_data.append(values)

for pair in range(400):
    temp = []
    for data in dipep_data:
        temp.append(data[pair])
    test_data[all_pairs[pair]] = temp
#####################

    
test = pd.DataFrame(data=test_data)
# print(test)
cls = RandomForestClassifier(n_estimators = 250)
cls.fit(f1f2, Lable)
pred = cls.predict(test)

# print(pred)
output = {}
output["ID"] = ID
output["Label"] = pred

final = pd.DataFrame(data=output)
final.to_csv(r'output.csv', index = False)

###########
# print("done")
# Any results you write to the current directory are saved as output.