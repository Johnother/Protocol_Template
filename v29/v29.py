import numpy as np
import pandas as pd

# BB88 protocol experiment without decoys


df = pd.read_csv(
    "data/bb84_1.csv",
    header=None,
    dtype={
        0: "int64",
        1: "string",
        2: "int64",
        3: "string",
        4: "int64",
        5: "string",
        6: "int64",
    },
)

df.columns = [
    "bit", "polA",
    "A", "polB",
    "B", "polE",
    "BI"
]

print(df.head())

bit = df["bit"].tolist()
polA = df["polA"].tolist()
A = df["A"].tolist()
polB = df["polB"].tolist()
B = df["B"].tolist()
polE = df["polE"].tolist()
BI = df["BI"].tolist()

print(bit)



# Print the sifted bits and their number
sifted = []

for i in range(len(bit)):
    if polA[i] == polB[i]:
        sifted.append([bit[i], A[i], B[i]])
    else: continue
    
print("-------------- Sifted bits for BB84 without Eve ----------------------------------")
print(f"{len(sifted)} sifted bits. Ratio: {len(sifted)/52:.4f}")
for i in range(len(sifted)):
    print(sifted[i])

# With Eve
sifted = []

for i in range(len(bit)):
    if polA[i] == polB[i]:
        sifted.append([bit[i], A[i], BI[i]])
    else: continue
    
print("-------------- Sifted bits for BB84 with Eve ----------------------------------")
print(f"{len(sifted)} sifted bits. Ratio: {len(sifted)/52:.4f}")
for i in range(len(sifted)):
    print(sifted[i])

N_err = 0
for i in range(len(sifted)):
    if sifted[i][1] != sifted[i][2]:
        N_err += 1
    else: continue

print(f"{N_err} wrong bits. QBER ratio: {N_err/len(sifted):.4f}")

# BB88 protocol experiment with decoys


df = pd.read_csv(
    "data/bb84_2.csv",
    header=None,
    dtype={
        0: "int64",
        1: "string",
        2: "int64",
        3: "string",
        4: "int64",
        5: "string",
        6: "int64",
    },
)

df.columns = [
    "bit", "polA",
    "A", "polB",
    "B", "polE",
    "BI"
]

print(df.head())

bit = df["bit"].tolist()
polA = df["polA"].tolist()
A = df["A"].tolist()
polB = df["polB"].tolist()
B = df["B"].tolist()
polE = df["polE"].tolist()
BI = df["BI"].tolist()

print(bit)



# Print the sifted bits and their number
sifted = []

for i in range(len(bit)):
    if polA[i] == polB[i]:
        sifted.append([bit[i], A[i], B[i]])
    else: continue
    
print("-------------- Sifted bits for BB84 without Eve ----------------------------------")
print(f"{len(sifted)} sifted bits. Ratio: {len(sifted)/52:.4f}")
for i in range(len(sifted)):
    print(sifted[i])

# With Eve
sifted = []

for i in range(len(bit)):
    if polA[i] == polB[i]:
        sifted.append([bit[i], A[i], BI[i]])
    else: continue
    
print("-------------- Sifted bits for BB84 with Eve ----------------------------------")
print(f"{len(sifted)} sifted bits. Ratio: {len(sifted)/52:.4f}")
for i in range(len(sifted)):
    print(sifted[i])

N_err = 0
N_Decoy = 0
for i in range(len(sifted)):
    if sifted[i][2] == 3:
        N_Decoy += 1
    elif sifted[i][1] != sifted[i][2]:
        N_err += 1
    else: continue

print(f"{N_err} wrong bits. QBER ratio: {N_err/len(sifted):.4f}\t{N_Decoy} decoy states\tDecoy ratio: {N_Decoy/len(sifted):.4f}")
print(f"{N_err+N_Decoy} total errors\tTotal error rate: {(N_err/len(sifted) + N_Decoy/len(sifted)):.4f}")