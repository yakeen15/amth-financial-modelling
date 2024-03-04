from lib.tools import *
import pandas as pd

S = 930.5
K_l = [785, 800, 850, 900, 975, 1000]
T_l = [1, 5, 10, 20, 60, 100]
r_l = [0.01, 0.02, 0.05, 0.08, 0.1]
vol_l = [0.23, 0.3, 0.45, 0.6]
tree_steps_l = [2, 5, 10, 30, 60]

#K = K_l[0]
T = T_l[0]
r = r_l[0]
vol = vol_l[0]
tree_steps = tree_steps_l[0]

ck = []
for K in K_l:
    ck1 = binomial_tree_option_price(S,K,T,r,vol,tree_steps)
    ck2 = trinomial_tree_option_price(S,K,T,r,vol,tree_steps)
    ck3 = monte_carlo_option_price(S,K,T,r,vol)
    ck4 = black_scholes_merton(S,K,T,r,vol)
    ck.append([ck1,ck2,ck3,ck4])

K = K_l[0]
#T = T_l[0]
r = r_l[0]
vol = vol_l[0]
tree_steps = tree_steps_l[0]

ct = []

for T in T_l:
    ct1 = binomial_tree_option_price(S,K,T,r,vol,tree_steps)
    ct2 = trinomial_tree_option_price(S,K,T,r,vol,tree_steps)
    ct3 = monte_carlo_option_price(S,K,T,r,vol)
    ct4 = black_scholes_merton(S,K,T,r,vol)
    ct.append([ct1,ct2,ct3,ct4])

K = K_l[0]
T = T_l[0]
#r = r_l[0]
vol = vol_l[0]
tree_steps = tree_steps_l[0]

cr = []

for r in r_l:
    cr1 = binomial_tree_option_price(S,K,T,r,vol,tree_steps)
    cr2 = trinomial_tree_option_price(S,K,T,r,vol,tree_steps)
    cr3 = monte_carlo_option_price(S,K,T,r,vol)
    cr4 = black_scholes_merton(S,K,T,r,vol)
    cr.append([cr1,cr2,cr3,cr4])

K = K_l[0]
T = T_l[0]
r = r_l[0]
#vol = vol_l[0]
tree_steps = tree_steps_l[0]

cv = []

for vol in vol_l:
    cv1 = binomial_tree_option_price(S,K,T,r,vol,tree_steps)
    cv2 = trinomial_tree_option_price(S,K,T,r,vol,tree_steps)
    cv3 = monte_carlo_option_price(S,K,T,r,vol)
    cv4 = black_scholes_merton(S,K,T,r,vol)
    cv.append([cv1,cv2,cv3,cv4])

K = K_l[0]
T = T_l[0]
r = r_l[0]
vol = vol_l[0]
#tree_steps = tree_steps_l[0]

cts = []

for tree_steps in tree_steps_l:
    cts1 = binomial_tree_option_price(S,K,T,r,vol,tree_steps)
    cts2 = trinomial_tree_option_price(S,K,T,r,vol,tree_steps)
    cts3 = monte_carlo_option_price(S,K,T,r,vol)
    cts4 = black_scholes_merton(S,K,T,r,vol)
    cts.append([cts1,cts2,cts3,cts4])


print("Base case: S0 = "+str(S)+", K = "+str(K)+", r = "+str(r)+", vol = "+str(vol)+", steps = "+str(tree_steps))
print("With K ="+str(K_l)+"\n")
data = {"K" : K_l,
        "BM" : [item[0] for item in ck],
        "TM" : [item[1] for item in ck],
        "MC" : [item[2] for item in ck],
        "BS" : [item[3] for item in ck]}
df = pd.DataFrame(data)
df.to_csv('K_variable.csv', index=False)
print("With T ="+str(T_l)+"\n")
data = {"T" : T_l,
        "BM" : [item[0] for item in ct],
        "TM" : [item[1] for item in ct],
        "MC" : [item[2] for item in ct],
        "BS" : [item[3] for item in ct]}
df = pd.DataFrame(data)
df.to_csv('T_variable.csv', index=False)
print("With vol ="+str(vol_l)+"\n")
data = {"vol" : vol_l,
        "BM" : [item[0] for item in cv],
        "TM" : [item[1] for item in cv],
        "MC" : [item[2] for item in cv],
        "BS" : [item[3] for item in cv]}
df = pd.DataFrame(data)
df.to_csv('vol_variable.csv', index=False)
print("With r ="+str(r_l)+"\n")
data = {"r" : r_l,
        "BM" : [item[0] for item in cr],
        "TM" : [item[1] for item in cr],
        "MC" : [item[2] for item in cr],
        "BS" : [item[3] for item in cr]}
df = pd.DataFrame(data)
df.to_csv('r_variable.csv', index=False)
print("With TreeSteps ="+str(tree_steps_l)+"\n")
data = {"TreeSteps" : tree_steps_l,
        "BM" : [item[0] for item in cts],
        "TM" : [item[1] for item in cts],
        "MC" : [item[2] for item in cts],
        "BS" : [item[3] for item in cts]}
df = pd.DataFrame(data)
df.to_csv('TreeSteps_variable.csv', index=False)
print("All data printed")