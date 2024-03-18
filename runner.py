from lib.tools import *
import pandas as pd
S = 930.5
K = 785
T = 100
r = 0.01
vol = 0.23
tree_steps = 1000

c1 = trinomial_tree_option_price(S,K,T/252,r,vol,tree_steps)
#c2 = trinomial_tree_option_price(S,K,T,r,vol,tree_steps)
#c3 = monte_carlo_option_price(S,K,T,r,vol)
#c4 = black_scholes_merton(S,K,T,r,vol)

#call = [c1,c2,c3,c4]

#p1 = binomial_tree_option_price(S,K,T,r,vol,tree_steps,option_type='put')
#p2 = trinomial_tree_option_price(S,K,T,r,vol,tree_steps,option_type='put')
#p3 = monte_carlo_option_price(S,K,T,r,vol,option_type='put')
#p4 = black_scholes_merton(S,K,T,r,vol,option_type='put')

#put = [p1,p2,p3,p4]
print("Data: S0 = "+str(S)+", K = "+str(K)+", r = "+str(r)+", vol = "+str(vol)+", steps = "+str(tree_steps))
#data = {"Method vs Type" : ["Binomial Tree","Trinomial Tree","Monte carlo simulation","Black-Scholes model"],
#        "Call" : call}
#df = pd.DataFrame(data)
print(c1)