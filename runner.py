from lib.tools import *
import pandas as pd
S = 100
K = 105
T = 1
r = 0.05
vol = 0.2
tree_steps = 2

c1 = binomial_tree_option_price(S,K,T,r,vol,tree_steps)
c2 = trinomial_tree_option_price(S,K,T,r,vol,tree_steps)
c3 = monte_carlo_option_price(S,K,T,r,vol)
c4 = black_scholes_merton(S,K,T,r,vol)

call = [c1,c2,c3,c4]

p1 = binomial_tree_option_price(S,K,T,r,vol,tree_steps,option_type='put')
p2 = trinomial_tree_option_price(S,K,T,r,vol,tree_steps,option_type='put')
p3 = monte_carlo_option_price(S,K,T,r,vol,option_type='put')
p4 = black_scholes_merton(S,K,T,r,vol,option_type='put')

put = [p1,p2,p3,p4]

data = {"Method vs Type" : ["Binomial Tree","Trinomial Tree","Monte carlo simulation","Black-Scholes model"],
        "Call" : call,
        "Put" : put}
df = pd.DataFrame(data)
print(df)