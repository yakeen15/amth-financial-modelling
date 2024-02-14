import numpy as np

def monte_carlo_option_price(S, K, T, r, sigma, option_type='call', num_simulations=100000):
    dt = 1 / 252
    num_steps = int(T / dt)
    S_T = np.zeros(num_simulations)

    for i in range(num_simulations):
        rand = np.random.normal(0, 1, num_steps)
        stock_path = S * np.exp(np.cumsum((r - 0.5 * sigma ** 2) * dt + sigma * np.sqrt(dt) * rand))
        S_T[i] = stock_path[-1]

    if option_type == 'call':
        payoff = np.maximum(0, S_T - K)
    elif option_type == 'put':
        payoff = np.maximum(0, K - S_T)
    else:
        raise ValueError("Invalid option type. Use 'call' or 'put'.")

    option_price = np.exp(-r * T) * np.mean(payoff)
    return option_price

S = 100
K = 105
T = 1
r = 0.05
sigma = 0.2

call_price = monte_carlo_option_price(S, K, T, r, sigma, option_type='call')
put_price = monte_carlo_option_price(S, K, T, r, sigma, option_type='put')

print("Monte Carlo Call option price:", call_price)
print("Monte Carlo Put option price:", put_price)
