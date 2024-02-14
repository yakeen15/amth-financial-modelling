import numpy as np
from scipy.stats import norm

def binomial_tree_option_price(S, K, T, r, sigma, n, option_type='call'):
    delta_t = T / n
    u = np.exp(sigma * np.sqrt(delta_t))
    d = 1 / u
    p = (np.exp(r * delta_t) - d) / (u - d)

    stock_prices = np.zeros((n+1, n+1))
    option_values = np.zeros((n+1, n+1))

    for i in range(n+1):
        for j in range(i+1):
            stock_prices[j, i] = S * (u ** (i - j)) * (d ** j)

    if option_type == 'call':
        option_values[:, -1] = np.maximum(0, stock_prices[:, -1] - K)
    elif option_type == 'put':
        option_values[:, -1] = np.maximum(0, K - stock_prices[:, -1])

    for i in range(n-1, -1, -1):
        for j in range(i+1):
            option_values[j, i] = np.exp(-r * delta_t) * (p * option_values[j, i+1] + (1 - p) * option_values[j+1, i+1])

    option_price = option_values[0, 0]
    return option_price

def trinomial_tree_option_price(S, K, T, r, sigma, n, option_type='call'):
    dt = T / n
    u = np.exp(sigma * np.sqrt(2 * dt))
    d = 1 / u
    p = (np.exp(r * dt / 2) - np.exp(-sigma * np.sqrt(dt / 2))) / (np.exp(sigma * np.sqrt(dt / 2)) - np.exp(-sigma * np.sqrt(dt / 2)))
    
    # Constructing the trinomial tree
    stock_price_tree = np.zeros((n + 1, 2 * n + 1))
    stock_price_tree[0, n] = S
    for i in range(1, n + 1):
        for j in range(-i, i + 1):
            if j == -i:
                stock_price_tree[i, j + n] = stock_price_tree[i - 1, j + 1 + n] * u
            elif j == i:
                stock_price_tree[i, j + n] = stock_price_tree[i - 1, j - 1 + n] * d
            else:
                stock_price_tree[i, j + n] = stock_price_tree[i - 1, j + n] * np.sqrt(u)
    
    option_price_tree = np.zeros((n + 1, 2 * n + 1))
    
    # Calculate option prices at maturity
    if option_type == 'call':
        option_price_tree[n, :] = np.maximum(0, stock_price_tree[n, :] - K)
    else:
        option_price_tree[n, :] = np.maximum(0, K - stock_price_tree[n, :])
    
    # Backward induction to calculate option prices at earlier nodes
    for i in range(n - 1, -1, -1):
        for j in range(-i, i + 1):
            option_price_tree[i, j + n] = np.exp(-r * dt) * (p * option_price_tree[i + 1, j + 1 + n] + (1 - p) * option_price_tree[i + 1, j + n] + p * option_price_tree[i + 1, j - 1 + n])
    
    return option_price_tree[0, n]


def black_scholes_merton(S, K, T, r, sigma, option_type='call'):
    d1 = (np.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)

    if option_type == 'call':
        option_price = S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    elif option_type == 'put':
        option_price = K * np.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)
    else:
        raise ValueError("Invalid option type. Use 'call' or 'put'.")

    return option_price

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