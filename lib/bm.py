import numpy as np

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
