import numpy as np

def trinomial_tree_option_price(S, K, T, r, sigma, n, option_type='call'):
    delta_t = T / n
    u = np.exp(sigma * np.sqrt(2 * delta_t))
    d = 1 / u
    m = 1
    p_u = ((np.exp(r * delta_t/2) - np.exp(-sigma * np.sqrt(delta_t/2))) /
           (np.exp(sigma * np.sqrt(delta_t/2)) - np.exp(-sigma * np.sqrt(delta_t/2))))**2
    p_d = ((np.exp(sigma * np.sqrt(delta_t/2)) - np.exp(r * delta_t/2)) /
           (np.exp(sigma * np.sqrt(delta_t/2)) - np.exp(-sigma * np.sqrt(delta_t/2))))**2
    p_m = 1 - p_u - p_d

    stock_prices = np.zeros((2*n+1, n+1))
    option_values = np.zeros((2*n+1, n+1))

    for j in range(n+1):
        for i in range(-j, j+1, 2):
            stock_prices[i+n, j] = S * (u ** (j - abs(i)))

    if option_type == 'call':
        option_values[:, -1] = np.maximum(0, stock_prices[:, -1] - K)
    elif option_type == 'put':
        option_values[:, -1] = np.maximum(0, K - stock_prices[:, -1])

    for j in range(n-1, -1, -1):
        for i in range(-j, j+1, 2):
            option_values[i+n, j] = np.exp(-r * delta_t) * (p_u * option_values[i+n+1, j+1] +
                                                            p_d * option_values[i+n-1, j+1] +
                                                            p_m * option_values[i+n, j+1])

    option_price = option_values[n, 0]
    return option_price

S = 100
K = 105
T = 1
r = 0.05
sigma = 0.2
n = 2

call_option_price = trinomial_tree_option_price(S, K, T, r, sigma, n, option_type='call')
put_option_price = trinomial_tree_option_price(S, K, T, r, sigma, n, option_type='put')

print("European Call Option Price:", call_option_price)
print("European Put Option Price:", put_option_price)
