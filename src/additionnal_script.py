import numpy as np
from itertools import product
from sklearn.metrics import mean_squared_error

def generate_weights(steps, num_matrices):
    # Generate all possible weights that add up to 1 with the given granularity
    step = 1 / steps
    return [w for w in product(np.arange(0, 1.01, step), repeat=num_matrices) if np.isclose(sum(w), 1)]

def calculate_rmse(Y, Y_hats, weights):
    # Calculate the weighted sum of Y_hats matrices
    weighted_sum = sum(w * Y_hat for w, Y_hat in zip(weights, Y_hats))
    # Calculate RMSE
    return np.sqrt(mean_squared_error(Y, weighted_sum))

def find_optimal_weights(Y, Y_hats, steps=10):
    num_matrices = len(Y_hats)
    weights_list = generate_weights(steps, num_matrices)
    best_rmse = float('inf')
    best_weights = None
    
    # Evaluate all possible combinations of weights
    for weights in weights_list:
        rmse = calculate_rmse(Y, Y_hats, weights)
        if rmse < best_rmse:
            best_rmse = rmse
            best_weights = weights
            
    return best_weights, best_rmse

def calculate_weighted_sum(Y_hats, weights):
    # Compute the weighted sum of the Y_hat matrices using the given weights
    weighted_sum = np.zeros_like(Y_hats[0])
    for Y_hat, weight in zip(Y_hats, weights):
        weighted_sum += weight * Y_hat
    return weighted_sum

def useless_function():
    print("This file only purpose is to show how to add modules or other scripts to your submissions.")