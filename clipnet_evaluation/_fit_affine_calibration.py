"""
This script fits an affine transformation between log10(obs) and log10(pred) on
training data and applies it to test data.
"""

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression


# Fit the transformation from source to target
def fit_affine(source_points, target_points):
    # Add a column of 1s to the source points for the intercept
    n_samples, n_features = source_points.shape
    X = np.hstack([source_points, np.ones((n_samples, 1))])  # (N, D+1)

    # Solve linear regression for each dimension
    reg = LinearRegression(fit_intercept=False).fit(X, target_points)
    affine_matrix = reg.coef_  # Shape: (D_target, D_source + 1)
    return affine_matrix


# Apply the affine transformation
def apply_affine(points, affine_matrix):
    n_samples = points.shape[0]
    X = np.hstack([points, np.ones((n_samples, 1))])  # Add bias term
    return X @ affine_matrix.T  # (N, D_target)


def total_least_squares(X, y):
    """
    Calculates the Total Least Squares solution for the equation Xb = y.

    Args:
      X (numpy.ndarray): Independent variable data, shape (n, m).
      Y (numpy.ndarray): Dependent variable data, shape (n, 1).

    Returns:
      numpy.ndarray: Solution matrix B, shape (m, 1).
    """
    # Combine X and y into a single matrix
    Z = np.hstack((X, y.reshape(-1, 1)))
    # Calculate the singular value decomposition of Z
    U, S, V = np.linalg.svd(Z)
    # The TLS solution is the last column of V, normalized
    tls_solution = V[-1, :-1] / V[-1, -1]
    return tls_solution


def total_least_squares(X: np.ndarray, y: np.ndarray):
    """
    Fits a Total Least Squares model to the data (X, y).

    Args:
        X (np.ndarray): A (n_samples, n_features) input matrix.
        y (np.ndarray): A (n_samples,) output vector.

    Returns:
        coef (np.ndarray): Estimated TLS coefficients of shape (n_features,)
        intercept (float): Estimated intercept
    """
    # Combine X and y into one matrix
    A = np.hstack([X, y.reshape(-1, 1)])
    # Perform SVD
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    V = Vt.T
    # Last right singular vector corresponds to smallest singular value
    v = V[:, -1]
    # Separate into coefficients and intercept
    coef = -v[:-1] / v[-1]
    # Estimate intercept
    intercept = -np.mean(X @ coef - y)
    return coef, intercept


X_ = np.expand_dims(
    np.log10(
        pd.read_csv("merged_procap_0.csv.gz", header=None, index_col=0)
        .to_numpy()[:, np.r_[250:750, 1250:1750]]
        .sum(axis=1)
    ),
    axis=-1,
)
y_ = np.log10(np.load("merged_clipnet_0.npz")["arr_1"])

X = []
y = []

for i in range(1, 10):
    X.append(
        np.expand_dims(
            np.log10(
                pd.read_csv(f"merged_procap_{i}.csv.gz", header=None, index_col=0)
                .to_numpy()[:, np.r_[250:750, 1250:1750]]
                .sum(axis=1)
            ),
            axis=-1,
        )
    )
    y.append(np.log10(np.load(f"merged_clipnet_{i}.npz")["arr_1"]))

X = np.concatenate(X)
y = np.concatenate(y)
mat = fit_affine(y, X)

LinearRegression(fit_intercept=True).fit(X_, apply_affine(y_, mat)).coef_
