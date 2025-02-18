import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from matplotlib.colors import LogNorm 
from sklearn.neighbors import KernelDensity


def kde_fit(data, k=2, eps=0.1):
    m = data.shape[1] # the data dimension
    # Fit and extract the minimum density
    bw = knn_bandwidth(data, k)
    kde = KernelDensity(kernel="gaussian", bandwidth=bw).fit(data)
    min_density = np.min(np.exp(kde.score_samples(data)))-eps
    # Define the bounding box
    bounding_box = np.zeros((m, 2))
    for i in range(m):
        min_i, max_i = data[:, i].min() - 1, data[:, i].max() + 1
        bounding_box[i, :] = [min_i, max_i]
    # Constants and prop density
    const_prop_density = 1 / np.prod(bounding_box[:, 1] - bounding_box[:, 0])
    # Define mesh grid
    if m > 3:
        grid_points = np.array([
            np.random.uniform(bounding_box[i, 0], bounding_box[i, 1], 100000) 
            for i in range(m)
        ]).T  # Shape (100000, m) --- this avoids memory issues
    else:
        grid = np.meshgrid(*[np.linspace(bounding_box[i, 0], bounding_box[i, 1], 100) for i in range(m)])
        grid_points = np.vstack([g.ravel() for g in grid]).T
    # Extract the KDE values
    kde_values = kde.score_samples(grid_points)
    scaling_constant = np.max(np.exp(kde_values))
    proposals = np.asarray([const_prop_density, scaling_constant, min_density, bw])
    
    return kde, proposals, bounding_box


def kde_sample(kde, proposals, bounding_box, num_samples):
    m = bounding_box.shape[0]  # number of dimensions
    const_prop_dens, scaling_constant, min_density, _ = proposals
    
    if m > 3:
        grid_points = np.array([
            np.random.uniform(bounding_box[i, 0], bounding_box[i, 1], 10000) 
            for i in range(m)
        ]).T
    else:
        grid = np.meshgrid(*[np.linspace(bounding_box[i, 0], bounding_box[i, 1], 100) for i in range(m)])
        grid_points = np.vstack([g.ravel() for g in grid]).T

    samples = []
    total_samples = 0
    accepted_samples = 0
    while len(samples) < num_samples:
        # Sample a point uniformly within the bounding box in all dimensions
        point = np.array([np.random.uniform(bounding_box[i, 0], bounding_box[i, 1]) for i in range(m)])
        # Evaluate KDE density at the sampled point
        dst = np.exp(kde.score_samples(point.reshape(1, -1)))[0]
        # Sample a uniform random value for rejection
        u = np.random.uniform(0, scaling_constant * const_prop_dens)
        # Accept the sample if it meets the acceptance criteria
        if u <= dst and dst >= min_density:
            samples.append(point)
            accepted_samples += 1
        total_samples += 1
    return np.array(samples), grid_points

def knn_bandwidth(data, k):    
    # Here since we are querying the data with a tree built from it, the 1st nn will be zeros.
    tree = cKDTree(data)
    # The k=k+1 corrects for the above statement, and the function is formulated properly
    bandwidths, _ = tree.query(data, k=k+1, p=2)
    return np.mean(bandwidths[:, k])/2

def plot_2dkde(kde, data, ax="", title = "KDE"):
    x_grid, y_grid = np.mgrid[data[:, 0].min():data[:, 0].max():100j, 
                          data[:, 1].min():data[:, 1].max():100j] 
    pos = np.vstack([x_grid.ravel(), y_grid.ravel()])
    den = np.exp(kde.score_samples(pos.T))
    den = np.reshape(den, x_grid.shape)
    if ax=="":
        _, ax = plt.subplots()
    ax.imshow(np.rot90(den), cmap='inferno', extent=[data[:, 0].min(), data[:, 0].max(),
                                                       data[:, 1].min(), data[:, 1].max()],
                                                       norm=LogNorm())
    ax.scatter(data[:, 0], data[:, 1], alpha=0.5, s=8)
    ax.set_title(title)
    plt.show()
