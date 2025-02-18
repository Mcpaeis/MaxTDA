import numpy as np
from DTM_filtrations import *
from kde_functions import *


def dtm_diagram(mat, m=0.1):
    simplextree_DTM = DTMFiltration(X=mat, m=m, p=2, dimension_max=2)  # creating a simplex tree
    diagram_DTM = simplextree_DTM.persistence()    
    
    # Convert to DataFrame format
    diagram_df = pd.DataFrame([{'dimension': point[0],'birth': point[1][0],'death': point[1][1]} for point in diagram_DTM])
    
    return diagram_DTM, diagram_df, simplextree_DTM

def dtm_bootstrap(base_mat, ref_mat, m, dim, B, kk=0, th=0, is_smooth=False):
    """
    The base_mat and ref_mat are different for the case where there is smoothing.
    The base_mat is the smoothed mat from the orginal ref_mat
    The ref_mat is the unsmoothed mat and should be the one bootstraped
    """
    _, dgm, splx_tree = dtm_diagram(base_mat, m)
    base_persistence = splx_tree.persistence_intervals_in_dimension(dim)
    #gd.bottleneck_distance(I0, I1)
    n_points = base_mat.shape[0]
    btlnck_distances = []
    for b in range(B):
        # Bootstrap sample
        indices = np.random.choice(n_points, size=n_points, replace=True)
        boot_sample = ref_mat[indices]
        # If we're dealing with smooth samples, then apply the kde sampling
        if is_smooth:            
            _, props, _ = kde_fit(boot_sample, k=kk, eps=0)
            min_dens = props[2]
            # Now use this min density to loop over the thresholds
            bt_kde_fit, bt_props, bt_bbox = kde_fit(ref_mat, k=kk, eps=min_dens-th)
            boot_sample, _= kde_sample(kde=bt_kde_fit, proposals=bt_props, bounding_box=bt_bbox, num_samples=n_points)  
        # Compute diagram for bootstrap sample
        _, _, boot_splx_tree = dtm_diagram(boot_sample, m)
        boot_persistence = boot_splx_tree.persistence_intervals_in_dimension(dim)
        # Add to the bottleneck distance
        btlnck_distances.append(gudhi.bottleneck_distance(base_persistence, boot_persistence))

        if (b + 1) % 100 == 0:
            print(f"Completed {b + 1}/{B} bootstrap iterations")
    return btlnck_distances, dgm

# Grid search for the optimal m for each dataspace:

def get_top_2_persistent(diagram, homology_dimension):
    # Filter diagram for given dimension
    dim_diagram = diagram[diagram['dimension'] == homology_dimension]
    if len(dim_diagram) == 0:
        # Return empty array with same number of columns as input
        return np.zeros((2, len(diagram.columns)))
    # Calculate persistence
    persistence = dim_diagram['death'] - dim_diagram['birth']
    # Sort by persistence and get top 2
    top_2 = dim_diagram.iloc[np.argsort(-persistence)[:2]]
    # Convert to numpy array
    top_2_array = top_2.to_numpy()
    # If only one feature exists, pad with zeros
    if len(top_2_array) < 2:
        padding = np.zeros((2 - len(top_2_array), top_2_array.shape[1]))
        top_2_array = np.vstack([top_2_array, padding])
    return top_2_array

def grid_search_dtm(planet, planetspot1, planetspot2, smoothed1, smoothed2, spot, m_grid, dim=1):
    """Perform grid search over DTM parameter m."""
    # Initialize DataFrame for persistence results
    pers_grid = pd.DataFrame({'m': pd.Series(dtype='float64'),'h1_1': pd.Series(dtype='float64'),'grp': pd.Series(dtype='str')})
    
    # Grid search over m
    for m in m_grid:
        # Compute diagrams for different matrices
        _,dgm_planet,_ = dtm_diagram(planet, m)
        _,dgm_planetspot1,_ = dtm_diagram(planetspot1, m)
        _,dgm_planetspot2,_ = dtm_diagram(planetspot2, m)
        _,dgm_smooth1,_ = dtm_diagram(smoothed1, m)
        _,dgm_smooth2,_ = dtm_diagram(smoothed2, m)
        _,dgm_spot,_ = dtm_diagram(spot, m)
        
        # Extract top most persistent H_1 features
        top_2_H1_planet = get_top_2_persistent(dgm_planet, dim)
        top_2_H1_planetspot1 = get_top_2_persistent(dgm_planetspot1, dim)
        top_2_H1_planetspot2 = get_top_2_persistent(dgm_planetspot2, dim)
        top_2_H1_smooth1 = get_top_2_persistent(dgm_smooth1, dim)
        top_2_H1_smooth2 = get_top_2_persistent(dgm_smooth2, dim)
        top_2_H1_spot = get_top_2_persistent(dgm_spot, dim)
        
        # Create new rows for each group
        new_rows = pd.DataFrame([
            {'m': m, 'h1_1': top_2_H1_planet[0, 2] - top_2_H1_planet[0, 1], 'grp': 'Planet'},
            {'m': m, 'h1_1': top_2_H1_planetspot1[0, 2] - top_2_H1_planetspot1[0, 1], 'grp': 'Planet+Spot(1)'},
            {'m': m, 'h1_1': top_2_H1_planetspot2[0, 2] - top_2_H1_planetspot2[0, 1], 'grp': 'Planet+Spot(2)'},
            {'m': m, 'h1_1': top_2_H1_spot[0, 2] - top_2_H1_spot[0, 1], 'grp': 'Spot'},
            {'m': m, 'h1_1': top_2_H1_smooth1[0, 2] - top_2_H1_smooth1[0, 1], 'grp': 'Smooth(1)'},
            {'m': m, 'h1_1': top_2_H1_smooth2[0, 2] - top_2_H1_smooth2[0, 1], 'grp': 'Smooth(2)'}
        ])
        
        # Append to results
        pers_grid = pd.concat([pers_grid, new_rows], ignore_index=True)
    return pers_grid