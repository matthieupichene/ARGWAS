import numpy as np

def slope(x, y):
    x = np.asarray(x)
    y = np.asarray(y)

    n = len(x)
    sum_x = np.sum(x)
    sum_y = np.sum(y)
    sum_xy = np.sum(x * y)
    sum_x2 = np.sum(x * x)

    numerator = n * sum_xy - sum_x * sum_y
    denominator = n * sum_x2 - sum_x ** 2

    return 0.0 if denominator == 0 else numerator / denominator

def r_squared(x, y):
    x = np.asarray(x)
    y = np.asarray(y)

    if x.size == 0 or x.size != y.size:
        return 0.0

    dx = x - np.mean(x)
    dy = y - np.mean(y)

    var_x = np.sum(dx ** 2)
    var_y = np.sum(dy ** 2)

    if var_x == 0 or var_y == 0:
        return 0.0

    cov_xy = np.sum(dx * dy)
    return (cov_xy ** 2) / (var_x * var_y)

def best_subtree_finder(tree, nodelist, ntotal, cut, phenotypes):
    results = [0.0, 0.0, ""]
    best_r_sq = 0.0
    half_ntotal = ntotal // 2

    for i in nodelist:
        subtree = list(tree.leaves(i))
        n = len(subtree)
        if n < cut or n > ntotal - cut:
            continue

        # Create boolean mask
        insubtree = np.zeros(ntotal, dtype=bool)
        insubtree[subtree] = True

        # Convert to diploid genotypes
        genotypes = insubtree[::2].astype(int) + insubtree[1::2].astype(int) - 1

        r_sq = r_squared(genotypes, phenotypes)
        if r_sq > best_r_sq:
            best_r_sq = r_sq
            results[0] = r_sq
            results[1] = slope(genotypes, phenotypes)
            results[2] = ' '.join(map(str,map(int,insubtree)))
            
    return results       
