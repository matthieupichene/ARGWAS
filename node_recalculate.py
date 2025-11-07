def node_recalculate(tree1, tree2):
    nodelist = set()

    # Get all nodes present in either tree
    nodes1 = set(tree1.nodes())
    nodes2 = set(tree2.nodes())
    all_nodes = nodes1.union(nodes2)

    # Build parent maps (with default -1 for missing nodes)
    max_node = max(max(nodes1, default=0), max(nodes2, default=0))
    parent1 = [-1] * (max_node + 1)
    parent2 = [-1] * (max_node + 1)

    for u in nodes1:
        parent1[u] = tree1.parent(u)
    for u in nodes2:
        parent2[u] = tree2.parent(u)

    # Nodes whose parent differs between the trees
    parents_t1 = {parent1[u] for u in all_nodes if parent1[u] != -1 and parent1[u] != parent2[u]}
    parents_t2 = {parent2[u] for u in all_nodes if parent2[u] != -1 and parent2[u] != parent1[u]}

    # Walk up tree2 from t2-exclusive parents
    for i in parents_t2:
        while i != -1:
            nodelist.add(i)
            i = parent2[i]

    # Walk up from t1-exclusive parents until shared ancestor or root
    for i in parents_t1:
        j = i
        while j != -1 and j in parents_t1:
            j = parent1[j]
        while j != -1:
            nodelist.add(j)
            j = parent2[j]

    return list(nodelist)
