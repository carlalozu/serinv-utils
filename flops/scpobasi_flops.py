
def scpobasi_flops(n_diagonals, n_offdiags, arrowhead_blocksize):


    FLOPS = 0
    counts = {
        'triangular_solve_nb3': 0,
        'vector_scaling_ns': 0,
        'vector_scaling_nb': 0,
        'element_scaling': 0,
        'div': 0,
        'dot_product_ns': 0,
        'dot_product_nb': 0,
        'matrix_vector_nsns': 0,
        'matrix_vector_nsnb': 0,
        'matrix_vector_nbnb': 0,
        'DGEMM_nb3': 0,
    }

    # Arrowhead inversion first
    counts['div'] += 1
    FLOPS += 1

    # X_{ndb+1, ndb+1}
    counts['DGEMM_nb3'] += 1
    FLOPS += 2 * arrowhead_blocksize**3

    # X_{ndb+1,ndb}
    counts['matrix_vector_nbnb'] += 1
    FLOPS += 2 * arrowhead_blocksize * arrowhead_blocksize
    counts['vector_scaling_nb'] += 1
    FLOPS += arrowhead_blocksize

    # X_{ndb,ndb}
    counts['dot_product_nb'] += 1
    FLOPS += 2 * arrowhead_blocksize
    counts['element_scaling'] += 1
    FLOPS += 1

    # Rest of the matrix
    for i in range(2, n_diagonals+1):

        # Adjust for the size of the block E, under the diagonal
        tail = min(i - 1, n_offdiags)

        # Inverse of the L diagonal value i, L_{i, i}^{-1}
        counts['div'] += 1
        FLOPS += 1

        # --- Off-diagonal slice part ---
        # X_{i+1, i} = (-X_{i+1, i+1} L_{i+1, i} -
        #              X_{ndb+1, i+1}^{T} L_{ndb+1, i}) L_{i, i}^{-1}
        counts['matrix_vector_nsns'] += 1
        FLOPS += 2 * tail**2
        counts['matrix_vector_nsnb'] += 1
        FLOPS += 2 * tail * arrowhead_blocksize
        counts['vector_scaling_ns'] += 1
        FLOPS += tail

        # --- Arrowhead part ---
        # X_{ndb+1, i} = (- X_{ndb+1, i+1} L_{i+1, i} -
        #                X_{ndb+1, ndb+1} L_{ndb+1, i}) L_{i, i}^{-1}
        counts['matrix_vector_nsnb'] += 1
        FLOPS += 2 * tail * arrowhead_blocksize
        counts['matrix_vector_nbnb'] += 1
        FLOPS += 2 * arrowhead_blocksize**2
        counts['vector_scaling_nb'] += 1
        FLOPS += arrowhead_blocksize

        # --- Diagonal value part ---
        # X_{i, i} = (L_{i, i}^{-T} - X_{i+1, i}^{T} L_{i+1, i} -
        #            X_{ndb+1, i}.conj().T L_{ndb+1, i}) L_{i, i}^{-1}
        counts['dot_product_ns'] += 1
        FLOPS += 2 * tail
        counts['dot_product_nb'] += 1
        FLOPS += 2 * arrowhead_blocksize
        counts['element_scaling'] += 1
        FLOPS += 1

    return FLOPS
