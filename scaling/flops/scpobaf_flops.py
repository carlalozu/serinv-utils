
def scpobaf_flops(n_diagonals, n_offdiags, arrowhead_blocksize):

    # Process banded part of the matrix
    FLOPS: int = 0

    counts = {
        'vector_scaling_ns': 0,
        'vector_scaling_nb': 0,
        'div': 0,
        'sqrt': 0,
        'dot_product_ns': 0,
        'matrix_vector_nsns': 0,
        'matrix_vector_nsnb': 0,
        'DGEMM_nb3': 0,
        'GER_nb2': 0
    }

    for i in range(n_diagonals-1):

        # L_{i, i} = chol(A_{i, i})
        counts['sqrt'] += 1
        FLOPS += 1  # 1 sqrt

        counts['div'] += 1
        FLOPS += 1  # 1 div

        # Update column i of the lower diagonals
        # matrix times vector (ns-1)*ns
        counts['matrix_vector_nsns'] += 1
        FLOPS += 2 * n_offdiags*(n_offdiags-1)

        # L_{i+1, i} = A_{i+1, i} @ L_{i, i}^{-T}
        counts['vector_scaling_ns'] += 1
        FLOPS += n_offdiags  # vector scaling ns

        # L_{ndb+1, i} = A_{ndb+1, i} @ L_{i, i}^{-T}
        counts['vector_scaling_nb'] += 1
        FLOPS += arrowhead_blocksize  # vector scaling nb

        # A_{ndb+1, i+1} = A_{ndb+1, i+1} - L_{ndb+1, i} @ L_{i+1, i}.conj().T
        b = i+1-max(0, i-n_offdiags+1)
        counts['matrix_vector_nsnb'] += 1
        FLOPS += 2 * arrowhead_blocksize * b

        # A_{ndb+1, ndb+1} = A_{ndb+1, ndb+1} - L_{ndb+1, i} @ L_{ndb+1,
        # i}.conj().T
        counts['GER_nb2'] += 1
        FLOPS += 2 * arrowhead_blocksize**2  # outer product

        # Update next diagonal
        # A_{i+1, i+1} = A_{i+1, i+1} - L_{i+1, i} @ L_{i+1, i}.conj().T
        counts['dot_product_ns'] += 1
        FLOPS += 2 * n_offdiags  # dot product

    # L_{ndb, ndb} = chol(A_{ndb, ndb})
    counts['sqrt'] += 1
    FLOPS += 1  # 1 sqrt

    # L_{ndb+1, ndb} = A_{ndb+1, ndb} @ L_{ndb, ndb}^{-T}
    counts['div'] += 1
    FLOPS += 1  # div
    counts['vector_scaling_nb'] += 1
    FLOPS += arrowhead_blocksize  # vector scaling nb

    # A_{ndb+1, ndb+1} = A_{ndb+1, ndb+1} - L_{ndb+1, ndb} @ L_{ndb+1, ndb}^{T}
    counts['GER_nb2'] += 1
    FLOPS += 2 * arrowhead_blocksize**2  # outer product, one mult one add

    # L_{ndb+1, ndb+1} = chol(A_{ndb+1, ndb+1})
    counts['cholesky_nb'] += 1
    FLOPS += 1/3*arrowhead_blocksize**3 + 1/2 * \
        arrowhead_blocksize**2 + 1/6*arrowhead_blocksize  # cholesky

    return int(FLOPS), counts
