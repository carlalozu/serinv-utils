# Copyright 2023-2024 ETH Zurich. All rights reserved.

def scpobbaf_flops(n_diag_blocks, diagonal_blocksize, arrowhead_blocksize, n_offdiags_blk):

    FLOPS = 0
    counts = {
        'cholesky_ns3': 0,
        'cholesky_nb3': 0,
        'triang_solve_ns3': 0,
        'triang_solve_ns2nb': 0,
        'dgemm_ns3': 0,
        'dgemm_ns2nb': 0,
        'dgemm_nb2ns': 0
    }
    for i in range(n_diag_blocks-1):
        # L_{i, i} = chol(A_{i, i})
        FLOPS += diagonal_blocksize**3/3 + diagonal_blocksize**2 / \
            2 + diagonal_blocksize/6  # cholesky ns
        counts['cholesky_ns3'] += 1

        for j in range(1, min(n_offdiags_blk + 1, n_diag_blocks - i)):
            # L_{i+j, i} = A_{i+j, i} @ L_{i, i}^{-T}
            FLOPS += diagonal_blocksize**3  # triangular solve ns3
            counts['triang_solve_ns3'] += 1

            for k in range(1, j):
                # L_{i+j, i+k} = A_{i+j, i+k} - L_{i+j, i} @ L_{i+k, i}^{T}
                FLOPS += 2*diagonal_blocksize**3  # dgemm ns3
                counts['dgemm_ns3'] += 1

            # Update next diagonal block
            # A_{i+j, i+j} = A_{i+j, i+j} - L_{i+j, i} @ L_{i+j, i}.conj().T
            FLOPS += 2*diagonal_blocksize**3  # dgemm ns3
            counts['dgemm_ns3'] += 1

        # Part of the decomposition for the arrowhead structure
        # L_{ndb+1, i} = A_{ndb+1, i} @ L_{i, i}^{-T}
        FLOPS += diagonal_blocksize**2*arrowhead_blocksize  # triangular solve ns2nb
        counts['triang_solve_ns2nb'] += 1

        for k in range(1, min(n_offdiags_blk + 1, n_diag_blocks - i)):
            # L_{ndb+1, i+k} = A_{ndb+1, i+k} - L_{ndb+1, i} @ L_{i+k, i}^{T}
            FLOPS += 2*diagonal_blocksize**2*arrowhead_blocksize  # dgemm ns2nb
            counts['dgemm_ns2nb'] += 1

        # L_{ndb+1, ndb+1} = A_{ndb+1, ndb+1} - L_{ndb+1, i} @ L_{ndb+1, i}^{T}
        FLOPS += 2*diagonal_blocksize*arrowhead_blocksize**2  # dgemm nsnb2
        counts['dgemm_nb2ns'] += 1

    # L_{ndb, ndb} = chol(A_{ndb, ndb})
    FLOPS += diagonal_blocksize**3/3 + diagonal_blocksize**2 / \
        2 + diagonal_blocksize/6  # cholesky ns
    counts['cholesky_ns3'] += 1

    # L_{ndb+1, nbd} = A_{ndb+1, nbd} @ L_{ndb, ndb}^{-T}
    FLOPS += diagonal_blocksize**2*arrowhead_blocksize  # triangular solve ns2nb
    counts['triang_solve_ns2nb'] += 1

    # A_{ndb+1, ndb+1} = A_{ndb+1, ndb+1} - L_{ndb+1, ndb} @ L_{ndb+1, ndb}^{T}
    FLOPS += 2*diagonal_blocksize*arrowhead_blocksize**2  # dgemm nsnb2
    counts['dgemm_nb2ns'] += 1

    # L_{ndb+1, ndb+1} = chol(A_{ndb+1, ndb+1})
    FLOPS += arrowhead_blocksize**3/3 + arrowhead_blocksize**2 / \
        2 + arrowhead_blocksize/6  # cholesky nb
    counts['cholesky_nb3'] += 1

    return int(FLOPS), counts
