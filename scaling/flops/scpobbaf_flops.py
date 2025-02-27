# Copyright 2023-2024 ETH Zurich. All rights reserved.

def scpobbaf_flops(n_diag_blocks, diagonal_blocksize, arrowhead_blocksize, n_offdiags_blk):

    FLOPS = 0
    counts = {
        'cholesky_ns': 0,
        'cholesky_nb': 0,
        'triangular_solve_ns3': 0,
        'triangular_solve_ns2nb': 0,
        'DGEMM_ns3': 0,
        'DGEMM_ns2nb': 0,
        'DGEMM_nsnb2': 0
    }
    for i in range(n_diag_blocks-1):
        # L_{i, i} = chol(A_{i, i})
        FLOPS += diagonal_blocksize**3/3 + diagonal_blocksize**2 / \
            2 + diagonal_blocksize/6  # cholesky ns
        counts['cholesky_ns'] += 1

        for j in range(1, min(n_offdiags_blk + 1, n_diag_blocks - i)):
            # L_{i+j, i} = A_{i+j, i} @ L_{i, i}^{-T}
            FLOPS += diagonal_blocksize**3  # triangular solve ns3
            counts['triangular_solve_ns3'] += 1

            for k in range(1, j):
                # L_{i+j, i+k} = A_{i+j, i+k} - L_{i+j, i} @ L_{i+k, i}^{T}
                FLOPS += 2*diagonal_blocksize**3  # DGEMM ns3
                counts['DGEMM_ns3'] += 1

            # Update next diagonal block
            # A_{i+j, i+j} = A_{i+j, i+j} - L_{i+j, i} @ L_{i+j, i}.conj().T
            FLOPS += 2*diagonal_blocksize**3  # DGEMM ns3
            counts['DGEMM_ns3'] += 1

        # Part of the decomposition for the arrowhead structure
        # L_{ndb+1, i} = A_{ndb+1, i} @ L_{i, i}^{-T}
        FLOPS += diagonal_blocksize**2*arrowhead_blocksize  # triangular solve ns2nb
        counts['triangular_solve_ns2nb'] += 1

        for k in range(1, min(n_offdiags_blk + 1, n_diag_blocks - i)):
            # L_{ndb+1, i+k} = A_{ndb+1, i+k} - L_{ndb+1, i} @ L_{i+k, i}^{T}
            FLOPS += 2*diagonal_blocksize**2*arrowhead_blocksize  # DGEMM ns2nb
            counts['DGEMM_ns2nb'] += 1

        # L_{ndb+1, ndb+1} = A_{ndb+1, ndb+1} - L_{ndb+1, i} @ L_{ndb+1, i}^{T}
        FLOPS += 2*diagonal_blocksize*arrowhead_blocksize**2  # DGEMM nsnb2
        counts['DGEMM_nsnb2'] += 1

    # L_{ndb, ndb} = chol(A_{ndb, ndb})
    FLOPS += diagonal_blocksize**3/3 + diagonal_blocksize**2 / \
        2 + diagonal_blocksize/6  # cholesky ns
    counts['cholesky_ns'] += 1

    # L_{ndb+1, nbd} = A_{ndb+1, nbd} @ L_{ndb, ndb}^{-T}
    FLOPS += diagonal_blocksize**2*arrowhead_blocksize  # triangular solve ns2nb
    counts['triangular_solve_ns2nb'] += 1

    # A_{ndb+1, ndb+1} = A_{ndb+1, ndb+1} - L_{ndb+1, ndb} @ L_{ndb+1, ndb}^{T}
    FLOPS += 2*diagonal_blocksize*arrowhead_blocksize**2  # DGEMM nsnb2
    counts['DGEMM_nsnb2'] += 1

    # L_{ndb+1, ndb+1} = chol(A_{ndb+1, ndb+1})
    FLOPS += arrowhead_blocksize**3/3 + arrowhead_blocksize**2 / \
        2 + arrowhead_blocksize/6  # cholesky nb
    counts['cholesky_nb'] += 1

    return int(FLOPS), counts
