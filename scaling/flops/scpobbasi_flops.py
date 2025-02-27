# Copyright 2023-2024 ETH Zurich. All rights reserved.


def scpobbasi_flops(n_diag_blocks, diagonal_blocksize, arrowhead_blocksize, n_offdiags_blk):
    FLOPS = 0
    counts = {
        'triangular_solve_ns3': 0,
        'triangular_solve_nb3': 0,
        'DGEMM_ns3': 0,
        'DGEMM_ns2nb': 0,
        'DGEMM_nsnb2': 0,
        'DGEMM_nb3': 0
    }

    counts['triangular_solve_nb3'] += 1
    FLOPS += arrowhead_blocksize**3

    counts['DGEMM_nb3'] += 1
    FLOPS += 2 * arrowhead_blocksize**3

    counts['triangular_solve_ns3'] += 1
    FLOPS += diagonal_blocksize**3

    # X_{ndb+1, ndb} = -X_{ndb+1, ndb+1} L_{ndb+1, ndb} L_{ndb, ndb}^{-1}
    counts['DGEMM_nsnb2'] += 1
    FLOPS += 2 * arrowhead_blocksize**2 * diagonal_blocksize
    counts['DGEMM_ns2nb'] += 1
    FLOPS += 2 * arrowhead_blocksize * diagonal_blocksize**2


    # X_{ndb, ndb} = (L_{ndb, ndb}^{-T} - X_{ndb+1, ndb}^{T} L_{ndb+1, ndb}) L_{ndb, ndb}^{-1}
    counts['DGEMM_ns2nb'] += 1
    FLOPS += 2 * diagonal_blocksize**2 * arrowhead_blocksize
    counts['DGEMM_ns3'] += 1
    FLOPS += 2 * diagonal_blocksize**3

    for i in range(n_diag_blocks - 2, -1, -1):

        # L_blk_inv = L_{i, i}^{-1}
        counts['triangular_solve_ns3'] += 1
        FLOPS += diagonal_blocksize**3

        # Arrowhead part
        # X_{ndb+1, i} = - X_{ndb+1, ndb+1} L_{ndb+1, i}
        counts['DGEMM_nsnb2'] += 1
        FLOPS += 2 * arrowhead_blocksize**2 * diagonal_blocksize

        for k in range(i + 1, min(i + n_offdiags_blk + 1, n_diag_blocks), 1):
            # X_{ndb+1, i} = X_{ndb+1, i} - X_{ndb+1, k} L_{k, i}
            counts['DGEMM_ns2nb'] += 1
            FLOPS += 2 * arrowhead_blocksize * diagonal_blocksize**2

        # X_{ndb+1, i} = X_{ndb+1, i} L_{i, i}^{-1}
        counts['DGEMM_ns2nb'] += 1
        FLOPS += 2 * arrowhead_blocksize * diagonal_blocksize**2

        # Off-diagonal block part
        for j in range(min(i + n_offdiags_blk, n_diag_blocks - 1), i, -1):
            # Take the effect of the arrowhead part into account
            # X_{j, i} = - X_{ndb+1, j}.T L_{ndb+1, i}
            counts['DGEMM_ns2nb'] += 1
            FLOPS += 2 * arrowhead_blocksize * diagonal_blocksize**2

            for k in range(i + 1, j):
                # X_{j, i} = X_{j, i} - X_{j, k} @ L_{k, i}, k<j
                counts['DGEMM_ns3'] += 1
                FLOPS += 2 * diagonal_blocksize**3

            # X_{j, i} = X_{j, i} - X_{j, j} @ L_{j, i}, k=j
            counts['DGEMM_ns3'] += 1
            FLOPS += 2 * diagonal_blocksize**3

            for k in range(j+1, min(i + n_offdiags_blk + 1, n_diag_blocks)):
                # X_{j, i} = X_{j, i} - X_{k, j}.T @ L_{k, i}, k>j
                counts['DGEMM_ns3'] += 1
                FLOPS += 2 * diagonal_blocksize**3

            # X_{j, i} = X_{j, i} L_{i, i}^{-1}

            counts['DGEMM_ns3'] += 1
            FLOPS += 2 * diagonal_blocksize**3

        # Diagonal block part
        # X_{i, i} = L_{i, i}^{-T} - X_{ndb+1, i}.T L_{ndb+1, i}
        counts['DGEMM_ns2nb'] += 1
        FLOPS += 2 * diagonal_blocksize**2 * arrowhead_blocksize

        for k in range(i + 1, min(i + n_offdiags_blk + 1, n_diag_blocks), 1):
            # X_{i, i} = X_{i, i} - X_{k, i}.T L_{k, i}
            counts['DGEMM_ns3'] += 1
            FLOPS += 2 * diagonal_blocksize**3

        # X_{i, i} = X_{i, i} L_{i, i}^{-1}
        counts['DGEMM_ns3'] += 1
        FLOPS += 2 * diagonal_blocksize**3

    return int(FLOPS), counts
