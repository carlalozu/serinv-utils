
def scpobaf_flops(n_diagonals, n_offdiags, arrowhead_blocksize):

    # Process banded part of the matrix
    FLOPS:int = 0
    for i in range(n_diagonals-1):

        # L_{i, i} = chol(A_{i, i})
        FLOPS += 1  # 1 sqrt

        FLOPS += 1 # 1 div

        # Update column i of the lower diagonals
        FLOPS += 2 * n_offdiags*(n_offdiags-1)  # matrix product, 1 mult, 1 add

        # L_{i+1, i} = A_{i+1, i} @ L_{i, i}^{-T}
        FLOPS += n_offdiags  # vector scaling ns

        # L_{ndb+1, i} = A_{ndb+1, i} @ L_{i, i}^{-T}
        FLOPS += arrowhead_blocksize  # vector scaling nb

        b =  i+1-max(0, i-n_offdiags+1)
        # FLOPS += 2 * arrowhead_blocksize * n_offdiags  # matrix times vector (ns-1)*ns
        FLOPS += 2 * arrowhead_blocksize * b  # matrix times vector (ns-1)*ns

        # A_{ndb+1, ndb+1} = A_{ndb+1, ndb+1} - L_{ndb+1, i} @ L_{ndb+1, i}.conj().T
        FLOPS += 2 * arrowhead_blocksize**2  # outer product, one mult one add

        # Update next diagonal
        # A_{i+1, i+1} = A_{i+1, i+1} - L_{i+1, i} @ L_{i+1, i}.conj().T
        FLOPS += 2 * n_offdiags  # dot product

    # L_{ndb, ndb} = chol(A_{ndb, ndb})
    FLOPS += 1  # 1 sqrt

    # L_{ndb+1, ndb} = A_{ndb+1, ndb} @ L_{ndb, ndb}^{-T}
    FLOPS += 1 # div
    FLOPS += arrowhead_blocksize  # vector scaling nb

    # A_{ndb+1, ndb+1} = A_{ndb+1, ndb+1} - L_{ndb+1, ndb} @ L_{ndb+1, ndb}^{T}
    FLOPS += 2 * arrowhead_blocksize**2  # outer product, one mult one add

    # L_{ndb+1, ndb+1} = chol(A_{ndb+1, ndb+1})
    FLOPS += 1/3*arrowhead_blocksize**3 + 1/2*arrowhead_blocksize**2 + 1/6*arrowhead_blocksize  # cholesky

    return FLOPS
