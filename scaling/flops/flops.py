# from const import OPERATIONS_FLOPS, ALG_OPERATION_COUNT


def T_flops_POBAF(nt: int, ns: int, nb: int):
    return int(0
               + (1/3) * nb**3 + (1/2) * nb**2 + (1/6) * nb  # cholesky nb
               + (1)*(nt)           # sqrt
               + (1)*(nt)           # div
               + (ns)*(nt-1)        # vector scaling ns
               + (nb)*(nt)          # vector scaling nb
               + (2*ns)*(nt-1)      # dot product
               + (2*(ns-1)*ns)*(nt-1)  # matrix product, ns-1*ns
               + (2*nb*ns)*(nt-1)   # matrix product, nb*ns
               + (2*nb**2)*(nt)     # outer product
               )


def T_flops_POBASI(nt: int, ns: int, nb: int):
    return int(0
               + (nb**3)*(1)        # triangular solve nb3
               + (ns)*(nt-1)        # vector scaling ns
               + (nb)*(nt)          # vector scaling nb
               + (1)*(nt)           # mult
               + (1)*(nt)           # div
               + (2*ns)*(nt-1)      # dot product ns
               + (2*nb)*(nt)        # dot product nb
               + (2*ns**2)*(nt-1)   # matrix_vector_nsns
               + (2*nb*ns)*2*(nt-1)  # matrix_vector_nsnb
               + (2*nb**2)*(nt)     # matrix_vector_nbnb
               + (2*nb**3)*(1)      # DGEMM nb3
               )


def T_flops_theory(nt: int, ns: int, nb: int):
    return int(0
               + (1/3) * nb**3 + (1/2) * nb**2 + (1/6) * nb  # cholesky nb
               + (1)*(nt)           # sqrt
               + (1)*(nt)           # div
               + (ns)*(nt-1)        # vector scaling ns
               + (nb)*(nt)          # vector scaling nb
               + (2*ns)*(nt-1)      # dot product
               + ((ns-1)*ns)*(nt-1)  # matrix product, ns-1*ns HALF
               + (2*nb*ns)*(nt-1)   # matrix product, nb*ns
               + (2*nb**2)*(nt)     # outer product
               )


def T_flops_POBTAF(nt: int, ns: int, nb: int):
    return int(0
               + (1/3 * nb**3 + 1/2 * nb**2 + 1/6 * nb)*(1)  # cholesky nb
               + (1/3 * ns**3 + 1/2 * ns**2 + 1/6 * ns)*(nt)  # cholesky ns
               + (ns**3)*(nt-1)        # triangular solve ns3
               + (ns**2*nb)*(nt)       # triangular solve ns2nb
               + (2*ns**3)*(nt-1)      # DGEMM ns3
               + (2*ns**2*nb)*(nt-1)   # DGEMM ns2nb
               + (2*ns*nb**2)*(nt)     # DGEMM nsnb2
               )


def T_flops_POBTASI(nt: int, ns: int, nb: int):
    return int(0
               + (ns**3)*(nt)           # triangular solve ns3
               + (nb**3)*(1)            # triangular solve ns2
               + (2*ns**3)*(4*(nt-1)+1)  # DGEMM ns3
               + (2*ns**2*nb)*(4*(nt-1)+2)  # DGEMM ns2nb
               + (2*ns*nb**2)*(nt)      # DGEMM nsnb2
               + (2*nb**3)*(1)          # DGEMM nb3
               )


def T_flops_POBBAF(nt, ns, nb, n):
    return int(0
               + (1/3 * ns**3 + 1/2 * ns**2 + 1/6 * ns)*(nt)  # cholesky ns
               + (1/3 * nb**3 + 1/2 * nb**2 + 1/6 * nb)*(1)  # cholesky nb
               + (ns**3)*(nt-1)*n       # triangular solve ns3
               + (ns**2*nb)*nt          # triangular solve ns2nb
               + (2*ns**3)*(nt-1)*n**2  # DGEMM ns3
               + (2*ns**2*nb)*(nt-1)*n  # DGEMM ns2nb
               + (2*ns*nb**2)*(nt)      # DGEMM nsnb2
               )


def T_flops_POBBASI(nt: int, ns: int, nb: int, n):
    return int(0
               + (ns**3)*(nt)           # triangular solve ns3
               + (nb**3)*(1)            # triangular solve nb3
               + (2*ns**3)*((nt-1)*(n+1)**2+1)  # DGEMM ns3
               + (2*ns**2*nb)*(2*(nt-1)*(n+1)+2)  # DGEMM ns2nb
               + (2*ns*nb**2)*(nt)      # DGEMM nsnb2
               + (2*nb**3)*(1)          # DGEMM nb3
               )


def T_FLOPS_unpack(alg: str, nt: int, ns: int, nb: int, n: int):
    flops: int = 0
    for op in ALG_OPERATION_COUNT[alg]:
        flops += OPERATIONS_FLOPS[op]['flops'](ns, nb) * \
            ALG_OPERATION_COUNT[alg][op](nt, n)
    return int(flops)
