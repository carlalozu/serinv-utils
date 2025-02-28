def calculate_parameters_banded(m: int, b: int, nb: int):

    parameters = {
        "m": m,
        "b": b,
        "nb": nb,
    }

    try:
        assert b % 2 == 1,  "bandwidth must be odd"

        ns_ = int((b - 1)/2)
        n_ = 1

        effective_b_ = ns_*(n_*2+1)
        parameters["effective_b"] = effective_b_

        assert b <= effective_b_, "Effective bandwidth is smaller than bandwidth"

        nt_ = m - nb
        parameters["nt"] = nt_

        m_ = nt_ + nb

        parameters["ns"] = ns_
        parameters["n"] = n_

        assert m_ == m, "Matrix size must be kept"

        return {
            'parameters': parameters,
            'flag': 1
        }
    except AssertionError as e:
        return {
            'parameters': parameters,
            'flag': 0,
            'error': str(e),
        }


def calculate_parameters_tri_diagonal(m: int, b: int, nb: int):

    parameters = {
        "m": m,
        "b": b,
        "nb": nb,
    }

    try:
        assert b % 2 == 1,  "bandwidth must be odd"

        n_ = 1
        inner_m = m - nb

        ns_ = int((b-1)/2)
        nt_ = inner_m//ns_

        iters = 0
        while inner_m % ns_:
            ns_ += 1
            nt_ = inner_m//ns_
            iters += 1

        parameters["nt"] = nt_
        assert nt_ > 1, "Matrix not divisible, might be prime"

        effective_b_ = ns_*(n_*2)+1
        parameters["effective_b"] = effective_b_
        assert b <= effective_b_, "Effective b is smaller than b"

        parameters["ns"] = ns_
        parameters["n"] = n_

        m_ = nt_*ns_ + nb
        assert m_ == m, "Matrix size must be kept"

        return {
            'parameters': parameters,
            'flag': 1
        }
    except AssertionError as e:
        return {
            'parameters': parameters,
            'flag': 0,
            'error': str(e),
        }


def calculate_parameters_n_diagonal(m: int, b: int, nb: int, n_: int = 1):

    parameters = {
        "m": m,
        "nb": nb,
        "b": b
    }
    try:
        assert b % 2 == 1,  "b must be odd"

        inner_m = m - nb

        ns_ = int((b-1+2*n_-1)/(2*n_))
        nt_ = inner_m//ns_

        iters = 0
        while inner_m % ns_:
            ns_ += 1
            nt_ = inner_m//ns_
            iters += 1

        parameters["nt"] = nt_
        assert nt_ > 1, "Matrix not divisible, might be prime"

        m_ = nt_*ns_ + nb

        parameters["ns"] = ns_
        parameters["n"] = n_

        assert m_ == m, "Matrix size must be kept"

        effective_b_ = ns_*n_*2+1
        parameters["effective_b"] = effective_b_
        assert b <= effective_b_, "Effective b is smaller than b"
        if b <= (ns_*(n_-1)*2+1):
            raise AssertionError('Extra block is being used for b')

        return {
            'parameters': parameters,
            'flag': 1
        }
    except AssertionError as e:
        return {
            'parameters': parameters,
            'flag': 0,
            'error': str(e),
        }
