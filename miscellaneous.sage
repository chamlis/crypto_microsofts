# Example: 2^13 (mod 101) -> square_and_multiply(2, 13, 101)
def square_and_multiply(base, exp, mod):
    pows = exp.digits(2)
    print(f"{base}^{exp} =", " * ".join(reversed(list(map(lambda x: f"{base}^(2^{x[0]})", filter(lambda p: p[1] != 0, enumerate(pows)))))))
    n = Integers(mod)(1)
    for i in reversed(pows):
        if i == 0:
            _n = n^2
            print(f"n = {n}^2 = {_n} (mod {mod})")
            n = _n
        if i == 1:
            _n = n^2 * base
            print(f"n = {n}^2 * {base} = {_n} (mod {mod})")
            n = _n
    return n
