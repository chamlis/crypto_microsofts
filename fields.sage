def is_finite_field(p, fx):
    r = PolynomialRing(FiniteField(p), 'x')
    x = r.gen()
    modulo = fx(x)
    factored = modulo.factor()
    print(factored)
    return len(factored) == 1 and factored[0][1] == 1

def make_finite_field(p, fx):
    r = PolynomialRing(FiniteField(p), 'x')
    x = r.gen()
    modulo = fx(x)

    s = r.quotient(modulo, 'a')
    a = s.gen()

    return s, a
