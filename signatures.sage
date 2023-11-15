def elgamal_sign(m, g, s_k, p, k=None):
    assert m.modulus() == p-1
    assert g.modulus() == p

    if k is None:
        k = FiniteField(p).random_element()

    assert k.modulus() == p

    r = g^k
    k_inv = Integers(p-1)(k)^-1
    sig = Integers(p-1)(k_inv * (ZZ(m) - ZZ(s_k)*ZZ(r)))

    return (r, sig)

def elgamal_verify(m, r, sig, g, p_k, p):
    assert m.modulus() == p-1
    assert r.modulus() == p
    assert sig.modulus() == p-1
    assert g.modulus() == p
    assert p_k.modulus() == p

    g_m = g^m
    other_side = p_k^r * r^sig

    return g_m == other_side

def rsa_sign(m, d, n):
    return (m^d) % n

def rsa_verify(m, sig, e, n):
    return (m % n) == ((sig^e) % n)
