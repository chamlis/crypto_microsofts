def elgamal_sign(m, g, s_k, p, k=None):
    assert m.base_ring().order() == p-1
    assert g.base_ring().order() == p

    if k is None:
        k = FiniteField(p).random_element()

    assert k.base_ring().order() == p

    r = g^k
    k_inv = Integers(p-1)(k)^-1
    sig = Integers(p-1)(k_inv * (ZZ(m) - ZZ(s_k)*ZZ(r)))

    return (r, sig)

def elgamal_verify(m, r, sig, g, p_k, p):
    assert m.base_ring().order() == p-1
    assert r.base_ring().order() == p
    assert sig.base_ring().order() == p-1
    assert g.base_ring().order() == p
    assert p_k.base_ring().order() == p

    g_m = g^m
    other_side = p_k^r * r^sig

    return g_m == other_side
