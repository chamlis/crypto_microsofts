# This is for questions where the nonce gets repeated

# Example:
# You ask Alice to sign two messages m1 = 12 and m2 = 23 using
# the ElGamal signature scheme. She sends you her public key
# (p, g, pkA) = (31, 3, 15) and two signatures (r1, sig1) = (24, 6)
# and (r2, sig2) = (24, 23)
#   -> elgamal_get_secret_key(31, 3, 15, 24, 12, 6, 23, 23)
def elgamal_get_secret_key(p, g, pk, r, m1, sig1, m2, sig2):
    mod = p-1
    R = Integers(mod)
    k = (R(m1) - R(m2))/(R(sig1) - R(sig2))
    print (f"nonce = ({m1}-{m2})/({sig1}-{sig2}) = {k} (mod {mod})")
    num = (m1 - sig1 * k)
    print(f"""
sk = (m1 - sig1 * k)/r (mod {mod})
   = ({m1} - {sig1} * {k})/{r} (mod {mod})
   = {num}/{r} (mod {mod})
""")

    cd = gcd(mod, r)
    _mod = mod/cd
    _R = Integers(_mod)
    a = _R(num) / (_R(r//cd))
    print(f"sk = {num}/{r/cd} = {a} (mod {_mod})")

    sk_pos = [int(a) + mult*int(_mod) for mult in range(cd)]
    print(f"candidate values of sk = {sk_pos} (mod {mod})")

    for sk in sk_pos:
        if (Integers(p)(g)^sk) == pk:
            print(f"sk = {sk} (mod {mod})")
            return sk

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

def check_rsa_public(e, n):
    factors = n.factor()
    coef_sum = sum([coef for _,coef in factors])
    assert len(factors) == 2 and coef_sum == 2, f"n is not of the form p*q\nn = {n} = {factors}"
    [(p,_),(q,_)] = factors
    print(f"n = {n}")
    print(f"p*q = {p}*{q}")
    phi = euler_phi(n)
    print(f"phi(n) = {phi}")
    assert gcd(phi, e) == 1, f"phi(n) = {phi} is not coprime to e = {e}"
    print(f"({e},{n}) is a valid key pair!")
