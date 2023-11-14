def pp_table(titles, rows, prefix=""):
    rows = list(rows)
    row_lens = [max((len(str(row[ix])) for row in ([titles] + rows))) for ix in range(len(rows[0]))]

    def print_row(r):
        s = f"{prefix}|"
        for ix, x in enumerate(r):
            s += f" {str(x).rjust(row_lens[ix])} |"
        print(s)

    print_row(titles)

    for row in rows:
        print_row(row)

    print()

def pohlig_hellman(a, b, p=None):
    p = a.modulus()
    assert p == b.modulus()

    factors = factor(p-1)

    print(f"Prime factorization of {p-1} is {' * '.join([f'{f[0]}^{f[1]}' for f in factors])}")
    print()

    answer_mod = {}
    for prime, max_power in factors:
        print(f"Computing result (mod {prime**max_power})")

        print(f"  Write result = {' + '.join([f'{prime**p}r_{p}' for p in range(max_power)])} (mod {prime**max_power})")
        print(f"  Where 0 <= r_i < {prime}")
        print()

        pbits = []

        for power in range(1, max_power+1):
            mod = prime**power
            f = (p-1)//mod

            rhs = a**f
            lhs = 1

            print(f"  Computing result (mod {mod})")

            for ix in range(power-1):
                lhs *= b**(f*pbits[ix]*prime**ix)

            print(f"  Find r_{power-1} st. {rhs} = {sum([x*(prime**ix) for ix, x in enumerate(pbits)])} + {prime**(power-1)}*r_{power-1} (mod {mod})")

            vals = []
            for candidate in range(prime):
                final = lhs * b**(f*candidate*prime**(power-1))
                vals.append(final)

                if final == rhs:
                    pbits.append(candidate)
                    break
            else:
                assert False

            pp_table((f"r_{power-1}", f"{sum([x*(prime**ix) for ix, x in enumerate(pbits[:-1])])} + {prime**(power-1)}*r_{power-1}"), enumerate(vals), prefix="  ")
            print(f"  r_{power-1} = {pbits[-1]}")
            print()

        overall = 0
        for ix, x in enumerate(pbits):
            overall += x * prime**ix
        answer_mod[prime**max_power] = overall

        print(f"  result = {' + '.join([f'{prime**p}*{pbits[p]}' for p in range(max_power)])} (mod {prime**max_power})")
        print(f"  result = {overall} (mod {prime**max_power})")
        print()


    print("Combining moduli with SRT")

    final_answer = 0
    current_mod = 1

    for base, modulus in answer_mod.items():
        m = current_mod
        n = base

        y, c, d = xgcd(m, n)
        assert y == 1

        a = final_answer
        b = modulus

        final_answer = b*c*m + a*d*n

        current_mod = m * n

        if m != 1:
            print(f"  1 = {c}*{m} + {d}*{n}")
            print(f"  a={a}, b={b}, c={c}, d={d}, m={m}, n={n}")
            print(f"  result = {b}*{c}*{m} + {a}*{d}*{n} (mod {current_mod})")

            print(f"  result = {final_answer % current_mod} (mod {current_mod})")
            print()

    return final_answer % (p-1)

def baby_step_giant_step(a, b):
    p = a.modulus()
    assert p == b.modulus()

    l = p-1
    bound = ceil(sqrt(l))

    print(f"Range is {0}--{bound}")
    print()

    x = b**(-bound)

    b = [b**i for i in range(bound)]

    pp_table(("i", "b_i"), enumerate(b))

    c = []

    for j in range(bound+1):
        value = a * x**j
        c.append(value)

        if value in b:
            i = b.index(value)
            result = i + (bound * j)

            pp_table(("j", "c_j"), enumerate(c))
            print(f"b_{i} == c_{j}")
            print(f"{i} + ({bound} * {j}) = {result}")

            return result

    return None

def pollard_rho(a, b):
    p = a.modulus()
    assert p == b.modulus()

    values = [(b, 1, 0)]

    while True:
        g_last, b_last, c_last = values[-1]

        mod = ZZ(g_last) % 3

        if mod == 0:
            next_values = (g_last * b, b_last + 1, c_last)
        elif mod == 1:
            next_values = (g_last * a, b_last, c_last + 1)
        elif mod == 2:
            next_values = (g_last^2, 2 * b_last, 2 * c_last)
        else:
            assert False

        values.append(next_values)

        for ix, value in enumerate(values[:-1]):
            if value[0] == next_values[0]:
                field = Integers(p-1)
                b, c = field(value[1]), field(value[2])
                b_prime, c_prime = field(next_values[1]), field(next_values[2])

                result = (b - b_prime)/(c_prime - c)

                pp_table(("G_i", "b_i", "c_i"), values)

                print(f"G_{ix} == G_{len(values)-1}")
                print(f"a = ({b} - {b_prime})/({c_prime} - {c})")
                print(f"a = {b - b_prime}/{c_prime - c}")
                print(f"a = {result}")

                return result
