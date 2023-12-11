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

def pohlig_hellman(a, b):
    p = b.modulus()
    group = Integers(p)
    order = b.multiplicative_order()

    factors = factor(order)

    print(f"Prime factorization of {order-1} is {' * '.join([f'{f[0]}^{f[1]}' for f in factors])}")
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
            f = (order)//mod

            rhs = group(a**f)
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

    return final_answer % order

def baby_step_giant_step(a, b, order=None):
    if order is None:
        order = b.multiplicative_order()

    bound = ceil(sqrt(order))

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
    order = b.multiplicative_order()

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
                field = Integers(order)
                b, c = field(value[1]), field(value[2])
                b_prime, c_prime = field(next_values[1]), field(next_values[2])

                denominator = c_prime - c

                if denominator == 0:
                    return None

                if gcd(denominator, order) != 1:
                    # TODO: proper cycle detector
                    return None

                result = (b - b_prime)/denominator

                pp_table(("G_i", "b_i", "c_i"), values)

                print(f"G_{ix} == G_{len(values)-1}")
                print(f"a = ({b} - {b_prime})/({c_prime} - {c})")
                print(f"a = {b - b_prime}/{c_prime - c}")
                print(f"a = {result}")

                return result

def index_calculus(a, b, factor_base):
    p = a.modulus()
    assert p == b.modulus()
    fp = FiniteField(p)
    fp1 = Integers(p-1)

    eqns = []
    ans = []

    sol = None
    for ix in range(p):
        factors = factor(ZZ(b^ix))

        if len(factors) == 0 or not all(map(lambda x: x[0] in factor_base, factors)):
            continue

        """if gcd(p-1, ix) != 1:
            continue"""
        """try:
            i = fp1(ix).inverse()
        except:
            continue"""

        eqn = [0 for _ in range(len(factor_base))]# + [fp1(ix)]
        for prime, power in factors:
            eqn[factor_base.index(prime)] = fp1(power)

        eqns.append(eqn)
        ans.append((fp1(ix),))

        if len(eqns) >= len(factor_base):
            m = Matrix(eqns)
            try:
                sol = m.solve_right(Matrix(ans))
                print(m)
                print(ans)
                break
            except:
                continue

    if sol is None:
        return None

    factor_base_logs = [sol[ix][0] for ix in range(len(factor_base))]

    print(factor_base_logs)

    # TODO: don't do this
    correct_factor_base_logs = [discrete_log(f, b) for f in factor_base]

    if factor_base_logs != correct_factor_base_logs:
        print("Wrong factor base logs :(")
        print(f"Should be {correct_factor_base_logs}")

    factor_base_logs = correct_factor_base_logs

    for ix in range(1, p+1):
        factors = factor(ZZ(b^ix * a))

        if len(factors) == 0 or not all(map(lambda x: x[0] in factor_base, factors)):
            continue

        result = 0
        for prime, power in factors:
            result += power * factor_base_logs[factor_base.index(prime)]

        result -= ix

        result = ZZ(result) % (p-1)

        return result


def _test_dlp():
    algs = (
        pohlig_hellman,
        baby_step_giant_step,
        pollard_rho,
        lambda a, b: index_calculus(a, b, (2,3,5))
    )

    for ix in range(100):
        p = random_prime(100000)
        fp = FiniteField(p)
        g = fp.multiplicative_generator()
        a = fp.random_element()


        try:
            correct = discrete_log(a, g)
            results = [alg(a, g) for alg in algs]
        except Exception as e:
            print(f"Error running {g} {a} {p}")
            raise e

        assert all(map(lambda r: r is None or r == correct, results))
