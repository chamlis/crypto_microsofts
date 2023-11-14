def pp_table(titles, rows):
    rows = list(rows)
    row_lens = [max((len(str(row[ix])) for row in ([titles] + rows))) for ix in range(len(rows[0]))]

    def print_row(r):
        s = "|"
        for ix, x in enumerate(r):
            s += f" {str(x).rjust(row_lens[ix])} |"
        print(s)

    print_row(titles)

    for row in rows:
        print_row(row)

    print()

def baby_step_giant_step(a, b):
    p = a.modulus()
    assert p == b.modulus()

    l = p-1
    bound = ceil(sqrt(l))

    print(f"Range is {0}--{bound}")
    print()

    x = b**(-sqrt(l))

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
