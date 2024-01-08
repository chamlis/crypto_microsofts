from itertools import compress

def print_through(f, v):
    if callable(f):
        f = f(v)

    print(f)
    return v

def iterate(f, x):
    while True:
        yield x
        x = f(x)

def shift_and_reduce(identity, shift_f, reduce_f, shift_op, reduce_op):
    def f(a, b):
        return reduce(
            reduce_f,
            print_through(
                lambda x: f" {reduce_op} ".join(map(str, x)),
                list(
                    compress(
                        map(
                            lambda x: x[1],
                            iterate(
                                lambda x: print_through(
                                    f"{a}{shift_op}{x[0]} = {x[1]}",
                                    (x[0]*2, shift_f(x[1]))
                                ),
                                (1, a)
                            )
                        ),
                        print_through(f"{b} in binary is {bin(b)}", b.digits(2))
                    )
                )
            ),
            identity
        )

    return f

double_and_add = shift_and_reduce(0, lambda a: a+a, lambda a, b: a+b, "*", "+")
square_and_multiply = shift_and_reduce(1, lambda a: a*a, lambda a, b: a*b, "^", "*")
