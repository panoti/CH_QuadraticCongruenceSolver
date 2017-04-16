#!/usr/bin/env python


def pow_mod(base, exponent, modulus):
    """
    Computing Modular exponentiation: base ^ exponent (mod modulus)

    Parameters
    ----------
    base : integer
    exponent : integer
    modulus : unsigned integer

    Returns
    -------
    int
        Result of calculation

    """
    return pow(base, exponent, modulus)


def inverse(a, p):
    """
    Computing multiplicative inverses in modular structures

    Parameters
    ----------
    a : integer
    p : unsigned integer

    Returns
    -------
    int
        Inverse of a (mod p)

    """
    t, newt = 0, 1
    r, newr = p, a

    while newr != 0:
        quotient = r // newr
        t, newt = newt, t - quotient * newt
        r, newr = newr, r - quotient * newr

    if r > 1:
        return None  # Vi p nguyen to nen truong hop nay khong xay ra

    if t < 0:
        t = t + p

    return t


def legendre(a, p):
    """
    Legendre function (a/p): (a/p) = a^((p-1)/2) (mod p)

    Parameters
    ----------
    a : integer
    p : unsigned integer

    Returns
    -------
    int
        Result of Legendre function

    """
    return pow_mod(a, (p - 1) // 2, p)


def tonelli_shanks(alpha, p):
    """
    Solve a congruence equation with Tonelli-Shanks algorithm:
        x^2 = alpha (mod p)
    With alpha \in Z_p, p is a prime number

    Parameters
    ----------
    alpha : integer
    p : unsigned integer

    Returns
    -------
    int
        One beta solution, the other can calculate by p - beta
    None
        No solution

    """
    if legendre(alpha, p) != 1:
        return None  # Vo nghiem

    q = p - 1
    s = 0

    while q % 2 == 0:
        q //= 2
        s += 1

    if s == 1:
        return pow_mod(alpha, (p + 1) // 4, p)

    for z in range(2, p):
        if p - 1 == legendre(z, p):
            break

    c = pow_mod(z, q, p)
    r = pow_mod(alpha, (q + 1) // 2, p)
    t = pow_mod(alpha, q, p)
    m = s
    t2 = 0

    while (t - 1) % p != 0:
        t2 = (t * t) % p

        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break

            t2 = (t2 * t2) % p

        b = pow_mod(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r


def solve_congruence(a, b, p):
    """
    Solve a congruence equation:
        a*x = b (mod p)
    With a, b \in Z_p, p is a prime number

    Parameters
    ----------
    a : integer
    b : integer
    p : unsigned integer

    Returns
    -------
    int
        Single unique solution
    []
        Infinitely many solutions
    None
        No solution

    """
    if a == 0:
        if b == 0:
            return []  # Vo so nghiem
        else:
            return None  # Vo nghiem
    else:
        return b * inverse(a, p) % p


def solve_quadratic_congruence(a, b, c, p):
    """
    Solve a quadratic congruence equation:
        a*x^2 + b*x + c = 0 (mod p)
    With a, b, c \in Z_p, p is a prime number less than 4 billion
    ...

    Parameters
    ----------
    a : integer
    b : integer
    c : integer
    p : unsigned integer

    Returns
    -------
    [int, int]
        Two solution of equation
    []
        Infinitely many solutions
    None
        No solution

    """
    if a == 0:
        return solve_congruence(b, -c, p)
    else:
        a_inv = inverse(a, p)
        ba = (b * a_inv) % p
        ca = (c * a_inv) % p
        b_div_2 = (ba * inverse(2, p)) % p
        alpha = (pow_mod(b_div_2, 2, p) - ca) % p
        y = tonelli_shanks(alpha, p)

        if y is None:
            return None  # Vo nghiem

        x1 = (y - b_div_2) % p
        x2 = (p - y - b_div_2) % p

        return [x1, x2]


if __name__ == '__main__':
    ttest = [(1, 1, -9, 11), (1, 6, 5, 7), (1, 6, 11, 31),
             (53212, 42124, 53321, 104395303)]

    for a, b, c, p in ttest:
        print(solve_quadratic_congruence(a, b, c, p))
