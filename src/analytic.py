"""Polynomial Taylor bounds; acb/arb arithmetic is used for certificates.

Coefficients are in ascending order. All exact inputs enter as integers or
rational strings; no binary floats enter the validated path.
"""
from math import comb, factorial
from flint import acb, arb


def mul(a, b):
    out = [acb(0) for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def derivative(a, order=1):
    for _ in range(order):
        a = [k*a[k] for k in range(1, len(a))]
    return a


def evaluate(a, z):
    out = acb(0)
    for x in reversed(a):
        out = out*z+x
    return out


def shift(a, z):
    return [sum((a[k]*comb(k, j)*z**(k-j) for k in range(j, len(a))), acb(0))
            for j in range(len(a))]


def family(name, epsilon):
    q = [acb(arb('2/5'), arb('1/5')), acb(-3), acb(0), acb(1)]
    t = acb(epsilon)/acb(3, 1)
    if name == 'aligned4':
        g = [acb(comb(4, k))*(-t)**k for k in range(5)]
    elif name == 'opposed2':
        g = [acb(1), acb(0), -2*t*t, acb(0), t**4]
    elif name == 'ring4':
        g = [acb(1), acb(0), acb(0), acb(0), -t**4]
    elif name == 'broken_ring4':
        g = mul([t**k for k in range(4)], [acb(1), -t/arb('1001/1000')])
    else:
        raise ValueError(name)
    return q, g, mul(q, g)


def stationary_certificate(a, z, radius):
    """Return a strict Rouche witness and stationary-value enclosure.

    f0, f1, f2 are actual derivatives at z, M3 bounds |f'''| on the
    complete closed complex disk. Every boolean uses a strict arb comparison.
    """
    r = arb(radius)
    assert r > 0
    t = shift(a, z)
    f0, f1, f2 = t[0], t[1], 2*t[2]
    m3 = sum((k*(k-1)*(k-2)*abs(t[k])*r**(k-3)
              for k in range(3, len(t))), arb(0))
    margin = abs(f2)*r-abs(f1)-m3*r*r/2
    if not (abs(f2) > 0 and margin > 0):
        return None
    v = f0-f1*f1/(2*f2)
    err = m3*r**3/6+m3*m3*r**4/(8*abs(f2))
    return dict(center=v, error=err, radius=r, margin=margin, m3=m3,
                f0=f0, f1=f1, f2=f2)


def find_certificate(a, z):
    """Search a fixed, rational radius list; no approximate-root assumption."""
    radii = ['1/1000000000000','1/300000000000','1/100000000000',
             '1/30000000000','1/10000000000','1/3000000000','1/1000000000',
             '1/300000000','1/100000000','1/30000000','1/10000000',
             '1/3000000','1/1000000','1/300000','1/100000','1/30000',
             '1/10000','1/3000','1/1000','1/300','1/100','1/30','1/10',
             '1/5','1/3','1/2']
    for radius in radii:
        cert = stationary_certificate(a, z, radius)
        if cert is not None:
            cert['radius_exact'] = radius
            return cert
    return None


def mp_coeffs(a, mp):
    # Conversion is exclusively for independent high-precision checks.
    return [mp.mpc(x.real.mid().str(110, radius=False),
                   x.imag.mid().str(110, radius=False)) for x in a]
