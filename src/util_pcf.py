from contextlib import contextmanager,redirect_stderr,redirect_stdout
from os import devnull


@contextmanager
def silence():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)

def argmin(d):
    return min(d, key=d.get)


def argmax(d):
    return max(d, key=d.get)


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def warning(msg):
    return '{}{}{}'.format('\033[93m', msg, '\033[0m')


def info(msg):
    return '{}{}{}'.format('\033[96m', msg, '\033[0m')


def debug(msg):
    return '{}{}{}'.format('\033[92m', msg, '\033[0m')

def safediv(v):
    return v if v > 0 else 1.0


def forward(f, i, g, limit):
    if limit is None:
        left = g
    else:
        left = min(g, limit)
    right = float(max(f[i] - f[i + 1], 0.0) / safediv(f[i]))
    return left - right


def estimate_forward(f, i, g, knw, limit):
    left = max(g, float(max(knw[i] - knw[i + 1], 0.0) / safediv(knw[i])))
    if limit is not None:
        left = min(left, limit)
    right = float(max(f[i] - f[i + 1], 0.0) / safediv(f[i]))
    return left - right


def central(f, i, g, limit):
    if limit is None:
        left = float(max(f[i - 1] - f[i], 0.0) / safediv(f[i - 1]))
    else:
        left = min(
            limit,
            float(max(f[i - 1] - f[i], 0.0) / safediv(f[i - 1])),
        )
    right = float(max(f[i] - f[i + 1], 0.0) / safediv(f[i]))
    return left - right


def backward(f, i, g, limit):
    if limit is None:
        left = float(max(f[i - 1] - f[i], 0.0) / safediv(f[i - 1]))
    else:
        left = min(
            limit,
            float(max(f[i - 1] - f[i], 0.0) / safediv(f[i - 1])),
        )
    right = g
    return left - right

