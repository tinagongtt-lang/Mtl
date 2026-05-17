"""
Mtl 数学库 - 高等数学与特殊函数综合模块 (advanced_math)
整合了正交多项式、解析数论、信号激活、矩阵核及特殊积分函数
"""

import time
from decimal import Decimal

# =====================================================================
# 0. 内部核心基础函数与常量注入 (确保当前文件完全闭环且无需外部交叉引用)
# =====================================================================

try:
    from .pi import calculate_pi_chudnovsky
    _PI = float(calculate_pi_chudnovsky(50))
except Exception:
    _PI = 3.14159265358979323846

_TWO_PI = _PI * 2.0
_HALF_PI = _PI / 2.0
_SQRT_2PI = 2.5066282746310005  # √(2π)

def _local_sin(x: float) -> float:
    """内部极速递推 sin 驱动"""
    x = float(x) % _TWO_PI
    if x > _PI: x -= _TWO_PI
    elif x < -_PI: x += _TWO_PI
    term = x
    sin_x = x
    x_squared = x * x
    n = 1
    while True:
        term *= -x_squared / ((2 * n) * (2 * n + 1))
        prev = sin_x
        sin_x += term
        if prev == sin_x: break
        n += 1
    return sin_x

def _local_cos(x: float) -> float:
    """内部极速递推 cos 驱动"""
    return _local_sin(float(x) + _HALF_PI)

def _local_exp(x: float) -> float:
    """内部极速递推 exp 驱动"""
    x = float(x)
    term, exp_x, n = 1.0, 1.0, 1
    while True:
        term *= x / n
        prev = exp_x
        exp_x += term
        if prev == exp_x: break
        n += 1
    return exp_x

def _local_log(x: float) -> float:
    """内部极速递推 log (ln) 驱动"""
    x = float(x)
    if x <= 0.0: raise ValueError("Math domain error")
    t = (x - 1.0) / (x + 1.0)
    t_squared = t * t
    term, log_x, n = t, t, 1
    while True:
        term *= t_squared * (2 * n - 1) / (2 * n + 1)
        prev = log_x
        log_x += term
        if prev == log_x: break
        n += 1
    return log_x * 2.0

def _local_sqrt(x: float) -> float:
    """内部开方驱动"""
    x = float(x)
    if x < 0.0: raise ValueError("Math domain error")
    if x == 0.0: return 0.0
    g = x / 2.0
    while True:
        prev = g
        g = 0.5 * (g + x / g)
        if abs(prev - g) < 1e-15: break
    return g


# =====================================================================
# 1. 正交多项式 (Polynomials)
# =====================================================================

def hermite_h(n: int, x: float) -> float:
    """埃尔米特多项式 Hermite Polynomials H_n(x)"""
    if n < 0: raise ValueError("n 必须为非负整数")
    if n == 0: return 1.0
    p0, p1 = 1.0, 2.0 * float(x)
    for i in range(2, n + 1):
        p2 = 2.0 * x * p1 - 2.0 * (i - 1) * p0
        p0, p1 = p1, p2
    return p1

def chebyshev_t(n: int, x: float) -> float:
    """第一类切比雪夫多项式 Chebyshev Polynomials T_n(x)"""
    if n < 0: raise ValueError("n 必须为非负整数")
    if n == 0: return 1.0
    p0, p1 = 1.0, float(x)
    for i in range(2, n + 1):
        p2 = 2.0 * x * p1 - p0
        p0, p1 = p1, p2
    return p1

def legendre_p(n: int, x: float) -> float:
    """勒让德多项式 Legendre Polynomials P_n(x)"""
    if n < 0: raise ValueError("n 必须为非负整数")
    if n == 0: return 1.0
    p0, p1 = 1.0, float(x)
    for i in range(2, n + 1):
        p2 = ((2.0 * i - 1.0) * x * p1 - (i - 1) * p0) / i
        p0, p1 = p1, p2
    return p1

def laguerre_l(n: int, x: float) -> float:
    """拉盖尔多项式 Laguerre Polynomials L_n(x)"""
    if n < 0: raise ValueError("n 必须为非负整数")
    if n == 0: return 1.0
    p0, p1 = 1.0, 1.0 - float(x)
    for i in range(2, n + 1):
        p2 = ((2.0 * i - 1.0 - x) * p1 - (i - 1) * p0) / i
        p0, p1 = p1, p2
    return p1


# =====================================================================
# 2. 数论与离散数学函数 (Number Theory)
# =====================================================================

def zeta(s: float, terms: int = 1000) -> float:
    """黎曼Zeta函数 Riemann Zeta Function ζ(s) [暂实现s>1的狄利克雷级数收敛]"""
    if s <= 1.0: raise NotImplementedError("s <= 1 时需要解析延拓暂未支持")
    return sum(1.0 / (n ** s) for n in range(1, terms + 1))

def moebius(n: int) -> int:
    """默比乌斯函数 Möbius Function μ(n)"""
    if n <= 0: raise ValueError("n 必须为正整数")
    if n == 1: return 1
    p, d = 0, 2
    while d * d <= n:
        if n % d == 0:
            p += 1
            n //= d
            if n % d == 0: return 0
        d += 1
    if n > 1: p += 1
    return -1 if p % 2 != 0 else 1

def totient(n: int) -> int:
    """欧拉Phi函数 Euler's Totient Function φ(n)"""
    if n <= 0: raise ValueError("n 必须为正整数")
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0: n //= p
            result -= result // p
        p += 1
    if n > 1: result -= result // n
    return result

def divisor_count(n: int) -> int:
    """因子个数函数 Divisor Function d(n)"""
    if n <= 0: return 0
    count, d = 0, 1
    while d * d <= n:
        if n % d == 0:
            count += 2 if d * d != n else 1
        d += 1
    return count

def divisor_sum(n: int) -> int:
    """因子和函数 Sigma Function σ(n)"""
    if n <= 0: return 0
    total, d = 0, 1
    while d * d <= n:
        if n % d == 0:
            total += d
            if d * d != n: total += n // d
        d += 1
    return total

def prime_pi(n: int) -> int:
    """黎曼质数计数函数 Prime-counting Function π(n)"""
    if n < 2: return 0
    prime = [True] * (n + 1)
    p = 2
    while p * p <= n:
        if prime[p]:
            for i in range(p * p, n + 1, p): prime[i] = False
        p += 1
    return sum(1 for x in range(2, n + 1) if prime[x])


# =====================================================================
# 3. 统计、概率、基础工具与激活函数 (Statistics & Utilities)
# =====================================================================

def sgn(x: float) -> float:
    """Sgn[x] 符号函数"""
    if x > 0.0: return 1.0
    if x < 0.0: return -1.0
    return 0.0

def clamp(x: float, min_val: float, max_val: float) -> float:
    """Clamp[x, min, max] 区间斩断锁死"""
    if x < min_val: return min_val
    if x > max_val: return max_val
    return x

def modulus(*args) -> float:
    """Modulus[...] 标量绝对值或多维向量模长 (L2 范数)"""
    if not args: return 0.0
    return _local_sqrt(sum(float(i) * float(i) for i in args))

def kronecker_delta(i: int, j: int) -> int:
    """克罗内克 Kronecker Delta 函数 δ_{ij}"""
    return 1 if i == j else 0

def heaviside_step(x: float) -> float:
    """亥维赛阶跃函数 Heaviside Step Function H(x)"""
    if x < 0.0: return 0.0
    if x == 0.0: return 0.5
    return 1.0

def relu(x: float) -> float:
    """线性整流激活函数 ReLU"""
    return max(0.0, x)

def leaky_relu(x: float, alpha: float = 0.01) -> float:
    """渗漏线性整流激活函数 Leaky ReLU"""
    return x if x > 0.0 else alpha * x

def softmax(vector: list) -> list:
    """归一化指数函数 Softmax"""
    max_v = max(vector)
    exps = [_local_exp(x - max_v) for x in vector]
    sum_exps = sum(exps)
    return [e / sum_exps for e in exps]


# =====================================================================
# 4. 高级特殊与超几何积分函数 (Special Functions)
# =====================================================================

def lambert_w(x: float) -> float:
    """兰伯特W函数 Lambert W Function W_0(x)"""
    if x < -0.36787944117144233: raise ValueError("Math domain error")
    if x == 0.0: return 0.0
    w = _local_log(x) if x > 1.0 else x
    for _ in range(100):
        ew = _local_exp(w)
        prev = w
        w -= (w * ew - x) / (ew * (w + 1.0) - (w + 2.0) * (w * ew - x) / (2.0 * w + 2.0))
        if abs(prev - w) < 1e-15: break
    return w

def sinc(x: float) -> float:
    """辛格函数 Sinc Function Sinc[x] = sin(x)/x"""
    if abs(x) < 1e-9: return 1.0
    return _local_sin(x) / x

def erf_val(x: float) -> float:
    """误差函数 Error Function Erf[x]"""
    term = x
    sum_erf = x
    x_squared = x * x
    n = 1
    while True:
        term *= -x_squared / n
        next_item = term / (2 * n + 1)
        prev = sum_erf
        sum_erf += next_item
        if prev == sum_erf: break
        n += 1
    return sum_erf * (2.0 / _local_sqrt(_PI))

def erfc(x: float) -> float:
    """互补误差函数 Complementary Error Function Erfc[x]"""
    return 1.0 - erf_val(x)

def beta(x: float, y: float) -> float:
    """贝塔函数 Beta Function B(x, y) (通过内部封装的极速 Lanczos 伽马转换)"""
    def _internal_gamma(v: float) -> float:
        if v < 0.5: return _PI / (_local_sin(_PI * v) * _internal_gamma(1.0 - v))
        p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
             771.32342877765313, -176.61502916214059, 12.507343278686905,
             -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
        v -= 1.0
        a = p[0]
        for i in range(1, len(p)): a += p[i] / (v + i)
        t = v + 7.5
        return _SQRT_2PI * (t ** (v + 0.5)) * _local_exp(-t) * a

    return (_internal_gamma(x) * _internal_gamma(y)) / _internal_gamma(x + y)

def exp_integral_ei(x: float, terms: int = 200) -> float:
    """指数积分 Exponential Integral Ei(x)"""
    if x <= 0.0: raise ValueError("Ei(x) 目前仅支持 x > 0")
    gamma_euler = 0.5772156649015328
    total = gamma_euler + _local_log(x)
    term = x
    for n in range(1, terms + 1):
        if n > 1: term *= (x * (n - 1)) / (n * n)
        else: term = x
        prev = total
        total += term
        if prev == total: break
    return total

def hypergeometric_2f1(a: float, b: float, c: float, x: float) -> float:
    """高斯超几何函数 Hypergeometric Functions 2F1"""
    if abs(x) >= 1.0: raise NotImplementedError("|x| >= 1 暂未实现解析延拓")
    term, total, n = 1.0, 1.0, 0
    while True:
        term *= ((a + n) * (b + n) / (c + n)) * x / (n + 1)
        prev = total
        total += term
        if prev == total: break
        n += 1
    return total


# =====================================================================
# 5. 矩阵核心与变换核算子 (Matrix & Kernels)
# =====================================================================

def matrix_trace(mat: list) -> float:
    """矩阵迹函数 Trace Function"""
    return sum(float(mat[i][i]) for i in range(len(mat)))

def matrix_det(mat: list) -> float:
    """矩阵行列式函数 Determinant Function (高斯消元法)"""
    n = len(mat)
    a = [[float(x) for x in row] for row in mat]
    det = 1.0
    for i in range(n):
        pivot = i
        for j in range(i + 1, n):
            if abs(a[j][i]) > abs(a[pivot][i]): pivot = j
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            det *= -1.0
        if abs(a[i][i]) < 1e-15: return 0.0
        det *= a[i][i]
        for j in range(i + 1, n):
            factor = a[j][i] / a[i][i]
            for k in range(i, n): a[j][k] -= factor * a[i][k]
    return det

def fourier_transform_kernel(w: float, t: float) -> complex:
    """傅里叶变换核 Fourier Transform Kernel e^{-iwt}"""
    return complex(_local_cos(w * t), -_local_sin(w * t))

def laplace_transform_kernel(s: float, t: float) -> float:
    """拉普拉斯变换核 Laplace Transform Kernel e^{-st}"""
    return _local_exp(-s * t)
# advanced_math.py 的最后一行
__all__ = [
    'hermite_h', 'chebyshev_t', 'legendre_p', 'laguerre_l',
    'zeta', 'moebius', 'totient', 'divisor_count', 'divisor_sum', 'prime_pi',
    'sgn', 'clamp', 'modulus', 'kronecker_delta', 'heaviside_step', 'relu', 'leaky_relu', 'softmax',
    'lambert_w', 'sinc', 'erf_val', 'erfc', 'beta', 'exp_integral_ei', 'hypergeometric_2f1',
    'matrix_trace', 'matrix_det', 'fourier_transform_kernel', 'laplace_transform_kernel'
]