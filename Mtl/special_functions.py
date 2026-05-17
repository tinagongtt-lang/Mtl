"""
Mtl 数学库 - 特殊函数模块 (Gamma, Erf, Bessel)
"""
from decimal import Decimal
from .exp_log import sqrt, sqrt_hp, log, log_hp, exp, exp_hp
from .pi import *

_SQRT_2PI = 2.5066282746310005  # √(2π)

def gamma(x: float) -> float:
    """Gamma[x] 伽马函数 (Lanczos 近似法，g=7, n=9)"""
    x = float(x)
    if x < 0.5:
        # 利用反射公式：Γ(x) * Γ(1-x) = π / sin(πx)
        from .sin import sin
        pi = 3.141592653589793
        return pi / (sin(pi * x) * gamma(1.0 - x))
    
    # Lanczos 系数
    p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
         771.32342877765313, -176.61502916214059, 12.507343278686905,
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
    
    x -= 1.0
    a = p[0]
    for i in range(1, len(p)):
        a += p[i] / (x + i)
    t = x + 7.5
    return _SQRT_2PI * (t ** (x + 0.5)) * exp(-t) * a


def erf(x: float) -> float:
    """ErrorFunction Erf[x] 误差函数 (泰勒级数展开递推)"""
    x = float(x)
    pi = 3.141592653589793
    term = x
    sum_erf = x
    x_squared = x * x
    n = 1
    while True:
        # 递推式
        term *= -x_squared / n
        next_item = term / (2 * n + 1)
        prev = sum_erf
        sum_erf += next_item
        if prev == sum_erf: break
        n += 1
    return sum_erf * (2.0 / sqrt(pi))


def bessel_j(n_order: int, x: float) -> float:
    """BesselJ[n, x] 第一类 n 阶贝塞尔函数 (仅支持整数阶 0 和 1)"""
    x = float(x)
    if n_order not in (0, 1):
        raise NotImplementedError("Only 0-th and 1-st order Bessel functions are currently supported.")
    
    term = 1.0 if n_order == 0 else x / 2.0
    sum_bessel = term
    x_squared_quarter = (x * x) / 4.0
    k = 1
    while True:
        # 递推式：分母结合了阶乘演进
        term *= -x_squared_quarter / (k * (k + n_order))
        prev = sum_bessel
        sum_bessel += term
        if prev == sum_bessel: break
        k += 1
    return sum_bessel