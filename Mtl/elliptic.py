"""
Mtl 数学库 - 椭圆积分模块 (Elliptic Integrals)
"""
from .exp_log import sqrt
from .pi import *

def elliptic_k(m: float) -> float:
    """EllipticK[m] 第一类完全椭圆积分 (AGM算法)"""
    m = float(m)
    if m >= 1.0 or m < 0.0: raise ValueError("Elliptic K requires 0 <= m < 1")
    
    a0 = 1.0
    b0 = sqrt(1.0 - m)
    
    while True:
        a1 = 0.5 * (a0 + b0)
        b1 = sqrt(a0 * b0)
        if abs(a1 - a0) < 1e-15:
            break
        a0, b0 = a1, b1
    return 3.141592653589793 / (2.0 * a1)


def elliptic_e(m: float) -> float:
    """EllipticE[m] 第二类完全椭圆积分 (AGM算法)"""
    m = float(m)
    if m == 1.0: return 1.0
    if m > 1.0 or m < 0.0: raise ValueError("Elliptic E requires 0 <= m <= 1")
    
    a0 = 1.0
    b0 = sqrt(1.0 - m)
    c0 = sqrt(m)
    
    sum_c = 0.5 * (c0 * c0)
    w = 1.0
    
    while True:
        w *= 2.0
        a1 = 0.5 * (a0 + b0)
        b1 = sqrt(a0 * b0)
        c1 = 0.5 * (a0 - b0)
        sum_c += w * (c1 * c1)
        if abs(a1 - a0) < 1e-15:
            break
        a0, b0 = a1, b1
        
    return (3.141592653589793 / (2.0 * a1)) * (1.0 - sum_c)