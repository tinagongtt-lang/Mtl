"""
Mtl 数学库 - 扩展工具模块
"""
from decimal import Decimal
from .exp_log import sqrt, sqrt_hp

def sgn(x):
    """Sgn[x] 符号函数：返回 1, -1 或 0"""
    if isinstance(x, Decimal):
        if x > Decimal('0'): return Decimal('1')
        if x < Decimal('0'): return Decimal('-1')
        return Decimal('0')
    if x > 0.0: return 1.0
    if x < 0.0: return -1.0
    return 0.0

def clamp(x, min_val, max_val):
    """Clamp[x, min, max] 将 x 限制在 [min, max] 区间内"""
    if x < min_val: return min_val
    if x > max_val: return max_val
    return x

def modulus(*args):
    """Modulus[...] 标量绝对值或多维向量的欧几里得模长 (L2范数)"""
    if not args: return 0.0
    if len(args) == 1:
        x = args[0]
        return x if x >= 0 else -x
        
    if isinstance(args[0], Decimal):
        sq_sum = sum(i * i for i in args)
        return sqrt_hp(sq_sum)
    else:
        sq_sum = sum(float(i) * float(i) for i in args)
        return sqrt(sq_sum)