"""
Mtl 数学库 - 扩展三角与双曲函数模块
"""
from decimal import Decimal
from .sin import sin, sin_hp
from .cos import cos, cos_hp
from .exp_log import exp, exp_hp

def cot(x: float) -> float:
    """Cot[x] 余切"""
    s = sin(x)
    if abs(s) < 1e-15: raise ValueError("Math domain error: cot(x) is undefined.")
    return cos(x) / s

def csc(x: float) -> float:
    """Csc[x] 余割"""
    s = sin(x)
    if abs(s) < 1e-15: raise ValueError("Math domain error: csc(x) is undefined.")
    return 1.0 / s

def sinh(x: float) -> float:
    """Sinh[x] 双曲正弦"""
    return (exp(x) - exp(-x)) / 2.0

def cosh(x: float) -> float:
    """Cosh[x] 双曲余弦"""
    return (exp(x) + exp(-x)) / 2.0

def tanh(x: float) -> float:
    """Tanh[x] 双曲正切"""
    e_pos, e_neg = exp(x), exp(-x)
    return (e_pos - e_neg) / (e_pos + e_neg)

# ================= High Precision =================
def sinh_hp(x) -> Decimal:
    x = Decimal(str(x)) if not isinstance(x, Decimal) else x
    return (exp_hp(x) - exp_hp(-x)) / Decimal('2')

def cosh_hp(x) -> Decimal:
    x = Decimal(str(x)) if not isinstance(x, Decimal) else x
    return (exp_hp(x) + exp_hp(-x)) / Decimal('2')