"""
Mtl 数学库 - 余弦模块 (cos)
"""

from decimal import Decimal
from .pi import *
from .sin import sin, sin_hp  # 假设你的 sin 函数在同级目录的 sin.py 中

# 初始化常量
_PI_DECIMAL = Decimal(calculate_pi_chudnovsky(100))
_HALF_PI_DECIMAL = _PI_DECIMAL / 2
_HALF_PI_FLOAT = float(_HALF_PI_DECIMAL)


def cos(x: float) -> float:
    """
    标准余弦函数。
    利用 cos(x) = sin(x + π/2) 完美复用 sin 的前向递推优化。
    """
    return sin(float(x) + _HALF_PI_FLOAT)


def cos_hp(x) -> Decimal:
    """
    高精度余弦函数（支持 100 位 Decimal 精度）。
    """
    x = Decimal(str(x)) if not isinstance(x, Decimal) else x
    return sin_hp(x + _HALF_PI_DECIMAL)