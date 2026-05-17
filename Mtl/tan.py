"""
Mtl 数学库 - 正切模块 (tan)
"""

from decimal import Decimal
from .sin import sin, sin_hp
from .cos import cos, cos_hp


def tan(x: float) -> float:
    """
    标准正切函数。
    :raises ValueError: 当 cos(x) 趋近于 0（即处于奇点位置）时抛出。
    """
    c = cos(x)
    # 64位浮点数下的极小值卡位防御
    if abs(c) < 1e-15:
        raise ValueError("Math domain error: tan(x) is undefined at this point (vertical asymptote).")
    return sin(x) / c


def tan_hp(x) -> Decimal:
    """
    高精度正切函数（支持 100 位 Decimal 精度）。
    """
    x = Decimal(str(x)) if not isinstance(x, Decimal) else x
    c = cos_hp(x)
    # 100位高精度环境下的极小值卡位防御
    if abs(c) < Decimal('1e-95'):
        raise ValueError("Math domain error: tan_hp(x) is undefined (vertical asymptote).")
    return sin_hp(x) / c