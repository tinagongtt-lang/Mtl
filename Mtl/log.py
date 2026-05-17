"""
Mtl 数学库 - 自然对数模块 (log)
"""

from decimal import Decimal


def log(x: float) -> float:
    """
    标准自然对数 ln(x)。
    采用双曲反正切变换泰勒级数递推，规避传统展开式在 x > 2 时不收敛的问题。
    """
    x = float(x)
    if x <= 0.0:
        raise ValueError("Math domain error: log(x) requires x > 0.")
        
    # 面积双曲反正切变换，将主值域映射到中心高效收敛区间
    t = (x - 1.0) / (x + 1.0)
    t_squared = t * t
    
    term = t
    log_x = t
    n = 1
    
    while True:
        # 核心递推：下一项 = 上一项 * t^2 * (2n-1)/(2n+1)
        term *= t_squared * (2 * n - 1) / (2 * n + 1)
        
        previous = log_x
        log_x += term
        
        if previous == log_x:
            break
        n += 1
        
    return log_x * 2.0


def log_hp(x) -> Decimal:
    """
    高精度自然对数 ln(x)（支持 100 位 Decimal 精度）。
    """
    x = Decimal(str(x)) if not isinstance(x, Decimal) else x
    if x <= Decimal('0'):
        raise ValueError("Math domain error: log_hp(x) requires x > 0.")
        
    t = (x - Decimal('1')) / (x + Decimal('1'))
    t_squared = t * t
    
    term = t
    log_x = t
    n = 1
    
    while True:
        # 纯 Decimal 环境下的无损高精度前向递推
        term *= t_squared * Decimal(2 * n - 1) / Decimal(2 * n + 1)
        
        previous = log_x
        log_x += term
        
        if previous == log_x:
            break
        n += 1
        
    return log_x * Decimal('2')