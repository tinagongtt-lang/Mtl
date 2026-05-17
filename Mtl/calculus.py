"""
Mtl 数学库 - 微积分模块 (calculus)
"""
from decimal import Decimal
from typing import Callable

# ==========================================
# 1. 数值微分 (Derivatives)
# ==========================================

def derivative(f: Callable[[float], float], x: float, h: float = 1e-5) -> float:
    """
    数值导数 f'(x) 
    采用五点中心差分公式 (Five-Point Center Difference)，误差级别为 O(h^4)。
    比普通的 (f(x+h)-f(x))/h 精度高出数个数量级。
    
    :param f: 目标函数，例如 sin, cos 或自定义的 lambda
    :param x: 导数求值点
    :param h: 步长步长
    """
    return (f(x - 2*h) - 8*f(x - h) + 8*f(x + h) - f(x + 2*h)) / (12 * h)


def derivative_hp(f: Callable[[Decimal], Decimal], x: Decimal, h: Decimal = Decimal('1e-15')) -> Decimal:
    """
    高精度数值导数 (支持 100 位 Decimal 精度)。
    """
    if not isinstance(x, Decimal): x = Decimal(str(x))
    if not isinstance(h, Decimal): h = Decimal(str(h))
    
    # 同样采用五点差分
    return (f(x - 2*h) - Decimal('8')*f(x - h) + Decimal('8')*f(x + h) - f(x + 2*h)) / (Decimal('12') * h)


# ==========================================
# 2. 数值积分 (Integrals)
# ==========================================

def integrate(f: Callable[[float], float], a: float, b: float, n: int = 1000) -> float:
    """
    定积分 (Definite Integral) ∫[a, b] f(x) dx
    采用 复合辛普森规则 (Composite Simpson's Rule)，比传统的矩形法、梯形法更精准。
    
    :param f: 被积函数
    :param a: 积分下限
    :param b: 积分上限
    :param n: 切分网格数（必须是偶数，n 越大越精准）
    """
    a, b = float(a), float(b)
    if n % 2 != 0:
        n += 1  # 强制转换为偶数
        
    h = (b - a) / n
    total = f(a) + f(b)
    
    # 奇数项系数为 4，偶数项系数为 2
    for i in range(1, n):
        x = a + i * h
        if i % 2 == 0:
            total += 2 * f(x)
        else:
            total += 4 * f(x)
            
    return total * h / 3.0


def integrate_hp(f: Callable[[Decimal], Decimal], a: Decimal, b: Decimal, n: int = 200) -> Decimal:
    """
    高精度定积分。
    由于 Decimal 计算开销较大，默认将 n 设为 200，在复合辛普森下已经能提供极高的精度。
    """
    a = Decimal(str(a)) if not isinstance(a, Decimal) else a
    b = Decimal(str(b)) if not isinstance(b, Decimal) else b
    
    if n % 2 != 0:
        n += 1
        
    n_dec = Decimal(n)
    h = (b - a) / n_dec
    total = f(a) + f(b)
    
    for i in range(1, n):
        x = a + Decimal(i) * h
        if i % 2 == 0:
            total += Decimal('2') * f(x)
        else:
            total += Decimal('4') * f(x)
            
    return total * h / Decimal('3')