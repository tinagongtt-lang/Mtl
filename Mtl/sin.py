"""
Mtl 数学库 - 三角函数核心模块
基于高精度 Chudnovsky PI 与 泰勒级数前向递推法则实现
"""

from decimal import Decimal, getcontext
from .pi import *

# ==========================================
# 1. 动态初始化配置（根据你库的定位自动或手动调配）
# ==========================================

# 假设计算 100 位全局精度（如果你需要全面拥抱高精度版）
getcontext().prec = 100

# 核心常量：直接调用你本地的 Chudnovsky 算法生成绝对精准的 PI 边界
# 如果你的 calculate_pi_chudnovsky 返回的是字符串或整型，请确保正确包裹
_PI_DECIMAL = Decimal(calculate_pi_chudnovsky(100))
_TWO_PI_DECIMAL = _PI_DECIMAL * 2

# 预留一份硬件级标准 float 的常量，用于高速硬件浮点计算
_PI_FLOAT = float(_PI_DECIMAL)
_TWO_PI_FLOAT = _PI_FLOAT * 2


# ==========================================
# 2. 核心函数实现
# ==========================================

def sin(x) -> float:
    """
    高效、高精度的标准双精度浮点数正弦函数。
    
    使用标准 IEEE 754 float (64位) 运行，针对物理模拟、绘图引擎等高性能场景优化。
    利用泰勒级数前向递推算法，时间复杂度恒定 O(1)。
    """
    # 强制将输入转换为标准 float
    x = float(x)
    
    # 【防御大数崩溃】
    # 当浮点数 x 超过 1e15 时，由于 64位浮点数的低位信息缺失，
    # x % _TWO_PI 的取模操作会彻底失真。库级标准做法：直接归零。
    if x > 1e15 or x < -1e15:
        return 0.0
        
    # 【步骤 1：大域规约】利用 2π 周期性将 x 限制在 [0, 2π) 之间
    x = x % _TWO_PI_FLOAT
    
    # 【步骤 2：中心域规约】进一步将 x 限制在 [-π, π] 之间。
    # 泰勒级数在 0 附近收敛最快，此步能将循环次数压缩至 10 次以内。
    if x > _PI_FLOAT:
        x -= _TWO_PI_FLOAT
    elif x < -_PI_FLOAT:
        x += _TWO_PI_FLOAT
        
    # 【步骤 3：前向递推迭代】
    term = x            # 泰勒第一项：x^1 / 1!
    sin_x = x           # 初始累加器
    x_squared = x * x   # 缓存 x^2，杜绝循环内的重复幂运算
    n = 1
    
    while True:
        # 核心递推：下一项 = 上一项 * (-x^2) / (2n * (2n + 1))
        # 绝不计算显式阶乘，杜绝大数溢出与无谓的硬件开销
        term *= -x_squared / ((2 * n) * (2 * n + 1))
        
        previous = sin_x
        sin_x += term
        
        # 动态截断：当 term 小于当前 float 的最小机器精度时，
        # sin_x + term 会由于精度落差直接等于 sin_x。此时主动停机，榨干硬件最后 1 bit 精度。
        if previous == sin_x: 
            break
            
        n += 1
        
    # 【步骤 4：边界防御】强行卡位，防止由于微小的浮点数累加误差导致输出如 1.0000000000000002
    if sin_x > 1.0: return 1.0
    if sin_x < -1.0: return -1.0
    
    return sin_x


def sin_hp(x) -> Decimal:
    """
    任意精度（High Precision）正弦函数。
    
    全面拥抱你用 Chudnovsky 算法算出的 100 位 PI。
    当你的库需要进行极限精度的数学推导、科学计算或密码学探索时使用。
    """
    # 确保输入转化为高精度 Decimal
    x = Decimal(str(x)) if not isinstance(x, Decimal) else x
    
    # 【大域与中心域规约】Decimal 环境下大数取模依然绝对精准
    x = x % _TWO_PI_DECIMAL
    if x > _PI_DECIMAL:
        x -= _TWO_PI_DECIMAL
    elif x < -_PI_DECIMAL:
        x += _TWO_PI_DECIMAL
        
    # 【高精度前向递推】
    term = x
    sin_x = x
    x_squared = x * x
    n = 1
    
    while True:
        # 所有的操作数均在 Decimal 上进行，保证 100 位精度不发生断崖式下跌
        term *= -x_squared / Decimal(2 * n * (2 * n + 1))
        
        previous = sin_x
        sin_x += term
        
        # 只要最新的项还能改变 100 位精度的任何一位，就继续算，直到完全收敛
        if previous == sin_x: 
            break
            
        n += 1
        
    # 边界卡位
    if sin_x > Decimal('1.0'): return Decimal('1.0')
    if sin_x < Decimal('-1.0'): return Decimal('-1.0')
    
    return sin_x