import math
from decimal import Decimal, getcontext


def calculate_pi_chudnovsky(precision):
    """使用丘德诺夫斯基算法计算指定位数的 pi"""
    # 设置 decimal 的全局精度（稍微多设几位防止四舍五入误差）
    getcontext().prec = precision + 10

    # 算法中的代数常数
    C = 426880 * Decimal(10005).sqrt()
    M = Decimal(1)
    L = Decimal(13591409)
    X = Decimal(1)
    K = Decimal(6)
    S = Decimal(13591409)

    # 计算需要迭代的步数 (每步大约 14 位精度)
    terms = math.ceil(precision / 14)

    print(f"设定目标精度: {precision} 位")
    print(f"预计需要迭代步数: {terms} 步\n" + "-" * 40)

    # 开始纯代数迭代
    for k in range(1, terms):
        # 优化项：利用前一项的状态递推计算大阶乘和高次幂，避免直接算 (6k)!
        M = M * (K**3 - 16 * K) // (k**3)
        L += 545140134
        X *= -262537412640768000
        S += (M * L) / X
        K += 12

    # 最终求出 pi
    pi = C / S

    # 截取到用户要求的精确位数
    getcontext().prec = precision
    return +pi  # 通过正号触发精度截断


