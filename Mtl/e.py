def calculate_e_pure_int(precision):
    """完全不使用 Decimal 库，纯靠 Python 整数大数运算计算 e"""
    # 放大倍数，多给 10 位防止截断误差
    scale = 10 ** (precision + 10)

    total = scale  # k = 0 时的项
    item = scale  # 当前项

    k = 1
    while item > 0:
        item //= k  # 纯整数地板除
        total += item
        k += 1

    # 转换为字符串并精确截取
    e_str = str(total)
    # 插入小数点：第一位是 '2'，后面接小数点
    return e_str[0] + "." + e_str[1:precision]


