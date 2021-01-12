import sys

class Alignment(object):
    def __init__(self):
        self.W1 = 0
        self.W2 = 0
        self.W3 = 0
        self.X = 0
        self.Y = 0
        self.M = 0
        self.O = 0


def max2(a, b):
    return a if a > b else b


def max3(a, b, c):
    f = a if a > b else b
    return f if f > c else c


# 替换矩阵：match分值为10，mismatch分值为-5 和 0
FMatrix = [10, 0, -5, 0, 10, 0, 0, 0, -5, 0,
    0, 0, 10, 0, 0, 0, 0, 0, 0, -5,
    0, 0, 0, 0, 0, -5, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 10]
dict = {'A': 0, 'C': 2, 'G': 6, 'T': 19}
listset=['A','C','G','T']

def get_float_score(a, b):
    if a in listset and b in listset:
        i = dict.get(a)
        j = dict.get(b)
        return FMatrix[i+j]
    else:
        #print("ilegal DNA",a,b)
        return -5
    

global optimal
optimal = 0
global global_a
global global_b

def print_align(aligns, i, j, s1, s2, s1align, s2align, n):
    global optimal
    if optimal == 1:
        return
    result = aligns[i][j]
    if i == 0 and j == 0:  # 到矩阵单元(0, 0)才算结束，这代表初始的两个字符串的每个字符都被比对到了
        a = ''
        b = ''
        k = n-1
        while k >=0:
            a += s1align[k]
            b += s2align[k]
            k -= 1
        global global_a
        global_a = a
        global global_b
        global_b = b
        optimal = 1
        return
    if result.W1 != 0:  # 向上回溯一格
        s1align[n] = s1[i-1]
        s2align[n] = '-'
        print_align(aligns, i-1, j, s1, s2, s1align, s2align, n+1)
    if result.W2 != 0:  # 向左上回溯一格
        s1align[n] = s1[i-1]
        s2align[n] = s2[j-1]
        print_align(aligns, i-1, j-1, s1, s2, s1align, s2align, n+1)
    if result.W3 != 0:  # 向左回溯一格
        s1align[n] = '-'
        s2align[n] = s2[j-1]
        print_align(aligns, i, j-1, s1, s2, s1align, s2align, n+1)


# return score, s1, s2 modified
def align(s1, s2):
    # print(s1, s2, 333333)
    global optimal
    optimal = 0
    # 两个空位的罚分
    opening_gap = -4
    extending_gap = -1
    m = len(s1)
    n = len(s2)
    aligns = []
    # 创建二维对象, 高m+1 长n+1
    for i in range(m+1):
        aligns.append([])
        for j in range(n+1):
            x = Alignment()
            aligns[i].append(x)
    # aligns = [[Alignment() for i in range(n+1)] for j in range(m+1)]
    # 初始化
    aligns[0][0].X = opening_gap
    aligns[0][0].Y = opening_gap
    aligns[0][0].M = 0
    aligns[0][0].O = max3(aligns[0][0].X, aligns[0][0].Y, aligns[0][0].M)
    # 循环
    for i in range(1, m+1):
        aligns[i][0].X = opening_gap + (i - 1) * extending_gap
        aligns[i][0].Y = 2 * opening_gap + (i - 1) * extending_gap
        aligns[i][0].M = opening_gap + (i - 1) * extending_gap
        aligns[i][0].O = max3(aligns[i][0].X, aligns[i][0].Y, aligns[i][0].M)
        aligns[i][0].W1 = 1

    for j in range(1, n+1):
        aligns[0][j].X = 2 * opening_gap + (j - 1) * extending_gap
        aligns[0][j].Y = opening_gap + (j - 1) * extending_gap
        aligns[0][j].M = opening_gap + (j - 1) * extending_gap
        aligns[0][j].O = max3(aligns[0][j].X, aligns[0][j].Y, aligns[0][j].M)
        aligns[0][j].W3 = 1

    # 动态规划算法计算得分矩阵每个单元的分值
    # i = 0
    # while i <= m:
    #     j = 0
    #     while j <= n:
    #         # print(aligns[i][j])
    #         print(aligns[i][j].W1,aligns[i][j].W2,aligns[i][j].W3,aligns[i][j].X,aligns[i][j].Y,aligns[i][j].M,aligns[i][j].O)
    #         # print(i, j)
    #         j+=1
    #     i+=1
    # print(aligns[1][2].X)
    for i in range(1, m+1):
        for j in range(1, n+1):
            aligns[i][j].X = max2(aligns[i - 1][j].X + extending_gap, aligns[i - 1][j].M + opening_gap)
            aligns[i][j].Y = max2(aligns[i][j - 1].Y + extending_gap, aligns[i][j - 1].M + opening_gap)
            f = get_float_score(s1[i - 1], s2[j - 1])
            aligns[i][j].M = max3(aligns[i - 1][j - 1].X + f, aligns[i - 1][j - 1].Y + f, aligns[i - 1][j - 1].M + f)
            aligns[i][j].O = max3(aligns[i][j].X, aligns[i][j].Y, aligns[i][j].M)
            if aligns[i][j].O == aligns[i][j].X:
                aligns[i][j].W1 = 1
            if aligns[i][j].O == aligns[i][j].M:
                aligns[i][j].W2 = 1
            if aligns[i][j].O == aligns[i][j].Y:
                aligns[i][j].W3 = 1

    # print("序列匹对得分:", aligns[m][n].O)
    salign = [0 for i in range(m+n+1)]
    ralign = [0 for i in range(m+n+1)]
    # print('alignment内存使用情况: 二维对象数组:', sys.getsizeof(aligns), '其他:', sys.getsizeof(salign)+sys.getsizeof(ralign))
    print_align(aligns, m, n, s1, s2, salign, ralign, 0)
    # print(global_a, global_b)
    return aligns[m][n].O, global_a, global_b

if __name__ == '__main__':
    print("hello word")
    a,b,c=align('ATGCTTAGTGCACTCACGC','ATGCTTAGTGCACTCACGC')
    print("a",a)
    print("b",b)
    print("c",c)