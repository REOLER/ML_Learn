import numpy as np
import matplotlib.pyplot as plt


class GA():
    """求出二进制编码的长度"""

    def __init__(self):
        self.boundsBegin = -2
        self.boundsEnd = 3
        # 运算精确度
        precision = 0.0001
        # %染色体长度
        self.BitLength = int(np.log2((self.boundsEnd - self.boundsBegin) / precision)) + 1
        # 初始种群大小
        self.initSize = 50
        # 最大进化代数
        self.GenerationMax = 12
        # 交叉概率
        self.pCrossOver = 0.90
        # 变异概率
        self.p_mutation = 0.2
        # 随机生成initSize行长度为BitLength的列表
        self.population = np.random.randint(0, 2, size=(self.initSize, self.BitLength))

    """计算出适应度"""

    def fitness(self, population):
        fit_value = []
        cumulative = []
        for i in population:
            # 个体的染色体由二进制对应的十进制
            x = self.transform2to10(i)
            # 十进制数映射到函数的定义域
            geneX = self.boundsBegin + x * (self.boundsEnd - self.boundsBegin) / (pow(2, self.BitLength) - 1)
            # 定义域上的数代入函数表达式求得适应值
            s = self.targetFun(geneX)
            # 将适应值依次记录到列表fit_value中
            fit_value.append(s)
        # 适应值的总和
        fit_sum = sum(fit_value)
        # 计算每个个体的适应性比例
        every_population = [k / fit_sum for k in fit_value]
        # 计算累计适应性比例,记录到列表cumulative中
        cumulative.append(every_population[0])
        every_population.remove(every_population[0])
        for j in every_population:
            p = cumulative[-1] + j
            cumulative.append(p)
        return fit_value, cumulative

    """选择两个基因，准备交叉"""

    def select(self, cumulative):
        sel_gene_num = []
        for i in range(2):
            j = 1
            r = np.random.uniform(0, 1)
            p_rank = [x - r for x in cumulative]
            while p_rank[j] < 0:
                j = j + 1
            sel_gene_num.append(j)
        return sel_gene_num

    """交叉"""

    def crossover(self, sel_gene_num, pc):
        d = self.population[sel_gene_num[0]].copy()
        f = self.population[sel_gene_num[1]].copy()
        crossRank = np.random.uniform()
        if crossRank < pc:
            print('Can exchange chromosomes')
            # 交叉点不能为首或尾
            c = np.random.randint(1, self.BitLength - 1)
            print('parent:')
            print(d)
            print(f)
            print('cross area:' + str(c))
            print('new individual:')
            a = self.population[sel_gene_num[0]][c:]
            b = self.population[sel_gene_num[1]][c:]
            d[c:] = b
            f[c:] = a
            print(d)
            print(f)
            g = d
            h = f
        else:
            g = self.population[sel_gene_num[0]]
            h = self.population[sel_gene_num[1]]
        return g, h

    """变异操作"""

    def mutation(self, sc_new, p_mutation):
        r = np.random.uniform(0, 1)
        if r < p_mutation:
            v = np.random.randint(0, self.BitLength)
            print('genetic mutation, location:' + str(v))
            sc_new[v] = abs(sc_new[v] - 1)
        else:
            sc_new = sc_new
        return sc_new

    """二进制转换为十进制"""

    def transform2to10(self, population):
        # x=population[-1]  #最后一位的值
        x = 0
        # n=len(population)
        n = self.BitLength
        p = population.copy()
        p = p.tolist()
        p.reverse()
        for j in range(n):
            x = x + p[j] * pow(2, j)
        return x  # 返回十进制的数

    """目标函数"""

    def targetFun(self, x):
        y = x * (np.sin(10 * np.pi * x)) + 2
        return y


if __name__ == '__main__':
    GenerationMax = 12
    ga = GA()
    sc_new = []
    yMax = []
    geneX = 0
    # print(ga.population)
    fit_value, cumulative = ga.fitness(ga.population)
    Generation = 1
    while Generation < GenerationMax + 1:
        fit_value, cumulative = ga.fitness(ga.population)
        for j in range(0, ga.initSize, 2):
            # 返回选中的2个个体的序号
            sel_gene_num = ga.select(cumulative)
            # 两条染色体交叉
            sCross = ga.crossover(sel_gene_num, ga.pCrossOver)
            # 两条染色体变异
            s1 = ga.mutation(sCross[0], ga.p_mutation)
            s2 = ga.mutation(sCross[1], ga.p_mutation)
            sc_new.append(s1)
            sc_new.append(s2)
        ga.population = sc_new
        fit_value, cumulative = ga.fitness(ga.population)
        fMax = max(fit_value)
        d = fit_value.index(fMax)
        yMax.append(fMax)
        x = ga.transform2to10(ga.population[d])
        geneX = ga.boundsBegin + x * (ga.boundsEnd - ga.boundsBegin) / (pow(2, ga.BitLength) - 1)
        Generation = Generation + 1
    bestPopulation = geneX
    targetMax = ga.targetFun(geneX)
    print(geneX)
    print(targetMax)

    x = np.linspace(-2, 3, 500)
    print(type(x))
    y = x * (np.sin(10 * np.pi * x)) + 2
    print(type(y))
    plt.scatter(geneX, targetMax, c='red')
    plt.xlim(0, 5)
    plt.ylim(0, 6)
    plt.plot(x, y)
    plt.annotate('local max', xy=(geneX, targetMax), xytext=(4, 5.5), arrowprops=dict(facecolor='black', shrink=0.05))
    plt.show()
