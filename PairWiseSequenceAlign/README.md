## 0. 文件结构介绍

- **data:**
  - sars-cov2.fasta:  SARS COV-2的30个菌株;
  - alignResult.txt: 存储对齐后的序列， 计算的得分以及程序运行情况；
  - result.txt: 对齐计算的得分以及程序运行情况；
  - resultMatrix.xls: 从result.txt 整理得到的 运行时间，得分，内存情况统计矩阵。
- **src:**
  - [STree4CS](https://github.com/pfmarteau/STree4CS):  提供后缀树构建算法，以及使用后缀树搜索最大公共子串；
  - aligment.py:  使用Needleman-Wunsch 算法进行局部对齐，并计算对齐后得分；
  - BinaryTree.py:  使用后缀树搜索最大公共子串，然后进行分割，构建二叉树，计算每个树节点局部序列匹配情况，然后通过中序遍历进行求和，得出最终的匹配序列以及对齐后序列的得分情况。
  - datavisual.py:  读取序列文件；
  - suffix_treesCon.py: 基于论文实现的代码，但是时间复杂度和空间占用比较大，作为临时文件；
  - datadistri.m:  使用matlab 绘图程序代码；

- 2020244002刘冬冬高级算法作业.doc：   project report;
- apt-2021.pdf： project assignment task requirement file;
- paper：相关文献记录
- ppt:  展示PPT;

## 1. 后缀树

### 1.1. 构建

> 后缀树是文本text所有后缀的压缩字典树。一个后缀树的生成经过一下几步：
> （1）生成给定文本text的所有后缀。
> （2）视所有后缀为有效单词，生成一个压缩字典树。

以”banana\0”（’\0’是结束字符）为例，字符串的所有后缀为：

```
banana\0
anana\0
nana\0
ana\0
na\0
a\0
\0
```

![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210110212147989.png)

![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210110212034737.png)

![Generalized Suffix Tree for Banana and Bonanza](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210112143124943.png)

### 1.2. 应用

（1）从目标串S中判断`是否包含`模式串T（Pattern Searching）

```
方案：用S构造后缀树，按在Trie中搜索子串的方法搜索T即可。
原理：若T在S中，则T必然是S的某个后缀的前缀。
例如：S = "leconte" T = "con"，查找T是否在S中,则T(con)必然是S(leconte)的后缀之一"conte"的前缀。
123
```

（2）从目标串S中查找串T`重复次数`

```
方案：用S+'$'构造后缀树，搜索T节点下的叶节点数目即为重复次数
原理：如果T在S中重复了两次，则S应有两个后缀以T为前缀，重复次数就自然统计出来了。
```

（3）从目标串S中查找`最长的重复子串`（Finding the longest repeated substring）

```
方案：原理同2，具体做法就是找到最深的非叶节点。
这个深是指从root所经历过的字符个数，最深非叶节点所经历的字符串起来就是最长重复子串。
为什么要非叶节点呢?因为既然是要重复，当然叶节点个数要>=2。
```

（4）从目标串T和S中查找`最长公共子串`（Finding the longest common substr ing）

```
方案：将S1#S2$作为字符串压入后缀树，找到最深的非叶节点，且该节点的叶节点既有#也有$(无#)。
```

（5）从目标串T中`查找最长的回文子串`（Finding the longest palindrome in a string）

### 1.3. 后缀树对齐算法

![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210111110749972.png)

1. build a suffix tree from one sequence, here S11 is assumed as the chosen one. and the tree's name is tree-S1;

2. pick out common strings:![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210111110921663.png)

   > The function named `select_prefix, which is a member function of suffix tree data structure, is applied to find the longest common prefix between a given input string and a suffix in the tree.` 

3. pick out common strings in S1 and S2 with sequence S2 and tree-S1 to form two common substring sets;

   - the starting positions of two matching substrings cannot be too far away;  `the start position should be less than their length`

4. aligns the remaining unmatched substrings between S1 &S2, and all substrings are concatenated to form the aligned sequence.

## 2. 序列比对算法

### 2.0. Fasta 文件格式

> 在生物信息学中，FASTA格式（又称为Pearson格式）是一种基于文本的、用于表示核苷酸序列或氨基酸序列的格式。在这种格式中碱基对或氨基酸用单个字母来表示，且允许在序列前添加序列名及注释。
>
> FASTA文件以序列表示和序列作为一个基本单元，各行记录信息如下：
>
> - 第一行是由大于号">"开头的任意文字说明，用于序列标记，为了保证后续分析软件能够区分每条序列，单个序列的标识必须具有唯一性。；
> - 从第二行开始为序列本身，只允许使用既定的核苷酸或氨基酸编码符号。通常核苷酸符号大小写均可，而氨基酸常用大写字母。使用时应注意有些程序对大小写有明确要求。文件每行的字母一般不应超过80个字符。

### 2.1. reads 截取

> NGS（二代测序）分析的起点往往是[fastq](https://en.wikipedia.org/wiki/FASTQ_format)文件。fastq文件其实就是一条条的记录，每个记录包含4行。其中比较重要的是第二行和第四行：第二行是测序得到的碱基序列，第四行是每个碱基相应的测序质量，测序质量越高代表该碱基被测错的概率越低，反之越高。在进行下游分析之前，常常要对fastq文件中的reads进行修剪（trim），将一条reads中测序质量不高的部分截掉。那么截取reads常用的策略有两种，**Fixed-length-trimming**以及**Phred-based-trimming**。
>

- **Fixed-length-trimming**： 截取固定长度的序列。一般来说，一条reads的头几个碱基和末尾几个碱基的测序质量比较差，所以你可以不加区分地将所有reads的前m个碱基以及后n个碱基去除。![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210110224420362.png)

- **Phred-based-trimming**: 将原始测序质量phred值都减去一个阈值（默认0.05）得到一系列新的数值（该新的序列有正值也有负值）。然后在该序列中找到和最大的子序列。 `最大子序列和问题`, 使用工具：[bwa](https://github.com/lh3/bwa)以及[seqtk](https://github.com/lh3/seqtk)；

### 2.2. 测试算法

- 蓄水池算法：fastq文件往往都很大，出于测试目的，我们经常要从fastq文件中随机抽取reads，生成一个小一点的fastq文件，以加快测试效率。![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210110232959528.png)

- **k-mer枚举![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210110233318087.png)**

### 2.3. 字典树

> Trie树，即字典树，又称单词查找树或键树，是一种树形结构，是一种哈希树的变种。典型应用是用于统计和排序大量的字符串（但不仅限于字符串），所以经常被搜索引擎系统用于文本词频统计。它的优点是：最大限度地减少无谓的字符串比较。

- **前缀树的创建**： ![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210110233700092.png)

> 1. 根节点不包含字符，除根节点外的每一个子节点都包含一个字符。
> 2. 从根节点到**某一个节点**，路径上经过的字符连接起来，为该节点对应的字符串。
> 3. 每个节点的所有子节点包含的字符互不相同。

> 插入和查询的效率很高，都为O(m)O(m)，其中 mm 是待插入/查询的字符串的长度。
>
> - 关于查询，会有人说 hash 表时间复杂度是O(1)不是更快？但是，哈希搜索的效率通常取决于 hash 函数的好坏，若一个坏的 hash 函数导致很多的冲突，效率并不一定比Trie树高。
> - Trie树中不同的关键字不会产生冲突。
> - Trie树只有在允许一个关键字关联多个值的情况下才有类似hash碰撞发生。
> - Trie树不用求 hash 值，对短字符串有更快的速度。通常，求hash值也是需要遍历字符串的。
> - Trie树可以对关键字按**字典序**排序。

- **应用---字符串检索**

> 从根节点开始一个一个字符进行比较：
>
> - 如果沿路比较，发现不同的字符，则表示该字符串在集合中不存在。
> - 如果所有的字符全部比较完并且全部相同，还需判断最后一个节点的标志位

- **词频统计**

> 对每一个关键字执行插入操作，若已存在，计数加1，若不存在，插入后`count`置1。

- **字符串排序**

> 对大量字符串按字典序进行排序，思路也很简单：遍历一次所有关键字，将它们全部插入trie树，树的每个结点的所有儿子很显然地按照字母表排序，然后先序遍历输出Trie树中所有关键字即可。

- **前缀匹配**

> trie树前缀匹配常用于搜索提示。如当输入一个网址，可以自动搜索出可能的选择。当没有完全匹配的搜索结果，可以返回前缀最相似的可能。

- **前缀树 增加删除操作**： https://juejin.cn/post/6844903750490914829

### 2.4.  全局比对Needleman-Wunsch

> 根据一个打分矩阵（替换矩阵）计算出两个序列比对最高得分的算法。
>
> **线性罚分**，即penalty＝g*d，其中g为空位长度，d为一个空位的罚分。
> **仿射型罚分**，即penalty=d+(g-1)*e， 其中g为连续空位的长度，d为连续空位中第一个空位的罚分，e为连续空位中第2个及以后空位的罚分。

- 打分矩阵：选用不同的打分矩阵或者罚分分值会导致比对结果不同，常用BLAST打分矩阵。
- 计算比对最高得分的算法：常用动态规划算法（Needleman-Wunsch算法）。![image-20210111084047244](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210111084047244.png)

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXSEQ 1000
#define GAP_CHAR '-'

// 对空位的罚分是仿射的
struct Unit {
    int W1;   // 是否往上回溯一格
    int W2;   // 是否往左上回溯一格
    int W3;   // 是否往左回溯一格
    float X;
    float Y;
    float M;
    float O;      // 得分矩阵第(i, j)这个单元的分值，即序列s(1,...,i)与序列r(1,...,j)比对的最高得分
};
typedef struct Unit *pUnit;

void strUpper(char *s);
float max2(float a, float b);
float max3(float a, float b, float c);
float getFScore(char a, char b);
void printAlign(pUnit** a, const int i, const int j, char* s, char* r, char* saln, char* raln, int n);
void align(char *s, char *r);

int main() {
    char s[MAXSEQ];
    char r[MAXSEQ];
    printf("The 1st seq: ");
    scanf("%s", s);
    printf("The 2nd seq: ");
    scanf("%s", r);
    align(s, r);
    return 0;
}

void strUpper(char *s) {
    while (*s != '\0') {
        if (*s >= 'a' && *s <= 'z') {
            *s -= 32;
        }
        s++;
    }
}

float max2(float a, float b) {
    return a > b ? a : b;
}

float max3(float a, float b, float c) {
    float f = a > b ? a : b;
    return f > c ? f : c;
}

// 替换矩阵：match分值为5，mismatch分值为-4
// 数组下标是两个字符的ascii码减去65之后的和
float FMatrix[] = {
    5, 0, -4, 0, 5, 0, -4, 0, -4, 0,
    0, 0, 5, 0, 0, 0, 0, 0, 0, -4,
    0, -4, 0, 0, 0, -4, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 5
};

float getFScore(char a, char b) {
    return FMatrix[a + b - 'A' - 'A'];
}

void printAlign(pUnit** a, const int i, const int j, char* s, char* r, char* saln, char* raln, int n) {
    int k;
    pUnit p = a[i][j];
    if (! (i || j)) {   // 到矩阵单元(0, 0)才算结束，这代表初始的两个字符串的每个字符都被比对到了
        for (k = n - 1; k >= 0; k--)
            printf("%c", saln[k]);
        printf("\n");
        for (k = n - 1; k >= 0; k--)
            printf("%c", raln[k]);
        printf("\n\n");
        return;
    }
    if (p->W1) {    // 向上回溯一格
        saln[n] = s[i - 1];
        raln[n] = GAP_CHAR;
        printAlign(a, i - 1, j, s, r, saln, raln, n + 1);
    }
    if (p->W2) {    // 向左上回溯一格
        saln[n] = s[i - 1];
        raln[n] = r[j - 1];
        printAlign(a, i - 1, j - 1, s, r, saln, raln, n + 1);
    }
    if (p->W3) {    // 向左回溯一格
        saln[n] = GAP_CHAR;
        raln[n] = r[j - 1];
        printAlign(a, i, j - 1, s, r, saln, raln, n + 1);
    }
}

void align(char *s, char *r) {
    int i, j;
    int m = strlen(s);
    int n = strlen(r);
    float d = -7;     // 对第一个空位的罚分
    float e = -2;     // 第二个及以后空位的罚分
    pUnit **aUnit;
    char* salign;
    char* ralign;
    float f;
    // 初始化
    if ((aUnit = (pUnit **) malloc(sizeof(pUnit*) * (m + 1))) == NULL) {
        fputs("Error: Out of space!\n", stderr);
        exit(1);
    }
    for (i = 0; i <= m; i++) {
        if ((aUnit[i] = (pUnit *) malloc(sizeof(pUnit) * (n + 1))) == NULL) {
            fputs("Error: Out of space!\n", stderr);
            exit(1);     
        }
        for (j = 0; j <= n; j++) {
            if ((aUnit[i][j] = (pUnit) malloc(sizeof(struct Unit))) == NULL) {
                fputs("Error: Out of space!\n", stderr);
                exit(1);     
            }
            aUnit[i][j]->W1 = 0;
            aUnit[i][j]->W2 = 0;
            aUnit[i][j]->W3 = 0;
        }
    }
    aUnit[0][0]->X = d;
    aUnit[0][0]->Y = d;
    aUnit[0][0]->M = 0;
    aUnit[0][0]->O = max3(aUnit[0][0]->X, aUnit[0][0]->Y, aUnit[0][0]->M);
    for (i = 1; i <= m; i++) {
        aUnit[i][0]->X = d + (i - 1) * e;
        aUnit[i][0]->Y = 2 * d + (i - 1) * e;
        aUnit[i][0]->M = d + (i - 1) * e;
        aUnit[i][0]->O = max3(aUnit[i][0]->X, aUnit[i][0]->Y, aUnit[i][0]->M);
        aUnit[i][0]->W1 = 1;
    }
    for (j = 1; j <= n; j++) {
        aUnit[0][j]->X = 2 * d + (j - 1) * e;
        aUnit[0][j]->Y = d + (j - 1) * e;
        aUnit[0][j]->M = d + (j - 1) * e;
        aUnit[0][j]->O = max3(aUnit[0][j]->X, aUnit[0][j]->Y, aUnit[0][j]->M);
        aUnit[0][j]->W3 = 1;
    }
    // 将字符串都变成大写
    strUpper(s);
    strUpper(r);
    // 动态规划算法计算得分矩阵每个单元的分值
    for (i = 1; i <= m; i++) {
        for (j = 1; j <= n; j++) {
            aUnit[i][j]->X = max2(aUnit[i - 1][j]->X + e, aUnit[i - 1][j]->M + d);
            aUnit[i][j]->Y = max2(aUnit[i][j - 1]->Y + e, aUnit[i][j - 1]->M + d);
            f = getFScore(s[i - 1], r[j - 1]);
            aUnit[i][j]->M = max3(aUnit[i - 1][j - 1]->X + f, aUnit[i - 1][j - 1]->Y + f, aUnit[i - 1][j - 1]->M + f);
            aUnit[i][j]->O = max3(aUnit[i][j]->X, aUnit[i][j]->Y, aUnit[i][j]->M);
            if (aUnit[i][j]->O == aUnit[i][j]->X) aUnit[i][j]->W1 = 1;
            if (aUnit[i][j]->O == aUnit[i][j]->M) aUnit[i][j]->W2 = 1;
            if (aUnit[i][j]->O == aUnit[i][j]->Y) aUnit[i][j]->W3 = 1;
        }
    }
/*
    // 打印得分矩阵
    for (i = 0; i <= m; i++) {
        for (j = 0; j <= n; j++)
            printf("%f ", aUnit[i][j]->O);
        printf("\n");
    }
*/
    printf("max score: %f\n", aUnit[m][n]->O);
    // 打印最优比对结果，如果有多个，全部打印
    // 递归法
    if ((salign = (char*) malloc(sizeof(char) * (m + n + 1))) == NULL) {
        fputs("Error: Out of space!\n", stderr);
        exit(1);
    }
    if ((ralign = (char*) malloc(sizeof(char) * (m + n + 1))) == NULL) {
        fputs("Error: Out of space!\n", stderr);
        exit(1);
    }
    printAlign(aUnit, m, n, s, r, salign, ralign, 0);
    // 释放内存
    free(salign);
    free(ralign);
    for (i = 0; i <= m; i++) {
        for (j = 0; j <= n; j++) {
            free(aUnit[i][j]);
        }
        free(aUnit[i]);
    }
    free(aUnit);
}
```

```python
a,b,c=align('AGTGCACTCACGC','ATGCTTAGTGCACTCACGC')
print("a",a)
print("b",b)
print("c",c)
```

![image-20210111164251774](C:/Users/dell/AppData/Roaming/Typora/typora-user-images/image-20210111164251774.png)

### 2.5. 局部匹配Smith-Waterman

![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210111084326579.png)

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXSEQ 1000
#define GAP_CHAR '-'
// 对空位的罚分是线性的
struct Unit {
    int W1;   // 是否往上回溯一格
    int W2;   // 是否往左上回溯一格
    int W3;   // 是否往左回溯一格
    float M;      // 得分矩阵第(i, j)这个单元的分值，即序列s(1,...,i)与序列r(1,...,j)比对的最高得分
};
typedef struct Unit *pUnit;
void strUpper(char *s);
float max4(float a, float b, float c, float d);
float getFScore(char a, char b);
void printAlign(pUnit** a, const int i, const int j, char* s, char* r, char* saln, char* raln, int n);
void align(char *s, char *r);
int main() {
    char s[MAXSEQ];
    char r[MAXSEQ];
    printf("The 1st seq: ");
    scanf("%s", s);
    printf("The 2nd seq: ");
    scanf("%s", r);
    align(s, r);
    return 0;
}
void strUpper(char *s) {
    while (*s != '\0') {
        if (*s >= 'a' && *s <= 'z') {
            *s -= 32;
        }
        s++;
    }
}
float max4(float a, float b, float c, float d) {
    float f = a > b ? a : b;
    float g = c > d ? c : d;
    return f > g ? f : g;
}
// 替换矩阵：match分值为5，mismatch分值为-4
// 数组下标是两个字符的ascii码减去65之后的和
float FMatrix[] = {
    5, 0, -4, 0, 5, 0, -4, 0, -4, 0,
    0, 0, 5, 0, 0, 0, 0, 0, 0, -4,
    0, -4, 0, 0, 0, -4, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 5
};
float getFScore(char a, char b) {
    return FMatrix[a + b - 'A' - 'A'];
}
void printAlign(pUnit** a, const int i, const int j, char* s, char* r, char* saln, char* raln, int n) {
    int k;
    pUnit p = a[i][j];
    if (p->M == 0) {
        for (k = n - 1; k >= 0; k--)
            printf("%c", saln[k]);
        printf("\n");
        for (k = n - 1; k >= 0; k--)
            printf("%c", raln[k]);
        printf("\n\n");
        return;
    }
    if (p->W1) {    // 向上回溯一格
        saln[n] = s[i - 1];
        raln[n] = GAP_CHAR;
        printAlign(a, i - 1, j, s, r, saln, raln, n + 1);
    }
    if (p->W2) {    // 向左上回溯一格
        saln[n] = s[i - 1];
        raln[n] = r[j - 1];
        printAlign(a, i - 1, j - 1, s, r, saln, raln, n + 1);
    }
    if (p->W3) {    // 向左回溯一格
        saln[n] = GAP_CHAR;
        raln[n] = r[j - 1];
        printAlign(a, i, j - 1, s, r, saln, raln, n + 1);
    }
}
void align(char *s, char *r) {
    int i, j;
    int m = strlen(s);
    int n = strlen(r);
    float gap = -2.5;     // 对空位的罚分
    float m1, m2, m3, maxm;
    float maxMatrix;    // 得分矩阵中的最高分
    pUnit **aUnit;
    char* salign;
    char* ralign;
    // 初始化
    if ((aUnit = (pUnit **) malloc(sizeof(pUnit*) * (m + 1))) == NULL) {
        fputs("Error: Out of space!\n", stderr);
        exit(1);
    }
    for (i = 0; i <= m; i++) {
        if ((aUnit[i] = (pUnit *) malloc(sizeof(pUnit) * (n + 1))) == NULL) {
            fputs("Error: Out of space!\n", stderr);
            exit(1);     
        }
        for (j = 0; j <= n; j++) {
            if ((aUnit[i][j] = (pUnit) malloc(sizeof(struct Unit))) == NULL) {
                fputs("Error: Out of space!\n", stderr);
                exit(1);     
            }
            aUnit[i][j]->W1 = 0;
            aUnit[i][j]->W2 = 0;
            aUnit[i][j]->W3 = 0;
        }
    }
    aUnit[0][0]->M = 0;
    for (i = 1; i <= m; i++) {
        aUnit[i][0]->M = 0;
    }
    for (j = 1; j <= n; j++) {
         aUnit[0][j]->M = 0;
    }
    // 将字符串都变成大写
    strUpper(s);
    strUpper(r);
    // 动态规划算法计算得分矩阵每个单元的分值
    for (i = 1; i <= m; i++) {
        for (j = 1; j <= n; j++) {
            m1 = aUnit[i - 1][j]->M + gap;
            m2 = aUnit[i - 1][j - 1]->M + getFScore(s[i - 1], r[j - 1]);
            m3 = aUnit[i][j - 1]->M + gap;
            maxm = max4(m1, m2, m3, 0);
            aUnit[i][j]->M = maxm;
            if (maxm != 0) { 
                if (m1 == maxm) aUnit[i][j]->W1 = 1;
                if (m2 == maxm) aUnit[i][j]->W2 = 1;
                if (m3 == maxm) aUnit[i][j]->W3 = 1;
            }
        }
    }
/*
    // 打印得分矩阵
    for (i = 0; i <= m; i++) {
        for (j = 0; j <= n; j++)
            printf("%f ", aUnit[i][j]->M);
        printf("\n");
    }
*/
    // 求取得分矩阵中的最高分
    maxMatrix = 0;
    for (i = 1; i <= m; i++) {
        for (j = 1; j <= n; j++) {
            if (aUnit[i][j]->M > maxMatrix)
                maxMatrix = aUnit[i][j]->M;
        }
    }
    printf("max score: %f\n", maxMatrix);
    // 打印最优比对结果，如果有多个，全部打印
    // 递归法
    if (maxMatrix == 0) {
        fputs("No seq aligned.\n", stdout);
    } else {
        if ((salign = (char*) malloc(sizeof(char) * (m + n + 1))) == NULL) {
            fputs("Error: Out of space!\n", stderr);
            exit(1);
        }
        if ((ralign = (char*) malloc(sizeof(char) * (m + n + 1))) == NULL) {
            fputs("Error: Out of space!\n", stderr);
            exit(1);
        }
        for (i = m; i >= 1; i--)
            for (j = n; j >= 1; j--)
                if (aUnit[i][j]->M == maxMatrix)
                    printAlign(aUnit, i, j, s, r, salign, ralign, 0);
        // 释放内存
        free(salign);
        free(ralign);
    }
    for (i = 0; i <= m; i++) {
        for (j = 0; j <= n; j++) {
            free(aUnit[i][j]);
        }
        free(aUnit[i]);
    }
    free(aUnit);
}
```

### 2.6. 交叉匹配

> 不同于全局匹配，交叉匹配中两端的序列可以不参与联配（或者说不罚分）。不同于局部匹配，交叉匹配中某一条序列的头部必须参与联配且某一条序列的尾部必须参加联配。

### 2.7. BWT算法

> BWT算法可以分为编码和解码两部分。编码后，原始字符串中的相似字符会处在比较相邻的位置；解码就是将编码后的字符串重新恢复成原始字符串的过程。BWT的一个特点就是经过编码后的字符串可以完全恢复成原始字符串。

![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210111090501197.png)

![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210111090610254.png)

![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210111090629729.png)

![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210111091208074.png)

![](https://gitee.com/github-25970295/blogImage/raw/master/img/image-20210111091006197.png)

## 3. coding

- 元素合法判断

```python
listset=['A','C','G','T']   # 实验中序列中含有N， 不在范围内
def get_float_score(a, b):
    if a in listset and b in listset:
        i = dict.get(a)
        j = dict.get(b)
        return FMatrix[i+j]
    else:
        print("ilegal DNA",a,b)
        return -5
```

- None 和“” 区别

> 是不同的一种数据类型；type(None) <class 'NoneType'>；type('') <class ''str'>；
>
> 判断的时候 均是False；ff=None  >>> if ff:    print('ff is define') ； 无输出
>
> 属性不同，dir(None) 查看；

- 内存资源查看

```python
import psutil
import os
print(u'当前进程的内存使用：%.4f GB' % (psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024) )
#当前进程的内存使用：0.6119 GB
info = psutil.virtual_memory()
info
#svmem(total=8506298368, available=3768840192, percent=55.7, used=4737458176, free=3768840192)
print( u'电脑总内存：%.4f GB' % (info.total / 1024 / 1024 / 1024) )
print(u'当前使用的总内存占比：',info.percent)
print(u'cpu个数：',psutil.cpu_count())

#查看所有程序内存
import psutil
# 查看所有进程
pid = psutil.pids()
for k, i in enumerate(pid):
    try:
        proc = psutil.Process(i)
        print(k, i, "%.2f%%" % (proc.memory_percent()), "%", proc.name(), proc.exe())
    except psutil.AccessDenied:
        print("无")
```

- 程序执行时间查看

```python
from time import clock
def _time_analyze_():
    start = clock()
    pass
    finish = clock()
```

- 二叉树建立

```python
class node():
	def __init__(self,k=None,l=None,r=None):
		self.key=k;
		self.left=l;
		self.right=r;
 
def create(root):
	a=raw_input('enter a key:');
	if a is '#':
		root=None;
	else:
		root=node(k=a);
		root.left=create(root.left);
		root.right=create(root.right);
	return root;
 
def preorder(root):      #前序遍历
	if root is None:
		return ;
	else :
		print root.key;
		preorder(root.left);
		preorder(root.right);
 
def inorder(root):     #中序遍历
	if root is None:
		return ;
	else:
		inorder(root.left);
		print root.key;
		inorder(root.right);
 
def postorder(root):   # 后序遍历
	if root is None:
		return ;
	else :
		postorder(root.left);
		postorder(root.right);
		print root.key;
		
root=None;     # 测试代码
root=create(root);
preorder(root);
inorder(root);
postorder(root);
```

```python
class PairsequenceItem:
    def __init__(self,seq1="",seq2="",start1=0,start2=0,common=False):
        self.start1=start1
        self.start2=start2
        if len(seq1)==0:
            self.seq1=""
        else:
            self.seq1=seq1
        if len(seq2)==0:
            self.seq2=""
        else:
            self.seq2=seq2
        self.common=common
    # return score, s1,s2 aligned
    def getPairScoreAlignSequence(self):
        if self.common:
            return len(self.seq1)*10,self.seq1,self.seq1
        else:
            return align("".join(self.seq1),"".join(self.seq2))


class Node:
    def __init__(self, data):
        self.left = None
        self.right = None
        self.data = data

def createTree(seq1,seq2):
    if len(seq1)==0 and len(seq2)==0:
        return None
    elif len(seq1)>200 and len(seq2)>200:  
        try:
            #print("seq1",len(seq1),"seq2",len(seq2))
            st=STree.STree4CS([seq1,seq2])
            #print("FindLSC")
            start,end=st.lcsEndStart()
            st1 = STree.STree4CS([seq2])
            start1=st1.find(seq1[start:end])
            #print("construct PairsequenceItem")
            pair=PairsequenceItem(seq1[start:end],seq2[start1:start1+end-start],start,start1,common=True)
            root=Node(data=pair)
            #print("construct left node")
            root.left=createTree(seq1[0:start],seq2[0:start1])
            #print("construct right node")
            root.right=createTree(seq1[end:-1],seq2[start1+end-start:-1])
            #print("return node")
            return root
        except :
            print("error seq1",len(seq1),"seq2",len(seq2))
    else:
        #print("uncommon sequence:\tseq1",seq1,"seq2",seq2)
        pa =PairsequenceItem(seq1,seq2,0,0,common=False)
        return Node(data=pa)


def inorderTraversal(root):   #先把根节点入栈，如果左子树一直不为空，就一直入栈，直到 
                                ##把所有左节点入栈，然后pop栈顶元素，指针指向栈顶元素的右子树
    """
    :type root: TreeNode
    :rtype: List[int]
    """
    stack=[]
    res=[]
    if not root:
        return []
    while root or stack:
        while root:
            stack.append(root)
            root=root.left
        if stack:
            a=stack.pop()
            root=a.right
            res.append(a.data)
    return res
```

## 4. 学习资源

- Fasta文件介绍：https://www.jianshu.com/p/cd232d34c408
- https://www.geeksforgeeks.org/pattern-searching-using-suffix-tree/
- 后缀树使用代码：https://github.com/pfmarteau/STree4CS
- 字典树：https://blog.csdn.net/lisonglisonglisong/article/details/45584721
- BWT算法： https://blog.csdn.net/biocity/article/details/102081350
- 二叉树： https://blog.csdn.net/u011608357/article/details/26075069
- pairwise-alignment-in-python： https://github.com/alevchuk/pairwise-alignment-in-python