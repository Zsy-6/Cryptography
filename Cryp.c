#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include <windows.h>
#include <wincrypt.h>
// 用于Windows加密API的库
#include <time.h>
#include <complex.h>
// 离散傅里叶变换测试用到的库

#define NUM_SAMPLES 1000 // 生成随机数的样本数量
#define RANGE 100         // 随机数范围 N-1
#define SubLen 5         // 重叠子序列长度

// 均匀分布的真随机数，N的范围内，生成一个
int generate_true_random(int N)
{
    // 使用CryptGenRandom生成一个随机整数。
    HCRYPTPROV hProvider;
    if (!CryptAcquireContext(&hProvider, NULL, NULL, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT))
    {
        fprintf(stderr, "CryptAcquireContext failed\n");
        exit(EXIT_FAILURE);
    }

    int random_number;
    if (!CryptGenRandom(hProvider, sizeof(random_number), (BYTE *)&random_number))
    {
        fprintf(stderr, "CryptGenRandom failed\n");
        CryptReleaseContext(hProvider, 0);
        exit(EXIT_FAILURE);
    }

    CryptReleaseContext(hProvider, 0);
    return abs(random_number) % N;
    // 将生成的整数取绝对值并对N取模，
    // 得到在 ([0, N-1]) 范围内的随机数
}

// 正态分布随机数
// 生成均值为mean，标准差为std_dev的正态分布随机数
double generateGaussian(double mean, double std_dev)
{
    static int hasSpare = 0;
    static double spare;
    if (hasSpare)
    {
        hasSpare = 0;
        return mean + std_dev * spare;
    }

    hasSpare = 1;
    double u, v, s;
    do
    {
        u = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
        v = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);

    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    return mean + std_dev * (u * s);
}
//生成[0,1]之间的均匀分布
double uniform_random() 
{
    return (double)rand() / (double)RAND_MAX;
}
// 生成标准正态分布的随机数
double normal_random() 
{
    double u1 = uniform_random();
    double u2 = uniform_random();
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159 * u2);
    return z0;
}
// 生成[0, N-1]之间的一个正态分布随机数
int normal_random_in_range(int N)
{
    double mean = (N - 1) / 2.0;
    double stddev = (N - 1) / 6.0;
    double z = normal_random();
    int result;
    do
    {
        result = (int)round(generateGaussian(mean, stddev));
    } while (result < 0 || result >= N);
    return result;
}
// 正态分布随机数生成代码结束

// 选择随机数函数，传入n,生成一个随机整数
int choice(int n)
{
    int num;
    switch (n)
    {
    case 1:
        num = rand() % RANGE;
        break;
    case 2:
        num = generate_true_random(RANGE);
        break;
    case 3:
        num = normal_random_in_range(RANGE);
        break;
    default:
        printf("Number is not 1, 2, or 3\n");
        break;
    }
    return num;
}

// 生成一个随机二进制数
int choice2(int n)
{
    int num;
    switch (n)
    {
    case 1:
        num = rand() % RANGE;
        break;
    case 2:
        num = generate_true_random(RANGE);
        break;
    case 3:
        num = normal_random_in_range(RANGE);
        break;
    default:
        printf("Number is not 1, 2, or 3\n");
        break;
    }
    return num % 2;
}

// 随机数的频率测试          1
// 传入n，选择随机数种类
void test_random_frequency(int n)
{
    printf("评率测试\n");
    int frequency[RANGE] = {0}; // 初始化频率数组为0
    int i, num;

    // 生成随机数并记录频率
    for (i = 0; i < NUM_SAMPLES; i++)
    {
        num = choice(n);
        frequency[num]++;
    }

    // 打印频率结果
    printf("Number | Frequency\n");
    printf("-----------------\n");
    for (i = 0; i < RANGE; i++)
    {
        printf("%6d | %9d\n", i, frequency[i]);
    }
    printf("评率测试结束\n");
}
// 频率测试代码结束

// 游程测试函数               0
// 传入二进制数组和n, n=sizeof(arr) / sizeof(arr[0])
void runsTest(int arr[], int n)
{
    printf("游程测试开始\n");
    int runs = 1; // 初始化游程数量
    int positive = 0, negative = 0;
    // 计算正数和负数的个数
    for (int i = 0; i < n; i++)
    {
        if (arr[i] == 1)
        {
            positive++;
        }
        else if (arr[i] == 0)
        {
            negative++;
        }
    }

    // 计算游程数量
    for (int i = 1; i < n; i++)
    {
        if (arr[i] != arr[i - 1])
        {
            runs++;
        }
    }

    // 期望的游程数量和标准差
    double expectedRuns = ((2 * positive * negative) / (double)n) + 1;
    double standardDeviation = sqrt((2 * positive * negative * (2 * positive * negative - n)) / (double)(n * n * (n - 1)));

    // Z值
    double z = (runs - expectedRuns) / standardDeviation;

    // 输出结果
    printf("Number of runs: %d\n", runs);
    printf("Expected number of runs: %f\n", expectedRuns);
    printf("Standard deviation: %f\n", standardDeviation);
    printf("Z value: %f\n", z);

    // 检验Z值是否在临界值范围内 (假设显著性水平为0.05，对应的Z值范围为-1.96到1.96)
    if (z > -1.96 && z < 1.96)
    {
        printf("The sequence is random.\n");
    }
    else
    {
        printf("The sequence is not random.\n");
    }
    printf("游程测试结束\n");
}
// 游程测试结束代码

// 离散傅里叶变换测试          1
// DFT函数
void DFT(double *inreal, double *inimag, double complex *out, int n)
{
    for (int k = 0; k < n; k++)
    { // 对于输出数组中的每个元素
        double complex sum = 0.0 + 0.0 * I;
        for (int t = 0; t < n; t++)
        { // 对于输入数组中的每个元素
            double angle = 2 * 3.14159 * t * k / n;
            sum += (inreal[t] + inimag[t] * I) * cexp(-I * angle);
        }
        out[k] = sum;
    }
}

// 时间复杂度高 数组选20规格
// 传入n，选择随机数类型
void fouri(int n)
{
    printf("离散傅里叶变换测试开始\n");
    int size = 20;
    double inreal[size];
    double inimag[size];
    double complex out[size];

    // 生成随机实部数组并初始化虚部为0
    for (int i = 0; i < size; i++)
    {
        inreal[i] = choice(n);
    }
    for (int i = 0; i < 20; i++)
    {
        inimag[i] = 0.0;
    }

    DFT(inreal, inimag, out, size); // 执行DFT

    // 打印输出数组
    printf("\nDFT output:\n");
    for (int i = 0; i < size; i++)
    {
        printf("%f + %fi\n", creal(out[i]), cimag(out[i]));
    }
    printf("离散傅里叶变换测试结束\n");
}
// 离散傅里叶变换测试函数结束

// 重叠子序列测试函数               1
//  计算卡方分布的p值
double chi2_p_value(double chi_square, int degrees_of_freedom)
{
    double k = degrees_of_freedom / 2.0;
    double x = chi_square / 2.0;
    double sum = exp(-x);
    double term = sum;

    for (int i = 1; i < k; i++)
    {
        term *= x / i;
        sum += term;
    }

    return 1.0 - sum;
}
// 重叠子序列测试函数,传入n选择 随机数种类
void Overlapping_Sub_Test(int n)
{
    printf("重叠子序列测试开始\n");
    // 生成随机二进制序列
    int sequence[NUM_SAMPLES];
    for (int i = 0; i < NUM_SAMPLES; i++)
    {
        sequence[i] = choice2(n);
    }

    int num_subsequences = (int)pow(2, SubLen);
    int subseq_counts[num_subsequences];
    for (int i = 0; i < num_subsequences; i++)
    {
        subseq_counts[i] = 0;
    }

    for (int i = 0; i <= NUM_SAMPLES - SubLen; i++)
    {
        int subseq = 0;
        for (int j = 0; j < SubLen; j++)
        {
            subseq = (subseq << 1) | sequence[i + j];
        }
        subseq_counts[subseq]++;
    }

    // 计算期望出现次数
    double expected_count = (double)(NUM_SAMPLES - SubLen + 1) / num_subsequences;

    // 计算卡方统计量
    double chi_square = 0.0;
    for (int i = 0; i < num_subsequences; i++)
    {
        double observed_count = subseq_counts[i];
        chi_square += (observed_count - expected_count) * (observed_count - expected_count) / expected_count;
    }

    // 计算p值
    int degrees_of_freedom = num_subsequences - 1;
    double p_value = chi2_p_value(chi_square, degrees_of_freedom);

    // 输出结果
    printf("Chi-square: %f\n", chi_square);
    printf("p-value: %f\n", p_value);
    if (p_value < 0.05)
    {
        printf("随机序列未通过重叠子序列测试。\n");
    }
    else
    {
        printf("随机序列通过了重叠子序列测试。\n");
    }
    printf("重叠子序列测试结束\n");
}
// 重叠子序列测试函数结束

// 随机数线性复杂度测试函数           1
// 传入n选择随机数种类
void linear_complexity(int n)
{
    printf("线性复杂度测试开始\n");
    int array[NUM_SAMPLES] = {0};
    for (int i = 0; i < NUM_SAMPLES; i++)
    {
        array[i] = choice2(n); // 生成0或1的随机数
    }
    int *b = (int *)calloc(NUM_SAMPLES, sizeof(int));
    int *c = (int *)calloc(NUM_SAMPLES, sizeof(int));
    b[0] = 1;
    c[0] = 1;
    int l = 0, m = -1, i = 0;

    for (int n = 0; n < NUM_SAMPLES; n++)
    {
        int d = array[n];
        for (int j = 1; j <= l; j++)
        {
            d ^= c[j] & array[n - j];
        }
        if (d == 1)
        {
            int *temp = (int *)calloc(NUM_SAMPLES, sizeof(int));
            for (i = 0; i < NUM_SAMPLES; i++)
            {
                temp[i] = c[i];
            }
            for (i = 0; (i + n - m) < NUM_SAMPLES; i++)
            {
                c[i + n - m] ^= b[i];
            }
            if (l <= n / 2)
            {
                l = n + 1 - l;
                m = n;
                for (i = 0; i < NUM_SAMPLES; i++)
                {
                    b[i] = temp[i];
                }
            }
            free(temp);
        }
    }
    free(b);
    free(c);
    printf("线性复杂度: %d\n", l);
    printf("线性复杂度测试结束\n");
}
// 线性复杂度测试函数结束
int main()
{
    int n = 0;
    int m;
    srand(time(NULL));
    printf("输入数字选择测试的随机数产生方式\n");
    printf("1为rand()产生的伪随机数\n");
    printf("2为均匀分布的真随机数\n");
    printf("3为正态分布的真随机数\n");
    scanf("%d", &n);
    if (n < 1 || n > 3)
    {
        printf("输入错误");
        return 0;
    }
    else
    {
        switch (n)
        {
            case 1:
            {   
                test_random_frequency(n);
                int a1[NUM_SAMPLES]={};
                for(int i=0;i<NUM_SAMPLES;i++)
                {
                    a1[i]=rand()%RANGE;
                }
                m=sizeof(a1) / sizeof(a1[0]);
                runsTest(a1, m);
                fouri(n);
                Overlapping_Sub_Test(n);
                linear_complexity(n);
                break;
            }
            case 2:
            {
                test_random_frequency(n);
                int a1[NUM_SAMPLES]={};
                for(int i=0;i<NUM_SAMPLES;i++)
                {
                    a1[i]=choice(2);
                }
                m=sizeof(a1) / sizeof(a1[0]);
                runsTest(a1, m);
                fouri(n);
                Overlapping_Sub_Test(n);
                linear_complexity(n);
                break;
            }
            default:
            {                
                test_random_frequency(n);
                int a1[NUM_SAMPLES]={};
                for(int i=0;i<NUM_SAMPLES;i++)
                {
                    a1[i]=choice2(3);
                }
                m=sizeof(a1) / sizeof(a1[0]);
                runsTest(a1, m);
                fouri(n);
                Overlapping_Sub_Test(n);
                linear_complexity(n);
                break;
            }
        }
    }
    return 0;
}