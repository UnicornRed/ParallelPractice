#include <iostream>
#include <cmath>
#include <omp.h>

using namespace std;

//Структура пары
struct My_pair
{
    double first;
    double second;
};

//Класс кубического полинома
class Cube_poly
{
private:
    //Коэффициенты кубического полинома
    double coeff[4];
    //Начало промежутка задания кубического полинома
    double a;
    //Конец промежутка задания кубического полинома
    double b;
public:
    //Кубический полином по-умолчанию
    Cube_poly() : a(0), b(0)
    {
        for (int i = 0; i < 4; i++)
            coeff[i] = 0;
    }

    //Кубический полином по коэффициентам
    Cube_poly(double coeff[]) : a(0), b(0)
    {
        for (int i = 0; i < 4; ++i)
            (*this).coeff[i] = coeff[i];
    }

    //Изменение коэффициентов и промежутка задания кубического полинома
    //по построению кубического полинома Эрмита
    //f1 - значение функции в начале промежутка
    //f2 - значение функции в конце промежутка
    //m1 - значение первой производной функции в начале промежутка
    //m2 - значение первой производной функции в конце промежутка
    //h - длина промежутка
    //x1 - начало промежутка
    //x2 - конец промежутка
    void constructor (double f1, double f2, double m1, double m2, double h, double x1, double x2)
    {
        double tCoeff[4] = { -x1 / h, 1.0 / h, 0, 0 };

        Cube_poly t = Cube_poly(tCoeff);
        Cube_poly poly;
        //Построение кубического полинома Эрмита
        poly = (t - 1.0) * (t - 1.0) * (t * 2.0 + 1) * f1 - t * t * (t * 2.0 - 3.0) * f2 +
               (t - 1) * (t - 1) * t * m1 * h + t * t * (t - 1) * m2 * h;

        for (int i = 0; i < 4; ++i)
            coeff[i] = poly.coeff[i];

        a = x1;
        b = x2;
    }

    //Переопределение сложения для кубических полиномов
    Cube_poly operator+(Cube_poly b)
    {
        double coeff[4];
        Cube_poly* a = this;

        for (int i = 0; i < 4; ++i)
            coeff[i] = a->get_coeff(i) + b.get_coeff(i);

        Cube_poly poly(coeff);

        return poly;
    }

    //Переопределение сложения для кубического полинома и числа
    Cube_poly operator+(double b)
    {
        double coeff[4];
        Cube_poly* a = this;

        for (int i = 0; i < 4; ++i)
            coeff[i] = a->get_coeff(i);

        coeff[0] += b;

        Cube_poly poly(coeff);

        return poly;
    }

    //Переопределение вычитания для кубических полиномов
    Cube_poly operator-(Cube_poly b)
    {
        double coeff[4];
        Cube_poly* a = this;

        for (int i = 0; i < 4; ++i)
            coeff[i] = a->get_coeff(i) - b.get_coeff(i);

        Cube_poly poly(coeff);

        return poly;
    }

    //Переопределение вычитания для кубического полинома и числа
    Cube_poly operator-(double b)
    {
        double coeff[4];
        Cube_poly* a = this;

        for (int i = 0; i < 4; ++i)
            coeff[i] = a->get_coeff(i);

        coeff[0] -= b;

        Cube_poly poly(coeff);

        return poly;
    }

    //Переопределение унарного минуса для кубического полинома
    Cube_poly operator-()
    {
        double coeff[4];
        Cube_poly* a = this;

        for (int i = 0; i < 4; ++i)
            coeff[i] = -a->get_coeff(i);

        Cube_poly poly(coeff);

        return poly;
    }

    //Переопределение умножения для кубических полиномов
    Cube_poly operator*(Cube_poly b)
    {
        double coeff[4] = { 0, 0, 0, 0 };
        Cube_poly* a = this;

        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j <= i; ++j)
                coeff[i] += a->get_coeff(j) * b.get_coeff(i - j);
        }

        Cube_poly poly(coeff);

        return poly;
    }

    //Переопределение умножения для кубического полинома и числа
    Cube_poly operator*(double b)
    {
        double coeff[4];
        Cube_poly* a = this;

        for (int i = 0; i < 4; ++i)
            coeff[i] = a->get_coeff(i) * b;

        Cube_poly poly(coeff);

        return poly;
    }

    //Получение коэффициента кубического полинома по номеру
    const double get_coeff(int i)
    {
        return coeff[i];
    }

    //Получение начала промежутка задания
    const double get_a()
    {
        return a;
    }

    //Получение конца промежутка задания
    const double get_b()
    {
        return b;
    }
};

//Класс кубических сплайнов
class Cubic_spline
{
private:
    //Кубические полиномы, образующие сплайн
    Cube_poly* poly;
    //Метод прогонки
    void Method_run_through(int n, double* a, double* b, double* c, double* d, double* x)
    {
        double* alpha = new double[n];
        double* beta = new double[n];

        alpha[0] = -c[0] / b[0], beta[0] = d[0] / b[0];

        //Прямой этап метода прогонки
        for (int i = 1; i < n; ++i)
        {
            double gamma = b[i] + a[i] * alpha[i - 1];
            if (i != n - 1)
                alpha[i] = -c[i] / gamma;
            beta[i] = (d[i] - a[i] * beta[i - 1]) / gamma;
        }

        x[n - 1] = beta[n - 1];

        //Обратный этап метода прогонки
        for (int i = n - 2; i >= 0; i--)
            x[i] = alpha[i] * x[i + 1] + beta[i];

        delete[] alpha;
        delete[] beta;
    }

    //Распараллеливание метод прогонки
    double* Method_run_through_parallel(int n, double* a_, double* b_, double* c_, double* f_)
    {
        int m, p;
        double* x = new double[n + 1];
        #pragma omp parallel
        {
            p = omp_get_num_threads();
            m = (n - 1) / p + 1;
        }

        double* a = new double[m * p + 1];
        double* b = new double[m * p + 1];
        double* c = new double[m * p + 1];
        double* f = new double[m * p + 1];

        for (int i = 0; i <= n; ++i)
        {
            a[i] = a_[i];
            b[i] = b_[i];
            c[i] = c_[i];
            f[i] = f_[i];
        }

        for (int i = n + 1; i <= m * p; ++i)
        {
            a[i] = 0;
            b[i] = 1;
            c[i] = 0;
            f[i] = 0;
        }
        
        double* xmu = new double[p + 1];
        double** v = new double*[p];
        double** z = new double* [p];
        double** w = new double* [p];

        for (int i = 0; i < p; ++i)
        {
            v[i] = new double[m + 1];
            z[i] = new double[m + 1];
            w[i] = new double[m + 1];
        }

        double* d;

        if (m > 1)
            d = new double[m - 1];

        //Вычисление вспомогательных коэффициентов v, z и w
        #pragma omp parallel for shared(m, p)
            for (int mu = 0; mu < p; ++mu)
            {
                if (m > 1)
                {
                    for (int i = 1; i < m; ++i)
                        d[i] = 0;
                    d[1] = -a[mu * m + 1];

                    Method_run_through(m - 1, a + mu * m + 1, b + mu * m + 1, c + mu * m + 1, d + 1, v[mu] + 1);
                }
                v[mu][0] = 1;
                v[mu][m] = 0;

                if (m > 1)
                {
                    d[1] = 0;
                    d[m - 1] = -c[mu * m + m - 1];

                    Method_run_through(m - 1, a + mu * m + 1, b + mu * m + 1, c + mu * m + 1, d + 1, z[mu] + 1);
                }
                z[mu][0] = 0;
                z[mu][m] = 1;

                if (m > 1)
                {
                    for (int i = 1; i < m; ++i)
                        d[i] = f[mu * m + i];

                    Method_run_through(m - 1, a + mu * m + 1, b + mu * m + 1, c + mu * m + 1, d + 1, w[mu] + 1);
                }
                w[mu][0] = 0;
                w[mu][m] = 0;
            }

        double* amu = new double[p + 1];
        double* bmu = new double[p + 1];
        double* cmu = new double[p + 1];
        double* fmu = new double[p + 1];

        amu[0] = cmu[p] = 0;

        #pragma omp parallel for shared(m, p)
            for (int mu = 0; mu <= p; ++mu)
            {
                if (mu == 0)
                {
                    amu[mu] = 0;
                    bmu[mu] = b[mu * m] + a[mu * m] * z[mu][m - 1] + c[mu * m] * v[mu][1];
                    cmu[mu] = c[mu * m] * z[mu][1];
                    fmu[mu] = f[mu * m] - c[mu * m] * w[mu][1];
                }
                else if (mu == p)
                {
                    amu[mu] = 0;
                    bmu[mu] = b[mu * m];
                    cmu[mu] = 0;
                    fmu[mu] = f[mu * m] - a[mu * m] * w[mu - 1][m - 1];
                }
                else
                {
                    amu[mu] = a[mu * m] * v[mu][m - 1];
                    bmu[mu] = b[mu * m] + a[mu * m] * z[mu][m - 1] + c[mu * m] * v[mu][1];
                    cmu[mu] = c[mu * m] * z[mu][1];
                    fmu[mu] = f[mu * m] - a[mu * m] * w[mu - 1][m - 1] - c[mu * m] * w[mu][1];
                }
            }

        //Вычисление x_mu - обобщающие неизвестные
        Method_run_through(p + 1, amu, bmu, cmu, fmu, xmu);

        delete[] amu;
        delete[] bmu;
        delete[] cmu;
        delete[] fmu;

        //Вычисление искомых неизвестных
        #pragma omp parallel for shared(m, p)
            for (int mu = 0; mu < p; ++mu)
                for (int j = 0; j <= m; ++j)
                    if (j + mu * m < n + 1)
                        x[j + mu * m] = v[mu][j] * xmu[mu] + z[mu][j] * xmu[mu + 1] + w[mu][j];

        delete[] a;
        delete[] b;
        delete[] c;
        delete[] f;
        delete[] xmu;
        delete[] v;
        delete[] z;
        delete[] w;

        return x;
    }
public:
    //Кубический сплайн по таблице аргумент-значение
    //n - количетво промежутков разбиения
    //f0 - значение производной функции в первом узле
    //fn - значение производной функции в последнем узле
    //argValue - таблица аргумент-значение
    //parallel - флаг, определяющий применение OpenMP
    Cubic_spline(int n, double f0, double fn, My_pair* argValue, bool parallel)
    {
        int i;
        poly = new Cube_poly[n + 1];
        //Значения первых производных функции в узлах сетки
        double* m;
        //Длины промежутков разбиения
        double* h = new double[n + 1];
        //Вспомогательные коэффициенты \mu
        double* mu = new double[n + 1];
        //Вспомогательные коэффициенты \lambda
        double* lambda = new double[n + 1];
        //Вспомогательные коэффициенты c
        double* c = new double[n + 1];
        //Вспомогательный массив, заполненый двойками.
        double* g = new double[n + 1];

        //Нахождение длин промежутков
        #pragma omp parallel for shared(argValue, n) private(i) if (parallel)
            for (i = 0; i < n; ++i)
            {
                h[i] = argValue[i + 1].first - argValue[i].first;
                g[i] = 2;
            }

        mu[0] = lambda[n] = 0;
        g[n] = 2;

        //Нахождение вспомогательных коэффициентов \mu и \lambda
        #pragma omp parallel for shared(n) private(i) if (parallel)
            for (i = 1; i < n; ++i)
                lambda[i] = 1 - (mu[i] = h[i - 1] / (h[i - 1] + h[i]));

        c[0] = 2 * f0, c[n] = 2 * fn;

        //Нахождение вспомогательных коэффициентов c
        #pragma omp parallel for shared(n) private(i) if (parallel)
            for (i = 1; i < n; ++i)
                c[i] = 3 * (mu[i] * (argValue[i + 1].second - argValue[i].second) / h[i] +
                            lambda[i] * (argValue[i].second - argValue[i - 1].second) / h[i - 1]);

        //Метод прогонки
        if (parallel)
            m = Method_run_through_parallel(n, lambda, g, mu, c);
        else
        {
            m = new double[n + 1];
            Method_run_through(n + 1, lambda, g, mu, c, m);
        }

        delete[] mu;
        delete[] lambda;

        //Построение сплайна
        #pragma omp parallel for shared(n, argValue) private(i) if (parallel)
            for (i = 0; i < n; ++i)
                poly[i].constructor(argValue[i].second, argValue[i + 1].second, m[i], m[i + 1], h[i],
                                    argValue[i].first, argValue[i + 1].first);

        delete[] h;
        delete[] m;
        delete[] c;
    }

    //Получение i-го кубического полинома сплайна
    const Cube_poly get_spline(int i)
    {
        return poly[i];
    }

    void deleteSpline()
    {
        delete[] poly;
    }
};

class Function
{
public:
    double Value(double arg)
    {
        return 1 / (1 + 25 * arg * arg) + sin(2 * arg);
    }
};

void output_spline(int n, Cubic_spline* spline1, Cubic_spline* spline2)
{
    for (int i = 0; i < n; ++i)
    {
        Cube_poly poly1 = spline1->get_spline(i);
        Cube_poly poly2 = spline2->get_spline(i);

        for (int j = 0; j < 4; ++j)
            cout << poly1.get_coeff(j) << " <-> " << poly2.get_coeff(j) << "\n";

        cout << "[" << poly1.get_a() << ", " << poly1.get_b() << "]\n\n";
    }
}

void test(double a, double b, int n)
{
    cout << n << "\n";

    Function f;
    double x = a;
    double h = (b - a) / n;

    My_pair* argValue = new My_pair[n + 1];

    for (int i = 0; i <= n; ++i)
    {
        argValue[i].first = x;
        argValue[i].second = f.Value(x);
        x += h;
    }

    double t1, t2;

    t1 = omp_get_wtime();
    Cubic_spline* spline1 = new Cubic_spline(n, 0, 0, argValue, true);
    t2 = omp_get_wtime();

    cout << "Time work (parallel): " << (t2 - t1) * 1000 << "ms\n";

    t1 = omp_get_wtime();
    Cubic_spline* spline2 = new Cubic_spline(n, 0, 0, argValue, false);
    t2 = omp_get_wtime();

    cout << "Time work: " << (t2 - t1) * 1000 << "ms\n";

    output_spline(n, spline1, spline2);

    delete[] argValue;
    spline1->deleteSpline();
    spline2->deleteSpline();
}

int main()
{
    int n;
    double a, b;

    cout << "Enter the number of split intervals: ";
    cin >> n;
    cout << "Enter the beginning of the interval: ";
    cin >> a;
    cout << "Enter the end of the interval: ";
    cin >> b;

    test(a, b, n);

    return 0;
}
