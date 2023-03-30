#include <iostream>
#include <vector>
#include <tuple>
#include <span>
#include <algorithm>
#include <numeric>
#include <functional>
#include <string>
#include <format>

#include <math.h> 
#include <conio.h>
#include <ctime>
#include <cassert>
using namespace std;


//-------------------//
int gsd(int a, int b) {
    if (a % b == 0) return b;
    if (b % a == 0) return a;
    return a > b ? gsd(a % b, b) : gsd(a, b % a);
}
//-------------------//

namespace exercise_1
{
    auto extended_gcd(int a, int b)
    {
        if (!a) return make_tuple(b, 0, 1);
        auto [gcd, x, y] = extended_gcd(b % a, a);
        return make_tuple(gcd, (y - (b / a) * x), x);
    }

    //---1. нахождение НОД и его линейного разложения по расширенному алгоритму Эвклида---//
    auto extended_gcd_2(int a, int b)
    {
        int first  = a > b ? a : b;
        int second = a > b ? b : a;

        vector<int> x{ 1, 0 };
        auto f = [&x](int cell) {return x[x.size() - 2] - cell * x[x.size() - 1]; };

        int cell = 0;
        while (a != 0 && b != 0)
        {
            if (a > b)
            {
                cell = a / b;
                a = a % b;
                x.push_back(f(cell));
            }
            else
            {
                cell = b / a;
                b = b % a;
                x.push_back(f(cell));
            }
        }

        int d = a + b;
        int u = x[x.size() - 2];
        int v = (d - first * u) / second;
       
        return make_tuple(first, second, u, v, d);
    }
}


namespace exercise_2
{
    //---2. поиск простых чисел до числа n---//
    auto search_for_prime_numbers(vector<bool> v, int N)
    {
        int i, j;

        for (i = 1; i <= N; i++)
        {
            v[i] = true;
        }

        i = 1; j = 1;
        while ((2 * i * j + i + j) <= N)
        {
            while (j <= (N - i) / (2 * i + 1))
            {
                v[2 * i * j + i + j] = false;
                j++;
            }
            i++;
            j = i;
        }

        cout << "2 ";
        for (i = 1; i <= N; i++)
        {
            if (v[i]) cout << 2 * i + 1 << " ";
        }
        cout << endl;
    }
}


namespace exercise_3
{
    //---3.	каноническое разложение числа на простые множители---//
    auto print_number_factors(int x)
    {
        cout << "Number factors " << x << ":\n";
        int divisor = 2;
        while (x != 1)
        {
            while (x % divisor == 0)
            {
                cout << divisor << endl;
                x /= divisor;
            }
            divisor += 1;
        }
    }
}


namespace exercise_4
{
    //---4.	расчет функции Эйлера для m---//
    int phi(int m)
    {
        int result = m;
        for (int i = 2; i * i <= m; ++i)
        {
            if (m % i == 0)
            {
                while (m % i == 0) m /= i;
                result -= result / i;
            }
        }
        if (m > 1) result -= result / m;
        return result;
    }
}


namespace exercise_5
{
    //---нахождение обратного элемента в Zm---//
    auto invert(int a, int m)
    {
        if (a < 1 or m < 2) return -1;

        int u1 = m;
        int u2 = 0;
        int v1 = a;
        int v2 = 1;

        while (v1 != 0)
        {
            int q  = u1 / v1;
            int t1 = u1 - q * v1;
            int t2 = u2 - q * v2;
            u1 = v1;
            u2 = v2;
            v1 = t1;
            v2 = t2;
        }

        return u1 == 1 ? (u2 + m) % m : -1;
    }
}


namespace exercise_6
{
    //---6. решение сравнений (для простого и составного m)---//
    auto solveCongruences(int A, int B, int M)
    {
        int x = 0, mod = 1;
        int r1 = A * mod, r2 = M, x1 = 1, x2 = 0, r = 1;

        while (r != 0)
        {
            int q = (int)floor((double)r1 / r2);
            int t = x1 - q * x2;
            x1 = x2;
            x2 = t;
            r = r1 - q * r2;
            r1 = r2;
            r2 = r;
        }

        int b = B - A * x;
        if (b % r1 != 0) return "no solutions!"s;

        x += mod * b * x1 / r1;
        mod *= M / r1;

        if (x < 0 || x >= mod) x -= mod * (int)floor((double)x / mod);

        return format("x = {}(mod {})", x, mod);
    }


    //---6. решение сравнений (для простого и составного m)---//
    int modexp(int x, int y, int N)
    {
        if (y == 0) return 1;
        int z = modexp(x, y / 2, N);
        return (y % 2 == 0) ? (z * z) % N : (x * z * z) % N;
    }
    
    auto solve_linear_congruence(int a, int b, int m)
    {
        int m_start = m;
        int g = gsd(a, m);
    
        if (g == 1)
        {
            if (b % g) throw "No solutions!";
            a = a / g, b = b / g, m = m / g;
            return make_tuple(modexp(a, -1, m) * b % m, m, g); 
        }
        else if (g > 1 && b % g == 0)
        {
            int n = 0;
            a = a / g, b = b / g, m = m / g;
            int x = modexp(a, -1, m) * b % m;
            int firt_x = x;
    
            for (size_t i = 0; i < g - 1; i++)
            {
                cout << format("x = {}+{}k\n", x + n, m_start);
                n = n + m;
            }
    
            return make_tuple(x + n, firt_x, g);
        }
        else throw "No solutions!";
    }
    
    auto result_processing(int a, int b, int m)
    {
        cout << format("{}x = {}(mod {})\nanswer:\n", a, b, m);
        
        try
        {
            auto [x, first_x, delitel] = solve_linear_congruence(a, b, m);
            cout << format("x = {}+{}k\n", x, m);
            cout << format("Common decision: x = {}(mod {})\n", first_x, m / delitel);
        }
        catch (const char* arg)
        {
            cout << arg << endl;
            return;
        }
    }
}


namespace exercise_7
{
    //---7.	решение системы сравнений китайским алгоритмом---//
    string chinese_algorithm(vector<int> v_b, vector<int> v_m)
    {
        //--------------------------//
        for (size_t i = 0; i < v_m.size(); i++)
        {
            for (size_t j = i + 1; j < v_m.size(); j++)
            {
                if (gsd(v_m[i], v_m[j]) != 1) return "no solutions!";
            }
        }
        //--------------------------//

        int M  = accumulate(v_m.begin(), v_m.end(), 1, multiplies());
        int x0 = 0;
        vector<int> vM(v_m.size(), 0);
        vector<int> vY(v_m.size(), 0);
        transform(v_m.begin(), v_m.end(), vM.begin(), [&M](auto e) { return M / e; });

        for (size_t i = 0; i < v_m.size(); i++)
        {
            for (size_t j = 1; j < v_m[i]; j++)
            {
                if ((vM[i] * (int)j - 1) % v_m[i] == 0) vY[i] = j;
            }
            x0 += vM[i] * vY[i] * v_b[i];
        }

        while (x0 > M) x0 -= M;

        return format("x = {}(mod {})", x0, M);
    }
}


namespace exercise_8
{
    //---8. нахождение вычета a^k(mod m) для простого и составного m---//
    auto power_residue_2(int a, int k, int m)
    {
        int temp = k % (m - 1);
        temp = pow(a, temp) - (m * (pow(a, temp) / m));
        temp = temp ? temp : 1;
        return format("{}^{}(mod{}) = {}(mod({}) = {}", a, k, m, temp, m, temp % m);
    }

    //-----------------------------------------------------------------//
    //auto power_residue(int a, int k, int p) 
    //{
    //    if (!(k > p - 1))      throw "error: k > p - 1";
    //    if (!(gsd(a, p) == 1)) throw "error: (a, p) = 1";
    //    return k % (p - 1);
    //}
    //
    //auto power_residue(int a, int k, vector<int> v_p) 
    //{
    //    string      result;
    //    vector<int> v_b(v_p.size(), 0);
    //
    //    for (size_t i = 0; i < v_p.size(); i++)
    //    {
    //        try   
    //        { 
    //            int r = power_residue(a, k, v_p[i]);
    //            result += format("[x = {}^{}(mod {})\n", a, r, v_p[i]);
    //            v_b[i]  = pow(a, r); 
    //        }
    //        catch (const char* msg) 
    //        {     
    //            return ""s + msg;     
    //        }
    //    }
    //
    //    return result + "answer: " + exercise_7::chinese_algorithm(v_b, v_p);
    //}
    
    //-----------------------------------------------------------------//
}


namespace exercise_9
{
    //---9. нахождение первообразного корня (образующего элемента) и формирование с его помощью приведенной системы вычетов---//


    //int powmod(int a, int b, int p) 
    //{
    //    int res = 1;
    //    while (b)
    //    {
    //        if (b & 1) res = int(res * 1ll * a % p), --b;
    //        else         a = int(a * 1ll * a % p), b >>= 1;
    //    }
    //    return res;
    //}
    //
    //int generator(int p) 
    //{
    //    vector<int> fact;
    //    int phi = p - 1;
    //    int n   = phi;
    //
    //    for (int i = 2; i * i <= n; ++i)
    //    {
    //        if (n % i == 0) {
    //            fact.push_back(i);
    //            while (n % i == 0) n /= i;
    //        }
    //    }
    //
    //    if (n > 1) fact.push_back(n);
    //
    //    for (int res = 2; res <= p; ++res) 
    //    {
    //        bool ok = true;
    //        for (size_t i = 0; i < fact.size() && ok; ++i)
    //        {
    //            ok &= powmod(res, phi / fact[i], p) != 1;
    //        }
    //        if (ok) return res;
    //    }
    //    return -1;
    //}
}


namespace exercise_11
{
    //---11. нахождение символа Лежандра и символа Якоби---//
    auto jacobi(int a, int n)
    {
        assert(n > a > 0 && n % 2 == 1);
        
        int temp = 1;

        while (a != NULL)
        {
            while (a % 2 == 0)
            {
                a /= 2;
                int r = n % 8;
                temp = r == 3 || r == 5 ? -temp : temp;
            }

            swap(a, n);
            temp = a % 4 == n % 4 == 3 ? -temp : temp;
            a %= n;
        }

        return n == 1 ? temp : 0;
    }

    auto factorize(int n)
    {
        vector<int> factors;

        int p = 2;
        while (true)
        {
            while (n % p == 0 && n > 0) 
            {
                factors.push_back(p);
                n = n / p;
            }
            p += 1;
            if (p > n / p) break;
        }

        if (n > 1) factors.push_back(n);

        return factors;
    }

    int calculateLegendre(int a, int p)
    {
        if      (a >= p || a < 0)        return calculateLegendre(a % p, p);
        else if (a == 0 || a == 1)       return a;
        else if (a == 2)                 return (p % 8 == 1 || p % 8 == 7) ? 1 : -1;
        else if (a == p - 1)             return p % 4 == 1 ? 1 : -1;
        else if (![](int a) {
                    for (size_t i = 2; i < a; i++)
                    {
                        if (!a % i) return false;
                    }
                    return true;
                }(a))  {
                auto factors = factorize(a);
                int product = 1;
                for (auto& i : factors) product = calculateLegendre(i, p);
                return product;
             }
        else                             return (((p - 1) / 2) % 2 == 0 || ((a - 1) / 2) % 2 == 0) ? calculateLegendre(p, a) : -1;
    }

    auto response_processing(int a, int n)
    {
        if ([](int n) {
            for (size_t i = 2; i < n / 2 + 1; i++)
            {
                if (n % i == 0) return false;
            }
            return true;
            }(n)) {
            cout << "Symbol of Legendre: " << calculateLegendre(a, n) << endl;
        }
        else if (n > 0 && n % 2 == 1) cout << "Jacobi symbol: " << jacobi(a, n) << endl;
        else cout << "Symbol of Legendre: " << calculateLegendre(a, n) << endl;
    }
}


namespace exercise_12
{
    //---12. нахождение порядка всех элементов в группе Z +m , Z *m---//
    void finding_order_plus(int m)
    {
        auto x = [=] {
            vector<int> temp;
            for (size_t i = 0; i < m; i++) temp.push_back(i);
            return temp;
        }();
        auto d = [=] {
            vector<int> temp;
            for (size_t i = 1; i < m / 2 + 1; i++)
            {
                if (m % i == 0) temp.push_back(i);
            }
            temp.push_back(m);
            return temp;
        }();
        for (auto& x_elem : x)
        {
            for (auto& d_elem : d)
            {
                if ((x_elem * d_elem) % m == 0) 
                {
                    cout << format("The order of the number {} is {}\n", x_elem, d_elem);
                    break;
                }
            }
        }
    }

    void finding_order_multiply(int m)
    {
        int Z_m = exercise_4::phi(m);
        auto x = [=] {
            vector<int> temp;
            for (size_t i = 1; i < m; i++) temp.push_back(i);
            return temp;
        }();
        auto d = [=] {
            vector<int> temp;
            for (size_t i = 1; i < m / 2 + 1; i++)
            {
                if (m % i == 0) temp.push_back(i);
            }
            temp.push_back(Z_m);
            return temp;
        }();

        for (auto& x_elem : x)
        {
            bool temp = false;
            for (auto& d_elem : d)
            {
                temp = false;
                if ((long long)pow((long long)x_elem, (long long)d_elem) % m == 1)
                {
                    cout << format("The order of the number {} is {}\n", x_elem, d_elem);
                    temp = true;
                    break;
                }
            }

            if (temp == false) cout << "There is no order for the number " << x_elem << endl;
        }
    }
}


namespace exercise_17
{
    //---17 построить псевдослучайную последовательность по одному из алгоритмов генерации---//
    class My_random
    {
        using uint64 = unsigned long long;
    public:
        My_random(uint64 seed = 0) : X(seed) {};

        auto& operator = (uint64 seed) 
        { 
            X = seed; 
            return *this; 
        }

        uint64 operator()(uint64 seed = uint64(-1))
        {
            const uint64 a = 3202034522624059733ULL;
            const uint64 c = 1ULL;

            if (seed != uint64(-1)) X = seed;
            uint64 Y = a * X + c;
            X = a * Y + c;
            Y = (Y & 0xFFFFFFFF00000000ULL) | (X >> 32);
            return Y;
        }

        uint64 operator()(uint64 min, uint64 max)
        {
            return (*this)() % (max - min) + min;
        }

    private:
        uint64 X;
    };
}


int main()
{
    string separator = string(30, '-') + "\n\n";
    //---1---//
    {
        using namespace exercise_1;
        cout << "N1:\n";

        auto [first, second, u, v, d] = extended_gcd_2(549, 387);
        cout << format("({}, {}) = {}*{} + {}*{} = {}", first, second, first, u, second, v, d) << endl;
        
        cout << separator;
    }
    //---2---//
    {
        using namespace exercise_2;
        cout << "N2:\n";

        int N = 1000;
        vector<bool> massiv(N);

        cout << "input: " << N << endl;
        search_for_prime_numbers(massiv, N / 2 - 1);

        cout << separator;
    }
    //---3---//
    {
        using namespace exercise_3;
        cout << "N3:\n";

        print_number_factors(2312);

        cout << separator;
    }
    //---4---// +
    {
        using namespace exercise_4;
        cout << "N4:\n";

        int input = 375;
        cout << "phi(" << input << ") = " << phi(input) << endl;
        
        cout << separator;
    }
    //---5---//
    {
        using namespace exercise_5;
        cout << "N5:\n";

        int a = 5, m = 11;
        cout << format("input(a, m): {}, {}\n", a, m);
        cout << "result: " << invert(5, 11) << endl;

        cout << separator;
    }
    //---6---//
    {
        using namespace exercise_6;
        cout << "N6:\n";

        int a = 4, b = 6, m = 2;
        result_processing(a, b, m);

        cout << separator;
    }
    //---7---// + 
    {
        using namespace exercise_7;
        cout << "N7:\n";

        vector<int> v_b{ 3, 2 };
        vector<int> v_m{ 8, 16 };

        for (size_t i = 0; i < v_m.size(); i++)
        {
            cout << "[x = " << v_b[i] << "(mod " << v_m[i] << ")\n";
        }
        cout << "answer: " << chinese_algorithm(v_b, v_m) << endl;

        cout << separator;
    }
    //---8---//
    {
        using namespace exercise_8;
        cout << "N8:\n";

        int a = 6, k = 12, m = 13;

        cout << format("a = {}\nb = {}\nm = {}\n", a, k, m);
        cout << power_residue_2(a, k, m) << endl;

        cout << separator;
    }
    //---11---//
    {
        //13, 15 якоби

        using namespace exercise_11;
        cout << "N11:\n";

        int a = 13, b = 15;

        cout << format("a = {}\nb = {}\n", a, b);
        response_processing(a, b);

        cout << separator;
    }
    //---12---//
    {
        using namespace exercise_12;
        cout << "N12:\n";

        int m = 12;
        cout << "Finding the order of all elements in a group Zm+, m = " << m << "\nresult Zm+:\n";
        finding_order_plus(m);
        cout << "Finding the order of all elements in a group Zm*, m = " << m << "\nresult Zm*:\n";
        finding_order_multiply(m);

        cout << separator;
    }
    //---17---//
    {
        using namespace exercise_17;
        cout << "N17:\n";

        int N = 10;
        cout << "count: " << N << endl;

        auto tm = std::tm{};
        auto time = static_cast<unsigned long long>(std::mktime(&tm));
        My_random rand(time);

        for (size_t i = 0; i < 10; i++)
        {
            cout << rand() << endl;
        }
    }
}












////Алгоритм Евклида
//int node(int a, int b) {
//    while (a && b) a > b ? a %= b : b %= a;
//    return (a + b);
//}
//
////Расширенный алгоритм Евклида 
//auto extended_gcd(int a, int b) {
//    if (!a) return make_tuple(b, 0, 1);
//    auto [gcd, x, y] = extended_gcd(b % a, a);
//    return make_tuple(gcd, (y - (b / a) * x), x);
//}
//
////Разложение на простые множители
//void prime_factorization(int n) {
//    int temp = 2;
//    while (n != 1) {
//        while (!(n % temp)) {
//            cout << temp << ' ';
//            n /= temp;
//        }
//        ++temp;
//    }
//    cout << endl;
//}
//
////Функция Эйлера для m
//int phi(int m) {
//    int result = m;
//    for (int i = 2; i * i <= m; ++i)
//        if (m % i == 0) {
//            while (m % i == 0) m /= i;
//            result -= result / i;
//        }
//    if (m > 1) result -= result / m;
//    return result;
//}
//
//
//
//int main() {
//    int a = 49104, b = 624960;
//
//    //1)==================================
//    //cout << node(49104, 624960) << endl;
//    //2)======================================================================
//    auto [gcd, x, y] = extended_gcd(22, 6);
//    cout << a << "*" << x << " + " << b << "*" << y << " = " << gcd << endl;
//
//    //3)?????????????????????
//    //4)=====================
//    //prime_factorization(a);
//
//    //5)=====================
//    //cout << phi(7) << endl;
//
//    //6)---------------------------------
//    //7)
//
//}