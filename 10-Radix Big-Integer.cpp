// g++ -std=c++11 -O2
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
typedef long long ll;
typedef unsigned long long ull;
typedef __int128_t lll;
namespace FFT_Solver
{
    struct complex
    {
        double a, b; // a : real number b : imaginary number
        double len() const { return a * a + b * b; }
        complex operator+(const complex &o) const { return {a + o.a, b + o.b}; }
        complex operator-(const complex &o) const { return {a - o.a, b - o.b}; }
        complex operator-() const { return {-a, -b}; }
        complex operator*(const double &o) const { return {a * o, b * o}; }
        complex operator*(const complex &o) const { return {a * o.a - b * o.b, b * o.a + a * o.b}; }
        complex operator/(const double &o) const { return {a / o, b / o}; }
        complex operator!() const { return {a, -b}; } // conjugate
        complex operator/(const complex &o) const { return ((*this) * (!o)) / o.len(); }
    };
    const int N = (1 << 16) | 3;
    const double PI = acos(-1.0);
    bool initialized;
    int L, brev[N]; // Butterfly operation
    complex w[N], v[N];
    complex com_a[N], com_b[N];
    void init(int _L)
    {
        L = _L, initialized = 1;
        for (int i = 0; i < (1 << L); ++i)
            brev[i] = (brev[i >> 1] >> 1) | ((i & 1) << (L - 1));
        for (int i = 0; i < (1 << L); ++i)
        {
            w[i] = {cos(i * PI * 2 / (1 << L)), sin(i * PI * 2 / (1 << L))};
            v[i] = {cos(i * PI * 2 / (1 << L)), -sin(i * PI * 2 / (1 << L))};
        }
    }
 
    struct initializer
    {
        // length is adjustable
        initializer() { init(16); }
    }fft_init;
 
    void fft(complex a[], int lgn, int flag)
    {
        int n = 1 << lgn;
        for (int i = 0; i < n; ++i)
        {
            int rv = brev[i] >> (L - lgn);
            if (rv < i)
            {
                // std::swap(a[rv], a[i]);
                complex tmp = a[rv];
                a[rv] = a[i], a[i] = tmp;
            }
        }
 
        int fa = L;
        complex *q = (flag == 1) ? w : v;
 
        for (int t = 1; t < n; t <<= 1)
        {
            --fa;
 
            for (int i = 0; i < n; i += t << 1)
            {
                complex *p = a + i;
 
                for (int j = 0; j < t; ++j)
                {
                    complex x = p[j + t] * q[j << fa];
                    p[j + t] = p[j] - x, p[j] = p[j] + x;
                }
            }
        }
 
        if (flag == -1)
            for (int i = 0; i < n; ++i)
                a[i] = {a[i].a / n, a[i].b / n};
    }
 
    void fft_two_seq(complex a[], complex b[], int lgn, int flag)
    {
        int n = 1 << lgn;
 
        for (int i = 0; i < n; ++i)
            a[i].b = b[i].a;
 
        fft(a, lgn, flag);
 
        b[0] = !a[0];
 
        for (int i = 1; i < n; ++i)
            b[i] = !a[n - i];
 
        for (int i = 0; i < n; ++i)
        {
            complex x = a[i], y = b[i];
            a[i] = {(x.a + y.a) / 2.0, (x.b + y.b) / 2.0};
            b[i] = (x - y) / (complex){0, 2};
        }
    }
 
    // a[0...n] * b[0...m] (assume that res is all zero)
    void mul(int a[], int n, int b[], int m, int res[], int mod, int *res_len)
    {
        // multiple case
        assert((*res_len) >= n + m + 1);
        memset(res, 0, sizeof(res[0]) * (*res_len));

        // brute force
        if (n < 100 / (m + 1) || n < 3 || m < 3)
            for (int i = 0; i <= n; ++i)
                for (int j = 0; j <= m; ++j)
                {
                    // res[i + j] += a[i] * b[j];
                    ll x = 1ll * res[i + j] + (1ll * a[i] * b[j]);
                    res[i + j] = 0;
                    int offset = 0;
                    while (x >= mod)
                        res[i + j + offset] += (x % mod), x /= mod, ++offset;
                    res[i + j + offset] += x;                   
                }   
        // FFT
        else
        {
            assert(initialized);
            int lgk = 0, k = 1, len = n + m;
            while ((1 << lgk) <= len)
                ++lgk, k <<= 1;
            for (int i = 0; i <= n; ++i)
                com_a[i].a = a[i], com_a[i].b = 0.0;
            for (int i = 0; i <= m; ++i)
                com_b[i].a = b[i], com_b[i].b = 0.0;
            // multiple_case
            memset(com_a + (n + 1), 0, sizeof(com_a[0]) * (k - n - 1));
            memset(com_b + (m + 1), 0, sizeof(com_b[0]) * (k - m - 1));
 
            fft_two_seq(com_a, com_b, lgk, 1);
            for (int i = 0; i < k; ++i)
                com_a[i] = com_a[i] * com_b[i];
            fft(com_a, lgk, -1);
            for (int i = 0; i <= n + m; ++i)
            {
                ll x = 1ll * res[i] + (ll)(com_a[i].a + 0.5);
                res[i] = 0;
                int offset = 0;
                while (x >= mod)
                    res[i + offset] += (x % mod), x /= mod, ++offset;
                res[i + offset] += x;
                // res[i] = (ll)(com_a[i].a + 0.5);
            }
        }
        // adjust by mod
        (*res_len) = n + m + 1;
        while (res[(*res_len)])
            ++(*res_len);
        while (res[(*res_len) - 1] >= mod)
            res[(*res_len)] = res[(*res_len) - 1] / mod, res[(*res_len) - 1] %= mod, ++(*res_len);
        while ((*res_len) && (!res[(*res_len) - 1]))
            --(*res_len);
    }
}


// base radix : 10
struct BigNumber
{
    const int digit[6] = {1, 10, 100, 1000, 10000, 100000}, MODD = 5;
    static const int MOD = 100000, base_radix = 10;

    int len;
    bool f;
    int *a;
    BigNumber() { a = NULL; }
    BigNumber(lll x)
    {
        a = new int[12];
        len = 0, f = true;
        if (x < 0)
            x = -x, f = false;
        int sum = 0, k = 0;
        while (x)
        {
            sum += x % base_radix * digit[k++];
            if (k == MODD)
                a[len++] = sum, sum = 0, k = 0;
            x /= base_radix;
        }
        if (sum > 0)
            a[len++] = sum;
        if (len == 0)
            len = 1, a[0] = 0;
    }
    BigNumber(const char *s)
    {
        int slen = strlen(s);
        a = new int[slen / MODD + 2];
        len = 0, f = true;
        int start = 0;
        if (s[0] == '+')
            start = 1;
        else if (s[0] == '-')
            start = 1, f = false;
        int sum = 0, k = 0;
        for (int i = slen - 1; i >= start; i--)
        {
            sum += (s[i] ^ '0') * digit[k++];
            if (k == MODD)
                a[len++] = sum, sum = 0, k = 0;
        }
        if (sum > 0)
            a[len++] = sum;
        while (len && (!a[len - 1]))
            --len;
        if (len == 0)
            len = 1, a[0] = 0, f = true;
    }
    BigNumber(const BigNumber &p)
    {
        len = p.len, f = p.f, a = new int[len];
        memcpy(a, p.a, len * sizeof(int));
    }
    void input()
    {
        static char s[100010]; // the length of the input string is adjustable
        scanf("%s", s);
        // printf("len is %d", strlen(s));
        (*this) = BigNumber(s);
    }
    bool is_zero() const { return (len == 1) && (a[0] == 0); }
    friend BigNumber abs(const BigNumber &A)
    {
        BigNumber B = A;
        B.f = true;
        return B;
    }
    int print_len() const
    {
        if (is_zero())
            return 1;
        int ret = (this->len - 1) * MODD;
        int tmp = this->a[this->len - 1], add = 0;
        while (digit[add] <= tmp && add < MODD)
            ++add;
        return ret + add;
    }
    void operator=(const BigNumber &p)
    {
        if (a != NULL)
            delete a;
        len = p.len, f = p.f;
        a = new int[len];
        memcpy(a, p.a, len * sizeof(int));
    }
    void trim()
    {
        while (len > 1 && !a[len - 1])
            --len;
    }
    int compare(const BigNumber &p) const
    {
        if (len == p.len)
        {
            for (int i = len - 1; i >= 0; i--)
                if (a[i] != p.a[i])
                    return (a[i] > p.a[i]) * 2 - 1;
            return 0;
        }
        return (len > p.len) * 2 - 1;
    }
    bool operator==(const BigNumber &p) const { return !compare(p); }
    bool operator!=(const BigNumber &p) const { return compare(p); }
    bool operator<(const BigNumber &p) const { return compare(p) < 0; }
    bool operator>(const BigNumber &p) const { return compare(p) > 0; }
    bool operator<=(const BigNumber &p) const { return compare(p) <= 0; }
    bool operator>=(const BigNumber &p) const { return compare(p) >= 0; }
    BigNumber operator-(const BigNumber &p) const
    {
        BigNumber tmp = p, tmp2 = (*this);
        tmp.f = !(tmp.f);
        return tmp + tmp2;
    }
    BigNumber operator+(const BigNumber &p) const
    {
        BigNumber tmp;
        BigNumber tmp2 = (*this);
        tmp.len = (tmp2.len > p.len) ? (tmp2.len + 1) : (p.len + 1);
        tmp.a = new int[tmp.len];
        if (tmp2.f == p.f)
        {
            tmp.a[0] = 0;
            for (int i = 0; i < tmp.len - 1; i++)
            {
                ll x = 0;
                x += i < tmp2.len ? tmp2.a[i] : 0, x += i < p.len ? p.a[i] : 0, x += tmp.a[i];
                tmp.a[i] = x % MOD, tmp.a[i + 1] = x / MOD;
            }
            tmp.f = f;
        }
        else
        {
            int com = compare(p);
            if (com == 0)
            {
                delete tmp.a;
                tmp.a = new int[1], tmp.f = true;
                tmp.a[0] = 0, tmp.len = 1;
            }
            else
            {
                BigNumber p1 = p;
                BigNumber *t1 = &tmp2, *t2 = &p1;
                if (com == -1)
                {
                    BigNumber *tmp = t1;
                    t1 = t2, t2 = tmp;
                }
                ll x = 0;
                for (int i = 0; i < t1->len; i++)
                {
                    bool y = false;
                    int t2num = i < t2->len ? t2->a[i] : 0;
                    if (x + t1->a[i] < t2num)
                        y = true, x += MOD;
                    tmp.a[i] = t1->a[i] + x - t2num;
                    x = y ? -1 : 0;
                }
                tmp.len = t1->len;
                tmp.f = t1->f;
            }
        }
        while (tmp.len > 1 && tmp.a[tmp.len - 1] == 0)
            tmp.len--;
        return tmp;
    }

    BigNumber operator*(const ll& _p) const
    {
        ll p = _p;
        if ((this->len == 1 && this->a[0] == 0) || !p)
            return BigNumber((lll)0);
        BigNumber tmp;
        tmp.f = !(f ^ (p > 0));
        if (p < 0)
            p = -p;
        tmp.len = 0;
        tmp.a = new int[this->len + 10];
        lll r = 0;
        for (int i = 0; i < len; ++i)
        {
            r += this->a[i] * p;
            tmp.a[tmp.len++] = r % MOD;
            r /= MOD;
        }
        while (r)
            tmp.a[tmp.len++] = r % MOD, r /= MOD;
        return tmp;
    }

    BigNumber operator*(const BigNumber &p) const
    {
        if ((this->len == 1 && this->a[0] == 0) || (p.len == 1 && p.a[0] == 0))
            return BigNumber((lll)0);
        BigNumber ret;
        ret.f = !(f ^ p.f);
        ret.len = len + p.len + 5;
        ret.a = new int[ret.len];
        FFT_Solver::mul(a, len - 1, p.a, p.len - 1, ret.a, MOD, &ret.len);
        return ret;
    }

    BigNumber operator/(const ll& _p) const
    {
        ll p = _p;
        assert(p);
        BigNumber tmp;
        lll r = 0;
        tmp.f = !(f ^ (p > 0)), tmp.len = len;
        tmp.a = new int[tmp.len];
        if (p < 0)
            p = -p;
        memset(tmp.a, 0, sizeof(int) * tmp.len);
        for (int i = len - 1; ~i; --i)
            r = r * MOD + a[i], tmp.a[i] = r / p, r %= p;
        while (tmp.len > 1 && tmp.a[tmp.len - 1] == 0)
            --tmp.len;
        return tmp;
    }

    ll operator%(const ll& _p) const
    {
        ll p = _p;
        assert(p);
        int is_pos = !(f ^ (p > 0));
        if (p < 0)
            p = -p;
        lll r = 0;
        for (int i = len - 1; ~i; --i)
            r = r * MOD + a[i], r %= p;
        return is_pos ? r : -r;
    }

    // (*this) * (MOD)^n
    BigNumber left_shift(int n) const
    {
        if (!n)
    		return (*this);
        BigNumber ret;
        ret.a = new int[len + n + 2];
        ret.len = len + n, ret.f = f;
        memset(ret.a, 0, sizeof(int) * n);
        memcpy(ret.a + n, a, sizeof(int) * len);
        return ret;
    }

    // (*this) / (mod)^n
    BigNumber right_shift(int n) const
    {
        if (!n)
    		return (*this);
        if (len <= n)
            return BigNumber((lll)0);
        BigNumber ret;
        ret.a = new int[len - n + 2];
        ret.len = len - n, ret.f = f;
        memcpy(ret.a, a + n, sizeof(int) * (ret.len));
        return ret;
    }
    // return a[len - n,.....,len - 1]
    BigNumber ahead_digit(int n) const
    {
        assert(n <= len);
        BigNumber ret;
        ret.f = f, ret.len = n;
        ret.a = new int[n + 2];
        memcpy(ret.a, a + len - n, sizeof(int) * n);
        return ret;
    }
    friend BigNumber inv(const BigNumber &A)
    {
        if (A.len == 1)
            return BigNumber(((lll)MOD * MOD) / (lll)A.a[0]);
        else if (A.len == 2)
            return BigNumber(((lll)MOD * (lll)MOD * (lll)MOD * (lll)MOD) / ((lll)A.a[0] + (lll)A.a[1] * MOD));
        else
        {
            int k = (A.len >> 1) + 1;
            BigNumber B = inv(A.ahead_digit(k));
            return ((B + B).left_shift(A.len - k)) - ((A * B) * B).right_shift(k << 1);
        }
    }
    // faster mode : newton raphson
    // abs round down (e.g. -7 / 3 = -2)

    BigNumber operator/(const BigNumber &o) const
    {
        assert(o.len > 1 || (o.len == 1 && o.a[0] != 0));
        bool flag = !(f ^ o.f);
        BigNumber A = abs((*this)), B = abs(o);

        if (A < B)
            return BigNumber((lll)0);
        if (A == B)
        {
            BigNumber ret = BigNumber((lll)1);
            ret.f = flag;
            return ret;
        }
        int n = A.len, m = B.len;
        if (n > (m << 1))
        {
            A = A.left_shift(n - (m << 1));
            B = B.left_shift(n - (m << 1));
            m = n - m, n = m << 1;
        }
        BigNumber D = inv(B), ans = BigNumber((lll)0), one = BigNumber((lll)1);
        BigNumber DB = D * B;
        if (DB.len >= (m << 1))
        {
            DB = DB - B;
            if (!D.a[0])
                D = D - one;
            else
                --D.a[0];
        }
        while (B <= A)
        {
            BigNumber C = (A * D).right_shift(m << 1);
            if (C.len == 1 && C.a[0] == 0)
                break;
            ans = ans + C;
            A = A - (B * C);
        }
        while (B <= A)
        {
            A = A - B;
            if (ans.a[0] == MOD - 1)
                ans = ans + one;
            else
                ++ans.a[0];
        }
        ans.f = flag;
        return ans;
    }
    // reminder
    // e.g. -7 % 3 = -1

    BigNumber operator%(const BigNumber &o) const
    {
        assert(o.len > 1 || (o.len == 1 && o.a[0] != 0));
        bool flag = !(f ^ o.f);
        BigNumber A = abs((*this)), B = abs(o);

        if (A < B)
        {
            A.f = flag;
            return A;
        }
        if (A == B)
            return BigNumber((lll)0);
        int n = A.len, m = B.len, right = 0;
        if (n > (m << 1))
        {
            A = A.left_shift(n - (m << 1));
            B = B.left_shift(n - (m << 1));
            right = n - (m << 1);
            m = n - m, n = m << 1;
        }
        BigNumber D = inv(B), one = BigNumber((lll)1);
        BigNumber DB = D * B;
        if (DB.len >= (m << 1))
        {
            DB = DB - B;
            if (!D.a[0])
                D = D - one;
            else
                --D.a[0];
        }
        while (B <= A)
        {
            BigNumber C = (A * D).right_shift(m << 1);
            if (C.len == 1 && C.a[0] == 0)
                break;
            A = A - (B * C);
        }
        while (B <= A)
            A = A - B;
        A.f = flag;
        return right ? A.right_shift(right) : A;
    }

    // update answer to ret_div and ret_rem
    void div_rem(const BigNumber& o, BigNumber& ret_div, BigNumber& ret_rem) 
    {
        assert(o.len > 1 || (o.len == 1 && o.a[0] != 0));
        bool flag = !(f ^ o.f);
        BigNumber A = abs((*this)), B = abs(o);

        if (A < B)
        {
            A.f = flag;
            ret_div = BigNumber((lll)0), ret_rem = A;
            return;
        }
        if (A == B)
        {
            BigNumber ret = BigNumber((lll)1);
            ret.f = flag;
            ret_div = ret, ret_rem = BigNumber((lll)0);
            return;
        }
        int n = A.len, m = B.len, right = 0;
        if (n > (m << 1))
        {
            A = A.left_shift(n - (m << 1));
            B = B.left_shift(n - (m << 1));
            right = n - (m << 1);
            m = n - m, n = m << 1;
        }
        BigNumber D = inv(B), ans = BigNumber((lll)0), one = BigNumber((lll)1);
        BigNumber DB = D * B;
        if (DB.len >= (m << 1))
        {
            DB = DB - B;
            if (!D.a[0])
                D = D - one;
            else
                --D.a[0];
        }
        while (B <= A)
        {
            BigNumber C = (A * D).right_shift(m << 1);
            if (C.len == 1 && C.a[0] == 0)
                break;
            ans = ans + C;
            A = A - (B * C);
        }
        while (B <= A)
        {
            A = A - B;
            if (ans.a[0] == MOD - 1)
                ans = ans + one;
            else
                ++ans.a[0];
        }
        ans.f = A.f = flag;
        ret_div = ans, ret_rem = right ? A.right_shift(right) : A;
        return;
    }

    BigNumber operator^(int p) const
    {
        assert(p >= 0);
        BigNumber tmp = *this, ans = BigNumber((lll)1);
        while (p)
        {
            if (p & 1)
                ans = ans * tmp;
            p >>= 1;
            tmp = tmp * tmp;
        }
        return ans;
    }

    BigNumber operator^(BigNumber p) const
    {
        assert(p.f);
        BigNumber tmp = *this, ans = BigNumber((lll)1);
        while (!(p.len == 1 && p.a[0] == 0))
        {
            if (p.a[0] & 1)
                ans = ans * tmp;
            p = p / 2ll;
            tmp = tmp * tmp;
        }
        return ans;
    }
    
    // \lfloor log_k (*this) \rfloor
    int log(int k) const
    {
        assert(k > 1);
        int idx = MODD, A = digit[idx], B = print_len() - 1;
        while (!(a[len - 1] / A))
            A = digit[--idx];
        A = a[len - 1] / A; // (*this) = A * 10^B
        int tmp = (log10(A) + B) / log10(k);
        BigNumber target = BigNumber((lll)k) ^ (tmp + 1);
        return (target > (*this)) ? tmp : (tmp + 1);
    }

    BigNumber root(int k) const
    {
        assert(!((!this->f) && (!(k & 1))));
        if (this->is_zero())
            return BigNumber((lll)0);
        else if (k == 1)
            return (*this);
        BigNumber abs_x = abs(*this), one = BigNumber((lll)1);
        // original method : using bit_length, but slower.
        // int bit_length = abs_x.log(2) + 1;
        // BigNumber x = BigNumber((lll)2) ^ (bit_length / k + 1);
        
        // optimize : binary search for top 1 or 2 element
        int len_ans = (print_len() + k - 1) / k;
        int target_top = MOD - 1, tail = (len_ans - 1) / MODD; // target_top * (MOD^(tail))
        int lo = 0, hi = MOD - 1, mi = 0;
        BigNumber x = BigNumber((lll)0), target = (BigNumber((lll)target_top) ^ k).left_shift(tail * k);;
        // 1 element
        if (target < (*this))
            tail++, lo = 1, x = (BigNumber((lll)lo)).left_shift(tail);
        else
        {
            while (lo < hi)
            {
                mi = (lo + hi) >> 1;
                target = (BigNumber((lll)mi) ^ k).left_shift(tail * k);
                if (target <= (*this))
                    lo = mi + 1;
                else
                    hi = mi;
            }
            if (!tail)
                x = (BigNumber((lll)lo)).left_shift(tail);
            else
            {
                // 2 element
                lll LO = 1ll * (lo - 1) * MOD, HI = 1ll * lo * MOD, MI = 0;
                --tail;
                while (LO < HI)
                {
                    MI = (LO + HI) >> 1;
                    target = (BigNumber((lll)MI) ^ k).left_shift(tail * k);
                    if (target <= (*this))
                        LO = MI + 1;
                    else
                        HI = MI;
                }
                x = (BigNumber((lll)LO)).left_shift(tail);
            }
        }
        // newton raphson from a proper initial value
        int k1 = k - 1;
        while (1)
        {
            BigNumber xx = abs_x / (x ^ k1);
			if (x <= xx)
                break;
            x = ((x * k1) + xx) / k;
        }
        x.f = this->f;
        return x;
    }

    BigNumber& operator+=(const BigNumber& o) { return (*this) = (*this) + o, (*this); }
    BigNumber& operator-=(const BigNumber& o) { return (*this) = (*this) - o, (*this); }
    BigNumber& operator*=(const BigNumber& o) { return (*this) = (*this) * o, (*this); }
    BigNumber& operator/=(const BigNumber& o) { return (*this) = (*this) / o, (*this); }
    BigNumber& operator%=(const BigNumber& o) { return (*this) = (*this) % o, (*this); }
    // pre-autoincrement
    
    BigNumber& operator++()
    {
        if (a[0] < MOD - 1)
            a[0]++;
        else
        {
            BigNumber one = BigNumber("1");
            (*this) += one;
        }
        return (*this);
    }

    BigNumber operator++(int)
    {
        BigNumber ret = (*this);
        if (a[0] < MOD - 1)
            a[0]++;
        else
        {
            BigNumber one = BigNumber("1");
            (*this) += one;
        }
        return ret;
    }

    BigNumber& operator--()
    {
        if (a[0] > 0)
            a[0]--;
        else
        {
            BigNumber one = BigNumber("1");
            (*this) -= one;
        }
        return (*this);
    }

    BigNumber operator--(int)
    {
        BigNumber ret = (*this);
        if (a[0] > 0)
            a[0]--;
        else
        {
            BigNumber one = BigNumber("1");
            (*this) -= one;
        }
        return (*this);
    }
    

    // bignum37 : optimize for gcd(a, b)

    struct bignum37
    {
        static const ll _base = 1000000000000000000ll;
        static const lll base = (lll)_base * (lll)_base * (lll)10;
        static const int bit_width = 37, base_radix = 10;
        lll *a;
        int len;
        bignum37() {}
        bignum37(const char *s)
        {
            int n = strlen(s), i = 0;
            a = new lll[n / bit_width + 5], len = 0;
            while (1)
            {
                a[len++] = 0;
                if (n - i <= bit_width)
                {
                    for (int j = 0; i + j < n; a[len - 1] = (a[len - 1] << 3) + (a[len - 1] << 1) + (s[j] ^ 48), ++j)
                        ;

                    break;
                }
                i += bit_width;
                for (int j = 0; j < bit_width; a[len - 1] = (a[len - 1] << 3) + (a[len - 1] << 1) + (s[n - i + j] ^ 48), ++j)
                    ;
            }
        }
        inline BigNumber print()
        {
            char *ch = new char[20500]; // the length is adjustable
            int print_len = 0, is_neg = 0;
            lll x = a[len - 1];

            while (x)
                ch[print_len++] = (x % base_radix) ^ '0', x /= base_radix;

            for (int i = is_neg; i < print_len + is_neg - 1 - i; i++)
                // std::swap(ch[i], ch[print_len + is_neg - 1 - i]);
                ch[i] ^= ch[print_len + is_neg - 1 - i], ch[print_len + is_neg - 1 - i] ^= ch[i], ch[i] ^= ch[print_len + is_neg - 1 - i];

            for (int i = len - 2; i >= 0; i--)
            {
                x = a[i];
                int k = bit_width - 1;

                while (x)
                    ch[print_len + k] = (x % base_radix) ^ '0', x /= base_radix, k--;

                while (k >= 0)
                    ch[print_len + k] = '0', k--;

                print_len += bit_width;
            }

            if (print_len == 0)
                ch[print_len++] = '0';

            ch[print_len] = '\0';
            // puts(ch);
            BigNumber ret = BigNumber(ch);
            delete ch;
            return ret;
        }
        inline int cmp(const bignum37 &b) const
        {
            if (len ^ b.len)
                return len < b.len ? -1 : 1;

            for (int i = len - 1; i >= 0; --i)
                if (a[i] ^ b.a[i])
                    return a[i] < b.a[i] ? -1 : 1;

            return 0;
        }
        inline void right_shift()
        {
            for (int i = len - 1; i >= 0; (a[i] & 1) && (a[i - 1] += base), a[i--] >>= 1);

            while (!a[len - 1])
                --len;
        }
        inline void left_shift()
        {
            for (int i = 0; i < len; a[i++] <<= 1);

            for (int i = 0; i < len - 1; a[i + 1] += a[i] / base, a[i] %= base, ++i);

            if (a[len - 1] >= base)
                a[len] = a[len - 1] / base, a[len - 1] %= base, ++len;
        }
        inline void operator-=(const bignum37 &b)
        {
            for (int i = 0; i < b.len; ++i)
            {
                a[i] -= b.a[i];

                for (; a[i] < (lll)0; a[i] += base, --a[i + 1]);
            }

            while (!a[len - 1])
                --len;
        }
        ~bignum37() { delete a; }
    };

    bignum37 trans_to_37() const
    {
        char *ch = new char[1000500];
        int print_len = 0, is_neg = 0;
        if (!f && (len > 1 || a[0] != 0))
            ch[print_len++] = '-', is_neg = 1;
        int x = a[len - 1];
        while (x)
            ch[print_len++] = x % base_radix + '0', x /= base_radix;
        for (int i = is_neg; i < print_len + is_neg - 1 - i; i++)
            ch[i] ^= ch[print_len + is_neg - 1 - i], ch[print_len + is_neg - 1 - i] ^= ch[i], ch[i] ^= ch[print_len + is_neg - 1 - i];
        for (int i = len - 2; i >= 0; i--)
        {
            x = a[i];
            int k = MODD - 1;
            while (x)
                ch[print_len + k] = x % base_radix + '0', x /= base_radix, k--;
            while (k >= 0)
                ch[print_len + k] = '0', k--;
            print_len += MODD;
        }
        if (print_len == 0)
            ch[print_len++] = '0';
        ch[print_len] = '\0';
        // puts(ch);
        bignum37 ret = bignum37(ch);
        delete ch;
        return ret;
    }

    friend BigNumber gcd_BF(BigNumber A, BigNumber B)
    {
        A.f = B.f = 1;
        if (A.is_zero())
            return B;
        else if (B.is_zero())
            return A;
        BigNumber* p1 = &A, *p2 = &B;
        while (!p2->is_zero())
        {
            (*p1) %= (*p2);
            BigNumber *tmp = p1;
            p1 = p2, p2 = tmp;
        }
        return (*p1);
    }

    friend BigNumber gcd(BigNumber A, BigNumber B)
    {
        A.f = B.f = 1;
        if (A.is_zero())
            return B;
        else if (B.is_zero())
            return A;
        bignum37 a = A.trans_to_37(), b = B.trans_to_37();
        int cnt = 0, tmp = 0;
        while (!(a.a[0] & 1) && !(b.a[0] & 1))
            a.right_shift(), b.right_shift(), ++cnt;
        while (1)
        {
            while (!(a.a[0] & 1))
                a.right_shift();

            while (!(b.a[0] & 1))
                b.right_shift();

            tmp = a.cmp(b);

            if (!tmp)
                break;

            if (tmp < 0)
                b -= a;
            else
                a -= b;
        }

        while (cnt--)
            a.left_shift();
        return a.print();
    }

    // (*this) = factorial(x)
    void factorial(int x)
    {
        assert(x >= 0);
        if (!x)
        {
            (*this) = BigNumber((lll)1);
            return;
        }

        int stack_sz = (int)(ceil(log2(x))) + 1; // level of the calling tree
        if (stack_sz <= 1)
            stack_sz = 2;
        BigNumber* s_big = new BigNumber[stack_sz];
        int* s_num = new int[stack_sz];
        int top = 0;
        for (int i = 1; i <= x; ++i)
        {
            BigNumber b = BigNumber((lll)i);
            int level = 1;
            while (top && s_num[top] == level)
                b = b * s_big[top--], ++level;
            s_big[++top] = b, s_num[top] = level;
        }
        (*this) = s_big[top];
        for (int i = top - 1; i; --i)
            (*this) *= s_big[i];
        delete[] s_big, delete s_num;
    }

    // (*this) = C(n, m) O(klog^2k) k = min(n, m - n)
    void combination(int n, int m)
    {
        if (m < 0 || m > n)
        {
            (*this) = BigNumber((lll)0);
            return;
        }
        else if (!m || m == n)
        {
            (*this) = BigNumber((lll)1);
            return;
        }    
        int number = ((n - m) > m) ? m : (n - m);
        int stack_sz = (int)(ceil(log2(number))) + 1; // level of the calling tree
        if (stack_sz <= 1)
            stack_sz = 2;
        BigNumber* s_big_1 = new BigNumber[stack_sz], * s_big_2 = new BigNumber[stack_sz];
        int* s_num = new int[stack_sz];
        int top = 0;
        for (int i = 1; i <= number; ++i)
        {
            BigNumber b = BigNumber((lll)i), c = BigNumber((lll)(n + 1 - i));
            int level = 1;
            while (top && s_num[top] == level)
                b = b * s_big_1[top], c = c * s_big_2[top--], ++level;
            s_big_1[++top] = b, s_big_2[top] = c, s_num[top] = level;
        }
        (*this) = s_big_2[top];
        BigNumber tmp = s_big_1[top];
        for (int i = top - 1; i; --i)
            (*this) *= s_big_2[i], tmp *= s_big_1[i];
        delete[] s_big_1, delete[] s_big_2, delete s_num;
        (*this) /= tmp;
    }

    // (*this) = catalan(n)
    void catalan(int n)
    {
        assert(n >= 0);
        if (!n)
            (*this) = BigNumber((lll)0);
        else
            combination(n << 1, n), (*this) /= (n + 1);
    }

    void output() const
    {
        char *ch = new char[500500];
        int print_len = 0, is_neg = 0;
        if (!f && (len > 1 || a[0] != 0))
            ch[print_len++] = '-', is_neg = 1;
        int x = a[len - 1];
        while (x)
            ch[print_len++] = (x % base_radix) ^ '0', x /= base_radix;
        for (int i = is_neg; i < print_len + is_neg - 1 - i; i++)
            // std::swap(ch[i], ch[print_len + is_neg - 1 - i]);
            ch[i] ^= ch[print_len + is_neg - 1 - i], ch[print_len + is_neg - 1 - i] ^= ch[i], ch[i] ^= ch[print_len + is_neg - 1 - i];
        for (int i = len - 2; i >= 0; i--)
        {
            x = a[i];
            int k = MODD - 1;
            while (x)
                ch[print_len + k] = (x % base_radix) ^ '0', x /= base_radix, k--;
            while (k >= 0)
                ch[print_len + k] = '0', k--;
            print_len += MODD;
        }
        if (print_len == 0)
            ch[print_len++] = '0';
        ch[print_len] = '\0';
        // printf("%s", ch);
        puts(ch);
        delete ch;
    }    
    ~BigNumber() { delete a; }
};
// above is the Library itself
// below is test area and sample for A / B
struct fastIO
{
    static const int BUFF_SZ = 1 << 18;
    char inbuf[BUFF_SZ], outbuf[BUFF_SZ];
    fastIO()
    {
        setvbuf(stdin, inbuf, _IOFBF, BUFF_SZ);
        setvbuf(stdout, outbuf, _IOFBF, BUFF_SZ);
    }
} IO;

char A[100010], B[100010];
int len_a, len_b;
BigNumber a, b, c, d;
ll _a, _b;
ll x;
int k;
int main() 
{
    // a.input(), b.input(), (a / b).output();
    scanf("%s%s", A, B);
    len_a = strlen(A), len_b = strlen(B);

    if (len_a <= 18 && len_b <= 18)
        sscanf(A, "%lld", &_a),  sscanf(B, "%lld", &_b), printf("%lld", _a / _b);
    else if (len_a < len_b)
        puts("0");
    else if (len_b <= 18)
        a = BigNumber(A), sscanf(B, "%lld", &_b), (a / _b).output();
    else
        a = BigNumber(A), b = BigNumber(B), (a / b).output();

}
