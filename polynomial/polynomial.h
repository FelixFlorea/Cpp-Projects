//
// Created by felix on 10/9/26.
//

#ifndef LIBRARIE_POLINOAME_POLYNOMIAL_REFINED_H
#define LIBRARIE_POLINOAME_POLYNOMIAL_REFINED_H


#include <concepts> //pentru definirea conceptelor algebrice riguroase si pentru a restrictiona template-urile din a face operatii absurde
#include <ostream>  //optimizare mica dar buna totusi
#include <istream>  //pentru a inlocui iostream caci nu folosesc cout, cin, ci doar <<, >>
#include <vector>   //pentru polinoame care au in mod constant dimensiuni variabile
#include <initializer_list> //pentru usurinta de a lucra cu vectori, utila la polinoame mai tarziu




template <typename T>
inline constexpr bool is_commutative = false;  //implicit, niciun tip nu e considerat comutativ, daca nu specific eu altfel

//lista tipurilor de baza pe care le declar comutative (eticheta nu mai sta in interiorul clasei, ca sa pot include si double)
template <> inline constexpr bool is_commutative<int> = true;
template <> inline constexpr bool is_commutative<long long> = true;
template <> inline constexpr bool is_commutative<float> = true;
template <> inline constexpr bool is_commutative<double> = true;

//inversul pentru tipurile de baza; trebuie declarat inaintea conceptelor, ca field sa le gaseasca
inline double inverse(double a) { return 1.0 / a; }
inline float  inverse(float a)  { return 1.0f / a; }

template <typename T>
concept ring=           //definesc conceptul de inel
    std::regular<T> and //elementul poate fi creat gol, copiat, distrus, si comparat
    requires(T a, T b)
{
    {-a} -> std::same_as<T>;
    {a + b} -> std::same_as<T>;
    {a - b} -> std::same_as<T>;
    {a * b} -> std::same_as<T>;
    {a += b} -> std::same_as<T&>;
    {a -= b} -> std::same_as<T&>;
    {a *= b} -> std::same_as<T&>;
    {T(0)} -> std::same_as<T>;
    {T(1)} -> std::same_as<T>;
};

template <typename T>   //definesc conceptul de inel comutativ: inel + tipul se afla pe lista is_commutative
concept commutative_ring = ring<T> and is_commutative<T>;

template <typename T>   //vom presupune ca toate corpurile sunt comutative(nu lucram cu quaternioni si octonioni etc)
concept field = commutative_ring<T> and
    requires(T a, T b)
{
    {inverse(a)} -> std::same_as<T>;    //functie separata
    {a / b} -> std::same_as<T>;
    {a /= b} -> std::same_as<T&>;
};

template <unsigned int mod>     //nu avem nevoie de concepts aici. template ul T sau in cazul actual unsigned int nu trebuie sa verifice conditie de corp/inel
class modular_element {         //consideram mereu mod nr prim
private:
    int v;                         //!!!voi modifica in curand din int v, intr un T v, iar tipul de date T va fi "integral", si de asemenea templateul clasei
    struct direct_tag {};
    modular_element(int val, direct_tag) {  // retin valoarea fara sa mai calculez %, foarte util!!!!
        v = val;
    }
public:
    modular_element() {     //constructor implicit, intializeaza cu 0 elementul modular
        v=0;
    }

    explicit modular_element(long long x) {      //constructor de conversie, lent dar sigur, obligatoriu explicit, ca sa nu faca conversii pe ascuns
        long long remainder = x % mod;
        if (remainder < 0)
            remainder += mod;
        v =static_cast<int>(remainder);     //elimina avertismentele compilatorului privind pierderea de precizie long long -> int
    }

    int val() const {   //metoda publica read-only care nu modifica v
        return v;
    }

    modular_element operator-() const {
        if (v == 0)
            return modular_element(0, direct_tag{});    //nu mai fac operatia % ci direct iau elementul modular care stiu ca se afla in Zp
        return modular_element(mod - v, direct_tag{});
    }

    friend modular_element operator+(const modular_element& x, const modular_element& y) {
        int res=x.v+y.v;
        if (res >= mod)
            res = res - mod;
        return modular_element(res, direct_tag{});
    }

    friend modular_element operator-(const modular_element& x, const modular_element& y) {
        int res=x.v-y.v;
        if (res < 0)
            res = res + mod;
        return modular_element(res, direct_tag{});
    }

    friend modular_element operator*(const modular_element& x, const modular_element& y) {
        int res = static_cast<int>((1LL * x.v * y.v) % mod);
        return modular_element(res, direct_tag{});
    }

    static modular_element pow_log(modular_element base, long long exp) {
        if (exp < 0)
            return pow_log(base.inverse(), -exp);       //daca am exponent negativ, pentru inversul inversului
        modular_element res(1);
        while (exp > 0) {
            if (exp & 1)
                res = res * base;
            base = base * base;
            exp >>= 1;
        }
        return res;
    }

    modular_element inverse() const {
        if (v == 0)
            throw "0 has no inverse";
        return pow_log(modular_element(v),mod-2);           // trebuie mentionat ca v trebuie adus la forma de element modular
    }

    friend modular_element operator/(const modular_element& x, const modular_element& y) {
        return x*y.inverse();
    }

    friend modular_element& operator+=(modular_element& x, const modular_element& y) {
        x.v += y.v;
        if (x.v >= mod)
            x.v -= mod;
        return x;
    }

    friend modular_element& operator-=(modular_element& x, const modular_element& y) {
        x.v -= y.v;
        if (x.v < 0)
            x.v += mod;
        return x;
    }

    friend modular_element& operator*=(modular_element& x, const modular_element& y) {
        x.v = static_cast<int>((1LL * x.v * y.v) % mod);
        return x;
    }

    friend modular_element& operator/=(modular_element& x, const modular_element& y) {
        x = x / y;
        return x;
    }

    friend bool operator==(const modular_element& x, const modular_element& y) {
        return x.v == y.v;
    }

    friend std::ostream& operator<<(std::ostream& out, const modular_element& p) {
        out << p.v;
        return out;
    }

    friend std::istream& operator>>(std::istream& in, modular_element& p) {
        long long x;
        if (!(in >> x)) return in;
        p = modular_element(x);
        return in;
    }
};

template <unsigned int mod>
inline constexpr bool is_commutative<modular_element<mod>> = true;  //trec elementele modulare pe lista tipurilor comutative

template <unsigned int mod>
modular_element<mod> inverse(const modular_element<mod>& a) {       //inversul ca functie separata, cerut de conceptul field
    return a.inverse();
}

template <typename T, unsigned int n>        //doar matrice patratice momentan, n>=1
requires ring<T>        //elementele din matrice trebuie sa fie dintr un inel pentru a putea efectua operatii
struct matrix {
    T m[n][n];

    matrix() {          //constructor implicit
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                m[i][j] = T(0);
    }
                                    //voi mai adauga un constructor ca sa pot declara matricile matrix<int> m={{0,1}, {0,1}} de exemplu si ca sa pot trece de conceptul ring, caci nu am constructor matrix(1)
    static matrix identity() {
        matrix I;
        for (int i = 0; i < n; i++)
            I.m[i][i] = T(1);
        return I;
    }

    matrix operator-() const {
        matrix result;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                result.m[i][j] = -m[i][j];
        return result;
    }

    friend matrix operator+(const matrix& m1, const matrix& m2) {
        matrix result;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                result.m[i][j] = m1.m[i][j] + m2.m[i][j];
        return result;
    }

    friend matrix operator-(const matrix& m1, const matrix& m2) {
        matrix result;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                result.m[i][j] = m1.m[i][j] - m2.m[i][j];
        return result;
    }

    friend matrix operator*(const matrix& m1, const matrix& m2) {
        matrix result;
        for (int i = 0; i < n; i++)
            for (int k = 0; k < n; k++)
                for (int j = 0; j < n; j++)
                    result.m[i][j] = result.m[i][j] + m1.m[i][k] * m2.m[k][j];
        return result;
    }

    friend matrix& operator+=(matrix& m1, const matrix& m2) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                m1.m[i][j] = m1.m[i][j] + m2.m[i][j];
        return m1;
    }

    friend matrix& operator-=(matrix& m1, const matrix& m2) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                m1.m[i][j] = m1.m[i][j] - m2.m[i][j];
        return m1;
    }

    friend matrix& operator*=(matrix& m1, const matrix& m2) {
        m1 = m1 * m2;      
        return m1;
    }

    friend bool operator==(const matrix& m1, const matrix& m2) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (m1.m[i][j] != m2.m[i][j])
                    return false;
        return true;
    }

    friend matrix transpose(const matrix& mat) {
        matrix result;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                result.m[i][j] = mat.m[j][i];
        return result;
    }

    friend T trace(const matrix& mat) {
        T result = T(0);
        for (int i = 0; i < n; i++)
            result = result + mat.m[i][i];
        return result;
    }

    static matrix pow_log(matrix base, long long exp) {
        matrix result = matrix::identity();
        while (exp > 0) {
            if (exp & 1)
                result = result * base;
            base = base * base;
            exp >>= 1;
        }
        return result;
    }

    friend std::ostream& operator<<(std::ostream& out, const matrix& mat) {
        for (int i = 0; i < n; i++) {
            out << "[ ";
            for (int j = 0; j < n; j++) {
                out << mat.m[i][j];
                if (j < n - 1) out << " ";
            }
            out << " ]\n";
        }
        return out;
    }

    friend std::istream& operator>>(std::istream& in, matrix& mat) {    //citeste cele n*n elemente linie cu linie, fara paranteze
        matrix res;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (!(in >> res.m[i][j]))
                    return in;
        mat = res;
        return in;
    }
};

template <typename T>
requires ring<T>
class polynomial {
private:
    std::vector<T> a;

    void normalize() {      //e privata, deoarece polinomul trebuie sa fie mereu curat, fiind mereu chemata functia de constructori
        while (a.size() > 1 && a.back() == T(0)) {      //ultimul coeficient nu e niciodată 0, cu excepția polinomului nul
            a.pop_back();
        }
    }       //daca nu ar fi privat, unele functii,(exemplul degree) ar da rezultate false

public:
    polynomial() {
        a.push_back(T(0));  //constructor implicit, mereu am polinomul nul
    }

    polynomial(std::initializer_list<T> list) {     //imi permite sa scriu mai usor polinomul ex: polynomial<int> p={1,2,3};
        for (auto x : list) {                     //!!! atentie, daca vreau elemente modulare voi scrie polynomial<Z7> p={Z7(1), Z7(2)}; datorita constructorului explicit la elemente modulare
            a.push_back(x);
        }
        if (a.empty()) {
            a.push_back(T(0));
        }
        normalize();
    }

    explicit polynomial(const std::vector<T>& vec) {    //constructor de conversie, converteste vectori in polinoame
        a = vec;                                        //explicit, nimic pe ascuns sau nedorit
        if (a.empty()) {
            a.push_back(T(0));
        }
        normalize();
    }

    int degree() const {
        return (int)a.size() - 1;           //!!!a.size() returneaza numar fara semn, ma ajuta mai incolo la T() T point si la [] si la multe altele
    }

    T operator[](int i) const {
        if (i < 0)
            throw "negative index";
        if (i >= (int)a.size())
            return T(0);
        return a[i];
    }

    void set_coef(int i, T val) {
        if (i < 0)
            throw "negative index";
        if (i >= (int)a.size()) {
            if (val == T(0))
                return;
            a.resize(i + 1, T(0));
        }
        a[i] = val;
        normalize();
    }

    T operator()(T point) const {
        int deg = degree();
        T val = a[deg];
        for (int i = deg - 1; i >= 0; i--) {
            val = val * point + a[i];
        }
        return val;
    }

    const std::vector<T>& get_coefficients() const {
        return a;
    }
};

template <typename T>
requires ring<T>
int degree(const polynomial<T>& P) {
    return P.degree();
}

template <typename T>
requires ring<T>
polynomial<T> shift_left(int k, const polynomial<T>& P) {
    if (degree(P) == 0 && P[0] == T(0))
        return P;
    std::vector<T> vec;
    for (int i = 0; i < k; i++) {
        vec.push_back(T(0));
    }
    for (int i = 0; i <= degree(P); i++) {
        vec.push_back(P[i]);
    }
    return polynomial<T>(vec);
}

template <typename T>
requires ring<T>
polynomial<T> shift_right(int k, const polynomial<T>& P) {
    if (k > degree(P))
        return polynomial<T>();
    std::vector<T> vec;
    for (int i = k; i <= degree(P); i++) {
        vec.push_back(P[i]);
    }
    return polynomial<T>(vec);
}

template <typename T>
requires ring<T>
polynomial<T> subpolynomial(int start, int stop, const polynomial<T>& P) {
    std::vector<T> vec;
    for (int i = start; i <= stop; i++) {
        vec.push_back(P[i]);
    }
    return polynomial<T>(vec);
}

template <typename T>
requires ring<T>
polynomial<T> add(const polynomial<T>& P1, const polynomial<T>& P2) {
    int deg1 = degree(P1), deg2 = degree(P2);
    int degp = (deg1 > deg2 ? deg1 : deg2);
    std::vector<T> vec(degp + 1, T(0));
    for (int i = 0; i <= degp; i++) {
        vec[i] = P1[i] + P2[i];
    }
    return polynomial<T>(vec);
}

template <typename T>
requires ring<T>
polynomial<T> subtract(const polynomial<T>& P1, const polynomial<T>& P2) {
    int deg1 = degree(P1), deg2 = degree(P2);
    int degp = (deg1 > deg2 ? deg1 : deg2);
    std::vector<T> vec(degp + 1, T(0));
    for (int i = 0; i <= degp; i++) {
        vec[i] = P1[i] - P2[i];
    }
    return polynomial<T>(vec);
}

template <typename T>
requires ring<T>
polynomial<T> scalar_multiply(T c, const polynomial<T>& P) {
    if (c == T(0))
        return polynomial<T>();
    int degp = degree(P);
    std::vector<T> vec(degp + 1, T(0));
    for (int i = 0; i <= degp; i++) {
        vec[i] = c * P[i];
    }
    return polynomial<T>(vec);
}

template <typename T>
requires ring<T>
polynomial<T> multiply_brute(const polynomial<T>& P1, const polynomial<T>& P2) {
    int deg1 = degree(P1), deg2 = degree(P2), degp = deg1 + deg2;
    std::vector<T> vec(degp + 1, T(0));
    for (int i = 0; i <= deg1; i++) {
        for (int j = 0; j <= deg2; j++) {
            vec[i + j] += P1[i] * P2[j];
        }
    }
    return polynomial<T>(vec);
}

template <typename T>
requires commutative_ring<T>
polynomial<T> multiply_karatsuba(const polynomial<T>& P1, const polynomial<T>& P2) {
    int deg1 = degree(P1);
    int deg2 = degree(P2);
    if (deg1 == 0)
        return scalar_multiply(P1[0], P2);
    if (deg2 == 0)
        return scalar_multiply(P2[0], P1);
    if (deg1 < 32 && deg2 < 32)
        return multiply_brute(P1, P2);

    int max_deg = (deg1 > deg2 ? deg1 : deg2);
    int m = (max_deg + 1) / 2;

    polynomial<T> P1_low = subpolynomial(0, m - 1, P1);
    polynomial<T> P1_high = shift_right(m, P1);
    polynomial<T> P2_low = subpolynomial(0, m - 1, P2);
    polynomial<T> P2_high = shift_right(m, P2);

    polynomial<T> A = multiply_karatsuba(P1_high, P2_high);
    polynomial<T> B = multiply_karatsuba(P1_low, P2_low);
    polynomial<T> C = multiply_karatsuba(add(P1_high, P1_low), add(P2_high, P2_low));

    polynomial<T> middle = subtract(subtract(C, A), B);
    polynomial<T> res1 = shift_left(2 * m, A);
    polynomial<T> res2 = shift_left(m, middle);

    return add(add(res1, res2), B);
}

template <typename T>
requires field<T>
polynomial<T> divide_quotient(const polynomial<T>& P1, const polynomial<T>& P2) {
    if (degree(P2) == 0 && P2[0] == T(0)) {
        throw "Division by zero polynomial";
    }
    std::vector<T> r_vec = P1.get_coefficients();
    int deg2 = degree(P2);
    if ((int)r_vec.size() - 1 < deg2) {
        return polynomial<T>();
    }
    int deg_q = (int)r_vec.size() - 1 - deg2;
    std::vector<T> q_vec(deg_q + 1, T(0));

    while ((int)r_vec.size() - 1 >= deg2 && !(r_vec.size() == 1 && r_vec[0] == T(0))) {
        int curr_deg_r = (int)r_vec.size() - 1;
        T factor = r_vec[curr_deg_r] / P2[deg2];
        int pow = curr_deg_r - deg2;
        q_vec[pow] = factor;

        for (int i = 0; i <= deg2; i++) {
            r_vec[pow + i] -= factor * P2[i];
        }
        r_vec.pop_back();   //coeficientul de sus e 0 matematic; il scot fara sa verific, caci la double poate ramane ceva de genul 1e-16 din rotunjiri
        while (r_vec.size() > 1 && r_vec.back() == T(0)) {
            r_vec.pop_back();
        }
    }
    return polynomial<T>(q_vec);
}

template <typename T>
requires field<T>
polynomial<T> divide_remainder(const polynomial<T>& P1, const polynomial<T>& P2) {
    if (degree(P2) == 0 && P2[0] == T(0)) {
        throw "Division by zero polynomial";
    }
    std::vector<T> r_vec = P1.get_coefficients();
    int deg2 = degree(P2);

    while ((int)r_vec.size() - 1 >= deg2 && !(r_vec.size() == 1 && r_vec[0] == T(0))) {
        int curr_deg_r = (int)r_vec.size() - 1;
        T factor = r_vec[curr_deg_r] / P2[deg2];
        int pow = curr_deg_r - deg2;

        for (int i = 0; i <= deg2; i++) {
            r_vec[pow + i] -= factor * P2[i];
        }
        r_vec.pop_back();   //coeficientul de sus e 0 matematic; il scot fara sa verific, caci la double poate ramane ceva de genul 1e-16 din rotunjiri
        while (r_vec.size() > 1 && r_vec.back() == T(0)) {
            r_vec.pop_back();
        }
    }
    return polynomial<T>(r_vec);
}

template <typename T>
requires ring<T>
polynomial<T> operator+(const polynomial<T>& P1, const polynomial<T>& P2) {
    return add(P1, P2);
}

template <typename T>
requires ring<T>
polynomial<T> operator-(const polynomial<T>& P1, const polynomial<T>& P2) {
    return subtract(P1, P2);
}

template <typename T>
requires ring<T>
polynomial<T> operator*(T c, const polynomial<T>& P) {
    return scalar_multiply(c, P);
}

template <typename T>
requires ring<T>
polynomial<T> operator*(const polynomial<T>& P, T c) {
    return multiply_brute(P, polynomial<T>{c});
}

template <typename T>
requires ring<T> && (!commutative_ring<T>)          //momentan nu ajuta la nimic
polynomial<T> operator*(const polynomial<T>& P1, const polynomial<T>& P2) {
    return multiply_brute(P1, P2);
}

template <typename T>
requires commutative_ring<T>
polynomial<T> operator*(const polynomial<T>& P1, const polynomial<T>& P2) {
    int max_deg = (degree(P1) > degree(P2) ? degree(P1) : degree(P2));
    if (max_deg < 64)
        return multiply_brute(P1, P2);
    else
        return multiply_karatsuba(P1, P2);
}

template <typename T>
requires field<T>
polynomial<T> operator/(const polynomial<T>& P1, const polynomial<T>& P2) {
    return divide_quotient(P1, P2);
}

template <typename T>
requires field<T>
polynomial<T> operator%(const polynomial<T>& P1, const polynomial<T>& P2) {
    return divide_remainder(P1, P2);
}

template <typename T>
requires ring<T>
polynomial<T>& operator+=(polynomial<T>& P1, const polynomial<T>& P2) {
    P1 = P1 + P2;
    return P1;
}

template <typename T>
requires ring<T>
polynomial<T>& operator-=(polynomial<T>& P1, const polynomial<T>& P2) {
    P1 = P1 - P2;
    return P1;
}

template <typename T>
requires ring<T>
polynomial<T>& operator*=(polynomial<T>& P1, const polynomial<T>& P2) {
    P1 = P1 * P2;
    return P1;
}

template <typename T>
requires ring<T>
bool operator==(const polynomial<T>& P1, const polynomial<T>& P2) {
    return P1.get_coefficients() == P2.get_coefficients();
}

template <typename T>
requires ring<T>
polynomial<T> derivative(const polynomial<T>& P) {
    int degp = degree(P);
    std::vector<T> vec(degp, T(0));
    for (int i = 1; i <= degp; i++) {
        vec[i - 1] = P[i] * T(i);
    }
    return polynomial<T>(vec);
}

template <typename T>
requires field<T>
polynomial<T> antiderivative(const polynomial<T>& P) {     //integral face parte din namespace ul std, asa ca am modificat din integral in antiderivative
    int degp = degree(P);
    std::vector<T> vec(degp + 2, T(0));
    for (int i = 0; i <= degp; i++) {
        if (P[i] == T(0))
            continue;
        if (T(i + 1) == T(0))
            throw "Integral undefined in this characteristic";      //eroare pentru corpuri cu caracteristica p, nu ai cum sa integrezi
        vec[i + 1] = P[i] / T(i + 1);
    }
    return polynomial<T>(vec);
}

template <typename T>
requires ring<T>
std::ostream& operator<<(std::ostream& out, const polynomial<T>& P) {
    int deg_p = degree(P);
    bool first = true;
    for (int i = deg_p; i >= 0; i--) {
        T coef = P[i];
        if (coef == T(0))
            continue;

        if (!first) out << "+";
        first = false;

        if (i == 0) {
            out << coef;
        } else if (i == 1) {
            if (coef == T(1)) out << "X";
            else out << coef << "X";
        } else {
            if (coef == T(1)) out << "X^" << i;
            else out << coef << "X^" << i;
        }
    }
    if (first) {
        out << T(0);
    }
    return out;
}

template <typename T>   //citeste gradul, apoi cei grad+1 coeficienti in ordine crescatoare (de la X^0 la X^grad)
requires ring<T>
std::istream& operator>>(std::istream& in, polynomial<T>& P) {
    int deg;
    if (!(in >> deg) || deg < 0) {
        in.setstate(std::ios::failbit);
        return in;
    }
    std::vector<T> vec;
    for (int i = 0; i <= deg; i++) {
        T coef;
        if (!(in >> coef))
            return in;
        vec.push_back(coef);
    }
    P = polynomial<T>(vec);
    return in;
}


#endif //LIBRARIE_POLINOAME_POLYNOMIAL_REFINED_H
