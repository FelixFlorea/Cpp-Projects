#ifndef POLYNOM_H_INCLUDED
#define POLYNOM_H_INCLUDED

#include <iostream>
#include <vector>
#include <concepts>
#include <initializer_list>

template <typename T>
concept ring=requires(T a, T b){
    {a+b}->std::same_as<T>;
    {a-b}->std::same_as<T>;
    {a*b}->std::same_as<T>;
    {T(0)}->std::same_as<T>;
};

template <typename T>
concept field=ring<T> && requires(T a, T b){
    {a/b}->std::same_as<T>;
    {T(1)}->std::same_as<T>;
};

template <int mod>
struct modular_element{
    int v;
    modular_element():v(0){}
    modular_element(long long x){
        long long remainder=x%mod;
        if (remainder<0)
            remainder+=mod;
        v=static_cast<int>(remainder);
    }
    friend modular_element operator+(const modular_element& p1, const modular_element& p2){
        return modular_element(p1.v+p2.v);
    }
    friend modular_element operator-(const modular_element& p1, const modular_element& p2){
        return modular_element(p1.v-p2.v+mod);
    }
    friend modular_element operator*(const modular_element& p1, const modular_element& p2){
        return modular_element(1LL*p1.v*p2.v);
    }
    static long long pow_log(long long base, long long exp) {
        long long result=1;
        base%=mod;
        while (exp>0) {
            if (exp%2==1) 
                result=(result*base)%mod;
            base=(base*base)%mod;
            exp/=2;
        }
        return result;
    }
    modular_element modular_inverse() const {
        return modular_element(pow_log(v, mod-2));
    }
    friend modular_element operator/(const modular_element& p1, const modular_element& p2){
        return p1*p2.modular_inverse();
    }
    modular_element& operator+=(const modular_element& obj){
        v=(v+obj.v)%mod;
        return *this;
    }
    modular_element& operator-=(const modular_element& obj){
        v=(v-obj.v+mod)%mod;
        return *this;
    }
    modular_element& operator*=(const modular_element& obj){
        v=(1LL*v*obj.v)%mod;
        return *this;
    }
    modular_element& operator/=(const modular_element& obj){
        *this=*this*obj.modular_inverse();
        return *this;
    }
    friend bool operator==(const modular_element& p1, const modular_element& p2) {
        return p1.v==p2.v;
    }
    friend bool operator!=(const modular_element& p1, const modular_element& p2) {
        return p1.v!=p2.v;
    }
    friend bool operator==(const modular_element& p1, long long val){
        return p1.v==((val%mod+mod))%mod;
    }
    friend bool operator!=(const modular_element& p1, long long val){
        return !(p1==val);
    }
    friend std::ostream& operator<<(std::ostream& out, const modular_element& p) {
        out<<p.v;
        return out;
    }
    friend std::istream& operator>>(std::istream& in, modular_element& p) {
        long long x;
        in>>x;
        long long remainder=x%mod;
        if (remainder<0) 
            remainder+=mod;
        p.v=static_cast<int>(remainder);
        return in;
    }
};

template <typename T>
requires ring<T>
struct polynomial{
    std::vector<T> a;
    polynomial(){
        a.push_back(T(0));
    }
    polynomial(int degree){
        a.resize(degree+1, T(0));
    }
    polynomial(std::initializer_list<T> list){
        for(auto x : list){
            a.push_back(x);
        }
        if(a.empty()){
            a.push_back(T(0));
        }
    }
    void set_dimension(int n){
        a.resize(n, T(0));
        if(a.empty()){
            a.push_back(T(0));
        }
    }
    T operator[](int i) const {
        if (i<0)
            throw std::out_of_range("negative index");
        if (i>=a.size())
            return T(0);
        return a[i];
    }
    T& operator[](int i) {
        if (i<0)
            throw std::out_of_range("negative index");
        if(i>=a.size())
            a.resize(i+1, T(0));
        return a[i];
    }
    T operator()(T point) const {
        int deg=a.size()-1;
        if (deg<0)
            return T(0);
        T val=a[deg];
        for (int i=deg-1;i>=0;i--) {
            val=val*point+a[i];
        }
        return val;
    }
};

template <typename T>
requires ring<T>
int degree(const polynomial<T>& P){
    if(P.a.empty()) return 0;
    return P.a.size()-1;
}

template <typename T>
requires ring<T>
T coefficient(int i, const polynomial<T>& P){
    if (i<0 or i>degree(P))
        return T(0);
    else
        return P.a[i];
}

template <typename T>
requires ring<T>
polynomial<T> trim(const polynomial<T>& P){
    polynomial<T> copy=P;
    while(copy.a.size()>1 and copy.a.back()==T(0)){
        copy.a.pop_back();
    }
    return copy;
}

template <typename T>
requires ring<T>
polynomial<T> shift_left(int k, const polynomial<T>& P){
    if (P.a.size()==1 && P.a[0]==T(0)) 
        return P;
    polynomial<T> copy=P;
    copy.a.insert(copy.a.begin(), k, T(0));
    return copy;
}

template <typename T>
requires ring<T>
polynomial<T> shift_right(int k, const polynomial<T>& P) {
    if (k>=P.a.size())
        return polynomial<T>{{T(0)}};
    polynomial<T> copy=P;
    copy.a.erase(copy.a.begin(), copy.a.begin()+k);
    return copy;
}

template <typename T>
requires ring<T>
polynomial<T> subpolynomial(int start, int stop, const polynomial<T>& P){
    polynomial<T> result;
    result.a.clear();
    for (int i=start; i<=stop; i++){
        result.a.push_back(coefficient(i, P));
    }
    return result;
}

template <typename T>
requires ring<T>
polynomial<T> add(const polynomial<T>& P1, const polynomial<T>& P2){
    int deg1=degree(P1),deg2=degree(P2),degp=std::max(deg1,deg2);
    polynomial<T> P;
    P.a.resize(degp+1,T(0));
    for (int i=0;i<=degp;i++){
        P.a[i]=coefficient(i,P1)+coefficient(i,P2);
    }
    return P;
}

template <typename T>
requires ring<T>
polynomial<T> subtract(const polynomial<T>& P1, const polynomial<T>& P2){
    int deg1=degree(P1),deg2=degree(P2),degp=std::max(deg1,deg2);
    polynomial<T> P;
    P.a.resize(degp+1,T(0));
    for (int i=0;i<=degp;i++){
        P.a[i]=coefficient(i,P1)-coefficient(i,P2);
    }
    return P;
}

template <typename T>
requires ring<T>
polynomial<T> scalar_multiply(T c, const polynomial<T>& P){
    polynomial<T> copy=P;
    int degp=degree(copy);
    for (int i=0;i<=degp;i++)
        copy.a[i]=copy.a[i]*c;
    return copy;
}

template <typename T>
requires ring<T>
polynomial<T> multiply_brute(const polynomial<T>& P1, const polynomial<T>& P2){
    polynomial<T> P;
    int deg1=degree(P1),deg2=degree(P2),degp=deg1+deg2;
    P.a.resize(degp+1,T(0));
    for (int i=0;i<=deg1;i++){
        for (int j=0;j<=deg2;j++){
            P.a[i+j]+=P1.a[i]*P2.a[j];
        }
    }
    return P;
}

template <typename T>
requires ring<T>
polynomial<T> multiply_karatsuba(const polynomial<T>& P1, const polynomial<T>& P2){
    polynomial<T> clean1=trim(P1);
    polynomial<T> clean2=trim(P2);
    int deg1=degree(clean1);
    int deg2=degree(clean2);
    if (deg1<0 or deg2<0) {
        return polynomial<T>();
    }
    if (deg1==0) {
        return scalar_multiply(clean1.a[0], clean2);
    }
    if (deg2==0) {
        return scalar_multiply(clean2.a[0], clean1);
    }
    if (deg1<32 and deg2<32) {
        return multiply_brute(clean1, clean2);
    }
    int m=(std::max(deg1, deg2)+1)/2;
    polynomial<T> P1_low=subpolynomial(0, m-1, clean1);
    polynomial<T> P1_high=shift_right(m, clean1);
    polynomial<T> P2_low=subpolynomial(0, m-1, clean2);
    polynomial<T> P2_high=shift_right(m, clean2);
    polynomial<T> A=multiply_karatsuba(P1_high, P2_high);
    polynomial<T> B=multiply_karatsuba(P1_low, P2_low);
    polynomial<T> C=multiply_karatsuba(add(P1_high, P1_low), add(P2_high, P2_low));
    polynomial<T> middle=subtract(subtract(C, A), B);
    polynomial<T> res1=shift_left(2*m, A);
    polynomial<T> res2=shift_left(m, middle);
    return trim(add(add(res1, res2), B));
}

template <typename T>
requires field<T>
polynomial<T> divide_quotient(const polynomial<T>& P1, const polynomial<T>& P2){
    polynomial<T> clean_p2=trim(P2);
    if (clean_p2.a.size()==1 && clean_p2.a[0]==T(0)){
        return polynomial<T>();
    }
    polynomial<T> R=trim(P1);
    polynomial<T> Q;
    int deg2=degree(clean_p2);
    if (degree(R)<deg2){
        return polynomial<T>();
    }
    int deg_q=degree(R)-deg2;
    Q.a.resize(deg_q+1,T(0));
    while (degree(R)>=deg2&&!(R.a.size()==1&&R.a[0]==T(0))){
        int curr_deg_r=degree(R);
        T factor=R.a[curr_deg_r]/clean_p2.a[deg2];
        int pow=curr_deg_r-deg2;
        Q.a[pow]=factor;
        polynomial<T> t;
        t.a.clear();
        t.a.resize(pow+1,T(0));
        t.a[pow]=factor;
        R=trim(R-(t*clean_p2));
    }
    return trim(Q);
}

template <typename T>
requires field<T>
polynomial<T> divide_remainder(const polynomial<T>& P1, const polynomial<T>& P2){
    polynomial<T> clean_p2=trim(P2);
    if (clean_p2.a.size()==1 && clean_p2.a[0]==T(0)){
        return trim(P1);
    }
    polynomial<T> R=trim(P1);
    int deg2=degree(clean_p2);
    while (degree(R)>=deg2 && !(R.a.size()==1 && R.a[0]==T(0))){
        int curr_deg_r=degree(R);
        T factor=R.a[curr_deg_r]/clean_p2.a[deg2];
        int pow=curr_deg_r-deg2;
        polynomial<T> t;
        t.a.clear();
        t.a.resize(pow+1, T(0));
        t.a[pow]=factor;
        R=trim(R-(t*clean_p2));
    }
    return R;
}

template <typename T>
requires ring<T>
polynomial<T> operator+(const polynomial<T>& P1, const polynomial<T>& P2){
    return add(P1,P2);
}

template <typename T>
requires ring<T>
polynomial<T> operator-(const polynomial<T>& P1, const polynomial<T>& P2){
    return subtract(P1,P2);
}

template <typename T>
requires ring<T>
polynomial<T> operator*(T c, const polynomial<T>& P){
    return scalar_multiply(c,P);
}

template <typename T>
requires ring<T>
polynomial<T> operator*(const polynomial<T>& P, T c){
    return scalar_multiply(c,P);
}

template <typename T>
requires ring<T>
polynomial<T> operator*(const polynomial<T>& P1, const polynomial<T>& P2){
    int max_deg=std::max(degree(P1), degree(P2));
    if (max_deg<64) 
        return multiply_brute(P1, P2);
    else 
        return multiply_karatsuba(P1, P2);
}

template <typename T>
requires ring<T>
polynomial<T>& operator+=(polynomial<T>& P1, const polynomial<T>& P2){
    P1=add(P1, P2);
    return P1;
}

template <typename T>
requires ring<T>
polynomial<T>& operator-=(polynomial<T>& P1, const polynomial<T>& P2){
    P1=subtract(P1, P2);
    return P1;
}

template <typename T>
requires ring<T>
polynomial<T>& operator*=(polynomial<T>& P1, const polynomial<T>& P2){
    P1=P1*P2;
    return P1;
}

template <typename T>
requires ring<T>
bool operator==(const polynomial<T>& P1, const polynomial<T>& P2){
    polynomial<T> clean1=trim(P1);
    polynomial<T> clean2=trim(P2);
    return clean1.a==clean2.a;
}

template <typename T>
requires ring<T>
bool operator!=(const polynomial<T>& P1, const polynomial<T>& P2){
    polynomial<T> clean1=trim(P1);
    polynomial<T> clean2=trim(P2);
    return clean1.a!=clean2.a;
}

template <typename T>
requires field<T>
polynomial<T> operator/(const polynomial<T>& P1, const polynomial<T>& P2){
    return divide_quotient(P1,P2);
}

template <typename T>
requires field<T>
polynomial<T> operator%(const polynomial<T>& P1, const polynomial<T>& P2){
    return divide_remainder(P1,P2);
}

template <typename T>
requires ring<T>
polynomial<T> operator--(const polynomial<T>& P){
    if (P.a.size() <= 1)
        return polynomial<T>{{T(0)}};
    polynomial<T> result;
    int degp=degree(P);
    result.a.resize(degp);
    for (int i=1;i<=degp;i++){
        result.a[i-1]=P.a[i]*i;
    }
    return result;
}

template <typename T>
requires field<T>
polynomial<T> operator++(const polynomial<T>& P){
    polynomial<T> result;
    int degp=degree(P);
    result.a.resize(degp+2, T(0));
    for (int i=0; i<=degp; i++){
        result.a[i+1]=P.a[i]/(i+1);
    }
    return result;
}

template <typename T>
requires ring<T>
std::ostream& operator<<(std::ostream& out, const polynomial<T>& P){
    polynomial<T> clean=trim(P);
    int deg_p=degree(clean);
    bool first=true;
    for (int i=deg_p;i>=0;i--){
        T coef=clean.a[i];
        if(coef==T(0)){
            continue;
        }
        if(!first){
            out<<"+";
        }
        first=false;
        if(i==0){
            out<<coef;
        }else if(i==1){
            if(coef==T(1)){
                out<<"X";
            }else{
                out<<coef<<"X";
            }
        }else{
            if(coef==T(1)){
                out<<"X^"<<i;
            }else{
                out<<coef<<"X^"<<i;
            }
        }
    }
    if(first){
        out<<T(0);
    }
    return out;
}

#endif // POLYNOM_H_INCLUDED
