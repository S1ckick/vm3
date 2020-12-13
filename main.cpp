#include <iostream>
#include <cmath>

#define SPLITTER 134217729.0               // = 2^27 + 1
#define SPLIT_THRESH 6.69692879491417e+299 // = 2^996


using namespace std;

void d_sum_d(double a, double b, double &s, double &err);
void d_sub_d(double a, double b, double &s, double &err);
void nint(double d, double &res);
void dd_sum_dd(double x_h, double x_l, double y_h, double y_l, double &t_h, double &t_l);
void split(double a, double &hi, double &lo);
void d_mul_d(double a, double b, double &s, double &err);
void dd_mul_dd(double x_h, double x_l, double y_h, double y_l, double &t_h, double &t_l);
void nint_dd(double &x_h, double &x_l);
void dd_mul_d(double x_h, double x_l, double y, double &res_h, double &res_l);
void d_qsum_d(double a, double b, double &s, double &err);
void dd_sub_dd(double x_h, double x_l, double y_h, double y_l, double &res_h, double &res_l);
void dd_add_d(double x_h, double x_l, double b, double &res_h, double &res_l);
void dd_div_dd(double x_h, double x_l, double y_h, double y_l, double &res_h, double &res_l);
void sin_taylor(double a_h, double a_l, double &res_h, double &res_l);
void cos_taylor(double a_h, double a_l, double &res_h, double &res_l);
bool dd_less_or_equal_dd(double x_h, double x_l, double y_h, double y_l);
void sin(double a_h, double a_l, double &res_h, double &res_l);

void d_sum_d(double a, double b, double &s, double &err) {
    s = a + b;
    double bb = s - a;
    err = (a - (s - bb)) + (b - bb);
}


void d_qsum_d(double a, double b, double &s, double &err){
    s = a + b;
    err = b - (s - a);
}


void d_sub_d(double a, double b, double &s, double &err){
    s = a - b;
    double bb = s - a;
    err = (a - (s - bb)) - (b + bb);
}

void nint(double d, double &res){
    if(d == floor(d)){
        res = d;
        return;
    }
    res = floor(d+0.5);
}


void split(double a, double &hi, double &lo) {
    double temp;
    if (a > SPLIT_THRESH || a < -SPLIT_THRESH) {
        a *= 3.7252902984619140625e-09;  // 2^-28
        temp = SPLITTER * a;
        hi = temp - (temp - a);
        lo = a - hi;
        hi *= 268435456.0;          // 2^28
        lo *= 268435456.0;          // 2^28
    } else {
        temp = SPLITTER * a;
        hi = temp - (temp - a);
        lo = a - hi;
    }
}


void d_mul_d(double a, double b, double &s, double &err) {
    double a_hi, a_lo, b_hi, b_lo;
    double pd = a*b;
    split(a,a_hi,a_lo);
    split(b, b_hi, b_lo);
    err = ((a_hi * b_hi - pd) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
    s = pd;
}


void dd_sum_dd(double x_h, double x_l, double y_h, double y_l, double &t_h, double &t_l){
    double s1,s2,t1,t2;
    d_sum_d(x_h,y_h,s1,s2);
    d_sum_d(x_l,y_l,t1,t2);
    s2+=t1;
    d_qsum_d(s1,s2,s1,s2);
    s2+=t2;
    d_qsum_d(s1,s2,s1,s2);
    t_h = s1;
    t_l = s2;

}


void dd_mul_dd(double x_h, double x_l, double y_h, double y_l, double &t_h, double &t_l){
    double p1, p2;
    d_mul_d(x_h,y_h,p1,p2);
    p2+=(x_h*y_l + x_l*y_h);
    d_qsum_d(p1,p2,p1,p2);
    t_h = p1;
    t_l = p2;
}


void nint_dd(double &x_h, double &x_l){
    double hi;
    nint(x_h,hi);
    double lo;
    if(hi == x_h){
        nint(x_l,lo);
        d_qsum_d(hi,lo,hi,lo);
    } else {
        lo = 0.0;
        if(abs(hi - x_h) == 0.5 && x_l < 0.0){
            hi-=1.0;
        }
    }
    x_h = hi;
    x_l = lo;
}


void dd_mul_d(double x_h, double x_l, double y, double &res_h, double &res_l){
    double p1, p2;
    d_mul_d(x_h, y, p1, p2);
    p2+=(x_l * y);
    d_qsum_d(p1,p2,p1,p2);
    res_h = p1;
    res_l = p2;
}


void dd_sub_dd(double x_h, double x_l, double y_h, double y_l, double &res_h, double &res_l){
    double s1,s2,t1,t2;
    d_sub_d(x_h,y_h,s1,s2);
    d_sub_d(x_l,y_l,t1,t2);
    s2+=t1;
    d_qsum_d(s1,s2,s1,s2);
    s2+=t2;
    d_qsum_d(s1,s2,s1,s2);
    res_h = s1;
    res_l = s2;
}


void dd_add_d(double x_h, double x_l, double b, double &res_h, double &res_l){
    double s1, s2;
    d_sum_d(x_h, b,s1,s2);
    s2 += x_l;
    d_qsum_d(s1,s2,s1,s2);
    res_h = s1;
    res_l = s2;
}


void dd_div_dd(double x_h, double x_l, double y_h, double y_l, double &res_h, double &res_l){
    double q1, q2, q3;
    q1 = x_h / y_h;

    double r_h, r_l;
    double temp_h, temp_l;

    dd_mul_d(y_h,y_l,q1,temp_h,temp_l);
    dd_sub_dd(x_h,x_l,temp_h,temp_l,r_h,r_l);

    q2 = r_h / y_h;
    dd_mul_d(y_h,y_l,q2,temp_h,temp_l);
    dd_sub_dd(r_h,r_l,temp_h,temp_l,r_h,r_l);

    q3 = r_h / y_h;
    d_qsum_d(q1,q2,q1,q2);

    dd_add_d(q1,q2,q3,res_h,res_l);
}



static const int n_inv_fact = 15;
static const double inv_fact[n_inv_fact][2] = {
        { 1.66666666666666657e-01,  9.25185853854297066e-18},
        { 4.16666666666666644e-02,  2.31296463463574266e-18},
        { 8.33333333333333322e-03,  1.15648231731787138e-19},
        { 1.38888888888888894e-03, -5.30054395437357706e-20},
        { 1.98412698412698413e-04,  1.72095582934207053e-22},
        { 2.48015873015873016e-05,  2.15119478667758816e-23},
        { 2.75573192239858925e-06, -1.85839327404647208e-22},
        { 2.75573192239858883e-07,  2.37677146222502973e-23},
        { 2.50521083854417202e-08, -1.44881407093591197e-24},
        { 2.08767569878681002e-09, -1.20734505911325997e-25},
        { 1.60590438368216133e-10,  1.25852945887520981e-26},
        { 1.14707455977297245e-11,  2.06555127528307454e-28},
        { 7.64716373181981641e-13,  7.03872877733453001e-30},
        { 4.77947733238738525e-14,  4.39920548583408126e-31},
        { 2.81145725434552060e-15,  1.65088427308614326e-31}
};

void sin_taylor(double a_h, double a_l, double &res_h, double &res_l){
    if(a_h == 0.0 && a_l == 0.0)
    {
        res_h = 0.0;
        res_l = 0.0;
        return;
    }
    double r_h=a_h, r_l=a_l;
    res_h = a_h;
    res_l = a_l;
    double t_h,t_l;
    double x_h,x_l;
    dd_mul_dd(a_h,a_l,a_h,a_l,x_h,x_l);
    int i = 0;
    int j = 1;
    do{
        dd_mul_dd(r_h,r_l,x_h,x_l,r_h,r_l);
        double fac_h=inv_fact[i][0];
        double fac_l=inv_fact[i][1];
        dd_mul_dd(r_h,r_l,fac_h,fac_l,t_h,t_l);
        if(j % 2 !=0){
            dd_sub_dd(res_h,res_l,t_h,t_l,res_h,res_l);
         }
        else
            dd_sum_dd(res_h,res_l,t_h,t_l,res_h,res_l);
        i+=2;
        j++;
    } while(i < n_inv_fact);
}

void cos_taylor(double a_h, double a_l, double &res_h, double &res_l){
    if(a_h == 0.0 && a_l == 0.0)
    {
        res_h = 1.0;
        res_l = 0.0;
        return;
    }

    double t_h,t_l;
    double x_h,x_l;
    dd_mul_dd(a_h,a_l,a_h,a_l,x_h,x_l);
    x_h = -x_h;
    x_l = -x_l;
    double r_h=x_h, r_l=x_l;
    res_h = r_h*0.5;
    res_l = r_l*0.5;
    dd_add_d(res_h,res_l,1.0,res_h,res_l);
    int i = 1;
    do {
        dd_mul_dd(r_h,r_l,x_h,x_l,r_h,r_l);
        double fac_h=inv_fact[i][0];
        double fac_l=inv_fact[i][1];
        dd_mul_dd(r_h,r_l,fac_h,fac_l,t_h,t_l);
        dd_sum_dd(res_h,res_l,t_h,t_l,res_h,res_l);
        i+=2;
    } while (i < n_inv_fact);

}

bool dd_less_or_equal_dd(double x_h, double x_l, double y_h, double y_l){
    return (x_h < y_h || (x_h == y_h && x_l <= y_l));
}

void sin(double a_h, double a_l, double &res_h, double &res_l){
    if(a_h == 0.0 && a_l == 0.0)
    {
        res_h = 0.0;
        res_l = 0.0;
        return;
    }

    double pi2_h = 1.570796326794896558e+00;
    double pi2_l = 6.123233995736766036e-17;

    if(dd_less_or_equal_dd(a_h,a_l,pi2_h,pi2_l) && dd_less_or_equal_dd(-pi2_h,-pi2_l,a_h,a_l)){
        sin_taylor(a_h,a_l,res_h,res_l);
        return;
    }

    double _2pi_h = 6.283185307179586232e+00;
    double _2pi_l = 2.449293598294706414e-16;

    double z_h,z_l;
    double temp_h,temp_l;
    dd_div_dd(a_h,a_l,_2pi_h,_2pi_l,z_h,z_l);
    nint_dd(z_h,z_l);
    dd_mul_dd(_2pi_h,_2pi_l,z_h,z_l,temp_h,temp_l);

    dd_sub_dd(a_h,a_l,temp_h,temp_l,res_h,res_l);

    double t_h,t_l;
    double q = floor(res_h / pi2_h + 0.5);
    dd_mul_d(pi2_h,pi2_l,q,temp_h,temp_l);
    dd_sub_dd(res_h,res_l,temp_h,temp_l,t_h,t_l);

    int j = static_cast<int>(q);

    if (j < -2 || j > 2) {
        cout << "Can't reduce modulo pi/2 \n";
        return;
    }

    switch (j) {
        case 0: {
            sin_taylor(t_h,t_l,res_h,res_l);
            return;
        }

        case 1:{
            cos_taylor(t_h,t_l,res_h,res_l);
            return;
        }
        case -1:{
            cos_taylor(t_h,t_l,res_h,res_l);
            res_h = -res_h;
            res_l = -res_l;
            return;
        }
        default:{
            sin_taylor(t_h,t_l,res_h,res_l);
            res_h = -res_h;
            res_l = -res_l;
            return;
        }
    }

}

/*
 * 2*pi 6.283185307179586232e+00 2.449293598294706414e-16
 * pi 3.141592653589793116e+00 1.224646799147353207e-16
 * pi/2 1.570796326794896558e+00 6.123233995736766036e-17
 * pi/4 7.853981633974482790e-01 3.061616997868383018e-17
 * 3/4*pi 2.356194490192344837e+00 9.1848509936051484375e-17
 * pi/16 1.963495408493620697e-01 7.654042494670957545e-18
 */

/* sin(pi/16) = 1.950903220161282758e-01, -7.991079068461731263e-18
 * sin(pi/8) = 3.826834323650897818e-01, -1.005077269646158761e-17
 * sin(3/16pi) = 5.555702330196021776e-01,  4.709410940561676821e-17
 * sin(pi/4) = 7.071067811865475727e-01, -4.833646656726456726e-17
 */

int main() {


    double some_pi_h = 0.0;
    double some_pi_l = 1e-20;

    double pi_h = 3.141592653589793116e+00;
    double pi_l = 1.224646799147353207e-16;

    //dd_mul_d(pi_h, pi_l, 100, pi_h, pi_l);
    //dd_sum_dd(some_pi_h,some_pi_l,pi_h,pi_l,some_pi_h,some_pi_l);

    double res_sin_h,res_sin_l;
    sin(some_pi_h,some_pi_l,res_sin_h,res_sin_l);
    printf("res_sin_h: %.18le \n",res_sin_h);

    int res_h_power;
    double res_h_mantissa = frexp(res_sin_h, &res_h_power);
    printf("mantissa_h : %f, power_h: %d, sign_h: %s\n\n",res_h_mantissa, res_h_power, signbit(res_sin_h) ? "-" : "+");

    printf("res_sin_l: %.18le \n",res_sin_l);

    int res_l_power;
    double res_l_mantissa = frexp(res_sin_l, &res_l_power);
    printf("mantissa_l : %f, power_l: %d, sign_l: %s\n\n",res_l_mantissa, res_l_power, signbit(res_sin_l) ? "-" : "+");

    printf("sin:    %.64le \n",res_sin_h + res_sin_l);

    return 0;
}
