#include <iostream>
#include <cfenv>
#include <cmath>
#include <iomanip>
#define SHIFT_POW 27
using namespace std;

void fast2Sum(double a, double b, double *s, double *t) {
    double dum;
    double z;
/* Branching below may hinder performance */
/* Suppress if we know in advance that a >= b */
    if (b > a) { dum = a; a = b; b = dum; }
    *s= a+b;
    z = *s - a;
    *t= b-z;
}

/*
void TwoSum(double a, double b, double *s, double *t){
    *s = a + b;
    double a_ = *s - b;
    double b_ = *s - a_;
    double del_a = a - a_;
    double del_b = b - b_;
    *t = del_a + del_b;
}
 */

void dekkerSum(double x_h, double x_l, double y_h, double y_l, double *t_h, double *t_l){
    double r_h;
    double r_l;
    double s;
    if(fabs(x_h) > fabs(y_h)){
        fast2Sum(x_h,y_h,&r_h,&r_l);
        s = (r_l + y_l) + x_l;
    }
    else {
        fast2Sum(y_h,x_h,&r_h,&r_l);
        s = (r_l + x_l) + y_l;
    }
    fast2Sum(r_h,s,t_h,t_l);
}

void veltkampSplit(double x, int sp, double *x_high, double *x_low) {
    unsigned long C = (1UL << sp) + 1;
    double gamma = C * x;
    double delta = x - gamma;
    *x_high = gamma + delta;
    *x_low = x - *x_high;
}

void dekkerMult(double x, double y, double *r_1, double *r_2) {
    double x_high, x_low; double y_high, y_low; double t_1;
    double t_2;
    double t_3;
    veltkampSplit(x, SHIFT_POW, &x_high, &x_low);
    veltkampSplit(y, SHIFT_POW, &y_high, &y_low);
    *r_1 = x * y;
    t_1  = -*r_1 + x_high * y_high;
    t_2  =   t_1 + x_high * y_low;
    t_3  =   t_2 + x_low  * y_high;
    *r_2 =   t_3 + x_low  * y_low;
}

void dekkerDWORDmult(double x_h, double x_l, double y_h, double y_l, double *t_h, double *t_l){
    double c_h,c_l;
    dekkerMult(x_h,y_h, &c_h, &c_l);
    double p1 = x_h*y_l;
    double p2 = x_l*y_h;
    c_l = c_l + (p1+p2);
    fast2Sum(c_h,c_l,t_h,t_l);
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

void sin_taylor(double a_h, double a_l, double *res_h, double *res_l){
    if(a_h == 0.0)
    {
        *res_h = 0.0;
        *res_l = 0.0;
        return;
    }
    double r_h=a_h, r_l=a_l;
    *res_h = a_h;
    *res_l = a_l;
    double t_h,t_l;
    int i = 0;
    int j = 1;
    do{
        dekkerDWORDmult(r_h,r_l,a_h,a_l,&r_h,&r_l);
        dekkerDWORDmult(r_h,r_l,a_h,a_l,&r_h,&r_l);
        double fac_h=inv_fact[i][0];
        double fac_l=inv_fact[i][1];
        dekkerDWORDmult(r_h,r_l,fac_h,fac_l,&t_h,&t_l);
        if(j % 2 !=0){
            t_h*=-1;
            t_l*=-1;
        }
        dekkerSum(*res_h,*res_l,t_h,t_l,res_h,res_l);
        i+=2;
        j++;
    } while(i < n_inv_fact);
}

/*
 * 2*pi 6.283185307179586232e+00 2.449293598294706414e-16
 * pi 3.141592653589793116e+00 1.224646799147353207e-16
 * pi/2 1.570796326794896558e+00 6.123233995736766036e-17
 * pi/4 7.853981633974482790e-01 3.061616997868383018e-17
 * 3/4*pi 2.356194490192344837e+00 9.1848509936051484375e-17
 * pi/16 1.963495408493620697e-01 7.654042494670957545e-18
 */

/* Table of sin(k * pi/16) and cos(k * pi/16). */
static const double sin_table [4][2] = {
        {1.950903220161282758e-01, -7.991079068461731263e-18},
        {3.826834323650897818e-01, -1.005077269646158761e-17},
        {5.555702330196021776e-01,  4.709410940561676821e-17},
        {7.071067811865475727e-01, -4.833646656726456726e-17}
};

int main() {
    std::fesetround(FE_TONEAREST);

    double some_pi_h = 1.963495408493620697e-01;
    double some_pi_l = 7.654042494670957545e-18;

    double res_sin_h,res_sin_l;
    sin_taylor(some_pi_h,some_pi_l,&res_sin_h,&res_sin_l);
    printf("res_sin_h: %.18le \n",res_sin_h);

    int res_h_power;
    double res_h_mantissa = frexp(res_sin_h, &res_h_power);
    printf("mantissa_h : %f, power_h: %d, sign_h: %s\n\n",res_h_mantissa, res_h_power, signbit(res_sin_h) ? "-" : "+");

    printf("res_sin_l: %.18le \n",res_sin_l);

    int res_l_power;
    double res_l_mantissa = frexp(res_sin_l, &res_l_power);
    printf("mantissa_l : %f, power_l: %d, sign_l: %s\n\n",res_l_mantissa, res_l_power, signbit(res_sin_l) ? "-" : "+");

    return 0;
}