#include "gmp.h"
#include "mpfr.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NN 8192

struct cmplx {
 mpfr_t re;
 mpfr_t im;
};

static int bitr(int j, int nu){
static int j1,i,j2,bit;
j1=j;
bit=0;
for(i=1;i<=nu;i++){
 j2=j1/2;
 bit=2*bit+(j1-2*j2);
 j1=j2;
}
return bit;
}

static void fft(int fl, int sign, mpfr_t xr[],mpfr_t xi[], int nu, mpfr_t pi, mp_rnd_t rmode){
//   tiny FFT program
//  sign: FFT, when sign>0; Inverse FFT, when sign<0
//  xr,xi: input & output; xr=array real part, xi=array of imaginary part
//  nu: N=^nu; N is length of xr and xi
//  pi: ratio of diameter vs. circumference
//  rmode: rounding mode
static int n2,nu1,l,i,k,k1,k1n2,p,n, b;
//static int flag=0;
static mp_prec_t prec;
static mpfr_t arg,c[NN],s[NN],tr,ti,w,d,f,pi2,w1,w2;
//  get the precision of variables
prec=mpfr_get_prec(pi);
n=pow(2,nu);
if(fl==0){
 mpfr_init2(arg,prec);
 for(i=0;i<n;i++){mpfr_init2(c[i],prec);mpfr_init2(s[i],prec);}
 mpfr_init2(tr,prec);mpfr_init2(ti,prec);mpfr_init2(d,prec);
 mpfr_init2(w,prec);mpfr_init2(f,prec);mpfr_init2(pi2,prec);
 mpfr_init2(w1,prec);mpfr_init2(w2,prec);
mpfr_set_d(pi2,2.0,rmode);
mpfr_set_si(w,n,rmode);
mpfr_div(pi2,pi2,w,rmode);
mpfr_mul(pi2,pi2,pi,rmode);
}

if(sign<0) for(i=0;i<n;i++) mpfr_neg(xi[i],xi[i],rmode);
n2=n/2;
nu1=nu-1;
k=0;

for(l=1;l<=nu;l++){
label1: for(i=1;i<=n2;i++){
  p=pow(2,nu1); p=k/p;p=bitr(p,nu);
  if(fl==0){
    mpfr_set_si(w,p,rmode);
    mpfr_mul(arg,pi2,w,rmode);
    mpfr_sin_cos(s[p],c[p],arg,rmode);
  }
  k1=k+1;
  k1n2=k1+n2;
  mpfr_mul(tr,xr[k1n2-1],c[p],rmode);
  mpfr_fma(tr,xi[k1n2-1],s[p],tr,rmode);
  mpfr_mul(ti,xi[k1n2-1],c[p],rmode);
  mpfr_mul(w1,xr[k1n2-1],s[p],rmode);
  mpfr_sub(ti,ti,w1,rmode);
  mpfr_sub(xr[k1n2-1],xr[k1-1],tr,rmode);
  mpfr_sub(xi[k1n2-1],xi[k1-1],ti,rmode);
  mpfr_add(xr[k1-1],xr[k1-1],tr,rmode);
  mpfr_add(xi[k1-1],xi[k1-1],ti,rmode);
  k=k+1;
 }
 k=k+n2;
 if(k<n) goto label1;
 k=0;
 nu1=nu1-1;
 n2=n2/2;
}
for(k=1;k<=n;k++){
 i=bitr(k-1,nu)+1;
 if(i<=k) goto label2;
 mpfr_swap(xr[k-1],xr[i-1]);
 mpfr_swap(xi[k-1],xi[i-1]);
label2: continue;
}
/*
mpfr_clear(arg);m
for(i=0;i<n;i++){pfr_clear(c[i]);mpfr_clear(s[i]);}
mpfr_clear(tr);mpfr_clear(ti);mpfr_clear(d);
mpfr_clear(w);mpfr_clear(f);mpfr_clear(pi2);
mpfr_clear(w1); mpfr_clear(w2);
*/
if(sign<0) for(i=0;i<n;i++) mpfr_neg(xi[i],xi[i],rmode);
return;
}
