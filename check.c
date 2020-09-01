#include <stdlib.h>
#include <stdio.h>
#include "gmp.h"
#include "mpfr.h"

#define NN 64
#define np 6   /* np = log2(NN) */

void fft(int flag, int sign, mpfr_t xr[],mpfr_t xi[], int nu, mpfr_t pi, mp_rnd_t rmode);

int main(){
mpfr_t xr[NN],xi[NN],y[NN]; mpfr_t pi; mp_rnd_t rmode=GMP_RNDN;
mp_prec_t prec;
mpfr_t x,theta,w,dt;

unsigned int i;
unsigned int nu,n;
int sign, flag;
nu=np;n=NN;prec=107; /* 30 digits precision */

mpfr_init2(pi,prec);
for(i=0;i<=NN-1;i++)
{mpfr_init2(xr[i],prec);mpfr_init2(xi[i],prec);mpfr_init2(y[i],prec);}

mpfr_init2(x,prec); mpfr_init2(theta,prec); mpfr_init2(w,prec); 
mpfr_init2(dt,prec);
mpfr_const_pi(pi,rmode);
mpfr_set_si(w,NN,rmode);
mpfr_div(dt,pi,w,rmode);

for(i=1;i<=NN;i++){
mpfr_mul_ui(theta,dt,i,rmode);
mpfr_mul_ui(theta,theta,i,rmode);
mpfr_cos(xr[i-1],theta,rmode);
mpfr_sin(xi[i-1],theta,rmode);
}

sign=1;flag=0;
fft(flag, sign,xr,xi,nu,pi,rmode);
for(i=0;i<=NN-1;i++){
mpfr_mul(x,xr[i],xr[i],rmode);
mpfr_fma(x,xi[i],xi[i],x,rmode);
mpfr_sqrt(x,x,rmode);
printf("%d: ",i+1); mpfr_out_str(stdout, 10, 0, x, rmode); putchar('\n');
}
mpfr_clear(pi);
for(i=0;i<=NN-1;i++){mpfr_clear(xr[i]);mpfr_clear(xi[i]);mpfr_clear(y[i]);}
mpfr_clear(x); mpfr_clear(theta); mpfr_clear(w); mpfr_clear(dt);
return 0;
}
