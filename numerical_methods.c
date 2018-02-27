#include <stdio.h>
#include <math.h>
#define f(x) (x*x*x)-x-1
#define g(x) x-((x*x*x)-x-1)/((3*x*x)-1)
#define fdiff(x) (3*x*x)-1
#define NO 10
#define TOL 1e-3

void bisectionIteration(double a, double b){
    int i=1;
    double FA=f(a);
    printf("n       a          b          p         f(p)\n");
    while(i<=NO){
        double p=((a+b)/2);
        double FP=f(p);
        printf("%d   %f   %f   %f   %f   \n",i,a,b,p,FP);
        if(FP==0 || ((b-a)/2)<TOL){
            printf("\t%f is a root of function\n",p);
            break;
        }
        i=i+1;
        if(FA*FP>0){
            a=p;
            FA=FP;
        }
        else{
            b=p;
        }
    }
    printf("Method failed after %d iteration\n\n\n",i);
}

void fixedPointIteration(double p0){
    int i=1;
    printf("n       p0         p\n");
    while(i<=NO){
        double p=g(p0);
        printf("%d   %f   %f   \n",i,p0,p);
        if(fabs(p-p0)<TOL){
            printf("\t%f is a root of function\n",p);
            break;
        }
        i=i+1;
        p0=p;
    }
    printf("Method failed after %d iteration\n\n\n",i);
}

void newtonRaphsonIteration(double p0){
    int i=1;
    printf("n       p         f(p)\n");
    while(i<=NO){
        double p=p0-((f(p0))/(fdiff(p0)));
        printf("%d   %f   %f   \n",i,p,f(p));
        if(fabs(p-p0)<TOL){
            printf("\t%f is a root of functions\n",p);
            break;
        }
        i=i+1;
        p0=p;
    }
    printf("Method failed after %d iteration\n\n\n",i);
}

void secantIteration(double p0,double p1){
    int i=2;
    double q0=f(p0);
    double q1=f(p1);
    printf("n       p0         p1         p\n");
    while(i<=NO){
        double p=p1-(q1*(p1-p0))/(q1-q0);
        printf("%d   %f   %f   %f   \n",i-1,p0,p1,p);
        if(fabs(p-p1)<TOL){
            printf("\t%f is a root of function\n",p);
            break;
        }
        i=i+1;
        p0=p1;
        q0=q1;
        p1=p;
        q1=f(p);
    }
    printf("Method failed after %d iteration\n\n\n",i-1);
}

void regulaPalsiIteration(double p0,double p1){
    int i=2;
    double q0=f(p0);
    double q1=f(p1);
    printf("n       p0         p1         p       \n");
    while(i<=NO){
        double p=p1-(q1*(p1-p0))/(q1-q0);
        double q=f(p);
        printf("%d   %f   %f   %f   \n",i-1,p0,p1,p);
        if(fabs(p-p1)<TOL){
            printf("\t%f is a root of function\n",p);
            break;
        }
        i=i+1;
        if((q*q1)<0){
            p0=p1;
            q0=q1;
        }
        p1=p;
        q1=q;
    }
    printf("Method failed after %d iteration\n\n\n",i-1);
}

int main(){
    printf("\t With Bisection Iteration\n");
    bisectionIteration(1,2);
    printf("\t With Fixed Point Iteration\n");
    fixedPointIteration(1.5);
    printf("\t With Newton-Raphson Iteration\n");
    newtonRaphsonIteration(1.5);
    printf("\t With Secant Iteration\n");
    secantIteration(1,1.5);
    printf("\t With Regula-Palsi Iteration\n");
    regulaPalsiIteration(1,1.5);

    return 0;
}
/*Emre ERÝNÇ*/
