#include <stdio.h>
#include "lapacke.h"
#include "cblas.h"
#include <math.h>
#include "rdtsc.h"
typedef unsigned long long u64;

void set_GB_operator_colMajor_poisson1D(double* AB, int lab, int la, int kv){
    int ii, jj, kk;
    for (jj=0;jj<(la);jj++){
        kk = jj*(lab);
        AB[3*jj]=-1.0;
        AB[3*jj+1]=2.0;
        AB[3*jj+2]=-1.0;
    }
    AB[0]=0.0;

    AB[(lab)*(la)-1]=0.0;
}

void print(double* A,double* x,double* b,u64 n,u64 m,int test)
{
    if(test==0)
    {
    printf("A=\n");
    for (int i = 0; i < n; ++i) {
        printf("[");
        for (int j = 0; j < m-1; ++j) {
            printf("%f,",A[i*n+j]);
        }
        printf("%f]\n",A[i*n+m-1]);
    }
    }
    printf("x=[\n");
    for (int i = 0; i < n-1; ++i) {
        printf("%f,",x[i]);
    }
    printf("%f]\n",x[n-1]);
    printf("b=[\n");
    for (int i = 0; i < n-1; ++i) {
        printf("%f,",b[i]);
    }
    printf("%f]\n",b[n-1]);

}
double * allocation(u64 n)
{
    double *t;
    t=(double *) malloc(n*sizeof(double));
    return t;
}

void init(double* A,double* x,double* b,double* sol,u64 n,double a)
{


    //x[0]=(double)rand()/(double )(RAND_MAX/a);
    x[0]=(double)rand()/(double )(RAND_MAX/a);
    sol[0]=0;
    b[0]=0;
    for (int j = 0; j < 2; ++j) {
        //A[j]=(double)rand()/(double )(RAND_MAX/a);
        A[j]=-1;
    }
    A[0]=abs(A[1])*2;

    for (int j = 3; j < n; ++j) {
        A[j]=0;

    }

    u64 deb=0;
    for (int i = 1; i < n-1; ++i) {
        for (int j = deb; j < deb+3; ++j) {
            //A[i*n+j]=(double)rand()/(double )(RAND_MAX/a);
            A[i*n+j]=-1;
        }
        b[i]=0;
        x[i]=(double)rand()/(double )(RAND_MAX/a);
        sol[i]=0;
        //A[i*n+i]=abs(A[i*n + i-1]+A[i*n + i+1])+1;
        //x[i]=i+1;
        A[i*n+i]=abs(A[i*n + i-1]+A[i*n + i+1]);



        for (int j = deb+3; j < n; ++j) {
            A[i*n+j]=0;
        }
        deb++;


    }


    for (int j = n-2; j < n; ++j) {
        //A[(n-1)*n+j]=(double)rand()/(double )(RAND_MAX/a);
        A[(n-1)*n+j]=-1;
    }
    A[(n-1)*n+n-1]=abs(A[(n-1)*n+n-1])*2;
    for (int j = 0; j < n-2; ++j) {
        A[(n-1)*n+j]=0;

    }
    x[n-1]=(double)rand()/(double )(RAND_MAX/a);
    //x[n-1]=n;
    b[n-1]=0;
    sol[n-1]=0;
    b[0]=0;
}
int gauss_seidel(double* A,double* AB,double* x,double *b,double* resles,u64 n,double eps,u64 itmax)
{
    double ra,*d,*b1,r,*temp;
    b1=allocation(n);

    d=allocation(n);
    temp=allocation(n);
    for (int i = 0; i < n; ++i) {
        d[i]=0;
        temp[i]=0;
    }
    r=1;
    cblas_daxpy(n,0,x,1,d,1);
    ra= cblas_ddot(n*3,AB, 1,AB,1);
    ra=sqrt(ra);
    cblas_dcopy(n,b,1,b1,1);

    u64 it=0;
    while (r>eps && it<itmax)
    {

        temp[0]=(1/A[0])*(-A[1]*x[1]+b[0]);
        for (int i = 1; i < n-1; ++i) {
            temp[i]=(1/A[i*n + i])*(-A[i*n + i-1]*temp[i-1]-A[i*n + i+1]*x[i+1]+b[i]);
        }


        temp[n-1]=(1/A[(n-1)*n+n-1])*(-A[n*(n-1)+(n-2)]*temp[n-2]+b[n-1]);
        cblas_dgbmv(CblasColMajor,CblasNoTrans,n,n,1,1,-1,AB,3,temp,1,1,b,1);
        //cblas_dgbmv(CblasRowMajor,CblasNoTrans,n,n,1,1,-1,AB,n,temp,1,1,b,1);
        r = cblas_ddot(n,b, 1,b,1);

        r = sqrt(r)/ra;
        cblas_dcopy(n,temp,1,x,1);
        cblas_dcopy(n,d,1,temp,1);
        cblas_dcopy(n,b1,1,b,1);
        it++;
        resles[it]=r;

    }

    return it;

}
int jacobi(double* A,double* AB,double* x,double *b,double* resles,u64 n,double eps,u64 itmax)
{
    double ra,*d,*b1,r,*temp;
    b1=allocation(n);

    d=allocation(n);
    temp=allocation(n);
    for (int i = 0; i < n; ++i) {
        d[i]=0;
        temp[i]=0;
    }
    r=1;
    cblas_daxpy(n,0,x,1,d,1);
    ra= cblas_ddot(n*3,AB, 1,AB,1);
    ra=sqrt(ra);
    cblas_dcopy(n,b,1,b1,1);

    int it=0;
    while (r>eps && it<itmax)
    {

        for (int i = 1; i < n-1; ++i) {
            temp[i]=(1/A[i*n + i])*(-A[i*n + i-1]*x[i-1]-A[i*n + i+1]*x[i+1]+b[i]);
        }
        temp[0]=(1/A[0])*(-A[1]*x[1]+b[0]);

        temp[n-1]=(1/A[(n-1)*n+n-1])*(-A[n*(n-1)+(n-2)]*x[n-2]+b[n-1]);
        cblas_dgbmv(CblasColMajor,CblasNoTrans,n,n,1,1,-1,AB,3,temp,1,1,b,1);
        //cblas_dgbmv(CblasRowMajor,CblasNoTrans,n,n,1,1,-1,AB,n,temp,1,1,b,1);
        //printf(" %f %f \n",r,temp[1]);
        r = cblas_ddot(n,b, 1,b,1);

        r = sqrt(r)/ra;
        cblas_dcopy(n,temp,1,x,1);
        //cblas_dcopy(n,d,1,temp,1);
        cblas_dcopy(n,b1,1,b,1);
        resles[it]=r;
        it++;

    }


    return it;
}
void test_rdtsc_perf(double* A,double* AB,double* sol,double* b,double* resles,double* n)
{
    double* t;
    t=allocation(n);
    cblas_dcopy(n,sol,1,t,1);

    unsigned long long d,f,s;
    printf("n:perf \n");
    double eps=0.5;
    for (int i = 0; i <10; i+=1) {
        s=0;
        d=rdtsc();
        for (int j = 0; j <40 ; ++j) {
            jacobi(A,AB,sol,b,resles,n,eps,10000);

            cblas_dcopy(n,t,1,sol,1);


        }
        f=rdtsc();
        s=s+(f-d);
        s=s/40;
        printf("%llu,",eps,s);
        eps=eps/10;


    }
}
void test_rdtsc_conv(double* A,double* AB,double* sol,double* b,double* resles,double* n)
{
    double* t;
    t=allocation(n);
    cblas_dcopy(n,sol,1,t,1);

    unsigned long long d,f,s;
    printf("n:perf \n");
    double eps=0.000005;
        s=gauss_seidel(A,AB,sol,b,resles,n,eps,1000000);
        printf("gauss%llu",s);
        cblas_dcopy(n,t,1,sol,1);
        s=jacobi(A,AB,sol,b,resles,n,eps,1000000);
        printf("jacobi %llu;",s);


    eps=eps/10;


    }



int main(int argc, char *argv[]) {
    u64 n=strtol(argv[1],NULL,10);
    double a=strtod(argv[2],NULL);
    double *A,*AB,*x,*b,*resles,*sol;
    double epsilon=0.000000005;
    u64 itmax=1000000,it;
    A=allocation(n*n);
    x=allocation(n);
    resles=allocation(itmax);
    b=allocation(n);
    sol=allocation(n);
    init(A,x,b,sol,n,a);

    AB=allocation(n*3);
    set_GB_operator_colMajor_poisson1D(AB,3,n,0);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,n,n,1,1,1,AB,3,x,1,1,b,1);
    //it=gauss_seidel(A,AB,sol,b,resles,n,epsilon,itmax);
    //it=jacobi(A,AB,sol,b,resles,n,epsilon,itmax);
    test_rdtsc_conv(A,AB,sol,b,resles,n);
    //printf("%i ",it);
    //print(A,x,b,n,n,1);
    //print(A,sol,b,n,n,1);



    return 0;
}