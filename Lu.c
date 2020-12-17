#include <stdio.h>
#include "lapacke.h"
#include "cblas.h"
#include <math.h>
typedef unsigned long long u64;
double * allocation(u64 n)
{
    double *t;
    t=(double *) malloc(n*sizeof(double));
    return t;
}
void init(double* A,double* A1,u64 n,u64 m,double a)
{



    for (int j = 0; j < 2; ++j) {
        A[j]=(double)rand()/(double )(RAND_MAX/a);
        A1[j]=A[j];
        //A[j]=j;

        //A[i*n+j]=i*n+j;
    }
    for (int j = 3; j < m; ++j) {
        A[j]=0;
        A1[j]=A[j];
        }
    u64 deb=0;
    for (int i = 1; i < n-1; ++i) {
        for (int j = deb; j < deb+3; ++j) {
            A[i*m+j]=(double)rand()/(double )(RAND_MAX/a);
            A1[i*m+j]=A[i*m+j];
            //A[i*m+j]=i*m+j;
        }

        for (int j = deb+3; j < m; ++j) {
            A[i*m+j]=0;
            A1[i*m+j]=A[i*m+j];

        }
        deb++;

    }
    for (int j = m-2; j < m; ++j) {
        A[(n-1)*m+j]=(double)rand()/(double )(RAND_MAX/a);
        A1[(n-1)*m+j]=A[(n-1)*m+j];
        //A[i*n+j]=i*n+j;
    }
    for (int j = 0; j < m-2; ++j) {
        A[(n-1)*m+j]=0;

    }
}
void LU(double* A,double* L,double* U,u64 n,u64 m)
{
    A[m]=A[m]/A[0];
    A[m+1]=A[m+1]-A[m]*A[1];
    for (int i = 1; i < n-2; ++i) {
        A[(i+1)*m + i]=A[(i+1)*m + i]/A[(i)*m +i ];
        A[(i+1)*m + i+1]=A[(i+1)*m + i+1]-A[(i+1)*m + i]*A[(i)*m + i+1];
    }
    A[(n-1)*m + m-2]=A[(n-1)*m + m-2]/A[(n-2)*m + m-2];
    A[(n-1)*m + m-1]=A[(n-1)*m + m-1]-A[(n-1)*m + m-2]*A[(n-2)*m + m-1];



    for (int i = 0; i < n; ++i) {
        L[i*m+i]=1;
        for (int j = 0; j < i; ++j) {
         L[i*m+j]=A[i*m+j];
         U[i*m+j]=0;
        }

    }
    for (int i = 0; i < n; ++i) {
        U[i*m+i]=A[i*m+i];
        for (int j = i+1; j < m; ++j) {
            U[i*m+j]=A[i*m+j];
            L[i*m+j]=0;
        }
    }
}
void print(double* A,u64 n,u64 m)
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
void prod(double* A,double* L,double* U,u64 n,u64 m)
{
    for (int i = 0; i <n ; ++i) {
        for (int j = 0; j <m ; ++j) {
            A[i*m+j]=0;
            for (int k = 0; k < m; ++k) {
                A[i*m+j]=A[i*m+j]+L[i*m+k]*U[k*m+j];
            }
        }
    }
}
int main(int argc, char *argv[]) {
    u64 n=strtol(argv[1],NULL,10);
    u64 m=strtol(argv[2],NULL,10);
    double a=strtod(argv[3],NULL);
    double alpha=strtod(argv[5],NULL);
    double beta=strtod(argv[6],NULL);
    CBLAS_TRANSPOSE type=(CBLAS_TRANSPOSE)argv[4];
    double *A,*L,*U,*p,*A1;
    A=allocation(m*n);
    A1=allocation(m*n);

    L=allocation(m*(n));
    U=allocation(m*(n));
    p=allocation(n);
    init(A,A1,n,m,a);

    //set_GB_operator_colMajor_poisson1D(AB,m,n,0);
    //Gba_c(A,AB,n,m);

    //print(A,n,m);
    //cblas_dgbmv(CblasColMajor,CblasNoTrans,n,m,(m-1),(n-1),alpha,AB,n+m-1,x,1,beta,y,1);
    LU(A,L,U,n,m);
    //print(L,n,m);
    //print(U,n,m);
    prod(A,L,U,n,m);
    //print(A,n,m);
    double temp;
    temp = cblas_ddot(n*m,A, 1,A,1);
    temp = sqrt(temp);
    //print(A1,n,m);
    cblas_daxpy(n*m, -1, A, 1, A1, 1);
    double relres;
    //print(A1,n,m);
    relres = cblas_ddot(n*m, A1, 1, A1,1);
    relres = sqrt(relres);
    relres = relres / temp;

    printf("\nThe relative residual error is relres = %e \n",relres);

    //LAPACKE_dgetrf(LAPACK_ROW_MAJOR,m,n,A1,m,p);
    //print(A1,n,m);
    //printf("%ld %ld",n,m);
    //print(AB,x,y,n,m,1);
    //cblas_dgemv(CblasRowMajor,CblasNoTrans,n,m,alpha,A,m,x,1,beta,y1,1);
    //print(AB,x,y1,n,m,1);


    return 0;
}
