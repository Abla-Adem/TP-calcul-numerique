#include <stdio.h>
#include "lapacke.h"
#include "cblas.h"
typedef unsigned long long u64;
void Gba(double* A,double* AB,u64 n,u64 m)
{
    int ligne=0;
    int cpt=0;
    int col;
    for (int j = m-1; j >=0; j-=1) {

        for (int k = ligne+1; k <n ; ++k) {
            AB[cpt]=0;
            cpt++;
        }
        col=j;
        for (int i = 0; i < ligne+1; ++i) {
        AB[cpt]=A[i*n+col];
        col++;
        cpt++;
        }
        if(ligne!=(n-1))
        {
            ligne++;
        }


    }

    for (int i = 1; i < n; ++i) {
        col=0;
        ligne=i;
        for (int j = i; j <n ; ++j) {
          AB[cpt]=A[ligne*n+col];
          col++;
          ligne++;
          cpt++;
        }
        for (int k = 0; k < i; ++k) {
            AB[cpt]=0;
            cpt++;
        }
    }


}
void Gba_c(double* A,double* AB,u64 n,u64 m)
{
    int saut=m+n-1;
    int col_zero=n-1;
    int ligne=0;
    int ligne_A,col_A;

    for (int j = m-1; j >= 0; j-=1) {
        for (int i = 0; i < col_zero; ++i) {
            AB[ligne+i*saut]=0;
            }
        ligne_A=0;
        col_A=j;
        for (int i =col_zero; i<m; i++) {
            AB[ligne+i*saut]=A[ligne_A*m+col_A];
            col_A++;
            ligne_A++;
        }
        col_zero--;
        ligne++;
    }
    int col;
    for (int i = 1; i <n ; ++i) {
        col_A=0;
        col=0;
        for (int j = i; j <n ; ++j) {
            AB[ligne+saut*col]=A[col_A+j*m];
            col_A++;
            col++;
        }
        ligne++;
    }




}





void print(double* A,double* x,double *y,u64 n,u64 m,int b)
{
    if(b==0)
    {
    printf("A=\n");
    for (int i = 0; i < n; ++i) {
        printf("[");
        for (int j = 0; j < m-1; ++j) {
            printf("%f,",A[i*n+j]);
        }
        printf("%f]\n",A[i*n+m-1]);
    }

    printf("x=[\n");
    for (int i = 0; i < n-1; ++i) {
        printf("%f,",x[i]);
    }
    printf("%f]\n",y[n-1]);
    printf("y=[\n");
    for (int i = 0; i < n-1; ++i) {
        printf("%f,",y[i]);
    }
    printf("%f]\n",y[n-1]);

     }
    else
    {
        /*
        printf("AB=[\n");
        for (int i = 0; i < m*n-1; ++i) {

            if(i%7==0 & i!=0)
            {
                printf("\n",i);
            }
                printf("%f,",A[i]);


        }
        printf("%f]\n",A[(n+m-1)*n-1]);
        */
        printf("y=[\n");
        for (int i = 0; i < n-1; ++i) {
            printf("%f,",y[i]);
        }
        printf("%f]\n",y[n-1]);

    }


}
double * allocation(u64 n)
{
    double *t;
        t=(double *) malloc(n*sizeof(double));
    return t;
}
void init(double* A,double* x,double *y,double *y1,u64 n,u64 m,double a)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            A[i*n+j]=(double)rand()/(double )(RAND_MAX/a);
            //A[i*n+j]=i*n+j;
        }
        //x[i]=i;
        //y[i]=i;
        x[i]=(double)rand()/(double )(RAND_MAX/a);
        y[i]=(double)rand()/(double )(RAND_MAX/a);
        y1[i]=y[i];

    }
}
int main(int argc, char *argv[]) {
    u64 n=strtol(argv[1],NULL,10);
    u64 m=strtol(argv[2],NULL,10);
    double a=strtod(argv[3],NULL);
    double alpha=strtod(argv[5],NULL);
    double beta=strtod(argv[6],NULL);
    CBLAS_TRANSPOSE type=(CBLAS_TRANSPOSE)argv[4];
    double *A,*x,*y,*AB,*y1;
    A=allocation(m*n);
    AB=allocation((m+n-1)*n);
    x=allocation(n);
    y=allocation(n);
    y1=allocation(n);
    init(A,x,y,y1,n,m,a);
    //set_GB_operator_colMajor_poisson1D(AB,m,n,0);
    Gba_c(A,AB,n,m);
    print(A,x,y,n,m,0);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,n,m,(m-1),(n-1),alpha,AB,n+m-1,x,1,beta,y,1);

    //printf("%ld %ld",n,m);
    print(AB,x,y,n,m,1);
    cblas_dgemv(CblasRowMajor,CblasNoTrans,n,m,alpha,A,m,x,1,beta,y1,1);
    print(AB,x,y1,n,m,1);


    return 0;
}
