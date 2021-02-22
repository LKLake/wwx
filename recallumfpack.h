//Data:2013-2-24
//修改了Ai Ax的类型 利用最大维数n*n来保存，可以调用正确结果 不过不知系统随机分配的值 函数没有用
//上面这句话后面修改了程序，利用二次循环，将Ai Ax的维数确定
//Data:2013-2-26
//调用UMFPACK包来实现求解方程组
//UMFPACK采用CSC（列压缩存储） matlab中的接口为A/b
#include <stdio.h>
#include <math.h>
#include "umffpack.h"
//传递的四个参数A b x n -----Data:2013-02-27
//意思为Ax=b n为矩阵维数
//umf(SparseMatAi,SparseMatAx,SparseMatAp,B,testx,(N+1)*(M+1));
void umff(int *Ai,double *Ax,int *Ap,double  *b,double *x,int n)
{
	//2013-4-17
	double *nul =NULL ;// 此处修改了定义的形式。之前是double *null=(double *)NULL 这种有可能编译器将null与NULL等同起来。重命名改为nul
    void *Symbolic, *Numeric ;

	umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, nul, nul);
    umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, nul, nul) ; 
    umfpack_di_free_symbolic (&Symbolic);
    umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, nul, nul);
    umfpack_di_free_numeric (&Numeric) ;
	//return 0;
}