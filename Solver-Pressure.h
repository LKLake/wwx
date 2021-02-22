#include <iomanip>
#include <fstream>
#include "matvec.h"
#include"GS.h"
#include"GMRES.h"
#include"ggje.h"
#include "callumfpack.h"
//2013-2-22 核查了A和B是否于Poisson离散出来形式一致
//加入了左边界条件为强制边界(Dirichlet条件)
//#define Temp(i,j) Temp[(i)*(N+1)*(M+1)+j] //使用宏定义将二维转换成一维
//2013-4-2 添加了umfpack计算压力
//2013-4-8 修改了在边界处离散的错误，以前为在每个边界按照一个方向差分 pian_u/pian_x=(u_n+1-u_n-1)/delta_x
//在边界处的差分都是由边界内部节点到外部节点
//2013-4-9 修改了pian_ustar/pian_x的迎风格式 按照具体物理模型进行离散
//2013-4-11 修改了压力输出的错误 以前导致有两条斜着的带子出现，错位导致
void Sol_P(int N,int M,double delta_x,double delta_y,double delta_t ,double **ustar ,double **vstar ,double **upie ,double **vpie,double **Pres,int *SparseMatAi,double *SparseMatAx,int *SparseMatAp)
{

    double *B =new double [(N+1)*(M+1)];
    for(int bk=0;bk<(N+1)*(M+1);bk++)
    {
        B[bk]=0;
        //             cout<<B[bk]<<" ";
    }

    /*---------------------------------------------------------------------------------------------*/

    B[0*(M+1)+0]=1/delta_t*(ustar[0+1][0]-ustar[0][0])/delta_x+1/delta_t*(vstar[0][0+1]-vstar[0][0])/delta_y+2/delta_x*(-(upie[0][0]-ustar[0][0])/delta_t)+2/delta_y*(-(vpie[0][0]-vstar[0][0])/delta_t);
    //                                    向后差分
    //cout<<"B[0]="<<B[0];
    /*i=0上的点（除去(0,0),(0,M)）*/
    for(int j=1;j<=M-1;j++)
    {
        int i=0;//第一行
        //A[i*(M+1)+j][i*(M+1)+j-1]=1/pow(delta_y,2);

        B[i*(M+1)+j]=1/delta_t*(ustar[i+1][j]-ustar[i][j])/delta_x+1/delta_t*(vstar[i][j+1]-vstar[i][j-1])/(2*delta_y)+2/delta_x*(-(upie[i][j]-ustar[i][j])/delta_t);
        //                                   向后差分和中心差分（能用中心差分就用中心差分）
    }

    B[0*(M+1)+M]=1/delta_t*(ustar[0+1][M]-ustar[0][M])/delta_x+1/delta_t*(vstar[0][M]-vstar[0][M-1])/delta_y+2/delta_x*(-(upie[0][M]-ustar[0][M])/delta_t)-2/delta_y*(-(vpie[0][M]-vstar[0][M])/delta_t);
    /*i=1~N-1,j=1~M-1 之间的离散方程*/
    for(int i=1;i<N;i++)
    {
        for(int j=1;j<M;j++)
        {

            B[i*(M+1)+j]=1/delta_t*((ustar[i+1][j]-ustar[i-1][j])/(2*delta_x)+(vstar[i][j+1]-vstar[i][j-1])/(2*delta_y));//对一阶导数采用中心差分
        }
    }
    /*------------------------------------------------------------------------------------------------------*/

    B[N*(M+1)+0]=1/delta_t*(ustar[N][0]-ustar[N-1][0])/delta_x+1/delta_t*(vstar[N][0+1]-vstar[N][0])/delta_y-2/delta_x*(-(upie[N][0]-ustar[N][0])/delta_t)+2/delta_y*(-(vpie[N][0]-vstar[N][0])/delta_t);//改变下边界符号
    /*i=N上的点（除去(0,0),(0,M)）*/
    for(int j=1;j<=M-1;j++)
    {
        int i=N;//第N行

        B[i*(M+1)+j]=1/delta_t*(ustar[i][j]-ustar[i-1][j])/delta_x+1/delta_t*(vstar[i][j+1]-vstar[i][j-1])/(2*delta_y)-2/delta_x*(-(upie[i][j]-ustar[i][j])/delta_t);
    }

    B[N*(M+1)+M]=1/delta_t*(ustar[N][M]-ustar[N-1][M])/delta_x+1/delta_t*(vstar[N][M]-vstar[N][M-1])/delta_y-2/delta_x*(-(upie[N][M]-ustar[N][M])/delta_t)-2/delta_y*(-(vpie[N][M]-vstar[N][M])/delta_t);

    /*---------------------------------------------------------------------------------------------------------------------*/
    /*j=M上的点（除去(0,M),(N,M)）*/
    for(int i=1;i<=N-1;i++)
    {
        int j=M;//第N行
        B[i*(M+1)+j]=1/delta_t*(ustar[i+1][j]-ustar[i-1][j])/(2*delta_x)+1/delta_t*(vstar[i][j]-vstar[i][j-1])/delta_y-2/delta_y*(-(vpie[i][j]-vstar[i][j])/delta_t);
    }
    /*---------------------------------------------------------------------------------------------------------------------*/
    /*j=0上的点（除去(0,0),(N,0)）*/
    for(int i=1;i<=N-1;i++)
    {
        int j=0;//第0行
        B[i*(M+1)+j]=1/delta_t*(ustar[i+1][j]-ustar[i-1][j])/(2*delta_x)+1/delta_t*(vstar[i][j+1]-vstar[i][j])/delta_y+2/delta_y*(-(vpie[i][j]-vstar[i][j])/delta_t);//改变下边界
    }
    /*---------------------------------------------------------------------------------------------------------------------*/

    //注：row为行号， col为列号
    //下面八行需要注意---需满足delta_x=delta_y
    for(int i=0;i<(N+1)*(M+1);i++)
    {

        B[i]=-1.0*pow(delta_x,2.0)*B[i];

    }

    /*-----------------------处理Ax=b中A的形式让其利用Gauss或者稀疏可解-----------------*/
    //????????????????????????????????
    for(int i=0;i<M+1;i++)
    {
        B[i]=10;//左边界赋值为常数10
    }

    double *testx =new double [(N+1)*(M+1)];
    for(int bk=0;bk<(N+1)*(M+1);bk++)
    {
        testx[bk]=0;
        //cout<<"B[bk]="<<B[bk]<<" ";
    }
    //cout<<endl;
    /*---------------------调用UMFPACK求解 A*testx=B ---------------------------*/
    umf(SparseMatAi,SparseMatAx,SparseMatAp,B,testx,(N+1)*(M+1));
    for(int row=0;row<=N;row++)
    {
        for(int col=0;col<=M;col++)
        {
            //Pres[row][col]=sol_P[(row-1)*M+col-1];  //Pres represents Pressure
            Pres[row][col]=testx[row*(M+1)+col];// 注意是M+1，非M
            //cout<<"testx[row*(M+1)+col]="<<testx[row*(M+1)+col]<<endl;
        }
    }

    /*撤销u*空间*/

    delete [] B;
    delete [] testx;

}