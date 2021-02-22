//2014.8.12完成，computers & fluids 33(2004)375-404
  //用抛物线（固体边界离散上连续三点）近似边界

#include "gjdnn.h"
//#include "gaus.h"
#include <iostream>
#include <math.h> 
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <cmath>

using namespace std;


void VelPlusTwo(double **u,double **X,double **v,double **inisc,double **inisvel,
				double **inifc,int N,int M,int NN,int MM,int N_sec,int M_sec,double delta_x,double delta_y,
				double ss,double delta_t,int ni,int nj)
{
    int i,j,k,nk=2,mm=3,nn=2;
    double PI=3.1415926,rho=1.0,usolid=0.0,vsolid=0.0,scorxn,scoryn,fcorx,fcory,d1,d2,s,s1,s2;	
	void (*Sol_X)(double **,double **,int, int); 
	/*将流体单元中心点的物理量通过光滑函数转换，求得固体节点的物理量*/
	double **D=new double *[ni];
	for(int k=0;k<ni;k++)
		D[k]= new double [nj];

	int Q=0;   //矩阵D中非零元素总和

	int *NumR=new int [ni];   //分别存储前1 - ni行非0元素个数总和，NumR[ni-1]=Q;

	for(i=0;i<ni;i++)
    {
		scorxn=inisc[i][0];
		scoryn=inisc[i][1];
		for(j=0;j<nj;j++)
		{
			fcorx=inifc[j][0];            
            fcory=inifc[j][1];            
            d1=fabs((fcorx-scorxn)/delta_x);
            if(d1>2.5) {s1=0.0;}
			else if(d1>1.5) {s1=1.0625-PI/64-3*d1/4+d1*d1/8+(d1-2)/16*(sqrt(-14+16*d1-4*d1*d1))+0.0625*asin(sqrt(2.0)*(d1-2));}
			else if(d1>0.5) {s1=0.25+(1-d1)/8*(sqrt(-2+8*d1-4*d1*d1))-0.125*asin(sqrt(2.0)*(d1-1));}
			else {s1=0.375+PI/32-d1*d1/4;}
		    d2=fabs((fcory-scoryn)/delta_y);
            if(d2>2.5) {s2=0.0;}
			else if(d2>1.5) {s2=1.0625-PI/64-3*d2/4+d2*d2/8+(d2-2)/16*(sqrt(-14+16*d2-4*d2*d2))+0.0625*asin(sqrt(2.0)*(d2-2));}
			else if(d2>0.5) {s2=0.25+(1-d2)/8*(sqrt(-2+8*d2-4*d2*d2))-0.125*asin(sqrt(2.0)*(d2-1));}
			else {s2=0.375+PI/32-d2*d2/4;} 
            s=s1*s2/delta_x/delta_y;
            D[i][j]=s;
			if(s!=0)  Q=Q+1;
		} 
		NumR[i]=Q;
	}
	cout<<"计算D[i][j]完成   Q="<<Q<<endl;

	int *NumL=new int [Q];   //分别存储Q个非0元素的列号
	int R=0;
	for(i=0;i<ni;i++)
    {
		for(j=0;j<nj;j++)
		{
			if(D[i][j]!=0)
			{
				NumL[R]=j;
				R=R+1;
			}
		} 
	}

	/* 求矩阵A，A[i][i]=A[i][i]+E[i][j]*C[j][i] */
	double **A=new double *[ni];
	for(int k=0;k<ni;k++)
		A[k]= new double [ni];
	//计算A矩阵，直接法，速度慢
	/*for(i=0;i<ni;i++)
	{
		for(j=0;j<ni;j++)
		{
			A[i][j]=0.0;
			for(k=0;k<nj;k++)
			{
				A[i][j]=A[i][j]+D[j][k]*D[i][k]*delta_s[j]*delta_t*delta_x*delta_y/rho;
			}
		}
	}*/

	//计算A矩阵，由于矩阵内部含有零元素，在矩阵相乘时避免零元素的乘运算

	for(i=0;i<ni;i++)   //由于[D][D]T得到的是对称矩阵，因此只需要上三角矩阵
	{
		for(j=i;j<ni;j++)
		{
			A[i][j]=0.0;
			
			int a=NumR[j],b,c=NumR[i],d;

			if(j==0) b=0;
			else b=NumR[j-1];

			if(i==0) d=0;
			else d=NumR[i-1];

			for(int l=b;l<a;l++)
			{
				int m=NumL[l];
				for(int n=d;n<c;n++)
				{
					int p=NumL[n];
					if(m==p)
						A[i][j]=A[i][j]+D[j][m]*D[i][p];
				}
			}
		}
	}
	for(j=0;j<ni-1;j++)   //由于[D][D]T得到的是对称矩阵，下三角元素是上三角元素的转置
	{
		for(i=j+1;i<ni;i++)
		{
			A[i][j]=A[j][i];
		}
	}
	for(i=0;i<ni;i++)  //矩阵[A]最终结果，乘上常系数和固体边界信息
	{
		for(j=0;j<ni;j++)
		{
			A[i][j]=A[i][j]*ss*delta_t*delta_x*delta_y/rho;
		}
	}
	cout<<"计算A[i][j]完成"<<endl;

	/* 求矩阵B */
	double **B=new double *[ni];
	for(int k=0;k<ni;k++)
		B[k]= new double [nk];

	for(i=0;i<ni;i++)
	{
		B[i][0]=0.0;
		B[i][1]=0.0;
		for(j=N_sec;j<N_sec+NN+1;j++)
		{
			for(k=M_sec;k<M_sec+MM+1;k++)
			{
				int l=(j-N_sec)*(MM+1)+k-M_sec;
				B[i][0]=B[i][0]+delta_x*delta_y*D[i][l]*u[j][k];
				B[i][1]=B[i][1]+delta_x*delta_y*D[i][l]*v[j][k];
			}
		}
		B[i][0]=inisvel[i][0]-B[i][0];     //固体节点速度：inisvel[i][k]，从主函数中传入
		B[i][1]=inisvel[i][1]-B[i][1];
	}
	cout<<"计算B[i][j]完成"<<endl;
	
	/* 计算固体节点力X[ni][1]*/
	Sol_X = gjdn;
	(*Sol_X)((double **)A,(double **)B,ni,nk);     //调用gidn函数,求线性代数方程组AX=B
  	      
    for(j=0;j<1;j++)                     //如何计算：还没弄清楚
		for(i=0;i<nk;i++)
		{X[j][i]=(B[ni-1][i]+2*B[j][i]+B[j+1][i])/4;}
	
	for(j=1;j<ni-1;j++)
		for (i=0;i<nk;i++)
		{X[j][i]=(B[j-1][i]+2*B[j][i]+B[j+1][i])/4;}

	for(j=ni-1;j<ni;j++)
		for(i=0;i<nk;i++)
		{X[j][i]=(B[j-1][i]+2*B[j][i]+B[0][i])/4;} 
    cout<<"计算固体节点力X[i][j]完成"<<endl;


    /* 将固体节点收到的力X[i][k]回代（光滑函数转换），求流体次区域节点修正速度delt_u[nj][nk]=plusfu[j][k]，
	   进而求得下一时间步速度u[nj][nk]=fu[j][k]=inifu[j][k]+plusfu[j][k]                                  */
	
	for(j=N_sec;j<N_sec+NN+1;j++)
	{
		for(k=M_sec;k<M_sec+MM+1;k++)
		{
			int l=(j-N_sec)*(MM+1)+k-M_sec;
			for(i=0;i<ni;i++)
			{
				u[j][k]=u[j][k]+delta_t/rho*X[i][0]*D[i][l]*ss;
				v[j][k]=v[j][k]+delta_t/rho*X[i][1]*D[i][l]*ss;
			}
		}
	}
    

	/* ******释放动态内存****** */
	for(i=0;i<ni;i++)
		delete []D[i];
	delete []D;

	for(i=0;i<ni;i++)
		delete []A[i];
	delete []A;

	for(i=0;i<ni;i++)
		delete []B[i];
	delete []B;

	delete []NumR;

	delete []NumL;

	//int *NumR=new int *[ni];   //分别存储前1 - ni行非0元素个数总和，NumR[ni-1]=Q;

	//int *NumL=new int *[Q];   //分别存储Q个非0元素的列号
	/* ******释放动态内存 结束****** */		
}