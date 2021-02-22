/*利用Gauss-Saidel迭代算法求解AX=b*/
using namespace std;
#define A(i,j)  A[(i)*N+j]
#define TOL 0.00001   //TOL represent the accuracy
void GS(double *A, int N,double b[],double x[])
{
	//double A[n][n]={{8,-3,2},{4,11,-1},{6,3,12}}; //A represent coefficient matrix
	//double b[n]={20,33,36};
	double *x0=new double [N];
	for(int i=0;i<N;i++)
		x0[i]=0;
	double *x1=new double [N];
	for(int i=0;i<N;i++)
		x1[i]=0;
	int counter=0; //counter represent the step of Computing
	double *err_x=new double [N];
	for(int i=0;i<N;i++)
		err_x[i]=0;
	double max_x=0;
	do
	{
		for (int i=0;i<N;i++)
		{
			double s=0;
			for(int j=0;j<N;j++)
			{
				if(j!=i)
				{
					s=s+A(i,j)*x0[j];
				}
			}
			x1[i]=(b[i]-s)/A(i,i);
			err_x[i]=x1[i]-x0[i];
			max_x=fabs(err_x[0]);
			if(fabs(err_x[i])>max_x)
			{
				max_x=fabs(err_x[i]);
			}
			x0[i]=x1[i];
		} 
    	counter++;
	}
	while(max_x>TOL);
	for (int i=0;i<N;i++)
		printf("%lf ",x1[i]);
	printf("\n");
	printf("%d",counter);
	delete [] x0;
	delete [] x1;
	delete [] err_x;
}
