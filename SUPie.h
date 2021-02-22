using namespace std;

void SolUPie(double **ustar,double **upie,double **Pres,int N,int M )
{
	double Length=15.0;   //求解区域的长
	double Wide=10.0;      //求解区域的宽
	double delta_x=Length/N;
	double delta_y=Wide/M;
	double delta_t=0.001;
	//double Re=10000;
    double uboundary=1;
	/*for(int row=0;row<(N+1);row++)
	{
		for(int col=0;col<(M+1);col++)
		{
			
			cout<<u[row][col]<<" ";
		}
		cout<<endl;
	}*/
	/*---------------------------------------i=0-------------------------------------------------------*/
	for(int i=0;i<1;i++)
	{
		for(int j=0;j<1;j++)
		{
			upie[i][j]=uboundary;
		}
	}
	for(int i=0;i<1;i++)
	{
		for(int j=1;j<M;j++)
		{
			upie[i][j]=uboundary;
		}
	}
	for(int i=0;i<1;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			upie[i][j]=uboundary;
		}
	}
	/*-----------------------------------------------------------------------------------------------------------------*/
	/*--------------------------------------i=1~N-1-----------------------------------------------------------------*/
	for(int i=1;i<=N-1;i++)
	{
		for(int j=0;j<=M;j++)
		{
			upie[i][j]=ustar[i][j]-delta_t*((Pres[i+1][j]-Pres[i-1][j])/(2*delta_x));
		}
	}
	/*-----------------------------------------------------------------------------------------------*/
	/*---------------------------------------i=N----------------------------------------------------*/
	for(int i=N;i<N+1;i++)
	{
		for(int j=0;j<1;j++)
		{
			upie[i][j]=uboundary;
		}
	}
	for(int i=N;i<N+1;i++)
	{
		for(int j=1;j<M;j++)
		{
			upie[i][j]=uboundary;
		}
	}
	for(int i=N;i<N+1;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			upie[i][j]=uboundary;
		}
	}
	/*-------------------------------------------------------------------------------------*/
	/*---------------------------------j=0------------------------------------------------*/
	for(int i=0;i<N+1;i++)
	{
		for(int j=0;j<1;j++)
		{
			upie[i][j]=uboundary;
		}
	}
	/*------------------------------------------------------------------------------------*/
	/*---------------------------------j=M----------------------------------------------*/
	for(int i=0;i<N+1;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			upie[i][j]=uboundary;
		}
	}
	/*------------------------------------------------------------------------------------*/
	//cout<<"Output the upie!"<<endl;
	//for(int i=N;i>=0;i--)
	//{
	//	for(int j=0;j<=M;j++)
	//	{
	//		cout<<" "<<upie[i][j];
	//	//	out<<setprecision(5)<<upie[i][j]<<" ";
	//	}
	//cout<<endl;
	////out<<setprecision(5)<<endl;
	//}
	}