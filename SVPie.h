using namespace std;

void SolVPie(double **vstar,double **vpie,double **Pres,int N,int M )
{
	double Length=15.0;   //求解区域的长
	double Wide=10.0;      //求解区域的宽
	double delta_x=Length/N;
	double delta_y=Wide/M;
	double delta_t=0.001;
	double vboundary=0.0;
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
			vpie[i][j]=vboundary;
		}
	}
	for(int i=0;i<1;i++)
	{
		for(int j=1;j<M;j++)
		{
			vpie[i][j]=vboundary;
		}
	}
	for(int i=0;i<1;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			vpie[i][j]=vboundary;
		}
	}
	/*-----------------------------------------------------------------------------------------------------------------*/
	/*--------------------------------------i=1~N-1-----------------------------------------------------------------*/
	for(int i=1;i<=N-1;i++)
	{
		for(int j=1;j<M;j++)
		{
			vpie[i][j]=vstar[i][j]-delta_t*((Pres[i][j+1]-Pres[i][j-1])/(2*delta_x));
		}
	}
	/*-----------------------------------------------------------------------------------------------*/
	/*---------------------------------------i=N----------------------------------------------------*/
	for(int i=N;i<N+1;i++)
	{
		for(int j=0;j<1;j++)
		{
			vpie[i][j]=vboundary;
		}
	}
	for(int i=N;i<N+1;i++)
	{
		for(int j=1;j<M;j++)
		{
			vpie[i][j]=vboundary;
		}
	}
	for(int i=N;i<N+1;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			vpie[i][j]=vboundary;
		}
	}
	/*-------------------------------------------------------------------------------------*/
	/*---------------------------------j=0------------------------------------------------*/
	for(int i=0;i<N+1;i++)
	{
		for(int j=0;j<1;j++)
		{
			vpie[i][j]=vboundary;
		}
	}
	/*------------------------------------------------------------------------------------*/
	/*---------------------------------j=M----------------------------------------------*/
	for(int i=0;i<N+1;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			vpie[i][j]=vboundary;
		}
	}
	/*------------------------------------------------------------------------------------*/
	//cout<<"Output the vpie!"<<endl;
	//for(int i=N;i>=0;i--)
	//{
	//	for(int j=0;j<=M;j++)
	//	{
	//		cout<<" "<<vpie[i][j];
	//		//	out<<setprecision(5)<<vpie[i][j]<<" ";
	//	}
	//	cout<<endl;
	//	//out<<setprecision(5)<<endl;
	//}
}