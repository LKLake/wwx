
using namespace std;
#define A(i,j)  A[(i)*N+j]
void LU(double *A, int N,double b[],double x[])
{

//     double U[N][N]={0};
         double **U=new double *[N];
         for(int k=0;k<N;k++)
                U[k]= new double [N];
         for(int i=0;i<N;i++)
                 for(int j=0;j<N;j++)
                        U[i][j]=0;
//     double L[N][N]={0};
         double **L=new double *[N];
         for(int k=0;k<N;k++)
                L[k]= new double [N];
         for(int i=0;i<N;i++)
                 for(int j=0;j<N;j++)
                        L[i][j]=0;
         for(int i=0;i<N;i++)
        {
                U[0][i]=A(0,i);
                L[i][0]=A(i,0)/U[0][0];
        }
         for(int r=1;r<N;r++)
        {
                 for(int i=r;i<N;i++)
                {
                         double s=0;
                         for(int k=0;k<=r-1;k++)
                                s=s+L[r][k]*U[k][i];
                        U[r][i]=A(r,i)-s;
                         double t=0;
                         for(int k=0;k<=r-1;k++)
                                t=t+L[i][k]*U[k][r];
                        L[i][r]=(A(i,r)-t)/U[r][r];
                }
        }

//     double y[N];
         double *y =new double [N];
        y[0]=b[0];

         for(int   i=1;i<N;i++)
        {
                 double ss=0;
                 for(int k=0;k<=i-1;k++)
                        ss=ss+L[i][k]*y[k];
                y[i]=b[i]-ss;
//             printf("%lf ",y[i]);
        }
        printf( "\n");
//     double x[N];
        cout<< "Output the current X" <<endl;
        x[N-1]=y[N-1]/U[N-1][N-1];
        printf( "%lf  " ,x[N-1]);
         for(int i=N-2;i>=0;i--)
        {
                 double ss=0;
                 for(int k=i+1;k<N;k++)
                        ss=ss+U[i][k]*x[k];
                x[i]=(y[i]-ss)/U[i][i];
                printf( "%lf " ,x[i]);
        }
        cout<<endl;
         for(int t=0;t<N;t++)
                 delete []U[t];
         delete []U;
         for(int t=0;t<N;t++)
                 delete []L[t];
         delete []L;
         delete [] y;
}
