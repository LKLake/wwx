  //3GJDN.CPP
  //全选主元Gauss-Jordan消去法求解实系数方程组
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  
  void gjdn(double  **a, double **b, int n, int m)        //执行Gauss-Jordan消去法
  { 
	 	  
	  int *js,l,k,i,j,is;
      double d,t;
      js = new int[n];
      l=1;
      for (k=0;k<=n-1;k++)
      { 
		  d=0.0;
          for (i=k;i<=n-1;i++)
          for (j=k;j<=n-1;j++)
          { 
			  t=fabs(a[i][j]);
              if (t>d) { d=t; js[k]=j; is=i;}
          }
          if (d+1.0==1.0) l=0;
          else
          { 
			  if (js[k]!=k)
                for (i=0;i<=n-1;i++)
                { 
                    t=a[i][k]; 
					a[i][k]=a[i][js[k]]; 
					a[i][js[k]]=t;
                }
              if (is!=k)
              { 
				  for (j=k;j<=n-1;j++)
                  { 
                      t=a[k][j]; 
					  a[k][j]=a[is][j]; 
					  a[is][j]=t;
                  }
                  for (j=0;j<=m-1;j++)
                  { 
                      t=b[k][j]; 
					  b[k][j]=b[is][j]; 
					  b[is][j]=t;
                  }
              }
          }
          if (l==0)
          { 
			  delete []js;
			  cout <<"\n系数矩阵奇异！无解. " <<endl;
              return;
          }
          d=a[k][k];
          for (j=k+1;j<=n-1;j++)
             a[k][j]=a[k][j]/d;
          for (j=0;j<=m-1;j++)
             b[k][j]=b[k][j]/d;
          for (j=k+1;j<=n-1;j++)
          for (i=0;i<=n-1;i++)
          { 
              if (i!=k)
                a[i][j]=a[i][j]-a[i][k]*a[k][j];
          }
          for (j=0;j<=m-1;j++)
          for (i=0;i<=n-1;i++)
          { 
              if (i!=k)
              b[i][j]=b[i][j]-a[i][k]*b[k][j];
          }
      }
      for (k=n-1;k>=0;k--)
        if (js[k]!=k)
          for (j=0;j<=m-1;j++)
          { 
              t=b[k][j]; b[k][j]=b[js[k]][j]; b[js[k]][j]=t;
          }
      delete [] js;
	  
  }

  