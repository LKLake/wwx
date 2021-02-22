   //1RINV.CPP
  //实矩阵求逆
  #include  <iostream>
  #include  <cmath>
  #include  <fstream>
  using namespace std;
  void rinv(double  **a,int n)        //执行实矩阵求逆
  { 
	  int *is,*js,i,j,k;
      double d,p;
      is = new int[n];
      js = new int[n];
      for (k=0; k<=n-1; k++)
      { 
		  d=0.0;
          for (i=k; i<=n-1; i++)
          for (j=k; j<=n-1; j++)
          { 
			  p=fabs(a[i][j]);
              if (p>d) { d=p; is[k]=i; js[k]=j;}
          }
          if (d+1.0==1.0)
          { 
			  delete [] is,js;
              cout <<"\nA为奇异矩阵！没有逆矩阵. " <<endl;
			  exit(1);
          }
          if (is[k]!=k)
            for (j=0; j<=n-1; j++)
            { 
                p=a[k][j]; a[k][j]=a[is[k]][j]; a[is[k]][j]=p;
            }
          if (js[k]!=k)
            for (i=0; i<=n-1; i++)
            { 
                p=a[i][k]; a[i][k]=a[i][js[k]]; a[i][js[k]]=p;
            }
          a[k][k]=1.0/a[k][k];
          for (j=0; j<=n-1; j++)
            if (j!=k)  a[k][j]=a[k][j]*a[k][k];
          for (i=0; i<=n-1; i++)
            if (i!=k)
              for (j=0; j<=n-1; j++)
                if (j!=k) a[i][j]=a[i][j]-a[i][k]*a[k][j];
          for (i=0; i<=n-1; i++)
            if (i!=k)  a[i][k]=-a[i][k]*a[k][k];
      }
      for (k=n-1; k>=0; k--)
      { 
		  if (js[k]!=k)
            for (j=0; j<=n-1; j++)
            { 
                p=a[k][j]; a[k][j]=a[js[k]][j]; a[js[k]][j]=p;
            }
          if (is[k]!=k)
            for (i=0; i<=n-1; i++)
            { 
                p=a[i][k]; a[i][k]=a[i][is[k]]; a[i][is[k]]=p;
            }
      }
      delete [] is, js;
  }

  