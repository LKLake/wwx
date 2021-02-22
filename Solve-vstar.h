using namespace std;

void CVEqution(double **vstar,double **u,double **v,double CF,double delta_x,double delta_y,
			   double delta_t,double Re,double Turb,int N,int M)
{
	/*中间节点*/
	double uplus=0.0,vplus=0.0;
	double usub=0.0,vsub=0.0;
	double uMvTox=0.0,vToTox=0.0;
	double vMvToy=0.0,vToToy=0.0;
	double MvTox=0.0,MvToy=0.0,MuToy=0.0,S=0.0;
	/*----------------------------------------------------------------------------------------------*/
	for(int i=0;i<1;i++)
	{
		for(int j=0;j<1;j++)
		{
			/*对流项的离散  */
			uMvTox=u[i][j]/(1*delta_x)*(v[i+1][j]-v[i][j]); 
			vMvToy=v[i][j]/(1*delta_y)*(v[i][j+1]-v[i][j]);
			/*扩散项的离散 */
			vToTox=(v[i][j]-2*v[i+1][j]+v[i+2][j])/(delta_x*delta_x);
			vToToy=(v[i][j]-2*v[i][j+1]+v[i][j+2])/(delta_y*delta_y);		
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i][j])/(1*delta_x); 
				MvToy=(v[i][j+1]-v[i][j])/(1*delta_y);
				MuToy=(u[i][j+1]-u[i][j])/(1*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=0;i<1;i++)
	{
		for(int j=1;j<M;j++)
		{
			
			/*对流项的离散  */
			uMvTox=u[i][j]/(1*delta_x)*(v[i+1][j]-v[i][j]);  
			vMvToy=v[i][j]/(2*delta_y)*(v[i][j+1]-v[i][j-1]); 
			/*扩散项的离散 */
			vToTox=(v[i][j]-2*v[i+1][j]+v[i+2][j])/(delta_x*delta_x);
			vToToy=(v[i][j+1]-2*v[i][j]+v[i][j-1])/(delta_y*delta_y);			
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i][j])/(1*delta_x); 
				MvToy=(v[i][j+1]-v[i][j-1])/(2*delta_y);
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}

	for(int i=0;i<1;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			
			/*对流项的离散  */
			uMvTox=u[i][j]/(1*delta_x)*(v[i+1][j]-v[i][j]);
			vMvToy=v[i][j]/(1*delta_y)*(v[i][j]-v[i][j-1]);
			/*扩散项的离散 */
			vToTox=(v[i][j]-2*v[i+1][j]+v[i+2][j])/(delta_x*delta_x);
			vToToy=(v[i][j-2]-2*v[i][j-1]+v[i][j])/(delta_y*delta_y);			
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i][j])/(1*delta_x); 
				MvToy=(v[i][j]-v[i][j-1])/(1*delta_y);
				MuToy=(u[i][j]-u[i][j-1])/(1*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	/*-------------------------------------------------------------------------------------------*/
	/*-------------------------------------------------------------------------------------------*/
	for(int i=1;i<2;i++)
	{
		for(int j=0;j<1;j++)
		{
			
			/*对流项的离散  */
			uMvTox=u[i][j]/(2*delta_x)*(v[i+1][j]-v[i-1][j]); 
			vMvToy=v[i][j]/(1*delta_y)*(v[i][j+1]-v[i][j]);
			/*扩散项的离散 */
			vToTox=(v[i+1][j]-2*v[i][j]+v[i-1][j])/(delta_x*delta_x);
			vToToy=(v[i][j]-2*v[i][j+1]+v[i][j+2])/(delta_y*delta_y);			
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x); 
				MvToy=(v[i][j+1]-v[i][j])/(1*delta_y);
				MuToy=(u[i][j+1]-u[i][j])/(1*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=1;i<2;i++)
	{
		for(int j=1;j<M;j++)
		{
			
			/*对流项的离散  */
			uMvTox=u[i][j]/(2*delta_x)*(v[i+1][j]-v[i-1][j]);  
			vMvToy=v[i][j]/(2*delta_y)*(v[i][j+1]-v[i][j-1]); 
			/*扩散项的离散 */
			vToTox=(v[i+1][j]-2*v[i][j]+v[i-1][j])/(delta_x*delta_x);
			vToToy=(v[i][j+1]-2*v[i][j]+v[i][j-1])/(delta_y*delta_y);			
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x); 
				MvToy=(v[i][j+1]-v[i][j-1])/(2*delta_y);
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=1;i<2;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			
			/*对流项的离散  */
			uMvTox=u[i][j]/(2*delta_x)*(v[i+1][j]-v[i-1][j]); 
			vMvToy=v[i][j]/(1*delta_y)*(v[i][j]-v[i][j-1]);
			/*扩散项的离散 */
			vToTox=(v[i+1][j]-2*v[i][j]+v[i-1][j])/(delta_x*delta_x);
			vToToy=(v[i][j-2]-2*v[i][j-1]+v[i][j])/(delta_y*delta_y);			
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x); 
				MvToy=(v[i][j]-v[i][j-1])/(1*delta_y);
				MuToy=(u[i][j]-u[i][j-1])/(1*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	/*-----------------------------------------------------------------------------------------------------------------*/
	/*-----------------------------------------------------------------------------------------------------------------*/
	for(int i=2;i<=N-2;i++)
	{
		for(int j=0;j<1;j++)
		{
			uplus=1.0/2*(u[i][j]+abs(u[i][j]));   //uplus 代表 u+
			usub=1.0/2*(u[i][j]-abs(u[i][j]));     //usub 代表  u-

			vplus=1.0/2*(v[i][j]+abs(v[i][j]));   //vplus 代表 v+
			vsub=1.0/2*(v[i][j]-abs(v[i][j]));     //vsub 代表  v-

			/*对流项的离散  */
			uMvTox=u[i][j]/(2*delta_x)*(v[i+1][j]-v[i-1][j])
				-uplus*CF/delta_x*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])
				-usub*CF/delta_x*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j]);
			vMvToy=v[i][j]/(1*delta_y)*(v[i][j+1]-v[i][j]);
			/*扩散项的离散 */
			vToTox=(v[i+1][j]-2*v[i][j]+v[i-1][j])/(delta_x*delta_x);
			vToToy=(v[i][j]-2*v[i][j+1]+v[i][j+2])/(delta_y*delta_y);			
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x)
					-CF*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])/delta_x
					-CF*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j])/delta_x;
				MvToy=(v[i][j+1]-v[i][j])/(1*delta_y);
				MuToy=(u[i][j+1]-u[i][j])/(1*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=2;i<=N-2;i++)
	{
		for(int j=1;j<2;j++)
		{
			uplus=1.0/2*(u[i][j]+abs(u[i][j]));   //uplus 代表 u+
			usub=1.0/2*(u[i][j]-abs(u[i][j]));     //usub 代表  u-

			vplus=1.0/2*(v[i][j]+abs(v[i][j]));   //vplus 代表 v+
			vsub=1.0/2*(v[i][j]-abs(v[i][j]));     //vsub 代表  v-

			/*对流项的离散  */
			uMvTox=u[i][j]/(2*delta_x)*(v[i+1][j]-v[i-1][j])
				-uplus*CF/delta_x*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])
				-usub*CF/delta_x*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j]); 
			vMvToy=v[i][j]/(2*delta_y)*(v[i][j+1]-v[i][j-1]);
			/*扩散项的离散 */
			vToTox=(v[i+1][j]-2*v[i][j]+v[i-1][j])/(delta_x*delta_x);
			vToToy=(v[i][j+1]-2*v[i][j]+v[i][j-1])/(delta_y*delta_y);			
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x)
					-CF*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])/delta_x
					-CF*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j])/delta_x;
				MvToy=(v[i][j+1]-v[i][j-1])/(2*delta_y);
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=2;i<=N-2;i++)
	{
		for(int j=2;j<M-1;j++)
		{
			uplus=1.0/2*(u[i][j]+abs(u[i][j]));   //uplus 代表 u+
			usub=1.0/2*(u[i][j]-abs(u[i][j]));     //usub 代表  u-

			vplus=1.0/2*(v[i][j]+abs(v[i][j]));   //vplus 代表 v+
			vsub=1.0/2*(v[i][j]-abs(v[i][j]));     //vsub 代表  v-

			/*对流项的离散  */
			uMvTox=u[i][j]/(2*delta_x)*(v[i+1][j]-v[i-1][j])
				-uplus*CF/delta_x*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])
				-usub*CF/delta_x*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j]); 
			vMvToy=v[i][j]/(2*delta_y)*(v[i][j+1]-v[i][j-1])
				-vplus*CF/delta_y*(v[i][j+1]-3*v[i][j]+3*v[i][j-1]-v[i][j-2])
				-vsub*CF/delta_y*(v[i][j+2]-3*v[i][j+1]+3*v[i][j]-v[i][j-1]);
			/*扩散项的离散 */
			vToTox=(v[i+1][j]-2*v[i][j]+v[i-1][j])/(delta_x*delta_x);
			vToToy=(v[i][j+1]-2*v[i][j]+v[i][j-1])/(delta_y*delta_y);			
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x)
					-CF*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])/delta_x
					-CF*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j])/delta_x;
				MvToy=(v[i][j+1]-v[i][j-1])/(2*delta_y)
					-CF*(v[i][j+1]-3*v[i][j]+3*v[i][j-1]-v[i][j-2])/delta_y
					-CF*(v[i][j+2]-3*v[i][j+1]+3*v[i][j]-v[i][j-1])/delta_y;
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y)
					-CF*(u[i][j+1]-3*u[i][j]+3*u[i][j-1]-u[i][j-2])/delta_y
					-CF*(u[i][j+2]-3*u[i][j+1]+3*u[i][j]-u[i][j-1])/delta_y;
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=2;i<=N-2;i++)
	{
		for(int j=M-1;j<M;j++)
		{
			uplus=1.0/2*(u[i][j]+abs(u[i][j]));   //uplus 代表 u+
			usub=1.0/2*(u[i][j]-abs(u[i][j]));     //usub 代表  u-

			vplus=1.0/2*(v[i][j]+abs(v[i][j]));   //vplus 代表 v+
			vsub=1.0/2*(v[i][j]-abs(v[i][j]));     //vsub 代表  v-

			/*对流项的离散  */
			uMvTox=u[i][j]/(2*delta_x)*(v[i+1][j]-v[i-1][j])
				-uplus*CF/delta_x*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])
				-usub*CF/delta_x*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j]); 
			vMvToy=v[i][j]/(2*delta_y)*(v[i][j+1]-v[i][j-1]);
			/*扩散项的离散 */
			vToTox=(v[i+1][j]-2*v[i][j]+v[i-1][j])/(delta_x*delta_x);
			vToToy=(v[i][j+1]-2*v[i][j]+v[i][j-1])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x)
					-CF*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])/delta_x
					-CF*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j])/delta_x;
				MvToy=(v[i][j+1]-v[i][j-1])/(2*delta_y);
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=2;i<=N-2;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			uplus=1.0/2*(u[i][j]+abs(u[i][j]));   //uplus 代表 u+
			usub=1.0/2*(u[i][j]-abs(u[i][j]));     //usub 代表  u-

			vplus=1.0/2*(v[i][j]+abs(v[i][j]));   //vplus 代表 v+
			vsub=1.0/2*(v[i][j]-abs(v[i][j]));     //vsub 代表  v-

			/*对流项的离散  */
			uMvTox=u[i][j]/(2*delta_x)*(v[i+1][j]-v[i-1][j])
				-uplus*CF/delta_x*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])
				-usub*CF/delta_x*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j]); 
			vMvToy=v[i][j]/(1*delta_y)*(v[i][j]-v[i][j-1]);
			/*扩散项的离散 */
			vToTox=(v[i+1][j]-2*v[i][j]+v[i-1][j])/(delta_x*delta_x);
			vToToy=(v[i][j]-2*v[i][j-1]+v[i][j-2])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x)
					-CF*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])/delta_x
					-CF*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j])/delta_x;
				MvToy=(v[i][j]-v[i][j-1])/(1*delta_y);
				MuToy=(u[i][j]-u[i][j-1])/(1*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	/*---------------------------------------------------------------------------------------------*/
	/*---------------------------------------------------------------------------------------------*/
	for(int i=N-1;i<N;i++)
	{
		for(int j=0;j<1;j++)
		{
			/*对流项的离散  */
			uMvTox=u[i][j]/(2*delta_x)*(v[i+1][j]-v[i-1][j]); 
			vMvToy=v[i][j]/(1*delta_y)*(v[i][j+1]-v[i][j]); 
			/*扩散项的离散 */
			vToTox=(v[i+1][j]-2*v[i][j]+v[i-1][j])/(delta_x*delta_x);
			vToToy=(v[i][j]-2*v[i][j+1]+v[i][j+2])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x); 
				MvToy=(v[i][j+1]-v[i][j])/(1*delta_y);
				MuToy=(u[i][j+1]-u[i][j])/(1*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=N-1;i<N;i++)
	{
		for(int j=1;j<M;j++)
		{
			/*对流项的离散  */
			uMvTox=u[i][j]/(2*delta_x)*(v[i+1][j]-v[i-1][j]);
			vMvToy=v[i][j]/(2*delta_y)*(v[i][j+1]-v[i][j-1]);
			/*扩散项的离散 */
			vToTox=(v[i+1][j]-2*v[i][j]+v[i-1][j])/(delta_x*delta_x);
			vToToy=(v[i][j+1]-2*v[i][j]+v[i][j-1])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x); 
				MvToy=(v[i][j+1]-v[i][j-1])/(2*delta_y);
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=N-1;i<N;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			/*对流项的离散  */
			uMvTox=u[i][j]/(2*delta_x)*(v[i+1][j]-v[i-1][j]);
			vMvToy=v[i][j]/(1*delta_y)*(v[i][j]-v[i][j-1]);
			/*扩散项的离散 */
			vToTox=(v[i+1][j]-2*v[i][j]+v[i-1][j])/(delta_x*delta_x);
			vToToy=(v[i][j-2]-2*v[i][j-1]+v[i][j])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x); 
				MvToy=(v[i][j]-v[i][j-1])/(1*delta_y);
				MuToy=(u[i][j]-u[i][j-1])/(1*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	/*-----------------------------------------------------------------------------------------------*/
	/*-----------------------------------------------------------------------------------------------*/
	for(int i=N;i<N+1;i++)
	{
		for(int j=0;j<1;j++)
		{
			/*对流项的离散  */
			uMvTox=u[i][j]/(1*delta_x)*(v[i][j]-v[i-1][j]); 
			vMvToy=v[i][j]/(1*delta_y)*(v[i][j+1]-v[i][j]);
			/*扩散项的离散 */
			vToTox=(v[i][j]-2*v[i-1][j]+v[i-2][j])/(delta_x*delta_x);
			vToToy=(v[i][j]-2*v[i][j+1]+v[i][j+2])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i][j]-v[i-1][j])/(1*delta_x); 
				MvToy=(v[i][j+1]-v[i][j])/(1*delta_y);
				MuToy=(u[i][j+1]-u[i][j])/(1*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=N;i<N+1;i++)
	{
		for(int j=1;j<M;j++)
		{
			/*对流项的离散  */
			uMvTox=u[i][j]/(1*delta_x)*(v[i][j]-v[i-1][j]); 
			vMvToy=v[i][j]/(2*delta_y)*(v[i][j+1]-v[i][j-1]);
			/*扩散项的离散 */
			vToTox=(v[i][j]-2*v[i-1][j]+v[i-2][j])/(delta_x*delta_x);
			vToToy=(v[i][j+1]-2*v[i][j]+v[i][j-1])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i][j]-v[i-1][j])/(1*delta_x); 
				MvToy=(v[i][j+1]-v[i][j-1])/(2*delta_y);
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=N;i<N+1;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			/*对流项的离散  */
			uMvTox=u[i][j]/(1*delta_x)*(v[i][j]-v[i-1][j]);
			vMvToy=v[i][j]/(1*delta_y)*(v[i][j]-v[i][j-1]);
			/*扩散项的离散 */
			vToTox=(v[i][j]-2*v[i-1][j]+v[i-2][j])/(delta_x*delta_x);
			vToToy=(v[i][j-2]-2*v[i][j-1]+v[i][j])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MvTox=(v[i][j]-v[i-1][j])/(1*delta_x);
				MvToy=(v[i][j]-v[i][j-1])/(1*delta_y);
				MuToy=(u[i][j]-u[i][j-1])/(1*delta_y); 
				S=MvTox+MvToy+MuToy+MvToy;
			}
			
			vstar[i][j]=v[i][j]-delta_t*(uMvTox+vMvToy-1/Re*(vToTox+vToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
}