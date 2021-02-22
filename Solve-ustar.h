#include "math.h"
using namespace std;
//2013-5-8修改了在节点i=2~N-2 j各点处将部分u[][]写为v[][]
void CUEqution(double **ustar,double **u,double **v,double CF,double delta_x,double delta_y,
			   double delta_t,double Re,double Turb,int N,int M)
{
	/*中间节点*/
	double uplus=0.0,vplus=0.0;
	double usub=0.0,vsub=0.0;
	double uMuTox=0.0,uToTox=0.0;
	double vMutoy=0.0,uToToy=0.0;
	double MuTox=0.0,MuToy=0.0,MvTox=0,S=0.0;
	
	/*----------------------------------------------------------------------------------------------*/
	for(int i=0;i<1;i++)
	{
		for(int j=0;j<1;j++)
		{
			
			/*对流项的离散  */
			uMuTox=u[i][j]/(1*delta_x)*(u[i+1][j]-u[i][j]); //Pian-u/Pian-x采用向后差分 uMuTox 代表 u.Pian-u/Pian-x
			vMutoy=v[i][j]/(1*delta_y)*(u[i][j+1]-u[i][j]);    //Pian-u/Pian-y采用向后差分 vMuToy 代表 v.Pian-u/Pian-y
			/*扩散项的离散 */
			uToTox=(u[i][j]-2*u[i+1][j]+u[i+2][j])/(delta_x*delta_x);
			uToToy=(u[i][j]-2*u[i][j+1]+u[i][j+2])/(delta_y*delta_y);

			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i][j])/(1*delta_x); 
				MuToy=(u[i][j+1]-u[i][j])/(1*delta_y);
				MvTox=(v[i+1][j]-v[i][j])/(1*delta_x); 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=0;i<1;i++)
	{
		for(int j=1;j<M;j++)
		{
			
			/*对流项的离散  */
			uMuTox=u[i][j]/(1*delta_x)*(u[i+1][j]-u[i][j]); //Pian-u/Pian-x采用向后差分 uMuTox 代表 u.Pian-u/Pian-x
			vMutoy=v[i][j]/(2*delta_y)*(u[i][j+1]-u[i][j-1]);    //Pian-u/Pian-y采用中心差分 vMuToy 代表 v.Pian-u/Pian-y
			/*扩散项的离散 */
			uToTox=(u[i][j]-2*u[i+1][j]+u[i+2][j])/(delta_x*delta_x);
			uToToy=(u[i][j+1]-2*u[i][j]+u[i][j-1])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i][j])/(1*delta_x); 
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y);
				MvTox=(v[i+1][j]-v[i][j])/(1*delta_x); 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=0;i<1;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			
			/*对流项的离散  */
			uMuTox=u[i][j]/(1*delta_x)*(u[i+1][j]-u[i][j]); //Pian-u/Pian-x采用向后差分 uMuTox 代表 u.Pian-u/Pian-x
			vMutoy=v[i][j]/(1*delta_y)*(u[i][j]-u[i][j-1]);    //Pian-u/Pian-y采用向前差分 vMuToy 代表 v.Pian-u/Pian-y
			/*扩散项的离散 */
			uToTox=(u[i][j]-2*u[i+1][j]+u[i+2][j])/(delta_x*delta_x);
			uToToy=(u[i][j-2]-2*u[i][j-1]+u[i][j])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i][j])/(1*delta_x); 
				MuToy=(u[i][j]-u[i][j-1])/(1*delta_y);
				MvTox=(v[i+1][j]-v[i][j])/(1*delta_x); 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
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
			uMuTox=u[i][j]/(2*delta_x)*(u[i+1][j]-u[i-1][j]); //Pian-u/Pian-x采用中心差分 uMuTox 代表 u.Pian-u/Pian-x
			vMutoy=v[i][j]/(1*delta_y)*(u[i][j+1]-u[i][j]);    //Pian-u/Pian-y采用向后差分 vMuToy 代表 v.Pian-u/Pian-y
			/*扩散项的离散 */
			uToTox=(u[i+1][j]-2*u[i][j]+u[i-1][j])/(delta_x*delta_x);
			uToToy=(u[i][j]-2*u[i][j+1]+u[i][j+2])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i-1][j])/(2*delta_x); 
				MuToy=(u[i][j+1]-u[i][j])/(1*delta_y);
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x); 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	
	for(int i=1;i<2;i++)
	{
		for(int j=1;j<M;j++)
		{
			
			/*对流项的离散  */
			uMuTox=u[i][j]/(2*delta_x)*(u[i+1][j]-u[i-1][j]); //Pian-u/Pian-x采用向后差分 uMuTox 代表 u.Pian-u/Pian-x
			vMutoy=v[i][j]/(2*delta_y)*(u[i][j+1]-u[i][j-1]);    //Pian-u/Pian-y采用中心差分 vMuToy 代表 v.Pian-u/Pian-y
			/*扩散项的离散 */
			uToTox=(u[i+1][j]-2*u[i][j]+u[i-1][j])/(delta_x*delta_x);
			uToToy=(u[i][j+1]-2*u[i][j]+u[i][j-1])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i-1][j])/(2*delta_x); 
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y);
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x); 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=1;i<2;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			
			/*对流项的离散  */
			uMuTox=u[i][j]/(2*delta_x)*(u[i+1][j]-u[i-1][j]); //Pian-u/Pian-x采用向后差分 uMuTox 代表 u.Pian-u/Pian-x
			vMutoy=v[i][j]/(1*delta_y)*(u[i][j]-u[i][j-1]);    //Pian-u/Pian-y采用向前差分 vMuToy 代表 v.Pian-u/Pian-y
			/*扩散项的离散 */
			uToTox=(u[i+1][j]-2*u[i][j]+u[i-1][j])/(delta_x*delta_x);
			uToToy=(u[i][j-2]-2*u[i][j-1]+u[i][j])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i-1][j])/(2*delta_x); 
				MuToy=(u[i][j]-u[i][j-1])/(1*delta_y);
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x); 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
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
			uMuTox=u[i][j]/(2*delta_x)*(u[i+1][j]-u[i-1][j])
				-uplus*CF/delta_x*(u[i+1][j]-3*u[i][j]+3*u[i-1][j]-u[i-2][j])
				-usub*CF/delta_x*(u[i+2][j]-3*u[i+1][j]+3*u[i][j]-u[i-1][j]); 
			vMutoy=v[i][j]/(1*delta_y)*(u[i][j+1]-u[i][j]); 
			/*扩散项的离散 */
			uToTox=(u[i+1][j]-2*u[i][j]+u[i-1][j])/(delta_x*delta_x);
			uToToy=(u[i][j]-2*u[i][j+1]+u[i][j+2])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i-1][j])/(2*delta_x)
					-CF*(u[i+1][j]-3*u[i][j]+3*u[i-1][j]-u[i-2][j])/delta_x
					-CF*(u[i+2][j]-3*u[i+1][j]+3*u[i][j]-u[i-1][j])/delta_x; 
				MuToy=(u[i][j+1]-u[i][j])/(1*delta_y);
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x)
					-CF*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])/delta_x
					-CF*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j])/delta_x; 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
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
			uMuTox=u[i][j]/(2*delta_x)*(u[i+1][j]-u[i-1][j])
				-uplus*CF/delta_x*(u[i+1][j]-3*u[i][j]+3*u[i-1][j]-u[i-2][j])
				-usub*CF/delta_x*(u[i+2][j]-3*u[i+1][j]+3*u[i][j]-u[i-1][j]);
			vMutoy=v[i][j]/(2*delta_y)*(u[i][j+1]-u[i][j-1]);
			/*扩散项的离散 */
			uToTox=(u[i+1][j]-2*u[i][j]+u[i-1][j])/(delta_x*delta_x);
			uToToy=(u[i][j+1]-2*u[i][j]+u[i][j-1])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i-1][j])/(2*delta_x)
					-CF*(u[i+1][j]-3*u[i][j]+3*u[i-1][j]-u[i-2][j])/delta_x
					-CF*(u[i+2][j]-3*u[i+1][j]+3*u[i][j]-u[i-1][j])/delta_x; 
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y);
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x)
					-CF*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])/delta_x
					-CF*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j])/delta_x; 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=2;i<=N-2;i++)
	{
		for(int j=2;j<=M-2;j++)
		{
			uplus=1.0/2*(u[i][j]+abs(u[i][j]));   //uplus 代表 u+
			usub=1.0/2*(u[i][j]-abs(u[i][j]));     //usub 代表  u-

			vplus=1.0/2*(v[i][j]+abs(v[i][j]));   //vplus 代表 v+
			vsub=1.0/2*(v[i][j]-abs(v[i][j]));     //vsub 代表  v-

			/*对流项的离散  */
			uMuTox=u[i][j]/(2*delta_x)*(u[i+1][j]-u[i-1][j])
				-uplus*CF/delta_x*(u[i+1][j]-3*u[i][j]+3*u[i-1][j]-u[i-2][j])
				-usub*CF/delta_x*(u[i+2][j]-3*u[i+1][j]+3*u[i][j]-u[i-1][j]); 
			vMutoy=v[i][j]/(2*delta_y)*(u[i][j+1]-u[i][j-1])
				-vplus*CF/delta_y*(u[i][j+1]-3*u[i][j]+3*u[i][j-1]-u[i][j-2])
				-vsub*CF/delta_y*(u[i][j+2]-3*u[i][j+1]+3*u[i][j]-u[i][j-1]);
			/*扩散项的离散 */
			uToTox=(u[i+1][j]-2*u[i][j]+u[i-1][j])/(delta_x*delta_x);
			uToToy=(u[i][j+1]-2*u[i][j]+u[i][j-1])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i-1][j])/(2*delta_x)
					-CF*(u[i+1][j]-3*u[i][j]+3*u[i-1][j]-u[i-2][j])/delta_x
					-CF*(u[i+2][j]-3*u[i+1][j]+3*u[i][j]-u[i-1][j])/delta_x; 
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y)
					-CF*(u[i][j+1]-3*u[i][j]+3*u[i][j-1]-u[i][j-2])/delta_y
					-CF*(u[i][j+2]-3*u[i][j+1]+3*u[i][j]-u[i][j-1])/delta_y;
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x)
					-CF*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])/delta_x
					-CF*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j])/delta_x; 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
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
			uMuTox=u[i][j]/(2*delta_x)*(u[i+1][j]-u[i-1][j])
				-uplus*CF/delta_x*(u[i+1][j]-3*u[i][j]+3*u[i-1][j]-u[i-2][j])
				-usub*CF/delta_x*(u[i+2][j]-3*u[i+1][j]+3*u[i][j]-u[i-1][j]); 
			vMutoy=v[i][j]/(2*delta_y)*(u[i][j+1]-u[i][j-1]);
			/*扩散项的离散 */
			uToTox=(u[i+1][j]-2*u[i][j]+u[i-1][j])/(delta_x*delta_x);
			uToToy=(u[i][j+1]-2*u[i][j]+u[i][j-1])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i-1][j])/(2*delta_x)
					-CF*(u[i+1][j]-3*u[i][j]+3*u[i-1][j]-u[i-2][j])/delta_x
					-CF*(u[i+2][j]-3*u[i+1][j]+3*u[i][j]-u[i-1][j])/delta_x; 
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y);
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x)
					-CF*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])/delta_x
					-CF*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j])/delta_x; 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
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
			uMuTox=u[i][j]/(2*delta_x)*(u[i+1][j]-u[i-1][j])
				-uplus*CF/delta_x*(u[i+1][j]-3*u[i][j]+3*u[i-1][j]-u[i-2][j])
				-usub*CF/delta_x*(u[i+2][j]-3*u[i+1][j]+3*u[i][j]-u[i-1][j]); 
			vMutoy=v[i][j]/(1*delta_y)*(u[i][j]-u[i][j-1]);
			/*扩散项的离散 */
			uToTox=(u[i+1][j]-2*u[i][j]+u[i-1][j])/(delta_x*delta_x);
			uToToy=(u[i][j]-2*u[i][j-1]+u[i][j-2])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i-1][j])/(2*delta_x)
					-CF*(u[i+1][j]-3*u[i][j]+3*u[i-1][j]-u[i-2][j])/delta_x
					-CF*(u[i+2][j]-3*u[i+1][j]+3*u[i][j]-u[i-1][j])/delta_x; 
				MuToy=(u[i][j]-u[i][j-1])/(1*delta_y);
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x)
					-CF*(v[i+1][j]-3*v[i][j]+3*v[i-1][j]-v[i-2][j])/delta_x
					-CF*(v[i+2][j]-3*v[i+1][j]+3*v[i][j]-v[i-1][j])/delta_x; 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
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
			uMuTox=u[i][j]/(2*delta_x)*(u[i+1][j]-u[i-1][j]); 
			vMutoy=v[i][j]/(1*delta_y)*(u[i][j+1]-u[i][j]);
			/*扩散项的离散 */
			uToTox=(u[i+1][j]-2*u[i][j]+u[i-1][j])/(delta_x*delta_x);
			uToToy=(u[i][j]-2*u[i][j+1]+u[i][j+2])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i-1][j])/(2*delta_x); 
				MuToy=(u[i][j+1]-u[i][j])/(1*delta_y);
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x); 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}

	for(int i=N-1;i<N;i++)
	{
		for(int j=1;j<M;j++)
		{
			
			/*对流项的离散  */
			uMuTox=u[i][j]/(2*delta_x)*(u[i+1][j]-u[i-1][j]); //Pian-u/Pian-x采用向后差分 uMuTox 代表 u.Pian-u/Pian-x
			vMutoy=v[i][j]/(2*delta_y)*(u[i][j+1]-u[i][j-1]);    //Pian-u/Pian-y采用中心差分 vMuToy 代表 v.Pian-u/Pian-y
			/*扩散项的离散 */
			uToTox=(u[i+1][j]-2*u[i][j]+u[i-1][j])/(delta_x*delta_x);
			uToToy=(u[i][j+1]-2*u[i][j]+u[i][j-1])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i-1][j])/(2*delta_x); 
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y);
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x); 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=N-1;i<N;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			
			/*对流项的离散  */
			uMuTox=u[i][j]/(2*delta_x)*(u[i+1][j]-u[i-1][j]); //Pian-u/Pian-x采用向后差分 uMuTox 代表 u.Pian-u/Pian-x
			vMutoy=v[i][j]/(1*delta_y)*(u[i][j]-u[i][j-1]);    //Pian-u/Pian-y采用向前差分 vMuToy 代表 v.Pian-u/Pian-y
			/*扩散项的离散 */
			uToTox=(u[i+1][j]-2*u[i][j]+u[i-1][j])/(delta_x*delta_x);
			uToToy=(u[i][j-2]-2*u[i][j-1]+u[i][j])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i+1][j]-u[i-1][j])/(2*delta_x); 
				MuToy=(u[i][j]-u[i][j-1])/(1*delta_y);
				MvTox=(v[i+1][j]-v[i-1][j])/(2*delta_x); 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
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
			uMuTox=u[i][j]/(1*delta_x)*(u[i][j]-u[i-1][j]); 
			vMutoy=v[i][j]/(1*delta_y)*(u[i][j+1]-u[i][j]); 
			/*扩散项的离散 */
			uToTox=(u[i][j]-2*u[i-1][j]+u[i-2][j])/(delta_x*delta_x);
			uToToy=(u[i][j]-2*u[i][j+1]+u[i][j+2])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i][j]-u[i-1][j])/(1*delta_x); 
				MuToy=(u[i][j+1]-u[i][j])/(1*delta_y);
				MvTox=(v[i][j]-v[i-1][j])/(1*delta_x); 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=N;i<N+1;i++)
	{
		for(int j=1;j<M;j++)
		{
			
			/*对流项的离散  */
			uMuTox=u[i][j]/(1*delta_x)*(u[i][j]-u[i-1][j]); 
			vMutoy=v[i][j]/(2*delta_y)*(u[i][j+1]-u[i][j-1]); 
			/*扩散项的离散 */
			uToTox=(u[i][j]-2*u[i-1][j]+u[i-2][j])/(delta_x*delta_x);
			uToToy=(u[i][j+1]-2*u[i][j]+u[i][j-1])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i][j]-u[i-1][j])/(1*delta_x); 
				MuToy=(u[i][j+1]-u[i][j-1])/(2*delta_y);
				MvTox=(v[i][j]-v[i-1][j])/(1*delta_x); 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
	for(int i=N;i<N+1;i++)
	{
		for(int j=M;j<M+1;j++)
		{
			
			/*对流项的离散  */
			uMuTox=u[i][j]/(1*delta_x)*(u[i][j]-u[i-1][j]); 
			vMutoy=v[i][j]/(1*delta_y)*(u[i][j]-u[i][j-1]);  
			/*扩散项的离散 */
			uToTox=(u[i][j]-2*u[i-1][j]+u[i-2][j])/(delta_x*delta_x);
			uToToy=(u[i][j-2]-2*u[i][j-1]+u[i][j])/(delta_y*delta_y);
			
			if(Turb==1.0)   //考虑湍流效应
			{
				/*亚格子应力项的离散 */
				MuTox=(u[i][j]-u[i-1][j])/(1*delta_x); 
				MuToy=(u[i][j]-u[i][j-1])/(1*delta_y);
				MvTox=(v[i][j]-v[i-1][j])/(1*delta_x); 
				S=MuTox+MuToy+MuTox+MvTox;
			}
			
			ustar[i][j]=u[i][j]-delta_t*(uMuTox+vMutoy-1/Re*(uToTox+uToToy)
				-sqrt(2.0)/200.0*pow(delta_x,2.0)*fabs(S)*S);
		}
	}
}