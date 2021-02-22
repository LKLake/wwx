#include <fstream>

#include"SBoundAndInt.h"

void SolEta(double **eta,int N,int M)
{
	double Length=30.0;   //求解区域的长
	double Wide=15.0;      //求解区域的宽
	double center_x=7.5;
	double center_y=7.5;
	double Radius=1;
	double SL=center_x-Radius;   //固体的最左边
	double SR=center_x+Radius;   //固体的最右边
	double SU=center_y+Radius;   //固体的最上边
	double SD=center_y-Radius;   //固体的最下边
	double RD=center_x-Radius;    //RD 区间的下限
	double RU=center_x+Radius;    //RU  区间的上限
	double *x=new double [N+1];
	double *y=new double [M+1];

	double delta_x=Length/N;
	double delta_y=Wide/M;
	for(int i=0;i<N+1;i++)
		x[i]=i*delta_x;
	for(int j=0;j<M+1;j++)
		y[j]=j*delta_y;

	/*-----------------------------------------------------------------------------------------------*/
	for(int i=1;i<N+1;i++)    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	{
		for(int j=1;j<M+1;j++)
		{
			if(x[i]<=SL) {eta[i-1][j-1]=0;continue;}
			if(x[i-1]>=SR){eta[i-1][j-1]=0;continue;}
			if(y[j]<=SD){eta[i-1][j-1]=0;continue;}
			if(y[j-1]>=SU){eta[i-1][j-1]=0;continue;}
			if(x[i-1]<=SL&&SL<x[i])
			{
				//cout<<"***"<<endl;
				//out<<setprecision(5)<<"***"<<" ";
				double f1,f2;//f1 is bigger than f2
				double XFD,XFU;
				double xp=x[i];
				if(Fun_Solve_y(SD,(SD+SU)/2,xp,center_x,center_y,Radius)<Fun_Solve_y((SD+SU)/2,SU,xp,center_x,center_y,Radius)) 
				{
					f1=Fun_Solve_y((SD+SU)/2,SU,xp,center_x,center_y,Radius);
					f2=Fun_Solve_y(SD,(SD+SU)/2,xp,center_x,center_y,Radius);
				}
				else {
					f1=Fun_Solve_y(SD,(SD+SU)/2,xp,center_x,center_y,Radius);
					f2=Fun_Solve_y((SD+SU)/2,SU,xp,center_x,center_y,Radius);
				}
				//cout<<"f1="<<f1<<" "<<"f2="<<f2<<endl;
				//out<<setprecision(5)<<"f1="<<f1<<" "<<"f2="<<f2<<endl;
				if(y[j]<f2) {eta[i-1][j-1]=0;continue;}
				if(y[j-1]>f1)     {eta[i-1][j-1]=0;continue;}
				if(y[j-1]>f2&&y[j]<f1) 
				{ 
					//	cout<<"RD="<<RD<<"(RD+RU)/2="<<(RD+RU)/2<<"y="<<y[j];
					XFD=Fun_Solve_x(RD,(RD+RU)/2,center_x,y[j-1],center_y,Radius);
					XFU=Fun_Solve_x(RD,(RD+RU)/2,center_x,y[j],center_y,Radius);
					//cout<<"XFD="<<XFD<<" "<<"XFU="<<XFU<<endl;
					//out<<setprecision(5)<<"XFD="<<XFD<<" "<<"XFU="<<XFU<<endl;
					double Temp_area,Area;
					if(y[j-1]>=center_y)  { Temp_area=RombergU(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
					Area=Temp_area+(x[i]-XFU)*(y[j]-y[j-1]);
					//cout<<"Temp_area="<<Temp_area<<"x[i]="<<x[i]<<XFU<<endl;
					eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
					if(y[j]<=center_y)  { Temp_area=RombergD(min(XFU,XFD),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
					Area=-Temp_area+(x[i]-XFU)*(y[j]-y[j-1]);
					//  cout<<"Temp_area="<<Temp_area<<"x[i]="<<x[i]<<XFU<<endl;
					eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
					if(y[j-1]<center_y&&y[j]>center_y){ double Temp1,Temp2;
					Temp1=RombergU(RD,XFU,center_x,center_y,Radius)-abs(XFU-RD)*center_y+(x[i]-XFU)*(y[j]-center_y);
					Temp2=-(RombergD(RD,XFD,center_x,center_y,Radius)-abs(XFD-RD)*y[j-1])+(x[i]-RD)*abs(y[j-1]-center_y);
					Area=Temp1+Temp2;
					eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
				}
				if(y[j-1]<f1&&y[j]>=f1) { double Area;
				XFD=Fun_Solve_x(RD,(RD+RU)/2,center_x,y[j-1],center_y,Radius);
				Area=RombergU(XFD,x[i],center_x,center_y,Radius)-abs(XFD-x[i])*y[j-1];
				eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
				if(y[j-1]<=f2&&y[j]>f2) { double Temp_area,Area;
				XFU=Fun_Solve_x(RD,(RD+RU)/2,center_x,y[j],center_y,Radius);
				Temp_area=RombergD(XFU,x[i],center_x,center_y,Radius)-abs(XFU-x[i])*y[j-1];
				Area=-Temp_area+abs(XFU-x[i])*(y[j]-y[j-1]);
				eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			}
			if(x[i-1]<SR&&SR<=x[i])
			{
				//cout<<"***"<<endl;
				double f1,f2;//f1 is bigger than f2
				double XFD,XFU;
				double xp=x[i-1];
				if(Fun_Solve_y(SD,(SD+SU)/2,xp,center_x,center_y,Radius)<Fun_Solve_y((SD+SU)/2,SU,xp,center_x,center_y,Radius)) 
				{
					f1=Fun_Solve_y((SD+SU)/2,SU,xp,center_x,center_y,Radius);
					f2=Fun_Solve_y(SD,(SD+SU)/2,xp,center_x,center_y,Radius);
				}
				else {
					f1=Fun_Solve_y(SD,(SD+SU)/2,xp,center_x,center_y,Radius);
					f2=Fun_Solve_y((SD+SU)/2,SU,xp,center_x,center_y,Radius);
				}
				//cout<<"f1="<<f1<<" "<<"f2="<<f2<<endl;
				if(y[j]<f2) {eta[i-1][j-1]=0;continue;} //添加=
				if(y[j-1]>f1)     {eta[i-1][j-1]=0;continue;} //添加=
				if(y[j-1]>=f2&&y[j]<f1) 
				{ 
					//	cout<<"RD="<<RD<<"(RD+RU)/2="<<(RD+RU)/2<<"y="<<y[j];
					XFD=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j-1],center_y,Radius);
					XFU=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j],center_y,Radius);
					//	    cout<<"XFD="<<XFD<<" "<<"XFU="<<XFU<<endl;
					double Temp_area,Area;
					if(y[j-1]>=center_y)  { Temp_area=RombergU(XFU,XFD,center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
					Area=Temp_area+abs(x[i-1]-XFU)*(y[j]-y[j-1]);
					//  cout<<"Temp_area'="<<Temp_area<<endl;
					eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
					if(y[j]<=center_y)  { Temp_area=RombergD(XFD,XFU,center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
					Area=-Temp_area+abs(x[i-1]-XFU)*(y[j]-y[j-1]);
					//  cout<<"Temp_area="<<Temp_area<<"x[i]="<<x[i]<<XFU<<endl;
					eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
					if(y[j-1]<center_y&&y[j]>center_y){ double Temp1,Temp2;
					Temp1=RombergU(XFU,RU,center_x,center_y,Radius)-abs(XFU-RU)*center_y+abs(x[i-1]-XFU)*(y[j]-center_y);
					Temp2=-(RombergD(XFD,RU,center_x,center_y,Radius)-abs(XFD-RU)*y[j-1])+abs(x[i-1]-RU)*abs(y[j-1]-center_y);
					Area=Temp1+Temp2;
					eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
				}
				if(y[j-1]<f1&&y[j]>=f1) { double Area;  //!!!!!!!!!!!!!!!!!
				XFD=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j-1],center_y,Radius);
				Area=RombergU(x[i-1],XFD,center_x,center_y,Radius)-abs(XFD-x[i-1])*y[j-1];
				eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
				if(y[j-1]<f2&&y[j]>f2) { double Temp_area,Area;
				XFU=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j],center_y,Radius);
				Temp_area=RombergD(x[i-1],XFU,center_x,center_y,Radius)-abs(XFU-x[i-1])*y[j-1];
				Area=-Temp_area+abs(XFU-x[i-1])*(y[j]-y[j-1]);
				eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			}
			/*----------------------------左端----------------------------*/
			if(x[i-1]>RD&&x[i]<=(RD+RU)/2)
			{   double YLU,YLD;
			YLU=Fun_Solve_y((SD+SU)/2,SU,x[i-1],center_x,center_y,Radius);
			YLD=Fun_Solve_y(SD,(SD+SU)/2,x[i-1],center_x,center_y,Radius);
			double YRU,YRD;
			YRU=Fun_Solve_y((SD+SU)/2,SU,x[i],center_x,center_y,Radius);
			YRD=Fun_Solve_y(SD,(SD+SU)/2,x[i],center_x,center_y,Radius);
			if(y[j]<min(YLD,YRD)){eta[i-1][j-1]=0;continue;}
			if(y[j-1]>=max(YLU,YRU)){eta[i-1][j-1]=0;continue;}
			if(y[j-1]>=max(YLD,YRD)&&y[j]<min(YLU,YRU)) {eta[i-1][j-1]=1;}
			if(y[j-1]<=YLU&&y[j]>YLU){ double XFU,XFD;
			if(y[j]<YRU)
			{
				XFU=Fun_Solve_x(RD,(RD+RU)/2,center_x,y[j],center_y,Radius);
				XFD=x[i-1];
				double Temp_area,Area;
				Temp_area=RombergU(XFD,XFU,center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
				Area=Temp_area+abs(x[i]-XFU)*(y[j]-y[j-1]);
				//  cout<<"Temp_area'="<<Temp_area<<endl;
				eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));
			}
			else{XFU=x[i];
			XFD=x[i-1];
			double Temp_area,Area;
			Temp_area=RombergU(XFD,XFU,center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=Temp_area+abs(x[i]-XFU)*(y[j]-y[j-1]);
			//  cout<<"Temp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}}
			if(y[j-1]>=YLU&&y[j]<=YRU){ double XFD,XFU;
			XFD=Fun_Solve_x(RD,(RD+RU)/2,center_x,y[j-1],center_y,Radius);
			XFU=Fun_Solve_x(RD,(RD+RU)/2,center_x,y[j],center_y,Radius);
			double Temp_area,Area;
			Temp_area=RombergU(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=Temp_area+abs(x[i]-XFU)*(y[j]-y[j-1]);
			// cout<<"Temp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			if(y[j-1]>YLU&&y[j]>=YRU){ double XFD,XFU;
			XFD=Fun_Solve_x(RD,(RD+RU)/2,center_x,y[j-1],center_y,Radius);
			XFU=x[i];
			double Temp_area,Area;
			Temp_area=RombergU(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=Temp_area+abs(x[i]-XFU)*(y[j]-y[j-1]);
			//  cout<<"Temp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			/*-----------------------------------------------------------------*/
			if(y[j-1]<=YLD&&y[j]>YLD){ double XFD,XFU;
			if(y[j-1]<YRD)
			{
				XFD=x[i];
				XFU=x[i-1];
				double Temp_area,Area;
				Temp_area=RombergD(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFD-XFU)*y[j-1];
				Area=-Temp_area+abs(x[i-1]-x[i])*(y[j]-y[j-1]);
				// cout<<"Temp_area'="<<Temp_area<<endl;
				eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));} 
			//cout<<"ext="<<eta[i-1][j-1];
			else { XFD=Fun_Solve_x(RD,(RD+RU)/2,center_x,y[j-1],center_y,Radius);
			double Temp_area,Area;
			Temp_area=RombergD(x[i-1],XFD,center_x,center_y,Radius)-abs(XFD-x[i-1])*y[j-1];
			Area=-Temp_area+abs(x[i-1]-x[i])*(y[j]-y[j-1]);
			//    cout<<"Temp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}}  
			if(y[j-1]>=YRD&&y[j]<=YLD){ double XFD,XFU;
			XFD=Fun_Solve_x(RD,(RD+RU)/2,center_x,y[j-1],center_y,Radius);
			XFU=Fun_Solve_x(RD,(RD+RU)/2,center_x,y[j],center_y,Radius);
			double Temp_area,Area;
			Temp_area=RombergD(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=-Temp_area+abs(x[i]-XFU)*(y[j]-y[j-1]);
			//    cout<<"Temp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			if(y[j-1]<=YRD&&y[j]<YLD&&y[j]>YRD){ double XFD,XFU;
			XFU=Fun_Solve_x(RD,(RD+RU)/2,center_x,y[j],center_y,Radius);
			XFD=x[i];
			double Temp_area,Area;
			Temp_area=RombergD(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=-Temp_area+abs(XFD-XFU)*(y[j]-y[j-1]);
			//   cout<<"Temp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			}
			/*--------------------右端--------------------------------------*/
			if(x[i-1]>=(RD+RU)/2&&x[i]<RU)
			{   double YLU,YLD;
			YLU=Fun_Solve_y((SD+SU)/2,SU,x[i-1],center_x,center_y,Radius);
			YLD=Fun_Solve_y(SD,(SD+SU)/2,x[i-1],center_x,center_y,Radius);
			double YRU,YRD;
			YRU=Fun_Solve_y((SD+SU)/2,SU,x[i],center_x,center_y,Radius);
			YRD=Fun_Solve_y(SD,(SD+SU)/2,x[i],center_x,center_y,Radius);
			if(y[j]<min(YLD,YRD)) {eta[i-1][j-1]=0;}
			if(y[j-1]>=max(YLU,YRU)) {eta[i-1][j-1]=0;}
			if(y[j-1]>=max(YLD,YRD)&&y[j]<min(YLU,YRU)) {eta[i-1][j-1]=1;}
			if(y[j-1]<=YRU&&y[j]>YRU){ double XFU,XFD;
			if(y[j]<YLU){
				XFU=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j],center_y,Radius);
				double Temp_area,Area;
				Temp_area= RombergU(min(x[i],XFU),max(x[i],XFU),center_x,center_y,Radius)-abs(XFU-x[i])*y[j-1];
				Area=Temp_area+abs(x[i-1]-XFU)*(y[j]-y[j-1]);
				// cout<<"TTTemp_area'="<<Temp_area<<endl;
				eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			else { XFU=x[i-1];
			XFD=x[i];
			double Temp_area,Area;
			Temp_area= RombergU(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=Temp_area+abs(x[i-1]-XFU)*(y[j]-y[j-1]);
			//    cout<<"TTTemp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			}
			if(y[j-1]>=YRU&&y[j]<=YLU){ double XFD,XFU;
			XFD=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j-1],center_y,Radius);
			XFU=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j],center_y,Radius);
			double Temp_area,Area;
			Temp_area=RombergU(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=Temp_area+abs(x[i-1]-XFU)*(y[j]-y[j-1]);
			// cout<<"Temp_area="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			if(y[j-1]>YRU&&y[j-1]<YLU&&y[j]>=YLU){ double XFD,XFU;
			XFD=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j-1],center_y,Radius);
			XFU=x[i-1];
			double Temp_area,Area;
			Temp_area=RombergU(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=Temp_area;
			//   cout<<"ATemp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			/*-----------------------------------------------------------------*/
			if(y[j-1]<=YLD&&y[j]>YLD){ double XFD,XFU;
			if(y[j]<YRD) {
				XFU=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j],center_y,Radius);
				XFD=x[i-1];
				double Temp_area,Area;
				Temp_area=RombergD(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFD-XFU)*y[j-1];
				Area=-Temp_area+abs(XFD-XFU)*(y[j]-y[j-1]);
				//     cout<<"Temp_area'="<<Temp_area<<endl;
				eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			else { XFU=x[i];
			XFD=x[i-1];
			double Temp_area,Area;
			Temp_area=RombergD(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFD-XFU)*y[j-1];
			Area=-Temp_area+abs(XFD-XFU)*(y[j]-y[j-1]);
			//    cout<<"Temp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}}
			if(y[j-1]>=YLD&&y[j]<=YRD){ double XFD,XFU;
			XFD=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j-1],center_y,Radius);
			XFU=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j],center_y,Radius);
			double Temp_area,Area;
			Temp_area=RombergD(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=-Temp_area+abs(x[i-1]-XFU)*(y[j]-y[j-1]);
			//    cout<<"Temp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			if(y[j-1]>YLD&&y[j-1]<YRD&&y[j]>=YRD){ double XFD,XFU;
			XFU=x[i];
			XFD=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j-1],center_y,Radius);
			double Temp_area,Area;
			Temp_area=RombergD(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=-Temp_area+abs(x[i-1]-XFU)*(y[j]-y[j-1]);
			//    cout<<"Temp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			}
		/*-------x[i-1]分布在左边-x[i]分布在右边--------*/
		if(x[i-1]<(RD+RU)/2&&x[i]>(RD+RU)/2)
			{   
		    double YLU,YLD;
			YLU=Fun_Solve_y((SD+SU)/2,SU,x[i-1],center_x,center_y,Radius);
			YLD=Fun_Solve_y(SD,(SD+SU)/2,x[i-1],center_x,center_y,Radius);
			double YRU,YRD;
			YRU=Fun_Solve_y((SD+SU)/2,SU,x[i],center_x,center_y,Radius);
			YRD=Fun_Solve_y(SD,(SD+SU)/2,x[i],center_x,center_y,Radius);
			if(y[j]<min(YLD,YRD)) {eta[i-1][j-1]=0;}
			if(y[j-1]>=max(YLU,YRU)) {eta[i-1][j-1]=0;}
			if(y[j-1]>=max(YLD,YRD)&&y[j]<min(YLU,YRU)) {eta[i-1][j-1]=1;}
			if(y[j-1]<=YRU&&y[j]>YRU){ double XFU,XFD;
			if(y[j]<YLU){
				XFU=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j],center_y,Radius);
				double Temp_area,Area;
				Temp_area= RombergU(min(x[i],XFU),max(x[i],XFU),center_x,center_y,Radius)-abs(XFU-x[i])*y[j-1];
				Area=Temp_area+abs(x[i-1]-XFU)*(y[j]-y[j-1]);
				// cout<<"TTTemp_area'="<<Temp_area<<endl;
				eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			else { XFU=x[i-1];
			XFD=x[i];
			double Temp_area,Area;
			Temp_area= RombergU(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=Temp_area+abs(x[i-1]-XFU)*(y[j]-y[j-1]);
			//    cout<<"TTTemp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			}
			if(y[j-1]>=YRU&&y[j]<=YLU){ double XFD,XFU;
			XFD=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j-1],center_y,Radius);
			XFU=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j],center_y,Radius);
			double Temp_area,Area;
			Temp_area=RombergU(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=Temp_area+abs(x[i-1]-XFU)*(y[j]-y[j-1]);
			// cout<<"Temp_area="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			if(y[j-1]>YRU&&y[j-1]<YLU&&y[j]>=YLU){ double XFD,XFU;
			XFD=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j-1],center_y,Radius);
			XFU=x[i-1];
			double Temp_area,Area;
			Temp_area=RombergU(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=Temp_area;
			//   cout<<"ATemp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			/*-----------------------------------------------------------------*/
			if(y[j-1]<=YLD&&y[j]>YLD){ double XFD,XFU;
			if(y[j]<YRD) {
				XFU=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j],center_y,Radius);
				XFD=x[i-1];
				double Temp_area,Area;
				Temp_area=RombergD(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFD-XFU)*y[j-1];
				Area=-Temp_area+abs(XFD-XFU)*(y[j]-y[j-1]);
				//     cout<<"Temp_area'="<<Temp_area<<endl;
				eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			else { XFU=x[i];
			XFD=x[i-1];
			double Temp_area,Area;
			Temp_area=RombergD(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFD-XFU)*y[j-1];
			Area=-Temp_area+abs(XFD-XFU)*(y[j]-y[j-1]);
			//    cout<<"Temp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}}
			if(y[j-1]>=YLD&&y[j]<=YRD){ double XFD,XFU;
			XFD=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j-1],center_y,Radius);
			XFU=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j],center_y,Radius);
			double Temp_area,Area;
			Temp_area=RombergD(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=-Temp_area+abs(x[i-1]-XFU)*(y[j]-y[j-1]);
			//    cout<<"Temp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			if(y[j-1]>YLD&&y[j-1]<YRD&&y[j]>=YRD){ double XFD,XFU;
			XFU=x[i];
			XFD=Fun_Solve_x((RD+RU)/2,RU,center_x,y[j-1],center_y,Radius);
			double Temp_area,Area;
			Temp_area=RombergD(min(XFD,XFU),max(XFD,XFU),center_x,center_y,Radius)-abs(XFU-XFD)*y[j-1];
			Area=-Temp_area+abs(x[i-1]-XFU)*(y[j]-y[j-1]);
			//    cout<<"Temp_area'="<<Temp_area<<endl;
			eta[i-1][j-1]=Area/((x[i]-x[i-1])*(y[j]-y[j-1]));}
			
			}
		//cout<<"eta[i-1][j-1]="<<eta[i-1][j-1]<<endl;
		}
	
	}
	
	/*ofstream outeta;
	outeta.open("F:\\fluid\\C-V-IBM\\eta.txt");
	cout<<"---the eta---"<<endl;
	for(int j=M-1;j>=0;j--)
	{
		for(int  i=0;i<N;i++)
		{
			cout<<" "<<eta[i][j]; 
			outeta<<setprecision(5)<<eta[i][j]<<" ";
		}
		cout<<endl;
		outeta<<endl;
	}*/
	delete [] x;
	delete [] y;
}

	