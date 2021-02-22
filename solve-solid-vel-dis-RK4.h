#include <iostream>
#include <math.h> 
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include<iostream>
#include<math.h> 
#include <iomanip>
#include <fstream>
using namespace std;


void solrigid(double **X,double **scc,double **svel,double **dis,double **vel,double**inisc,double mass,double damp,double rigid,double rot_iner,double rot_damp,double rot_rigid,double ss,double delta_t,double delta_x,int ni,int kk)
{

    int i,j,k,nk=2;
	double pi=3.14159;
    double K1,K2,K3,K4,L1,L2,L3,L4;
	double disxxx=0,disyyy=0,angleaccc=0,anglevel,theta;
	double aef,angle,angleve,angleac,accex,accey,accexx,acceyy,velx,vely,velangle,fx=0,fy=0,rem=0,disx,disy,disangle,rx,ry,temp,temp1,xx,rr;
   
	 /*  char buffff[256];
	   sprintf(buffff, "C:\\movement\\Re300-2\\solidmov-%d.txt",kk-1); 
	   ifstream infilee;
	   infilee.open(buffff, ios::out);
   		if(!infilee)
		{
			cout<<"read is ok\n";
			exit(1);
		}		
			//infilee >>angleac>>" ">>angleve>>" ">>accex>>" ">>accey>>" ">>velx>>" ">>vely>>" ">>disx>>" ">>disy;
		    infilee >>angleac;
			infilee >>angleve;
			infilee >>accex;
			infilee >>accey;
			infilee >>velx;
            infilee >>vely;
			infilee >>disx;
			infilee >>disy;
	    infilee.close();		
        cout<<angleac<<"   "<<angleve<<endl;
        cout<<accex<<"   "<<accey<<endl;   
        cout<<ss<<"   "<<delta_x<<"   "<<delta_t<<" "<<rot_iner<<endl; */
	/*-----------------求合力矩和合力------------------*/
	    
	    for(i=0;i<ni;i++)
	    {
		    rx=scc[i][0]-dis[0][1];
		    ry=scc[i][1]-dis[1][1];
           // rem+=ss*delta_x*(X[i][0]*ry-X[i][1]*rx);		  
		   // fx+=ss*delta_x*X[i][0];
		    fy+=-2*2*ss*delta_x*X[i][1]/10/pi;
	    }		
		cout<<"fy="<<fy<<endl; 
    /*-----------------求位移和旋转角速度------------------*/
        /* //X方向
		    K1=vel[0][1];
	        L1=-damp/mass*vel[0][1]-rigid/mass*dis[0][1]+fy/mass;
			K2=vel[0][1]+0.5*delta_t*L1;
			L2=-damp/mass*(vel[0][1]+0.5*delta_t*L1)-rigid/mass*(dis[0][1]+0.5*delta_t*K1)+fy/mass;
			K3=vel[0][1]+0.5*delta_t*L2;
		    L3=-damp/mass*(vel[0][1]+0.5*delta_t*L2)-rigid/mass*(dis[0][1]+0.5*delta_t*K2)+fy/mass;
			K4=vel[0][1]+delta_t*L3;
            L4=-damp/mass*(vel[0][1]+delta_t*L3)-rigid/mass*(dis[0][1]+delta_t*K3)+fy/mass;
            disx=dis[0][1]+delta_t*(K1+2*K2+2*K3+K4)/6;
			velx=vel[0][1]+delta_t*(L1+2*L2+2*L3+L4)/6;*/
		
		 //y方向
		
		
		    K1=vel[1][1];
	        L1=-damp/mass*vel[1][1]-rigid/mass*dis[1][1]+fy/mass;
			K2=vel[1][1]+0.5*delta_t*L1;
			L2=-damp/mass*(vel[1][1]+0.5*delta_t*L1)-rigid/mass*(dis[1][1]+0.5*delta_t*K1)+fy/mass;
			K3=vel[1][1]+0.5*delta_t*L2;
		    L3=-damp/mass*(vel[1][1]+0.5*delta_t*L2)-rigid/mass*(dis[1][1]+0.5*delta_t*K2)+fy/mass;
			K4=vel[1][1]+delta_t*L3;
            L4=-damp/mass*(vel[1][1]+delta_t*L3)-rigid/mass*(dis[1][1]+delta_t*K3)+fy/mass;
            disy=dis[1][1]+delta_t*(K1+2*K2+2*K3+K4)/6;
			vely=vel[1][1]+delta_t*(L1+2*L2+2*L3+L4)/6;

           //旋转


		/*	K1=vel[2][1];
	        L1=-rot_damp/rot_iner*vel[2][1]-rot_rigid/rot_iner*dis[2][1]+fy/rot_iner;
			K2=vel[2][1]+0.5*delta_t*L1;
			L2=-rot_damp/rot_iner*(vel[2][1]+0.5*delta_t*L1)-rot_rigid/rot_iner*(dis[2][1]+0.5*delta_t*K1)+fy/rot_iner;
			K3=vel[2][1]+0.5*delta_t*L2;
		    L3=-rot_damp/rot_iner*(vel[2][1]+0.5*delta_t*L2)-rot_rigid/rot_iner*(dis[2][1]+0.5*delta_t*K2)+fy/rot_iner;
			K4=vel[2][1]+delta_t*L3;
            L4=-rot_damp/rot_iner*(vel[2][1]+delta_t*L3)-rot_rigid/rot_iner*(dis[2][1]+delta_t*K3)+fy/rot_iner;
            disangle=dis[2][1]+delta_t*(K1+2*K2+2*K3+K4)/6;
			velangle=vel[2][1]+delta_t*(L1+2*L2+2*L3+L4)/6;*/



		
	/*-----------------求位移和旋转角速度------------------*/	  
		
		   // dis[0][0]=dis[0][1];
           // dis[0][1]=disx;
           //   vel[0][0]=vel[0][1];
           //   vel[0][1]=velx;

              dis[1][0]=dis[1][1];
              dis[1][1]=disy;
			  vel[1][0]=vel[1][1];
              vel[1][1]=vely;

          //  dis[2][0]=dis[2][1];
           // dis[2][1]=disangle;		
		  //    vel[2][0]=vel[2][1];
          //    vel[2][1]=velangle;
            cout<<"y位移dy="<<dis[1][0]<<endl; 
			char buff[256];
			sprintf(buff, "E:\\wwq\\programm-13\\005\\results\\rigidmov-%d.txt", kk); //每隔MemIntStep步，写入到文本中，文件名做相应改变，不能覆盖

			ofstream outfor;
			outfor.open(buff, ios::out);			
            outfor<<dis[1][0]<<"   "<<vel[1][0];
			outfor<<endl;	
	        outfor.close ();	
            cout<<"输出线位移、线速度完成"<<endl; 

				
			/*char buf[256];
			sprintf(buf, "C:\\movement\\Re300-2\\angle-%d.txt", kk); //每隔MemIntStep步，写入到文本中，文件名做相应改变，不能覆盖
			ofstream outf;
			outf.open(buf, ios::out);           		 
		   	outf<<remm<<"   "<<angleaccc<<"   "<<anglevell<<"   "<<thetaa<<"  "<<disxxx;  		     
	        outf<<endl;			
	        outf.close ();	
            cout<<"输出固体节点力完成"<<endl;*/
	       
	
           for(i=0;i<ni;i++)
	       {
            /*rx=inisc[i][0];
		    ry=inisc[i][1];	    				
			xx=ry/rx;
			rr=sqrt(rx*rx+ry*ry);			
			if (rx>0)
	       {
		    aef=atan(xx);
		   }
		   else
		   {
           aef=pi+atan(xx);
	       }
		   angle=dis[2][1]+aef;		  
	       scc[i][0]=dis[0][1]+rr*cos(angle);
		   scc[i][1]=dis[1][1]+rr*sin(angle);		  
		   svel[i][0]=(dis[0][1]-dis[0][0])/delta_t-rr*dis[2][1]*sin(angle);
		   svel[i][1]=(dis[1][1]-dis[1][0])/delta_t+rr*dis[2][1]*cos(angle);*/
		   scc[i][0]=dis[0][1]+inisc[i][0];
		   scc[i][1]=dis[1][1]+inisc[i][1];		  
		   svel[i][0]=vel[0][1];
		   svel[i][1]=vel[1][1];
	     }		

	
}
