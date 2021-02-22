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


void solrigid(double **X,double **scc,double **svel,double **dis,double**inisc,double mass,double damp,double rigid,double rot_iner,double rot_damp,double rot_rigid,double ss,double delta_t,double delta_x,int ni,int kk)
{

    int i,j,k,nk=2;
	double pi=3.14159;
    double con1=0.5,con2=0.25;
	double disxxx=0,disyyy=0,rem=0,angleaccc=0,anglevel,theta;
	double aef,angle,angleve,angleac,accex,accey,accexx,acceyy,velx,vely,velxx,velyy,fx=0,fy=0,disx,disy,rx,ry,temp,temp1,xx,rr;
   
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
		    fy+=2*2*ss*delta_x*X[i][1]/10/pi;
	    }		
		cout<<rem<<"   "<<fx<<"   "<<fy<<endl; 
    /*-----------------求位移和旋转角速度------------------*/
	       // disx=dis[0][1]*(2-damp/mass*delta_t-rigid/mass*delta_t*delta_t)+dis[0][0]*(damp/mass*delta_t-1)+delta_t*delta_t*fx/mass;//求得下一时刻x方向位移disx	
           // velx=(disx-dis[0][1])/delta_t;//求得下一时刻x方向线速度
			disy=dis[1][1]*(2-damp/mass*delta_t-rigid/mass*delta_t*delta_t)+dis[1][0]*(damp/mass*delta_t-1)+delta_t*delta_t*fy/mass;//求得下一时刻y方向位移disx	
            vely=(disy-dis[1][1])/delta_t;//求得下一时刻y方向线速度
           // theta=dis[2][1]*(2-rot_damp/rot_iner*delta_t-rot_rigid/rot_iner*delta_t*delta_t)+dis[2][0]*(rot_damp/rot_iner*delta_t-1)+delta_t*delta_t*rem/rot_iner;///求得下一时刻转过的角度theta
          //  anglevel=(theta-dis[2][1])/delta_t;//求得下一时刻角速度
            //angleacc=(thetaaa-2*thetaa+theta)/delta_t/delta_t;//求得下一时刻角加速度
	/*-----------------求位移和旋转角速度------------------*/	  
		
		   // dis[0][0]=dis[0][1];
           // dis[0][1]=disx;
            dis[1][0]=dis[1][1];
            dis[1][1]=disy;
          //  dis[2][0]=dis[2][1];
           // dis[2][1]=anglevel;
		
		   /* accexx=fx;
	        acceyy=fy;           
	        angleaccc=remm/rot_iner;//求得下一时刻角加速度
	        anglevell=angleve+((1-con1)*angleac+con1*angleaccc)*delta_t;//求得下一时刻角速度
	        velxx=velx+((1-con1)*accex+con1*accexx)*delta_t;
	        velyy=vely+((1-con1)*accey+con1*acceyy)*delta_t;
	        thetaa=angleve*delta_t+((0.5-con2)*angleac+con2*thetaa)*delta_t*delta_t;       //求得转过的角度theta
            disxxx=velx*delta_t+((0.5-con2)*accex+con2*accexx)*delta_t*delta_t;
	        disyyy=vely*delta_t+((0.5-con2)*accey+con2*acceyy)*delta_t*delta_t;

			char bufff[256];
			sprintf(bufff, "C:\\movement\\Re300-2\\solidmov-%d.txt",kk); //每隔MemIntStep步，写入到文本中，文件名做相应改变，不能覆盖

			ofstream outtforcee;
			outtforcee.open(bufff, ios::out);
			outtforcee <<angleaccc;
			outtforcee <<anglevell;
			outtforcee <<accexx;
			outtforcee <<acceyy;
			outtforcee <<velxx;
            outtforcee <<velyy;
			outtforcee <<disxxx;
			outtforcee <<disyyy;
             //outtforcee<<angleaccc<<"   "<<anglevell<<"   "<<accexx<<"   "<<acceyy<<"   "<<velxx<<"   "<<velyy<<"   "<<disxxx<<"   "<<disyyy; 
			        
	        outtforcee.close ();	
            cout<<"输出角加速度、角速度完成"<<endl; */

			char buff[256];
			sprintf(buff, "f:\\wwq\\programm-4\\cylinder\\Re200-5\\rigidmov-%d.txt", kk); //每隔MemIntStep步，写入到文本中，文件名做相应改变，不能覆盖

			ofstream outfor;
			outfor.open(buff, ios::out);			
            outfor<<dis[1][0]<<"   "<<vely;
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
		   svel[i][0]=(dis[0][1]-dis[0][0])/delta_t;
		   svel[i][1]=(dis[1][1]-dis[1][0])/delta_t;
	     }		

	
}
