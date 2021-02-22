
using namespace std;

double fun_Solve_y(double x,double  center_x,double y,double  center_y,double R )
{
	double z;
	z=(x-center_x)*(x-center_x)+(y-center_y)*(y-center_y)-R*R; //圆的方程
	return z;
}

double Fun_Solve_y(double a,double b,double xp,double center_x,double center_y,double R)    
{
	double ya;
	double yb;
	//double fun_Solve_y(double x,double  center_x,double,double y,double  center_y,double R )
	ya=fun_Solve_y(xp,center_x,a,center_y,R);
	yb=fun_Solve_y(xp,center_x,b,center_y,R);
	double c;
	double INF=0.00001; //容许误差
	int K=1000;   //最大的迭代次数
	double yc;
	int NUM=0;  //计算迭代的次数
	//cout<<"fun-s-y="<<fun_Solve_y(a,center_x,xp,center_y,R);  这个地方xp和a b放错了位置 花了不少时间
	if(abs(fun_Solve_y(xp,center_x,b,center_y,R))<INF||abs(fun_Solve_y(xp,center_x,a,center_y,R))<INF) 
	{
		if(abs(fun_Solve_y(xp,center_x,b,center_y,R))<INF) {c=b;}
		else c=a;
	}
	else {
	for(int i=0;i<K;i++)
	{
		c=(a+b)/2;
		yc=fun_Solve_y(xp,center_x,c,center_y,R);
		NUM++;
		if(abs(yc)<INF) { break;}
		else if(ya*yc<0) { b=c; yb=yc;}
		else { a=c ; ya=yc;}
	}
	if(abs(yc)>INF) {cout<<"No root in this Region !"<<endl;}
	//cout<<"NUM="<<NUM<<endl;
	
	}
	//cout<<"c="<<c<<endl;
	return c;
}
double fun_Solve_x(double x,double  center_x,double y,double  center_y,double R )
{
	double z;
	z=(x-center_x)*(x-center_x)+(y-center_y)*(y-center_y)-R*R; //圆的方程
	return z;
}

double Fun_Solve_x(double a,double b,double center_x,double yp,double center_y,double R)    
{
	double xa;
	double xb;
	//double fun_Solve_x(double x,double  center_x,double,double y,double  center_y,double R )
	xa=fun_Solve_x(a,center_x,yp,center_y,R);
	xb=fun_Solve_x(b,center_x,yp,center_y,R);
	double c;
	double INF=0.00001; //容许误差
	int K=1000;   //最大的迭代次数
	double xc;
	int NUM=0;  //计算迭代的次数
	if(abs(fun_Solve_x(b,center_x,yp,center_y,R))<INF||abs(fun_Solve_x(a,center_x,yp,center_y,R))<INF) 
	{
		if(abs(fun_Solve_x(b,center_x,yp,center_y,R))<INF) {c=b;}
		else c=a;
	}
	else {
		for(int i=0;i<K;i++)
		{
			c=(a+b)/2;
			xc=fun_Solve_x(c,center_x,yp,center_y,R);
			NUM++;
			if(abs(xc)<INF) { break;}
			else if(xa*xc<0) { b=c; xb=xc;}
			else { a=c ; xa=xc;}
		}
		if(abs(xc)>INF) {cout<<"No root in this Region !"<<endl;}
		//cout<<"NUM="<<NUM<<endl;
		//cout<<"c="<<c<<endl;
	}
	return c;
}
double fun_U_semicircle(double x,double center_x,double center_y,double R)      //定义函数f(x)
{
	double r;
	r=center_y+sqrt(R*R-(x-center_x)*(x-center_x));
	return r;
}

double fun_D_semicircle(double x,double center_x,double center_y,double R )   
{
	double r;
	r=center_y-sqrt(R*R-(x-center_x)*(x-center_x));
	return r;
}

double RombergU(double aa,double bb,double center_x,double center_y,double R)            //Romberg（龙贝格)公式
{
	double re;
	double a[4][4];
	
	//double aa,bb;         //aa,bb为端点
	double h,sum;       //h为区间的长度
	double p;            //p为分点
	int i,j,k;
	//printf("请输入端点aa bb\n(注：格式为aa bb)\n");
	//scanf("%lf %lf",&aa,&bb);
	//if(aa=0){}
    sum=fun_U_semicircle(aa,center_x,center_y,R)+fun_U_semicircle(bb,center_x,center_y,R);
	a[0][0]=((bb-aa)/2)*sum;
	for(k=1;k<4;k++)         //k为二分的次数
	{
		sum=fun_U_semicircle(aa,center_x,center_y,R)+fun_U_semicircle(bb,center_x,center_y,R);
		h=(bb-aa)/pow(2.0,k);
	    for(i=1;i<pow(2.0,k);i++)      //计算除端点外的分点和的2倍
		{
	    	p=aa+i*h;                //p为分点
	    	sum=sum+2*fun_U_semicircle(p,center_x,center_y,R);
		}
            a[k][0]=(h/2.0)*sum;     //计算第一列       
	}
			 
	for(j=1;j<4;j++)              //i为行，j为列
		for(i=j;i<4;i++)
	        a[i][j]=(pow(4.0,j)/(pow(4.0,j)-1))*a[i][j-1]-(1/(pow(4.0,j)-1))*a[i-1][j-1];  //计算后面三列
       
		return re=a[3][3];
}
double RombergD(double aa,double bb,double center_x,double center_y,double R)            //Romberg（龙贝格)公式
{
	double re;
	double a[4][4];
	
	//double aa,bb;         //aa,bb为端点
	double h,sum;       //h为区间的长度
	double p;            //p为分点
	int i,j,k;
	//printf("请输入端点aa bb\n(注：格式为aa bb)\n");
	//scanf("%lf %lf",&aa,&bb);
	//if(aa=0){}
    sum=fun_D_semicircle(aa,center_x,center_y,R)+fun_D_semicircle(bb,center_x,center_y,R);
	a[0][0]=((bb-aa)/2)*sum;
	for(k=1;k<4;k++)         //k为二分的次数
	{
		sum=fun_D_semicircle(aa,center_x,center_y,R)+fun_D_semicircle(bb,center_x,center_y,R);
		h=(bb-aa)/pow(2.0,k);
	    for(i=1;i<pow(2.0,k);i++)      //计算除端点外的分点和的2倍
		{
	    	p=aa+i*h;                //p为分点
	    	sum=sum+2*fun_D_semicircle(p,center_x,center_y,R);
		}
            a[k][0]=(h/2.0)*sum;     //计算第一列       
	}
			 
	for(j=1;j<4;j++)              //i为行，j为列
		for(i=j;i<4;i++)
	        a[i][j]=(pow(4.0,j)/(pow(4.0,j)-1))*a[i][j-1]-(1/(pow(4.0,j)-1))*a[i-1][j-1];  //计算后面三列
       
		return re=a[3][3];
}