#include<iostream>
#include<math.h>
#include <iomanip>
#include <fstream>
#include <Eigen/Sparse>
//#include "umfpack.h"
//#include <Eigen/src/UmfPackSupport/UmfPackSupport.h>
#include "Solve-ustar.h"
#include "Solve-vstar.h"
#include "Solver-Pressure.h"
#include "SUPie.h"
#include "SVPie.h"
#include "SVelPlusOne-1.h"
#include "SVelPlusOne-2.h"
//#include "solve-solid-vel-dis.h"
#include "solve-solid-vel-dis-RK4.h"


/*===========================add====添加开始====add================================*/
#include "dirent.h"
/*添加原因：这个头文件用来实现目录操作，需要自己添加*/
/*==================================添加结束=======================================*/


using namespace Eigen;

using namespace std;

//void (*Sol_Press)(int ,int ,double  ,double ,double ,double ** ,double ** ,double ** ,double **,double **); //定义一个函数指针变量 Sol_Press

//void (*Sol_ustar)(double ** ,double ** ,double ** ,int , int ); //定义一个函数指针变量 Sol_ustar
void (*Sol_ustar)(double ** ,double ** ,double ** ,double ,double ,double ,double ,double ,double ,int ,int );
//void (*Sol_vstar)(double ** ,double ** ,double ** ,int , int ); //定义一个函数指针变量 Sol_vstar
void (*Sol_vstar)(double ** ,double ** ,double ** ,double  ,double  ,double  ,double  ,double  ,double ,int ,int );

void (*Sol_upie)(double **ustar,double **upie,double **Pres,int N,int M ); //定义一个函数指针变量 Sol_upie
void (*Sol_vpie)(double **ustar,double **upie,double **Pres,int N,int M ); //定义一个函数指针变量 Sol_vpie

void (*Sol_VelPlusOne)(double **u,double **X,double **upie,double **v,double **vpie,double **iniscc_1,double **inisvell_1,double **inifcc_1,int N,int M,int NN_1,int MM_1,int N_sec_1,int M_sec_2,double delta_xx,double delta_yy,double ss_1,double delta_t,int nii_1,int njj_1);
void (*Sol_VelPlusTwo)(double **u,double **X,double **v,double **iniscc_2,double **inisvell_2,double **inifcc_2,int N,int M,int NN_2,int MM_2,int N_sec_2,int M_sec_2,double delta_xx,double delta_yy,double ss_2,double delta_t,int nii_2,int njj_2);

void (*Sol_rigid)(double **X2,double **scc,double **svel,double **dis,double **vel,double **inisc_2,double mass,double damp,double rigid,double rot_iner,double rot_damp,double rot_rigid,double ss_2,double delta_t,double delta_x,int nii_2,int k);

//iniscc－固体节点初始位置坐标，inisvell-固体节点初始速度，inifcc－流体次区域节点初始坐标（流体次区域,NN－流体在X方向单元数，其在X方向节点为NN+1;MM－流体在Y方向单元数，其在Y方向节点为MM+1）
//N－流体在X方向单元数，其在X方向节点为N+1；M－流体在Y方向单元数，其在Y方向节点为M+1;delta_x 为X方向流体网格间距,delta_y为Y方向流体网格间距,ttot为计算时间步长,nii为固体边界节点数,njj为流体次区域流体节点总数
typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double


static string __base_address = "/content/wwx/resource/";
static string base_address = __base_address;
/*===========================add====添加开始====add================================*/

/*用来在state目录查找最新的状态文件的序号*/
static int find_current_index(){
    DIR *dir;
    base_address = __base_address;
    if ( !(dir=opendir(base_address.append("state").data()))){
        return -1;
    }
    char current_index_buf[256];
    char state_file_name_buf[256];
    int final_index = -1;
    int current_index = -1;
    int in_top_two_new_state = false;
    struct dirent *stdinfo;
    int top_two_new_state[2]={-1,-1};
    while (1) {
        if ((stdinfo = readdir(dir)) == 0) {
            break;
        }
        if(stdinfo->d_type == DT_REG && (strncmp(stdinfo->d_name,"state_",6) == 0)){
            strcpy(current_index_buf,stdinfo->d_name + 6);
            current_index = (int)strtol(current_index_buf,NULL,10);
            /*找到最新的两个state*/
            if(current_index >= top_two_new_state[0]){
                if(top_two_new_state[0] > top_two_new_state[1]){
                    top_two_new_state[1] = top_two_new_state[0];
                }
                top_two_new_state[0] = current_index;
                continue;
            }
            if(current_index >= top_two_new_state[1]){
                top_two_new_state[1] = current_index;
            }
        }
    }
    base_address = __base_address;
    if ( !(dir=opendir(base_address.append("state").data()))){
        return -1;
    }
    while (1) {
        if ((stdinfo = readdir(dir)) == 0) {
            break;
        }
        if(stdinfo->d_type == DT_REG && (strncmp(stdinfo->d_name,"state_",6) == 0)){
            in_top_two_new_state = false;
            strcpy(current_index_buf,stdinfo->d_name + 6);
            current_index = (int)strtol(current_index_buf,NULL,10);
            /*是否是最新的个stat之一e*/
            for(int i = 0;i < 2;i++){
                if(current_index == top_two_new_state[i]){
                    in_top_two_new_state = true;
                    break;
                }
            }
            /*如果不是最新的两个state之一，删除它*/
            if(!in_top_two_new_state){
                base_address = __base_address;
                sprintf(state_file_name_buf,base_address.append("state/%s").data(),stdinfo->d_name);
                cout<<"删除旧文件："<<stdinfo->d_name<<endl;
                remove(state_file_name_buf);
            }
        }
    }
    return top_two_new_state[0];
}


static void save_state(int current_index,double** ustar,double **vstar,
                       double **u,double **v,double **upie,double **vpie,double **Pres,
                       int *SparseMatAi,double *SparseMatAx,int *SparseMatAp,
                       double ** X2,double **scc,double **svel,double **inifcc_2,
                       int N,int M,int SMAnZ,int SMADi,int nii_2,int njj_2,int nkk,
                       double delta_t,double delta_x,double delta_y,
                       double CF,double Re,double Turb,
                       int NN_2 ,int MM_2,int N_sec_2,int M_sec_2,double ss_2){
    ofstream outtforce;
    char name_buff[256];
    cout<<"开始保存状态:state_"<<current_index<<endl;
    base_address = __base_address;
    sprintf(name_buff,base_address.append("state/state_%d").data(),current_index);
    outtforce.open(name_buff, ios::out);
    if(!outtforce){
        cout<<"创建状态文件失败"<<name_buff<<endl;
        exit(-1);
    }
    //预留标志文件完整性的魔数
    outtforce<<"_____"<<endl;
    //保存小变量(依照传参顺序)
    outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1) <<N<<" "<<M<<" "<<SMAnZ<<" "<<SMADi<<" "<<nii_2<<" "<<njj_2<<" "<<nkk<<" "<<
             delta_t<<" "<<delta_x<<" "<<delta_y<<" "<<
             CF<<" "<<Re<<" "<<Turb<<" "<<
             NN_2<<" "<<MM_2<<" "<<N_sec_2<<" "<<M_sec_2<<" "<<ss_2<<" "<<endl;

    //保存ustar
    for (int i = 0; i < N + 1; i++)
        for (int j = 0; j < M + 1; j++)
            outtforce<<ustar[i][j] << " ";
    outtforce << endl;
    //保存vstar
    for (int i = 0; i < N + 1; i++)
        for (int j = 0; j < M + 1; j++)
            outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1)<<vstar[i][j] << " ";
    outtforce << endl;
    //保存u
    for (int i = 0; i < N + 1; i++)
        for (int j = 0; j < M + 1; j++)
            outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1)<<u[i][j] << " ";
    outtforce << endl;
    //保存v
    for (int i = 0; i < N + 1; i++)
        for (int j = 0; j < M + 1; j++)
            outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1)<<v[i][j] << " ";
    outtforce << endl;
    //保存upie
    for (int i = 0; i < N + 1; i++)
        for (int j = 0; j < M + 1; j++)
            outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1)<<upie[i][j] << " ";
    outtforce << endl;
    //保存vpie
    for (int i = 0; i < N + 1; i++)
        for (int j = 0; j < M + 1; j++)
            outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1)<<vpie[i][j] << " ";
    outtforce << endl;
    //保存Pres
    for (int i = 0; i < N + 1; i++)
        for (int j = 0; j < M + 1; j++)
            outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1)<<Pres[i][j] << " ";
    outtforce << endl;
    //保存SparseMatAi
    for (int i=0;i<SMAnZ;i++){
        outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1)<<SparseMatAi[i]<< " ";
    }
    outtforce<<endl;
    //保存SparseMatAx
    for (int i=0;i<SMAnZ;i++){
        outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1)<<SparseMatAx[i]<< " ";
    }
    outtforce<<endl;
    //保存SparseMatAp
    for (int i=0;i<SMADi+1;i++){
        outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1)<<SparseMatAp[i]<< " ";
    }
    outtforce<<endl;

    //保存X2
    for (int i = 0; i < nii_2; i++)
        for (int j = 0; j < nkk; j++)
            outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1)<<X2[i][j] << " ";
    outtforce << endl;
    //保存scc
    for (int i = 0; i < nii_2; i++)
        for (int j = 0; j < nkk; j++)
            outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1)<<scc[i][j] << " ";
    outtforce << endl;
    //保存svel
    for (int i = 0; i < nii_2; i++)
        for (int j = 0; j < nkk; j++)
            outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1)<<svel[i][j] << " ";
    outtforce << endl;
    //保存inifcc_2
    for (int i = 0; i < njj_2; i++)
        for (int j = 0; j < nkk; j++)
            outtforce<<std::setprecision (std::numeric_limits<double>::digits10 + 1)<<inifcc_2[i][j] << " ";
    outtforce << endl;
    //回到文件头，写标志文件完整的魔数
    outtforce.seekp(0,ios::beg);
    outtforce << "magic";
    outtforce.close();
    cout<<"保存状态结束"<<endl;
    //调用一下这个函数，可以清理旧的state文件
    find_current_index();
}

static int load_state(double** ustar,double **vstar,
                      double **u,double **v,double **upie,double **vpie,double **Pres,
                      int *SparseMatAi,double *SparseMatAx,int *SparseMatAp,
                      double ** X2,double **scc,double **svel,double **inifcc_2,
                      int *N,int *M,int *SMAnZ,int *SMADi,int *nii_2,int *njj_2,int *nkk,
                      double *delta_t,double *delta_x,double *delta_y,
                      double *CF,double *Re,double *Turb,
                      int *NN_2 ,int *MM_2,int *N_sec_2,int *M_sec_2,double *ss_2){
    ifstream infile;
    char name_buff[256];
    char magic_buff[256];
    int current_index;
    current_index = find_current_index();
    cout<<"开始恢复状态:state_"<<current_index<<endl;
    if(current_index == -1){
        cout<<"错误:未找到保存的状态"<<endl;
        return -1;
    }
    base_address = __base_address;
    sprintf(name_buff,base_address.append("state/state_%d").data(),current_index);
    infile.open(name_buff, ios::in);
    infile >> magic_buff;
    if(strncmp(magic_buff,"magic",5) != 0){
        infile.close();
        cout<<"错误:最新的保存状态不完整"<<endl;
        return -1;
    }
    //检查恢复小变量
    infile>>std::setprecision (std::numeric_limits<double>::digits10 + 1)>>*N>>*M>>*SMAnZ>>*SMADi>>*nii_2>>*njj_2>>*nkk>>
          *delta_t>>*delta_x>>*delta_y>>
          *CF>>*Re>>*Turb>>
          *NN_2>>*MM_2>>*N_sec_2>>*M_sec_2>>*ss_2;

    //加载ustar
    for (int i = 0; i < *N + 1; i++)
        for (int j = 0; j < *M + 1; j++)
            infile >>std::setprecision (std::numeric_limits<double>::digits10 + 1) >> ustar[i][j];

    //加载vstar
    for (int i = 0; i < *N + 1; i++)
        for (int j = 0; j < *M + 1; j++)
            infile >> std::setprecision (std::numeric_limits<double>::digits10 + 1) >> vstar[i][j];

    //加载u
    for (int i = 0; i < *N + 1; i++)
        for (int j = 0; j < *M + 1; j++)
            infile >> std::setprecision (std::numeric_limits<double>::digits10 + 1) >> u[i][j];

    //加载v
    for (int i = 0; i < *N + 1; i++)
        for (int j = 0; j < *M + 1; j++)
            infile >> std::setprecision (std::numeric_limits<double>::digits10 + 1) >> v[i][j];

    //加载upie
    for (int i = 0; i < *N + 1; i++)
        for (int j = 0; j < *M + 1; j++)
            infile >> std::setprecision (std::numeric_limits<double>::digits10 + 1) >> upie[i][j];

    //加载vpie
    for (int i = 0; i < *N + 1; i++)
        for (int j = 0; j < *M + 1; j++)
            infile >> std::setprecision (std::numeric_limits<double>::digits10 + 1) >> vpie[i][j];

    //加载Pres
    for (int i = 0; i < *N + 1; i++)
        for (int j = 0; j < *M + 1; j++)
            infile >> std::setprecision (std::numeric_limits<double>::digits10 + 1) >> Pres[i][j];

    //加载SparseMatAi
    for (int i=0;i < *SMAnZ;i++){
        infile >> std::setprecision (std::numeric_limits<double>::digits10 + 1) >> SparseMatAi[i];
    }

    //加载SparseMatAx
    for (int i=0;i < *SMAnZ;i++){
        infile >> std::setprecision (std::numeric_limits<double>::digits10 + 1) >> SparseMatAx[i];
    }

    //加载SparseMatAp
    for (int i = 0;i < *SMADi + 1;i++){
        infile >> std::setprecision (std::numeric_limits<double>::digits10 + 1) >> SparseMatAp[i];
    }


    //加载X2
    for (int i = 0; i < *nii_2; i++)
        for (int j = 0; j < *nkk; j++)
            infile >> std::setprecision (std::numeric_limits<double>::digits10 + 1) >> X2[i][j];

    //加载scc
    for (int i = 0; i < *nii_2; i++)
        for (int j = 0; j < *nkk; j++)
            infile >> std::setprecision (std::numeric_limits<double>::digits10 + 1) >> scc[i][j];

    //加载svel
    for (int i = 0; i < *nii_2; i++)
        for (int j = 0; j < *nkk; j++)
            infile >> std::setprecision (std::numeric_limits<double>::digits10 + 1) >> svel[i][j];

    //加载inifcc_2
    for (int i = 0; i < *njj_2; i++)
        for (int j = 0; j < *nkk; j++)
            infile >> std::setprecision (std::numeric_limits<double>::digits10 + 1) >> inifcc_2[i][j];
    infile.close();
    cout<<"恢复状态结束"<<endl;
    return current_index;
}
/*添加原因：这部分是用来实现保存与恢复的函数*/
/*==================================添加结束=======================================*/


int main(int argc,char* argv[])
{
    base_address = argv[1];
    std::ios::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);


/*==================================建议开始=======================================*/
    int N=1200;
    cout<<"Please the row"<<endl;//输入行（X方向）
    //cin>>N;
    int M=800;
/*这里的M、N好像是和fluidfield-20000那个文件的“长”“宽”相匹配的，如果改成不带序号的
fluidfield，这两个值应该要改变*/
/*==================================建议结束=======================================*/


    cout<<"Please the lie"<<endl;//输入列（Y方向）
    //cin>>M;
    //int N_sec_1=344;
    int N_sec_2=344;
    cout<<"Please the second row index"<<endl;//输入次流体区域所在起始行（X方向）
    //cin>>N_sec;
    //int M_sec_1=296;
    int M_sec_2=344;
    cout<<"Please the second lie index"<<endl;//输入次流体区域所在起始列（Y方向）
    //cin>>M_sec;
    double Re=100;
    cout<<"Please the Re"<<endl;//输入Re
    //cin>>Re;
    int nii_1;
    cout<<"Please the first solid node "<<endl;//输入固体边界节点总数
    cin>>nii_1;
    int nii_2;
    cout<<"Please the second solid node "<<endl;//输入固体边界节点总数
    cin>>nii_2;
    double Length=15;   //求解区域的长
    double Wide=10;      //求解区域的宽
    //double Seclength_1=2.6;   //求解次区域的长
    //double Secwide_1=2.6;      //求解次区域的宽
    double Seclength_2=1.4;   //求解次区域的长
    double Secwide_2=1.4;      //求解次区域的宽
    double delta_x=Length/N;
    double delta_y=Wide/M;
    //int MM_1=Seclength_1/delta_x;
    //int NN_1=Secwide_1/delta_y;
    int MM_2=Seclength_2/delta_x;
    int NN_2=Secwide_2/delta_y;
    double delta_t=0.001;
    //int njj_1=(MM_1+1)*(NN_1+1);
    int njj_2=(MM_2+1)*(NN_2+1);
    int nkk=2,ntt=2,dof=3;//ntt为存储时间步数，dof固体自由度，平面为3
    int i,j,k,kk,kkk;
    double Spressure=0.0,U=1.0,ss_1=0.05,ss_2=0.024933275/2,pi=3.14159,ps=1.0,fcorx,fcory,ddd,solidcorx,solidcory,solidforcex,solidforcey;
    double currenttime,tt,rr,seta,xx,aef,angle,ft,fn;
    double mass=1,damp=0.02*2*pi/5.1,rigid=(2*pi/5.1)*(2*pi/5.1),rot_iner=0.005081548,rot_damp=0.0004,rot_rigid=1;//ps为固体密度
    double disxx=0,disyy=0;	//ss=0.024933275
    double anglevel,density=1,L=1.1928,L0=0.5,L1=1.75,seta_0=2*pi*45/360,beta=1.0,omegaa=2*pi*0.10,T=2*pi/omegaa, aero=(L0*L0+L1*L1)*seta_0,moment,cm,power,cp;
    double Turb=1.0, CF=0.125;   //若考虑湍流效应，Turb=1.0，否则为Turb=0.0


    cout<<"读入固体节点初始坐标inisc[i][j]"<<endl;
    double **inisc_1=new double *[nii_1];
    for(int k=0;k<nii_1;k++)
        inisc_1[k]= new double [nkk];

    double **inisc_2=new double *[nii_2];
    for(int k=0;k<nii_2;k++)
        inisc_2[k]= new double [nkk];


    ifstream infile;
    base_address = __base_address;
    infile.open(base_address.append("guidevane-node-cor-ini-2D-103.txt").data(),ios::in);

    if(!infile){

        cout<<"不能打开solid-node-cor-ini.txt数据文件!\n";
        exit(1);
    }
    for (i=0; i<nii_1; i++)                      //读入初始构型下first固体节点初始坐标inisc[i][j]
    {
        for (j=0; j<nkk; j++)
        {
            infile>>inisc_1[i][j];
        }
    }
    infile.close();
    base_address = __base_address;
    infile.open(base_address.append("solid-node-cor-ini-2D-252-0.txt").data(),ios::in);

    if(!infile){

        cout<<"不能打开solid-node-cor-ini.txt数据文件!\n";
        exit(1);
    }
    for (i=0; i<nii_2; i++)                      //读入初始构型下second固体节点初始坐标inisc[i][j]
    {
        for (j=0; j<nkk; j++)
        {
            infile>>inisc_2[i][j];
        }
    }
    infile.close();
    cout<<"读入固体节点初始坐标inisc[i][j]结束"<<endl;

    /*ofstream outfile;
    outfile.open("f:\\wwq\\programm-13\\C-V-IBM-60\\solid-node-cor.dat",ios::out);
    if(!outfile){

            cout<<"不能打开solid-node-cor.txt数据文件!\n";
            exit(1);
            }

     for (i=0; i<nii_2; i++)
      {
          for (j=0; j<nkk; j++)
          {

              outfile<<inisc_2[i][j] <<"   ";

          }
          outfile<<endl;
      }
      outfile <<endl;
      outfile.close ();	*/



    cout<<"生成流体次区域节点初始坐标inifcc[i][j]"<<endl;
    double temp1,temp2;
    int n,m;
    /*double **inifcc_1=new double *[njj_1];
  for(int k=0;k<njj_1;k++)
  inifcc_1[k]= new double [nkk];

  for(int i=0;i<NN_1+1;i++)
  {

      for(int j=0;j<MM_1+1;j++)
      {

          temp1=-0.7+i*delta_x;
          temp2=-1.3+j*delta_y;
          n=i*(MM_1+1)+j;
          inifcc_1[n][0]=temp1;
          inifcc_1[n][1]=temp2;
      }
  }*/


    double **inifcc_2=new double *[njj_2];
    for(int k=0;k<njj_2;k++)
        inifcc_2[k]= new double [nkk];

    for(int i=0;i<NN_2+1;i++)
    {

        for(int j=0;j<MM_2+1;j++)
        {

            temp1=-0.7+i*delta_x;
            temp2=-0.7+j*delta_y;
            n=i*(MM_2+1)+j;
            inifcc_2[n][0]=temp1;
            inifcc_2[n][1]=temp2;
        }
    }
    cout<<"计算流体次区域节点初始坐标inifcc[i][j]完成"<<endl;
    double **iniscc_1=new double *[nii_1];//每一时刻固体节点位置
    for(int k=0;k<nii_1;k++)
        iniscc_1[k]= new double [nkk];


    double **scc=new double *[nii_2];//每一时刻固体节点位置
    for(int k=0;k<nii_2;k++)
        scc[k]= new double [nkk];

    double **inisvell_1=new double *[nii_1];//每一时刻固体节点速度
    for(int k=0;k<nii_1;k++)
        inisvell_1[k]= new double [nkk];

    double **svel=new double *[nii_2];//每一时刻固体节点速度
    for(int k=0;k<nii_2;k++)
        svel[k]= new double [nkk];

    double **X1=new double *[nii_1];
    for(int k=0;k<nii_1;k++)
        X1[k]= new double [nkk];

    double **X2=new double *[nii_2];
    for(int k=0;k<nii_2;k++)
        X2[k]= new double [nkk];

    double **dis=new double *[dof];//每一时刻固体位置，ntt为存储的时间步数
    for(int k=0;k<dof;k++)
        dis[k]= new double [ntt];
    for(int i=0;i<dof;i++)
        for(int j=0;j<ntt;j++)
            dis[i][j]=0;

    double **vel=new double *[dof];//每一时刻固体速度，ntt为存储的时间步数
    for(int k=0;k<dof;k++)
        vel[k]= new double [ntt];
    for(int i=0;i<dof;i++)
        for(int j=0;j<ntt;j++)
            vel[i][j]=0;

    double **Pres=new double *[N+1];//Pres代表压力
    for(int k=0;k<N+1;k++)
        Pres[k]= new double [M+1];
    for(int i=0;i<N+1;i++)
        for(int j=0;j<M+1;j++)
            Pres[i][j]=0;

    double **ustar=new double *[N+1];
    for(int k=0;k<N+1;k++)
        ustar[k]=new double [M+1];
    for(int i=0;i<N+1;i++)
        for(int j=0;j<M+1;j++)
            ustar[i][j]=0;
    cout<<endl;
    double **vstar=new double *[N+1];
    for(int k=0;k<N+1;k++)
        vstar[k]=new double [M+1];
    for(int i=0;i<N+1;i++)
        for(int j=0;j<M+1;j++)
            vstar[i][j]=0;
    cout<<endl;
    double **u=new double *[N+1];
    for(int k=0;k<N+1;k++)
        u[k]=new double [M+1];
    for(int i=0;i<N+1;i++)
        for(int j=0;j<M+1;j++)
            u[i][j]=1;
    double **v=new double *[N+1];
    for(int k=0;k<N+1;k++)
        v[k]=new double [M+1];
    for(int i=0;i<N+1;i++)
        for(int j=0;j<M+1;j++)
            v[i][j]=0;

    //cout<<endl;
    /*---------边值条件--------------*/
    for(int i=0;i<1;i++)
    {
        for(int j=0;j<=M;j++)
        {
            u[i][j]=U;
            v[i][j]=0.0;//初始的扰动项
        }
    }
    for(int i=N;i<N+1;i++)
    {
        for(int j=0;j<=M;j++)
        {
            u[i][j]=U;
            v[i][j]=0.0;
        }
    }
    for(int i=0;i<N+1;i++)
    {
        for(int j=0;j<1;j++)
        {
            u[i][j]=U;
            v[i][j]=0.0;
        }
    }
    for(int i=0;i<N+1;i++)
    {
        for(int j=M;j<M+1;j++)
        {
            u[i][j]=U;
            v[i][j]=0.0;
        }
    }
    /*------------------------------------------------------------------------------------------------------*/
    /*cout<<"-----Output the u ------"<<endl;
    for(int i=0;i<N+1;i++)
    {
        for(int j=0;j<M+1;j++)
        {
            cout<<" "<<u[i][j];
        }
        cout<<endl;
    }*/
    /*---------------------给v分配空间并赋初值---------------------------*/

    /*给u'开辟空间并赋初值*/
    double **upie=new double *[N+1];
    for(int k=0;k<N+1;k++)
        upie[k]= new double [M+1];
    for(int i=0;i<N+1;i++)
        for(int j=0;j<M+1;j++)
            upie[i][j]=U;

    /*给v'开辟空间并赋初值*/
    double **vpie=new double *[N+1];
    for(int k=0;k<N+1;k++)
        vpie[k]= new double [M+1];
    for(int i=0;i<N+1;i++)
        for(int j=0;j<M+1;j++)
            vpie[i][j]=0;
    /*---------边值条件---------u'------*/
    for(int i=0;i<1;i++)
    {
        for(int j=0;j<=M;j++)
        {
            upie[i][j]=U;
            vpie[i][j]=0;
        }
    }
    for(int i=N;i<N+1;i++)
    {
        for(int j=0;j<=M;j++)
        {
            upie[i][j]=U;
            vpie[i][j]=0;
        }
    }
    for(int i=0;i<N+1;i++)
    {
        for(int j=0;j<1;j++)
        {
            upie[i][j]=U;
            vpie[i][j]=0;
        }
    }
    for(int i=0;i<N+1;i++)
    {
        for(int j=M;j<M+1;j++)
        {
            upie[i][j]=U;
            vpie[i][j]=0;
        }
    }


    // cout<<"读入初始流场u[i][j],v[i][j],p[i][j]"<<endl;


/*==================================建议开始=======================================*/
    base_address = __base_address;
    infile.open(base_address.append("fluidfield-20000.txt").data(),ios::in);
    /*建议原因：打开的文件可能需要是不带序号的那个fluidfield*/
/*==================================建议结束=======================================*/

    if(!infile){

        cout<<"不能打开solid-node-cor-ini.txt数据文件!\n";
        exit(1);
    }

    for(int j=M;j>=0;j--)
    {
        for(int  i=0;i<N+1;i++)
        {
            infile>>u[i][j];
            infile>>v[i][j];
            infile>>Pres[i][j];
        }
    }
    infile.close();
    cout<<"读入初始流场u[i][j],v[i][j],p[i][j]结束"<<endl;


    /*double **eta=new double *[N];
    for(int k=0;k<N;k++)
        eta[k]= new double [M];
    for(int i=0;i<N;i++)
        for(int j=0;j<M;j++)
            eta[i][j]=10;//本来eta[i][j]在0~1之间，这里取10为了防止错误，便于检查
    /*----将Ax=b中A用Eigen中的数据格式保存起来，同时将UMFPACK中的数据接口的三个变量实现出来------*/
    int SMADi=(N+1)*(M+1);
    SpMat SparseMatA(SMADi,SMADi); //定义稀疏矩阵A--SparseMatA，维数为SMADi
    typedef Eigen::Triplet<double> Tri;
    vector<Tri> coefficients;
    //coefficients.push_back(Tri(0,0,2.0));
    for(int i=0;i<(M+1);i++)
        coefficients.push_back(Tri(i,i,1.0));
    for(int i=M+1;i<SMADi-(M+1);i++)
    {
        if(i%(M+1)==0) //B矩阵中的首行
        {
            coefficients.push_back(Tri(i,i,4.0));
            coefficients.push_back(Tri(i,i+1,-2.0));
            coefficients.push_back(Tri(i,i-(M+1),-1.0));
            coefficients.push_back(Tri(i,i+(M+1),-1.0));
        }
        else if((i+1)%(M+1)==0) //B矩阵中的末行
        {
            coefficients.push_back(Tri(i,i,4.0));
            coefficients.push_back(Tri(i,i-1,-2.0));
            coefficients.push_back(Tri(i,i-(M+1),-1.0));
            coefficients.push_back(Tri(i,i+(M+1),-1.0));
        }
        else //B矩阵的中间行
        {
            coefficients.push_back(Tri(i,i,4.0));
            coefficients.push_back(Tri(i,i-1,-1.0));
            coefficients.push_back(Tri(i,i+1,-1.0));
            coefficients.push_back(Tri(i,i-(M+1),-1.0));
            coefficients.push_back(Tri(i,i+(M+1),-1.0));
        }
    }
    for(int i=SMADi-(M+1);i<SMADi;i++)
    {
        if(i==SMADi-(M+1)) //B矩阵中的首行
        {
            coefficients.push_back(Tri(i,i,4.0));
            coefficients.push_back(Tri(i,i+1,-2.0));
            coefficients.push_back(Tri(i,i-(M+1),-2.0));
        }
        else if(i==SMADi-1) //B矩阵中的末行
        {
            coefficients.push_back(Tri(i,i,4.0));
            coefficients.push_back(Tri(i,i-1,-2.0));
            coefficients.push_back(Tri(i,i-(M+1),-2.0));
        }
        else //B矩阵的中间行
        {
            coefficients.push_back(Tri(i,i,4.0));
            coefficients.push_back(Tri(i,i-1,-1.0));
            coefficients.push_back(Tri(i,i+1,-1.0));
            coefficients.push_back(Tri(i,i-(M+1),-2.0));
        }
    }

    SparseMatA.setFromTriplets(coefficients.begin(),coefficients.end());
    cout<<endl;
    int *SparseMatAp=new int[SMADi+1];
    SparseMatAp[0]=0;
    int SMAnZ=SparseMatA.nonZeros();//稀疏矩阵非零元的个数
    int *SparseMatAi=new int[SMAnZ];
    double *SparseMatAx=new double[SMAnZ];
    SMAnZ=0;
    cout<<"------SparseMatAx-------"<<endl;
    cout<<"SparseMatA.outerSize()="<<SparseMatA.outerSize()<<endl;
    //TODO 我修改了这里
    for(int i=0;i<SparseMatA.outerSize();i++)
    {
        SparseMatAp[i+1]=SparseMatAp[i];
        for (Eigen::SparseMatrix<double>::InnerIterator it(SparseMatA,i); it; ++it)
        {
            SparseMatAx[SMAnZ]=it.value();

            SparseMatAi[SMAnZ]=it.row();
            //cout<<SparseMatAx[SMAnZ]<<" ";
            SMAnZ++;
            SparseMatAp[i+1]++;
        }
        //cout<<SparseMatAp[i]<<" ";
        //cout<<endl;
    }
    /*----------对稀疏矩阵(Ax=b中的A)操作结束-----------------*/
    //SolEta((double **)eta, N, M);//更新eta

    int ComputeStep=100000;//总共计算的步数


/*========================modify====修改开始====modify=============================*/
    int MemIntStep=100;//MemIntStep 代表Memory interval step 存储的间隔步数
    int forestep=100;
    /*修改原因：这两个参数原本用来控制保存文件的间隔，但是使用新的状态保存方法后，就没有必要一直保存了，所以可以把间隔调大*/
/*==================================修改结束=======================================*/



/*===========================add====添加开始====add================================*/
    int save_state_setp = 10; //代表状态存储的间隔步数
    /*添加原因：这个参数用来控制每隔多少步进行一次状态保存，定得太小会因为保存频繁而浪费时间和磁盘空间*/
/*==================================添加结束=======================================*/


/*===========================add====添加开始====add================================*/
    /*
     * 恢复状态开始*/
    int current_index;
    char recover_state;
    cout<<"是否从恢复状态？y/n"<<endl;
    cin>>recover_state;
    //从头计算
    if(recover_state == 'y'){
        current_index = load_state(ustar,vstar,u,v,upie,vpie,Pres,
                                   SparseMatAi,SparseMatAx,SparseMatAp,
                                   X2,scc,svel,inifcc_2,
                                   &N,&M,&SMAnZ,&SMADi,&nii_2,&njj_2,&nkk,
                                   &delta_t,&delta_x,&delta_y,
                                   &CF,&Re,&Turb,&NN_2,&MM_2,&N_sec_2,&M_sec_2,&ss_2);
        //如果恢复状态不成功
        if(current_index == -1){
            char start_from_begin;
            cout<<"是否从头计算？y/n"<<endl;
            cin>>start_from_begin;
            //从头计算
            if(start_from_begin == 'y'){
                current_index = -1;
            }else{//不从头计算
                exit(-1);
            }
        }
    }else{
        current_index = -1;
    }

    /*
     * 恢复状态结束*/
    /*添加原因：这部分的作用是在每次重新开始运行前，恢复之前的状态*/
/*==================================添加结束=======================================*/


/*========================modify====修改开始====modify=============================*/
    for(int k=current_index+1;k<=ComputeStep;k++)
    {
        /*修改原因：把k=0修改为k=current_index+1,目的是从恢复后的序号开始运行，而不是从头开始运行，current_index就是状态目录中最新的state文件的序号*/
/*==================================修改结束=======================================*/


        currenttime =k*delta_t;
        if(k<=20000)
        {
            for (i = 0; i<nii_2; i++)
            {
                scc[i][0] = inisc_2[i][0];
                scc[i][1] = inisc_2[i][1];
                svel[i][0] = 0.0;
                svel[i][1] = 0.0;
            }
            cout<<"---the "<<k<<"   step---"<<endl;
            Sol_ustar = CUEqution; //将函数指针指向函数CU;
            (*Sol_ustar)((double **)ustar, (double **)u, (double **)v, CF, delta_x, delta_y, delta_t, Re, Turb, N, M); //调用CU函数

            Sol_vstar = CVEqution; //将函数指针指向函数CV;
            (*Sol_vstar)((double **)vstar, (double **)u, (double **)v, CF, delta_x, delta_y, delta_t, Re, Turb, N, M); //调用CV函数

            Sol_P(N, M, delta_x, delta_y, delta_t, (double **)ustar, (double **)vstar, (double **)upie, (double **)vpie, (double **)Pres, SparseMatAi, SparseMatAx, SparseMatAp); //调用CV函数

            Sol_upie = SolUPie; //将函数指针指向函数SolUPie;
            (*Sol_upie)((double **)ustar, (double **)upie, (double **)Pres, N, M); //调用CU函数

            Sol_vpie = SolVPie; //将函数指针指向函数SolVPie;
            (*Sol_vpie)((double **)vstar, (double **)vpie, (double **)Pres, N, M); //调用CV函数

            //Sol_VelPlusOne = VelPlusOne; //将函数指针指向函数Sol-Vel在下一步的值;
            //(*Sol_VelPlusOne)((double **)u, (double **)X1, (double **)upie, (double **)v, (double **)vpie, (double **)iniscc_1, (double **)inisvell_1, (double **)inifcc_1, N, M, NN_1, MM_1, N_sec_1, M_sec_1, delta_x, delta_y, ss_1, delta_t, nii_1, njj_1);

            Sol_VelPlusTwo = VelPlusTwo; //将函数指针指向函数Sol-Vel在下一步的值;
            (*Sol_VelPlusTwo)((double **)u, (double **)X2, (double **)v, (double **)scc, (double **)svel, (double **)inifcc_2, N, M, NN_2, MM_2, N_sec_2, M_sec_2, delta_x, delta_y, ss_2, delta_t, nii_2, njj_2);

        }

        else
        {

            for(int i=0;i<NN_2+1;i++)
            {
                for(int j=0;j<MM_2+1;j++)
                {
                    temp1=-0.7+int(dis[0][1]/delta_x)*delta_x+i*delta_x;
                    temp2=-0.7+int(dis[1][1]/delta_y)*delta_y+j*delta_y;
                    n=i*(MM_2+1)+j;
                    inifcc_2[n][0]=temp1;
                    inifcc_2[n][1]=temp2;
                }
            }

            cout<<"---the "<<k<<"   step---"<<endl;
            Sol_ustar=CUEqution; //将函数指针指向函数CU;
            (*Sol_ustar)((double **)ustar,(double **)u,(double **)v,CF,delta_x,delta_y,delta_t,Re,Turb,N,M); //调用CU函数

            Sol_vstar=CVEqution; //将函数指针指向函数CV;
            (*Sol_vstar)((double **)vstar,(double **)u,(double **)v,CF,delta_x,delta_y,delta_t,Re,Turb,N,M); //调用CV函数

            Sol_P(N,M,delta_x,delta_y,delta_t,(double **) ustar,(double **) vstar,(double **) upie,(double **) vpie,(double **) Pres,SparseMatAi,SparseMatAx,SparseMatAp); //调用CV函数

            Sol_upie=SolUPie; //将函数指针指向函数SolUPie;
            (*Sol_upie)((double **)ustar,(double **)upie,(double **)Pres, N, M ); //调用CU函数

            Sol_vpie=SolVPie; //将函数指针指向函数SolVPie;
            (*Sol_vpie)((double **)vstar,(double **)vpie,(double **)Pres, N, M ); //调用CV函数

            //Sol_VelPlusOne=VelPlusOne; //将函数指针指向函数Sol-Vel在下一步的值;
            //(*Sol_VelPlusOne)((double **)u,(double **)X1,(double **)upie,(double **)v,(double **)vpie,(double **)iniscc_1,(double **)inisvell_1,(double **)inifcc_1, N, M,NN_1,MM_1,N_sec_1,M_sec_1,delta_x,delta_y,ss_1,delta_t,nii_1,njj_1);


            Sol_VelPlusTwo=VelPlusTwo; //将函数指针指向函数Sol-Vel在下一步的值;
            (*Sol_VelPlusTwo)((double **)u,(double **)X2,(double **)v,(double **)scc,(double **)svel,(double **)inifcc_2, N, M,NN_2,MM_2,N_sec_2,M_sec_2,delta_x,delta_y,ss_2,delta_t,nii_2,njj_2);

            Sol_rigid=solrigid;
            (*Sol_rigid)((double **)X2,(double **)scc,(double **)svel,(double **)dis,(double **)vel,(double **)inisc_2,mass,damp,rigid,rot_iner,rot_damp,rot_rigid,ss_2,delta_t,delta_x,nii_2,k);
        }


/*===========================add====添加开始====add================================*/
        /*保存状态开始*/
        if(k% save_state_setp == 0){
            save_state(k,ustar,vstar,u,v,upie,vpie,Pres,
                       SparseMatAi,SparseMatAx,SparseMatAp,
                       X2,scc,svel,inifcc_2,
                       N,M,SMAnZ,SMADi,nii_2,njj_2,nkk,delta_t,delta_x,delta_y,CF,Re,Turb,NN_2,MM_2,N_sec_2,M_sec_2,ss_2);
        }
        /*保存状态结束*/
/*添加原因：这部分用来保存状态*/
/*==================================添加结束=======================================*/


        if(k%MemIntStep==0)
        {
            //printf("%d.txt",i);
            char buff[256];
            base_address = __base_address;
            sprintf(buff, base_address.append("results/vel-%d.plt").data(), k); //每隔MemIntStep步，写入到文本中，文件名做相应改变，不能覆盖
            ofstream outVelandPre;
            outVelandPre.open(buff, ios::out);
            outVelandPre<<"Title=\"meth"<<N<<"and"<<ComputeStep<<"step\""<<endl;
            outVelandPre<<"Variables=\"x\"\"y\"\"u\"\"v\"\"p\""<<endl;
            outVelandPre<<"Zone I="<<N+1<<" j="<<M+1<<" f=point"<<endl;
            for(int j=M;j>=0;j--)
            {
                for(int  i=0;i<N+1;i++)
                {

                    outVelandPre<<-5+i*delta_x<<" "<<-5+j*delta_y<<" "<<u[i][j]<<" "<<v[i][j]<<" "<<Pres[i][j];
                    outVelandPre<<endl;
                }

                //outVelandPre<<endl;
            }
            outVelandPre.close();
        }
        cout<<"输出流场信息完成"<<endl;

        if(k%forestep==0)
        {
            //printf("%d.txt",i);
            solidforcex=0.0;
            solidforcey=0.0;
            char bufff[256];
            base_address = __base_address;
            sprintf(bufff, base_address.append("results/liftdrag2-%d.txt").data(), k); //每隔MemIntStep步，写入到文本中，文件名做相应改变，不能覆盖
            ofstream outtforce;
            outtforce.open(bufff, ios::out);
            for (i=0; i<nii_2; i++)
            {
                solidforcex=solidforcex-2*X2[i][0]*ss_2*delta_x;
                solidforcey=solidforcey-2*X2[i][1]*ss_2*delta_x;
            }
            outtforce<<solidforcex<<"   "<<solidforcey;
            outtforce <<endl;
            outtforce.close ();
            cout<<"输出升力、阻力完成"<<endl;
        }

        if(k%forestep==0)
        {
            //printf("%d.txt",i);
            double solidforcexx=0.0;
            double solidforceyy=0.0;
            char bufff[256];
            base_address = __base_address;
            sprintf(bufff, base_address.append("results/solidforcee-%d.txt").data(), k); //每隔MemIntStep步，写入到文本中，文件名做相应改变，不能覆盖
            ofstream outtforce;
            outtforce.open(bufff, ios::out);
            for (i=0; i<nii_2; i++)
            {
                solidforcexx=X2[i][0];
                solidforceyy=X2[i][1];
                solidcorx=inisc_2[i][0];
                solidcory=inisc_2[i][1];
                xx=solidcory/solidcorx;
                if (solidcorx>0)
                {
                    aef=atan(xx);
                }
                else
                {
                    aef=pi+atan(xx);
                }
                ft=delta_x*(X2[i][0]*sin(aef)-X2[i][1]*cos(aef));
                fn=delta_x*(X2[i][0]*cos(aef)+X2[i][1]*sin(aef));
                outtforce<<solidcorx<<"   "<<solidcory<<"   "<<solidforcexx<<"   "<<solidforceyy<<"   "<<ft<<"   "<<fn;
                outtforce <<endl;
            }
            outtforce.close ();
            cout<<"输出固体节点力完成"<<endl;
        }


        /*if(k%forestep==0&k>=0)
        {
            //printf("%d.txt",i);
            double solidforcexx=0.0;
            double solidforceyy=0.0;
            char bufff[256];
            sprintf(bufff, "f:\\wwq\\programm-4\\cylinder\\Re200-20\\move-%d.txt", k); //每隔MemIntStep步，写入到文本中，文件名做相应改变，不能覆盖
            ofstream outtforce;
            outtforce.open(bufff, ios::out);
           for (i=0; i<nii; i++)
           {
           // solidforcexx=X[i][0];
            //solidforceyy=X[i][1];
            //solidcorx=scc[i][0];
            //solidcory=scc[i][1];
               //outtforce<<solidcorx<<"   "<<solidcory<<"   "<<solidforcexx <<"   "<<solidforceyy<<"   "<<svel[i][0]<<"   "<<svel[i][1];
            outtforce<<scc[i][0]<<"   "<<scc[i][1]<<"   "<<X[i][0] <<"   "<<X[i][1]<<"   "<<svel[i][0]<<"   "<<svel[i][1];
            outtforce <<endl;
            }
            outtforce.close ();
            cout<<"输出固体节点力完成"<<endl;
         }	*/
    }

    //for(int t=0;t<N;t++)
    //delete [] eta[t];
    //delete [] eta;

    for(int t=0;t<N+1;t++)
        delete [] Pres[t];
    delete [] Pres;

    for(int t=0;t<N+1;t++)
        delete []u[t];
    delete []u;
    for(int t=0;t<N+1;t++)
        delete []v[t];
    delete []v;
    for(int t=0;t<N+1;t++)
        delete []ustar[t];
    delete []ustar;
    for(int t=0;t<N+1;t++)
        delete []vstar[t];
    delete []vstar;

    /*撤销u'空间*/
    for(int t=0;t<N+1;t++)
        delete [] upie[t];
    delete [] upie;

    /*撤销v'空间*/
    for(int t=0;t<N+1;t++)
        delete [] vpie[t];
    delete [] vpie;

    return 0;
}
