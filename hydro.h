#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

using namespace std;

// trim from start
static inline string &ltrim(string &s) {
 s.erase(s.begin(), std::find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
 return s;
}

// trim from end
static inline std::string &rtrim(string &s) {
 s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
 return s;
}

// trim from both ends
static inline string &trim(string &s) {
 return ltrim(rtrim(s));
}
struct Hydro
{
 private:
  static const int maxx=120;
  static const int maxy=60;
  static const int maxh=20;
  static const double tau0=0.6;
  static const double tau1=24.6;
  static const double deltat=0.3;
  static const double deltax=0.3;
  static const double deltay=0.3;
  static const double deltah=0.3;
  int quench_method;
  double power;
  double power_t;
  double bmin_bin;
  double bmax_bin;
  double max_ncoll;
 public:
  double hydro_T[101][121][61][21];
  double hydro_E[101][121][61][21];
  double hydro_vx[101][121][61][21];
  double hydro_vy[101][121][61][21];
  double hydro_vh[101][121][61][21];
  
  static const int  nx=500;
  static const int  ny=500;
  static const int  nz=500;
  static const int  nb=100;
  static const int  nrho2=1000;
  
  static const double bmax=13.05;
  static const double xmax=10;
  static const double ymax=10;
  static const double zmax=10;
  static const double rho2max=150;
  double Ncoll[nx][ny][nb];
  double TAA[nrho2];
  
  Hydro(){};
  Hydro(bool is_txt, int centrality_bin, const char *file_Ncoll, const char *file_TAA)
  {
  
   std::string min_centrality_bins[]={"00","05","10","20","30","40","50","60","70"};
   std::string max_centrality_bins[]={"05","10","20","30","40","50","60","70","80"};

   std::stringstream fhydro;
   fhydro <<"../HiranoHydro/PbPb2760_"<<min_centrality_bins[centrality_bin]<<"-"<<max_centrality_bins[centrality_bin];   
   
   if(is_txt){  
    fhydro<<".dat";
    cout<<"start"<<endl; 
    ifstream input_file(fhydro.str().c_str());
    string line;
    int state=1;
    quench_method=0;
    power=1;
    power=0;
    while(getline(input_file,line))
    {
     trim(line);
     if(line.length()==0) continue;
     if(state==1)
     {
      double tau;
      double x;
      double y;
      double eta;
      double E;
      double T;
      double vx;
      double vy;
      double veta;
      double phase;
      stringstream line_s;
      line_s << line;
      line_s >> tau >> x >> y >> eta >> E >> T >> phase >> vx >> vy >> veta;
  	 if(tau==-999)break;
      int it = int((tau+deltat/2.-tau0)/deltat);
      int ix = int((x+deltax*maxx/2.+deltax/2.)/deltax);
      int iy = int((y+deltay/2.)/deltay);
      int ih = int((eta+deltah/2.)/deltah);
      hydro_T[it][ix][iy][ih]=T;
      hydro_E[it][ix][iy][ih]=E;
      hydro_vx[it][ix][iy][ih]=vx;
      hydro_vy[it][ix][iy][ih]=vy;
      hydro_vh[it][ix][iy][ih]=veta;
     }
     else
     {
      cout <<"Unexpected else. Exiting..." << endl;
      exit(1);
     }
    }
    input_file.close();
    cout<<"end"<<endl;
   }else{
    fhydro<<".bin";
    ifstream myFile;
    myFile.open(fhydro.str().c_str(), ios::in | ios::binary);
 
    myFile.read((char *)(hydro_T), sizeof(hydro_T));
    myFile.read((char *)(hydro_E), sizeof(hydro_E));
    myFile.read((char *)(hydro_vx), sizeof(hydro_vx));
    myFile.read((char *)(hydro_vy), sizeof(hydro_vy));
    myFile.read((char *)(hydro_vh), sizeof(hydro_vh));

    myFile.close();   
   }
   
   double bmin_values[]={0.,3.5,4.94,6.98,8.55,9.88,11.04,12.09};
   double bmax_values[]={3.5,4.94,6.98,8.55,9.88,11.04,12.09,13.05};
  
   this->set_bmin_bin(bmin_values[centrality_bin]);
   this->set_bmax_bin(bmax_values[centrality_bin]);
   
   ifstream fileNcoll;
   fileNcoll.open(file_Ncoll, ios::in | ios::binary);
   fileNcoll.read((char *)(Ncoll), sizeof(Ncoll));
   
   ifstream fileTAA;
   fileTAA.open(file_TAA, ios::in | ios::binary);
   fileTAA.read((char *)(TAA), sizeof(TAA));
   
   max_ncoll=0;
   double max_taa=0;
   for(int i=0;i<nrho2;i++){
    if(TAA[i]>max_taa) max_taa=TAA[i];
   }
   
   for(int i=0;i<nb;i++){
    for(int j=0;j<nx;j++){
	 for(int k=0;k<ny;k++){
	  if(Ncoll[j][k][i]>max_ncoll) max_ncoll=Ncoll[j][k][i];
	 }
	}
   }
   
   if((max_taa*max_taa)>max_ncoll) max_ncoll=max_taa*max_taa;
   set_max_ncoll(max_ncoll);
   
  }
  double get_maxx(){return maxx;}
  double get_maxy(){return maxy;}
  double get_maxh(){return maxh;}
  
  void set_bmin_bin(double bmin){
   this->bmin_bin=bmin;
  }
  
  void set_bmax_bin(double bmax){
   this->bmax_bin=bmax;
  }
  
  double get_bmin_bin(){
   return bmin_bin;
  }
  
  double get_bmax_bin(){
   return bmax_bin;
  }
  
  void set_max_ncoll(double max_ncoll){
   this->max_ncoll=max_ncoll;
  }
  
  double get_max_ncoll(){
   return max_ncoll;
  }
  
  int get_it(double tau){
   int it=int((tau-tau0)/deltat);
   return it;
  }
 
  double get_dt(int it,double tau){
   double dt=(tau-tau0-double(it)*deltat)/deltat;
   return dt;
  }
 
  int get_ix(double x){
   int ix;
   if(x>0.) ix=((int)(x/deltax)+maxx/2.);
   else ix=((int)(-x/deltax)+maxx/2.-1);
   return ix;
  }
 
  double get_dx(int ix,double x){
   double dx=(fabs(x)-double(ix-maxx/2.)*deltax)/deltax;
   return dx;
  }
 
  int get_iy(double y){
   int iy=(int)(fabs(y)/deltay);
   return iy;
  }
 
  double get_dy(int iy,double y){
   double dy=(fabs(y)-double(iy)*deltay)/deltay;
   return dy;
  }
 
  int get_ih(double eta){
   int ih=int(fabs(eta/deltah));
   return ih;
  }
 
  double get_dh(int ih,double eta){
   double dh=(fabs(eta)-double(ih)*deltah)/deltah;
   return dh;  
  }
  
  double get_hydro_E(double tau, double x, double y, double eta){ 
   int it=get_it(tau);
   int ix=get_ix(x);
   int iy=get_iy(y);
   int ih=get_ih(eta);
   double dt=get_dt(it,tau);
   double dx=get_dx(ix,x);
   double dy=get_dy(iy,y);
   double dh=get_dh(ih,eta);
   double E=0;
   
   if (tau>tau1) {
    return E;
   }
  
   if (ix<0 || ix>maxx || iy>maxy || ih>maxh) return E;
   
   E=hydro_E[it][ix][iy][ih]*(1.-dt)*(1.-dx)*(1.-dy)*(1.-dh);
   E+=hydro_E[it+1][ix][iy][ih]*(1.-dx)*(1.-dy)*(1.-dh)*dt;
   E+=hydro_E[it][ix+1][iy][ih]*(1.-dt)*(1.-dy)*(1.-dh)*dx;
   E+=hydro_E[it][ix][iy+1][ih]*(1.-dx)*(1.-dt)*(1.-dh)*dy;
   E+=hydro_E[it][ix][iy][ih+1]*(1.-dx)*(1.-dy)*(1.-dt)*dh;
   E+=hydro_E[it+1][ix+1][iy][ih]*dx*(1.-dy)*(1.-dh)*dt;
   E+=hydro_E[it][ix+1][iy+1][ih]*(1.-dt)*dy*(1.-dh)*dx;
   E+=hydro_E[it][ix][iy+1][ih+1]*(1.-dx)*(1.-dt)*dh*dy;
   E+=hydro_E[it][ix+1][iy][ih+1]*dx*(1.-dy)*(1.-dt)*dh;
   E+=hydro_E[it+1][ix][iy+1][ih]*(1.-dx)*dt*(1.-dh)*dy;
   E+=hydro_E[it+1][ix][iy][ih+1]*(1.-dx)*(1.-dy)*dt*dh;
   E+=hydro_E[it][ix+1][iy+1][ih+1]*dx*(1.-dt)*dh*dy;
   E+=hydro_E[it+1][ix+1][iy][ih+1]*dx*(1.-dy)*dt*dh;
   E+=hydro_E[it+1][ix+1][iy+1][ih]*dx*dt*(1.-dh)*dy;
   E+=hydro_E[it+1][ix][iy+1][ih+1]*(1.-dx)*dy*dt*dh;
   E+=hydro_E[it+1][ix+1][iy+1][ih+1]*dx*dy*dt*dh;
  
   return E;  
  }
 
  double get_hydro_T(double tau, double x, double y, double eta){ 
   int it=get_it(tau);
   int ix=get_ix(x);
   int iy=get_iy(y);
   int ih=get_ih(eta);
   double dt=get_dt(it,tau);
   double dx=get_dx(ix,x);
   double dy=get_dy(iy,y);
   double dh=get_dh(ih,eta);
   double T=0;
  
   if (tau>tau1) {
    return T;
   }
   if (ix<0 || ix>maxx || iy>maxy || ih>maxh) return T;
  
   T=hydro_T[it][ix][iy][ih]*(1.-dt)*(1.-dx)*(1.-dy)*(1.-dh);
   T+=hydro_T[it+1][ix][iy][ih]*(1.-dx)*(1.-dy)*(1.-dh)*dt;
   T+=hydro_T[it][ix+1][iy][ih]*(1.-dt)*(1.-dy)*(1.-dh)*dx;
   T+=hydro_T[it][ix][iy+1][ih]*(1.-dx)*(1.-dt)*(1.-dh)*dy;
   T+=hydro_T[it][ix][iy][ih+1]*(1.-dx)*(1.-dy)*(1.-dt)*dh;
   T+=hydro_T[it+1][ix+1][iy][ih]*dx*(1.-dy)*(1.-dh)*dt;
   T+=hydro_T[it][ix+1][iy+1][ih]*(1.-dt)*dy*(1.-dh)*dx;
   T+=hydro_T[it][ix][iy+1][ih+1]*(1.-dx)*(1.-dt)*dh*dy;
   T+=hydro_T[it][ix+1][iy][ih+1]*dx*(1.-dy)*(1.-dt)*dh;
   T+=hydro_T[it+1][ix][iy+1][ih]*(1.-dx)*dt*(1.-dh)*dy;
   T+=hydro_T[it+1][ix][iy][ih+1]*(1.-dx)*(1.-dy)*dt*dh;
   T+=hydro_T[it][ix+1][iy+1][ih+1]*dx*(1.-dt)*dh*dy;
   T+=hydro_T[it+1][ix+1][iy][ih+1]*dx*(1.-dy)*dt*dh;
   T+=hydro_T[it+1][ix+1][iy+1][ih]*dx*dt*(1.-dh)*dy;
   T+=hydro_T[it+1][ix][iy+1][ih+1]*(1.-dx)*dy*dt*dh;
   T+=hydro_T[it+1][ix+1][iy+1][ih+1]*dx*dy*dt*dh;
   
   return T/1000.; 
  }
 
  double get_hydro_vx(double tau, double x, double y, double eta){ 
   int it=get_it(tau);
   int ix=get_ix(x);
   int iy=get_iy(y);
   int ih=get_ih(eta);
   double dt=get_dt(it,tau);
   double dx=get_dx(ix,x);
   double dy=get_dy(iy,y);
   double dh=get_dh(ih,eta);
   double vx=0;
   
   if (tau>tau1) {
    return vx;
   }
   
   if (ix<0 || ix>maxx || iy>maxy || ih>maxh) return vx;
  
   vx=hydro_vx[it][ix][iy][ih]*(1.-dt)*(1.-dx)*(1.-dy)*(1.-dh);
   vx+=hydro_vx[it+1][ix][iy][ih]*(1.-dx)*(1.-dy)*(1.-dh)*dt;
   vx+=hydro_vx[it][ix+1][iy][ih]*(1.-dt)*(1.-dy)*(1.-dh)*dx;
   vx+=hydro_vx[it][ix][iy+1][ih]*(1.-dx)*(1.-dt)*(1.-dh)*dy;
   vx+=hydro_vx[it][ix][iy][ih+1]*(1.-dx)*(1.-dy)*(1.-dt)*dh;
   vx+=hydro_vx[it+1][ix+1][iy][ih]*dx*(1.-dy)*(1.-dh)*dt;
   vx+=hydro_vx[it][ix+1][iy+1][ih]*(1.-dt)*dy*(1.-dh)*dx;
   vx+=hydro_vx[it][ix][iy+1][ih+1]*(1.-dx)*(1.-dt)*dh*dy;
   vx+=hydro_vx[it][ix+1][iy][ih+1]*dx*(1.-dy)*(1.-dt)*dh;
   vx+=hydro_vx[it+1][ix][iy+1][ih]*(1.-dx)*dt*(1.-dh)*dy;
   vx+=hydro_vx[it+1][ix][iy][ih+1]*(1.-dx)*(1.-dy)*dt*dh;
   vx+=hydro_vx[it][ix+1][iy+1][ih+1]*dx*(1.-dt)*dh*dy;
   vx+=hydro_vx[it+1][ix+1][iy][ih+1]*dx*(1.-dy)*dt*dh;
   vx+=hydro_vx[it+1][ix+1][iy+1][ih]*dx*dt*(1.-dh)*dy;
   vx+=hydro_vx[it+1][ix][iy+1][ih+1]*(1.-dx)*dy*dt*dh;
   vx+=hydro_vx[it+1][ix+1][iy+1][ih+1]*dx*dy*dt*dh;
   
   if(eta<=0) vx=-vx;
   return vx; 
  }
  
  double get_hydro_vy(double tau, double x, double y, double eta){ 
   int it=get_it(tau);
   int ix=get_ix(x);
   int iy=get_iy(y);
   int ih=get_ih(eta);
   double dt=get_dt(it,tau);
   double dx=get_dx(ix,x);
   double dy=get_dy(iy,y);
   double dh=get_dh(ih,eta);
   double vy=0;
   
   if (tau>tau1) {
    return vy;
   }
   
   if (ix<0 || ix>maxx || iy>maxy || ih>maxh) return vy;
   
   vy=hydro_vy[it][ix][iy][ih]*(1.-dt)*(1.-dx)*(1.-dy)*(1.-dh);
   vy+=hydro_vy[it+1][ix][iy][ih]*(1.-dx)*(1.-dy)*(1.-dh)*dt;
   vy+=hydro_vy[it][ix+1][iy][ih]*(1.-dt)*(1.-dy)*(1.-dh)*dx;
   vy+=hydro_vy[it][ix][iy+1][ih]*(1.-dx)*(1.-dt)*(1.-dh)*dy;
   vy+=hydro_vy[it][ix][iy][ih+1]*(1.-dx)*(1.-dy)*(1.-dt)*dh;
   vy+=hydro_vy[it+1][ix+1][iy][ih]*dx*(1.-dy)*(1.-dh)*dt;
   vy+=hydro_vy[it][ix+1][iy+1][ih]*(1.-dt)*dy*(1.-dh)*dx;
   vy+=hydro_vy[it][ix][iy+1][ih+1]*(1.-dx)*(1.-dt)*dh*dy;
   vy+=hydro_vy[it][ix+1][iy][ih+1]*dx*(1.-dy)*(1.-dt)*dh;
   vy+=hydro_vy[it+1][ix][iy+1][ih]*(1.-dx)*dt*(1.-dh)*dy;
   vy+=hydro_vy[it+1][ix][iy][ih+1]*(1.-dx)*(1.-dy)*dt*dh;
   vy+=hydro_vy[it][ix+1][iy+1][ih+1]*dx*(1.-dt)*dh*dy;
   vy+=hydro_vy[it+1][ix+1][iy][ih+1]*dx*(1.-dy)*dt*dh;
   vy+=hydro_vy[it+1][ix+1][iy+1][ih]*dx*dt*(1.-dh)*dy;
   vy+=hydro_vy[it+1][ix][iy+1][ih+1]*(1.-dx)*dy*dt*dh;
   vy+=hydro_vy[it+1][ix+1][iy+1][ih+1]*dx*dy*dt*dh;
   
   if(y<=0) vy=-vy;
   return vy; 
 } 
  
  double get_hydro_vh(double tau, double x, double y, double eta){ 
   int it=get_it(tau);
   int ix=get_ix(x);
   int iy=get_iy(y);
   int ih=get_ih(eta);
   double dt=get_dt(it,tau);
   double dx=get_dx(ix,x);
   double dy=get_dy(iy,y);
   double dh=get_dh(ih,eta);
   double vh=0;
   
   if (tau>tau1) {
    return vh;
   }
   
   if (ix<0 || ix>maxx || iy>maxy || ih>maxh) return vh;
   
   vh=hydro_vh[it][ix][iy][ih]*(1.-dt)*(1.-dx)*(1.-dy)*(1.-dh);
   vh+=hydro_vh[it+1][ix][iy][ih]*(1.-dx)*(1.-dy)*(1.-dh)*dt;
   vh+=hydro_vh[it][ix+1][iy][ih]*(1.-dt)*(1.-dy)*(1.-dh)*dx;
   vh+=hydro_vh[it][ix][iy+1][ih]*(1.-dx)*(1.-dt)*(1.-dh)*dy;
   vh+=hydro_vh[it][ix][iy][ih+1]*(1.-dx)*(1.-dy)*(1.-dt)*dh;
   vh+=hydro_vh[it+1][ix+1][iy][ih]*dx*(1.-dy)*(1.-dh)*dt;
   vh+=hydro_vh[it][ix+1][iy+1][ih]*(1.-dt)*dy*(1.-dh)*dx;
   vh+=hydro_vh[it][ix][iy+1][ih+1]*(1.-dx)*(1.-dt)*dh*dy;
   vh+=hydro_vh[it][ix+1][iy][ih+1]*dx*(1.-dy)*(1.-dt)*dh;
   vh+=hydro_vh[it+1][ix][iy+1][ih]*(1.-dx)*dt*(1.-dh)*dy;
   vh+=hydro_vh[it+1][ix][iy][ih+1]*(1.-dx)*(1.-dy)*dt*dh;
   vh+=hydro_vh[it][ix+1][iy+1][ih+1]*dx*(1.-dt)*dh*dy;
   vh+=hydro_vh[it+1][ix+1][iy][ih+1]*dx*(1.-dy)*dt*dh;
   vh+=hydro_vh[it+1][ix+1][iy+1][ih]*dx*dt*(1.-dh)*dy;
   vh+=hydro_vh[it+1][ix][iy+1][ih+1]*(1.-dx)*dy*dt*dh;
   vh+=hydro_vh[it+1][ix+1][iy+1][ih+1]*dx*dy*dt*dh;
   
   if(eta<=0) vh=-vh;
   return vh; 
  }
  void set_quench_method(int quench_method){
    //0=collisional,1=radiative,2=light
    this->quench_method=quench_method;
    if(quench_method==0) power=2.;
    else if(quench_method==1){
     power=3.;
     power_t=1.;
    }
   }
   int get_quench_method(){
    return quench_method;
   }
   double get_power(){return power;}
   double get_power_t(){return power_t;}
   
};