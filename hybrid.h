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
// #include "Glauber/IntegrateT_v4.C"
// #include "Glauber/LQNumericIntegral.C"
// #include "Glauber/LQNumericIntegral_v2.C"
// #include "Glauber/jetPosition.C"
// #include "Glauber/fTnucleon.h"

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

struct Vector3
{
 private:
  double x;
  double y;
  double z;
 public:
  void reset(double x, double y, double z){
   this->x=x;
   this->y=y;
   this->z=z;
  }
  Vector3(double x, double y, double z){
   reset(x,y,z);
  }
  double get_x(){return x;}
  double get_y(){return y;}
  double get_z(){return z;}
  double dot(Vector3 *v2){
   return x*v2->get_x()+y*v2->get_y()+z*v2->get_z();
  }
};

struct Hydro
{
 private:
 double hydro_T[101][121][61][21];
 double hydro_E[101][121][61][21];
 double hydro_vx[101][121][61][21];
 double hydro_vy[101][121][61][21];
 double hydro_vh[101][121][61][21];
 
 static const int maxx=120;
 static const int maxy=60;
 static const int maxh=20;
 static const double tau0=0.6;
 static const double tau1=24.6;
 static const double deltat=0.3;
 static const double deltax=0.3;
 static const double deltay=0.3;
 static const double deltah=0.3;
 
 public:
 Hydro(const char *file_path)
 { 
  cout<<"start"<<endl;
  int ipoint=0;
  ifstream input_file(file_path);
  string line;
  int state=1;
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
 }
 
 double get_maxx(){return maxx;}
 double get_maxy(){return maxy;}
 double get_maxh(){return maxh;}
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
  
  return E/1000;  
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
  // cout << "it= "<<it<<" ix="<<ix<<" iy= "<<iy<<" ih="<<ih <<" dx="<<dx<<" dy= "<<dy<<" dh="<<dh<<" hydro_E[it][ix][iy][ih]="<<hydro_E[it][ix][iy][ih]<<endl;
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
  
  return vx/1000; 
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
  
  return vy/1000; 
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
  
  return vh/1000; 
 }
};


//----------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

struct Fragment 
{  
 private:
  int id;
  double px;
  double py;
  double pz;
  double E;
  double m;
  int QG;
  int initialIQG;
  // int dinfo;
  int momID;
  Fragment *mother;
  int number_of_daughters;
  Fragment **daughters;
  int number_of_ancestors;
  Fragment **ancestors;
  bool bool_incone;
  bool incone_ancestor;
  double taus;
  double initial_x;
  double initial_y;
  double initial_z;
  double final_x;
  double final_y;
  double final_z;
  double phi;
  double f;
  double quenched_E;
  double path_lenght;
  double T;
  double tstop;
  int jet_index_incone;
  int jet_index_incone_ancestor;
  int QG_jet;
  int quench_method;
  double power;
  double power_t;
 public:
  Fragment(int id, double px, double py, double pz, double E, double m, int momID, int QG, int initialIQG)
  {
   this->id = id;
   this->px = px;
   this->py = py;
   this->pz = pz;
   this->E = E;
   this->m = m;
   this->QG = QG;
   this->initialIQG = initialIQG;
   this->momID = momID;
   bool_incone = false;
   incone_ancestor=false;
   daughters = NULL;
   ancestors = NULL;
   mother = NULL;
   number_of_daughters = 0;
   number_of_ancestors = 0;
   taus = 0;
   initial_x=0;
   initial_y=0;
   initial_z=0;
   final_x=0;
   final_y=0;
   final_z=0;
   phi=0;
   f=0;
   quenched_E = E;
   path_lenght=0.5;
   T = 0.3;
   tstop = 0;
   jet_index_incone=-1;
   jet_index_incone_ancestor=-1;
   QG_jet=0;
   int quench_method=0;
   double power=1.;
   double power_t=0.;
  } 
//!PARTICLE PROPERTIES
  int get_id(){return id;}
  double get_px(){return px;}
  double get_py(){return py;}
  double get_pz(){return pz;}
  void set_px(double px){this->px=px;}
  void set_py(double py){this->py=py;}
  double get_E(){return E;}
  double get_m(){return m;}
  int get_QG(){return QG;}
  double get_initialIQG(){return initialIQG;}
  double get_tauf(){
   if(m==0) return 10000;
   else return 0.2*2*E/pow(m,2);
  }
  double get_tauf_geometry(){return get_gammaL()*0.2*2*E/pow(m,2);}
  double get_gammaL(){return E/sqrt(pow(px,2)+pow(py,2));}
  void set_taus(){
   if(mother==NULL) taus=get_tauf();
   if(mother!=NULL) taus = mother->get_taus()+get_tauf();
  }

  void set_taus_geometry(){
   if(id==3){
	   taus=0;
   } 
   if(mother!=NULL){
     taus = mother->get_taus()+get_tauf();   
   }
  }
  double get_taus(){return taus;} 
  int get_QG_jet(){
   if(id!=3 && mother->get_id()==3) QG_jet=get_QG();
   if(id!=3 && mother->get_id()!=3) QG_jet=mother->get_QG_jet();
   return QG_jet;
  }
  void set_initial_coordinates(){
   if(id>3 && mother->get_id()==3){
    final_x=initial_x+get_tauf()*px/E;
    final_y=initial_y+get_tauf()*py/E;
    final_z=initial_z+get_tauf()*pz/E;
    
    phi=atan2(get_py(),get_px());
   }
   else if(mother->get_id()!=3){
    initial_x = mother->get_final_x();
    initial_y = mother->get_final_y();
    initial_z = mother->get_final_z();
    final_x=initial_x+get_tauf()*px/E;
    final_y=initial_y+get_tauf()*py/E;
    final_z=initial_z+get_tauf()*pz/E;
    phi = atan2(get_py(),get_px());
   }
   else cout << "problem"<<endl;
  }

  void set_initial_x(double initial_x){this->initial_x=initial_x;}
  void set_initial_y(double initial_y){this->initial_y=initial_y;}
  void set_initial_z(double initial_z){this->initial_z=initial_z;}
  double get_initial_x(){return initial_x;}
  double get_initial_y(){return initial_y;}
  double get_initial_z(){return initial_z;}
  double get_final_x(){return final_x;}
  double get_final_y(){return final_y;}
  double get_final_z(){return final_z;}
  double get_phi(){return phi;}
  void set_f(){f=E/this->mother->E;}
  double get_f(){return f;}
  double get_quenched_E(){return quenched_E;}
  double get_initial_E(){
   if(id>3) return f*mother->get_quenched_E();
   if(id<=3) return get_E();
  } 
  
  void set_quench_method(int quench_method){
   this->quench_method=quench_method;
   if(quench_method==0) power=1.;
   else if(quench_method==1) power=2.;
   else if(quench_method==2) power=1.33333;
   else if(quench_method==4){
    power=3.;
    power_t=1.;
   }
  }
  int get_quench_method(){
   return quench_method;
  }
  double get_power(){return power;}
  double get_power_t(){return power_t;}
  
  void set_quenched_E(double quenched_E){
  //!quench_method = 0 -> Light quark quenching 
  //!quench_method = 1 -> Heavy quark quenching 
  //!quench_method = 2 -> Energy independent quenching 
   this->quenched_E=quenched_E;
  } 
  
  double get_tstop(){
   if(mother!=NULL){
    double initial_E = f*(mother->get_quenched_E());
    tstop = 0.2*0.0016*sqrt(initial_E/m)/T;
    return tstop;
   }
   if(mother==NULL){
    return 0;
   }
  }
//!MOTHER RELATED FUNCTIONS
  int get_momID(){return momID;}
  void set_mother(Fragment *mother){this->mother = mother;}
  Fragment * get_mother(){return mother;}
//!DAUGTHER RELATED FUNCTIONS
  void add_daughter(Fragment *daughter)
  {
   if(daughters == NULL)
   {
    daughters = (Fragment **)calloc(2,sizeof(Fragment *));
    daughters[0] = daughter;
	  daughters[1] = NULL;
   }
   else
   {
    daughters = (Fragment **)realloc(daughters, (number_of_daughters+2)*sizeof(Fragment *));
  	daughters[number_of_daughters] = daughter;
	  daughters[number_of_daughters+1] = NULL;
   }
   number_of_daughters++;
  }
  int get_number_of_daughters(){return number_of_daughters;}
  Fragment * get_ieth_daughter(int i){return daughters[i];}
//!ANCESTOR RELATED FUNCTIONS
  void add_ancestor(Fragment *ancestor)
  {
   if(ancestors == NULL)
   {
    ancestors = (Fragment **)calloc(2,sizeof(Fragment *));    
    ancestors[0] = ancestor;    
  	ancestors[1] = NULL;    
   }
   else
   {    
    ancestors = (Fragment **)realloc(ancestors, (number_of_ancestors+2)*sizeof(Fragment *));    
	  ancestors[number_of_ancestors] = ancestor;    
  	ancestors[number_of_ancestors+1] = NULL;    
   }    
   number_of_ancestors++;    
  }
  int get_number_of_ancestors(){return number_of_ancestors;}
  Fragment * get_ieth_ancestor(int i){return ancestors[i];}
  bool is_descendant(Fragment *ancestor_candidate){
   for(int i=0; i<number_of_ancestors; i++){
    if(get_ieth_ancestor(i)->get_id() == ancestor_candidate->get_id()) return true;
   }
  }
//!JET CONE RELATED FUNCTIONS
  void set_incone(bool val){bool_incone=val;}
  void set_incone_of_ieth_jet(int i){jet_index_incone=i;}
  bool is_incone(){return bool_incone;}
  bool is_incone_of_ieth_jet(int i)
  {
  if(jet_index_incone==i) return true;
  else return false;
  }
  void set_incone_ancestor(bool val){incone_ancestor=val;}
  void set_incone_ancestor_of_ieth_jet(int i){jet_index_incone_ancestor=i;}
  bool is_incone_ancestor(){return incone_ancestor;}
  bool is_incone_ancestor_of_ieth_jet(int i){
   if(jet_index_incone_ancestor==i) return true;
   else return false;
  }
  int get_jet_index_incone_ancestor(){return jet_index_incone_ancestor;}
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
struct Incone
{
 private:
  int ind;
  double y;
  double phi;
  double pt;
  double E;
 public:
  Incone(int ind, double y, double phi, double pt, double E)
  {
   this->ind = ind;
   this->y = y;
   this->phi = phi;
   this->pt = pt;
   this->E = E;
  }
  int get_ind(){return ind;}
  double get_y(){return y;}
  double get_phi(){return phi;}
  double get_pt(){return pt;}
  double get_E(){return E;}
};

struct Jet
{
 private:
  double pt;
  int number_of_incones;
  Incone **incones;
  double phi;
 public:
   void reset()
  { 
   pt=0;
   incones=NULL;
   number_of_incones=0;
  }
  Jet()
  {
   reset();
  }
  void set_pt(double pt){this->pt=pt;}
  double get_pt(){return pt;}
  void set_phi(double phi){this->phi=phi;}
  double get_phi(){return phi;}
  int get_number_of_incones(){return number_of_incones;}
  void add_incone(int ind, double y, double phi, double pt, double E)
  {
   Incone *incone = new Incone(ind, y, phi, pt, E);
   if(incones == NULL)
   {
    incones = (Incone **)calloc(2,sizeof(Incone *));
    incones[0] = incone;
  	incones[1] = NULL;
   }
   else
   {
    incones = (Incone **)realloc(incones, (number_of_incones+2)*sizeof(Incone *));
	  incones[number_of_incones] = incone;
	  incones[number_of_incones+1] = NULL;
   }
   number_of_incones++;
  }
  Incone * get_ieth_incone(int i){return incones[i];}
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

class utilities {
 public:
  static void integrate_T(Fragment *fragment, Hydro *hydro, double factor, double Tc){
   int id=fragment->get_id();

   double Ei=fragment->get_initial_E();

   double px=fragment->get_px();
   double py=fragment->get_py();
   double pz=fragment->get_pz();

   Vector3 *w=new Vector3(px/Ei,py/Ei,pz/Ei);
   double wdotw=w->dot(w);
   double E=Ei;
   if(id<=3 || fragment->get_taus()<0.6){
    fragment->set_quenched_E(E);
	  return;
   }
   double integral=0;
   double initial_x=fragment->get_initial_x();
   double initial_y=fragment->get_initial_y();
   double initial_z=fragment->get_initial_z();
   double initial_t=fragment->get_mother()->get_taus();
   double initial_r=sqrt(pow(initial_x,2)+pow(initial_y,2)+pow(initial_z,2));
   double initial_eta=0.5*log((initial_r+initial_z)/(initial_r-initial_z));

   double final_x=fragment->get_final_x();
   double final_y=fragment->get_final_y();
   double final_z=fragment->get_final_z();
   double final_t=fragment->get_taus();
   
   double power=fragment->get_power();
   double power_t=fragment->get_power_t();
   double T = 0;
   double tstop=0;
   int quench_method=fragment->get_quench_method();
   
   double C;
   int QG=fragment->get_QG();
   if(QG==1) C=1; 
   if(QG==2) C=9/4;
   

   double dt=0.01; 
   double t=initial_t;
   if(initial_t<0.6) t=0.6;
   double x=initial_x+(final_x-initial_x)*(t-initial_t)/(final_t-initial_t);
   double y=initial_y+(final_y-initial_y)*(t-initial_t)/(final_t-initial_t);
   double z=initial_z+(final_z-initial_z)*(t-initial_t)/(final_t-initial_t);
   double r=sqrt(pow(x,2)+pow(y,2)+pow(z,2));
   double eta=0.5*log((r+z)/(r-z));
   cout<<"id"<<fragment->get_id()<<" mother id="<<fragment->get_mother()->get_id()<<" mother initial_x"<<fragment->get_mother()->get_initial_x()<<" mother final_x"<<fragment->get_mother()->get_initial_x()<<" t="<<t<<" x="<<x<<" initial_x="<<initial_x<<" final_x="<<final_x<<" initial_y="<<initial_y<<" y="<<y<<" eta="<<eta<<endl;
   T = hydro->get_hydro_T(t,x,y,eta);
   double vx = hydro->get_hydro_vx(t,x,y,eta);
   double vy = hydro->get_hydro_vy(t,x,y,eta);
   double vh = hydro->get_hydro_vh(t,x,y,eta);
   double gamma=cosh(vh);
   double vz=sqrt(1-1/pow(gamma,2.)-pow(vx,2.)-pow(vy,2.));
   
   Vector3 *v= new Vector3(vx,vy,vz);
   double wdotv=v->Vector3::dot(w);
   double vdotv=v->Vector3::dot(v);
   double coef=sqrt(wdotw+pow(gamma,2)*(vdotv-2*wdotv+pow(wdotv,2)));
    
   double initial_tplasma=initial_t*coef;
   double tplasma=initial_tplasma;
   while(t<final_t){    
    if(T<Tc) break;
    if(quench_method!=3){
     integral += pow(T,power)*pow((tplasma-initial_tplasma)*5,power_t)*dt*coef*5;
    }else{
     double EFi=Ei*gamma*(1-wdotv);
     tstop=0.2*pow(EFi/C,1./3.)/(2*factor*pow(T,4./3.)); 
     if((tplasma-initial_tplasma)>tstop){
      E=0;
      break;
     } 
     E += (-((4*EFi/3.141592)*pow(tplasma-initial_tplasma,2))/(pow(tstop,2)*sqrt(pow(tstop,2)-pow(tplasma-initial_tplasma,2))))*dt*coef;
     if(E<0){
      E=0;
      break;
     }
    }
    t=t+dt;
    tplasma=tplasma+dt*coef;
    x=initial_x+(final_x-initial_x)*(t-initial_t)/(final_t-initial_t);
    y=initial_y+(final_y-initial_y)*(t-initial_t)/(final_t-initial_t);
    z=initial_z+(final_z-initial_z)*(t-initial_t)/(final_t-initial_t);
    r=sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    eta=0.5*log((r+z)/(r-z));
     
    T = hydro->get_hydro_T(t,x,y,eta);

    vx = hydro->get_hydro_vx(t,x,y,eta);
    vy = hydro->get_hydro_vy(t,x,y,eta);
    vh = hydro->get_hydro_vh(t,x,y,eta);
       
    gamma=cosh(vh);
    vz=sqrt(1-1/pow(gamma,2.)-pow(vx,2.)-pow(vy,2.));
   
    v->reset(vx,vy,vz);
    wdotv=v->Vector3::dot(w);
    vdotv=v->Vector3::dot(v);
    coef=sqrt(wdotw+pow(gamma,2)*(vdotv-2*wdotv+pow(wdotv,2)));
   }
                         cout<<"integrate_T.17"<<endl;

   if(quench_method == 0){
    E=C*Ei*exp(-factor*integral); 
   }
   else if(quench_method == 1){
    E=Ei-C*factor*integral;
   }
   else if(quench_method==2 && Ei>0.0001){
    double deltaE= 1-2*factor*pow(C,1.33333)*(1/pow(Ei,0.33333))*integral;
	  if(deltaE>=0)E=Ei*sqrt(deltaE);
	  else E=0; 
   }
   else if(quench_method==4){
    E=Ei-C*factor*integral;
   }
   if(E<0) E=0;
   fragment->set_quenched_E(E);
  }
  
  static double embed(Fragment *fragment, Hydro *hydro){

   double t=0.6;
   double maxx=hydro->get_maxx();
   double maxy=hydro->get_maxy();
   double maxh=hydro->get_maxh();
   
   bool passed_value=false;
   double maxT6=0.005;
   double x;
   double y;
   double h;
   double T6,randT6;

   while(!passed_value){
	// cout<<"time"<<time(NULL)<<endl;
    x=2*maxx*((double) (rand()%10000) / 10000.)-maxx;
    y=2*maxy*((double) (rand()%10000) / 10000.)-maxy;
    h=2*maxh*((double) (rand()%10000) / 10000.)-maxh;
    randT6=maxT6*((double) (rand()%10000) / 10000.);
    T6 = pow(hydro->get_hydro_T(t,x,y,h),6);
    if(T6>randT6)passed_value=true;    
   }
   double z=sqrt(pow(x,2)+pow(y,2))*sinh(h);
   fragment->set_initial_x(x);
   fragment->set_initial_y(y);
   fragment->set_initial_z(z);
  }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
struct Event
{
 private:
  Jet **jets;//!SET BY READING FROM THE DATA FILE
  Fragment **fragments;
  Fragment *root;
  int number_of_jets;
  int number_of_fragments;//!SAME AS NF BUT IS AUTOMATICALLY SET WHEN ADDING THE FRAGMENTS (see add_fragment)
 public:
  void reset()
  {
   fragments = NULL;
   jets = NULL;
   root = NULL;
   number_of_fragments=0;
   number_of_jets=0;
  }
  
  Event()
  {
   reset();
  }
  
  ~Event()
  {
   
  }
  
  //!FUNCTIONS TO CONSTRUCT AN EVENT
  void add_fragment(int id, double px, double py, double pz, double E, double m, int momID, int QG, int initialIQG)
  { 
   Fragment *fragment = new Fragment(id, px, py, pz, E, m, momID, QG, initialIQG);
   if(fragments == NULL)
   {
    fragments = (Fragment **)calloc(2,sizeof(Fragment *));
    fragments[0] = fragment;
	  fragments[1] = NULL;
   }
   else
   {
    fragments = (Fragment **)realloc(fragments, (number_of_fragments+2)*sizeof(Fragment *));
	  fragments[number_of_fragments] = fragment;
	  fragments[number_of_fragments+1] = NULL;
   }
   number_of_fragments++;
  }
  
  void add_jet(Jet *jet)
  {
   if(jets == NULL)
   {
    jets = (Jet **)calloc(2,sizeof(Jet *));
    jets[0] = jet;
	  jets[1] = NULL;
   }
   else
   {
    jets = (Jet **)realloc(jets, (number_of_jets+2)*sizeof(Jet *));
	  jets[number_of_jets] = jet;
  	jets[number_of_jets+1] = NULL;
   }
   number_of_jets++;
  }
  
  //!FUNCTIONS TO GET EVENT INFO
  int get_number_of_fragments(){return number_of_fragments;}
  int get_number_of_jets(){return number_of_jets;}
  Fragment * get_ieth_fragment(int i){return fragments[i];}  
  Jet * get_ieth_jet(int i){return jets[i];}  
  Fragment *findFragment_byId(int id){
   for(int i=0; i<number_of_fragments; i++)
   {    
    if(fragments[i]->get_id() == id)
	  {    
     return fragments[i];
    }    
   }    
   return NULL; 
  }
 
  //!FUNCTIONS THAT MAKE THE CONNECTIONS IN THE FAMILY TREE
  
  void findInconeFragments() //!PAIR THE PARTICLES IN JET CONE WITH THE PARTICLES IN THE LONG LIST
  {
   for(int i=0; i<number_of_jets; i++)
   {
    for(int j=0; j<jets[i]->get_number_of_incones(); j++)
	{
     for(int k=0; k<number_of_fragments; k++)
	 {
     if(abs(jets[i]->get_ieth_incone(j)->get_E()-fragments[k]->get_E())<0.001) fragments[k]->set_incone(true);	
     if(abs(jets[i]->get_ieth_incone(j)->get_E()-fragments[k]->get_E())<0.001) fragments[k]->set_incone_of_ieth_jet(i);	
     }
    }
   }
  }
  
  void build_shower()//!SET THE MOTHER AND DAUGHTERS OF EACH PARTICLE IN THE LONG LIST
  {
   for(int i=0; i<number_of_fragments; i++)
   {
   	Fragment* mother = findFragment_byId(fragments[i]->get_momID());
	if(mother !=NULL){
	 fragments[i]->set_mother(mother);
	 fragments[i]->set_f();
	}
    for(int j=0; j<number_of_fragments;j++){
	   if(fragments[j]->get_momID()==fragments[i]->get_id()){
	     fragments[i]->add_daughter(fragments[j]);
	    }
	  }
    }
  }  
  
  void set_ancestors(){//!SET ALL THE ANCESTORS OF EACH PARTICLE IN THE LONG LIST
   for(int i=0; i<number_of_fragments; i++){
    Fragment *ancestor = fragments[i];
    Fragment *elem = fragments[i];    
    while(ancestor->get_id()>3 && ancestor->get_mother()!=NULL){
 	   elem->add_ancestor(ancestor);
     ancestor = ancestor->get_mother();     
     if(ancestor->get_id()==3)
     {elem->add_ancestor(ancestor);          
     }
    }
   }
  }
  
  Fragment* find_root_of_ieth_jet(int ijet){//!FIND THE COMMON ANCESTOR OF ALL THE PARTICLES IN THE LEADING JET CONE
   Fragment * root_candidate = NULL;
   int i;
   for(i=0;i<number_of_fragments;i++){
    if(fragments[i]->is_incone_of_ieth_jet(i)){
     root_candidate=fragments[i];
     break;
    }
   }
   for(int j=i;j<number_of_fragments;j++){   
    if(!(fragments[j]->is_incone_of_ieth_jet(ijet))) continue;
    bool goon=true;
    for(int k=0; k<root_candidate->get_number_of_ancestors() and goon;k++){
     for(int l=0; l<fragments[j]->get_number_of_ancestors() and goon; l++){
      if(root_candidate->get_ieth_ancestor(k)->get_id() == fragments[j]->get_ieth_ancestor(l)->get_id()){
       root_candidate=root_candidate->get_ieth_ancestor(k);
       goon = false;
      }
     }
    }
   }
   return root_candidate;
  }
  
  void find_incone_ancestors(){//!FIND WHICH PARTICLES ARE ANCESTORS OF THE PARTICLES IN THE JET CONES
   for(int l=0;l<number_of_jets;l++){
    for(int i=0; i<number_of_fragments; i++){
     if(fragments[i]->is_incone_of_ieth_jet(l)){
      for(int k=0; k<fragments[i]->get_number_of_ancestors(); k++){
       for(int j=0; j<number_of_fragments; j++){
        if(fragments[i]->get_ieth_ancestor(k)->get_id()==fragments[j]->get_id()){
         fragments[j]->set_incone_ancestor(true);
         fragments[j]->set_incone_ancestor_of_ieth_jet(l);
         continue;
        }
       }
      }
     }    
    }
   }
  }
  
  void set_taus_until_common_ancestors(){//!SET TAUS OF EACH PARTICLE UP TO THE TWO COMMON ANCESTORS OF THE TWO JET CONES 
   for(int i=0; i<number_of_jets;i++){
    Fragment * root = find_root_of_ieth_jet(i);
    if(root->get_mother() != NULL) root->get_mother()->set_taus();
    for(int j=0; j<number_of_fragments; j++){
     Fragment * elem = fragments[j];
     if(elem->is_descendant(root) and elem->is_incone_ancestor()){
      if(fragments[j]->get_m()>1) fragments[j]->set_taus();
     }
    }
   }
  }
  

   void set_all_taus_coordinates_geometry(Hydro *hydro){//!SET TAUS OF EACH PARTICLE UP TO THE HARD SCATTERING 
    Fragment * hard_scattering=findFragment_byId(3);
	utilities::embed(hard_scattering,hydro);

    for(int i=0; i<number_of_fragments; i++){
     Fragment * elem = fragments[i];
     if(elem->is_descendant(hard_scattering) && elem->get_id()!=3){
      fragments[i]->set_taus_geometry();
      if(fragments[i]->get_mother()->get_id()==3){
       if(fragments[i]->get_jet_index_incone_ancestor()>=0) jets[fragments[i]->get_jet_index_incone_ancestor()]->set_phi(fragments[i]->get_phi());
      }
      fragments[i]->set_initial_coordinates();
     }
    }
   }

   void quench_geometry(Hydro* hydro, double factor,int quench_method=0,double Tc=200){
    Fragment * hard_scattering=findFragment_byId(3);
    for(int i=0; i<number_of_fragments; i++){
     Fragment * elem = fragments[i];
     if(elem->is_descendant(hard_scattering) and elem->get_id()!=3){
      fragments[i]->set_quench_method(quench_method);
      utilities::integrate_T(fragments[i], hydro, factor,Tc);
     }
    }
   }
 
  
  bool is_good_event(){//!AN EVENT IS GOOD IF THE PARTICLES IN JET CONE CAN BE FOUND IN THE LONG LIST
   bool good_event=true;
   for(int i=0;i<number_of_jets;i++){
    int nincone=0;
    for(int j=0; j<number_of_fragments; j++){
     if(fragments[j]->is_incone_of_ieth_jet(i)) nincone++;
    }
	  if(jets[i]->get_number_of_incones()!=nincone) good_event=false;
   }
   return good_event;
  }
  
  void print_event(ofstream *myfile) //!PRINTS THE EVENT AS IN THE .DAT FILE
  { 
   for(int i=0; i<number_of_fragments; i++){
    *myfile << fragments[i]->get_id() << " " << fragments[i]->get_px() << " " << fragments[i]->get_py() << " " << fragments[i]->get_pz() << " " << fragments[i]->get_E() << " " << fragments[i]->get_m() << " " << fragments[i]->get_momID() << " " << fragments[i]->get_QG() << " " << fragments[i]->get_initialIQG()<<"\n";
   }
   *myfile << "JETS\n";
   for(int i=0;i<number_of_jets;i++){
    *myfile<<i<<endl;
    *myfile<<jets[i]->get_number_of_incones()<<" "<<jets[i]->get_pt()<<endl;
	  for(int j=0;j<jets[i]->get_number_of_incones();j++){
	   *myfile <<jets[i]->get_ieth_incone(j)->get_ind() << " "<< jets[i]->get_ieth_incone(j)->get_y() << " " << jets[i]->get_ieth_incone(j)->get_phi() << " " << jets[i]->get_ieth_incone(j)->get_pt() << " "<< jets[i]->get_ieth_incone(j)->get_E() << "\n"; ;
	  }
   }
   *myfile << "END\n\n";
  }
};

struct DataFile_Parser
{
 private:
 vector <Event> events;

 public:
 vector <Event> get_event_vector()
 {
  return events;
 } 

 DataFile_Parser(const char *file_path)
 { 
  int ievent=0;
  bool hasjet=false;
  // cout<<"file: "<<file_path<<endl;
  ifstream input_file(file_path);
  int state=0;
  Event current;
  string line;
  state=1;
  while(getline(input_file,line))
  {
   // cout<<"inside the while loop"<<endl;
   trim(line);
   if(line.length()==0) continue;
   // if(ievent==200000) break;
   if(state==1 and string::npos == line.find("JETS"))
   {
    int id;
    double px;
    double py;
    double pz;
    double E;
    double m;
    int momID;
	  int QG;
	  int initialIQG;
    stringstream line_s;
    line_s << line;
	  // cout << "Line length: " << line.length() << endl;
	  // cout << "LINE: " << line << endl;
    line_s >> id >> px >> py >> pz >> E >> m >> momID >> QG >> initialIQG;
    // cout << "Adding fragment with the following params:" << " " << id << " " << px << " " << py << " " << pz << " " << E << " " << m << " " << momID << " " << QG << " " << initialIQG << endl;
    current.add_fragment(id, px, py, pz, E, m, momID,QG,initialIQG);
   }
   else if(state==1 and string::npos != line.find("JETS"))
   {
    state=2;
   }
   else if(state==2 and string::npos == line.find("END"))
   {
    hasjet=true;
    getline(input_file,line);
	  trim(line);
    if(line.length()==0) continue;
    int nparts;
    double pt;
    stringstream line_s;
    line_s << line;
    line_s >> nparts >> pt;
	  Jet *jet = new Jet();
  	jet->set_pt(pt);
  	for(int i=0; i<nparts; i++)
	{
	 getline(input_file,line);
	 trim(line);
   if(line.length()==0) continue;
	 int ind;
   double y;
   double phi;
   double pt;
   double E;
	 stringstream line_s;
   line_s << line; 
   line_s >> ind >> y >> phi >> pt >> E;
	 jet->add_incone(ind, y, phi, pt, E);
	}
	current.add_jet(jet);
   }
   else if(state==2 and string::npos != line.find("END"))
   {
    state=1;
	  events.push_back(current);
    current.reset();
	  if(ievent%100==0) cout<<"number of parsed events = "<<ievent<<endl;
	  ievent++;
   }
   else
   {
    cout <<"Unexpected else. Exiting..." << endl;
    exit(1);
   }
  }
  input_file.close();
 }
};
