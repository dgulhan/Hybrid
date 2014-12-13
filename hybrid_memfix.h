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
#include "hydro.h"
#include "incone.h"
#include "vector3.h"

using namespace std;

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
  double initial_T;
  double final_T;
  double phi;
  double f;
  double quenched_E;
  double path_lenght;
  double T;
  double tstop;
  int jet_index_incone;
  int jet_index_incone_ancestor;
  int QG_jet;
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
   initial_T=0;
   final_T=0;
   phi=0;
   f=0;
   quenched_E = E;
   path_lenght=0.5;
   T = 0.3;
   tstop = 0;
   jet_index_incone=-1;
   jet_index_incone_ancestor=-1;
   QG_jet=0;
  }

  ~Fragment()
  {
   free(daughters);
   free(ancestors);
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
   if(id>3){
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
  void set_final_x(double final_x){this->final_x=final_x;}
  void set_final_y(double final_y){this->final_y=final_y;}
  void set_final_z(double final_z){this->final_z=final_z;}
  double get_final_x(){return final_x;}
  double get_final_y(){return final_y;}
  double get_final_z(){return final_z;}
  double get_phi(){return phi;}
  void set_f(){f=E/this->mother->E;}
  double get_f(){return f;}
  double get_quenched_E(){return quenched_E;}
  double get_initial_E(){
   if(id>3 && mother->get_id()!=3){
    return f*mother->get_quenched_E();
   }
   else{
    return get_E();
   }
   if(id<=3){ 
    return get_E();
   }
  } 
  void set_quenched_E(double quenched_E){
  //!quench_method = 0 -> Light quark quenching 
  //!quench_method = 1 -> Heavy quark quenching 
  //!quench_method = 2 -> Energy independent quenching 
   this->quenched_E=quenched_E;
  } 
  
  void set_initial_T(double initial_T){
   this->initial_T=initial_T;
  } 
  void set_final_T(double final_T){
   this->final_T=final_T;
  } 
  double get_initial_T(){
   return initial_T;
  } 
  double get_final_T(){
   return final_T;
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
   return false;
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
  void delete_ancestors(){
   for(int i=0;i<number_of_ancestors;i++){
    free(ancestors);
   }
  }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------


struct Jet
{
 private:
  double pt;
  int number_of_incones;
  Incone **incones;
  double phi;
 public:
  Jet()
  { 
   pt=0;
   incones=NULL;
   number_of_incones=0;
  }

  ~Jet()
  {
   for(int i=0; i<number_of_incones; i++)
   {
    delete incones[i];
   }
   free(incones);
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
  Event()
  {
   fragments = NULL;
   jets = NULL;
   root = NULL;
   number_of_fragments=0;
   number_of_jets=0;
  }
  
  ~Event()
  {
   for(int i=0; i<number_of_fragments; i++)
   {
    delete fragments[i];
   }
   free(fragments);
   free(jets);
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
    if(mother !=NULL)
    {
     fragments[i]->set_mother(mother);
     fragments[i]->set_f();
    }
    for(int j=0; j<number_of_fragments;j++)
    {
     if(fragments[j]->get_momID()==fragments[i]->get_id())
     {
      fragments[i]->add_daughter(fragments[j]);
     }
    }
   }
  }  
  
  void set_ancestors(){//!SET ALL THE ANCESTORS OF EACH PARTICLE IN THE LONG LIST
   for(int i=0; i<number_of_fragments; i++){
    Fragment *ancestor = fragments[i];
    Fragment *elem = fragments[i];    
    while(ancestor->get_id()>3 && ancestor->get_mother()!=NULL)
    {
     elem->add_ancestor(ancestor);
     ancestor = ancestor->get_mother();     
     if(ancestor->get_id()==3)
     {
      elem->add_ancestor(ancestor);          
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
  
  void integrate_T(Fragment *fragment, Hydro *hydro, double factor, double Tc){
   int id=fragment->get_id();
   // cout<<"id="<<fragment->get_id()<<endl;
   double Ei=fragment->get_initial_E();
   // cout<<"id="<<id<<" Ei="<<Ei<<" mother Ei"<<fragment->get_mother()->get_quenched_E()<<" mother id="<<fragment->get_mother()->get_id()<<endl;
   double px=fragment->get_px();
   double py=fragment->get_py();
   double pz=fragment->get_pz();
   double original_E=fragment->get_E();
   Vector3 w(px/original_E,py/original_E,pz/original_E);
   double wdotw=w.dot(w);
   double E=Ei;
   if(id<=3){
    fragment->set_quenched_E(E);
	  return;
   }
   //bool first=true;
   //double integral=0;
   double initial_x=fragment->get_initial_x();
   double initial_y=fragment->get_initial_y();
   double initial_z=fragment->get_initial_z();
   double initial_t=fragment->get_mother()->get_taus();
   //double initial_r=sqrt(pow(initial_x,2)+pow(initial_y,2)+pow(initial_z,2));
   //double initial_eta=0.5*log((initial_r+initial_z)/(initial_r-initial_z));

   double final_x=fragment->get_final_x();
   double final_y=fragment->get_final_y();
   double final_z=fragment->get_final_z();
   double final_t=fragment->get_taus();
   
   if(sqrt(pow(final_t,2)-pow(final_z,2))<0.6){
    fragment->set_quenched_E(E);
	  return;
   }
   
   double power=hydro->get_power();
   double power_t=hydro->get_power_t();
   double T = 0;
   double Tplus = 0;
   double tstop=0;
   int quench_method=hydro->get_quench_method();
   
   double C;
   int QG=fragment->get_QG();
   if(QG==1) C=1; 
   if(QG==2) C=9./4.;
   

   double dt=0.01; 
   double t=initial_t;
   //double tau_initial=sqrt(pow(initial_t,2)-pow(initial_z,2));
   double tplasma=0;
   double final_T=0;
   double x,y,z,tau,eta,xplus,yplus,zplus,tauplus,etaplus;
   while(t<final_t && QG<3){    
    if((final_t-t)<dt){
	  dt=final_t-t;
	 }
   x=initial_x+(final_x-initial_x)*(t-initial_t)/(final_t-initial_t);
   y=initial_y+(final_y-initial_y)*(t-initial_t)/(final_t-initial_t);
   z=initial_z+(final_z-initial_z)*(t-initial_t)/(final_t-initial_t);
   tau=sqrt(pow(t,2)-pow(z,2));
   eta=0.5*log((t+z)/(t-z));
 
   xplus=initial_x+(final_x-initial_x)*(t+dt-initial_t)/(final_t-initial_t);
   yplus=initial_y+(final_y-initial_y)*(t+dt-initial_t)/(final_t-initial_t);
   zplus=initial_z+(final_z-initial_z)*(t+dt-initial_t)/(final_t-initial_t);
   tauplus=sqrt(pow(t+dt,2)-pow(zplus,2));
   etaplus=0.5*log((t+dt+zplus)/(t+dt-zplus));
	
   if((tau+tauplus)/2<0.6){
    t=t+dt;
    continue;
   }
   else{
    T = hydro->get_hydro_T(tau,x,y,eta);
    Tplus = hydro->get_hydro_T(tauplus,xplus,yplus,etaplus);
	  T=(T+Tplus)/2;
    if(T<(Tc/1000)) break;

    double vx = hydro->get_hydro_vx(tau,x,y,eta);
    double vy = hydro->get_hydro_vy(tau,x,y,eta);
    double vh = hydro->get_hydro_vh(tau,x,y,eta);
    double vz=tanh(vh);

    Vector3 v(vx,vy,vz);
    double wdotv=v.dot(w);
    double vdotv=v.dot(v);
    
    if(vdotv>=1) vdotv=0.9999999;
    double gamma=1/sqrt(1-vdotv);

    double coef=sqrt(wdotw+pow(gamma,2)*(vdotv-2*wdotv+pow(wdotv,2)));
    if(quench_method==0){
     E += -C*factor*pow(T,power)*dt*5;
    }
    else if(quench_method==1){
     E += -C*factor*pow(T,power)*pow((tplasma)*5,power_t)*dt*5;
	  }else if(quench_method==2){
     double EFi=Ei*gamma*(1-wdotv);
     tstop=0.2*pow(EFi/C,1./3.)/(2*factor*pow(T,4./3.)); 
     if(tplasma>tstop){
      E=0;
      break;
     } 
     E += (-((4*EFi/3.141592)*pow(tplasma,2))/(pow(tstop,2)*sqrt(pow(tstop,2)-pow(tplasma,2))))*dt;
    }
    t+=dt;
    tplasma+=dt*coef;
    final_T=T;
   }
   if(E<0){
    E=0;
    break;
   }
   }
   fragment->set_final_T(final_T);
   fragment->set_quenched_E(E);
  }
  
  double embed(Fragment *fragment, Hydro *hydro){
   bool do_minimum_bias_b=false;
   bool do_embed_ncoll=false;
   bool do_fixed_b=false;
   if(!do_minimum_bias_b) do_fixed_b=false;
   bool passed_value=false;
   double maxNcoll=hydro->get_max_ncoll();
   double x=0;
   double y=0;
   //double h=0; 
   double z=0;
   double b=0;
   //double ncoll,randNcoll;
   
   double db=hydro->bmax/((double)(hydro->nb));
   double dx=hydro->xmax/((double)(hydro->nx));
   //double dy=hydro->ymax/((double)(hydro->ny));
   //double drho2=hydro->rho2max/((double)(hydro->nrho2));
   if(do_minimum_bias_b){
    double b2=(pow(hydro->get_bmax_bin(),2)-pow(hydro->get_bmin_bin(),2))*((double) (rand()%10000) / 10000.)+pow(hydro->get_bmin_bin(),2);
    b=sqrt(b2);
    if(do_fixed_b) b=2.5;
   }    
   int ib = (int)(b/db);

   if(do_embed_ncoll){
    while(!passed_value){
     y=2*hydro->ymax*((double) (rand()%10000) / 10000.)-hydro->ymax;
     x=(2*hydro->xmax-b)*((double) (rand()%10000) / 10000.)-(2*hydro->xmax-b)/2.;
     if(!do_minimum_bias_b){
      b=(hydro->get_bmax_bin()-hydro->get_bmin_bin())*((double) (rand()%10000) / 10000.)+hydro->get_bmin_bin();
      ib = (int)(b/db);
     }
     double randNcoll=maxNcoll*((double) (rand()%10000) / 10000.);
    
     int ix = (int)(fabs(x)/dx);
     int iy = (int)(fabs(y)/dx);
	  
     double ncoll=0;
     if(ix<hydro->nx && iy<hydro->ny && ib<hydro->nb) ncoll=hydro->Ncoll[ix][iy][ib];
     if(ncoll>randNcoll)passed_value=true;    
    } 
   }else{
    while(!passed_value){
     double rho2=hydro->rho2max*((double) (rand()%10000) / 10000.);
     double phi=2*3.14159265359*((double) (rand()%10000) / 10000.);
	  
     if(!do_minimum_bias_b){
      double b2=(pow(hydro->get_bmax_bin(),2)-pow(hydro->get_bmin_bin(),2))*((double) (rand()%10000) / 10000.)+pow(hydro->get_bmin_bin(),2);
      b=sqrt(b2);
     }
	  
     x=cos(phi)*sqrt(rho2);
     y=sin(phi)*sqrt(rho2);
     double xm=x-b/2;
     double xp=x+b/2;
	 	  
     double rhom=pow(xm,2)+pow(y,2);
     double rhop=pow(xp,2)+pow(y,2);
	
     if((rhom)>(hydro->rho2max)) continue;
     if((rhop)>(hydro->rho2max)) continue;

     double drho2=hydro->rho2max/(double)hydro->nrho2;
	
     int irhom=(int)(rhom/drho2);
     int irhop=(int)(rhop/drho2);
	
     double ncoll=hydro->TAA[irhom]*hydro->TAA[irhop];
     double randNcoll=maxNcoll*((double) (rand()%10000) / 10000.);
     if(ncoll>randNcoll) passed_value=true;
    }
   }
   fragment->set_final_x(x);
   fragment->set_final_y(y);
   fragment->set_final_z(z);
   
   return b;
  }

   void set_all_taus_coordinates_geometry(Hydro *hydro){//!SET TAUS OF EACH PARTICLE UP TO THE HARD SCATTERING 
    Fragment * hard_scattering=findFragment_byId(3);
    embed(hard_scattering,hydro);
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
    hydro->set_quench_method(quench_method);
    for(int i=0; i<number_of_fragments; i++){
     Fragment * elem = fragments[i];
     if(elem->is_descendant(hard_scattering) and elem->get_id()!=3){
      integrate_T(fragments[i], hydro, factor, Tc);
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
    *myfile << fragments[i]->get_id() << " " 
            << fragments[i]->get_px() << " " 
            << fragments[i]->get_py() << " " 
            << fragments[i]->get_pz() << " " 
            << fragments[i]->get_E() << " " 
            << fragments[i]->get_m() << " " 
            << fragments[i]->get_momID() << " " 
            << fragments[i]->get_QG() << " " 
            << fragments[i]->get_initialIQG()<<"\n";
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
 ifstream input_file;
 
public:
 Event *get_next_event()
 {
  Event *current = new Event;
  string line;
  while(getline(input_file,line))
  {
   trim(line);
   if(line.length()==0) continue;
   if(string::npos == line.find("0.000123"))
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
    line_s >> id >> px >> py >> pz >> E >> m >> momID >> QG >> initialIQG;
    current->add_fragment(id, px, py, pz, E, m, momID,QG,initialIQG);
   }
   else if(string::npos != line.find("0.000123"))
   {
    return current;
   }
   else
   {
    cout <<"Unexpected else. Exiting..." << endl;
    exit(1);
   }
  }
  return NULL;
 }

 DataFile_Parser(const char *file_path)
 {
  input_file.open(file_path);  
 }
 
 ~DataFile_Parser()
 {
  input_file.close();
 }
};
