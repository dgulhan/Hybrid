#include "hybrid.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

// g++ `root-config --cflags` runQR.cc -o runQR `fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` `root-config --libs`
      
int main(int argc, char* argv[]) {
 double quench_method=0;
 int min_centrality=10;
 int max_centrality=20;
 double factor=0.1;
 int ifile=1; 
 std::stringstream fhydro;
 fhydro <<"../HiranoHydro/PbPb2760_"<<"00"<<"-"<<"05"<<".dat";   
 Hydro *hydro = new Hydro(fhydro.str().c_str());
 cout<<fhydro.str().c_str()<<endl;

 std::stringstream fpythia;
 fpythia<<"/afs/cern.ch/work/d/dgulhan/dataQG/maindata_p"<<ifile<<"_1000evts.txt"<<endl;
 DataFile_Parser *file = new DataFile_Parser(fpythia.str().c_str());//pythia file
 vector <Event> event_vector = file->get_event_vector();

 int event=0;
 for(vector <Event>::iterator it = event_vector.begin(); it != event_vector.end(); ++it){
  if(it->get_number_of_jets()==0) continue;
  cout<<"number of jets = "<<it->get_number_of_jets()<<endl;
  it->findInconeFragments();
  cout<<"in event loop at event number"<<event<<"\n";
  it->build_shower();
  cout << "showers built " << endl;
  it->set_ancestors();
  cout << "ancestors set" << endl;
  it->find_incone_ancestors();
  cout << "incone ancestors found" << endl; 
  it->set_all_taus_coordinates_geometry(hydro);//here
  cout<<"taus and coordinates set"<<endl;
  it->quench_geometry(hydro,factor,quench_method);
  cout<<"quenched"<<endl;
  event++;

  int nfrags = it->get_number_of_fragments();
  for(int i=0;i<nfrags;i++){
   Fragment * elem=it->get_ieth_fragment(i);
   
   double qE = elem->get_quenched_E();
   double p_quench_factor = elem->get_quenched_E()/elem->get_E();
   double qpx = p_quench_factor*elem->get_px();
   double qpy = p_quench_factor*elem->get_py();
   double qpz = p_quench_factor*elem->get_pz(); 
  }   
 }
}