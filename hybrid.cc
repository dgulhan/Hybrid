#include "hybrid.h"
// #include "TH1.h"
// #include "TH2D.h"
// #include "TH1D.h"
// #include "TFile.h"

// g++ `root-config --cflags` runQR.cc -o runQR `fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` `root-config --libs`
      
int main(int argc, char* argv[]) {
 srand((unsigned)time(NULL));

 int quench_method=atoi(argv[1]);
 if(quench_method>4){
  cout<<"argv[1] cannot be larger than 4"<<endl;
  return 0;
 }
 int centrality_bin=atoi(argv[2]);
 if(centrality_bin>8){
  cout<<"argv[2] cannot be larger than 8"<<endl;
  return 0;
 }
 
 cout<<"centrality_bin="<<centrality_bin<<endl;
 cout<<"quench_method="<<quench_method<<endl;
 double factor=atof(argv[3]);
 cout<<"factor="<<factor<<endl;
 double Tc=atof(argv[4]);
  cout<<"Tc="<<Tc<<endl;

 std::string min_centrality_bins[]={"00","05","10","20","30","40","50","60","70"};
 std::string max_centrality_bins[]={"05","10","20","30","40","50","60","70","80"};

 std::stringstream fhydro;
 fhydro <<"../HiranoHydro/PbPb2760_"<<"00"<<"-"<<"05"<<".dat";   
 Hydro *hydro = new Hydro(fhydro.str().c_str());
 cout<<fhydro.str().c_str()<<endl;

 for(int ifile=0; ifile<1;ifile++){
  std::stringstream fpythia;
  // fpythia<<"/afs/cern.ch/work/d/dgulhan/dataQG/maindata_p"<<ifile<<"_1000evts.txt";
  fpythia<<"oneevt.dat";
  DataFile_Parser *file = new DataFile_Parser(fpythia.str().c_str());//pythia file
  vector <Event> event_vector = file->get_event_vector();
  cout <<fpythia.str().c_str()<<endl;
 
  std::stringstream fout;
  fout << "Outfile/datafile_jetfile_" << ifile << "_centrality_" << min_centrality_bins[centrality_bin] << "_"<< max_centrality_bins[centrality_bin] << "_factor" << (factor) << "_Tc"<<Tc<<".txt";   
  ofstream myfile;
  myfile.open(fout.str().c_str());

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
   it->quench_geometry(hydro,factor,quench_method,Tc);
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
    myfile <<elem->get_id()<<" " <<elem->get_px() << " "<< elem->get_py() << " " << elem->get_pz() << " " << elem->get_E() <<" "<<elem->get_px()*p_quench_factor << " "<< elem->get_py()*p_quench_factor << " " << elem->get_pz()*p_quench_factor << " " << elem->get_quenched_E() <<" "<<elem->get_QG()<<"\n";

   }   
  }
  myfile.close();
 }
 // fhydro.close();
 // fpythia.close();
}