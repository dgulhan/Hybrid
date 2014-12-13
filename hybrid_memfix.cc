#include "hybrid_memfix.h"
// #include "TH1.h"

// g++ `root-config --cflags` runQR.cc -o runQR `fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` `root-config --libs`
      
int main(int argc, char* argv[]) {
 srand((unsigned)time(NULL));

 int quench_method=atoi(argv[1]);
 cout<<"quench_method="<<quench_method<<endl;

 if(quench_method>2){
  cout<<"argv[1] cannot be larger than 2"<<endl;
  return 0;
 }
 int centrality_bin=atoi(argv[2]);
 cout<<"centrality_bin="<<centrality_bin<<endl;
 if(centrality_bin>8){
  cout<<"argv[2] cannot be larger than 8"<<endl;
  return 0;
 }
 
 double factor=atof(argv[3]);
 cout<<"factor="<<factor<<endl;
 double Tc=atof(argv[4]);
 cout<<"Tc="<<Tc<<endl;
 
 char * file_number=argv[5];
 
 Hydro *hydro = new Hydro(0,centrality_bin,"Ncoll_matrix.bin","TAA_matrix.bin");
 
 
 std::string min_centrality_bins[]={"00","05","10","20","30","40","50","60","70"};
 std::string max_centrality_bins[]={"05","10","20","30","40","50","60","70","80"};

 for(int ifile=0; ifile<1;ifile++){ 
  std::stringstream fpythia;
  fpythia<<"/data/dgulhan/Adsdata_pthat50/TREES/full";
  DataFile_Parser *file = new DataFile_Parser(fpythia.str().c_str());//pythia file
  cout <<fpythia.str().c_str()<<endl;

  std::stringstream fout;
  fout << "Outfile_m"<<quench_method<<"_dijet_pthat50/"<< "datafile_jetfile_" << file_number <<"_method_"<<quench_method<< "_centrality_" << min_centrality_bins[centrality_bin] << "_"<< max_centrality_bins[centrality_bin] << "_factor" << (factor) << "_Tc"<<Tc<<".txt";   
  ofstream myfile;
  myfile.open(fout.str().c_str());

  int event=0;
  while(true) {
   Event *it = file->get_next_event();
   if(!it) break;

   event++;  
   //cout<<"event="<<event<<endl;

   //cout<<"number of jets = "<<it->get_number_of_jets()<<endl;
   it->findInconeFragments();
   it->build_shower();
   //cout << "showers built " << endl;
   it->set_ancestors();
   //cout << "ancestors set" << endl;
   it->find_incone_ancestors();
   //cout << "incone ancestors found" << endl; 
   it->set_all_taus_coordinates_geometry(hydro);//here
   //cout<<"taus and coordinates set"<<endl;
   it->quench_geometry(hydro,factor,quench_method,Tc);
   //cout<<"quenched"<<endl;
      
   int nfrags = it->get_number_of_fragments();
   
   for(int i=0;i<nfrags;i++) {
    Fragment * elem=it->get_ieth_fragment(i);
    if(elem->get_taus()<50){
     continue;
    }
    double p_quench_factor = elem->get_quenched_E()/elem->get_E();
    myfile << elem->get_id() << " " 
           << elem->get_px() << " "
           << elem->get_py() << " " 
           << elem->get_pz() << " " 
           << elem->get_E()  << " "
           << elem->get_px()*p_quench_factor << " "
           << elem->get_py()*p_quench_factor << " " 
           << elem->get_pz()*p_quench_factor << " " 
           << elem->get_quenched_E() << " "
           << elem->get_QG()*elem->get_initialIQG() << "\n";
    if(i==0)cout<<"event="<<event<<endl;
   }
   
   //cout<<"event="<<event<<endl;
   myfile<<"END\n";
   delete it;
  }
  myfile.close();
  delete file;
 }
 delete hydro;
}
