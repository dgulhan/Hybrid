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
     

int main(int argc, char* argv[]){

 const int n_centrality=9;

 std::string min_centrality_bins[]={"00","05","10","20","30","40","50","60","70"};
 std::string max_centrality_bins[]={"05","10","20","30","40","50","60","70","80"};

 for(int i_cent=0;i_cent<n_centrality;i_cent++){
  std::stringstream fhydro;
  fhydro <<"../HiranoHydro/PbPb2760_"<<min_centrality_bins[i_cent]<<"-"<<max_centrality_bins[i_cent]<<".dat";   
  cout<<"Hydro file:"<<fhydro.str().c_str()<<endl;
  Hydro *hydro = new Hydro(1,fhydro.str().c_str());
    
  std::stringstream fout;
  fout <<"../HiranoHydro/PbPb2760_"<<min_centrality_bins[i_cent]<<"-"<<max_centrality_bins[i_cent]<<".bin"; 
  ofstream myFile;
  myFile.open(fout.str().c_str(), ios::out | ios::binary);

  myFile.write((char *)(hydro->hydro_T), sizeof(hydro->hydro_T));
  myFile.write((char *)(hydro->hydro_E), sizeof(hydro->hydro_E));
  myFile.write((char *)(hydro->hydro_vx), sizeof(hydro->hydro_vx));
  myFile.write((char *)(hydro->hydro_vy), sizeof(hydro->hydro_vy));
  myFile.write((char *)(hydro->hydro_vh), sizeof(hydro->hydro_vh));

  myFile.close();

  delete hydro;
 }
}