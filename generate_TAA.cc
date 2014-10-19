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

using namespace std;

int main(int argc, char* argv[]){
 

 Hydro* hydro = new Hydro();
 double R = 6.62; 
 double del = 0.556;
 double drho2=hydro->rho2max/(double)hydro->nrho2;
 double dz=hydro->zmax/(double)hydro->nz;
 double maxTAA=0;
 for(int irho=0;irho<hydro->nrho2;irho++){ 
  double T=0;
  for(int k=0;k<hydro->nz; k++){
   double z=dz*((double)k)+dz/2;
   double rho2=drho2*((double)irho)+drho2/2;
   double r=sqrt(rho2+pow(z,2));
   T+=dz/(exp((r-R)/del)+1);
  }
  if(T>maxTAA) maxTAA=T;
  hydro->TAA[irho]=T;
 }
 
 cout<<maxTAA<<endl;
 stringstream fout;
 fout <<"TAA_matrix.bin"; 
 ofstream myFile;
 myFile.open(fout.str().c_str(), ios::out | ios::binary);
 myFile.write((char *)(hydro->TAA), sizeof(hydro->TAA));
 myFile.close();
 
}