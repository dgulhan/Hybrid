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
 

 Hydro* hydro=new Hydro();
 double R=7.11; 
 double del = 0.5;
 double sigma = 7.2;//fm^2
  cout<<3<<endl;
 double Ncollmax=0;
 double dx=hydro->xmax/(double)hydro->nx;
 double dy=hydro->ymax/(double)hydro->ny;
 double db=hydro->bmax/(double)hydro->nb;
 double dz=hydro->zmax/(double)hydro->nz;
 for(int ib=0;ib<hydro->nb;ib++){ 
  for(int i=0; i<hydro->nx; i++){
   for(int j=0; j<hydro->ny; j++){
    double Tm=0;   
    double Tp=0;
    double xm=dx*((double)i)-db*((double)ib)/2.;
    double xp=dx*((double)i)+db*((double)ib)/2.;
    double y=dy*((double)j);    
    for(int k=0;k<hydro->nz; k++){
     double z=dz*((double)k);
     double rm=sqrt(pow(xm,2)+pow(y,2)+pow(z,2));
     double rp=sqrt(pow(xp,2)+pow(y,2)+pow(z,2));
     Tm=1/(exp((rm-R)/del)+1);
     Tp=1/(exp((rp-R)/del)+1);
    }
    hydro->Ncoll[i][j][ib]=sigma*208*208*Tm*Tp;
    if(Ncollmax<sigma*208*208*Tm*Tp) Ncollmax=sigma*208*208*Tm*Tp;
   }
  }
 }
  cout<<Ncollmax<<endl;

 stringstream fout;
 fout <<"Ncoll_matrix.bin"; 
 ofstream myFile;
 myFile.open(fout.str().c_str(), ios::out | ios::binary);
 myFile.write((char *)(hydro->Ncoll), sizeof(hydro->Ncoll));
 myFile.close();
 
}