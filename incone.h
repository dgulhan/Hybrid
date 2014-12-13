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
