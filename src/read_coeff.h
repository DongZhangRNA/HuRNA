#ifndef CBS_READ_COEFF_H
#define CBS_READ_COEFF_H

#include "pointers.h"

namespace CBS_NS {

class ReadCoeff : protected Pointers {
 public:
  vector<string> Vmass,Vbondcoeff,Vanglecoeff,Vdihedralcoeff,Vimpropercoeff,Vpaircoeff;
  map<string,string> Mtype_charge;
  map<string,string> Mpair_style_args;
  string default_pair_style;
  
  ReadCoeff(Cobmms *);
  ReadCoeff(const ReadCoeff &);
  ReadCoeff& operator=(const ReadCoeff &);
  ~ReadCoeff();

  void read_RNA_coeff();
  void add_coeff(string);
  void reset_paircoeff_lj();
  void reset_paircoeff_soft(double =10.0, double = 0.0);
  void get_interaction_style();
  void add_extra_RNA_basepair_coeff(string);
  vector<string> get_weighted_RNA_basepair_coeff(string, double) const;
  vector<string> get_pair_coeff() const;
  void clear();
};

}

#endif