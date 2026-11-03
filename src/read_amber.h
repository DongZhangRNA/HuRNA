#ifndef CBS_READ_AMBER_H
#define CBS_READ_AMBER_H

#include "pointers.h"
#include "molecule.h"
#include "struct_definition.h"

namespace CBS_NS {

class ReadAmber : protected Pointers {
 public:
  vector<Atom> Vatoms;
  vector<Bond> Vbonds;
  vector<Angle> Vangles;
  vector<Dihedral> Vdihedrals;
  vector<Dihedral> Vimpropers;
  vector<string> Vmass;
  vector<string> Vbondcoeff,Vanglecoeff,Vdihedralcoeff,Vimpropercoeff,Vpaircoeff;
  vector<double> Vboxsize;
  
  ReadAmber(Cobmms *);
  ReadAmber(const ReadAmber &);
  ReadAmber& operator=(const ReadAmber &);
  ~ReadAmber();

  static bool check_tleap();
  void prepare_amber_files(const Molecule &, string, int =0);
  void prepare_amber_files(string, string, int =0);
  void prepare_gromacs_files(const Molecule &, string, const GmxBoxSet &, int =0);
  void prepare_gromacs_files(string, string, const GmxBoxSet &, int =0);
  void read_prmtop(string);
  void read_inpcrd(string);
  void clear();

 private:
  pair<int,int> parse_format(string);
  vector<string> parse_line(string, int);
};

}

#endif