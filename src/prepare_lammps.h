#ifndef CBS_PREPARE_LAMMPS_H
#define CBS_PREPARE_LAMMPS_H

#include "struct_definition.h"
#include "molecule.h"
#include "read_coeff.h"

namespace CBS_NS {

class PrepareLAMMPS : protected Pointers {
 public:
  vector<double> Vtemperature;
  int md_label;
  static const int remd,md,em;
  ReadCoeff Rcoeff;
  map<string,string>Mrestraint_type;
  vector<string> Vhead,Vio,Vgroup,Vem,Vfix,Vrun;  //in_script
  vector<string> Voutput_logs,Voutput_trjs;
  vector<double> Vboxsize;
  vector<string> Vbondstyle,Vanglestyle,Vdihedralstyle,Vimproperstyle;
  vector<string> Vatomtype,Vbondtype,Vangletype,Vdihedraltype,Vimpropertype;
  
  PrepareLAMMPS(Cobmms *);
  PrepareLAMMPS(const PrepareLAMMPS &);
  PrepareLAMMPS& operator=(const PrepareLAMMPS &);
  ~PrepareLAMMPS();

  void set_temperature(double, double);
  void set_remd_temperature();
  void set_coeff(const ReadCoeff &);
  vector<string> prepare_data(const Molecule &, const vector<Bond> &, const vector<Angle> &, 
   const vector<Dihedral> &, const vector<Dihedral> &, double =100.0);
  vector<string> prepare_data(const vector<Atom> &, const vector<Bond> &, const vector<Angle> &, 
   const vector<Dihedral> &, const vector<Dihedral> &, double =100.0);
  vector<string> prepare_data_from_amber(const Molecule &, string);
  void update_data_xyz(vector<string> &, const vector<Atom> &);
  void prepare_in_script(string, string, string, int, int);
  void prepare_in_script_allatom(string, string, string, int =20000, int =100);
  void write_in_script(string, int =0);
  vector<string> write_in_script(string, string, size_t, int =20000, int =100);
  vector<string> write_run_script(const vector<string> &, string, int =0);
  void clear();
};

}

#endif