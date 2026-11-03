#ifndef CBS_READ_LAMMPS_H
#define CBS_READ_LAMMPS_H

#include "atom.h"
#include "pointers.h"
#include "struct_definition.h"

namespace CBS_NS {

class ReadLAMMPS : protected Pointers {
 public:
  vector<vector<vector<Atom> > >VVatoms;
  vector<vector<int> > VVsteps;
  vector<vector<string> > VVenergies;  //step+" "+energy
  map<string,string> atomlabel;
  vector<Bond> Vbonds;
  vector<Angle> Vangles;
  vector<Dihedral> Vdihedrals; 
  vector<Dihedral> Vimpropers;
  int step_begin,step_end;
  
  ReadLAMMPS(Cobmms *);
  ReadLAMMPS(const ReadLAMMPS &);
  ReadLAMMPS& operator=(const ReadLAMMPS &);
  ~ReadLAMMPS();

  void read_log(string, int =0, int =-1);
  void read_log(const vector<string> &, int =0, int =-1);
  void read_lammpstrj(string, int =0, int =-1);
  void read_lammpstrj(const vector<string> &, int =0, int =-1);
  void recover_residue_info(string);
  void recover_residue_info(string, vector<Atom> &);
  vector<vector<Atom> > extract_top_structures(double, double =1.0, int =1) const;
  vector<Atom> extract_lowest_energy_structure() const;
  void read_data(string);
  void clear();
};

}

#endif