#ifndef CBS_RNA_RECOVERY_H
#define CBS_RNA_RECOVERY_H

#include "molecule.h"
#include "superpose.h"

namespace CBS_NS {

class RNARecovery : protected Pointers{
 public:
  Molecule rna_aa;

  RNARecovery(Cobmms *);
  RNARecovery(const RNARecovery &);
  RNARecovery& operator=(const RNARecovery &);
  ~RNARecovery();

  void recover_RNA_to_AA(const Molecule &);
  void recover_RNA_to_AA_v2(const Molecule &);
  void clear();

 private:
  vector<vector<Atom> > Vbackbones_AA_atoms;
  vector<vector<double> > Vbackbones_CG_xyz;
  vector<vector<double> > Vbackbones_CG_N_xyz; // N9/N1 atom for extra alignment
  vector<double> Vetas;
  vector<double> Vdistance_SpPns,Vdistance_SBNs;
  vector<double> Vthetas;
  map<string,vector<Atom> > Mbase_AA_atoms;
  map<string,vector<double> > Mbase_CG_xyz;

  void read_backbone_conformer();
  void set_base_conformer();
  static int calculate_RMSD(Superpose &, const vector<double> &, vector<double> &, int); 
};

}

#endif