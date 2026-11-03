#ifndef CBS_STRUCT_DEFINITION_H
#define CBS_STRUCT_DEFINITION_H

namespace CBS_NS {

struct Bond {
  string type;
  int atom_i,atom_j;
};

struct Angle {
  string type;
  int atom_i,atom_j,atom_k;
};

struct Dihedral {
  string type;
  int atom_i,atom_j,atom_k,atom_l;
};

struct ResiduePair {
  string residue_i,residue_j;
};

struct DistanceRestraint {
  string atom_i,atom_j;
  string type,coeff;
};

struct InteractionPair{
  string residue_i,residue_j;
  string type;
};

struct ConfusionMatrixScore {
  size_t TP = 0;
  size_t FP = 0;
  size_t FN = 0;

  double PPV = 0.0;
  double STY = 0.0;
  double INF = 0.0;
  double F1 = 0.0;
};

struct GmxTopMole {
  vector<string> Vheads;
  vector<vector<string> > Vresidues_info;
  vector<string> Vbonds_info,Vpairs_info,Vangles_info,Vdihedrals_info,Vimpropers_info;
  vector<string> Vresnames_info;
  vector<string> Vposres_info;
  vector<string> Vothers;
};

struct GmxBoxSet {
  double distance,concentration;
  int neutralization;
  string positivename,negativename;
  string boxshape;
};

}

#endif