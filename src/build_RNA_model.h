#ifndef CBS_BUILD_RNA_MODEL_H
#define CBS_BUILD_RNA_MODEL_H

#include "RNA_motif.h"
#include "RNA_molecule.h"
#include "read_restraint.h"

namespace CBS_NS {class BuildRNAModel;}
ostream& operator <<(ostream &, CBS_NS::BuildRNAModel &);

namespace CBS_NS {

class BuildRNAModel : protected Pointers {
 public:
  vector<vector<RNAMotif> > models;
  vector<vector<RNAMotif> > adjusted_motifs;
  vector<vector<string> > Vtemplate_information;

  BuildRNAModel(Cobmms *);
  BuildRNAModel(const BuildRNAModel &);
  BuildRNAModel& operator=(const BuildRNAModel &);
  ~BuildRNAModel();

  double build_RNA_model(const vector<RNAMotif> &, ReadRestraint &);
  void clear();
  friend ostream& ::operator <<(ostream &, BuildRNAModel &);

 private:
  RNAMolecule RNAHelix, RNALoop;
  size_t max_helix_length,max_loop_length;

  void read_RNA_helix_loop();
  void build_Helix(RNAMotif &);
  void build_Long_loop(RNAMotif &);
  int build_Hairpin(RNAMotif &);
  int build_Tail(RNAMotif &);
  int build_Junction(RNAMotif &);
  void build_Loop(RNAMotif &);
};

}

#endif