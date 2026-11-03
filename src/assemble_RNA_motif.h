#ifndef CBS_ASSEMBLE_RNA_MOTIF_H
#define CBS_ASSEMBLE_RNA_MOTIF_H

#include "RNA_motif.h"

namespace CBS_NS {

class AssembleRNAMotif : protected Pointers {
 public:
  RNAMotif motif_adjust;
  vector<string> overlapped_resenames;
  double moved_distance;

  AssembleRNAMotif(Cobmms *);
  AssembleRNAMotif(const AssembleRNAMotif &);
  AssembleRNAMotif& operator=(const AssembleRNAMotif &);
  ~AssembleRNAMotif();

  void assemble_RNA_motif(const vector<RNAMotif> &, const RNAMotif &, int =1); //with overlapped residues
  void place_RNA_motif(const vector<RNAMotif> &, vector<RNAMotif> &); //without overlapped residues
  void place_RNA_motif(const vector<RNAMotif> &, vector<vector<RNAMotif> > &);
  static int calculate_number_of_clash_atoms(const vector<Atom> &, const vector<Atom> &);
  void clear();
};

}

#endif