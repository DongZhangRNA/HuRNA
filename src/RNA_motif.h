#ifndef CBS_RNA_MOTIF_H
#define CBS_RNA_MOTIF_H

#include "molecule.h"

namespace CBS_NS {class RNAMotif;}
ostream& operator <<(ostream &, CBS_NS::RNAMotif &);

namespace CBS_NS {

class RNAMotif : public Molecule {
 public:
  string type;
  vector<string> sequences;
  vector<string> ends;
  vector<int> lengths;
  vector<string> secstrs;
  int flag_rigid;

  RNAMotif(Cobmms *);
  RNAMotif(const RNAMotif &);
  RNAMotif& operator=(const RNAMotif &);
  ~RNAMotif();

  void set_RNA_motif(const RNAMotif &);
  void insert_RNA_segment(const Chain &, int =-2);
  void set_RNA_motif(string, int =0);
  void get_RNA_motif_sequence();
  string get_RNA_motif_name() const;
  string get_RNA_motif_info() const;
  void set_RNA_motif_secondary_structure();
  void set_RNA_motif_xyz(RNAMotif &);
  bool is_identical_junction(RNAMotif &);
  void clear();
  size_t size() const;
  Chain& operator [](size_t);
  friend ostream& ::operator <<(ostream &, RNAMotif &);
};

}

#endif