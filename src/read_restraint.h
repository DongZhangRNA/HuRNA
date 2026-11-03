#ifndef CBS_READ_RESTRAINT_H
#define CBS_READ_RESTRAINT_H

#include "struct_definition.h"
#include "RNA_motif.h"

namespace CBS_NS {class ReadRestraint;}
ostream& operator <<(ostream &, CBS_NS::ReadRestraint &);

namespace CBS_NS {

class ReadRestraint : protected Pointers {
 public:
  vector<RNAMotif> RNA_templates;
  vector<ResiduePair> RNA_basepairs;
  vector<ResiduePair> RNA_basestacks;
  vector<vector<ResiduePair> > rigid_bodies;
  vector<DistanceRestraint> distance_restraints;

  ReadRestraint(Cobmms *);
  ReadRestraint(const ReadRestraint &);
  ReadRestraint& operator=(const ReadRestraint &);
  ~ReadRestraint();

  void read_restraint(string);
  void add_RNA_template(const vector<string> &);
  void add_rigid_body(string);
  void add_RNA_basepair(string);
  void add_RNA_basestack(string);
  void add_distance_restraint(string);
  void clear();
  friend ostream& ::operator <<(ostream &, ReadRestraint &);
};

}

#endif