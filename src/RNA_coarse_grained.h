#ifndef CBS_RNA_COARSE_GRAINED_H
#define CBS_RNA_COARSE_GRAINED_H

#include "struct_definition.h"
#include "RNA_molecule.h"

namespace CBS_NS {

class RNACoarseGrained : protected Pointers{
 public:
  Molecule rna_cg;
  vector<Bond> bonds;
  vector<Angle> angles;
  vector<Dihedral> dihedrals;
  vector<Dihedral> impropers;

  RNACoarseGrained(Cobmms *);
  RNACoarseGrained(const RNACoarseGrained &);
  RNACoarseGrained& operator=(const RNACoarseGrained &);
  ~RNACoarseGrained();

  void convert_RNA_to_CG(const Molecule &, size_t =0);
  void convert_RNA_to_CG(const RNAMolecule &, size_t =0);
  void convert_RNA_to_CG_v2(const Molecule &, size_t =0);
  void convert_RNA_to_CG_v2(const RNAMolecule &, size_t =0);
  void reset_RNA_CG(const vector<Atom> &, size_t =0);
  void reset_RNA_CG(const vector<Atom> &, const vector<string> &, size_t =0);
  void reset_RNA_CG(const Molecule &);
  void build_RNA_topology();
  void build_RNA_topology_v2();
  void add_circular_RNA_topology();
  void add_circular_RNA_topology_v2();
  void add_RNA_basepairs(const map<string,string> &, string ="");
  void add_RNA_basepairs(const vector<ResiduePair> &, string ="");
  void add_RNA_basepairs_v2(const map<string,string> &, string ="");
  void add_RNA_basepairs_v2(const vector<ResiduePair> &, string ="");
  void add_RNA_basestacks(const vector<string> &);
  void add_RNA_basestacks(const vector<ResiduePair> &);
  void add_distance_restraint(string, string, string);
  void add_distance_restraint(DistanceRestraint &);
  void add_distance_restraint(const vector<DistanceRestraint> &);
  size_t number_of_CG_atoms() const;
  void clear();

 private:
  map<int,string> atomid_name;
};

}

#endif