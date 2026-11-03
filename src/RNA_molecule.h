#ifndef CBS_RNA_MOLECULE_H
#define CBS_RNA_MOLECULE_H

#include "struct_definition.h"
#include "molecule.h"

namespace CBS_NS {

class RNAMolecule : public Molecule {
 public:
  vector<string> sequences;
  
  vector<InteractionPair> basepairs_WC;
  vector<InteractionPair> basepairs_nWC;
  vector<InteractionPair> basestacks;

  RNAMolecule(Cobmms *);
  RNAMolecule(Cobmms *, string);
  RNAMolecule(const RNAMolecule &);
  RNAMolecule& operator=(const RNAMolecule &);
  ~RNAMolecule();

  void set_RNA_molecule(const RNAMolecule &);
  void set_RNA_molecule(string, int =0);
  void set_RNA_molecule(const vector<Atom> &);
  string get_sequences() const;
  void replace_sequence(string &);
  string extract_secondary_structure() const;
  void extract_interactions();
  string extract_secondary_structure_MCAnnotate() const;
  void extract_interactions_MCAnnotate();
  void clear();
  bool is_basepaired(string, string) const;
  bool is_basepaired(Residue, Residue) const;
  static void swap_interactionpair_type(InteractionPair &);
  vector<string> read_sequence(string);
  vector<string> get_resnames_from_sequence(string);

 protected:
  inline bool is_WC_basepair(string) const;
  inline string change_dssr_residue_name(string);
};

}

#endif