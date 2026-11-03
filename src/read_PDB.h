#ifndef CBS_READ_PDB_H
#define CBS_READ_PDB_H

#include "atom.h"
#include "pointers.h"

namespace CBS_NS {

class ReadPDB : protected Pointers {
 public:
  vector<vector<Atom> > Full_atoms;
  vector<vector<Atom> > RNA_atoms;
  vector<vector<Atom> > Protein_atoms;
  vector<vector<Atom> > DNA_atoms;
  string PDBid;
  int num_modified;
  string expt_method;
  string resolution;
  
  ReadPDB(Cobmms *);
  ReadPDB(Cobmms *, string);
  ReadPDB(const ReadPDB &);
  ReadPDB& operator=(const ReadPDB &);
  ~ReadPDB();

  void read_PDB(string);
  void parse_PDB(string);
  void parse_PDB(const vector<string> &);
  void parse_CIF(string);
  void extract_RNA_atoms();
  void extract_protein_atoms();
  void extract_DNA_atoms();
  void clear();
  size_t size() const;
  void calculate_NMR_resolution();

 private:
  void check_residue_name();
  string standardize_RNA_atom_name(string);
  string standardize_RNA_nucleotide_name(string);
};

}

#endif