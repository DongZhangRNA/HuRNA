#ifndef CBS_MOLECULE_H
#define CBS_MOLECULE_H

#include "chain.h"

namespace CBS_NS {class Molecule;}
ostream& operator <<(ostream &, CBS_NS::Molecule);

namespace CBS_NS {

class Molecule : protected Pointers {
 public:
  vector<vector<string> > chain_resnames; //only for sequences
  vector<Chain> chains;
  vector<string> chain_ids;
  vector<int> chain_labels;
  vector<string> Vresnames;  // for full resname: "chainid/resid/resname"

  explicit Molecule(Cobmms *);
  explicit Molecule(Cobmms *, vector<vector<string> > &, int);
  Molecule(const Molecule &);
  Molecule& operator=(const Molecule &);
  ~Molecule();

  void set_molecule(const vector<vector<string> > &, int);
  void set_molecule(const vector<Atom> &, int);
  void set_molecule(const vector<Chain> &);
  void set_molecule(const Molecule &);
  Molecule copy_molecule() const;
  void standardize_molecule(int =0);
  void simplify_sidechain();
  void recovery_sidechain();
  void replace_sequence(vector<vector<string> > &);
  void update_residue_name();
  vector<vector<string> > get_residue_name() const;
  void change_chain_name(vector<vector<string> > &, vector<vector<string> > &);
  vector<double> get_molecule_xyz();
  vector<double> get_molecule_xyz() const;
  void reset_molecule_xyz(const vector<double> &);
  void reset_molecule_xyz(const vector<Atom> &);
  void set_molecule_xyz(Molecule &);
  void move_molecule(double, double, double);
  void move_molecule(const vector<double> &);
  string next_residue(string) const;
  string previous_residue(string) const;
  Residue find_residue(string) const;
  vector<Atom> get_atom_list() const;
  Atom find_atom(string) const;
  int find_atomid(string) const;
  string get_atomid_list(string) const;
  string get_atomid_list(string, string) const;
  bool is_residue_exist(string) const;
  int get_residue_index(string, size_t =0) const;
  void insert_chain(const Chain &);
  void clear();
  size_t size() const;
  vector<string> get_residue_list(string,string) const;
  size_t number_of_residues() const;
  //size_t number_of_residues(string) const;
  size_t number_of_residues(string,string) const;
  size_t number_of_atoms() const;
  vector<double> center_of_geometry() const;
  string get_sequence() const;
  Chain& operator [](int);
  void set_atomid(int =0);
  friend ostream& ::operator <<(ostream &, Molecule);

 protected:
  vector<int> Vchainid;
};

}

#endif