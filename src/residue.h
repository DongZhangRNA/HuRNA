#ifndef CBS_RESIDUE_H
#define CBS_RESIDUE_H

#include "atom.h"
#include "pointers.h"

namespace CBS_NS {class Residue;}
ostream& operator <<(ostream &, CBS_NS::Residue &);

namespace CBS_NS {

class Residue : protected Pointers {
 public:
  int resid;
  string chainid,resname;
  vector<Atom> atoms;

  Residue(Cobmms *);
  Residue(Cobmms *, string, int, string);
  Residue(const Residue &);
  Residue& operator=(const Residue &);
  ~Residue();

  void set_residue(string, int, string);
  void set_residue(string);
  void set_residue(const Residue &);
  void set_residue(const vector<Atom> &);
  //void set_RNA_residue(vector<Atom> &);
  string get_residue_name() const;
  void reset_residue_name(string);
  void reset_residue_name(string,int,string);
  void recovery_RNA_backbone_atoms();
  void simplify_RNA_base();
  void recovery_RNA_base();
  vector<double> get_residue_xyz() const;
  vector<double> get_RNA_residue_xyz();
  void reset_residue_xyz(const vector<double> &);
  void reset_residue_xyz(const vector<Atom> &);
  void set_residue_xyz(const Residue &);
  void set_RNA_residue_xyz(Residue);
  bool check_atom(string) const;
  bool is_atom_exist(string) const;
  Atom find_atom(string) const;
  Atom& operator [](string);
  void delete_atom(string);
  void delete_atom(vector<string> &);
  void delete_atom_hydrogen();
  size_t size() const;
  vector<double> center_of_geometry() const;
  void clear();
  int set_atomid(int);
  int set_RNA_atomid(int,string ="std");
  int set_atomid_empty(int);
  friend ostream& ::operator <<(ostream &, Residue &);

 private:
  map<string,vector<double> > atomxyz;

  void add_simplified_RNA_base_atoms();
  inline bool check_atoms(const vector<string> &);
  inline bool add_atom(string, const vector<double> &);
  inline void remove_atom(string);
};

}

#endif