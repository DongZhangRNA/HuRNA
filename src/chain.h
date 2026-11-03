#ifndef CBS_CHAIN_H
#define CBS_CHAIN_H

#include "residue.h"

namespace CBS_NS {class Chain;}
ostream& operator <<(ostream &, CBS_NS::Chain &);

namespace CBS_NS {

class Chain : protected Pointers {
 public:
  string chainid;
  vector<string> resnames; //only for sequence
  vector<Residue> residues;
  int chain_label;
  static const int begin,end;
  static const int rna,dna,protein,ligand,cg,other,empty;

  explicit Chain(Cobmms *);
  explicit Chain(Cobmms *, vector<string> &, int, string ="A", int =0);
  Chain(const Chain &);
  Chain& operator=(const Chain &);
  ~Chain();

  void set_chain(const vector<string> &, int, string, int);
  void set_chain(const vector<string> &, int);
  void set_chain(const vector<Residue> &, int);
  void set_chain(const Chain &);
  void standardize_chain();
  void standardize_chain(string, int =0);
  void simplify_sidechain();
  void recovery_sidechain();
  void replace_sequence(vector<string> &);
  vector<string> get_residue_name() const;
  void reset_chain_name(string);
  void reset_residues_name(string,int);
  void change_chain_name(vector<string> &, vector<string> &);
  vector<double> get_chain_xyz();
  vector<double> get_chain_xyz() const;
  void reset_chain_xyz(const vector<double> &);
  void set_chain_xyz(Chain &);
  void move_chain(double, double, double);
  void move_chain(const vector<double> &);
  void insert_residue(Residue &, int);
  void insert_residue(string, int);
  void delete_residue(int);
  void delete_residue(string);
  string first_of_residue() const;
  string last_of_residue() const;
  void clear();
  size_t size() const;
  size_t number_of_atoms() const;
  vector<double> center_of_geometry() const;
  Residue& operator [](int);
  int set_atomid(int);
  string get_sequence() const;
  friend ostream& ::operator <<(ostream &, Chain &);
};

}

#endif