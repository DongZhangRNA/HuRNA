#ifndef CBS_ATOM_H
#define CBS_ATOM_H

namespace CBS_NS {class Atom;}
ostream& operator <<(ostream &, CBS_NS::Atom &);

namespace CBS_NS {

class Atom {
 public:
  int id,resid;
  string chainid,resname,atomname,atomalt;
  char iCode;
  double q;
  double x,y,z;

  Atom();
  Atom(string);
  Atom(const Atom &);
  ~Atom();

  void set_atom(int, string, int, string, string, double, double, double);
  void set_atom(string, double, double, double);
  void set_atom(size_t, string, int, string, string, vector<double> &);
  void set_atom(string);
  void set_atom(const Atom &);
  string get_residue_name() const;
  string get_atom_name() const;
  void reset_residue_name(string);
  void reset_residue_name(string,int,string);
  void reset_atom_xyz(const vector<double> &);
  void reset_atom_xyz(const Atom &);
  void set_atom_xyz(const Atom &);
  void move_atom(double, double, double);
  void move_atom(const vector<double> &);
  vector<double> get_atom_xyz() const;
  double distance(const Atom &);
  void clear();
  friend ostream& ::operator <<(ostream &, Atom &);
};

}

#endif