#ifndef CBS_RNA_BASEPAIR_H
#define CBS_RNA_BASEPAIR_H

#include "pointers.h"

namespace CBS_NS {class RNABasePair;}
ostream& operator <<(ostream &, CBS_NS::RNABasePair &);

namespace CBS_NS {

class RNABasePair : protected Pointers {
 public:
  string secstr;  //without chain breaks
  string secstr_out;  //with chain breaks;
  map<string, string> mapBP;
  vector<string> Vresnames;

  RNABasePair(Cobmms *);
  RNABasePair(Cobmms *, const vector<vector<string> > &, string);
  RNABasePair(const RNABasePair &);
  RNABasePair& operator=(const RNABasePair &);
  ~RNABasePair();

  void set_RNA_basepairs(const RNABasePair &);
  void parse_RNA_basepairs(const vector<vector<string> > &, string);
  void set_RNA_basepairs(const vector<vector<string> > &, const map<string,string> &);
  void add_RNA_basepair(string, string);
  bool is_basepaired(string, string) const;
  bool is_basepaired(string) const;
  size_t size() const;
  void clear();
  size_t get_residue_index(string, size_t =0) const;
  int get_basepair_level(string, string) const;
  string get_segment_secondary_structure(string, string) const;
  string previous_residue(string) const;
  string next_residue(string) const;
  bool check_noncanonical_basepair(string, string) const;
  friend ostream& ::operator <<(ostream &, RNABasePair &);

 protected:
  vector<int> chainlength;
  map<string,vector<size_t> > Mchainid_indexs; //index_begin, index_end
  map<string,int> BPlevel;

  void set_basepair_level();
};

}

#endif