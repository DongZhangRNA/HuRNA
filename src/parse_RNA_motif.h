#ifndef CBS_PARSE_RNA_MOTIF_H
#define CBS_PARSE_RNA_MOTIF_H

#include "RNA_motif.h"
#include "RNA_basepair.h"

namespace CBS_NS {class ParseRNAMotif;}
ostream& operator <<(ostream &, CBS_NS::ParseRNAMotif &);

namespace CBS_NS {

class ParseRNAMotif : protected Pointers {
 public:
  vector<RNAMotif> motifs;
  static const int sl,ti,mw,slfl,tifl,mwfl,coil; //for RNA topology
  
  ParseRNAMotif(Cobmms *);
  ParseRNAMotif(const ParseRNAMotif &);
  ParseRNAMotif& operator=(const ParseRNAMotif &);
  ~ParseRNAMotif();

  void parse_RNA_motifs(const RNABasePair &);  //for the structure prediction
  int get_RNA_topology() const;
  string get_RNA_topology_string() const;
  bool is_flexible_loop(const RNAMotif &) const;
  size_t count_flexible_loop() const;
  size_t count_unpaired_chain() const;
  size_t size() const;
  void clear();
  friend ostream& ::operator <<(ostream &, ParseRNAMotif &);

 private:
  RNABasePair basepair;
  vector<Molecule> helices;
  map<string,string> hx_ends;
  vector<Molecule> loops;

  void get_helices();
  void get_loops();
  void extract_helix();
  void extract_tail();
  void extract_hairpin();
  void extract_kiss();
  void extract_junction();
  void extract_loop();
};

}

#endif