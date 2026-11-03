#ifndef CBS_EXTRACT_RNA_TEMPLATE_H
#define CBS_EXTRACT_RNA_TEMPLATE_H

#include "RNA_motif.h"
#include "RNA_molecule.h"
#include "RNA_basepair.h"

namespace CBS_NS {class ExtractRNATemplate;}
ostream& operator <<(ostream &, CBS_NS::ExtractRNATemplate &);

namespace CBS_NS {

class ExtractRNATemplate : protected Pointers {
 public:
  vector<RNAMotif> templates;
  vector<RNAMotif> extended_templates;
  vector<string> Vinfo;  //record log info
  
  ExtractRNATemplate(Cobmms *);
  ExtractRNATemplate(const ExtractRNATemplate &);
  ExtractRNATemplate& operator=(const ExtractRNATemplate &);
  ~ExtractRNATemplate();

  void extract_RNA_templates(const RNAMolecule &, const RNABasePair &, int =0);
  size_t size() const;
  void clear();
  friend ostream& ::operator <<(ostream &, ExtractRNATemplate &);

 private:
  RNABasePair *basepair;
  vector<Molecule> helices;
  map<string,string> hx_ends;
  map<string,string> pseudo_basepairs;
  vector<Molecule> loops;

  void get_helices();
  void get_loops();
  void extract_helix();
  void extract_tail();
  void extract_hairpin();
  void extract_loop_template();
  void extract_pseudoknot_template();
  void extract_kiss_template();
  void extract_junction_template();
  vector<RNAMotif> rotate_junction(RNAMotif &);
  inline bool is_pseudo_basepairs(string &);
  inline bool is_pseudo_basepairs(string &, string &);
  void extend_loop_kiss_loop(RNAMotif &);
  void extend_hairpin_kiss(RNAMotif &);
  void extend_junction(RNAMotif &);
};

}

#endif