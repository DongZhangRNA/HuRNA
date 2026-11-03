#ifndef CBS_FIND_RNA_TEMPLATE_H
#define CBS_FIND_RNA_TEMPLATE_H

#include "match_RNA_motif.h"
#include "read_RNA_template.h"

namespace CBS_NS {

class FindRNATemplate : protected Pointers {
 public:
  vector<MatchRNAMotif> Vtemplate_scores;
  vector<RNAMotif> Vtemplates;

  FindRNATemplate(Cobmms *);
  FindRNATemplate(const FindRNATemplate &);
  FindRNATemplate& operator=(const FindRNATemplate &);
  ~FindRNATemplate();

  void find_template(const RNAMotif &);
  void find_template_DM(const RNAMotif &);
  void find_Hairpin(const RNAMotif &);
  void find_Tail(const RNAMotif &);
  void find_Hairpin_kiss(const RNAMotif &);
  void find_Junction(const RNAMotif &);
  void find_Junction_kiss(const RNAMotif &);
  void find_Loop(const RNAMotif &);
  void find_Loop(const RNAMotif &, double);
  void find_Loop_kiss_loop(const RNAMotif &);
  void find_Loop_kiss_loop(const RNAMotif &, double);
  void clear();

 private:
  void select_template(vector<MatchRNAMotif> &, ReadRNATemplate &);
};

}

#endif