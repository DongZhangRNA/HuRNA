#ifndef CBS_MATCH_RNA_MOTIF_H
#define CBS_MATCH_RNA_MOTIF_H

#include "RNA_motif.h"

namespace CBS_NS {

class MatchRNAMotif : protected Pointers {
 public:
  int match_score;   //match_score=0, perfect match; =1/2/3/..., one/two/three/... segments are not match 
  //double length_difference;
  double sequence_similarity;
  string match_code;
  string template_info;
  static const int MATCH_NULL;

  MatchRNAMotif(Cobmms *);
  MatchRNAMotif(const MatchRNAMotif &);
  MatchRNAMotif& operator=(const MatchRNAMotif &);
  ~MatchRNAMotif();

  void get_match_score(const RNAMotif &, string &);
  string get_match_info() const;
  void clear();

 private:
  const RNAMotif *motif_ptr;
  RNAMotif *template_ptr;
  
  void match_Loop();
  void match_Hairpin();
  void match_Tail();
  void match_Loop_kiss(); 
  void match_Junction();
  void modify_paired_nucleotide(string &, size_t);
  double calculate_sequence_similarity(const string &, const string &);
  
};

}

#endif