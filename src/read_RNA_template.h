#ifndef CBS_READ_RNA_TEMPLATE_H
#define CBS_READ_RNA_TEMPLATE_H

#include "RNA_motif.h"

namespace CBS_NS {

class ReadRNATemplate : protected Pointers {
 public:
  vector<vector<string> > template_list;
  vector<string> full_template_list;
  string type,template_name;

  ReadRNATemplate(Cobmms *, string);
  ReadRNATemplate(const ReadRNATemplate &);
  ReadRNATemplate& operator=(const ReadRNATemplate &);
  ~ReadRNATemplate();
  
  bool is_template_exist();
  void read_templates();
  void read_full_template_list();
  RNAMotif read_template_full(string) const;
  vector<double> read_template_xyz(string) const;
  //vector<double> read_template_xyz(Molecule &, string);
  void write_templates();
  void write_full_template_list();
  void clear();
};

}

#endif