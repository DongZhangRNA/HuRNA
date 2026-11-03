#pragma once
#ifndef CBS_UPDATE_RNA_TEMPLATE_H
#define CBS_UPDATE_RNA_TEMPLATE_H

#include "read_PDB.h"
#include "extract_RNA_template.h"
#include "read_RNA_template.h"

namespace CBS_NS {

  Cobmms *cbs;
  string _version_;
  vector<string> PDB_list;
  vector<string> Vpdbs;
  vector<vector<string> > VVresnames_old;
  vector<RNAMolecule> Vrnas;
  vector<ExtractRNATemplate> Vtemplates;

  int extract_template(const RNAMolecule &, const RNABasePair &, vector<ExtractRNATemplate> &, int);
  void write_info(string &);
  void check_template_list (string &);
  int read_pdb_xyz(string, ReadPDB &, vector<Molecule> &, int);
  int read_template_xyz(const Molecule &, string, vector<vector<double> > &, int);
  //int read_template_xyz(ReadRNATemplate &, string, vector<vector<double> > &, int);
}

#endif