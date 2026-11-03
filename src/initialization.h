#ifndef CBS_INITIALIZATION_H
#define CBS_INITIALIZATION_H

#include "atom.h"
#include "pointers.h"

namespace CBS_NS {

class Initialization : protected Pointers {
 public:
  int noH;
  double RNA_distance_O3p_P, RNA_distance_O3p_P_min, RNA_distance_O3p_P_max;
  int RNA_check_chain_break,RNA_use_DM;
  string HOME_DIR, INPUT_DIR, RNA_TEMPLATE_DIR;
  string WORK_DIR;
  string DSSR,MCAnnotate,LMP,RNA_DM,Contrafold,USalign;
  string MPIRUN_OPT;
  double RNA_BP_RMSD_cutoff, RNA_template_RMSD_cutoff;
  double RNA_sequence_similarity_high, RNA_sequence_similarity_medium;
  int RNA_template_nmax, Num_model_max;
  int user_template_rigid;
  double assembly_delta_distance;
  int nthread;
  int Num_remd_dup;
  string output_format;
  string RNA_version;

  vector<string> Excluded_PDBs;
  vector<string> CircRNA_chains;
  
  map<string,vector<Atom> > RESIDUES_XYZ;
  map<string,vector<string> > RESIDUES_ATOMS;
  map<string,string> NUCLEOTIDES_LIST;
  map<string,string> RESIDUES_LIST;
  vector<string> RNA_STANDARD_NUCLEOTIDES;
  vector<string> PROTEIN_STANDARD_NUCLEOTIDES;

  map<string,string> SEQUENCE_MAP;

  vector<string> CHAIN_ID_LIST;

  vector<char> RNA_BP_NOTATION;

  Initialization(Cobmms *);
  ~Initialization();

  void set_parameters();
  void reset_parameters(string, string);
  vector<string> read_water_ions_name();
  map<string,string> read_mutation_atommap();
  void load_RNA_DM(int);

 private:
  void set_residues_list();
  void set_standard_residues();
  void set_residue_name_map();
  void set_RNA_basepair_notation();
};

}

#endif