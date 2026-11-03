#include "basic.h"
#include "readfile.h"
#include "initialization.h"
#include "error.h"
#include "superpose.h"
#include "cluster.h"
#include "read_PDB.h"

#include "CifParserBase.h"
#include "CifScannerInt.h"
#include "CifDefs.h"
#include "CifExcept.h"
#include "CifDataInfo.h"
#include "regex.h"
#include "rcsb_math.h"
#include "DICScannerInt.h"
#include "DICParserBase.h"
#include "GenCont.h"
#include "RcsbFile.h"
#include "Exceptions.h"
#include "mapped_vector.h"
#include "RcsbPlatform.h"

using namespace CBS_NS;

const string snull="NULL";

enum{PDB,CIF};

ReadPDB::ReadPDB(Cobmms *cbs) : Pointers(cbs) {
  num_modified=0;
  resolution="None";
}

ReadPDB::ReadPDB(Cobmms *cbs, string fstr) : Pointers(cbs) {
  num_modified=0;
  resolution="None";
  ReadPDB::read_PDB(fstr);
}

ReadPDB::ReadPDB(const ReadPDB &other) : Pointers(other.cbs) {
  PDBid=other.PDBid;
  num_modified=other.num_modified;
  Full_atoms=other.Full_atoms;
  RNA_atoms=other.RNA_atoms;
  Protein_atoms=other.Protein_atoms;
  DNA_atoms=other.DNA_atoms;
  expt_method=other.expt_method;
  resolution=other.resolution;
}

ReadPDB& ReadPDB::operator=(const ReadPDB &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    PDBid=other.PDBid;
    num_modified=other.num_modified;
    Full_atoms=other.Full_atoms;
    RNA_atoms=other.RNA_atoms;
    Protein_atoms=other.Protein_atoms;
    DNA_atoms=other.DNA_atoms;
    expt_method=other.expt_method;
    resolution=other.resolution;
  }

  return *this;
}

ReadPDB::~ReadPDB() {
  PDBid="";
  num_modified=-1;
  expt_method="";
  resolution="";
  Full_atoms.clear();
  Full_atoms.shrink_to_fit();
  RNA_atoms.clear();
  RNA_atoms.shrink_to_fit();
  Protein_atoms.clear();
  Protein_atoms.shrink_to_fit();
  DNA_atoms.clear();
  DNA_atoms.shrink_to_fit();
}

void ReadPDB::read_PDB(string filename) {
  if(!is_file_exist(filename)) {
    string str="read_PDB.cpp: Cannot find file "+filename;
    error->one(str);
  }
  vector<string> stoken=split_string(filename,"/\\");
  string pdbname=stoken.back();
  if(pdbname.find(".pdb")==string::npos && pdbname.find(".cif")==string::npos)  {
    string str="read_PDB.cpp: Unknown file type --> "+pdbname;
    error->one(str);
  }
  string pdbtype=pdbname.substr(pdbname.size()-4);
  int type=-1;
  if(pdbtype==".pdb") type=PDB;
  else if(pdbtype==".cif") type=CIF;
  else {
    string str="read_PDB.cpp: Unknown file type --> "+pdbtype;
    error->one(str);
  }

  PDBid=pdbname.substr(0,pdbname.size()-4);
  for(size_t ii=0;ii<PDBid.size();ii++)
    PDBid[ii]=tolower(PDBid[ii]);

  if(type==PDB) ReadPDB::parse_PDB(filename);
  else if(type==CIF) ReadPDB::parse_CIF(filename);
}

void ReadPDB::parse_PDB(string fstr) {
  Full_atoms.clear();

  vector<string> Vstr=read_file(fstr);
  vector<Atom> atoms;
  for(size_t ii=0;ii<Vstr.size();ii++) {
    vector<string> stoken=split_string(Vstr[ii]);
    if(stoken[0]=="ENDMDL") {
      if(atoms.size()>0) Full_atoms.push_back(atoms);
      atoms.clear();
      continue;
    }
    else if(stoken[0]=="END") break;
    else if(stoken[0].find("ATOM")==string::npos && stoken[0].find("HETATM")==string::npos) continue;
    Atom atom_one(Vstr[ii]);
    string atomname_one=atom_one.atomname;
    if(isdigit(atomname_one[0]) && atomname_one[1]=='H')
      atomname_one=atomname_one.substr(1)+atom_one.atomname.substr(0,1);
    atom_one.atomname=atomname_one;
    if(init->noH && atom_one.atomname[0]=='H') continue;
    atoms.push_back(atom_one);
  }
  if(atoms.size()>0) Full_atoms.push_back(atoms);

  ReadPDB::check_residue_name();
  ReadPDB::extract_RNA_atoms();
  ReadPDB::extract_protein_atoms();
  ReadPDB::extract_DNA_atoms();
}

void ReadPDB::parse_PDB(const vector<string> &Vstr) {
  Full_atoms.clear();

  vector<Atom> atoms;
  for(size_t ii=0;ii<Vstr.size();ii++) {
    vector<string> stoken=split_string(Vstr[ii]);
    if(stoken[0]=="ENDMDL") {
      if(atoms.size()>0) Full_atoms.push_back(atoms);
      atoms.clear();
      continue;
    }
    else if(stoken[0]=="END") break;
    else if(stoken[0].find("ATOM")==string::npos && stoken[0].find("HETATM")==string::npos) continue;
    Atom atom_one(Vstr[ii]);
    string atomname_one=atom_one.atomname;
    if(isdigit(atomname_one[0]) && atomname_one[1]=='H')
      atomname_one=atomname_one.substr(1)+atom_one.atomname.substr(0,1);
    atom_one.atomname=atomname_one;
    if(init->noH && atom_one.atomname[0]=='H') continue;
    atoms.push_back(atom_one);
  }
  if(atoms.size()>0) Full_atoms.push_back(atoms);

  ReadPDB::check_residue_name();
  ReadPDB::extract_RNA_atoms();
  ReadPDB::extract_protein_atoms();
  ReadPDB::extract_DNA_atoms();
}

void ReadPDB::parse_CIF(string cifFileName) {
  Full_atoms.clear();
  string diags;
  CifFile cifFile;

  CifParser *P=new CifParser(&cifFile, cifFile.GetVerbose());
  P->Parse(cifFileName, diags);
  delete(P);

  if(!diags.empty()) {
    string str="read_PDB.cpp: Diagnostics --> \n"+diags;
    error->message(str);
  }

  vector<string> block_names;
  cifFile.GetBlockNames(block_names);

  if(block_names.size()==0) {
    string str="read_PDB.cpp: No block found.";
    error->one(str);
  }

  string general_info_block_name=block_names[0];
  Block& general_info_block=cifFile.GetBlock(general_info_block_name);
  ISTable& entity_table=general_info_block.GetTable("entry");

  vector<string> entity_id_values;
  entity_table.GetColumn(entity_id_values, "id");

  if(block_names.size()==0) {
    string str="read_PDB.cpp: No entry found.";
    error->one(str);
  }
  
  for(size_t mm=0;mm<entity_id_values.size();mm++) {
    string molecule_structure_id=entity_id_values[mm];

    Block& block=cifFile.GetBlock(molecule_structure_id);
    if(block.IsTablePresent("exptl")) {
      ISTable& exptl_table=block.GetTable("exptl");
      vector<string> exptl_column_names=exptl_table.GetColumnNames();
      bool EXPTL_METHOD = (find(exptl_column_names.begin(),exptl_column_names.end(),"method") != exptl_column_names.end()) ? true : false;
      vector<string> V_exptl_method;
      if(EXPTL_METHOD==true) {
        exptl_table.GetColumn(V_exptl_method, "method");
      }
      expt_method=V_exptl_method[0];
    }
    else
      expt_method="None";
    
    vector<string> V_resolution;
    if(expt_method.find("DIFFRACTION")!=string::npos) {
      if(block.IsTablePresent("refine_hist")) {
        ISTable& refine_hist_table=block.GetTable("refine_hist");
        vector<string> refine_hist_column_names=refine_hist_table.GetColumnNames();
        bool RESOLUTION = (find(refine_hist_column_names.begin(),refine_hist_column_names.end(),"d_res_high") != refine_hist_column_names.end()) ? true : false;
        if(RESOLUTION==true) {
          refine_hist_table.GetColumn(V_resolution, "d_res_high");
          resolution=V_resolution[0];
        }
      }
      if((resolution=="." || resolution=="?" || resolution=="None" || resolution=="") 
          && block.IsTablePresent("reflns")) {
        ISTable& reflns_table=block.GetTable("reflns");
        vector<string> reflns_column_names=reflns_table.GetColumnNames();
        bool RESOLUTION_HIGH = (find(reflns_column_names.begin(),reflns_column_names.end(),"d_resolution_high") != reflns_column_names.end()) ? true : false;
        if(RESOLUTION_HIGH==true) {
          reflns_table.GetColumn(V_resolution, "d_resolution_high");
          resolution=V_resolution[0];
        }
      }
    }
    else if(expt_method.find("MICROSCOPY")!=string::npos) {
      if(block.IsTablePresent("em_3d_reconstruction")) {
        ISTable& reconstruction_table=block.GetTable("em_3d_reconstruction");
        vector<string> reconstruction_column_names=reconstruction_table.GetColumnNames();
        bool RESOLUTION = (find(reconstruction_column_names.begin(),reconstruction_column_names.end(),"resolution") != reconstruction_column_names.end()) ? true : false;
        if(RESOLUTION==true) {
          reconstruction_table.GetColumn(V_resolution, "resolution");
          resolution=V_resolution[0];
        }
      }
      if((resolution=="." || resolution=="?" || resolution=="None" || resolution=="")
          && block.IsTablePresent("refine_hist")) {
        ISTable& refine_hist_table=block.GetTable("refine_hist");
        vector<string> refine_hist_column_names=refine_hist_table.GetColumnNames();
        bool RESOLUTION = (find(refine_hist_column_names.begin(),refine_hist_column_names.end(),"d_res_high") != refine_hist_column_names.end()) ? true : false;
        if(RESOLUTION==true) {
          refine_hist_table.GetColumn(V_resolution, "d_res_high");
          resolution=V_resolution[0];
        }
      }
    }
    else
      resolution="None";
    
    if(!block.IsTablePresent("atom_site")) {
      string str="read_PDB.cpp: parse CIF file error --> "+cifFileName;
      error->warning(str);
      return; 
    }

    ISTable& atom_site_table=block.GetTable("atom_site");
    vector<string> column_names=atom_site_table.GetColumnNames();

    bool GROUP_PDB = (find(column_names.begin(),column_names.end(),"group_PDB") != column_names.end()) ? true : false;
    bool ID = (find(column_names.begin(),column_names.end(),"id") != column_names.end()) ? true : false;
    bool TYPE_SYMBOL = (find(column_names.begin(),column_names.end(),"type_symbol") != column_names.end()) ? true : false;
    bool LABEL_ATOM_ID = (find(column_names.begin(),column_names.end(),"label_atom_id") != column_names.end()) ? true : false;
    bool LABEL_ALT_ID = (find(column_names.begin(),column_names.end(),"label_alt_id") != column_names.end()) ? true : false;
    bool LABEL_COMP_ID = (find(column_names.begin(),column_names.end(),"label_comp_id") != column_names.end()) ? true : false;
    bool LABEL_ASYM_ID = (find(column_names.begin(),column_names.end(),"label_asym_id") != column_names.end()) ? true : false;
    bool LABEL_ENTITY_ID = (find(column_names.begin(),column_names.end(),"label_entity_id") != column_names.end()) ? true : false;
    bool LABEL_SEQ_ID = (find(column_names.begin(),column_names.end(),"label_seq_id") != column_names.end()) ? true : false;
    bool PDBX_PDB_INS_CODE = (find(column_names.begin(),column_names.end(),"pdbx_PDB_ins_code") != column_names.end()) ? true : false;
    bool CARTN_X = (find(column_names.begin(),column_names.end(),"Cartn_x") != column_names.end()) ? true : false;
    bool CARTN_Y = (find(column_names.begin(),column_names.end(),"Cartn_y") != column_names.end()) ? true : false;
    bool CARTN_Z = (find(column_names.begin(),column_names.end(),"Cartn_z") != column_names.end()) ? true : false;
    bool OCCUPANCY = (find(column_names.begin(),column_names.end(),"occupancy") != column_names.end()) ? true : false;
    bool B_ISO_OR_EQUIV = (find(column_names.begin(),column_names.end(),"B_iso_or_equiv") != column_names.end()) ? true : false;
    bool CARTN_X_ESD = (find(column_names.begin(),column_names.end(),"Cartn_x_esd") != column_names.end()) ? true : false;
    bool CARTN_Y_ESD = (find(column_names.begin(),column_names.end(),"Cartn_y_esd") != column_names.end()) ? true : false;
    bool CARTN_Z_ESD = (find(column_names.begin(),column_names.end(),"Cartn_z_esd") != column_names.end()) ? true : false;
    bool OCCUPANCY_ESD = (find(column_names.begin(),column_names.end(),"occupancy_esd") != column_names.end()) ? true : false;
    bool B_ISO_OR_EQUIV_ESD = (find(column_names.begin(),column_names.end(),"B_iso_or_equiv_esd") != column_names.end()) ? true : false;
    bool PDBX_FORMAL_CHARGE = (find(column_names.begin(),column_names.end(),"pdbx_formal_charge") != column_names.end()) ? true : false;
    bool AUTH_SEQ_ID = (find(column_names.begin(),column_names.end(),"auth_seq_id") != column_names.end()) ? true : false;
    bool AUTH_COMP_ID = (find(column_names.begin(),column_names.end(),"auth_comp_id") != column_names.end()) ? true : false;
    bool AUTH_ASYM_ID = (find(column_names.begin(),column_names.end(),"auth_asym_id") != column_names.end()) ? true : false;
    bool AUTH_ATOM_ID = (find(column_names.begin(),column_names.end(),"auth_atom_id") != column_names.end()) ? true : false;
    bool PDBX_PDB_MODEL_NUM = (find(column_names.begin(),column_names.end(),"pdbx_PDB_model_num") != column_names.end()) ? true : false;

    vector<string> V_group_PDB;
    if(GROUP_PDB==true) {atom_site_table.GetColumn(V_group_PDB, "group_PDB");}
    vector<string> V_id;
    if(ID==true) {atom_site_table.GetColumn(V_id, "id");}
    vector<string> V_type_symbol;
    if(TYPE_SYMBOL==true) {atom_site_table.GetColumn(V_type_symbol, "type_symbol");}
    vector<string> V_label_atom_id;
    if(LABEL_ATOM_ID==true) {atom_site_table.GetColumn(V_label_atom_id, "label_atom_id");}
    vector<string> V_label_alt_id;
    if(LABEL_ALT_ID==true) {atom_site_table.GetColumn(V_label_alt_id, "label_alt_id");}
    vector<string> V_label_comp_id;
    if(LABEL_COMP_ID==true) {atom_site_table.GetColumn(V_label_comp_id, "label_comp_id");}
    vector<string> V_label_asym_id;
    if(LABEL_ASYM_ID==true) {atom_site_table.GetColumn(V_label_asym_id, "label_asym_id");}
    vector<string> V_label_entity_id;
    if(LABEL_ENTITY_ID==true) {atom_site_table.GetColumn(V_label_entity_id, "label_entity_id");}
    vector<string> V_label_seq_id;
    if(LABEL_SEQ_ID==true) {atom_site_table.GetColumn(V_label_seq_id, "label_seq_id");}
    vector<string> V_pdbx_PDB_ins_code;
    if(PDBX_PDB_INS_CODE==true) {atom_site_table.GetColumn(V_pdbx_PDB_ins_code, "pdbx_PDB_ins_code");}
    vector<string> V_Cartn_x;
    if(CARTN_X==true) {atom_site_table.GetColumn(V_Cartn_x, "Cartn_x");}
    vector<string> V_Cartn_y;
    if(CARTN_Y==true) {atom_site_table.GetColumn(V_Cartn_y, "Cartn_y");}
    vector<string> V_Cartn_z;
    if(CARTN_Z==true) {atom_site_table.GetColumn(V_Cartn_z, "Cartn_z");}
    vector<string> V_occupancy;
    if(OCCUPANCY==true) {atom_site_table.GetColumn(V_occupancy, "occupancy");}
    vector<string> V_B_iso_or_equiv;
    if(B_ISO_OR_EQUIV==true) {atom_site_table.GetColumn(V_B_iso_or_equiv, "B_iso_or_equiv");}
    vector<string> V_Cartn_x_esd;
    if(CARTN_X_ESD==true) {atom_site_table.GetColumn(V_Cartn_x_esd, "Cartn_x_esd");}
    vector<string> V_Cartn_y_esd;
    if(CARTN_Y_ESD==true) {atom_site_table.GetColumn(V_Cartn_y_esd, "Cartn_y_esd");}
    vector<string> V_Cartn_z_esd;
    if(CARTN_Z_ESD==true) {atom_site_table.GetColumn(V_Cartn_z_esd, "Cartn_z_esd");}
    vector<string> V_occupancy_esd;
    if(OCCUPANCY_ESD==true) {atom_site_table.GetColumn(V_occupancy_esd, "occupancy_esd");}
    vector<string> V_B_iso_or_equiv_esd;
    if(B_ISO_OR_EQUIV_ESD==true) {atom_site_table.GetColumn(V_B_iso_or_equiv_esd, "B_iso_or_equiv_esd");}
    vector<string> V_pdbx_formal_charge;
    if(PDBX_FORMAL_CHARGE==true) {atom_site_table.GetColumn(V_pdbx_formal_charge, "pdbx_formal_charge");}
    vector<string> V_auth_seq_id;
    if(AUTH_SEQ_ID==true) {atom_site_table.GetColumn(V_auth_seq_id, "auth_seq_id");}
    vector<string> V_auth_comp_id;
    if(AUTH_COMP_ID==true) {atom_site_table.GetColumn(V_auth_comp_id, "auth_comp_id");}
    vector<string> V_auth_asym_id;
    if(AUTH_ASYM_ID==true) {atom_site_table.GetColumn(V_auth_asym_id, "auth_asym_id");}
    vector<string> V_auth_atom_id;
    if(AUTH_ATOM_ID==true) {atom_site_table.GetColumn(V_auth_atom_id, "auth_atom_id");}
    vector<string> V_pdbx_PDB_model_num;
    if(PDBX_PDB_MODEL_NUM==true) {atom_site_table.GetColumn(V_pdbx_PDB_model_num, "pdbx_PDB_model_num");}

    int Nmodel=(int)mm+1;
    vector<Atom> atoms;
    for(size_t atom_size=0; atom_size!=V_group_PDB.size();++atom_size ) {
      int modelSerial=0;
      string chainID="";
      int resSeq=0;
      string atomName="";
      string resName="";
      int atomSerial=0;
      double x=0.0;
      double y=0.0;
      double z=0.0;
      string altID="";
      string atom_type_symbol="";
      string insertion_code="";

      if(ID==true && TYPE_SYMBOL==true && LABEL_ATOM_ID==true && LABEL_ALT_ID==true &&
        LABEL_COMP_ID==true && LABEL_ASYM_ID==true && LABEL_SEQ_ID==true && CARTN_X==true &&
        CARTN_Y==true && CARTN_Z==true && PDBX_PDB_MODEL_NUM==true && PDBX_PDB_INS_CODE==true ) 
      {
        atomSerial=(V_id[atom_size] == "?" || V_id[atom_size] == ".") ? -1 : stoi(V_id[atom_size]);
        atom_type_symbol=V_type_symbol[atom_size];
        atomName=V_label_atom_id[atom_size];
        altID=V_label_alt_id[atom_size];
        resName=V_label_comp_id[atom_size];
        if(AUTH_ASYM_ID==true) chainID=V_auth_asym_id[atom_size];
        else chainID=V_label_asym_id[atom_size];
        if(AUTH_SEQ_ID==true) resSeq=stoi(V_auth_seq_id[atom_size]);
        else resSeq=(V_label_seq_id[atom_size] == "?" || V_label_seq_id[atom_size] == ".") ? -1 : stoi(V_label_seq_id[atom_size]);
        insertion_code=V_pdbx_PDB_ins_code[atom_size];
        x=stod(V_Cartn_x[atom_size]);
        y=stod(V_Cartn_y[atom_size]);
        z=stod(V_Cartn_z[atom_size]);
        modelSerial=(V_pdbx_PDB_model_num[atom_size] == "?" || V_pdbx_PDB_model_num[atom_size] == ".") ? -1 : stoi(V_pdbx_PDB_model_num[atom_size]);
      }
      else {
        string str="read_PDB.cpp: parse CIF file error --> "+cifFileName;
        error->warning(str);
        return; 
      }

      Atom atom_one;
      atom_one.set_atom(atomSerial,chainID,resSeq,resName,atomName,x,y,z);
      if(altID=="?" || altID==".") altID="";
      atom_one.atomalt=altID;
      if(insertion_code=="?" || insertion_code==".") insertion_code=" ";
      if(insertion_code.size()>0)
        atom_one.iCode=insertion_code[0];
      if(modelSerial!=Nmodel) {
        Full_atoms.push_back(atoms);
        atoms.clear();
        Nmodel=modelSerial;
      }
      string atomname_one=atom_one.atomname;
      if(isdigit(atomname_one[0]) && atomname_one[1]=='H')
        atomname_one=atomname_one.substr(1)+atom_one.atomname.substr(0,1);
      atom_one.atomname=atomname_one;
      if(init->noH && atom_one.atomname[0]=='H') continue;
      atoms.push_back(atom_one);
    }
    if(atoms.size()>0) Full_atoms.push_back(atoms);
  }

  ReadPDB::check_residue_name();
  ReadPDB::extract_RNA_atoms();
  ReadPDB::extract_protein_atoms();
  ReadPDB::extract_DNA_atoms();
}

void ReadPDB::extract_RNA_atoms() {
  RNA_atoms.clear();
  
  for(size_t ii=0;ii<Full_atoms.size();ii++) {
    vector<Atom> atoms;
    for(size_t jj=0;jj<Full_atoms[ii].size();jj++) {
      Atom atom_one=Full_atoms[ii][jj];
      atom_one.resname=ReadPDB::standardize_RNA_nucleotide_name(atom_one.resname);
      atom_one.atomname=ReadPDB::standardize_RNA_atom_name(atom_one.atomname);
      if(atom_one.resname==snull || atom_one.atomname==snull) continue;
      atoms.push_back(atom_one);
    }
    if(atoms.size()>0) RNA_atoms.push_back(atoms);
  }
}

void ReadPDB::extract_protein_atoms() {
  Protein_atoms.clear();
  
  for(size_t ii=0;ii<Full_atoms.size();ii++) {
    vector<Atom> atoms;
    for(size_t jj=0;jj<Full_atoms[ii].size();jj++) {
      Atom atom_one=Full_atoms[ii][jj];
    }
    if(atoms.size()>0) Protein_atoms.push_back(atoms);
  }
}

void ReadPDB::extract_DNA_atoms() {
  DNA_atoms.clear();
  
  for(size_t ii=0;ii<Full_atoms.size();ii++) {
    vector<Atom> atoms;
    for(size_t jj=0;jj<Full_atoms[ii].size();jj++) {
      Atom atom_one=Full_atoms[ii][jj];
    }
    if(atoms.size()>0) DNA_atoms.push_back(atoms);
  }
}

void ReadPDB::clear() {
  PDBid="";
  num_modified=-1;
  Full_atoms.clear();
  Full_atoms.shrink_to_fit();
  RNA_atoms.clear();
  RNA_atoms.shrink_to_fit();
  Protein_atoms.clear();
  Protein_atoms.shrink_to_fit();
  DNA_atoms.clear();
  DNA_atoms.shrink_to_fit();
  expt_method="";
  resolution="None";
}

size_t ReadPDB::size() const {
  return Full_atoms.size();
}

void ReadPDB::calculate_NMR_resolution() {
  if(expt_method.find("NMR")==string::npos) return;

  if(Full_atoms.size()==1) {
    resolution="None";
    return;
  }

  Cluster cluset_one(cbs);
  cluset_one.input_data(Full_atoms);
  cluset_one.calculate_RMSD_matrix();
  cluset_one.get_RMSD_cutoff(1.0);
  double rmsd_cutoff=cluset_one.RMSD_cutoff;

  ostringstream oss;
  oss<<std::fixed<<std::setprecision(2)<<rmsd_cutoff;
  resolution=oss.str();
  //resolution+=" ("+to_string(Full_atoms.size())+")";
}

string ReadPDB::standardize_RNA_atom_name(string name_one) {
  vector<string> stoken=split_string(name_one);
  string name=stoken[0];
  
  while(name.find("*")!=string::npos)
    name.replace(name.find_first_of("*"),1,"'");
  if(name.find("O")!=string::npos && name.find("P")!=string::npos && name.find("1")!=string::npos)
    name="OP1";
  if(name.find("O")!=string::npos && name.find("P")!=string::npos && name.find("2")!=string::npos)
    name="OP2";
  if(count(init->RESIDUES_ATOMS["RNA"].begin(),init->RESIDUES_ATOMS["RNA"].end(),name)==0) {
    name=snull;
    //string str="read_PDB.cpp: Unknown RNA atom name --> "+name_one;
    //error->warning(str);
  }

  return name;
}

void ReadPDB::check_residue_name() {
  vector<string> Vresnames_one;
  for(size_t ii=0;ii<Full_atoms.size();ii++) {
    for(size_t jj=0;jj<Full_atoms[ii].size();jj++) {
      string resname_one=Full_atoms[ii][jj].resname;
      if(count(Vresnames_one.begin(),Vresnames_one.end(),resname_one)==0)
        Vresnames_one.push_back(resname_one);
    }
  }
  
  vector<string> Vsol_names=init->read_water_ions_name();
  for(size_t ii=0;ii<Vresnames_one.size();ii++) {
    string resname_one=Vresnames_one[ii];
    if(count(Vsol_names.begin(),Vsol_names.end(),resname_one)) continue;
    int flag_nt=0;
    if(init->NUCLEOTIDES_LIST.count(resname_one)) flag_nt=1;
    int flag_aa=0;
    if(init->RESIDUES_LIST.count(resname_one)) flag_aa=1;

    if(flag_aa==0 && flag_nt==0) {
      string str="read_PDB.cpp: Unknown residue name --> "+resname_one;
      error->warning(str);
    }
  }
}

string ReadPDB::standardize_RNA_nucleotide_name(string seq_one) {
  string seq="";
  num_modified=0;
  vector<string> stoken=split_string(seq_one);
  //for(size_t ii=0;ii<stoken[0].size();ii++) {
    //char ch=toupper(stoken[0][ii]);
    //seq+=to_string(ch);
  //}
  seq=stoken[0];
  if(init->NUCLEOTIDES_LIST.count(seq)) {
    if(count(init->RNA_STANDARD_NUCLEOTIDES.begin(),init->RNA_STANDARD_NUCLEOTIDES.end(),seq)==0)
      num_modified++;
    seq=init->NUCLEOTIDES_LIST[seq];
  }
  if(count(init->RNA_STANDARD_NUCLEOTIDES.begin(),init->RNA_STANDARD_NUCLEOTIDES.end(),seq)==0) {
    seq=snull;
    //string str="read_PDB.cpp: Unknown RNA nucleotide name --> "+seq_one;
    //error->warning(str);
  }

  return seq;
}