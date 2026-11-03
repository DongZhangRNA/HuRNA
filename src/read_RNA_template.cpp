#include "basic.h"
#include "readfile.h"
#include "initialization.h"
#include "error.h"
#include "read_PDB.h"
#include "read_RNA_template.h"

using namespace CBS_NS;

ReadRNATemplate::ReadRNATemplate(Cobmms *cbs, string name_one) : Pointers(cbs) {
  template_name=name_one;
  vector<string> stoken=split_string(name_one,"#");
  type=stoken[0];
}

ReadRNATemplate::ReadRNATemplate(const ReadRNATemplate &other) : Pointers(other.cbs) {
  type=other.type;
  template_name=other.template_name;
  template_list=other.template_list;
  full_template_list=other.full_template_list;
}

ReadRNATemplate& ReadRNATemplate::operator=(const ReadRNATemplate &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    type=other.type;
    template_name=other.template_name;
    template_list=other.template_list;
    full_template_list=other.full_template_list;
  }

  return *this;

}

ReadRNATemplate::~ReadRNATemplate() {
  type="";
  template_name="";
  template_list.clear();
  template_list.shrink_to_fit();
  full_template_list.clear();
  full_template_list.shrink_to_fit();
}

bool ReadRNATemplate::is_template_exist() {
  string fstr=init->RNA_TEMPLATE_DIR+type+"/"+template_name+".dat";
  return is_file_exist(fstr);
}

void ReadRNATemplate::read_templates() {
  template_list.clear();
  if(!ReadRNATemplate::is_template_exist()) return;

  string fstr=init->RNA_TEMPLATE_DIR+type+"/"+template_name+".dat";
  vector<string> Vlines=read_file(fstr);
  vector<string> Vclu;
  for(size_t ii=0;ii<Vlines.size();ii++) {
    if(Vlines.size()<2) continue;
    if(Vlines[ii].find("Template")!=string::npos) {
      if(Vclu.size()>0) template_list.push_back(Vclu);
      Vclu.clear();
      continue;
    }
    Vclu.push_back(Vlines[ii]);
  }
  if(Vclu.size()>0) template_list.push_back(Vclu);
}

void ReadRNATemplate::read_full_template_list() {
  full_template_list.clear();
  if(!ReadRNATemplate::is_template_exist()) return;
  string fstr=init->RNA_TEMPLATE_DIR+"template_list/"+template_name+"_list.dat";
  vector<string> Vlines=read_file(fstr);
  for(auto sline : Vlines) {
    if(sline.size()<2) continue;
    full_template_list.push_back(sline);
  }
}

RNAMotif ReadRNATemplate::read_template_full(string temp_info) const {
  vector<string> stoken=split_string(temp_info);
  string PDBid=stoken[0];
  string type_one=stoken[1];
  string str_size=stoken[2];
  string str_end=stoken[3];
  string str_seq=stoken[4];
  string str_secstr=stoken[5];
  vector<string> Vends;
  stoken=split_string(str_end,"-&");
  for(size_t ii=0;ii<stoken.size();ii++)
    Vends.push_back(stoken[ii]);
  vector<int> lengths;
  stoken=split_string(str_size,"-&");
  for(size_t ii=0;ii<stoken.size();ii++)
    lengths.push_back(stoi(stoken[ii]));
  if(type_one=="Hairpin_kiss" || type_one=="Loop_kiss_loop")
    lengths[0]=lengths[0]+lengths[1]+lengths[2];
  else if(type=="Tail5" || type=="Tail3")
    lengths[0]=lengths[0]-1;

  RNAMotif motif_one(cbs);
  motif_one.type=type;
  stoken=split_string(str_secstr,"-&");
  for(size_t ii=0;ii<stoken.size();ii++)
    motif_one.secstrs.push_back(stoken[ii]);
  
  string fstr=init->RNA_TEMPLATE_DIR+"pdb/"+PDBid+".pdb";
  ReadPDB readpdb_one(cbs,fstr);
  Molecule mole_one(cbs);
  mole_one.set_molecule(readpdb_one.RNA_atoms[0],Chain::rna);
  size_t nchain=Vends.size()/2;
  for(size_t ii=0;ii<nchain;ii++) {
    string residue_begin=Vends[ii*2];
    string residue_end=Vends[ii*2+1];
    vector<string> stoken=split_string(residue_begin,"/");
    string chainid_one=stoken[0];
    Chain loop_i(cbs);
    int check=0;
    for(size_t jj=0;jj<mole_one.size();jj++) {
      for(int kk=0;kk<(int)mole_one.chains[jj].size();kk++) {
        if(mole_one.chains[jj].chainid!=chainid_one) continue;
        if(mole_one.chains[jj][kk].get_residue_name()==residue_begin) {
          if(kk+lengths[ii]+1>=(int)mole_one.chains[jj].size()) break;
          else if(mole_one.chains[jj][kk+lengths[ii]+1].get_residue_name()!=residue_end) continue;
          else {
            loop_i.insert_residue(mole_one.chains[jj][kk],Chain::end);
            check=1;
            continue;
          } 
        }
        if(check==1 && mole_one.chains[jj][kk].get_residue_name()==residue_end) {
          loop_i.insert_residue(mole_one.chains[jj][kk],Chain::end);
          check=2;
          break;
        }
        if(check==1) loop_i.insert_residue(mole_one.chains[jj][kk],Chain::end);
      }
      if(check==2) break;
    }
    loop_i.chain_label=Chain::rna;
    motif_one.insert_RNA_segment(loop_i,-2);
    if(check==0) {
      string str="read_RNA_template.cpp: Cannot find residues "+residue_begin+"-"+residue_end+" in file "+fstr;
      error->one(str);
    }
  }
  
  return motif_one;
}

vector<double> ReadRNATemplate::read_template_xyz(string temp_info) const {
  vector<double> Vxyz;
  vector<string> stoken=split_string(temp_info);
  string PDBid=stoken[0];
  string type_one=stoken[1];
  string str_size=stoken[2];
  string str_end=stoken[3];
  vector<string> Vends;
  stoken=split_string(str_end,"-&");
  for(size_t ii=0;ii<stoken.size();ii++)
    Vends.push_back(stoken[ii]);
  vector<int> lengths;
  stoken=split_string(str_size,"-&");
  for(size_t ii=0;ii<stoken.size();ii++)
    lengths.push_back(stoi(stoken[ii]));
  if(type_one=="Hairpin_kiss" || type_one=="Loop_kiss_loop")
    lengths[0]=lengths[0]+lengths[1]+lengths[2];
  else if(type=="Tail5" || type=="Tail3")
    lengths[0]=lengths[0]-1;
  
  string fstr=init->RNA_TEMPLATE_DIR+"pdb/"+PDBid+".pdb";
  ReadPDB readpdb_one(cbs,fstr);
  Molecule mole_one(cbs);
  mole_one.set_molecule(readpdb_one.RNA_atoms[0],Chain::rna);
  size_t nchain=Vends.size()/2;
  for(size_t ii=0;ii<nchain;ii++) {
    string residue_begin=Vends[ii*2];
    string residue_end=Vends[ii*2+1];
    vector<string> stoken=split_string(residue_begin,"/");
    string chainid_one=stoken[0];
    int check=0;
    for(size_t jj=0;jj<mole_one.size();jj++) {
      for(int kk=0;kk<(int)mole_one.chains[jj].size();kk++) {
        if(mole_one.chains[jj].chainid!=chainid_one) continue;
        if(mole_one.chains[jj][kk].get_residue_name()==residue_begin) {
          if(kk+lengths[ii]+1>=(int)mole_one.chains[jj].size()) break;
          else if(mole_one.chains[jj][kk+lengths[ii]+1].get_residue_name()!=residue_end) continue;
          else {
            vector<double> xyz=mole_one.chains[jj][kk].get_residue_xyz();
            Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
            check=1;
            continue;
          } 
        }
        if(check==1 && mole_one.chains[jj][kk].get_residue_name()==residue_end) {
          vector<double> xyz=mole_one.chains[jj][kk].get_residue_xyz();
          Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
          check=2;
          break;
        }
        if(check==1) {
          vector<double> xyz=mole_one.chains[jj][kk].get_residue_xyz();
          Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
        }
      }
      if(check==2) break;
    }
    if(check==0) {
      string str="read_RNA_template.cpp: Cannot find residues "+residue_begin+"-"+residue_end+" in file "+fstr;
      error->one(str);
    }
  }

  //string str="read_RNA_template.cpp: "+temp_info+" --> "+to_string(Vxyz.size());
  //error->message(str);
  
  return Vxyz;
}

/*
vector<double> ReadRNATemplate::read_template_xyz(Molecule &mole_one, string temp_info) {
  vector<double> Vxyz;
  vector<string> stoken=split_string(temp_info);
  string PDBid=stoken[0];
  string str_size=stoken[2];
  string str_end=stoken[3];
  vector<string> Vends;
  stoken=split_string(str_end,"-&");
  for(size_t ii=0;ii<stoken.size();ii++)
    Vends.push_back(stoken[ii]);
  vector<int> lengths;
  stoken=split_string(str_size,"-&");
  for(size_t ii=0;ii<stoken.size();ii++)
    lengths.push_back(stoi(stoken[ii]));
  if(type=="Hairpin_kiss" || type=="Loop_kiss_loop")
    lengths[0]=lengths[0]+lengths[1]+lengths[2];
  else if(type=="Tail5" || type=="Tail3")
    lengths[0]=lengths[0]-1;
  
  size_t nchain=Vends.size()/2;
  for(size_t ii=0;ii<nchain;ii++) {
    string residue_begin=Vends[ii*2];
    string residue_end=Vends[ii*2+1];
    vector<string> stoken=split_string(residue_begin,"/");
    string chainid_one=stoken[0];
    int check=0;
    for(size_t jj=0;jj<mole_one.size();jj++) {
      for(int kk=0;kk<(int)mole_one.chains[jj].size();kk++) {
        if(mole_one.chains[jj].chainid!=chainid_one) continue;
        if(mole_one.chains[jj][kk].get_residue_name()==residue_begin) {
          if(kk+lengths[ii]+1>=(int)mole_one.chains[jj].size()) break;
          else if(mole_one.chains[jj][kk+lengths[ii]+1].get_residue_name()!=residue_end) continue;
          else {
            vector<double> xyz=mole_one.chains[jj][kk].get_residue_xyz();
            Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
            check=1;
            continue;
          } 
        }
        if(check==1 && mole_one.chains[jj][kk].get_residue_name()==residue_end) {
          vector<double> xyz=mole_one.chains[jj][kk].get_residue_xyz();
          Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
          check=2;
          break;
        }
        if(check==1) {
          vector<double> xyz=mole_one.chains[jj][kk].get_residue_xyz();
          Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
        }
      }
      if(check==2) break;
    }
    if(check==0) {
      string str="read_RNA_template.cpp: Cannot find residues "+residue_begin+"-"+residue_end+" in PDB "+PDBid;
      error->one(str);
    }
  }
  
  return Vxyz;
}
*/

void ReadRNATemplate::write_templates() {
  string dir=init->RNA_TEMPLATE_DIR+type;
  if(!is_directory_exist(dir)) {
    string cmd="mkdir "+dir;
    cbs->execute_system_command(cmd);
  }
  string fstr=dir+"/"+template_name+".dat";
  ofstream outT (fstr.c_str());
  if(!outT.is_open()) {
    string str="read_RNA_template.cpp: Cannot open file "+fstr+" to write templates!";
    error->one(str);
  }
  for(size_t ii=0;ii<template_list.size();ii++) {
    outT<<"Template "<<ii+1<<endl;
    for(size_t jj=0;jj<template_list[ii].size();jj++)
      outT<<template_list[ii][jj]<<endl;
    outT<<endl;
  }
  outT.close(); outT.clear();
}

void ReadRNATemplate::write_full_template_list() {
  string dir=init->RNA_TEMPLATE_DIR+"template_list";
  if(!is_directory_exist(dir)) {
    string cmd="mkdir "+dir;
    cbs->execute_system_command(cmd);
  }
  string fstr=dir+"/"+template_name+"_list.dat";
  ofstream outT (fstr.c_str());
  if(!outT.is_open()) {
    string str="read_RNA_template.cpp: Cannot open file "+fstr+" to write template list!";
    error->one(str);
  }
  for(auto sline : full_template_list)
    outT<<sline<<endl;
  outT.close(); outT.clear();
}

void ReadRNATemplate::clear() {
  type="";
  template_name="";
  template_list.clear();
  template_list.shrink_to_fit();
  full_template_list.clear();
  full_template_list.shrink_to_fit();
}