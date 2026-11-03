#include "basic.h"
#include "error.h"
#include "RNA_motif.h"

using namespace CBS_NS;

RNAMotif::RNAMotif(Cobmms *cbs) : Molecule(cbs) {
  flag_rigid=0;
}

RNAMotif::RNAMotif(const RNAMotif &other) : Molecule(other) {
  type=other.type;
  sequences=other.sequences;
  ends=other.ends;
  lengths=other.lengths;
  secstrs=other.secstrs;
  flag_rigid=other.flag_rigid;
}

RNAMotif& RNAMotif::operator=(const RNAMotif &other) {
  if(this != &other) {
    Molecule::operator=(other);
    type=other.type;
    sequences=other.sequences;
    ends=other.ends;
    lengths=other.lengths;
    secstrs=other.secstrs;
    flag_rigid=other.flag_rigid;
  }

  return *this;
}

RNAMotif::~RNAMotif() {
  type="";
  chain_resnames.clear();
  chain_resnames.shrink_to_fit();
  chains.clear();
  chains.shrink_to_fit();
  chain_ids.clear();
  chain_ids.shrink_to_fit();
  chain_labels.clear();
  chain_labels.shrink_to_fit();
  Vchainid.clear();
  Vchainid.shrink_to_fit();
  Vresnames.clear();
  Vresnames.shrink_to_fit();
  sequences.clear();
  sequences.shrink_to_fit();
  ends.clear();
  ends.shrink_to_fit();
  lengths.clear();
  lengths.shrink_to_fit();
  secstrs.clear();
  secstrs.shrink_to_fit();
  flag_rigid=0;
}

void RNAMotif::set_RNA_motif(const RNAMotif &motif_one) {
  type=motif_one.type;
  chain_resnames=motif_one.chain_resnames;
  chains=motif_one.chains;
  chain_ids=motif_one.chain_ids;
  chain_labels=motif_one.chain_labels;
  Vchainid=motif_one.Vchainid;
  Vresnames=motif_one.Vresnames;
  sequences=motif_one.sequences;
  ends=motif_one.ends;
  lengths=motif_one.lengths;
  secstrs=motif_one.secstrs;
  flag_rigid=motif_one.flag_rigid;
}

void RNAMotif::insert_RNA_segment(const Chain &chain_one, int length_offset) {
  if(chain_one.chain_label!=Chain::rna) {
    string str="RNA_motif.cpp: Not RNA chain in insert_RNA_segment! --> ";
    str+=chain_one.first_of_residue()+"-"+chain_one.last_of_residue();
    error->one(str);
  }
  Molecule::insert_chain(chain_one);
  string res_one=chain_one.first_of_residue();
  ends.push_back(res_one);
  res_one=chain_one.last_of_residue();
  ends.push_back(res_one);
  int len=(int)chain_one.size()+length_offset;
  lengths.push_back(len);
  RNAMotif::get_RNA_motif_sequence();
}

void RNAMotif::set_RNA_motif(string info, int index_offset) {
  vector<string> stoken=split_string(info);
  if((int)stoken.size()<5+index_offset) {
    string str="RNA_motif.cpp: Wrong motif info --> "+info+" index_offset="+to_string(index_offset);
    error->one(str);
  }

  type=stoken[index_offset];

  string str_size=stoken[1+index_offset];
  string str_end=stoken[2+index_offset];
  string str_seq=stoken[3+index_offset];
  string str_secstr=stoken[4+index_offset];
  lengths.clear();
  size_t found=str_size.find("L");
  if(found!=string::npos) {
    str_size.erase(str_size.begin()+found);
  }
  stoken=split_string(str_size,"-&");
  for(size_t ii=0;ii<stoken.size();ii++)
    lengths.push_back(stoi(stoken[ii]));
  ends.clear();
  stoken=split_string(str_end,"-&");
  for(size_t ii=0;ii<stoken.size();ii++)
    ends.push_back(stoken[ii]);
  sequences.clear();
  stoken=split_string(str_seq,"-&");
  for(size_t ii=0;ii<stoken.size();ii++)
    sequences.push_back(stoken[ii]);
  secstrs.clear();
  stoken=split_string(str_secstr,"-&");
  for(size_t ii=0;ii<stoken.size();ii++)
    secstrs.push_back(stoken[ii]);
    
  vector<vector<string> > Vresname_one;
  for(size_t ii=0;ii<sequences.size();ii++) {
    string nt_begin=ends[ii*2];
    vector<string> stoken=split_string(nt_begin,"/");
    string chainid_one=stoken[0];
    int resid_one=stoi(stoken[1]);
    vector<string> Vtmp;
    for(size_t jj=0;jj<sequences[ii].size();jj++) {
      string res_one=chainid_one+"/"+to_string(jj+resid_one)+"/"+to_string((char)sequences[ii][jj]);
      Vtmp.push_back(res_one);
    }
    Vresname_one.push_back(Vtmp);
  }

  Molecule::set_molecule(Vresname_one,Chain::rna);
  Molecule::simplify_sidechain();

  flag_rigid=0;
}

void RNAMotif::get_RNA_motif_sequence() {
  sequences.clear();
  for(size_t ii=0;ii<chains.size();ii++) {
    string seq_one=chains[ii].get_sequence();
    sequences.push_back(seq_one);
  }
  /*
  for(size_t ii=0;ii<chain_resnames.size();ii++) {
    string seq_one="";
    for(size_t jj=0;jj<chain_resnames[ii].size();jj++)
      seq_one+=chain_resnames[ii][jj];
    sequences.push_back(seq_one);
  }
  */
}

string RNAMotif::get_RNA_motif_name() const {
  string name=RNAMotif::get_RNA_motif_info();
  vector<string> stoken=split_string(name);
  name=stoken[0]+"#";
  string str_size=stoken[1];
  /*
  if(type.find("way_junction_kiss")!=string::npos) {
    size_t found=str_size.find_first_of("-");
    if(found!=string::npos) 
      str_size=str_size.substr(0,found+1);
  }
  */
  name+=str_size;

  return name;
}

string RNAMotif::get_RNA_motif_info() const {
  string str_size=to_string(lengths[0]);
  if(type.find("way_junction_kiss")!=string::npos) {
    vector<string> stoken=split_string(type,"-");
    int nloop=stoi(stoken[0]);
    if(nloop==(int)lengths.size()) {
      for(size_t ii=1;ii<lengths.size()-1;ii++)
        str_size+="&"+to_string(lengths[ii]);
      str_size+="&"+to_string(lengths.back());
    }
    else if((nloop+3)==(int)lengths.size()) {
      for(size_t ii=1;ii<lengths.size()-4;ii++)
        str_size+="&"+to_string(lengths[ii]);
      str_size+="&"+to_string(lengths[lengths.size()-4]);
      //for(size_t ii=lengths.size()-3;ii<lengths.size();ii++) 
        //str_size+="-"+to_string(lengths[ii]);
    }
    else {
      for(size_t ii=1;ii<lengths.size();ii++)
        str_size+="&"+to_string(lengths[ii]);
    }
  }
  else if(type!="Helix" && type!="THelix") {
    for(size_t ii=1;ii<lengths.size();ii++)
      str_size+="&"+to_string(lengths[ii]);
  }
  string str_end=ends[0]+"-"+ends[1];
  for(size_t ii=2;ii<ends.size();ii++) {
    if(ii%2==0) str_end+="&"+ends[ii];
    else str_end+="-"+ends[ii];
  }
  string str_seq=sequences[0];
  for(size_t ii=1;ii<sequences.size();ii++)
    str_seq+="&"+sequences[ii];
  string str_secstr=secstrs[0];
  for(size_t ii=1;ii<secstrs.size();ii++)
    str_secstr+="&"+secstrs[ii];
  string info=type+" "+str_size+" "+str_end+" "+str_seq+" "+str_secstr;

  return info;
}

void RNAMotif::set_RNA_motif_secondary_structure() {
  secstrs.clear();
  if(type=="Helix") {
    string str="";
    for(int ii=0;ii<lengths[0];ii++)
      str+="(";
    secstrs.push_back(str);
    str="";
    for(int ii=0;ii<lengths[1];ii++)
      str+=")";
    secstrs.push_back(str);
  }
  else if(type=="Unpaired_chain") {
    string str="";
    for(int ii=0;ii<lengths[0];ii++)
      str+=".";
    secstrs.push_back(str);
  }
  else if(type=="Tail5") {
    string str="";
    for(int ii=0;ii<lengths[0];ii++)
      str+=".";
    str+="(";
    secstrs.push_back(str);
  }
  else if(type=="Tail3") {
    string str=")";
    for(int ii=0;ii<lengths[0];ii++)
      str+=".";
    secstrs.push_back(str);
  }
  else if(type=="Hairpin") {
    string str="(";
    for(int ii=0;ii<lengths[0];ii++)
      str+=".";
    str+=")";
    secstrs.push_back(str);
  }
  else if(type=="Hairpin_kiss") {
    string str="(";
    for(int ii=0;ii<lengths[0];ii++)
      str+=".";
    for(int ii=0;ii<lengths[1];ii++)
      str+="[";
    for(int ii=0;ii<lengths[2];ii++)
      str+=".";
    str+=")";
    secstrs.push_back(str);
  }
  else if(type=="Loop_kiss_loop") {
    string str=")";
    for(int ii=0;ii<lengths[0];ii++)
      str+=".";
    for(int ii=0;ii<lengths[1];ii++)
      str+="[";
    for(int ii=0;ii<lengths[2];ii++)
      str+=".";
    str+="(";
    secstrs.push_back(str);
  }
  else if(type=="Tertiary_helix") {
    string str="";
    for(int ii=0;ii<lengths[0];ii++)
      str+="[";
    secstrs.push_back(str);
  }
  else if(type=="Loop") {
    string str=")";
    for(int ii=0;ii<lengths[0];ii++)
      str+=".";
    str+="(";
    secstrs.push_back(str);
  }
  else if(type.find("-way_junction")!=string::npos && type.find("kiss")==string::npos) {
    for(size_t ii=0;ii<lengths.size();ii++) {
      string str=")";
      for(int jj=0;jj<lengths[ii];jj++)
        str+=".";
      str+="(";
      secstrs.push_back(str);
    }
  }
  else if(type.find("-way_junction_kiss")!=string::npos) {
    vector<string> stoken=split_string(type,"-");
    int nloop=stoi(stoken[0]);
    if(nloop==(int)lengths.size()) {
      for(size_t ii=0;ii<lengths.size();ii++) {
        string str=")";
        for(int jj=0;jj<lengths[ii];jj++)
          str+=".";
        str+="(";
        secstrs.push_back(str);
      }
    }
    else if((nloop+3)==(int)lengths.size()) {
      for(size_t ii=0;ii<lengths.size()-4;ii++) {
        string str=")";
        for(int jj=0;jj<lengths[ii];jj++)
          str+=".";
        str+="(";
        secstrs.push_back(str);
      }
      string str=")";
      for(int jj=0;jj<lengths[lengths.size()-3];jj++)
        str+=".";
      for(int jj=0;jj<lengths[lengths.size()-2];jj++)
        str+="[";
      for(int jj=0;jj<lengths.back();jj++)
        str+=".";
      str+="(";
      secstrs.push_back(str);
    }
    else {
      string str="RNA_motif.cpp: Wrong motif type! --> "+get_RNA_motif_info();
      error->one(str);
    }
  }
  else {
    string str="RNA_motif.cpp: Unknown motif type! --> "+type;
    error->one(str);
  }
}

void RNAMotif::set_RNA_motif_xyz(RNAMotif &motif_one) {
  if(chains.size()!=motif_one.chains.size()) {
    string str="RNA_motif.cpp: Number of chains is not matched in setting of xyz --> ";
    str+=to_string(chains.size())+"/"+to_string(motif_one.chains.size());
    error->one(str);
  }
  for(size_t ii=0;ii<chains.size();ii++) {
    chains[ii].set_chain_xyz(motif_one.chains[ii]);
  }
}

bool RNAMotif::is_identical_junction(RNAMotif &motif_one) {
  if(type!=motif_one.type) 
    return false;

  if(type.find("way_junction")==string::npos || type.find("way_junction_kiss")!=string::npos) {
    string str="RNA_motif.cpp: Not Junction --> "+type;
    error->warning(str);
    return false;
  }
  if(motif_one.type.find("way_junction")==string::npos || motif_one.type.find("way_junction_kiss")!=string::npos) {
    string str="RNA_motif.cpp: Not Junction --> "+motif_one.type;
    error->warning(str);
    return false;
  }

  string length_code="";
  for(size_t ii=0;ii<lengths.size();ii++) {
    length_code+=to_string(lengths[ii]);
  }
  string seq_full=sequences[0];
  for(size_t ii=1;ii<sequences.size();ii++) {
    seq_full+="&"+sequences[ii];
  }
  size_t index_i=0;
  int flag_check=0;
  for(size_t ii=0;ii<motif_one.lengths.size();ii++) {
    index_i=ii;
    string length_code_i="";
    for(size_t jj=index_i;jj<motif_one.lengths.size();jj++) {
      length_code_i+=to_string(motif_one.lengths[jj]);
    }
    for(size_t jj=0;jj<index_i;jj++) {
      length_code_i+=to_string(motif_one.lengths[jj]);
    }
    if(length_code!=length_code_i) continue;
    string seq_full_i=motif_one.sequences[index_i];
    for(size_t jj=index_i+1;jj<motif_one.sequences.size();jj++) {
      seq_full_i+="&"+motif_one.sequences[jj];
    }
    for(size_t jj=0;jj<index_i;jj++) {
      seq_full_i+="&"+motif_one.sequences[jj];
    }
    if(seq_full==seq_full_i) {
      flag_check=1;
      break;
    }
  }

  if(flag_check==0) 
    return false;
  else {
    RNAMotif mtf(cbs);
    mtf.type=motif_one.type;
    for(size_t jj=index_i;jj<motif_one.chains.size();jj++) {
      mtf.insert_RNA_segment(motif_one.chains[jj],-2);
      mtf.secstrs.push_back(motif_one.secstrs[jj]);
    }
    for(size_t jj=0;jj<index_i;jj++) {
      mtf.insert_RNA_segment(motif_one.chains[jj],-2);
      mtf.secstrs.push_back(motif_one.secstrs[jj]);
    }

    motif_one=mtf;

    return true;
  }
}

void RNAMotif::clear() {
  type="";
  Molecule::clear();
  sequences.clear();
  sequences.shrink_to_fit();
  ends.clear();
  ends.shrink_to_fit();
  lengths.clear();
  lengths.shrink_to_fit();
  secstrs.clear();
  secstrs.shrink_to_fit();
  flag_rigid=0;
}

size_t RNAMotif::size() const {
  return chains.size();
}

Chain & RNAMotif::operator [](size_t index) {
  return chains[index];
}

ostream& operator <<(ostream &output,RNAMotif &motif_one) {
  output<<"REMARK  "<<motif_one.get_RNA_motif_info()<<endl;
  /*
  output<<"REMARK  SEQ:";
  for(size_t ii=0;ii<motif_one.sequences.size();ii++) {
    output<<motif_one.sequences[ii];
    if(ii<motif_one.sequences.size()-1)
      output<<"&";
  }
  output<<endl;
  output<<"REMARK  SECSTR:";
  for(size_t ii=0;ii<motif_one.secstrs.size();ii++) {
    output<<motif_one.secstrs[ii];
    if(ii<motif_one.secstrs.size()-1)
      output<<"&";
  }
  output<<endl<<endl;
  */

  int index=0;
  for(size_t ii=0;ii<motif_one.chains.size();ii++) {
    index=motif_one.chains[ii].set_atomid(index);
    output<<motif_one.chains[ii];
  }

  return output;
}