#include "basic.h"
#include "constant.h"
#include "initialization.h"
#include "error.h"
#include "RNA_basepair.h"

using namespace CBS_NS;

RNABasePair::RNABasePair(Cobmms *cbs) : Pointers(cbs) {}

RNABasePair::RNABasePair(Cobmms *cbs, const vector<vector<string> > &Vresname_one, string ss_one) : Pointers(cbs) {
  RNABasePair::parse_RNA_basepairs(Vresname_one,ss_one);
}

RNABasePair::RNABasePair(const RNABasePair &other) : Pointers(other.cbs) {
  secstr=other.secstr;
  secstr_out=other.secstr_out;
  mapBP=other.mapBP;
  Vresnames=other.Vresnames;
  chainlength=other.chainlength;
  Mchainid_indexs=other.Mchainid_indexs;
  BPlevel=other.BPlevel;
}

RNABasePair& RNABasePair::operator=(const RNABasePair &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    secstr=other.secstr;
    secstr_out=other.secstr_out;
    mapBP=other.mapBP;
    Vresnames=other.Vresnames;
    chainlength=other.chainlength;
    Mchainid_indexs=other.Mchainid_indexs;
    BPlevel=other.BPlevel;
  }

  return *this;
}

RNABasePair::~RNABasePair() {
  secstr="";
  secstr_out="";
  mapBP.clear();
  Vresnames.clear();
  Vresnames.shrink_to_fit();
  chainlength.clear();
  chainlength.shrink_to_fit();
  Mchainid_indexs.clear();
  BPlevel.clear();
}

void RNABasePair::set_RNA_basepairs(const RNABasePair &basepair_one) {
  secstr=basepair_one.secstr;
  secstr_out=basepair_one.secstr_out;
  mapBP=basepair_one.mapBP;
  Vresnames=basepair_one.Vresnames;
  chainlength=basepair_one.chainlength;
  Mchainid_indexs=basepair_one.Mchainid_indexs;
  BPlevel=basepair_one.BPlevel;
}

void RNABasePair::parse_RNA_basepairs(const vector<vector<string> > &Vresname_one, string ss_one) {
  secstr="";
  mapBP.clear();
  Vresnames.clear();
  chainlength.clear();
  Mchainid_indexs.clear();
  
  size_t nres=0;
  for(size_t ii=0;ii<Vresname_one.size();ii++) {
    string chainid_one="";
    for(size_t jj=0;jj<Vresname_one[ii].size();jj++) {
      string res_one=Vresname_one[ii][jj];
      vector<string> stoken=split_string(res_one,"/");
      if(stoken.size()!=3) {
        string str="RNA_basepair.cpp (parse_RNA_basepairs): Wrong nucleotide name(chainid/resid/resname) --> "+res_one;
        error->one(str);
      }
      if(jj==0)
        chainid_one=stoken[0];
      Vresnames.push_back(res_one);
    }
    vector<size_t> Vindex;
    Vindex.push_back(nres);
    nres+=Vresname_one[ii].size();
    Vindex.push_back(nres-1);
    if(Mchainid_indexs.count(chainid_one)==0)
      Mchainid_indexs.emplace(chainid_one,Vindex);
    else
      Mchainid_indexs[chainid_one][1]=Vindex[1];
    chainlength.push_back((int)Vresname_one[ii].size());
  }
  
  if(Vresnames.size()==0) {
    string str="RNA_basepair.cpp: Please provide nucleotides' name first!";
    error->one(str);
  }

  vector<string> stoken=split_string(ss_one);
  ss_one=stoken[0];
  for(size_t ii=1;ii<stoken.size();ii++)
    ss_one+=" "+stoken[ii];
  for(size_t ii=0;ii<ss_one.size();ii++) {
    char ss=ss_one[ii];
    if(ss==' ' || ss=='-' || ss=='&') continue;
    else if(count(init->RNA_BP_NOTATION.begin(),init->RNA_BP_NOTATION.end(),ss)) secstr+=to_string(ss);
    else {
      string str="RNA_basepair.cpp: Unknown base pair notation --> "+to_string(ss);
      error->one(str);
    }
  }

  int nt=((int)init->RNA_BP_NOTATION.size()-1)/2;
  for(int ii=0;ii<nt;ii++) {
    char Ln=init->RNA_BP_NOTATION[ii*2+1];
    char Rn=init->RNA_BP_NOTATION[ii*2+2];
    int nL=0,nR=0;
    for(size_t jj=0;jj<secstr.size();jj++) {
      if(secstr[jj]==Ln) nL++;
      else if(secstr[jj]==Rn) nR++;
    }
    if(nL!=nR) {
      string str="RNA_basepair.cpp: Wrong secondary structure --> (n('"+to_string(Ln)+"')="+to_string(nL)+"/";
      str+="(n('"+to_string(Rn)+"')="+to_string(nR)+"\n"+ss_one;
      error->one(str);
    }
  }

  if(Vresnames.size()!=secstr.size()) {
    string str="RNA_basepair.cpp: Sequence and 2D structure is not matched --> \n";
    str+=to_string(Vresnames.size())+" nucleotides : "+to_string(secstr.size())+" symbols ("+secstr+")";
    error->one(str);
  }

  string top=secstr;
  for(int ii=0;ii<nt;ii++) {
    string Ln=to_string(init->RNA_BP_NOTATION[ii*2+1]);
    string Rn=to_string(init->RNA_BP_NOTATION[ii*2+2]);
    int viab=0;
    if(top.find(Ln)!=string::npos) {
      do{
        int a=top.find_last_of(Ln);
        string tmp=top.substr(a+1);
        int b=tmp.find_first_of(Rn);
        if(a<0 || b<0) {
          string str="RNA_basepair.cpp: Wrong secondary structure --> \n"+top;
          error->one(str);
        }
        b=a+1+b;
        string bpi=Vresnames[a];
        string bpj=Vresnames[b];
        mapBP.insert(pair<string,string>(bpi,bpj));
        top[a]='.';
        top[b]='.';
        if(top.find(Ln)==string::npos) viab=1;
      }while(viab==0);
    }
  }

  for(auto& it : mapBP) {
    if(RNABasePair::check_noncanonical_basepair(it.first,it.second)) {
      string str="RNA_basepair.cpp: Noncanonical base pair --> "+it.first+"-"+it.second;
      error->warning(str);
    }
  }

  RNABasePair::set_basepair_level();
}

void RNABasePair::set_RNA_basepairs(const vector<vector<string> > &Vresname_one, const map<string,string> &BP_one) {
  secstr="";
  mapBP.clear();
  Vresnames.clear();
  chainlength.clear();
  Mchainid_indexs.clear();
  
  size_t nres=0;
  for(size_t ii=0;ii<Vresname_one.size();ii++) {
    string chainid_one="";
    for(size_t jj=0;jj<Vresname_one[ii].size();jj++) {
      string res_one=Vresname_one[ii][jj];
      vector<string> stoken=split_string(res_one,"/");
      if(stoken.size()!=3) {
        string str="RNA_basepair.cpp (set_RNA_basepairs): Wrong nucleotide name(chainid/resid/resname) --> "+res_one;
        error->one(str);
      }
      if(jj==0)
        chainid_one=stoken[0];
      Vresnames.push_back(res_one);
    }
    vector<size_t> Vindex;
    Vindex.push_back(nres);
    nres+=Vresname_one[ii].size();
    Vindex.push_back(nres-1);
    if(Mchainid_indexs.count(chainid_one)==0)
      Mchainid_indexs.emplace(chainid_one,Vindex);
    else 
      Mchainid_indexs[chainid_one][1]=Vindex[1];
    chainlength.push_back((int)Vresname_one[ii].size());
  }

  if(Vresnames.size()==0) {
    string str="RNA_basepair.cpp: Please provide nucleotides' name first!";
    error->one(str);
  }

  for(auto& it : BP_one) {
    size_t ni=RNABasePair::get_residue_index(it.first);
    size_t nj=RNABasePair::get_residue_index(it.second);
    if(ni<nj)
      mapBP.emplace(it.first,it.second);
    else
      mapBP.emplace(it.second,it.first);

    if(RNABasePair::check_noncanonical_basepair(it.first,it.second)) {
      string str="RNA_basepair.cpp: Noncanonical base pair --> "+it.first+"-"+it.second;
      error->warning(str);
    }
  }

  secstr="";
  for(size_t ii=0;ii<Vresnames.size();ii++)
    secstr+=".";

  RNABasePair::set_basepair_level();
}

void RNABasePair::add_RNA_basepair(string bpi, string bpj) {
  vector<string> stoken=split_string(bpi,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (add_RNA_basepair): Wrong nucleotide name(chainid/resid/resname) --> "+bpi;
    error->one(str);
  }
  stoken=split_string(bpj,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (add_RNA_basepair): Wrong nucleotide name(chainid/resid/resname) --> "+bpj;
    error->one(str);
  }

  if(RNABasePair::is_basepaired(bpi)) {
    string str="RNA_basepair.cpp: Nucleotide "+bpi+" has been paired.";
    error->one(str);
  }
  if(RNABasePair::is_basepaired(bpj)) {
    string str="RNA_basepair.cpp: Nucleotide "+bpj+" has been paired.";
    error->one(str);
  }

  if(RNABasePair::check_noncanonical_basepair(bpi,bpj)) {
    string str="RNA_basepair.cpp: Noncanonical base pair --> "+bpi+"-"+bpj;
    error->warning(str);
  }

  size_t ni=RNABasePair::get_residue_index(bpi);
  size_t nj=RNABasePair::get_residue_index(bpj);

  if(ni<nj)
    mapBP.insert(pair<string,string>(bpi,bpj));
  else
    mapBP.insert(pair<string,string>(bpj,bpi));

  RNABasePair::set_basepair_level();
}

bool RNABasePair::is_basepaired(string bpi, string bpj) const {
  vector<string> stoken=split_string(bpi,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (is_basepaired): Wrong nucleotide name(chainid/resid/resname) --> "+bpi;
    error->one(str);
  }
  stoken=split_string(bpj,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (is_basepaired): Wrong nucleotide name(chainid/resid/resname) --> "+bpj;
    error->one(str);
  }

  int check=0;
  for(auto& it : mapBP) {
    if(it.first==bpi && it.second==bpj) {check=1; break;}
    else if(it.first==bpj && it.second==bpi) {check=1; break;}
  }

  return check;
}

bool RNABasePair::is_basepaired(string bpi) const {
  vector<string> stoken=split_string(bpi,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (is_basepaired): Wrong nucleotide name(chainid/resid/resname) --> "+bpi;
    error->one(str);
  }

  int check=0;
  for(auto& it : mapBP) {
    if(it.first==bpi || it.second==bpi) {check=1; break;}
  }

  return check;
}

size_t RNABasePair::size() const {
  return mapBP.size();
}

void RNABasePair::clear() {
  secstr="";
  secstr_out="";
  mapBP.clear();
  Vresnames.clear();
  Vresnames.shrink_to_fit();
  chainlength.clear();
  chainlength.shrink_to_fit();
  Mchainid_indexs.clear();
  BPlevel.clear();
}

size_t RNABasePair::get_residue_index(string res_one, size_t index_offset) const {
  vector<string> stoken=split_string(res_one,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (get_residue_index): Wrong nucleotide name(chainid/resid/resname) --> "+res_one;
    error->one(str);
  }
  
  size_t ni=0;
  int flag_check=0;
  for(size_t ii=index_offset; ii<Vresnames.size();ii++) {
    if(Vresnames[ii]==res_one) {
      flag_check=1;
      ni=ii;
      break;
    }
  }

  if(flag_check==0) {
    string str="RNA_basepair.cpp: Cannot find nucleotide "+res_one;
    error->one(str);
  }

  return ni;
}

int RNABasePair::get_basepair_level(string bpi, string bpj) const {
  vector<string> stoken=split_string(bpi,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (get_basepair_level): Wrong nucleotide name(chainid/resid/resname) --> "+bpi;
    error->one(str);
  }
  stoken=split_string(bpj,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (get_basepair_level): Wrong nucleotide name(chainid/resid/resname) --> "+bpj;
    error->one(str);
  }

  int nlevel=0;
  string bp=bpi+" "+bpj;
  if(BPlevel.count(bp)) 
    nlevel=BPlevel.at(bp);
  else {
    bp=bpj+" "+bpi;
    if(BPlevel.count(bp)) nlevel=BPlevel.at(bp);
  }

  return nlevel;
}

string RNABasePair::get_segment_secondary_structure(string res_begin, string res_end) const{
  vector<string> stoken=split_string(res_begin,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (get_segment_secondary_structure): Wrong nucleotide name(chainid/resid/resname) --> "+res_begin;
    error->one(str);
  }
  stoken=split_string(res_end,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (get_segment_secondary_structure): Wrong nucleotide name(chainid/resid/resname) --> "+res_end;
    error->one(str);
  }
  
  string chainid_one=stoken[0];
  size_t index_begin0=0, index_end0=0;
  if(Mchainid_indexs.count(chainid_one)==0) {
    string str="RNA_basepair.cpp: Cannot find nucleotide "+res_end;
    error->one(str);
  }
  else {
    vector<size_t> Vindex=Mchainid_indexs.at(chainid_one);
    index_begin0=Vindex[0];
    index_end0=Vindex[1];
  }

  size_t index_begin=RNABasePair::get_residue_index(res_begin);
  size_t index_end=0;
  if(init->CircRNA_chains.size()>0 && count(init->CircRNA_chains.begin(),init->CircRNA_chains.end(),chainid_one))
    index_end=RNABasePair::get_residue_index(res_end);
  else
    index_end=RNABasePair::get_residue_index(res_end,index_begin);

  string secstr_one="";
  if(index_begin>index_end) {
    if(init->CircRNA_chains.size()>0 && count(init->CircRNA_chains.begin(),init->CircRNA_chains.end(),chainid_one)) {
      secstr_one=secstr.substr(index_begin,index_end0-index_begin+1);
      secstr_one+=secstr.substr(index_begin0,index_end-index_begin0+1);
    }
    else {
      string str="RNA_basepair.cpp: Invalid residue_begin & residue_end --> "+res_begin+"("+to_string(index_begin)+") ";
      str+=res_end+"("+to_string(index_end)+")";
      error->one(str);
    }
  }
  else {
    size_t len=index_end-index_begin+1;
    secstr_one=secstr.substr(index_begin,len);
  }
  
  return secstr_one;
}

string RNABasePair::previous_residue(string res_one) const {
  string res_pre=ResidueConst::RESIDUE_NULL;

  vector<string> stoken=split_string(res_one,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (previous_residue): Wrong nucleotide name(chainid/resid/resname) --> "+res_one;
    error->one(str);
  }
  string chainid_one=stoken[0];

  size_t index_begin=0, index_end=0;
  if(Mchainid_indexs.count(chainid_one)==0) {
    string str="RNA_basepair.cpp: Cannot find residue "+res_one;
    error->warning(str);
    return res_pre;
  }
  else {
    vector<size_t> Vindex=Mchainid_indexs.at(chainid_one);
    index_begin=Vindex[0];
    index_end=Vindex[1];
  }

  size_t ni=0;
  int flag_check=0;
  for(size_t ii=0;ii<Vresnames.size();ii++) {
    if(Vresnames[ii]==res_one) {
      ni=ii;
      flag_check=1;
      break;
    }
  }
  if(flag_check==0) {
    string str="RNA_basepair.cpp: Cannot find residue "+res_one;
    error->warning(str);
    return res_pre;
  }
  
  if(ni>index_begin) {
    res_pre=Vresnames[ni-1];
  }
  else if(ni==index_begin && init->CircRNA_chains.size()>0 && count(init->CircRNA_chains.begin(),init->CircRNA_chains.end(),chainid_one)) {
    res_pre=Vresnames[index_end];
  }
  else {
    return res_pre;
  }

  stoken=split_string(res_pre,"/");
  if(stoken[0]!=chainid_one) res_pre=ResidueConst::RESIDUE_NULL;

  return res_pre;
}

string RNABasePair::next_residue(string res_one) const {
  string res_next=ResidueConst::RESIDUE_NULL;

  vector<string> stoken=split_string(res_one,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (next_residue): Wrong nucleotide name(chainid/resid/resname) --> "+res_one;
    error->one(str);
  }
  string chainid_one=stoken[0];

  size_t index_begin=0, index_end=0;
  if(Mchainid_indexs.count(chainid_one)==0) {
    string str="RNA_basepair.cpp: Cannot find residue "+res_one;
    error->warning(str);
    return res_next;
  }
  else {
    vector<size_t> Vindex=Mchainid_indexs.at(chainid_one);
    index_begin=Vindex[0];
    index_end=Vindex[1];
  }

  size_t ni=0;
  int flag_check=0;
  for(size_t ii=0;ii<Vresnames.size();ii++) {
    if(Vresnames[ii]==res_one) {
      ni=ii;
      flag_check=1;
      break;
    }
  }
  if(flag_check==0) {
    string str="RNA_basepair.cpp: Cannot find residue "+res_one;
    error->warning(str);
    return res_next;
  }
  
  if(ni<index_end) {
    res_next=Vresnames[ni+1];
  }
  else if(ni==index_end && init->CircRNA_chains.size()>0 && count(init->CircRNA_chains.begin(),init->CircRNA_chains.end(),chainid_one)) {
    res_next=Vresnames[index_begin];
  }
  else {
    return res_next;
  }

  stoken=split_string(res_next,"/");
  if(stoken[0]!=chainid_one) res_next=ResidueConst::RESIDUE_NULL;

  return res_next;
}

bool RNABasePair::check_noncanonical_basepair(string bpi, string bpj) const{
  string BPs[]={"GC","CG","AU","UA","GU","UG"};
  vector<string> stoken=split_string(bpi,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (check_noncanonical_basepair): Wrong nucleotide name(chainid/resid/resname) --> "+bpi;
    error->one(str);
  }
  string resname_i=stoken[2];
  stoken=split_string(bpj,"/");
  if(stoken.size()!=3) {
    string str="RNA_basepair.cpp (check_noncanonical_basepair): Wrong nucleotide name(chainid/resid/resname) --> "+bpj;
    error->one(str);
  }
  string resname_j=stoken[2];
  string bp=resname_i+resname_j;
  if(count(begin(BPs),end(BPs),bp)) return false;
  else return true;
}

//(,)=1, [,]=2, {,}=3, <,>=4, a,A=5......
void RNABasePair::set_basepair_level() {
  BPlevel.clear();

  map<string,int> nt_indexs;
  for(auto& it : mapBP) {
    int ni=static_cast<int>(RNABasePair::get_residue_index(it.first));
    int nj=static_cast<int>(RNABasePair::get_residue_index(it.second));
    nt_indexs.insert(pair<string,int>(it.first,ni));
    nt_indexs.insert(pair<string,int>(it.second,nj));
  }

  map<string,string> crossBP;
  int nlevel=1;
  for(auto& it : mapBP) {
    int check=0;
    string bpi=it.first;
    string bpj=it.second;
    int a1=nt_indexs[bpi];
    int b1=nt_indexs[bpj];
    if(a1>b1) {swap(a1,b1); swap(bpi,bpj);}
    for(auto& it1 : mapBP) {
      int a2=nt_indexs[it1.first];
      int b2=nt_indexs[it1.second];
      if(a2>b2) swap(a2,b2);
      if((a1<a2 && a2<b1 && b1<b2) || (a2<a1 && a1<b2 && b2<b1))
        {check=1; break;}
    }
    if(check)
      crossBP.insert(pair<string,string>(bpi,bpj));
    else {
      string bp=bpi+" "+bpj;
      BPlevel.insert(pair<string,int>(bp,nlevel));
    }
  }

  vector<map<string,string> > Vhx;
  while(crossBP.size()>0) {
    map<string,string> hx;
    hx.emplace(crossBP.begin()->first,crossBP.begin()->second);
    crossBP.erase(crossBP.begin());
    int check=1;
    while(check) {
      check=0;
      for(auto& it : crossBP) {
        int a1=nt_indexs[it.first];
        int b1=nt_indexs[it.second];
        for(auto& it1 : hx) {
          int a2=nt_indexs[it1.first];
          int b2=nt_indexs[it1.second];
          if(abs(a1-a2)<=2 && abs(b1-b2)<=2) {
            hx.insert(pair<string,string>(it.first,it.second));
            crossBP.erase(it.first);
            check=1;
            break;
          }
        }
        if(check) break;
      }
    }
    Vhx.push_back(hx);
  }
  
  while(Vhx.size()>0) {
    vector<map<string,string> > cross_hx;
    while(1) {
      int nmax=0;
      int index_max=-1;
      for(size_t ii=0;ii<Vhx.size();ii++) {
        int nc=0;
        int a1=nt_indexs[Vhx[ii].begin()->first];
        int b1=nt_indexs[Vhx[ii].begin()->second];
        for(size_t jj=0;jj<Vhx.size();jj++) {
          if(ii==jj) continue;
          int a2=nt_indexs[Vhx[jj].begin()->first];
          int b2=nt_indexs[Vhx[jj].begin()->second];
          if((a1<a2 && a2<b1 && b1<b2) || (a2<a1 && a1<b2 && b2<b1))
            nc+=Vhx[jj].size();
        }
        if(nc>nmax) {
          nmax=nc;
          index_max=ii;
        }
      }
      if(nmax==0) break;
      cross_hx.push_back(Vhx[index_max]);
      Vhx.erase(Vhx.begin()+index_max);
    }

    for(size_t ii=0;ii<Vhx.size();ii++) {
      for(auto& it : Vhx[ii]) {
        string bp=it.first+" "+it.second;
        BPlevel.insert(pair<string,int>(bp,nlevel));
      }
    }

    Vhx=cross_hx;
    nlevel++;
  }

  for(auto& it : BPlevel) {
    vector<string> stoken=split_string(it.first);
    int a=nt_indexs[stoken[0]];
    int b=nt_indexs[stoken[1]];
    int nlevel=it.second;
    if(a>b) swap(a,b);
    secstr[a]=init->RNA_BP_NOTATION[2*nlevel-1];
    secstr[b]=init->RNA_BP_NOTATION[2*nlevel];
  }

  secstr_out=secstr;
  int nc=chainlength[0];
  for(size_t ii=1;ii<chainlength.size();ii++) {
    secstr_out.insert(nc,1,'&');
    nc+=chainlength[ii];
    nc++;
  }
}

ostream& operator <<(ostream &output, RNABasePair &bp_one) {
  output<<bp_one.secstr_out<<endl;
  output<<"List of "<<bp_one.mapBP.size()<<" base pairs:"<<endl;
  output<<"#Nt_i  Nt_j  BP_level"<<endl;
  for(auto& it : bp_one.BPlevel)
    output<<it.first<<" "<<it.second<<endl;

  return output;
}