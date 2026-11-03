#include "basic.h"
#include "initialization.h"
#include "error.h"
#include "chain.h"

using namespace CBS_NS;

const int RESID_MAX=9999;

enum{BEGIN,END,RNA,DNA,PROTEIN,LIGAND,CG,OTHER,EMPTY};
const int Chain::begin=BEGIN;
const int Chain::end=END;

const int Chain::rna=RNA;
const int Chain::dna=DNA;
const int Chain::protein=PROTEIN;
const int Chain::ligand=LIGAND;
const int Chain::cg=CG;
const int Chain::other=OTHER;
const int Chain::empty=EMPTY;

Chain::Chain(Cobmms *cbs) : Pointers(cbs) {
  chainid="A";
  chain_label=-1;
}

Chain::Chain(Cobmms *cbs, vector<string> &Vseq_one, int label_one, string chainid_one, int resid_offset) 
 : Pointers(cbs) 
{
  Chain::set_chain(Vseq_one,label_one,chainid_one,resid_offset);
}

Chain::Chain(const Chain &other) : Pointers(other.cbs) {
  chainid=other.chainid;
  resnames=other.resnames;
  residues=other.residues;
  chain_label=other.chain_label;
}

Chain& Chain::operator=(const Chain &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    chainid=other.chainid;
    resnames=other.resnames;
    residues=other.residues;
    chain_label=other.chain_label;
  }

  return *this;
}

Chain::~Chain() {
  chainid="";
  chain_label=-1;
  resnames.clear();
  resnames.shrink_to_fit();
  residues.clear();
  residues.shrink_to_fit();
}

void Chain::set_chain(const vector<string> &Vseq_one, int label_one, string chainid_one, int resid_offset) {
  chain_label=label_one;
  chainid=chainid_one;
  resnames.clear();
  residues.clear();

  for(int ii=0; ii<(int)Vseq_one.size();ii++) {
    int resid_one=(ii+resid_offset)%RESID_MAX+1;
    string resname_one=Vseq_one[ii];
    Residue res_one(cbs);
    res_one.set_residue(chainid,resid_one,resname_one);
    residues.push_back(res_one);
    resnames.push_back(resname_one);
  }
}

void Chain::set_chain(const vector<string> &Vresname_one, int label_one) {
  chain_label=label_one;
  resnames.clear();
  residues.clear();

  for(size_t ii=0;ii<Vresname_one.size();ii++) {
    vector<string> stoken=split_string(Vresname_one[ii],"/");
    if(stoken.size()!=3) {
      string str="chain.cpp: Wrong residue name(chainid/resid/resname) in set_chain--> "+Vresname_one[ii];
      error->one(str);
    }

    string chainid_one=stoken[0];
    int resid_one=stoi(stoken[1]);
    if(resid_one>RESID_MAX)
      resid_one=(resid_one-1)%RESID_MAX+1;
    string resname_one=stoken[2];

    if(ii==0) chainid=chainid_one;
    if(chainid_one!=chainid) {
      string str="chain.cpp: Chain ids are not matched in set_chain--> "+chainid+" : "+chainid_one;
      error->one(str);
    }

    Residue res_one(cbs,chainid_one,resid_one,resname_one);
    residues.push_back(res_one);
    resnames.push_back(Vresname_one[ii]);
  }
}

void Chain::set_chain(const vector<Residue> &Vres, int label_one) {
  chainid=Vres[0].chainid;
  chain_label=label_one;
  residues.clear();
  resnames.clear();
  for(size_t ii=0;ii<Vres.size();ii++) {
    if(Vres[ii].chainid!=chainid) {
      string str="chain.cpp: Failed to set chain --> different chain id "+chainid+"/"+Vres[ii].chainid;
      error->one(str);
    }
    residues.push_back(Vres[ii]);
    resnames.push_back(Vres[ii].resname);
  }
}

void Chain::set_chain(const Chain &chain_one) {
  Chain::clear();
  chainid=chain_one.chainid;
  resnames=chain_one.resnames;
  for(size_t ii=0;ii<chain_one.residues.size();ii++)
    residues.push_back(chain_one.residues[ii]);
  chain_label=chain_one.chain_label;
}

void Chain::standardize_chain() {
  for(size_t ii=0;ii<residues.size();ii++) {
    int resid_one=ii+1;
    if(resid_one>RESID_MAX)
      resid_one=(resid_one-1)%RESID_MAX+1;
    string resname_one=residues[ii].resname;
    residues[ii].reset_residue_name(chainid,resid_one,resname_one);
  }
}

void Chain::standardize_chain(string chainid_one, int resid_offset) {
  chainid=chainid_one;
  for(size_t ii=0;ii<residues.size();ii++) {
    int resid_one=(ii+resid_offset)+1;
    if(resid_one>RESID_MAX)
      resid_one=(resid_one-1)%RESID_MAX+1;
    string resname_one=residues[ii].resname;
    residues[ii].reset_residue_name(chainid_one,resid_one,resname_one);
  }
}

void Chain::simplify_sidechain() {
  if(chain_label==Chain::rna) {
    for(size_t ii=0;ii<residues.size();ii++) 
      residues[ii].simplify_RNA_base();
  }
}

void Chain::recovery_sidechain() {
  if(chain_label==Chain::rna) {
    for(size_t ii=0;ii<residues.size();ii++)
      residues[ii].recovery_RNA_base();
  }
}

void Chain::replace_sequence(vector<string> &Vseq_one) {
  if(residues.size()!=Vseq_one.size()) {
    string str="chain.cpp: Number of residues is not matched for sequence replacement --> ";
    str+=residues.size()+"/"+Vseq_one.size();
    error->one(str);
  }

  resnames.clear();
  for(size_t ii=0;ii<Vseq_one.size(); ii++)
    resnames.push_back(Vseq_one[ii]);
  
  if(chain_label==Chain::rna) {
    for(size_t ii=0;ii<residues.size();ii++) {
      if(residues[ii].resname==Vseq_one[ii]) continue;
      residues[ii].simplify_RNA_base();
      residues[ii].reset_residue_name(chainid,residues[ii].resid,Vseq_one[ii]);
      residues[ii].recovery_RNA_base();
    }
  }
  else {
    for(size_t ii=0;ii<residues.size();ii++) {
      if(residues[ii].resname==Vseq_one[ii]) continue;
      residues[ii].reset_residue_name(chainid,residues[ii].resid,Vseq_one[ii]);
    }
  }
}

vector<string> Chain::get_residue_name() const{
  vector<string> Vres;
  for(size_t ii=0;ii<residues.size();ii++) {
    string res_name=residues[ii].get_residue_name();
    Vres.push_back(res_name);
  }

  return Vres;
}

void Chain::reset_chain_name(string cid_one) {
  chainid=cid_one;
  for(size_t ii=0;ii<residues.size();ii++) {
    residues[ii].reset_residue_name(cid_one,residues[ii].resid,residues[ii].resname);
  }
} 

void Chain::reset_residues_name(string cid_one, int resid_begin) {
  chainid=cid_one;
  for(size_t ii=0;ii<residues.size();ii++) {
    residues[ii].reset_residue_name(cid_one,resid_begin+ii,residues[ii].resname);
  }
}

void Chain::change_chain_name(vector<string> &Vold, vector<string> &Vnew) {
  int resid_begin=residues[0].resid;
  int resid_end=residues.back().resid;
  int resid_begin_new=residues[0].resid;
  chainid=Vnew[0];
  if(Vold.size()==3) {
    resid_begin=stoi(Vold[1]);
    resid_end=stoi(Vold[2]);
    resid_begin_new=stoi(Vnew[1]);
  }
  int check=0;
  int num=0;
  for(size_t ii=0;ii<residues.size();ii++) {
    if(residues[ii].resid==resid_begin) check=1;
    else if(residues[ii].resid==resid_end) check=2;
    
    if(check==1 || residues[ii].resid==resid_end) {
      residues[ii].reset_residue_name(chainid,resid_begin_new+num,residues[ii].resname);
      num++;
    }
  }
}

vector<double> Chain::get_chain_xyz() {
  vector<double> Vxyz_c;
  for(size_t ii=0;ii<residues.size();ii++) {
    vector<double> Vxyz_res;
    if(chain_label==Chain::rna)
      Vxyz_res=residues[ii].get_RNA_residue_xyz();
    else
      Vxyz_res=residues[ii].get_residue_xyz();
    Vxyz_c.insert(Vxyz_c.end(),Vxyz_res.begin(),Vxyz_res.end());
  }

  return Vxyz_c;
}

vector<double> Chain::get_chain_xyz() const {
  vector<double> Vxyz_c;
  for(size_t ii=0;ii<residues.size();ii++) {
    vector<double> Vxyz_res=residues[ii].get_residue_xyz();
    Vxyz_c.insert(Vxyz_c.end(),Vxyz_res.begin(),Vxyz_res.end());
  }

  return Vxyz_c;
}

void Chain::reset_chain_xyz(const vector<double> &Vxyz) {
  int n0=0;
  for(size_t ii=0;ii<residues.size();ii++) {
    vector<double> Vxyz_res;
    int nxyz=3*residues[ii].size();
    if(n0+nxyz>(int)Vxyz.size()) {
      string str="chain.cpp: Size of Vxyz is not matched in setting of xyz --> ";
      str+=to_string(n0+nxyz)+"/"+to_string(Vxyz.size());
      error->one(str);
    }
    Vxyz_res.insert(Vxyz_res.end(),Vxyz.begin()+n0,Vxyz.begin()+n0+nxyz);
    residues[ii].reset_residue_xyz(Vxyz_res);
    n0+=nxyz;
  }

  if(n0!=(int)Vxyz.size()) {
    string str="chain.cpp: Size of Vxyz is not matched in setting of xyz --> ";
    str+=to_string(n0)+"/"+to_string(Vxyz.size());
    error->one(str);
  }
}

void Chain::set_chain_xyz(Chain &chain_one) {
  if(residues.size()!=chain_one.size()) {
    string str="chain.cpp: Number of residues is not matched in setting of xyz --> ";
    str+=to_string(residues.size())+"/"+to_string(chain_one.size());
    error->one(str);
  }
  if(chain_label!=chain_one.chain_label) {
    string str="chain.cpp: Chain labels are not matched in setting of xyz --> ";
    str+=to_string(chain_label)+"/"+to_string(chain_one.chain_label);
    error->one(str);
  }
  for(size_t ii=0;ii<residues.size();ii++) {
    if(chain_label==Chain::rna) 
      residues[ii].set_RNA_residue_xyz(chain_one.residues[ii]);
    else {
      residues[ii].set_residue_xyz(chain_one.residues[ii]);
    }
  }
}

void Chain::move_chain(double dx, double dy, double dz) {
  for(size_t ii=0;ii<residues.size();ii++) {
    for(size_t jj=0;jj<residues[ii].atoms.size();jj++) {
      residues[ii].atoms[jj].x+=dx;
      residues[ii].atoms[jj].y+=dy;
      residues[ii].atoms[jj].z+=dz;
    }
  }
}

void Chain::move_chain(const vector<double> & Vxyz) {
  if(Vxyz.size()!=3) {
    string str="chain.cpp: Vxyz.size()!=3 in moving chain --> "+to_string(Vxyz.size());
    error->one(str);
  }

  for(size_t ii=0;ii<residues.size();ii++) {
    for(size_t jj=0;jj<residues[ii].atoms.size();jj++) {
      residues[ii].atoms[jj].x+=Vxyz[0];
      residues[ii].atoms[jj].y+=Vxyz[1];
      residues[ii].atoms[jj].z+=Vxyz[2];
    }
  }
}

void Chain::insert_residue(Residue &res_one, int pos) {
  string cid_one=res_one.chainid;
  string resname_one=res_one.resname;

  if(residues.size()==0) {
    chainid=cid_one;
    resnames.push_back(resname_one);
    residues.push_back(res_one);
  }
  else {
    if(chainid != cid_one) {
      string str="chain.cpp: Cannot insert residue --> different chain name "+chainid+"/"+cid_one;
      error->one(str);
    }

    if(pos==Chain::begin) {
      resnames.insert(resnames.begin(),resname_one);
      residues.insert(residues.begin(),res_one);
    }
    else if(pos==Chain::end) {
      resnames.push_back(resname_one);
      residues.push_back(res_one);
    }
    else {
      string str="chain.cpp: pos=0/1 only for insert_residue";
      error->one(str);
    }
  }
}

void Chain::insert_residue(string res_one, int pos) {
  Residue res(cbs);
  res.set_residue(res_one);

  if(residues.size()==0) {
    chainid=res.chainid;
    resnames.push_back(res.resname);
    residues.push_back(res);
  }
  else {
    if(chainid != res.chainid) {
      string str="chain.cpp: Cannot insert residue --> different chain name "+chainid+"/"+res.chainid;
      error->one(str);
    }

    if(pos==Chain::begin) {
      resnames.insert(resnames.begin(),res.resname);
      residues.insert(residues.begin(),res);
    }
    else if(pos==Chain::end) {
      resnames.push_back(res.resname);
      residues.push_back(res);
    }
    else {
      string str="chain.cpp: pos=0/1 only for insert_residue";
      error->one(str);
    }
  }
}

void Chain::delete_residue(int pos) {
  if(residues.size()==0) return;
  if(pos==Chain::begin) {
    resnames.erase(resnames.begin());
    residues.erase(residues.begin());
  }
  else if(pos==Chain::end) {
    resnames.pop_back();
    residues.pop_back();
  }
  else {
    string str="chain.cpp: pos=0/1 only for delete_residue";
    error->one(str);
  }
}

void Chain::delete_residue(string res_one) {
  int check=0;
  for(size_t ii=0;ii<residues.size();ii++) {
    if(residues[ii].get_residue_name()==res_one) {
      resnames.erase(resnames.begin()+ii);
      residues.erase(residues.begin()+ii);
      check=1;
      break;
    }
  }

  if(!check) {
    string str="chain.cpp: Failed to delete residue "+res_one;
    error->warning(str);
  }
}

string Chain::first_of_residue() const {
  return residues[0].get_residue_name();
}

string Chain::last_of_residue() const {
  return residues.back().get_residue_name();
}

void Chain::clear() {
  chainid="";
  chain_label=-1;
  resnames.clear();
  resnames.shrink_to_fit();
  residues.clear();
  residues.shrink_to_fit();
}

size_t Chain::size() const{
  return residues.size();
}

size_t Chain::number_of_atoms() const {
  size_t natoms=0;
  for(size_t ii=0;ii<residues.size();ii++) {
    natoms+=residues[ii].atoms.size();
  }

  return natoms;
}

vector<double> Chain::center_of_geometry() const {
  vector<double> Vxyz;
  for(size_t ii=0;ii<3;ii++)
    Vxyz.push_back(0.0);
  
  size_t natoms=0;
  for(size_t ii=0;ii<residues.size();ii++) {
    natoms+=residues[ii].atoms.size();
    for(size_t jj=0;jj<residues[ii].atoms.size();jj++) {
      Vxyz[0]+=residues[ii].atoms[jj].x;
      Vxyz[1]+=residues[ii].atoms[jj].y;
      Vxyz[2]+=residues[ii].atoms[jj].z;
    }
  }
  if(natoms>0) {
    Vxyz[0]/=natoms;
    Vxyz[1]/=natoms;
    Vxyz[2]/=natoms;
  }

  return Vxyz;
}

Residue & Chain::operator[](int index_one) {
  if(index_one>=0) 
    return residues[index_one];
  else 
    return residues[residues.size()+index_one];
}

int Chain::set_atomid(int index_one) {
  if(init->noH) {
    for(size_t ii=0;ii<residues.size();ii++)
      residues[ii].delete_atom_hydrogen();
  }

  vector<string> Vatom_Pgroup;
  Vatom_Pgroup.push_back("P");
  Vatom_Pgroup.push_back("OP1");
  Vatom_Pgroup.push_back("OP2");

  if(chain_label==Chain::rna) {
    if(init->output_format=="amb") {
      residues[0].delete_atom(Vatom_Pgroup);
    }
    for(size_t ii=0;ii<residues.size();ii++) 
      index_one=residues[ii].set_RNA_atomid(index_one,init->output_format);
  }
  else if(chain_label==Chain::cg) {
    for(size_t ii=0;ii<residues.size();ii++) {
      index_one=residues[ii].set_atomid(index_one);
    }
  }
  else if(chain_label==Chain::other || chain_label==Chain::empty) {
    for(size_t ii=0;ii<residues.size();ii++)
      index_one=residues[ii].set_atomid_empty(index_one);
  }
  else {
    string str="chain.cpp: Please set the chain label first!";
    error->one(str);
  }

  return index_one;
}

string Chain::get_sequence() const {
  string seq="";
  for(size_t ii=0;ii<residues.size();ii++) {
    string resname_one=residues[ii].resname;
    string sname;
    if(chain_label==Chain::rna) {
      sname="R_"+resname_one;
    }
    else if(chain_label==Chain::protein) {
      sname="P_"+resname_one;
    }
    else if(chain_label==Chain::dna) {
      sname="D_"+resname_one;
    }
    else if(chain_label==Chain::cg) {
      if(resname_one[0]=='c')
        resname_one=resname_one.substr(1);
    }
    else 
      sname="";

    if(init->SEQUENCE_MAP.count(sname))
      resname_one=init->SEQUENCE_MAP[sname];
    
    if(resname_one.size()==1) {
      seq+=resname_one;
    }
    else if(resname_one.size()>1) {
      seq+="("+resname_one+")";
    }
  }

  return seq;
}

ostream & operator <<(ostream &output, Chain &chain_one) {
  string resname_last;
  int resid_last;
  int nindex;
  for(size_t ii=0;ii<chain_one.size();ii++) {
    output<<chain_one.residues[ii];
    if(ii==chain_one.size()-1) {
      resid_last=chain_one[ii].resid;
      resname_last=chain_one[ii].resname;
      nindex=chain_one[ii].atoms.back().id+1;
    }
  }
  char str[128];
  char resseq[4]="   ";
  if(resname_last.size()>3) 
    resname_last=resname_last.substr(resname_last.size()-3);
  size_t len=resname_last.size();
  for(int ii=len-1;ii>=0;ii--)
    resseq[3-(len-ii)]=resname_last[ii];
  sprintf(str,"TER   %5d      %s %c%4d",nindex,resseq,chain_one.chainid[0],resid_last);
  output<<str<<endl;

  return output;
}