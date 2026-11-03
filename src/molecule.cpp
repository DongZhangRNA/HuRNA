#include "basic.h"
#include "constant.h"
#include "initialization.h"
#include "error.h"
#include "molecule.h"

using namespace CBS_NS;

Molecule::Molecule(Cobmms *cbs) : Pointers(cbs) {}

Molecule::Molecule(Cobmms *cbs, vector<vector<string> > &Vseq_one, int label_one) : Pointers(cbs) {
  Molecule::set_molecule(Vseq_one,label_one);
}

Molecule::Molecule(const Molecule &other) : Pointers(other.cbs) {
  Molecule::clear();
  chain_resnames=other.chain_resnames;
  chains=other.chains;
  chain_ids=other.chain_ids;
  chain_labels=other.chain_labels;
  Vresnames=other.Vresnames;
  Vchainid=other.Vchainid;
}

Molecule& Molecule::operator=(const Molecule &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    chain_resnames=other.chain_resnames;
    chains=other.chains;
    chain_ids=other.chain_ids;
    chain_labels=other.chain_labels;
    Vresnames=other.Vresnames;
    Vchainid=other.Vchainid;
  }

  return *this;
}

Molecule::~Molecule() {
  chain_resnames.clear();
  chain_resnames.shrink_to_fit();
  chains.clear();
  chains.shrink_to_fit();
  chain_ids.clear();
  chain_ids.shrink_to_fit();
  chain_labels.clear();
  chain_labels.shrink_to_fit();
  Vresnames.clear();
  Vresnames.shrink_to_fit();
  Vchainid.clear();
  Vchainid.shrink_to_fit();
}

void Molecule::set_molecule(const vector<vector<string> > &Vresname_one, int label_one) {
  chain_resnames.clear();
  chains.clear();
  chain_ids.clear();
  chain_labels.clear();
  Vresnames.clear();
  Vchainid.clear();

  if(label_one!=Chain::rna && label_one!=Chain::dna && label_one!=Chain::protein && label_one!=Chain::ligand && label_one!=Chain::cg) {
    string str="molecule.cpp: Molecule label is invalid! --> "+label_one;
    error->one(str);
  }

  for(size_t ii=0;ii<Vresname_one.size();ii++) {
    Chain chain_one(cbs);
    chain_one.set_chain(Vresname_one[ii],label_one);
    chains.push_back(chain_one);
    chain_resnames.push_back(chain_one.resnames);
    chain_ids.push_back(chain_one.chainid);
    chain_labels.push_back(label_one);
    for(size_t jj=0;jj<Vresname_one[ii].size();jj++) {
      Vchainid.push_back(ii);
      Vresnames.push_back(Vresname_one[ii][jj]);
    }
  }
}

void Molecule::set_molecule(const vector<Atom> &Vatom_one, int label_one) {
  chain_resnames.clear();
  chains.clear();
  chain_ids.clear();
  chain_labels.clear();
  Vresnames.clear();
  Vchainid.clear();

  if(label_one!=Chain::rna && label_one!=Chain::dna && label_one!=Chain::protein && label_one!=Chain::ligand 
     && label_one!=Chain::cg && label_one!=Chain::other && label_one!=Chain::empty) {
    string str="molecule.cpp: Molecule label is invalid! --> "+label_one;
    error->one(str);
  }
  
  string chainid_p=Vatom_one[0].chainid;
  string res_p=Vatom_one[0].get_residue_name();
  Chain chain_one(cbs);
  chain_one.chain_label=label_one;
  Residue res_one(cbs);
  vector<Atom> Tatoms;
  Tatoms.push_back(Vatom_one[0]);

  for(size_t ii=1;ii<Vatom_one.size();ii++) {
    string res_i=Vatom_one[ii].get_residue_name();
    if(res_i!=res_p) {
      res_one.clear();
      res_one.set_residue(Tatoms);
      chain_one.insert_residue(res_one,Chain::end);
      Tatoms.clear();
    }
    if(Vatom_one[ii].chainid!=chainid_p) {
      if(Tatoms.size()>0) {
        res_one.clear();
        res_one.set_residue(Tatoms);
        chain_one.insert_residue(res_one,Chain::end);
        Tatoms.clear();
      }
      if(chain_one.size()>0)  
        Molecule::insert_chain(chain_one);
      chain_one.clear();
      chain_one.chain_label=label_one;
    }
    Tatoms.push_back(Vatom_one[ii]);
    chainid_p=Vatom_one[ii].chainid;
    res_p=res_i;
  }
  if(Tatoms.size()>0) {
    res_one.clear();
    res_one.set_residue(Tatoms);
    chain_one.insert_residue(res_one,Chain::end);
    Molecule::insert_chain(chain_one);
    Tatoms.clear();
    chain_one.clear();
  }
}

void Molecule::set_molecule(const vector<Chain> &Vchains_one) {
  Molecule::clear();
  for(size_t ii=0;ii<Vchains_one.size();ii++) {
    Molecule::insert_chain(Vchains_one[ii]);
  }
}

void Molecule::set_molecule(const Molecule &molecule_one) {
  Molecule::clear();
  chain_resnames=molecule_one.chain_resnames;
  for(size_t ii=0;ii<molecule_one.chains.size();ii++)
    chains.push_back(molecule_one.chains[ii]);
  chain_ids=molecule_one.chain_ids;
  chain_labels=molecule_one.chain_labels;
  Vresnames=molecule_one.Vresnames;
  Vchainid=molecule_one.Vchainid;
}

Molecule Molecule::copy_molecule() const{
  Molecule mol_one(cbs);
  mol_one.chain_resnames=chain_resnames;
  for(size_t ii=0;ii<chains.size();ii++)
    mol_one.chains.push_back(chains[ii]);
  mol_one.chain_ids=chain_ids;
  mol_one.chain_labels=chain_labels;
  mol_one.Vresnames=Vresnames;
  mol_one.Vchainid=Vchainid;

  return mol_one;
}

void Molecule::standardize_molecule(int cid_offset) {
  Vresnames.clear();
  Vchainid.clear();
  
  for(int ii=0;ii<(int)chains.size();ii++) {
    int cid=(ii+cid_offset)%((int)init->CHAIN_ID_LIST.size());
    int nround=(ii+cid_offset)/((int)init->CHAIN_ID_LIST.size());
    string chainid_one=init->CHAIN_ID_LIST[cid];
    if(nround>0) chainid_one+=init->CHAIN_ID_LIST[nround-1];
    chains[ii].standardize_chain(chainid_one);
    for(size_t jj=0;jj<chains[ii].size();jj++) {
      string res_name=chains[ii][jj].get_residue_name();
      Vchainid.push_back(ii);
      Vresnames.push_back(res_name);
    }
  }
}

void Molecule::simplify_sidechain() {
  for(size_t ii=0;ii<chains.size();ii++)
    chains[ii].simplify_sidechain();
}

void Molecule::recovery_sidechain() {
  for(size_t ii=0;ii<chains.size();ii++)
    chains[ii].recovery_sidechain();
}

void Molecule::replace_sequence(vector<vector<string> > &Vseq_one) {
  if(Vseq_one.size()!=chains.size()) {
    string str="molecule.cpp: the new sequences do not match the number of chains in the replaced molecule --> ";
    str+=to_string((int) Vseq_one.size())+"/"+to_string((int) chains.size());
    error->one(str);
  }

  vector<vector<string> > Tresnames=chain_resnames;
  chain_resnames.clear();
  for(size_t ii=0;ii<Vseq_one.size();ii++) {
    if(Vseq_one[ii].size()==1 && Vseq_one[ii][0]=="N") {
      chain_resnames.push_back(Tresnames[ii]);
    }
    else {
      vector<string> Vres;
      for(size_t jj=0;jj<Vseq_one[ii].size();jj++)
        Vres.push_back(Vseq_one[ii][jj]);
      chain_resnames.push_back(Vres);
    }
  }

  Vresnames.clear();
  for(size_t ii=0;ii<chains.size();ii++) {
    chains[ii].replace_sequence(chain_resnames[ii]);
    for(size_t jj=0;jj<chains[ii].size();jj++) {
      string res_name=chains[ii][jj].get_residue_name();
      Vresnames.push_back(res_name);
    }
  }
}

void Molecule::update_residue_name() {
  Vresnames.clear();
  for(size_t ii=0;ii<chains.size();ii++) {
    vector<string> Tres=chains[ii].get_residue_name();
    Vresnames.insert(Vresnames.end(),Tres.begin(),Tres.end());
  }
}

vector<vector<string> > Molecule::get_residue_name() const{
  vector<vector<string> > Vres;
  for(size_t ii=0;ii<chains.size();ii++) {
    vector<string> Tres=chains[ii].get_residue_name();
    Vres.push_back(Tres);
  }

  return Vres;
}

void Molecule::change_chain_name(vector<vector<string> > &VVolds, vector<vector<string> > &VVnews) {
  for(size_t ii=0;ii<chains.size();ii++) {
    for(size_t jj=0;jj<VVolds.size();jj++) {
      if(chains[ii].chainid==VVolds[jj][0])
        chains[ii].change_chain_name(VVolds[jj],VVnews[jj]);
    }
  }

  Vresnames.clear();
  for(size_t ii=0;ii<chains.size();ii++)
    for(size_t jj=0;jj<chains[ii].size();jj++)
      Vresnames.push_back(chains[ii][jj].get_residue_name());
}

vector<double> Molecule::get_molecule_xyz() {
  vector<double> Vxyz_m;
  for(size_t ii=0;ii<chains.size();ii++) {
    vector<double> Vxyz_c=chains[ii].get_chain_xyz();
    Vxyz_m.insert(Vxyz_m.end(),Vxyz_c.begin(),Vxyz_c.end());
  }

  return Vxyz_m;
}

vector<double> Molecule::get_molecule_xyz() const {
  vector<double> Vxyz_m;
  for(size_t ii=0;ii<chains.size();ii++) {
    vector<double> Vxyz_c=chains[ii].get_chain_xyz();
    Vxyz_m.insert(Vxyz_m.end(),Vxyz_c.begin(),Vxyz_c.end());
  }

  return Vxyz_m;
}

void Molecule::reset_molecule_xyz(const vector<double> &Vxyz) {
  int n0=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    vector<double> Vxyz_chain;
    for(size_t jj=0;jj<chains[ii].size();jj++) {
      int nxyz=3*chains[ii][jj].size();
      if(n0+nxyz>(int)Vxyz.size()) {
        string str="molecule.cpp: Size of Vxyz is not matched in setting of xyz --> ";
        str+=to_string(n0+nxyz)+"/"+to_string(Vxyz.size());
        error->one(str);
      }
      Vxyz_chain.insert(Vxyz_chain.end(),Vxyz.begin()+n0,Vxyz.begin()+n0+nxyz);
      n0+=nxyz;
    }
    chains[ii].reset_chain_xyz(Vxyz_chain);
  }

  if(n0!=(int)Vxyz.size()) {
    string str="molecule.cpp: Size of Vxyz is not matched in setting of xyz --> ";
    str+=to_string(n0)+"/"+to_string(Vxyz.size());
    error->one(str);
  }
}

void Molecule::reset_molecule_xyz(const vector<Atom> &Vatoms_one) {
  vector<vector<Atom> > Vresidue_atoms;
  string res_p=Vatoms_one[0].get_residue_name();
  vector<Atom> Tatoms;
  Tatoms.push_back(Vatoms_one[0]);
  for(size_t ii=1;ii<Vatoms_one.size();ii++) {
    string res_i=Vatoms_one[ii].get_residue_name();
    if(res_i!=res_p) {
      Vresidue_atoms.push_back(Tatoms);
      Tatoms.clear();
    }
    Tatoms.push_back(Vatoms_one[ii]);
    res_p=res_i;
  }
  if(Tatoms.size()>0) {
    Vresidue_atoms.push_back(Tatoms);
    Tatoms.clear();
  }
  
  if(Vresidue_atoms.size()!=Molecule::number_of_residues()) {
    string str="molecule.cpp: Number of residues is not matched in setting of xyz --> ";
    str+=to_string(Vresidue_atoms.size())+"/"+to_string(Molecule::number_of_residues());
    error->one(str);
  }

  size_t nres=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    for(size_t jj=0;jj<chains[ii].size();jj++) {
      chains[ii].residues[jj].reset_residue_xyz(Vresidue_atoms[nres]);
      nres++;
    }
  }
  
  /*
  map<string,vector<double> > atomname_xyz;
  for(size_t ii=0;ii<Vatoms_one.size();ii++) {
    string atomname_one=Vatoms_one[ii].get_atom_name();
    vector<double> xyz_one=Vatoms_one[ii].get_atom_xyz();
    if(atomname_xyz.count(atomname_one)==0) {
      atomname_xyz.emplace(atomname_one,xyz_one);
    }
  }
  for(size_t ii=0;ii<chains.size();ii++) {
    for(size_t jj=0;jj<chains[ii].residues.size();jj++) {
      for(size_t kk=0;kk<chains[ii].residues[jj].atoms.size();kk++) {
        string atomname_one=chains[ii].residues[jj].atoms[kk].get_atom_name();
        if(atomname_xyz.count(atomname_one)) {
          chains[ii].residues[jj].atoms[kk].reset_atom_xyz(atomname_xyz[atomname_one]);
        }
        else {
          string str="molecule.cpp: Cannot find xyz for atom "+atomname_one;
          error->one(str);
        }
      }
    }
  }
  */
}

void Molecule::set_molecule_xyz(Molecule &mole_one) {
  if(chains.size()!=mole_one.size()) {
    string str="molecule.cpp: Number of chains is not matched in setting of xyz --> ";
    str+=to_string(chains.size())+"/"+to_string(mole_one.size());
    error->one(str);
  }
  for(size_t ii=0;ii<chains.size();ii++) {
    chains[ii].set_chain_xyz(mole_one.chains[ii]);
  }
}

void Molecule::move_molecule(double delx, double dely, double delz) {
  for(size_t ii=0;ii<chains.size();ii++) {
    for(size_t jj=0;jj<chains[ii].residues.size();jj++) {
      for(size_t kk=0;kk<chains[ii].residues[jj].atoms.size();kk++) {
        chains[ii].residues[jj].atoms[kk].x+=delx;
        chains[ii].residues[jj].atoms[kk].y+=dely;
        chains[ii].residues[jj].atoms[kk].z+=delz;
      }
    }
  }
}

void Molecule::move_molecule(const vector<double> & Vxyz) {
  if(Vxyz.size()!=3) {
    string str="molecule.cpp: Vxyz.size()!=3 in moving molecule --> "+to_string(Vxyz.size());
    error->one(str);
  }

  for(size_t ii=0;ii<chains.size();ii++) {
    for(size_t jj=0;jj<chains[ii].residues.size();jj++) {
      for(size_t kk=0;kk<chains[ii].residues[jj].atoms.size();kk++) {
        chains[ii].residues[jj].atoms[kk].x+=Vxyz[0];
        chains[ii].residues[jj].atoms[kk].y+=Vxyz[1];
        chains[ii].residues[jj].atoms[kk].z+=Vxyz[2];
      }
    }
  }
}

string Molecule::previous_residue(string res_one) const {
  vector<string> stoken=split_string(res_one,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+res_one;
    error->one(str);
  }
  
  string res_p=ResidueConst::RESIDUE_NULL;
  int ni=-1;
  for(size_t ii=0;ii<Vresnames.size();ii++) {
    if(Vresnames[ii]==res_one)
      {ni=ii; break;}
  }
  if(ni<0) {
    string str="molecule.cpp: Cannot find residue "+res_one;
    error->warning(str);
  }
  else if(ni>0) {
    res_p=Vresnames[ni-1];
    if(Vchainid[ni]!=Vchainid[ni-1]) res_p=ResidueConst::RESIDUE_NULL;
  }

  return res_p;
}

string Molecule::next_residue(string res_one) const{
  vector<string> stoken=split_string(res_one,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+res_one;
    error->one(str);
  }
  
  string res_n=ResidueConst::RESIDUE_NULL;
  int ni=-1;
  for(size_t ii=0;ii<Vresnames.size();ii++) {
    if(Vresnames[ii]==res_one)  
      {ni=ii; break;}
  }
  if(ni<0) {
    string str="molecule.cpp: Cannot find residue "+res_one;
    error->warning(str);
  }
  else if(ni<(int)Vresnames.size()-1) {
    res_n=Vresnames[ni+1];
    if(Vchainid[ni]!=Vchainid[ni+1]) res_n=ResidueConst::RESIDUE_NULL;
  }

  return res_n;
}

Residue Molecule::find_residue(string res_one) const {
  vector<string> stoken=split_string(res_one,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+res_one;
    error->one(str);
  }
  string chainid_one=stoken[0];
  int resid_one=stoi(stoken[1]);
  string resname_one=stoken[2];
  int nc=count(Vresnames.begin(),Vresnames.end(),res_one);
  if(nc<1) {
    string str="molecule.cpp: Cannot find residue "+res_one;
    error->one(str);
  }
  else if(nc>1) {
    string str="molecule.cpp: More than one residue was found for "+res_one+" --> "+to_string(nc);
    error->one(str);
  }

  size_t index_i=0;
  size_t index_j=0;
  int flag_check=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    if(chains[ii].chainid!=chainid_one) continue;
    index_i=ii;
    for(size_t jj=0;jj<chains[ii].residues.size();jj++) {
      if(chains[ii].residues[jj].resid==resid_one && chains[ii].residues[jj].resname==resname_one) {
        index_j=jj;
        flag_check=1;
        break;
      }
    }
    if(flag_check) break;
  }

  if(flag_check==0) {
    string str="molecule.cpp: Cannot find residue "+res_one+" in molecule";
    error->one(str);
  }
  
  string res_check=chains[index_i].residues[index_j].get_residue_name();
  if(res_one!=res_check) {
    string str="molecule.cpp: find residue error --> "+res_one+" : "+res_check;
    error->one(str);
  }

  return chains[index_i].residues[index_j];
}

vector<Atom> Molecule::get_atom_list() const {
  vector<Atom> Vatoms_list;
  for(size_t ii=0;ii<chains.size();ii++) {
    for(size_t jj=0;jj<chains[ii].residues.size();jj++) {
      for(size_t kk=0;kk<chains[ii].residues[jj].atoms.size();kk++) {
        Vatoms_list.push_back(chains[ii].residues[jj].atoms[kk]);
      }
    }
  }

  return Vatoms_list;
}

Atom Molecule::find_atom(string atom_one) const {
  vector<string> stoken=split_string(atom_one,"/");
  if(stoken.size()!=4) {
    string str="molecule.cpp: Wrong atom name(chainid/resid/resname/atomname) --> "+atom_one;
    error->one(str);
  }
  string chainid_one=stoken[0];
  int resid_one=stoi(stoken[1]);
  string resname_one=stoken[2];
  string atomname_one=stoken[4];
  
  string res_one=stoken[0]+"/"+stoken[1]+"/"+stoken[2];
  int nc=count(Vresnames.begin(),Vresnames.end(),res_one);
  if(nc<1) {
    string str="molecule.cpp: Cannot find residue "+res_one;
    error->one(str);
  }
  else if(nc>1) {
    string str="molecule.cpp: More than one residue was found for "+res_one+" --> "+to_string(nc);
    error->one(str);
  }

  size_t index_i=0;
  size_t index_j=0;
  int flag_check=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    if(chains[ii].chainid!=chainid_one) continue;
    index_i=ii;
    for(size_t jj=0;jj<chains[ii].residues.size();jj++) {
      if(chains[ii].residues[jj].resid==resid_one && chains[ii].residues[jj].resname==resname_one) {
        index_j=jj;
        flag_check=1;
        break;
      }
    }
    if(flag_check) break;
  }
  
  string res_check=chains[index_i].residues[index_j].get_residue_name();
  if(res_one!=res_check) {
    string str="molecule.cpp: find residue error --> "+res_one+" : "+res_check;
    error->one(str);
  }

  Residue residue_one=chains[index_i].residues[index_j];
  Atom atom_find=residue_one.find_atom(atomname_one);

  return atom_find;
}

int Molecule::find_atomid(string atom_one) const {
  vector<string> stoken=split_string(atom_one,"/");
  if(stoken.size()!=4) {
    string str="molecule.cpp: Wrong atom name(chainid/resid/resname/atomname) --> "+atom_one;
    error->one(str);
  }
  string chainid_one=stoken[0];
  int resid_one=stoi(stoken[1]);
  string resname_one=stoken[2];
  string atomname_one=stoken[4];
  
  string res_one=stoken[0]+"/"+stoken[1]+"/"+stoken[2];
  int nc=count(Vresnames.begin(),Vresnames.end(),res_one);
  if(nc<1) {
    string str="molecule.cpp: Cannot find residue "+res_one;
    error->one(str);
  }
  else if(nc>1) {
    string str="molecule.cpp: More than one residue was found for "+res_one+" --> "+to_string(nc);
    error->one(str);
  }

  size_t index_i=0;
  size_t index_j=0;
  int flag_check=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    if(chains[ii].chainid!=chainid_one) continue;
    index_i=ii;
    for(size_t jj=0;jj<chains[ii].residues.size();jj++) {
      if(chains[ii].residues[jj].resid==resid_one && chains[ii].residues[jj].resname==resname_one) {
        index_j=jj;
        flag_check=1;
        break;
      }
    }
    if(flag_check) break;
  }
  
  string res_check=chains[index_i].residues[index_j].get_residue_name();
  if(res_one!=res_check) {
    string str="molecule.cpp: find residue error --> "+res_one+" : "+res_check;
    error->one(str);
  }

  Residue residue_one=chains[index_i].residues[index_j];
  Atom atom_find=residue_one.find_atom(atomname_one);

  return atom_find.id;
}

string Molecule::get_atomid_list(string residue_range) const {  // residue_range=resname_begin-resname_end, return atomid_list=id_begin:id_end
  string chainid_one;
  int resid_begin,resid_end;
  /*
  vector<string> stoken=split_string(resid_range,":");
  size_t pos=resid_range.find_first_of(":");
  if(pos==0) {
    chainid_one="";
    if(stoken.size()!=2) {
      string str="molecule.cpp: Wrong resid range (chainid:resid_begin:resid_end) --> "+resid_range;
      error->one(str);
    }
    resid_begin=stoi(stoken[0]);
    resid_end=stoi(stoken[1]);
  }
  else if(pos>0 && stoken.size()==3) {
    chainid_one=stoken[0];
    resid_begin=stoi(stoken[1]);
    resid_end=stoi(stoken[2]);
  }
  else {
    string str="molecule.cpp: Wrong resid range (chainid:resid_begin:resid_end) --> "+resid_range;
    error->one(str);
  }
  */
  vector<string> stoken=split_string(residue_range,"-");
  if(stoken.size()!=2) {
    string str="molecule.cpp: Wrong residue range (residue_begin-residue_end) --> "+residue_range;
    error->one(str);
  }
  string resname_begin=stoken[0];
  string resname_end=stoken[1];
  stoken.clear();
  stoken=split_string(resname_begin,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_begin;
    error->one(str);
  }
  chainid_one=stoken[0];
  resid_begin=stoi(stoken[1]);
  stoken.clear();
  stoken=split_string(resname_end,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_end;
    error->one(str);
  }
  if(stoken[0]!=chainid_one) {
    string str="molecule.cpp: Chain id should be the same --> "+residue_range;
    error->one(str);
  }
  resid_end=stoi(stoken[1]);

  int atomid_begin=0,atomid_end=0;
  int flag_check1=0,flag_check2=0;
  size_t index_chain=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    if(chainid_one==chains[ii].chainid || (chainid_one.size()==0 && chains[ii].chainid.size()==0)) {
      index_chain=ii;
      for(size_t jj=0;jj<chains[ii].residues.size();jj++) {
        int resid_one=chains[ii].residues[jj].resid;
        if(resid_one==resid_begin) {
          flag_check1=1;
          atomid_begin=chains[ii].residues[jj].atoms[0].id;
          for(size_t kk=1;kk<chains[ii].residues[jj].atoms.size();kk++) {
            int atomid_one=chains[ii].residues[jj].atoms[kk].id;
            if(atomid_begin>atomid_one) 
              atomid_begin=atomid_one;
          }
        }
        if(resid_one==resid_end) {
          flag_check2=1;
          atomid_end=chains[ii].residues[jj].atoms[0].id;
          for(size_t kk=1;kk<chains[ii].residues[jj].atoms.size();kk++) {
            int atomid_one=chains[ii].residues[jj].atoms[kk].id;
            if(atomid_end<atomid_one)
              atomid_end=atomid_one;
          }
        }
      }
    }
  }

  string atomid_list;
  if(flag_check1 && flag_check2) {
    if(atomid_begin<=atomid_end)
      atomid_list=to_string(atomid_begin)+":"+to_string(atomid_end);
    else if(atomid_begin>atomid_end && count(init->CircRNA_chains.begin(),init->CircRNA_chains.end(),chainid_one)) {
      int chain_atomid_begin;
      for(size_t kk=0;kk<chains[index_chain].residues[0].atoms.size();kk++) {
        int id_one=chains[index_chain].residues[0].atoms[kk].id;
        if(kk==0) chain_atomid_begin=id_one;
        else if(id_one<chain_atomid_begin) chain_atomid_begin=id_one;
      }
      int chain_atomid_end;
      for(size_t kk=0;kk<chains[index_chain].residues.back().atoms.size();kk++) {
        int id_one=chains[index_chain].residues.back().atoms[kk].id;
        if(kk==0) chain_atomid_end=id_one;
        else if(id_one>chain_atomid_end) chain_atomid_end=id_one;
      }
      atomid_list=to_string(chain_atomid_begin)+":"+to_string(atomid_end)+" "+to_string(atomid_begin)+":"+to_string(chain_atomid_end);
    }
    else {
      string str="molecule.cpp: Wrong atomid list for residue range "+resname_begin+"-" +resname_end;
      error->one(str);
    }
  }
  else {
    string str="molecule.cpp: Cannot get atomid list for residue range "+residue_range;
    error->one(str);
  }

  return atomid_list;
}

string Molecule::get_atomid_list(string resname_begin, string resname_end) const {  // residue_range=resname_begin-resname_end, return atomid_list=id_begin:id_end
  string chainid_one;
  int resid_begin,resid_end;
  vector<string> stoken=split_string(resname_begin,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_begin;
    error->one(str);
  }
  chainid_one=stoken[0];
  resid_begin=stoi(stoken[1]);
  stoken.clear();
  stoken=split_string(resname_end,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_end;
    error->one(str);
  }
  if(stoken[0]!=chainid_one) {
    string str="molecule.cpp: Chain id should be the same --> "+resname_begin+" / " +resname_end;
    error->one(str);
  }
  resid_end=stoi(stoken[1]);

  int atomid_begin=0,atomid_end=0;
  int flag_check1=0,flag_check2=0;
  size_t index_chain=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    if(chainid_one==chains[ii].chainid || (chainid_one.size()==0 && chains[ii].chainid.size()==0)) {
      index_chain=ii;
      for(size_t jj=0;jj<chains[ii].residues.size();jj++) {
        int resid_one=chains[ii].residues[jj].resid;
        if(resid_one==resid_begin) {
          flag_check1=1;
          atomid_begin=chains[ii].residues[jj].atoms[0].id;
          for(size_t kk=1;kk<chains[ii].residues[jj].atoms.size();kk++) {
            int atomid_one=chains[ii].residues[jj].atoms[kk].id;
            if(atomid_begin>atomid_one) 
              atomid_begin=atomid_one;
          }
        }
        if(resid_one==resid_end) {
          flag_check2=1;
          atomid_end=chains[ii].residues[jj].atoms[0].id;
          for(size_t kk=1;kk<chains[ii].residues[jj].atoms.size();kk++) {
            int atomid_one=chains[ii].residues[jj].atoms[kk].id;
            if(atomid_end<atomid_one)
              atomid_end=atomid_one;
          }
        }
      }
    }
  }

  string atomid_list;
  if(flag_check1 && flag_check2) {
    if(atomid_begin<=atomid_end)
      atomid_list=to_string(atomid_begin)+":"+to_string(atomid_end);
    else if(atomid_begin>atomid_end && count(init->CircRNA_chains.begin(),init->CircRNA_chains.end(),chainid_one)) {
      int chain_atomid_begin;
      for(size_t kk=0;kk<chains[index_chain].residues[0].atoms.size();kk++) {
        int id_one=chains[index_chain].residues[0].atoms[kk].id;
        if(kk==0) chain_atomid_begin=id_one;
        else if(id_one<chain_atomid_begin) chain_atomid_begin=id_one;
      }
      int chain_atomid_end;
      for(size_t kk=0;kk<chains[index_chain].residues.back().atoms.size();kk++) {
        int id_one=chains[index_chain].residues.back().atoms[kk].id;
        if(kk==0) chain_atomid_end=id_one;
        else if(id_one>chain_atomid_end) chain_atomid_end=id_one;
      }
      atomid_list=to_string(chain_atomid_begin)+":"+to_string(atomid_end)+" "+to_string(atomid_begin)+":"+to_string(chain_atomid_end);
    }
    else {
      string str="molecule.cpp: Wrong atomid list for residue range "+resname_begin+"-" +resname_end;
      error->one(str);
    }
  }
  else {
    string str="molecule.cpp: Cannot get atomid list for residue range "+resname_begin+"-" +resname_end;
    error->one(str);
  }

  return atomid_list;
}

bool Molecule::is_residue_exist(string res_one) const {
  vector<string> stoken=split_string(res_one,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+res_one;
    error->one(str);
  }

  if(count(Vresnames.begin(),Vresnames.end(),res_one)) 
    return true;
  else 
    return false;
}

int Molecule::get_residue_index(string res_one, size_t index_offset) const {
  vector<string> stoken=split_string(res_one,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+res_one;
    error->one(str);
  }
  
  int ni=-1;
  for(size_t ii=index_offset; ii<Vresnames.size();ii++) {
    if(Vresnames[ii]==res_one) {
      ni=ii;
      break;
    }
  }

  if(ni<0) {
    string str="molecule.cpp: Cannot find residue "+res_one;
    error->one(str);
  }

  return ni;
}

void Molecule::insert_chain(const Chain &chain_one) {
  chain_resnames.push_back(chain_one.resnames);
  chains.push_back(chain_one);
  chain_ids.push_back(chain_one.chainid);
  chain_labels.push_back(chain_one.chain_label);
  for(size_t jj=0;jj<chain_one.size();jj++) {
    string res_name=chain_one.residues[jj].get_residue_name();
    Vresnames.push_back(res_name);
    Vchainid.push_back((int)chains.size()-1);
  }
}

void Molecule::clear() {
  chain_resnames.clear();
  chain_resnames.shrink_to_fit();
  chains.clear();
  chains.shrink_to_fit();
  chain_ids.clear();
  chain_ids.shrink_to_fit();
  chain_labels.clear();
  chain_labels.shrink_to_fit();
  Vresnames.clear();
  Vresnames.shrink_to_fit();
  Vchainid.clear();
  Vchainid.shrink_to_fit();
}

size_t Molecule::size() const {
  return chains.size();
}

vector<string> Molecule::get_residue_list(string resname_begin, string resname_end) const {
  string chainid_one;
  int resid_begin,resid_end;
  vector<string> stoken=split_string(resname_begin,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_begin;
    error->one(str);
  }
  chainid_one=stoken[0];
  resid_begin=stoi(stoken[1]);
  stoken.clear();
  stoken=split_string(resname_end,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_end;
    error->one(str);
  }
  if(stoken[0]!=chainid_one) {
    string str="molecule.cpp: Chain id should be the same --> "+resname_begin+" / " +resname_end;
    error->one(str);
  }
  resid_end=stoi(stoken[1]);

  int flag_check1=0,flag_check2=0;
  size_t index_begin=0,index_end=0, index_chain=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    if(chainid_one==chains[ii].chainid || (chainid_one.size()==0 && chains[ii].chainid.size()==0)) {
      index_chain=ii;
      for(size_t jj=0;jj<chains[ii].residues.size();jj++) {
        int resid_one=chains[ii].residues[jj].resid;
        if(resid_one==resid_begin) {
          flag_check1=1;
          index_begin=jj;
        }
        if(resid_one==resid_end) {
          flag_check2=1;
          index_end=jj;
        }
      }
    }
  }

  if(flag_check1==0) {
    string str="molecule.cpp: Cannot find residue "+resname_begin;
    error->one(str);
  }
  if(flag_check2==0) {
    string str="molecule.cpp: Cannot find residue "+resname_end;
    error->one(str);
  }

  vector<string> Vresname_one;
  if(index_begin<=index_end) {
    for(size_t jj=index_begin;jj<=index_end;jj++) {
      string resname_one=chains[index_chain].residues[jj].get_residue_name();
      Vresname_one.push_back(resname_one);
    }
  }
  else {
    if(count(init->CircRNA_chains.begin(),init->CircRNA_chains.end(),chainid_one)) {
      for(size_t jj=index_begin;jj<chains[index_chain].residues.size();jj++) {
        string resname_one=chains[index_chain].residues[jj].get_residue_name();
        Vresname_one.push_back(resname_one);
      }
      for(size_t jj=0;jj<=index_end;jj++) {
        string resname_one=chains[index_chain].residues[jj].get_residue_name();
        Vresname_one.push_back(resname_one);
      }
    }
    else {
      string str="molecule.cpp: Wrong residue range --> "+resname_begin+"-"+resname_end;
      error->one(str);
    }
  }

  return Vresname_one;
}

size_t Molecule::number_of_residues() const {
  return Vresnames.size();
}

/*
size_t Molecule::number_of_residues(string residue_range) const {  // residue_range=resname_begin-resname_end
  vector<string> stoken=split_string(residue_range,"-");
  if(stoken.size()!=2) {
    string str="molecule.cpp: Wrong residue range (residue_begin-residue_end) --> "+residue_range;
    error->one(str);
  }
  string resname_begin=stoken[0];
  string resname_end=stoken[1];

  string chainid_one;
  int resid_begin,resid_end;
  stoken.clear();
  stoken=split_string(resname_begin,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_begin;
    error->one(str);
  }
  chainid_one=stoken[0];
  resid_begin=stoi(stoken[1]);
  stoken.clear();
  stoken=split_string(resname_end,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_end;
    error->one(str);
  }
  if(stoken[0]!=chainid_one) {
    string str="molecule.cpp: Chain id should be the same --> "+residue_range;
    error->one(str);
  }
  resid_end=stoi(stoken[1]);
  
  int flag_check1=0,flag_check2=0;
  size_t index_begin=0,index_end=0, nres_chain=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    if(chainid_one==chains[ii].chainid || (chainid_one.size()==0 && chains[ii].chainid.size()==0)) {
      nres_chain=chains[ii].residues.size();
      for(size_t jj=0;jj<nres_chain;jj++) {
        int resid_one=chains[ii].residues[jj].resid;
        if(resid_one==resid_begin) {
          flag_check1=1;
          index_begin=jj;
        }
        if(resid_one==resid_end) {
          flag_check2=1;
          index_end=jj;
        }
      }
    }
  }

  if(flag_check1==0) {
    string str="molecule.cpp: Cannot find residue "+resname_begin;
    error->one(str);
  }
  if(flag_check2==0) {
    string str="molecule.cpp: Cannot find residue "+resname_end;
    error->one(str);
  }

  size_t nres=0;
  if(index_begin<=index_end)
    nres=index_end-index_begin+1;
  else {
    if(count(init->CircRNA_chains.begin(),init->CircRNA_chains.end(),chainid_one)) {
      nres=nres_chain-index_begin+index_end+1;
    }
    else {
      string str="molecule.cpp: Wrong residue range --> "+residue_range;
      error->one(str);
    }
  }

  return nres;
}
*/

size_t Molecule::number_of_residues(string resname_begin, string resname_end) const {  // residue_range=resname_begin-resname_end
  string chainid_one;
  int resid_begin,resid_end;
  vector<string> stoken=split_string(resname_begin,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_begin;
    error->one(str);
  }
  chainid_one=stoken[0];
  resid_begin=stoi(stoken[1]);
  stoken.clear();
  stoken=split_string(resname_end,"/");
  if(stoken.size()!=3) {
    string str="molecule.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_end;
    error->one(str);
  }
  if(stoken[0]!=chainid_one) {
    string str="molecule.cpp: Chain id should be the same --> "+resname_begin+" / " +resname_end;
    error->one(str);
  }
  resid_end=stoi(stoken[1]);

  int flag_check1=0,flag_check2=0;
  size_t index_begin=0,index_end=0, nres_chain=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    if(chainid_one==chains[ii].chainid || (chainid_one.size()==0 && chains[ii].chainid.size()==0)) {
      nres_chain=chains[ii].residues.size();
      for(size_t jj=0;jj<nres_chain;jj++) {
        int resid_one=chains[ii].residues[jj].resid;
        if(resid_one==resid_begin) {
          flag_check1=1;
          index_begin=jj;
        }
        if(resid_one==resid_end) {
          flag_check2=1;
          index_end=jj;
        }
      }
    }
  }

  if(flag_check1==0) {
    string str="molecule.cpp: Cannot find residue "+resname_begin;
    error->one(str);
  }
  if(flag_check2==0) {
    string str="molecule.cpp: Cannot find residue "+resname_end;
    error->one(str);
  }

  size_t nres=0;
  if(index_begin<=index_end)
    nres=index_end-index_begin+1;
  else {
    if(count(init->CircRNA_chains.begin(),init->CircRNA_chains.end(),chainid_one)) {
      nres=nres_chain-index_begin+index_end+1;
    }
    else {
      string str="molecule.cpp: Wrong residue range --> "+resname_begin+"-"+resname_end;
      error->one(str);
    }
  }

  return nres;
}

size_t Molecule::number_of_atoms() const {
  size_t natoms=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    for(size_t jj=0;jj<chains[ii].residues.size();jj++) {
      natoms+=chains[ii].residues[jj].atoms.size();
    }
  }

  return natoms;
}

vector<double> Molecule::center_of_geometry() const {
  vector<double> Vxyz;
  for(size_t ii=0;ii<3;ii++)
    Vxyz.push_back(0.0);

  size_t natoms=0;
  for(size_t ii=0;ii<chains.size();ii++) {
    for(size_t jj=0;jj<chains[ii].residues.size();jj++) {
      natoms+=chains[ii].residues[jj].atoms.size();
      for(size_t kk=0;kk<chains[ii].residues[jj].atoms.size();kk++) {
        Vxyz[0]+=chains[ii].residues[jj].atoms[kk].x;
        Vxyz[1]+=chains[ii].residues[jj].atoms[kk].y;
        Vxyz[2]+=chains[ii].residues[jj].atoms[kk].z;
      }
    }
  }
  if(natoms>0) {
    Vxyz[0]/=natoms;
    Vxyz[1]/=natoms;
    Vxyz[2]/=natoms;
  }

  return Vxyz;
}

string Molecule::get_sequence() const {
  string seq=chains[0].get_sequence();
  for(size_t ii=1;ii<chains.size();ii++)
    seq+="&"+chains[ii].get_sequence();

  return seq;
}

Chain & Molecule::operator[](int index_one) {
  if(index_one>=0)
    return chains[index_one];
  else
    return chains[chains.size()+index_one];
}

void Molecule::set_atomid(int index_offset) {
  for(size_t ii=0;ii<chains.size();ii++)
    index_offset=chains[ii].set_atomid(index_offset);
}

ostream & operator <<(ostream &output, Molecule molecule_one) {
  int index=0;
  for(size_t ii=0;ii<molecule_one.chains.size();ii++) {
    index=molecule_one.chains[ii].set_atomid(index)+1;
    output<<molecule_one.chains[ii];
  }

  return output;
}