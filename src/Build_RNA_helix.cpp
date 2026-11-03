#include "basic.h"
#include "readfile.h"
#include "error.h"
#include "initialization.h"
#include "read_PDB.h"
#include "superpose.h"
#include "RNA_molecule.h"

using namespace CBS_NS;

int main(int narg, char **arg) {
  map<string,vector<string> > Mcmd;
  string style="";
  vector<string> Vargs;
  for(int ii=1;ii<narg;ii++) {
    if(arg[ii][0]=='-') {
      if(style!="") {
        if(Mcmd.count(style)) {
          cout<<"Wrong commands!"<<endl;
          exit(0);
        }
        Mcmd.emplace(style,Vargs);
      }
      style=arg[ii];
      Vargs.clear();
      continue;
    }
    Vargs.push_back(arg[ii]);
  }
  if(Mcmd.count(style)) {
    cout<<"Wrong commands!"<<endl;
    exit(0);
  }
  Mcmd.emplace(style,Vargs);

  string seq="";
  if(Mcmd.count("-seq")) {
    seq=Mcmd["-seq"][0];
    for(size_t jj=1;jj<Mcmd["-seq"].size();jj++)
      seq+="&"+Mcmd["-seq"][jj];
    Mcmd.erase("-seq");
  }
  else {
    cout<<"Wrong commands: No input sequence (-seq)!"<<endl;
    exit(0);
  }
  
  string fout="helix.pdb";
  if(Mcmd.count("-o")) {
    fout=Mcmd["-o"][0];
    Mcmd.erase("-o");
  }
  else {
    cout<<"Wrong commands: No output (-o)!"<<endl;
    exit(0);
  }

  vector<string> Vres_info;
  if(Mcmd.count("-res")) {
    Vres_info=Mcmd["-res"];
    Mcmd.erase("-res");
  }

  string flog="./Cobmms.log";
  if(Mcmd.count("-log")) {
    vector<string> Vlog=Mcmd["-log"];
    if(Vlog.size()!=1) {
      string str="Unknown command option -->  -log";
      for(auto ss : Vlog)
        str+=" "+ss;
      cout<<str<<endl;
      exit(0);
    }
    flog=Vlog[0];
    Mcmd.erase("-log");
  } 
  if(flog.find(".log")==string::npos)
    flog+=".log";

  Cobmms *cbs= new Cobmms(flog);
  if(Mcmd.count("-homedir")==0) {
    cbs->parse_path_to_Cobmms(arg[0]);
  }
  else {
    cbs->set_path_to_Cobmms(Mcmd["-homedir"][0]);
    Mcmd.erase("-homedir");
  }
  cbs->initialization();
  cbs->parse_command(Mcmd);

  RNAMolecule rna_one(cbs);
  rna_one.set_RNA_molecule(seq);

  if(rna_one.chains.size()!=2 || rna_one.chains[0].residues.size()!=rna_one.chains[1].residues.size()) {
    string str="Wrong input sequence --> "+seq;
    cbs->error->one(str);
  }

  if(Vres_info.size()!=0 && Vres_info.size()!=2) {
    string str="Wrong residue information --> ";
    for(size_t jj=0;jj<Vres_info.size();jj++) {
      str+=" "+Vres_info[jj];
    }
    cbs->error->one(str);
  }

  if(Vres_info.size()==2) {
    for(size_t ii=0;ii<Vres_info.size();ii++) {
      vector<string> stoken=split_string(Vres_info[ii],":");
      string chainid_one=stoken[0];
      int resid_one=rna_one.chains[ii].residues[0].resid;
      if(stoken.size()>1)
        resid_one=stoi(stoken[1]);
      rna_one.chains[ii].reset_residues_name(chainid_one,resid_one);
    }
  } 

  string fstr=cbs->init->INPUT_DIR+"RNA_Aform_helix.pdb";
  ReadPDB readpdb_one(cbs);
  readpdb_one.read_PDB(fstr);
  RNAMolecule RNAHelix(cbs);
  RNAHelix.set_RNA_molecule(readpdb_one.RNA_atoms[0]);
  size_t max_helix_length=RNAHelix.chains[0].size();

  size_t index_begin=0;
  size_t length=rna_one.chains[0].size();
  size_t nres_overlapped=2;
  vector<RNAMolecule> Vhx;
  while(index_begin<length) {
    if((length-index_begin)<=max_helix_length) {
      Chain hx1(cbs),hx2(cbs);
      hx1.chain_label=Chain::rna;
      hx2.chain_label=Chain::rna;
      for(size_t ii=0;ii<length-index_begin;ii++) {
        hx1.insert_residue(rna_one.chains[0].residues[ii+index_begin],Chain::end);
        hx2.insert_residue(rna_one.chains[1].residues[length-ii-index_begin-1],Chain::begin);
      }
      Chain hx1_xyz(cbs),hx2_xyz(cbs);
      hx1_xyz.chain_label=Chain::rna;
      hx2_xyz.chain_label=Chain::rna;
      for(size_t ii=0;ii<length-index_begin;ii++) {
        hx1_xyz.insert_residue(RNAHelix.chains[0].residues[ii],Chain::end);
        hx2_xyz.insert_residue(RNAHelix.chains[1].residues[max_helix_length-ii-1],Chain::begin);
      }
      hx1.set_chain_xyz(hx1_xyz);
      hx2.set_chain_xyz(hx2_xyz);
      RNAMolecule hx_one(cbs);
      hx_one.insert_chain(hx1);
      hx_one.insert_chain(hx2);
      Vhx.push_back(hx_one);
      index_begin=length;
    }
    else {
      Chain hx1(cbs),hx2(cbs);
      hx1.chain_label=Chain::rna;
      hx2.chain_label=Chain::rna;
      for(size_t ii=0;ii<max_helix_length;ii++) {
        hx1.insert_residue(rna_one.chains[0].residues[ii+index_begin],Chain::end);
        hx2.insert_residue(rna_one.chains[1].residues[length-ii-index_begin-1],Chain::begin);
      }
      Chain hx1_xyz(cbs),hx2_xyz(cbs);
      hx1_xyz.chain_label=Chain::rna;
      hx2_xyz.chain_label=Chain::rna;
      for(size_t ii=0;ii<max_helix_length;ii++) {
        hx1_xyz.insert_residue(RNAHelix.chains[0].residues[ii],Chain::end);
        hx2_xyz.insert_residue(RNAHelix.chains[1].residues[max_helix_length-ii-1],Chain::begin);
      }
      hx1.set_chain_xyz(hx1_xyz);
      hx2.set_chain_xyz(hx2_xyz);
      RNAMolecule hx_one(cbs);
      hx_one.insert_chain(hx1);
      hx_one.insert_chain(hx2);
      Vhx.push_back(hx_one);
      index_begin+=max_helix_length-nres_overlapped;
    }
  }

  for(size_t ii=1;ii<Vhx.size();ii++) {
    vector<double> Vxyz_ref,Vxyz_tar;
    size_t nres_one=Vhx[ii-1].chains[0].residues.size();
    for(size_t jj=0;jj<nres_overlapped;jj++) {
      vector<double> xyz=Vhx[ii-1].chains[0].residues[nres_one-nres_overlapped+jj].get_residue_xyz();
      Vxyz_ref.insert(Vxyz_ref.end(),xyz.begin(),xyz.end());
      xyz.clear();
      xyz=Vhx[ii].chains[0].residues[jj].get_residue_xyz();
      Vxyz_tar.insert(Vxyz_tar.end(),xyz.begin(),xyz.end());
    }
    nres_one=Vhx[ii].chains[1].residues.size();
    for(size_t jj=0;jj<nres_overlapped;jj++) {
      vector<double> xyz=Vhx[ii-1].chains[1].residues[jj].get_residue_xyz();
      Vxyz_ref.insert(Vxyz_ref.end(),xyz.begin(),xyz.end());
      xyz.clear();
      xyz=Vhx[ii].chains[1].residues[nres_one-nres_overlapped+jj].get_residue_xyz();
      Vxyz_tar.insert(Vxyz_tar.end(),xyz.begin(),xyz.end());
    }
    Superpose sp_one(Vxyz_ref,Vxyz_tar);
    sp_one.calculate_rotation_matrix();
    vector<double> Vxyz_before=Vhx[ii].get_molecule_xyz();
    sp_one.reset_x2(Vxyz_before);
    vector<double> Vxyz_after=sp_one.rotate_x2();
    Vhx[ii].reset_molecule_xyz(Vxyz_after);
    for(size_t jj=nres_overlapped;jj<Vhx[ii].chains[0].residues.size();jj++) {
      Vhx[ii-1].chains[0].insert_residue(Vhx[ii].chains[0].residues[jj],Chain::end);
      Vhx[ii-1].chains[1].insert_residue(Vhx[ii].chains[1].residues[nres_one-jj-1],Chain::begin);
    } 
  }

  ofstream outP (fout.c_str());
  outP<<Vhx[0];
  outP.close(); outP.clear();  

  delete cbs;
  
  return 0;
}