#include "basic.h"
#include "error.h"
#include "readfile.h"
#include "initialization.h"
#include "RNA_molecule.h"
#include "read_PDB.h"

using namespace CBS_NS;

const int resid_null = -1000000;

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
  string inFile,outFile;
  string flog="./Cobmms.log";
  int label=Chain::rna;

  if(Mcmd.count("-i")==0) {
    cout<<"Wrong commands: No input (-i) !"<<endl;
    exit(0);
  }
  inFile=Mcmd["-i"][0];
  Mcmd.erase("-i");

  if(Mcmd.count("-o")==0) {
    cout<<"Wrong commands: No output (-o) !"<<endl;
    exit(0);
  }
  outFile=Mcmd["-o"][0];
  Mcmd.erase("-o");

  string arg_one="";
  if(Mcmd.count("-sel")) {
    arg_one=Mcmd["-sel"][0];
    for(size_t ii=1;ii<Mcmd["-sel"].size();ii++)
      arg_one+=" "+Mcmd["-sel"][ii];
    Mcmd.erase("-sel");
  }
  else {
    cout<<"Wrong commands: No chain information (-sel) !"<<endl;
    exit(0);
  }

  string seq_one="";
  if(Mcmd.count("-seq")) {
    seq_one=Mcmd["-seq"][0];
    for(size_t ii=1;ii<Mcmd["-seq"].size();ii++) 
      seq_one+="&"+Mcmd["-seq"][ii];
    Mcmd.erase("-seq");
  }
  else {
    cout<<"Wrong commands: No sequence information (-seq) !"<<endl;
    exit(0);
  }

  int nmodel=-1;
  if(Mcmd.count("-nmodel")) {
    nmodel=stoi(Mcmd["-nmodel"][0]);
    Mcmd.erase("-nmodel");
  }

  int flag_noH=0;
  if(Mcmd.count("-noH")) {
    flag_noH=1;
    Mcmd.erase("-noH");
  }

  int flag_osep=0;
  if(Mcmd.count("-os")) {
    flag_osep=1;
    Mcmd.erase("-os");
  }

  int flag_only=0;
  if(Mcmd.count("-only")) {
    flag_only=1;
    Mcmd.erase("-only");
  }

  if(Mcmd.count("-label")) {
    string slabel=Mcmd["-label"][0];
    if(slabel=="rna")
      label=Chain::rna;
    else if(slabel=="dna")
      label=Chain::dna;
    else if(slabel=="protein")
      label=Chain::protein;
    else {
      cout<<"Wrong commands: Unknown label (-label) --> "<<slabel<<endl;
      exit(0);
    }
  }
  Mcmd.erase("-label");

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

  cbs->init->noH=flag_noH;

  ReadPDB readpdb_one(cbs);
  readpdb_one.read_PDB(inFile);
  int nmax=static_cast<int>(readpdb_one.Full_atoms.size());
  if(nmodel<0 || nmodel>nmax)
    nmodel=nmax;

  vector<string> Vchain_args=split_string(arg_one,",");

  struct ResInfo {
    string chainid;
    int resid_begin, resid_end;
  };
  vector<ResInfo> Vsegments;
  for(size_t ii=0;ii<Vchain_args.size();ii++) {
    string sseg=Vchain_args[ii];
    ResInfo seg_i;
    size_t pos=sseg.find_first_of(":");
    if(pos==string::npos && sseg.size()>0) {
      seg_i.chainid=sseg;
      seg_i.resid_begin=resid_null;
      seg_i.resid_end=resid_null;
    }
    else if(pos==0) {
      seg_i.chainid="";
      vector<string> stoken=split_string(sseg,":");
      if(stoken.size()==1) {
        seg_i.resid_begin=stoi(stoken[0]);
        seg_i.resid_end=resid_null;
      }
      else if(stoken.size()==2) {
        seg_i.resid_begin=stoi(stoken[0]);
        seg_i.resid_end=stoi(stoken[1]);
      }
    }
    else {
      vector<string> stoken=split_string(sseg,":");
      if(stoken.size()==1) {
        seg_i.chainid=stoken[0];
        seg_i.resid_begin=resid_null;
        seg_i.resid_end=resid_null;
      }
      else if(stoken.size()==2) {
        seg_i.chainid=stoken[0];
        seg_i.resid_begin=stoi(stoken[1]);
        seg_i.resid_end=resid_null;
      }
      else if(stoken.size()==3) {
        seg_i.chainid=stoken[0];
        seg_i.resid_begin=stoi(stoken[1]);
        seg_i.resid_end=stoi(stoken[2]);
      }
    }
    Vsegments.push_back(seg_i);
  }

  RNAMolecule rna_one(cbs);
  vector<vector<string> > Vresnames_all;
  vector<string> Vseq;
  if(label==Chain::rna) {
    Vseq=rna_one.read_sequence(seq_one);
    for(size_t ii=0;ii<Vseq.size();ii++) {
      vector<string> Vresnames=rna_one.get_resnames_from_sequence(Vseq[ii]);
      Vresnames_all.push_back(Vresnames);
    }
  }
  if(Vsegments.size()!=Vresnames_all.size()) {
    string str="Number of segments (-arg) and sequence (-seq) is not matched  --> ";
    str+=to_string(Vsegments.size())+" / "+to_string(Vresnames_all.size());
    cbs->error->one(str);
  }

  vector<Molecule> Vmol_new;
  for(int ii=0;ii<nmodel;ii++) {
    Molecule mol_one(cbs);
    mol_one.set_molecule(readpdb_one.Full_atoms[ii],Chain::other);
    Molecule mol_new(cbs);
    vector<int> Vsegments_check;
    for(size_t kk=0;kk<Vsegments.size();kk++)
      Vsegments_check.push_back(0);
    for(size_t jj=0;jj<mol_one.chains.size();jj++) {
      string chainid_one=mol_one.chains[jj].chainid;
      Chain chain_modified=mol_one.chains[jj];
      int flag_check=0;
      vector<size_t> Vindex_modified;
      for(size_t kk=0;kk<Vsegments.size();kk++) {
        if(Vsegments[kk].chainid==chainid_one || (Vsegments[kk].chainid.size()==0 && chainid_one.size()==0)) {
          flag_check=1;
          Vsegments_check[kk]=1;
          if(Vsegments[kk].resid_begin==resid_null) 
            Vsegments[kk].resid_begin=mol_one.chains[jj].residues[0].resid;
          if(Vsegments[kk].resid_end==resid_null) 
            Vsegments[kk].resid_end=mol_one.chains[jj].residues.back().resid;
          
          int nres=Vsegments[kk].resid_end-Vsegments[kk].resid_begin+1;
          if(nres!=(int)Vresnames_all[kk].size()) {
            string str="Number of residues is not matched for segment "+Vsegments[kk].chainid;
            str+=":"+to_string(Vsegments[kk].resid_begin)+":"+to_string(Vsegments[kk].resid_end);
            str+=" --> "+to_string(nres)+" / "+to_string(Vresnames_all[kk].size())+" (";
            str+=Vseq[kk]+")";
            cbs->error->one(str);
          }
          size_t index_i=0;
          for(size_t nn=0;nn<mol_one.chains[jj].residues.size();nn++) {
            Residue &residue_one=mol_one.chains[jj].residues[nn];
            int resid_one=residue_one.resid;
            if(resid_one>=Vsegments[kk].resid_begin && resid_one<=Vsegments[kk].resid_end) {
              if(label==Chain::rna)
                chain_modified[nn].simplify_RNA_base();
              chain_modified[nn].reset_residue_name(chainid_one,resid_one,Vresnames_all[kk][index_i]);
              if(label==Chain::rna)
                chain_modified[nn].recovery_RNA_base();
              index_i++;
              Vindex_modified.push_back(nn);
            }
          }
        }
      }
      Chain chain_new(cbs);
      if(flag_check) {
        if(flag_only==0) {
          chain_new=chain_modified;
        }
        else {
          for(size_t nn=0;nn<chain_modified.residues.size();nn++) {
            if(count(Vindex_modified.begin(),Vindex_modified.end(),nn))
              chain_new.insert_residue(chain_modified.residues[nn],Chain::end);
          }
        }
        chain_new.chain_label=label;
      }
      else {
        if(flag_only==0)
          chain_new=mol_one.chains[jj];  
      }
      if(chain_new.residues.size()>0) mol_new.insert_chain(chain_new);
    }
    Vmol_new.push_back(mol_new);

    for(size_t kk=0;kk<Vsegments_check.size();kk++) {
      if(Vsegments_check[kk]==0) {
        string str="Cannot find segment "+Vsegments[kk].chainid+":"+to_string(Vsegments[kk].resid_begin)+":";
        str+=to_string(Vsegments[kk].resid_end)+" in "+inFile;
        cbs->error->warning(str);
      }
    }
  }

  size_t pos_fmt=max(0,(int)outFile.size()-4);
  string fout_fmt=outFile.substr(pos_fmt);
  if(fout_fmt!=".pdb") 
    outFile+=".pdb";
  string fout_sep=outFile.substr(0,outFile.size()-4)+"_";
  ofstream outP (outFile.c_str());
  for(size_t ii=0;ii<Vmol_new.size();ii++) {
    if(Vmol_new.size()>1)
      outP<<"MODEL     "<<ii+1<<endl;
    outP<<Vmol_new[ii];
    if(Vmol_new.size()>1)
      outP<<"ENDMDL"<<endl;

    if(flag_osep) {
      string fstr=fout_sep+to_string(ii+1)+".pdb";
      ofstream outF (fstr.c_str());
      outF<<Vmol_new[ii];
      outF.close(); outF.clear();
    }
  }
  if(Vmol_new.size()>1)
    outP<<"END"<<endl;
  outP.close(); outP.clear();

  delete cbs;
  
  return 0;
}