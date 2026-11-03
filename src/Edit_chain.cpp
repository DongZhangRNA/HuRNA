#include "basic.h"
#include "error.h"
#include "readfile.h"
#include "initialization.h"
#include "molecule.h"
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

  vector<string> Vedit_args;
  if(Mcmd.count("-arg")) {
    Vedit_args=Mcmd["-arg"];
    Mcmd.erase("-arg");
  }
  else {
    cout<<"Wrong commands: No edit information (-arg) !"<<endl;
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

  struct ResInfo {
    string chainid;
    int resid_begin, resid_end;
  };
  vector<ResInfo> Volds,Vnews;
  for(size_t ii=0;ii<Vedit_args.size();ii++) {
    vector<string> stoken=split_string(Vedit_args[ii],">");
    if(stoken.size()!=2) {
      string str="Wrong edit args: "+Vedit_args[ii];
      cbs->error->one(str);
    }
    string sold=stoken[0];
    string snew=stoken[1];
    ResInfo old_i,new_i;
    size_t pos=sold.find_first_of(":");
    if(pos==string::npos && sold.size()>0) {
      old_i.chainid=sold;
      old_i.resid_begin=resid_null;
      old_i.resid_end=resid_null;
    }
    else if(pos==0) {
      old_i.chainid="";
      stoken.clear();
      stoken=split_string(sold,":");
      if(stoken.size()==1) {
        old_i.resid_begin=stoi(stoken[0]);
        old_i.resid_end=resid_null;
      }
      else if(stoken.size()==2) {
        old_i.resid_begin=stoi(stoken[0]);
        old_i.resid_end=stoi(stoken[1]);
      }
    }
    else {
      stoken.clear();
      stoken=split_string(sold,":");
      if(stoken.size()==1) {
        old_i.chainid=stoken[0];
        old_i.resid_begin=resid_null;
        old_i.resid_end=resid_null;
      }
      else if(stoken.size()==2) {
        old_i.chainid=stoken[0];
        old_i.resid_begin=stoi(stoken[1]);
        old_i.resid_end=resid_null;
      }
      else if(stoken.size()==3) {
        old_i.chainid=stoken[0];
        old_i.resid_begin=stoi(stoken[1]);
        old_i.resid_end=stoi(stoken[2]);
      }
    }
    Volds.push_back(old_i);

    pos=snew.find_first_of(":");
    if(pos==string::npos && snew.size()>0) {
      new_i.chainid=snew;
      new_i.resid_begin=resid_null;
      new_i.resid_end=resid_null;
    }
    else if(pos==0) {
      new_i.chainid="";
      stoken.clear();
      stoken=split_string(snew,":");
      if(stoken.size()==1) {
        new_i.resid_begin=stoi(stoken[0]);
        new_i.resid_end=resid_null;
      }
      else if(stoken.size()==2) {
        new_i.resid_begin=stoi(stoken[0]);
        new_i.resid_end=stoi(stoken[1]);
      }
    }
    else {
      stoken.clear();
      stoken=split_string(snew,":");
      if(stoken.size()==1) {
        new_i.chainid=stoken[0];
        new_i.resid_begin=resid_null;
        new_i.resid_end=resid_null;
      }
      else if(stoken.size()==2) {
        new_i.chainid=stoken[0];
        new_i.resid_begin=stoi(stoken[1]);
        new_i.resid_end=resid_null;
      }
      else if(stoken.size()==3) {
        new_i.chainid=stoken[0];
        new_i.resid_begin=stoi(stoken[1]);
        new_i.resid_end=stoi(stoken[2]);
      }
    }
    Vnews.push_back(new_i);
  }
  vector<Molecule> Vmol_new;
  for(int ii=0;ii<nmodel;ii++) {
    Molecule mol_one(cbs);
    mol_one.set_molecule(readpdb_one.Full_atoms[ii],Chain::other);
    Molecule mol_new(cbs);
    for(size_t jj=0;jj<mol_one.chains.size();jj++) {
      string chainid_one=mol_one.chains[jj].chainid;
      Chain chain_modified=mol_one.chains[jj];
      int flag_check=0;
      vector<size_t> Vindex_modified;
      for(size_t kk=0;kk<Volds.size();kk++) {
        if(Volds[kk].chainid==chainid_one || (Volds[kk].chainid.size()==0 && chainid_one.size()==0)) {
          flag_check=1;
          if(Volds[kk].resid_begin==resid_null) 
            Volds[kk].resid_begin=mol_one.chains[jj].residues[0].resid;
          if(Volds[kk].resid_end==resid_null) 
            Volds[kk].resid_end=mol_one.chains[jj].residues.back().resid;
          if(Vnews[kk].resid_begin==resid_null)
            Vnews[kk].resid_begin=Volds[kk].resid_begin;
        
          int resid_new=Vnews[kk].resid_begin;
          for(size_t nn=0;nn<mol_one.chains[jj].residues.size();nn++) {
            Residue &residue_one=mol_one.chains[jj].residues[nn];
            int resid_one=residue_one.resid;
            if(resid_one>=Volds[kk].resid_begin && resid_one<=Volds[kk].resid_end) {
              chain_modified.residues[nn].reset_residue_name(Vnews[kk].chainid,resid_new,residue_one.resname);
              resid_new++;
              Vindex_modified.push_back(nn);
            }
          }
        }
      }
      Chain chain_new(cbs);
      chain_new.chain_label=chain_modified.chain_label;
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
      }
      else {
        if(flag_only==0)
          chain_new=mol_one.chains[jj];  
      }
      if(chain_new.residues.size()>0) mol_new.insert_chain(chain_new);
    }
    Vmol_new.push_back(mol_new);
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