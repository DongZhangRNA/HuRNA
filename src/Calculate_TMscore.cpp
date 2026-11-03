#include "basic.h"
#include "error.h"
#include "readfile.h"
#include "initialization.h"
#include "read_PDB.h"
#include "molecule.h"

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

  string fin;
  if(Mcmd.count("-i")) {
    fin=Mcmd["-i"][0];
    Mcmd.erase("-i");
  }
  else {
    cout<<"Wrong commands: No input (-i)!"<<endl;
    exit(0);
  }

  string tar_res="all";
  if(Mcmd.count("-sel_res")) {
    tar_res=Mcmd["-sel_res"][0];
    for(size_t jj=1;jj<Mcmd["-sel_res"].size();jj++)
      tar_res+=" "+Mcmd["-sel_res"][jj];
    Mcmd.erase("-sel_res");
  }

  string fref="";
  if(Mcmd.count("-ref")) {
    fref=Mcmd["-ref"][0];
    Mcmd.erase("-ref");
  }
  else {
    cout<<"Wrong commands: No reference structure (-ref)!"<<endl;
    exit(0);
  }

  string ref_res="all";
  if(Mcmd.count("-ref_res")) {
    ref_res=Mcmd["-ref_res"][0];
    for(size_t jj=1;jj<Mcmd["-ref_res"].size();jj++)
      ref_res+=" "+Mcmd["-ref_res"][jj];
    Mcmd.erase("-ref_res");
  }
  
  string fout="";
  if(Mcmd.count("-o")) {
    fout=Mcmd["-o"][0];
    Mcmd.erase("-o");
  }
  else {
    cout<<"Wrong commands: No output (-o)!"<<endl;
    exit(0);
  }

  int nmodel=-1;
  if(Mcmd.count("-nmodel")) {
    nmodel=stoi(Mcmd["-nmodel"][0]);
    Mcmd.erase("-nmodel");
  }

  string USalign="";
  int flag_USalign=0;
  if(Mcmd.count("-usalign")) {
    flag_USalign=1;
    USalign=Mcmd["-usalign"][0];
    Mcmd.erase("-usalign");
  }
  int flag_rm=1;
  if(Mcmd.count("-rm")) {
    flag_rm=stoi(Mcmd["-rm"][0]);
    Mcmd.erase("-rm");
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

  cbs->init->noH=1;  //no H atoms

  if(flag_USalign==0) USalign=cbs->init->USalign;

  if(!is_file_exist(USalign)) {
    string str="Cannot find USalign! Please compile USalign in PATH_to_Cobmms/external/TMscore first!";
    cbs->error->one(str);
  }

  ReadPDB readpdb_ref(cbs);
  readpdb_ref.read_PDB(fref);
  if(readpdb_ref.Full_atoms.size()==0 || readpdb_ref.Full_atoms[0].size()==0) {
    string str="No atom in reference structure "+fref;
    cbs->error->one(str);
  }
  Molecule mol_ref(cbs);
  mol_ref.set_molecule(readpdb_ref.Full_atoms[0],Chain::empty);

  ReadPDB readpdb_tar(cbs);
  readpdb_tar.read_PDB(fin);
  if(readpdb_tar.Full_atoms.size()==0 || readpdb_tar.Full_atoms[0].size()==0) {
    string str="No atom in input structure "+fin;
    cbs->error->one(str);
  }

  int nmax=static_cast<int>(readpdb_tar.Full_atoms.size());
  if(nmodel<0 || nmodel>nmax)
    nmodel=nmax;

  struct ResInfo {
    string chainid;
    int resid_begin, resid_end;
  };
  vector<Residue> Vresidues_ref;
  if(ref_res=="all") {
    for(size_t ii=0;ii<mol_ref.chains.size();ii++) {
      for(size_t jj=0;jj<mol_ref.chains[ii].residues.size();jj++) {
        Vresidues_ref.push_back(mol_ref.chains[ii].residues[jj]);
      }
    }
  }
  else {
    vector<ResInfo> Vrefs;
    vector<string> Varg=split_string(ref_res,",");
    for(size_t ii=0;ii<Varg.size();ii++) {
      string &seg_i=Varg[ii];
      ResInfo resinfo_i;
      size_t pos=seg_i.find_first_of(":");
      if(pos==string::npos && seg_i.size()>0) {
        resinfo_i.chainid=seg_i;
        resinfo_i.resid_begin=resid_null;
        resinfo_i.resid_end=resid_null;
      }
      else if(pos==0) {
        resinfo_i.chainid="";
        vector<string> stoken=split_string(seg_i,":");
        if(stoken.size()==1) {
          resinfo_i.resid_begin=stoi(stoken[0]);
          resinfo_i.resid_end=resid_null;
        }
        else if(stoken.size()==2) {
          resinfo_i.resid_begin=stoi(stoken[0]);
          resinfo_i.resid_end=stoi(stoken[1]);
        }
      }
      else {
        vector<string> stoken=split_string(seg_i,":");
        if(stoken.size()==1) {
          resinfo_i.chainid=stoken[0];
          resinfo_i.resid_begin=resid_null;
          resinfo_i.resid_end=resid_null;
        }
        else if(stoken.size()==2) {
          resinfo_i.chainid=stoken[0];
          resinfo_i.resid_begin=stoi(stoken[1]);
          resinfo_i.resid_end=resid_null;
        }
        else if(stoken.size()==3) {
          resinfo_i.chainid=stoken[0];
          resinfo_i.resid_begin=stoi(stoken[1]);
          resinfo_i.resid_end=stoi(stoken[2]);
        }
      }
      Vrefs.push_back(resinfo_i);
    }

    for(size_t ii=0;ii<Vrefs.size();ii++) {
      for(size_t jj=0;jj<mol_ref.chains.size();jj++) {
        string chainid_one=mol_ref.chains[jj].chainid;
        if(Vrefs[ii].chainid==chainid_one || (Vrefs[ii].chainid.length()==0 && chainid_one.length()==0)) {
          if(Vrefs[ii].resid_begin==resid_null) 
            Vrefs[ii].resid_begin=mol_ref.chains[jj].residues[0].resid;
          if(Vrefs[ii].resid_end==resid_null) 
            Vrefs[ii].resid_end=mol_ref.chains[jj].residues.back().resid;
          for(size_t kk=0;kk<mol_ref.chains[jj].residues.size();kk++) {
            Residue &residue_one=mol_ref.chains[jj].residues[kk];
            int resid_one=residue_one.resid;
            if(resid_one>=Vrefs[ii].resid_begin && resid_one<=Vrefs[ii].resid_end)
              Vresidues_ref.push_back(residue_one);
          }
          break;
        }
      }
    }
  }
  
  vector<ResInfo> Vtars;
  if(tar_res!="all") {
    vector<string> Varg=split_string(tar_res,",");
    for(size_t ii=0;ii<Varg.size();ii++) {
      string &seg_i=Varg[ii];
      ResInfo resinfo_i;
      size_t pos=seg_i.find_first_of(":");
      if(pos==string::npos && seg_i.size()>0) {
        resinfo_i.chainid=seg_i;
        resinfo_i.resid_begin=resid_null;
        resinfo_i.resid_end=resid_null;
      }
      else if(pos==0) {
        resinfo_i.chainid="";
        vector<string> stoken=split_string(seg_i,":");
        if(stoken.size()==1) {
          resinfo_i.resid_begin=stoi(stoken[0]);
          resinfo_i.resid_end=resid_null;
        }
        else if(stoken.size()==2) {
          resinfo_i.resid_begin=stoi(stoken[0]);
          resinfo_i.resid_end=stoi(stoken[1]);
        }
      }
      else {
        vector<string> stoken=split_string(seg_i,":");
        if(stoken.size()==1) {
          resinfo_i.chainid=stoken[0];
          resinfo_i.resid_begin=resid_null;
          resinfo_i.resid_end=resid_null;
        }
        else if(stoken.size()==2) {
          resinfo_i.chainid=stoken[0];
          resinfo_i.resid_begin=stoi(stoken[1]);
          resinfo_i.resid_end=resid_null;
        }
        else if(stoken.size()==3) {
          resinfo_i.chainid=stoken[0];
          resinfo_i.resid_begin=stoi(stoken[1]);
          resinfo_i.resid_end=stoi(stoken[2]);
        }
      }
      Vtars.push_back(resinfo_i);
    }
  }

  Molecule mol_ref_sel(cbs);
  Chain chain_one(cbs);
  string chainid0=Vresidues_ref[0].chainid;
  for(size_t ii=0;ii<Vresidues_ref.size();ii++) {
    string chainid_one=Vresidues_ref[ii].chainid;
    if(chainid_one!=chainid0) {
      chain_one.chain_label=Chain::empty;
      if(chain_one.size()>0) 
        mol_ref_sel.insert_chain(chain_one);
      chain_one.clear();
    }
    chain_one.insert_residue(Vresidues_ref[ii],Chain::end);
    chainid0=chainid_one;
  }
  if(chain_one.size()>0) {
    chain_one.chain_label=Chain::empty;
    mol_ref_sel.insert_chain(chain_one);
    chain_one.clear();
  }
  mol_ref_sel.standardize_molecule();
  string fref_sel="./TMscore_reference.pdb";
  ofstream outP (fref_sel.c_str());
  outP<<mol_ref_sel;
  outP.close(); outP.clear();

  string ftar_sel="./TMscore_target.pdb";
  string fout_TMscore="./TMscore_out.dat";
  ofstream outF (fout.c_str());
  outF<<"#Model TMscore"<<endl;
  for(int ii=0;ii<nmodel;ii++) {
    Molecule mol_tar(cbs);
    mol_tar.set_molecule(readpdb_tar.Full_atoms[ii],Chain::empty);
    vector<Residue> Vresidues_tar;
    if(tar_res=="all") {
      for(size_t jj=0;jj<mol_tar.chains.size();jj++) {
        for(size_t kk=0;kk<mol_tar.chains[jj].residues.size();kk++) {
          Vresidues_tar.push_back(mol_tar.chains[jj].residues[kk]);
        }
      }
    }
    else {
      for(size_t jj=0;jj<Vtars.size();jj++) {
        for(size_t kk=0;kk<mol_tar.chains.size();kk++) {
          string chainid_one=mol_tar.chains[kk].chainid;
          if(Vtars[jj].chainid==chainid_one || (Vtars[jj].chainid.length()==0 && chainid_one.length()==0)) {
            if(Vtars[jj].resid_begin==resid_null) 
              Vtars[jj].resid_begin=mol_tar.chains[kk].residues[0].resid;
            if(Vtars[jj].resid_end==resid_null) 
              Vtars[jj].resid_end=mol_tar.chains[kk].residues.back().resid;
            for(size_t nn=0;nn<mol_tar.chains[kk].residues.size();nn++) {
              Residue &residue_one=mol_tar.chains[kk].residues[nn];
              int resid_one=residue_one.resid;
              if(resid_one>=Vtars[jj].resid_begin && resid_one<=Vtars[jj].resid_end)
                Vresidues_tar.push_back(residue_one);
            }
            break;
          }
        }
      }
    }

    if(Vresidues_ref.size()!=Vresidues_tar.size()) {
      string str="Number of residues is not matched --> "+to_string(Vresidues_ref.size())+"(ref) / "+to_string(Vresidues_tar.size())+"(tar)";
      cbs->error->one(str);
    }

    Molecule mol_tar_sel(cbs);
    chain_one.clear();
    chainid0=Vresidues_tar[0].chainid;
    for(size_t jj=0;jj<Vresidues_tar.size();jj++) {
      string chainid_one=Vresidues_tar[jj].chainid;
      if(chainid_one!=chainid0) {
        chain_one.chain_label=Chain::empty;
        if(chain_one.size()>0) 
          mol_tar_sel.insert_chain(chain_one);
        chain_one.clear();
      }
      chain_one.insert_residue(Vresidues_tar[jj], Chain::end);
      chainid0=chainid_one;
    }
    if(chain_one.size()>0) {
      chain_one.chain_label=Chain::empty;
      mol_tar_sel.insert_chain(chain_one);
      chain_one.clear();
    }
    mol_tar_sel.standardize_molecule();
    outP.open(ftar_sel.c_str());
    outP<<mol_tar_sel;
    outP.close(); outP.clear();

    string cmd=USalign+" "+ftar_sel+" "+fref_sel;
    if(mol_ref_sel.size()>1)
      cmd+=" -mm 1 -ter 0";
    cmd+=" > "+fout_TMscore;

    cbs->execute_system_command(cmd);

    vector<string> Vlines=read_file(fout_TMscore);
    double TMscore;
    for(size_t kk=0;kk<Vlines.size();kk++) {
      string sline=Vlines[kk];
      if(sline.size()<10) continue;
      string shead=sline.substr(0,9);
      if(shead.find("TM-score=")!=string::npos) {
        sline=sline.substr(9);
        vector<string> stoken=split_string(sline);
        TMscore=stod(stoken[0]);
        break;
      }
    }
    if(nmodel>1) {
      outF<<fin<<"_Model#"<<ii+1<<" "<<TMscore<<endl;
    }
    else {
      outF<<fin<<" "<<TMscore<<endl;
    }
  }
  outF.close(); outF.clear();

  if(flag_rm) {
    string cmd="rm "+fref_sel+" "+ftar_sel+" "+fout_TMscore;
    cbs->execute_system_command(cmd);
  }

  delete cbs;
  
  return 0;
}