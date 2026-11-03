#include "basic.h"
#include "error.h"
#include "readfile.h"
#include "initialization.h"
#include "RNA_molecule.h"
#include "read_PDB.h"
#include "superpose.h"

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

  vector<string> pdblist;
  int nmodel=1;  //reserve #nmodel models in input PDB
  int flag_model=1; // flag_model=0, output moldes one by one; flag_model=1, output all models in single file
  int flag_out=1;  //flag_out=0, output to file "outname", changed filename; flag_out=1, output to directory "outDIR", unchanged filename 
  int flag_rna=1, flag_dna=0, flag_protein=0;
  string inDIR="";
  string outDIR="";
  string outname="";
  string flog="./Cobmms.log";

  if(Mcmd.count("-i")) {
    for(auto str : Mcmd["-i"])
      pdblist.push_back(str);
    Mcmd.erase("-i");
  }
  if(Mcmd.count("-inDIR")) {
    inDIR=Mcmd["-inDIR"][0];
    Mcmd.erase("-inDIR");
  }
  if(inDIR.size()!=0 && inDIR.back()!='/')
    inDIR+="/";
  if(Mcmd.count("-list")) {
    string flist=Mcmd["-list"][0];
    vector<string> Vlines=read_file(flist);
    for(auto sline : Vlines) {
      vector<string> stoken=split_string(sline);
      for(size_t ii=0;ii<stoken.size();ii++) {
        if(stoken[ii].find(".pdb")==string::npos && stoken[ii].find(".cif")==string::npos) continue;
        pdblist.push_back(stoken[ii]);
      }
    }
    Mcmd.erase("-list");
  }
  if(Mcmd.count("-outDIR")) {
    outDIR=Mcmd["-outDIR"][0];
    Mcmd.erase("-outDIR");
  }
  if(outDIR.size()!=0 && outDIR.back()!='/')
    outDIR+="/";
  if(Mcmd.count("-o")) {
    flag_out=0;
    string name=Mcmd["-o"][0];
    size_t found=name.find(".pdb");
    if(found!=string::npos)
      name=name.substr(0,found);
    outname=name;
    Mcmd.erase("-o");
  }
  if(Mcmd.count("-nmodel")) {
    nmodel=stoi(Mcmd["-nmodel"][0]);
    Mcmd.erase("-nmodel");
  }
  if(Mcmd.count("-merge")) {
    flag_model=stoi(Mcmd["-merge"][0]);
    Mcmd.erase("-merge");
    if(flag_model!=0 && flag_model!=1) {
      string str="Wrong value for -merge: 0/1 only!";
      cout<<str<<endl;
      exit(0);
    }
  }
  int flag_break_check=1;
  if(Mcmd.count("-chain_break_check")) {
    flag_break_check=stoi(Mcmd["-chain_break_check"][0]);
    Mcmd.erase("-chain_break_check");
    if(flag_break_check!=0 && flag_break_check!=1) {
      string str="Wrong value for -chain_break_check: 0/1 only!";
      cout<<str<<endl;
      exit(0);
    }
  }
  int flag_stand=1;
  if(Mcmd.count("-stand")) {
    flag_stand=stoi(Mcmd["-stand"][0]);
    Mcmd.erase("-stand");
    if(flag_stand!=0 && flag_stand!=1) {
      string str="Wrong value for -stand: 0/1 only!";
      cout<<str<<endl;
      exit(0);
    }
  }
  vector<string> Vtype;
  if(Mcmd.count("-type")) {
    Vtype=Mcmd["-type"];
    Mcmd.erase("-type");
  }
  if(Vtype.size()>0) {
    flag_rna=0;
    for(size_t ii=0;ii<Vtype.size();ii++) {
      if(Vtype[ii]=="rna") flag_rna=1;
      else if(Vtype[ii]=="dna") flag_dna=1;
      else if(Vtype[ii]=="protein") flag_protein=1;
      else if(Vtype[ii]=="all") {
        flag_rna=1;
        flag_dna=1;
        flag_protein=1;
      }
      else {
        string str="Unknown type: rna/dna/protein/all only!";
        cout<<str<<endl;
        exit(0);
      }
    }
  }

  vector<string> Vsel_chains;
  if(Mcmd.count("-sel")) {
    Vsel_chains=Mcmd["-sel"];
    Mcmd.erase("-sel");
  }
  if(Vsel_chains.size()==0)
    Vsel_chains.push_back("all");
  
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

  if(outDIR.size()>0) {
    if(!is_directory_exist(outDIR)) {
      string cmd="mkdir "+outDIR;
      cbs->execute_system_command(cmd);
    }

    cbs->init->WORK_DIR=outDIR;
  }

  cbs->init->RNA_check_chain_break=flag_break_check;

  if(pdblist.size()==0) {
    string str="No input!";
    cbs->error->one(str);
  }

  if(flag_out==0) {
    vector<Molecule> Vcomplex;
    for(size_t ii=0;ii<pdblist.size();ii++) {
      string fstr=inDIR+pdblist[ii];
      ReadPDB readpdb_one(cbs,fstr);
      int nm=(int)readpdb_one.Full_atoms.size();
      if(nm>=nmodel) nm=nmodel;

      if(Vsel_chains[0]!="all") {
        for(int jj=0;jj<nm;jj++) {
          vector<Atom> Vfull_tmp;
          for(size_t kk=0;kk<readpdb_one.Full_atoms[jj].size();kk++) {
            Atom& atom_one=readpdb_one.Full_atoms[jj][kk];
            if(count(Vsel_chains.begin(),Vsel_chains.end(),atom_one.chainid))
              Vfull_tmp.push_back(atom_one);
          }
          readpdb_one.Full_atoms[jj]=Vfull_tmp;

          vector<Atom> Vrna_tmp;
          for(size_t kk=0;kk<readpdb_one.RNA_atoms[jj].size();kk++) {
            Atom& atom_one=readpdb_one.RNA_atoms[jj][kk];
            if(count(Vsel_chains.begin(),Vsel_chains.end(),atom_one.chainid))
              Vrna_tmp.push_back(atom_one);
          }
          readpdb_one.RNA_atoms[jj]=Vrna_tmp;

          vector<Atom> Vdna_tmp;
          for(size_t kk=0;kk<readpdb_one.DNA_atoms[jj].size();kk++) {
            Atom& atom_one=readpdb_one.DNA_atoms[jj][kk];
            if(count(Vsel_chains.begin(),Vsel_chains.end(),atom_one.chainid))
              Vdna_tmp.push_back(atom_one);
          }
          readpdb_one.DNA_atoms[jj]=Vdna_tmp;

          vector<Atom> Vprotein_tmp;
          for(size_t kk=0;kk<readpdb_one.Protein_atoms[jj].size();kk++) {
            Atom& atom_one=readpdb_one.Protein_atoms[jj][kk];
            if(count(Vsel_chains.begin(),Vsel_chains.end(),atom_one.chainid))
              Vprotein_tmp.push_back(atom_one);
          }
          readpdb_one.Protein_atoms[jj]=Vprotein_tmp;
        }
      }

      if(nm==0 || readpdb_one.Full_atoms[0].size()==0) {
        string str="No atoms in file "+pdblist[ii]+", skipped!";
        cbs->error->message(str);
        continue;
      }
      for(int jj=0;jj<nm;jj++) {
        Molecule mol_one(cbs);
        if(flag_rna) {
          if(readpdb_one.RNA_atoms[jj].size()==0) {
            string str="No RNA atoms in file "+pdblist[ii]+", skipped!";
            cbs->error->message(str);
          }
          else {
            RNAMolecule rna_one(cbs);
            rna_one.set_RNA_molecule(readpdb_one.RNA_atoms[jj]);
            for(size_t kk=0;kk<rna_one.chains.size();kk++)
              mol_one.insert_chain(rna_one.chains[kk]);
          }
        }
        if(flag_protein) {
          //set_protein_molecule
        }
        if(flag_dna) {
          //set_dna_molecule
        }

        if(mol_one.chains.size()==0 || mol_one.chains[0].residues.size()==0) continue;

        if(flag_stand) mol_one.standardize_molecule();
        Vcomplex.push_back(mol_one);
      }
    }

    if(flag_model==0) {
      for(size_t jj=0;jj<Vcomplex.size();jj++) {
        string fstr=outDIR+outname+"-"+to_string(jj+1)+".pdb";
        ofstream outP (fstr.c_str());
        outP<<Vcomplex[jj];
        outP.close(); outP.clear();
      }
    }
    else if(flag_model==1) {
      string fstr=outDIR+outname+".pdb";
      ofstream outP (fstr.c_str());
      for(size_t jj=0;jj<Vcomplex.size();jj++) {
        if(Vcomplex.size()>1)
          outP<<"MODEL     "<<jj+1<<endl;
        outP<<Vcomplex[jj];
        if(Vcomplex.size()>1)
          outP<<"ENDMDL"<<endl;
      }
      outP<<"END"<<endl;
      outP.close(); outP.clear();
    }
  }
  else if(flag_out==1) {
    if(outDIR.size()==0) {
      string str="Please specify the output directory by -outDIR!";
      cbs->error->one(str);
    }
    for(size_t ii=0;ii<pdblist.size();ii++) {
      if(ii>0 && ii%100==0) {
        string str=to_string(ii)+"/"+to_string(pdblist.size())+" ...";
        cbs->error->message(str);
      }
      string fstr=inDIR+pdblist[ii];
      ReadPDB readpdb_one(cbs,fstr);
      int nm=(int)readpdb_one.Full_atoms.size();
      if(nm>=nmodel) nm=nmodel;

      if(Vsel_chains[0]!="all") {
        for(int jj=0;jj<nm;jj++) {
          vector<Atom> Vfull_tmp;
          for(size_t kk=0;kk<readpdb_one.Full_atoms[jj].size();kk++) {
            Atom& atom_one=readpdb_one.Full_atoms[jj][kk];
            if(count(Vsel_chains.begin(),Vsel_chains.end(),atom_one.chainid))
              Vfull_tmp.push_back(atom_one);
          }
          readpdb_one.Full_atoms[jj]=Vfull_tmp;

          vector<Atom> Vrna_tmp;
          for(size_t kk=0;kk<readpdb_one.RNA_atoms[jj].size();kk++) {
            Atom& atom_one=readpdb_one.RNA_atoms[jj][kk];
            if(count(Vsel_chains.begin(),Vsel_chains.end(),atom_one.chainid))
              Vrna_tmp.push_back(atom_one);
          }
          readpdb_one.RNA_atoms[jj]=Vrna_tmp;

          vector<Atom> Vdna_tmp;
          for(size_t kk=0;kk<readpdb_one.DNA_atoms[jj].size();kk++) {
            Atom& atom_one=readpdb_one.DNA_atoms[jj][kk];
            if(count(Vsel_chains.begin(),Vsel_chains.end(),atom_one.chainid))
              Vdna_tmp.push_back(atom_one);
          }
          readpdb_one.DNA_atoms[jj]=Vdna_tmp;

          vector<Atom> Vprotein_tmp;
          for(size_t kk=0;kk<readpdb_one.Protein_atoms[jj].size();kk++) {
            Atom& atom_one=readpdb_one.Protein_atoms[jj][kk];
            if(count(Vsel_chains.begin(),Vsel_chains.end(),atom_one.chainid))
              Vprotein_tmp.push_back(atom_one);
          }
          readpdb_one.Protein_atoms[jj]=Vprotein_tmp;
        }
      }
      
      if(nm==0 || readpdb_one.Full_atoms[0].size()==0) {
        string str="No atoms in file "+pdblist[ii]+", skipped!";
        cbs->error->message(str);
        continue;
      }
      vector<Molecule> Vcomplex;
      for(int jj=0;jj<nm;jj++) {
        Molecule mol_one(cbs);
        if(flag_rna) {
          if(readpdb_one.RNA_atoms[jj].size()==0) {
            string str="No RNA atoms in file "+pdblist[ii]+", skipped!";
            cbs->error->message(str);
          }
          else {
            RNAMolecule rna_one(cbs);
            rna_one.set_RNA_molecule(readpdb_one.RNA_atoms[jj]);
            for(size_t kk=0;kk<rna_one.chains.size();kk++)
              mol_one.insert_chain(rna_one.chains[kk]);
          }
        }
        if(flag_protein) {
          //set_protein_molecule
        }
        if(flag_dna) {
          //set_dna_molecule
        }

        if(mol_one.chains.size()==0 || mol_one.chains[0].residues.size()==0) continue;

        if(flag_stand) mol_one.standardize_molecule();
        Vcomplex.push_back(mol_one);
      }

      if(flag_model==0) {
        for(size_t jj=0;jj<Vcomplex.size();jj++) {
          string fstr=outDIR+readpdb_one.PDBid+"-"+to_string(jj+1)+".pdb";
          ofstream outP (fstr.c_str());
          outP<<Vcomplex[jj];
          outP.close(); outP.clear();
        }
      }
      else if(flag_model==1) {
        string fstr=outDIR+readpdb_one.PDBid+".pdb";
        ofstream outP (fstr.c_str());
        for(size_t jj=0;jj<Vcomplex.size();jj++) {
          if(Vcomplex.size()>1)
            outP<<"MODEL     "<<jj+1<<endl;
          outP<<Vcomplex[jj];
          if(Vcomplex.size()>1)
            outP<<"ENDMDL"<<endl;
        }
        outP<<"END"<<endl;
        outP.close(); outP.clear();
      }
    }
  }

  delete cbs;
  
  return 0;
}