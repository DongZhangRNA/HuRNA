#include "basic.h"
#include "error.h"
#include "readfile.h"
#include "initialization.h"
#include "molecule.h"
#include "read_PDB.h"
#include "read_lammps.h"

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

  string fintrj="";
  if(Mcmd.count("-i")) {
    fintrj=Mcmd["-i"][0];
    Mcmd.erase("-i");
  }
  if(fintrj.length()==0) {
    cout<<"Wrong commands: No input trj file (-i)!"<<endl;
    exit(0);
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
  if(fout.find(".pdb")==string::npos) {
    fout+=".pdb";
  }
  string fref;
  if(Mcmd.count("-ref")) {
    fref=Mcmd["-ref"][0];
    Mcmd.erase("-ref");
  }
  else {
    cout<<"Wrong commands: No input for residue info (-ref)!"<<endl;
    exit(0);
  }

  int step_begin=0;
  if(Mcmd.count("-b")) {
    step_begin=stoi(Mcmd["-b"][0]);
    Mcmd.erase("-b");
  }
  int step_end=-1;
  if(Mcmd.count("-e")) {
    step_end=stoi(Mcmd["-e"][0]);
    Mcmd.erase("-e");
  }

  int flag_sep=0;
  if(Mcmd.count("-sep"))  {
    flag_sep=1;
    Mcmd.erase("-sep");
  }

  int flag_center=0;
  if(Mcmd.count("-center")) {
    flag_center=1;
    Mcmd.erase("-center");
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

  ReadPDB readpdb_one(cbs);
  readpdb_one.read_PDB(fref);
  vector<Atom> Vatoms=readpdb_one.Full_atoms[0];

  Molecule mol_ref(cbs);
  mol_ref.set_molecule(Vatoms,Chain::cg);
  mol_ref.set_atomid(0);

  ReadLAMMPS readlmp_one(cbs);
  readlmp_one.read_lammpstrj(fintrj,step_begin,step_end);
  size_t nmodel=readlmp_one.VVatoms[0].size();

  if(flag_sep==0) {
    ofstream outR (fout.c_str());
    for(size_t ii=0;ii<nmodel;ii++) {
      vector<double> Vxyz;
      for(size_t jj=0;jj<readlmp_one.VVatoms[0][ii].size();jj++) {
        vector<double> xyz=readlmp_one.VVatoms[0][ii][jj].get_atom_xyz();
        Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
      }
      mol_ref.reset_molecule_xyz(Vxyz);
      if(flag_center) {
        vector<double> Vxyz_cog=mol_ref.center_of_geometry();
        mol_ref.move_molecule(-Vxyz_cog[0],-Vxyz_cog[1],-Vxyz_cog[2]);
      }
      if(nmodel>1)
        outR<<"MODEL     "<<ii+1<<endl;
      outR<<mol_ref;
      if(nmodel>1)
        outR<<"ENDMDL"<<endl;
    }
    outR<<"END"<<endl;
    outR.close(); outR.clear();  
  }
  else {
    for(size_t ii=0;ii<nmodel;ii++) {
      vector<double> Vxyz;
      for(size_t jj=0;jj<readlmp_one.VVatoms[0][ii].size();jj++) {
        vector<double> xyz=readlmp_one.VVatoms[0][ii][jj].get_atom_xyz();
        Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
      }
      mol_ref.reset_molecule_xyz(Vxyz);
      if(flag_center) {
        vector<double> Vxyz_cog=mol_ref.center_of_geometry();
        mol_ref.move_molecule(-Vxyz_cog[0],-Vxyz_cog[1],-Vxyz_cog[2]);
      }
      string fstr=fout.substr(0,fout.size()-4);
      fstr+="_"+to_string(readlmp_one.VVsteps[0][ii])+".pdb";
      ofstream outR (fstr.c_str());
      outR<<mol_ref;
      outR.close(); outR.clear();
    }
  }

  delete cbs;
  
  return 0;
}