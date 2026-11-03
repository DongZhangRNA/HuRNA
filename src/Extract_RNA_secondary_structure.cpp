#include "basic.h"
#include "error.h"
#include "readfile.h"
#include "initialization.h"
#include "RNA_molecule.h"
#include "RNA_basepair.h"
#include "read_PDB.h"

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
  string inFile,outFile;

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

  int flag_break_check=0;
  if(Mcmd.count("-break_check")) {
    flag_break_check=1;
    Mcmd.erase("-break_check");
  }

  string tool="dssr";
  if(Mcmd.count("-tool")) {
    tool=Mcmd["-tool"][0];
    Mcmd.erase("-tool");
  }
  if(tool!="dssr" && tool!="mcannotate") {
    cout<<"Wrong commands: only [dssr, mcannotate] allowed for -tool!"<<endl;
    exit(0);
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
  readpdb_one.read_PDB(inFile);
  RNAMolecule rna_one(cbs);

  cbs->init->RNA_check_chain_break=flag_break_check;

  if(readpdb_one.RNA_atoms.size()==0 || readpdb_one.RNA_atoms[0].size()==0) {
    string str="No RNA atoms in "+inFile+"!";
    cbs->error->one(str);
  }

  rna_one.set_RNA_molecule(readpdb_one.RNA_atoms[0]);
  
  string seq=rna_one.get_sequences();
  string secstr;
  if(tool=="dssr")
    secstr=rna_one.extract_secondary_structure();
  else if( tool=="mcannotate")
    secstr=rna_one.extract_secondary_structure_MCAnnotate();

  ofstream outF (outFile.c_str());
  outF<<seq<<endl;
  outF<<secstr<<endl;
  outF.close(); outF.clear();

  delete cbs;
  
  return 0;
}