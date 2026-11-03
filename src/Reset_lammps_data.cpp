#include "basic.h"
#include "error.h"
#include "initialization.h"
#include "read_lammps.h"
#include "prepare_lammps.h"

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

  string fdata="";
  if(Mcmd.count("-i")) {
    fdata=Mcmd["-i"][0];
    Mcmd.erase("-i");
  }
  else {
    cout<<"No input lammps data file (-i)!"<<endl;
    exit(0);
  }
  
  string ftrj="";
  if(Mcmd.count("-trj")) {
    ftrj=Mcmd["-trj"][0];
    Mcmd.erase("-trj");
  }
  else {
    cout<<"No input lammps trj file (-trj)!"<<endl;
    exit(0);
  }

  string fout="";
  if(Mcmd.count("-o")) {
    fout=Mcmd["-o"][0];
    Mcmd.erase("-o");
  }
  else {
    cout<<"No output file (-o)!"<<endl;
    exit(0);
  }

  int step=-1;
  if(Mcmd.count("-step")) {
    step=stoi(Mcmd["-step"][0]);
    Mcmd.erase("-step");
  }

  int fcenter=0;
  if(Mcmd.count("-center")) {
    fcenter=1;
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
  }
  cbs->initialization();

  vector<string> Vdata;
  string sline;
  ifstream inF (fdata.c_str());
  if(!inF.is_open()) {
    string str="Cannot find file "+fdata;
    cbs->error->one(str);
  }
  while(getline(inF,sline)) {
    Vdata.push_back(sline);
  }
  inF.close(); inF.clear();

  ReadLAMMPS readlmp_one(cbs);
  readlmp_one.read_lammpstrj(ftrj);

  vector<Atom> Vatoms;
  if(step==-1) {
    Vatoms=readlmp_one.VVatoms[0].back();
  }
  else {
    int flag_check=0;
    for(size_t ii=0;ii<readlmp_one.VVsteps[0].size();ii++) {
      if(readlmp_one.VVsteps[0][ii]==step) {
        flag_check=1;
        Vatoms=readlmp_one.VVatoms[0][ii];
        break;
      }  
    }
    if(flag_check==0) {
      string str="Cannot find Step = "+to_string(step)+" in trj file "+ftrj;
      cbs->error->one(str);
    }
  }

  if(fcenter) {
    double x0cm=0.0,y0cm=0.0,z0cm=0.0;
    size_t nline=Vdata.size();
    if(nline>25)
      nline=25;
    for(size_t ii=0;ii<20;ii++) {
      string sline=Vdata[ii];
      if(sline.size()<2) continue;
      vector<string> stoken=split_string(sline);
      if(stoken.size()!=4) continue;
      if(sline.find("xlo")!=string::npos && sline.find("xhi")!=string::npos) {
        x0cm=(stod(stoken[0])+stod(stoken[1]))/2.0;
      }
      else if(sline.find("ylo")!=string::npos && sline.find("yhi")!=string::npos) {
        y0cm=(stod(stoken[0])+stod(stoken[1]))/2.0;
      }
      else if(sline.find("zlo")!=string::npos && sline.find("zhi")!=string::npos) {
        z0cm=(stod(stoken[0])+stod(stoken[1]))/2.0;
      }
    }

    double x1cm=0.0, y1cm=0.0, z1cm=0.0;
    size_t natom=Vatoms.size();
    for(size_t ii=0;ii<natom;ii++) {
      x1cm+=Vatoms[ii].x;
      y1cm+=Vatoms[ii].y;
      z1cm+=Vatoms[ii].z;
    }
    x1cm/=natom;
    y1cm/=natom;
    z1cm/=natom;

    double delx=x0cm-x1cm;
    double dely=y0cm-y1cm;
    double delz=z0cm-z1cm;

    for(size_t ii=0;ii<natom;ii++) {
      Vatoms[ii].move_atom(delx,dely,delz);
    }
  }

  PrepareLAMMPS prelmp_one(cbs);
  prelmp_one.update_data_xyz(Vdata,Vatoms);

  ofstream outF (fout.c_str());
  for(size_t ii=0;ii<Vdata.size();ii++)
    outF<<Vdata[ii]<<endl;
  outF.close(); outF.clear();

  delete cbs;

  return 0;
}