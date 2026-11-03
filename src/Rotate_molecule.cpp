#include "basic.h"
#include "error.h"
#include "readfile.h"
#include "initialization.h"
#include "molecule.h"
#include "read_PDB.h"
#include "superpose.h"

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

  if(Mcmd.count("-ref")==0 || Mcmd.count("-tar")==0) {
    cout<<"Wrong commands: Please give reference/target file (-ref/-tar)!"<<endl;
    exit(0);
  }
  string refPDB=Mcmd["-ref"][0];
  string tarPDB=Mcmd["-tar"][0];
  Mcmd.erase("-ref");
  Mcmd.erase("-tar");

  if(Mcmd.count("-o")==0) {
    cout<<"Wrong commands: Please give output file (-o)!"<<endl;
    exit(0);
  }
  string outPDB=Mcmd["-o"][0];
  Mcmd.erase("-o");

  string ref_res="all";
  if(Mcmd.count("-ref_res")) {
    ref_res=Mcmd["-ref_res"][0];
    for(size_t jj=1;jj<Mcmd["-ref_res"].size();jj++)
      ref_res+=" "+Mcmd["-ref_res"][jj];
    Mcmd.erase("-ref_res");
  }
  string tar_res="all";
  if(Mcmd.count("-tar_res")) {
    tar_res=Mcmd["-tar_res"][0];
    for(size_t jj=1;jj<Mcmd["-tar_res"].size();jj++)
      tar_res+=" "+Mcmd["-tar_res"][jj];
    Mcmd.erase("-tar_res");
  }

  int nmodel=-1;
  if(Mcmd.count("-nmodel")) {
    nmodel=stoi(Mcmd["-nmodel"][0]);
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
  readpdb_one.read_PDB(refPDB);
  if(readpdb_one.Full_atoms.size()==0 || readpdb_one.Full_atoms[0].size()==0) {
    string str="No atoms in reference file "+refPDB;
    cbs->error->one(str);
  }
  Molecule mol_ref(cbs);
  mol_ref.set_molecule(readpdb_one.Full_atoms[0],Chain::other);
  
  readpdb_one.clear();
  readpdb_one.read_PDB(tarPDB);
  if(readpdb_one.Full_atoms.size()==0 || readpdb_one.Full_atoms[0].size()==0) {
    string str="No atoms in target file "+tarPDB;
    cbs->error->one(str);
  }
  Molecule mol_tar(cbs);
  mol_tar.set_molecule(readpdb_one.Full_atoms[0],Chain::other);

  if(nmodel<0)
    nmodel=static_cast<int>(readpdb_one.Full_atoms.size());
  else if(nmodel > (int)readpdb_one.Full_atoms.size())
    nmodel=static_cast<int>(readpdb_one.Full_atoms.size());

  struct ResInfo {
    string chainid;
    int resid_begin, resid_end;
  };
  vector<ResInfo> Vrefs;
  vector<string> stoken=split_string(ref_res,",");
  for(size_t ii=0;ii<stoken.size();ii++) {
    if(stoken[0]=="all") {
      for(size_t jj=0;jj<mol_ref.chains.size();jj++) {
        ResInfo seg_i;
        seg_i.chainid=mol_ref.chains[jj].chainid;
        seg_i.resid_begin=mol_ref.chains[jj].residues[0].resid;
        seg_i.resid_end=mol_ref.chains[jj].residues.back().resid;
        Vrefs.push_back(seg_i);
      }
      break;
    }
    else {
      string sseg=stoken[ii];
      ResInfo seg_i;
      size_t pos=sseg.find_first_of(":");
      if(pos==string::npos && sseg.size()>0) {
        seg_i.chainid=sseg;
        seg_i.resid_begin=resid_null;
        seg_i.resid_end=resid_null;
      }
      else if(pos==0) {
        seg_i.chainid="";
        vector<string> stoken1=split_string(sseg,":");
        if(stoken1.size()==1) {
          seg_i.resid_begin=stoi(stoken1[0]);
          seg_i.resid_end=resid_null;
        }
        else if(stoken1.size()==2) {
          seg_i.resid_begin=stoi(stoken1[0]);
          seg_i.resid_end=stoi(stoken1[1]);
        }
      }
      else {
        vector<string> stoken1=split_string(sseg,":");
        if(stoken1.size()==1) {
          seg_i.chainid=stoken1[0];
          seg_i.resid_begin=resid_null;
          seg_i.resid_end=resid_null;
        }
        else if(stoken1.size()==2) {
          seg_i.chainid=stoken1[0];
          seg_i.resid_begin=stoi(stoken1[1]);
          seg_i.resid_end=resid_null;
        }
        else if(stoken1.size()==3) {
          seg_i.chainid=stoken1[0];
          seg_i.resid_begin=stoi(stoken1[1]);
          seg_i.resid_end=stoi(stoken1[2]);
        }
      }
      Vrefs.push_back(seg_i);
    }
  }

  vector<ResInfo> Vtars;
  stoken.clear();
  stoken=split_string(tar_res,",");
  for(size_t ii=0;ii<stoken.size();ii++) {
    if(stoken[0]=="all") {
      for(size_t jj=0;jj<mol_tar.chains.size();jj++) {
        ResInfo seg_i;
        seg_i.chainid=mol_tar.chains[jj].chainid;
        seg_i.resid_begin=mol_tar.chains[jj].residues[0].resid;
        seg_i.resid_end=mol_tar.chains[jj].residues.back().resid;
        Vtars.push_back(seg_i);
      }
      break;
    }
    else {
      string sseg=stoken[ii];
      ResInfo seg_i;
      size_t pos=sseg.find_first_of(":");
      if(pos==string::npos && sseg.size()>0) {
        seg_i.chainid=sseg;
        seg_i.resid_begin=resid_null;
        seg_i.resid_end=resid_null;
      }
      else if(pos==0) {
        seg_i.chainid="";
        vector<string> stoken1=split_string(sseg,":");
        if(stoken1.size()==1) {
          seg_i.resid_begin=stoi(stoken1[0]);
          seg_i.resid_end=resid_null;
        }
        else if(stoken1.size()==2) {
          seg_i.resid_begin=stoi(stoken1[0]);
          seg_i.resid_end=stoi(stoken1[1]);
        }
      }
      else {
        vector<string> stoken1=split_string(sseg,":");
        if(stoken1.size()==1) {
          seg_i.chainid=stoken1[0];
          seg_i.resid_begin=resid_null;
          seg_i.resid_end=resid_null;
        }
        else if(stoken1.size()==2) {
          seg_i.chainid=stoken1[0];
          seg_i.resid_begin=stoi(stoken1[1]);
          seg_i.resid_end=resid_null;
        }
        else if(stoken1.size()==3) {
          seg_i.chainid=stoken1[0];
          seg_i.resid_begin=stoi(stoken1[1]);
          seg_i.resid_end=stoi(stoken1[2]);
        }
      }
      Vtars.push_back(seg_i);
    }
  }

  string ATOMNAME[]={"P","C5'","C4'","O4'","C1'","C2'","C3'","O3'"};
  int nres_ref=0;
  map<string,vector<double> > Matomname_xyz;
  for(size_t ii=0;ii<Vrefs.size();ii++) {
    for(size_t jj=0;jj<mol_ref.chains.size();jj++) {
      if(mol_ref.chains[jj].chainid==Vrefs[ii].chainid) {
        int resid_begin=Vrefs[ii].resid_begin;
        if(resid_begin==resid_null) resid_begin=mol_ref.chains[jj].residues[0].resid;
        int resid_end=Vrefs[ii].resid_end;
        if(resid_end==resid_null) resid_end=mol_ref.chains[jj].residues.back().resid;
        for(size_t kk=0;kk<mol_ref.chains[jj].residues.size();kk++) {
          int &resid_one=mol_ref.chains[jj].residues[kk].resid;
          if(resid_one>=resid_begin && resid_one<=resid_end) {
            nres_ref++;
            for(size_t nn=0;nn<mol_ref.chains[jj].residues[kk].atoms.size();nn++) {
              Atom &atom_one=mol_ref.chains[jj].residues[kk].atoms[nn];
              if(count(begin(ATOMNAME),end(ATOMNAME),atom_one.atomname)) {
                string aname=to_string(nres_ref)+"_"+atom_one.atomname;
                vector<double> xyz=atom_one.get_atom_xyz();
                if(Matomname_xyz.count(aname)==0)
                  Matomname_xyz.emplace(aname,xyz);
              }
            }
          }
        }
        break;
      }
    }
  }
  
  size_t pos_fmt=max(0,(int)outPDB.size()-4);
  string fout_fmt=outPDB.substr(pos_fmt);
  if(fout_fmt!=".pdb") 
    outPDB+=".pdb";
  ofstream outP (outPDB.c_str());
  for(int ii=0;ii<nmodel;ii++) {
    if(ii>0) {
      mol_tar.clear();
      mol_tar.set_molecule(readpdb_one.Full_atoms[ii],Chain::other);
    }
    vector<double> Vxyz_ref,Vxyz_tar;
    int nres_tar=0;
    for(size_t jj=0;jj<Vtars.size();jj++) {
      for(size_t kk=0;kk<mol_tar.chains.size();kk++) {
        if(mol_tar.chains[kk].chainid==Vtars[jj].chainid) {
          int resid_begin=Vtars[jj].resid_begin;
          if(resid_begin==resid_null) resid_begin=mol_tar.chains[kk].residues[0].resid;
          int resid_end=Vtars[jj].resid_end;
          if(resid_end==resid_null) resid_end=mol_tar.chains[kk].residues.back().resid;
          for(size_t nn=0;nn<mol_tar.chains[kk].residues.size();nn++) {
            int &resid_one=mol_tar.chains[kk].residues[nn].resid;
            if(resid_one>=resid_begin && resid_one<=resid_end) {
              nres_tar++;
              for(size_t mm=0;mm<mol_tar.chains[kk].residues[nn].atoms.size();mm++) {
                Atom &atom_one=mol_tar.chains[kk].residues[nn].atoms[mm];
                if(count(begin(ATOMNAME),end(ATOMNAME),atom_one.atomname)) {
                  string aname=to_string(nres_tar)+"_"+atom_one.atomname;
                  if(Matomname_xyz.count(aname)) {
                    vector<double> &xyz_ref=Matomname_xyz[aname];
                    Vxyz_ref.insert(Vxyz_ref.end(),xyz_ref.begin(),xyz_ref.end());
                    vector<double> xyz_tar=atom_one.get_atom_xyz();
                    Vxyz_tar.insert(Vxyz_tar.end(),xyz_tar.begin(),xyz_tar.end());
                  }
                }
              }
            }
          }
          break;
        }
      }
    }

    if(nres_ref!=nres_tar) {
      string str="Numbers of residues for matching are different --> ";
      str+=to_string(nres_ref)+"/"+to_string(nres_tar);
      cbs->error->one(str);
    }
    
    Superpose sp_one;
    sp_one.set_x1x2(Vxyz_ref,Vxyz_tar);
    sp_one.calculate_rotation_matrix();
    double rmsd_one=sp_one.get_RMSD();

    string str="RMSD = "+to_string(rmsd_one)+" for Model #"+to_string(ii+1);
    cbs->error->message(str);
    cbs->error->log(str);
    
    vector<double> Vxyz_mol=mol_tar.get_molecule_xyz();
    vector<double> Vxyz_mol_rotated=sp_one.rotate_x2(Vxyz_mol);
    mol_tar.reset_molecule_xyz(Vxyz_mol_rotated);

    if(nmodel>1)
      outP<<"MODEL     "<<ii+1<<endl;
    outP<<mol_tar;
    if(nmodel>1)
      outP<<"ENDMDL"<<endl;
  }
  if(nmodel>1)
    outP<<"END"<<endl;
  outP.close(); outP.clear();

  delete cbs;
  
  return 0;
}