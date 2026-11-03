#include "basic.h"
#include "error.h"
#include "readfile.h"
#include "initialization.h"
#include "read_PDB.h"
#include "RNA_molecule.h"
#include "RNA_coarse_grained.h"
#include "read_lammps.h"
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

  string fin;
  if(Mcmd.count("-i")) {
    fin=Mcmd["-i"][0];
    Mcmd.erase("-i");
  }
  else {
    cout<<"Wrong commands: No input CGMD trj (-i)!"<<endl;
    exit(0);
  }

  string tar_res="all";
  if(Mcmd.count("-sel_res")) {
    tar_res=Mcmd["-sel_res"][0];
    for(size_t jj=1;jj<Mcmd["-sel_res"].size();jj++)
      tar_res+=" "+Mcmd["-sel_res"][jj];
    Mcmd.erase("-sel_res");
  }
  string fres_info="";
  if(Mcmd.count("-res_info")) {
    fres_info=Mcmd["-res_info"][0];
    Mcmd.erase("-res_info");
  }
  else {
    if(tar_res!="all") {
      cout<<"Wrong commands: Please provide residue information file first (-res_info)!"<<endl;
      exit(0);
    }
  }

  string fref="";
  if(Mcmd.count("-ref")) {
    fref=Mcmd["-ref"][0];
    Mcmd.erase("-ref");
  }
  string ref_res="all";
  if(Mcmd.count("-ref_res")) {
    ref_res=Mcmd["-ref_res"][0];
    for(size_t jj=1;jj<Mcmd["-ref_res"].size();jj++)
      ref_res+=" "+Mcmd["-ref_res"][jj];
    Mcmd.erase("-ref_res");
  }
  string fref_format="aa";
  if(Mcmd.count("-rf")) {
    fref_format=Mcmd["-rf"][0];
    Mcmd.erase("-rf");
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

  ReadLAMMPS readlmp_one(cbs);
  readlmp_one.read_lammpstrj(fin,step_begin,step_end);

  ofstream outF (fout.c_str());
  outF<<"#Step RMSD"<<endl;

  if(fref.size()==0 && tar_res=="all") {
    vector<double> Vxyz_ref;
    for(size_t ii=0;ii<readlmp_one.VVatoms[0][0].size();ii++) {
      vector<double> xyz=readlmp_one.VVatoms[0][0][ii].get_atom_xyz();
      Vxyz_ref.insert(Vxyz_ref.end(),xyz.begin(),xyz.end());
    }
    for(size_t ii=0;ii<readlmp_one.VVatoms[0].size();ii++) {
      vector<double> Vxyz_tar;
      for(size_t jj=0;jj<readlmp_one.VVatoms[0][ii].size();jj++) {
        vector<double> xyz=readlmp_one.VVatoms[0][ii][jj].get_atom_xyz();
        Vxyz_tar.insert(Vxyz_tar.end(),xyz.begin(),xyz.end());
      }
      Superpose sp_one;
      sp_one.set_x1x2(Vxyz_ref,Vxyz_tar);
      double rmsd=sp_one.get_RMSD();
      outF<<readlmp_one.VVsteps[0][ii]<<" "<<rmsd<<endl;
    }
    outF.close(); outF.clear();

    delete cbs;
  
    return 0;
  }

  Molecule mol_ref(cbs);
  if(fref.size()!=0) {
    ReadPDB readpdb_one(cbs);
    readpdb_one.read_PDB(fref);
    if(readpdb_one.RNA_atoms.size()==0 || readpdb_one.RNA_atoms[0].size()==0) {
      string str="No RNA atoms in "+fref;
      cbs->error->one(str);
    }
    if(fref_format=="aa") {
      RNAMolecule rna_one(cbs);
      rna_one.set_RNA_molecule(readpdb_one.RNA_atoms[0]);
      RNACoarseGrained rnacg_one(cbs);
      if(cbs->init->RNA_version=="latest")
        rnacg_one.convert_RNA_to_CG(rna_one);
      else if(cbs->init->RNA_version=="v2")
        rnacg_one.convert_RNA_to_CG_v2(rna_one);
      mol_ref=rnacg_one.rna_cg;
    }
    else if(fref_format=="cg") {
      mol_ref.set_molecule(readpdb_one.RNA_atoms[0],Chain::cg);
    }
    else {
      string str="Unknown input format (-rf) for reference strucutre "+fref;
      cbs->error->one(str);
    }
  }

  Molecule mol_tar(cbs);
  vector<Atom> Vatoms_one=readlmp_one.VVatoms[0][0];
  if(fres_info.size()>0) 
    readlmp_one.recover_residue_info(fres_info,Vatoms_one);
  mol_tar.set_molecule(Vatoms_one,Chain::cg);
  
  if(fref.size()==0) {
    mol_ref=mol_tar;
    ref_res=tar_res;
  }

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

  if(fres_info.size()==0 && tar_res=="all") {
    vector<double> Vxyz_ref;
    size_t natom_ref=0;
    for(size_t ii=0;ii<Vresidues_ref.size();ii++) {
      for(size_t jj=0;jj<Vresidues_ref[ii].atoms.size();jj++) {
        natom_ref++;
        vector<double> xyz=Vresidues_ref[ii].atoms[jj].get_atom_xyz();
        Vxyz_ref.insert(Vxyz_ref.end(),xyz.begin(),xyz.end());
      }
    }
    for(size_t ii=0;ii<readlmp_one.VVatoms[0].size();ii++) {
      size_t natom_tar=readlmp_one.VVatoms[0][ii].size();
      if(natom_ref!=natom_tar) {
        string str="Number of atoms is not matched --> "+to_string(natom_ref)+"(ref) / "+to_string(natom_tar)+"(tar)";
        cbs->error->one(str);
      }
      vector<double> Vxyz_tar;
      for(size_t jj=0;jj<natom_tar;jj++) {
        vector<double> xyz=readlmp_one.VVatoms[0][ii][jj].get_atom_xyz();
        Vxyz_tar.insert(Vxyz_tar.end(),xyz.begin(),xyz.end());
      }
      Superpose sp_one;
      sp_one.set_x1x2(Vxyz_ref,Vxyz_tar);
      double rmsd=sp_one.get_RMSD();
      outF<<readlmp_one.VVsteps[0][ii]<<" "<<rmsd<<endl;
    }
    outF.close(); outF.clear();

    delete cbs;
  
    return 0;
  }
  
  vector<Residue> Vresidues_tar;
  if(tar_res=="all") {
    for(size_t ii=0;ii<mol_tar.chains.size();ii++) {
      for(size_t jj=0;jj<mol_tar.chains[ii].residues.size();jj++) {
        Vresidues_tar.push_back(mol_tar.chains[ii].residues[jj]);
      }
    }
  }
  else {
    vector<ResInfo> Vtars;
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

    for(size_t ii=0;ii<Vtars.size();ii++) {
      for(size_t jj=0;jj<mol_tar.chains.size();jj++) {
        string chainid_one=mol_tar.chains[jj].chainid;
        if(Vtars[ii].chainid==chainid_one || (Vtars[ii].chainid.length()==0 && chainid_one.length()==0)) {
          if(Vtars[ii].resid_begin==resid_null) 
            Vtars[ii].resid_begin=mol_tar.chains[jj].residues[0].resid;
          if(Vtars[ii].resid_end==resid_null) 
            Vtars[ii].resid_end=mol_tar.chains[jj].residues.back().resid;
          for(size_t kk=0;kk<mol_tar.chains[jj].residues.size();kk++) {
            Residue &residue_one=mol_tar.chains[jj].residues[kk];
            int resid_one=residue_one.resid;
            if(resid_one>=Vtars[ii].resid_begin && resid_one<=Vtars[ii].resid_end)
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

  vector<double> Vxyz_ref;
  vector<int> Vids_used;
  for(size_t ii=0;ii<Vresidues_ref.size();ii++) {
    map<string,vector<double> > Matomname_xyz;
    for(size_t jj=0;jj<Vresidues_ref[ii].atoms.size();jj++) {
      string atomname_one=Vresidues_ref[ii].atoms[jj].atomname;
      vector<double> xyz_one=Vresidues_ref[ii].atoms[jj].get_atom_xyz();
      Matomname_xyz.emplace(atomname_one,xyz_one);
    }
    for(size_t jj=0;jj<Vresidues_tar[ii].atoms.size();jj++) {
      string atomname_one=Vresidues_tar[ii].atoms[jj].atomname;
      if(Matomname_xyz.count(atomname_one)) {
        vector<double>& xyz_one=Matomname_xyz[atomname_one];
        Vxyz_ref.insert(Vxyz_ref.end(),xyz_one.begin(),xyz_one.end());
        Vids_used.push_back(Vresidues_tar[ii].atoms[jj].id);
      }
    }
  }
  vector<size_t> Vindexs_used;
  for(size_t ii=0;ii<Vids_used.size();ii++) {
    int flag_check=0;
    size_t index_i=0;
    for(size_t jj=0;jj<readlmp_one.VVatoms[0][0].size();jj++) {
      if(readlmp_one.VVatoms[0][0][jj].id==Vids_used[ii]) {
        index_i=jj;
        flag_check=1;
        break;
      }
    }
    if(flag_check==0) {
      string str="Cannot find atom id="+to_string(Vids_used[ii])+" in CGMD trj!";
      cbs->error->one(str);
    }
    Vindexs_used.push_back(index_i);
  }
  for(size_t ii=0;ii<readlmp_one.VVatoms[0].size();ii++) {
    vector<double> Vxyz_tar;
    for(size_t jj=0;jj<Vindexs_used.size();jj++) {
      vector<double> xyz=readlmp_one.VVatoms[0][ii][Vindexs_used[jj]].get_atom_xyz();
      Vxyz_tar.insert(Vxyz_tar.end(),xyz.begin(),xyz.end());
    }
    Superpose sp_one;
    sp_one.set_x1x2(Vxyz_ref,Vxyz_tar);
    double rmsd=sp_one.get_RMSD();
    outF<<readlmp_one.VVsteps[0][ii]<<" "<<rmsd<<endl;
  }
  outF.close(); outF.clear();

  delete cbs;
  
  return 0;
}