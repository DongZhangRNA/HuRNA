#include "basic.h"
#include "constant.h"
#include "initialization.h"
#include "error.h"
#include "residue.h"
#include "superpose.h"

using namespace CBS_NS;

const int ATOMID_MAX=99999;
const double Backbone_RMSD_CUT=0.5;
const double RNA_base_RMSD_CUT=0.2;

Residue::Residue(Cobmms *cbs) : Pointers(cbs) {}

Residue::Residue(Cobmms *cbs, string cid_one, int resid_one, string resname_one) : Pointers(cbs) 
{
  Residue::set_residue(cid_one,resid_one,resname_one);
}

Residue::Residue(const Residue &other) : Pointers(other.cbs) {
  chainid=other.chainid;
  resid=other.resid;
  resname=other.resname;
  atoms=other.atoms;
  atomxyz=other.atomxyz;
}

Residue& Residue::operator =(const Residue &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    chainid=other.chainid;
    resid=other.resid;
    resname=other.resname;
    atoms=other.atoms;
    atomxyz=other.atomxyz;
  }

  return *this;
}

Residue::~Residue() {
  chainid="";
  resid=-1;
  resname="";
  atoms.clear();
  atoms.shrink_to_fit();
  atomxyz.clear();
}

void Residue::set_residue(string cid_one, int resid_one, string resname_one) {
  chainid=cid_one;
  resid=resid_one;
  resname=resname_one;
  atoms.clear();
  atomxyz.clear();

  /*
  if(init->RESIDUES_ATOMS.count(resname)==0) {
    string str="residue.cpp: Undefined atoms list for residue "+resname;
    error->warning(str);
  }
  else {
    vector<string> Vatoms=init->RESIDUES_ATOMS[resname];
    for(size_t ii=0;ii<Vatoms.size();ii++) {
      vector<double> xyz;
      for(int jj=0;jj<3;jj++)
        xyz.push_back(0.0);
      add_atom(Vatoms[ii],xyz);
    }
  }
  */
}

void Residue::set_residue(string res_one) {
  vector<string> stoken=split_string(res_one,"/");
  if(stoken.size()!=3) {
    string str="residue.cpp: Wrong residue name(chainid/resid/resname) --> "+res_one;
    error->one(str);
  }

  chainid=stoken[0];
  resid=stoi(stoken[1]);
  resname=stoken[2];
  atoms.clear();
  atomxyz.clear();

  //string cid_one=stoken[0];
  //int resid_one=stoi(stoken[1]);
  //string resname_one=stoken[2];
  //set_residue(cid_one,resid_one,resname_one);
}

void Residue::set_residue(const Residue &res_one) {
  chainid=res_one.chainid;
  resid=res_one.resid;
  resname=res_one.resname;
  atoms=res_one.atoms;
  atomxyz=res_one.atomxyz;
}

void Residue::set_residue(const vector<Atom> &Vatom_one) {
  atoms.clear();
  atomxyz.clear();
  chainid=Vatom_one[0].chainid;
  resid=Vatom_one[0].resid;
  resname=Vatom_one[0].resname;
  
  for(size_t ii=0;ii<Vatom_one.size();ii++) {
    string name_one=Vatom_one[ii].atomname;
    if(name_one=="NULL") continue;
    if(Vatom_one[ii].atomalt!="" && Vatom_one[ii].atomalt!="A" && Vatom_one[ii].atomalt!=" ") continue;
    if(Vatom_one[ii].chainid!=chainid || Vatom_one[ii].resid!=resid || Vatom_one[ii].resname!=resname) {
      string str="residue.cpp: Residue name is not match in set_residue --> "+chainid+"/"+to_string(resid)+"/"+resname;
      str+=" : "+Vatom_one[ii].chainid+"/"+to_string(Vatom_one[ii].resid)+"/"+Vatom_one[ii].resname;
      error->one(str);
    }
    if(atomxyz.count(name_one)==0) {
      atoms.push_back(Vatom_one[ii]);
      vector<double> xyz=Vatom_one[ii].get_atom_xyz();
      atomxyz.emplace(name_one,xyz);
    }
    else {
      string str="residue.cpp: More than one "+name_one+" atom in set_residue --> Ignored!";
      error->warning(str);
    }
  }
}

/*
void Residue::set_RNA_residue(vector<Atom> &Vnt) {
  atoms.clear();
  atomxyz.clear();
  chainid=Vnt[0].chainid;
  resid=Vnt[0].resid;
  resname=Vnt[0].resname;

  for(size_t ii=0;ii<Vnt.size();ii++) {
    string name_one=Vnt[ii].atomname;
    if(name_one=="NULL") continue;
    if(Vnt[ii].atomalt!="" && Vnt[ii].atomalt!="A" && Vnt[ii].atomalt!=" ") continue;
    if(Vnt[ii].chainid!=chainid || Vnt[ii].resid!=resid || Vnt[ii].resname!=resname) {
      string str="residue.cpp: Residue name is not match in set_residue --> ";
      str+=chainid+"/"+to_string(resid)+"/"+resname+" : "+Vnt[ii].chainid+"/"+to_string(Vnt[ii].resid)+"/"+Vnt[ii].resname;
      error->one(str);
    }
    Atom atom1=Vnt[ii];
    vector<double> xyz=atom1.get_atom_xyz();
    atomxyz.insert(pair<string,vector<double> >(name_one,xyz));
  }

  for(auto& it : atomxyz) {
    string aname=it.first;
    vector<double> xyz=it.second;
    add_atom(aname,xyz);
  }
}
*/

string Residue::get_residue_name() const {
  string name=chainid+"/"+to_string(resid);
  if(atoms.size()>0 && !isspace(atoms[0].iCode)) name+=to_string((char)atoms[0].iCode);
  name+="/"+resname;
  if(chainid.size()==0 && resname.size()==0) name=ResidueConst::RESIDUE_NULL;

  return name;
}

void Residue::reset_residue_name(string name_one) {
  vector<string> stoken=split_string(name_one,"/");
  if(stoken.size()!=3) {
    string str="residue.cpp: Wrong residue name(chainid/resid/resname) --> "+name_one;
    error->one(str);
  }

  chainid=stoken[0];
  resid=stoi(stoken[1]);
  resname=stoken[2];
  for(size_t ii=0;ii<atoms.size();ii++) {
    atoms[ii].chainid=chainid;
    atoms[ii].resid=resid;
    atoms[ii].resname=resname;
    atoms[ii].atomalt="";
    atoms[ii].iCode=' ';
  }
}

void Residue::reset_residue_name(string chainid_one, int resid_one, string resname_one) {
  chainid=chainid_one;
  resid=resid_one;
  resname=resname_one;
  for(size_t ii=0;ii<atoms.size();ii++) {
    atoms[ii].chainid=chainid;
    atoms[ii].resid=resid;
    atoms[ii].resname=resname;
    atoms[ii].atomalt="";
    atoms[ii].iCode=' ';
  }
}

void Residue::recovery_RNA_backbone_atoms() {
  map<string,vector<double> > BB_XYZ;
  vector<Atom> VBB=init->RESIDUES_XYZ["rBB"];
  for(auto atom : VBB) {
    string name=atom.atomname;
    vector<double> xyz=atom.get_atom_xyz();
    BB_XYZ.emplace(name,xyz);
  }

  string res_name=Residue::get_residue_name();
  if((!Residue::check_atom("C5'") && !Residue::check_atom("O3'")) || !Residue::check_atom("C4'") || !Residue::check_atom("C3'")) { //Basic atoms C4', C3', C5'/O3'
    string str="residue.cpp: Failed to recovery backbone atoms --> One or more basic atoms (C4', C3', C5'/O3') is missing for residue "+res_name;
    error->one(str);
  }  

  enum {O3p,C5p};
  int flag_O3p=O3p;
  if(!Residue::check_atom("O3'")) flag_O3p=C5p;

  if(!Residue::check_atom("O3'")) { //add atom O3'
    vector<double> Vxyz,Vxyz0;
    vector<double> xyz;
    string ATOMR[]={"C5'","C4'","C3'"};
    for(auto atom : ATOMR) {
      if(atomxyz.count(atom)) {
        xyz=atomxyz[atom];
        Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
        xyz=BB_XYZ[atom];
        Vxyz0.insert(Vxyz0.end(),xyz.begin(),xyz.end());
      }
    }
    Superpose sp(Vxyz,Vxyz0);
    double rmsd=sp.get_RMSD();
    if(rmsd>Backbone_RMSD_CUT) {
      string str="residue.cpp: RMSD="+to_string(rmsd)+" Angstrom when recovery atom O3' for residue "+res_name;
      error->warning(str);
    }
    xyz=BB_XYZ["O3'"];
    sp.reset_x2(xyz);
    xyz=sp.rotate_x2();
    Residue::add_atom("O3'",xyz);
  }

  if(!Residue::check_atom("C5'")) { //add atom C5'
    vector<double> Vxyz,Vxyz0;
    vector<double> xyz;
    string ATOMR[]={"C4'","C3'","O3'","O5'"};
    for(auto atom : ATOMR) {
      if(atomxyz.count(atom)) {
        xyz=atomxyz[atom];
        Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
        xyz=BB_XYZ[atom];
        Vxyz0.insert(Vxyz0.end(),xyz.begin(),xyz.end());
      }
    }
    Superpose sp(Vxyz,Vxyz0);
    double rmsd=sp.get_RMSD();
    if(rmsd>Backbone_RMSD_CUT) {
      string str="residue.cpp: RMSD="+to_string(rmsd)+" Angstrom when recovery atom C5' for residue "+res_name;
      error->warning(str);
    }
    xyz=BB_XYZ["C5'"];
    sp.reset_x2(xyz);
    xyz=sp.rotate_x2();
    Residue::add_atom("C5'",xyz);
  }

  if(!Residue::check_atom("O5'")) { //add atom O5'
    vector<double> Vxyz,Vxyz0;
    vector<double> xyz;
    string ATOMR[]={"C5'","C4'","C3'","P"};
    for(auto atom : ATOMR) {
      if(atomxyz.count(atom)) {
        xyz=atomxyz[atom];
        Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
        xyz=BB_XYZ[atom];
        Vxyz0.insert(Vxyz0.end(),xyz.begin(),xyz.end());
      }
    }
    Superpose sp(Vxyz,Vxyz0);
    double rmsd=sp.get_RMSD();
    if(rmsd>Backbone_RMSD_CUT) {
      string str="residue.cpp: RMSD="+to_string(rmsd)+" Angstrom when recovery atom O5' for residue "+res_name;
      error->warning(str);
    }
    xyz=BB_XYZ["O5'"];
    sp.reset_x2(xyz);
    xyz=sp.rotate_x2();
    Residue::add_atom("O5'",xyz);
  }

  if(!Residue::check_atom("P") || !Residue::check_atom("OP1") || !Residue::check_atom("OP2")) { //add atom P,OP1,OP2
    vector<double> Vxyz,Vxyz0;
    vector<double> xyz;
    string ATOMR[]={"O5'","C5'","C4'","P","OP1","OP2"};
    for(auto atom : ATOMR) {
      if(atomxyz.count(atom)) {
        xyz=atomxyz[atom];
        Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
        xyz=BB_XYZ[atom];
        Vxyz0.insert(Vxyz0.end(),xyz.begin(),xyz.end());
      }
    }
    Superpose sp(Vxyz,Vxyz0);
    double rmsd=sp.get_RMSD();
    if(rmsd>Backbone_RMSD_CUT) {
      string str="residue.cpp: RMSD="+to_string(rmsd)+" Angstrom when recovery atoms P/OP1/OP2 for residue "+res_name;
      error->warning(str);
    }
    if(atomxyz.count("P")==0) {
      xyz=BB_XYZ["P"];
      sp.reset_x2(xyz);
      xyz=sp.rotate_x2();
      Residue::add_atom("P",xyz);
    }
    if(atomxyz.count("OP1")==0) {
      xyz=BB_XYZ["OP1"];
      sp.reset_x2(xyz);
      xyz=sp.rotate_x2();
      Residue::add_atom("OP1",xyz);
    }
    if(atomxyz.count("OP2")==0) {
      xyz=BB_XYZ["OP2"];
      sp.reset_x2(xyz);
      xyz=sp.rotate_x2();
      Residue::add_atom("OP2",xyz);
    }
  }

  if(!Residue::check_atom("C2'") && !Residue::check_atom("O4'") && !Residue::check_atom("C1'")) { //add atom C2'/O4'
    vector<double> Vxyz,Vxyz0;
    if(flag_O3p==O3p) { //add atom C2'
      vector<double> xyz;
      string ATOMR[]={"C4'","C3'","O3'"};
      for(auto atom : ATOMR) {
        xyz=atomxyz[atom];
        Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
        xyz=BB_XYZ[atom];
        Vxyz0.insert(Vxyz0.end(),xyz.begin(),xyz.end());
      }
      Superpose sp(Vxyz,Vxyz0);
      double rmsd=sp.get_RMSD();
      if(rmsd>Backbone_RMSD_CUT) {
        string str="residue.cpp: RMSD="+to_string(rmsd)+" Angstrom when recovery atom C2' for residue "+res_name;
        error->warning(str);
      }
      xyz=BB_XYZ["C2'"];
      sp.reset_x2(xyz);
      xyz=sp.rotate_x2();
      Residue::add_atom("C2'",xyz);
    }
    else if(flag_O3p==C5p) {  //add atom O4'
      vector<double> xyz;
      string ATOMR[]={"C4'","C3'","C5'"};
      for(auto atom : ATOMR) {
        xyz=atomxyz[atom];
        Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
        xyz=BB_XYZ[atom];
        Vxyz0.insert(Vxyz0.end(),xyz.begin(),xyz.end());
      }
      Superpose sp(Vxyz,Vxyz0);
      double rmsd=sp.get_RMSD();
      if(rmsd>Backbone_RMSD_CUT) {
        string str="residue.cpp: RMSD="+to_string(rmsd)+" Angstrom when recovery atom O4' for residue "+res_name;
        error->warning(str);
      }
      xyz=BB_XYZ["O4'"];
      sp.reset_x2(xyz);
      xyz=sp.rotate_x2();
      Residue::add_atom("O4'",xyz);
    }
  }

  if(!Residue::check_atom("C2'") || !Residue::check_atom("O4'") || !Residue::check_atom("C1'")) { //add atom O4'/C2'/C1' in sugar ring
    vector<double> Vxyz,Vxyz0;
    vector<double> xyz;
    string ATOMR[]={"C4'","C3'","O4'","C2'","C1'"};
    for(auto atom : ATOMR) {
      if(atomxyz.count(atom)) {
        xyz=atomxyz[atom];
        Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
        xyz=BB_XYZ[atom];
        Vxyz0.insert(Vxyz0.end(),xyz.begin(),xyz.end());
      }
    }
    Superpose sp(Vxyz,Vxyz0);
    double rmsd=sp.get_RMSD();
    if(rmsd>Backbone_RMSD_CUT) {
      string str="residue.cpp: RMSD="+to_string(rmsd)+" Angstrom when recovery atoms O4'/C2'/C1' for residue "+res_name;
      error->warning(str);
    }
    if(atomxyz.count("C1'")==0) {
      xyz=BB_XYZ["C1'"];
      sp.reset_x2(xyz);
      xyz=sp.rotate_x2();
      Residue::add_atom("C1'",xyz);
    }
    if(atomxyz.count("C2'")==0) {
      xyz=BB_XYZ["C2'"];
      sp.reset_x2(xyz);
      xyz=sp.rotate_x2();
      Residue::add_atom("C2'",xyz);
    }
    if(atomxyz.count("O4'")==0) {
      xyz=BB_XYZ["O4'"];
      sp.reset_x2(xyz);
      xyz=sp.rotate_x2();
      Residue::add_atom("O4'",xyz);
    }
  }

  string S_atoms[]={"C4'","O4'","C3'","C2'","C1'"};    //check sugar configuration
  vector<double> Sxyz,Sxyz0;
  for(auto atom : S_atoms) {
    vector<double> xyz=atomxyz[atom];
    Sxyz.insert(Sxyz.end(),xyz.begin(),xyz.end());
    xyz=BB_XYZ[atom];
    Sxyz0.insert(Sxyz0.end(),xyz.begin(),xyz.end());
  }
  Superpose sp(Sxyz,Sxyz0);
  double rmsd=sp.get_RMSD();
  if(rmsd>Backbone_RMSD_CUT) {
    string str="residue.cpp: Replace sugar ring for residue "+res_name+" --> RMSD=";
    str+=to_string(rmsd)+" Angstrom to standard C3'-endo configuration.";
    error->warning(str);

    vector<double> Vxyz=sp.rotate_x2();
    for(int ii=0;ii<get_array_size(S_atoms);ii++) {
      vector<double> xyz;
      for(int jj=0;jj<3;jj++)
        xyz.push_back(Vxyz[ii*3+jj]);
      string aname=S_atoms[ii];
      Residue::remove_atom(aname);
      atomxyz.erase(aname);
      Residue::add_atom(aname,xyz);
    }

    if(atomxyz.count("O2'")) {
      vector<double> xyz=BB_XYZ["O2'"];
      sp.reset_x2(xyz);
      xyz=sp.rotate_x2();
      Residue::remove_atom("O2'");
      atomxyz.erase("O2'");
      Residue::add_atom("O2'",xyz);
    }
  }

  if(!Residue::check_atom("O2'")) { //add atoms O2'
    vector<double> Vxyz,Vxyz0;
    vector<double> xyz;
    string ATOMR[]={"C3'","C2'","C1'"};
    for(auto atom : ATOMR) {
      if(atomxyz.count(atom)) {
        xyz=atomxyz[atom];
        Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
        xyz=BB_XYZ[atom];
        Vxyz0.insert(Vxyz0.end(),xyz.begin(),xyz.end());
      }
    }
    Superpose sp(Vxyz,Vxyz0);
    double rmsd=sp.get_RMSD();
    if(rmsd>Backbone_RMSD_CUT) {
      string str="residue.cpp: RMSD="+to_string(rmsd)+" Angstrom when recovery atom O2' for residue "+res_name;
      error->warning(str);
    }
    xyz=BB_XYZ["O2'"];
    sp.reset_x2(xyz);
    xyz=sp.rotate_x2();
    Residue::add_atom("O2'",xyz);
  }

  Residue::add_simplified_RNA_base_atoms();
  string B_atoms[]={"B1","B2","B3"}; //add atom B1/B2/B3 in base
  int nb=0;
  for(auto atom : B_atoms)
    if(Residue::check_atom(atom)) nb++;
  if(nb==2) {
    vector<double> Vxyz,Vxyz0;
    vector<double> xyz;
    string ATOMR[]={"C1'","B1","B2","B3"};
    for(auto atom : ATOMR) {
      if(atomxyz.count(atom)) {
        xyz=atomxyz[atom];
        Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
        xyz=BB_XYZ[atom];
        Vxyz0.insert(Vxyz0.end(),xyz.begin(),xyz.end());
      }
    }
    Superpose sp(Vxyz,Vxyz0);
    double rmsd=sp.get_RMSD();
    if(rmsd>Backbone_RMSD_CUT) {
      string str="residue.cpp: RMSD="+to_string(rmsd)+" Angstrom when recovery atom B1/B2/B3 for residue "+res_name;
      error->warning(str);
    }
    for(auto atom : B_atoms) {
      if(atomxyz.count(atom)==0) {
        xyz=BB_XYZ[atom];
        sp.reset_x2(xyz);
        xyz=sp.rotate_x2();
        Residue::add_atom(atom,xyz);
      }
    }
  }
  else if(nb<=1) {
    vector<double> Vxyz,Vxyz0;
    vector<double> xyz;
    string ATOMR[]={"C1'","C2'","O4'","B1","B2","B3"};
    for(auto atom : ATOMR) {
      if(atomxyz.count(atom)) {
        xyz=atomxyz[atom];
        Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
        xyz=BB_XYZ[atom];
        Vxyz0.insert(Vxyz0.end(),xyz.begin(),xyz.end());
      }
    }
    Superpose sp(Vxyz,Vxyz0);
    double rmsd=sp.get_RMSD();
    if(rmsd>Backbone_RMSD_CUT) {
      string str="residue.cpp: RMSD="+to_string(rmsd)+" Angstrom when recovery atom B1/B2/B3 for residue "+res_name;
      error->warning(str);
    }
    vector<double> Bxyz0;
    int Blen=get_array_size(B_atoms);
    for(int jj=0;jj<Blen;jj++) {
      xyz=BB_XYZ[B_atoms[jj]];
      Bxyz0.insert(Bxyz0.end(),xyz.begin(),xyz.end());
    }
    sp.reset_x2(Bxyz0);
    Bxyz0=sp.rotate_x2();
    for(int jj=0;jj<Blen;jj++) {
      string aname=B_atoms[jj];
      if(atomxyz.count(aname)) {
        atomxyz.erase(aname);
        Residue::remove_atom(aname);
      }
      xyz.clear();
      for(int kk=0;kk<3;kk++)
        xyz.push_back(Bxyz0[jj*3+kk]);
      Residue::add_atom(aname,xyz);
    }
  }
}

void Residue::simplify_RNA_base() {
  if(atoms.size()==0) return;
  if(atomxyz.count("B1")==0 || atomxyz.count("B2")==0 || atomxyz.count("B3")==0)
    Residue::add_simplified_RNA_base_atoms();
  if(atomxyz.count("B1") && !check_atom("B1")) 
    Residue::add_atom("B1",atomxyz["B1"]);
  if(atomxyz.count("B2") && !check_atom("B2"))
    Residue::add_atom("B2",atomxyz["B2"]); 
  if(atomxyz.count("B3") && !check_atom("B3"))
    Residue::add_atom("B3",atomxyz["B3"]);   

  vector<Atom> Tatoms=atoms;
  atoms.clear();
  for(auto it : Tatoms) {
    if(count(init->RESIDUES_ATOMS["rBB"].begin(),init->RESIDUES_ATOMS["rBB"].end(),it.atomname))
      atoms.push_back(it);
  }

  Residue::set_RNA_atomid(0,"std");
}

void Residue::recovery_RNA_base() {
  static double Base_RMSD_CUT=0.2;

  if(resname=="R" || resname=="Y" || resname=="N") {
    Residue::simplify_RNA_base();
    return;
  }

  map<string,vector<double> > base_x;
  vector<string> Vatoms;
  if(resname=="A" || resname=="C" || resname=="G" || resname=="U") {
    string basename="base"+resname;
    for(auto atom : init->RESIDUES_XYZ[basename]) {
      string name=atom.atomname;
      vector<double> xyz=atom.get_atom_xyz();
      base_x.emplace(name,xyz);
    }
    Vatoms=init->RESIDUES_ATOMS[basename];
  }

  Residue::remove_atom("B1");
  Residue::remove_atom("B2");
  Residue::remove_atom("B3");
  if(Residue::check_atoms(Vatoms)) {
    vector<double> Vxyz_tar,Vxyz_ref;
    for(size_t ii=0; ii<Vatoms.size();ii++) {
      vector<double> xyz=find_atom(Vatoms[ii]).get_atom_xyz();
      Vxyz_tar.insert(Vxyz_tar.end(),xyz.begin(),xyz.end());
      xyz=base_x[Vatoms[ii]];
      Vxyz_ref.insert(Vxyz_ref.end(),xyz.begin(),xyz.end());
    }
    Superpose sp(Vxyz_tar,Vxyz_ref);
    double rmsd=sp.get_RMSD();
    if(rmsd>Base_RMSD_CUT) {
      string str="residue.cpp: Replace base atoms for residue "+get_residue_name()+" --> RMSD=";
      str+=to_string(rmsd)+" Angstrom to standard configuration.";
      error->warning(str);
    }
    else return;
  } 

  Residue::simplify_RNA_base();
  Residue::remove_atom("B1");
  Residue::remove_atom("B2");
  Residue::remove_atom("B3");

  vector<double> Vxyz_B;
  vector<double> xyz=atomxyz["B1"];
  Vxyz_B.insert(Vxyz_B.end(),xyz.begin(),xyz.end());
  xyz=atomxyz["B2"];
  Vxyz_B.insert(Vxyz_B.end(),xyz.begin(),xyz.end());
  xyz=atomxyz["B3"];
  Vxyz_B.insert(Vxyz_B.end(),xyz.begin(),xyz.end());

  vector<double> Vxyz_B0;
  if(resname=="A" || resname=="G") {
    xyz=base_x["N9"];
    Vxyz_B0.insert(Vxyz_B0.end(),xyz.begin(),xyz.end());
    xyz=base_x["C8"];
    Vxyz_B0.insert(Vxyz_B0.end(),xyz.begin(),xyz.end());
    xyz=base_x["C4"];
    Vxyz_B0.insert(Vxyz_B0.end(),xyz.begin(),xyz.end());
  }
  else if(resname=="C" || resname=="U") {
    xyz=base_x["N1"];
    Vxyz_B0.insert(Vxyz_B0.end(),xyz.begin(),xyz.end());
    xyz=base_x["C6"];
    Vxyz_B0.insert(Vxyz_B0.end(),xyz.begin(),xyz.end());
    xyz=base_x["C2"];
    Vxyz_B0.insert(Vxyz_B0.end(),xyz.begin(),xyz.end());
  }
  vector<double> Vxyz0;
  for(auto& it : base_x) {
    xyz=it.second;
    Vxyz0.insert(Vxyz0.end(),xyz.begin(),xyz.end());
  }

  Superpose sp(Vxyz_B,Vxyz_B0);
  double rmsd=sp.get_RMSD();
  if(rmsd>0.5) {
    string str="residue.cpp: RMSD="+to_string(rmsd)+" Angstrom when recovery base atoms for residue "+get_residue_name();
    error->warning(str);
  }
  sp.reset_x2(Vxyz0);
  vector<double> Vxyz=sp.rotate_x2();

  int ni=0;
  for(auto& it : base_x) {
    xyz.clear();
    for(int jj=0;jj<3;jj++)
      xyz.push_back(Vxyz[ni*3+jj]);
    Residue::add_atom(it.first,xyz);
    ni++;
  }
}

vector<double> Residue::get_residue_xyz() const {
  vector<double> Vxyz;
  for(auto atom : atoms) {
    Vxyz.push_back(atom.x);
    Vxyz.push_back(atom.y);
    Vxyz.push_back(atom.z);
  }

  return Vxyz;
}

vector<double> Residue::get_RNA_residue_xyz() {
  Residue::set_RNA_atomid(0,"std");
  vector<double> Vxyz;
  for(auto atom : atoms) {
    Vxyz.push_back(atom.x);
    Vxyz.push_back(atom.y);
    Vxyz.push_back(atom.z);
  }

  return Vxyz;
}

void Residue::reset_residue_xyz(const vector<double> &Vxyz) {
  if(Vxyz.size()!=3*atoms.size()) {
    string str="residue.cpp: Number of atoms is not matched in setting of xyz --> ";
    str+=to_string(atoms.size())+" "+to_string(Vxyz.size()/3);
    error->one(str);
  }

  for(size_t ii=0;ii<atoms.size();ii++) {
    vector<double> xyz;
    for(int jj=0;jj<3;jj++)
      xyz.push_back(Vxyz[ii*3+jj]);
    atoms[ii].reset_atom_xyz(xyz);
    if(atomxyz.count(atoms[ii].atomname))
      atomxyz[atoms[ii].atomname]=xyz;
    else
      atomxyz.emplace(atoms[ii].atomname,xyz);
  }
}

void Residue::reset_residue_xyz(const vector<Atom> &Vatoms_one) {
  map<string,vector<double> >Matomname_xyz_one;
  for(size_t ii=0;ii<Vatoms_one.size();ii++) {
    string atomname_one=Vatoms_one[ii].atomname;
    vector<double> xyz=Vatoms_one[ii].get_atom_xyz();
    Matomname_xyz_one.emplace(atomname_one,xyz);
  } 

  for(size_t ii=0;ii<atoms.size();ii++) {
    string atomname_one=atoms[ii].atomname;
    if(Matomname_xyz_one.count(atomname_one)) {
      vector<double> &xyz=Matomname_xyz_one[atomname_one];
      atoms[ii].reset_atom_xyz(xyz);
      if(atomxyz.count(atomname_one))
        atomxyz[atomname_one]=xyz;
      else
        atomxyz.emplace(atomname_one,xyz);
    }
    else {
      string str="residue.cpp: Cannot find xyz for atom "+atoms[ii].get_atom_name();
      error->warning(str);
    }
  }
}

void Residue::set_residue_xyz(const Residue &res_one) {
  map<string,vector<double> > atomname_xyz;
  for(size_t ii=0;ii<res_one.atoms.size();ii++) {
    string atomname_one=res_one.atoms[ii].atomname;
    vector<double> xyz_one=res_one.atoms[ii].get_atom_xyz();
    if(atomname_xyz.count(atomname_one)==0) {
      atomname_xyz.emplace(atomname_one,xyz_one);
    }
  }

  for(size_t ii=0;ii<atoms.size();ii++) {
    string atomname_one=atoms[ii].atomname;
    if(atomname_xyz.count(atomname_one)) {
      atoms[ii].reset_atom_xyz(atomname_xyz[atomname_one]);
      if(atomxyz.count(atomname_one))
        atomxyz[atomname_one]=atomname_xyz[atomname_one];
      else
        atomxyz.emplace(atomname_one,atomname_xyz[atomname_one]);
    }
    else {
      string str="residue.cpp: Cannot find xyz for atom "+atoms[ii].atomname+" in residue "+get_residue_name();
      error->one(str);
    }
  }
}

void Residue::set_RNA_residue_xyz(Residue res_one) {
  atoms.clear();
  atomxyz.clear();
  if(res_one.resname==resname) {
    for(size_t ii=0;ii<res_one.atoms.size();ii++) {
      string atomname_one=res_one.atoms[ii].atomname;
      vector<double> xyz_one=res_one.atoms[ii].get_atom_xyz();
      Residue::add_atom(atomname_one,xyz_one);
    }
  }
  else {
    res_one.simplify_RNA_base();
    for(size_t ii=0;ii<res_one.atoms.size();ii++) {
      string atomname_one=res_one.atoms[ii].atomname;
      vector<double> xyz_one=res_one.atoms[ii].get_atom_xyz();
      Residue::add_atom(atomname_one,xyz_one);
    }
  }

  Residue::recovery_RNA_base();
  Residue::set_RNA_atomid(0,"std");
}

bool Residue::check_atom(string name_one) const {
  int check=0;
  for(size_t ii=0;ii<atoms.size();ii++) {
    if(atoms[ii].atomname==name_one) {
      check=1;
      break;
    }
  }
  
  if(check==1) return true;
  else return false;
}

bool Residue::is_atom_exist(string name_one) const {
  int check=0;
  for(size_t ii=0;ii<atoms.size();ii++) {
    if(atoms[ii].atomname==name_one) {
      check=1;
      break;
    }
  }
  
  if(check==1) return true;
  else return false;
}

Atom Residue::find_atom(string name_one) const {
  int check=0;
  int ni=-1;
  for(size_t ii=0;ii<atoms.size();ii++) {
    if(atoms[ii].atomname==name_one) {
      ni=ii;
      check=1;
      break;
    }
  }
  if(!check) {
    string str="residue.cpp: Can not find atom "+name_one+" in residue "+get_residue_name();
    error->one(str);
  }

  return atoms[ni];
}

Atom & Residue::operator[](string name_one) {
  int check=0;
  int ni=-1;
  for(size_t ii=0;ii<atoms.size();ii++) {
    if(atoms[ii].atomname==name_one) {
      ni=ii;
      check=1;
      break;
    }
  }
  if(!check) {
    string str="residue.cpp: Can not find atom "+name_one+" in residue "+get_residue_name();
    error->one(str);
  }

  return atoms[ni];
}

void Residue::delete_atom(string name_one) {
  vector<Atom> Tatoms=atoms;
  atoms.clear();
  for(size_t ii=0;ii<Tatoms.size();ii++) {
    if(Tatoms[ii].atomname==name_one) continue;
    atoms.push_back(Tatoms[ii]);
  }
  if(atomxyz.count(name_one))
    atomxyz.erase(name_one);
}

void Residue::delete_atom(vector<string> &Vname_one) {
  vector<Atom> Tatoms=atoms;
  atoms.clear();
  for(size_t ii=0;ii<Tatoms.size();ii++) {
    if(count(Vname_one.begin(),Vname_one.end(),Tatoms[ii].atomname)) continue;
    atoms.push_back(Tatoms[ii]);
  }
  for(size_t ii=0;ii<Vname_one.size();ii++) {
    if(atomxyz.count(Vname_one[ii]))
      atomxyz.erase(Vname_one[ii]);
  }
}

void Residue::delete_atom_hydrogen() {
  vector<Atom> Tatoms=atoms;
  atoms.clear();
  atomxyz.clear();
  for(size_t ii=0;ii<Tatoms.size();ii++) {
    string atomname_one=Tatoms[ii].atomname;
    if(atomname_one[0]=='H' || (atomname_one.size()>1 && isdigit(atomname_one[0]) && atomname_one[1]=='H')) continue;
    atoms.push_back(Tatoms[ii]);
    vector<double> xyz=Tatoms[ii].get_atom_xyz();
    atomxyz.emplace(atomname_one,xyz);
  }
}

size_t Residue::size() const {
  return atoms.size();
}

vector<double> Residue::center_of_geometry() const {
  vector<double> Vxyz;
  for(size_t ii=0;ii<3;ii++)
    Vxyz.push_back(0.0);

  size_t natoms=atoms.size();
  for(size_t ii=0;ii<natoms;ii++) {
    Vxyz[0]+=atoms[ii].x;
    Vxyz[1]+=atoms[ii].y;
    Vxyz[2]+=atoms[ii].z;
  }
  if(natoms>0) {
    Vxyz[0]/=natoms;
    Vxyz[1]/=natoms;
    Vxyz[2]/=natoms;
  }

  return Vxyz;
}

void Residue::clear() {
  chainid="";
  resid=-1;
  resname="";
  atoms.clear();
  atoms.shrink_to_fit();
  atomxyz.clear();
}

int Residue::set_atomid(int index_one) {
  int index_after=index_one+atoms.size();
  if(index_after>ATOMID_MAX) index_after-=ATOMID_MAX;
  int nid=0;
  vector<Atom> Tatoms=atoms;
  atoms.clear();

  if(init->RESIDUES_ATOMS.count(resname)==0) {
    string str="residue.cpp: Undefined atoms list for residue "+resname;
    error->one(str);
  }

  for(auto str : init->RESIDUES_ATOMS[resname]) {
    for(size_t ii=0;ii<Tatoms.size();ii++) {
      if(Tatoms[ii].atomname==str) {
        nid++;
        int index=nid+index_one;
        if(index>ATOMID_MAX) index-=ATOMID_MAX;
        Tatoms[ii].id=index;
        atoms.push_back(Tatoms[ii]);
        Tatoms.erase(Tatoms.begin()+ii);
        break;
      }
    }
  }

  for(size_t ii=0;ii<Tatoms.size();ii++) {
    nid++;
    int index=nid+index_one;
    if(index>ATOMID_MAX) index-=ATOMID_MAX;
    Tatoms[ii].id=index;
    atoms.push_back(Tatoms[ii]);
  }

  return index_after;
}

int Residue::set_RNA_atomid(int index_one, string format) {
  int index_after=index_one+atoms.size();
  if(index_after>ATOMID_MAX) index_after-=ATOMID_MAX;
  int nid=0;
  string resname_out=resname;
  vector<Atom> Tatoms=atoms;
  atoms.clear();
  for(auto str : init->RESIDUES_ATOMS["rBB"]) {
    for(size_t ii=0;ii<Tatoms.size();ii++) {
      if(Tatoms[ii].atomname==str) {
        nid++;
        int index=nid+index_one;
        if(index>ATOMID_MAX) index-=ATOMID_MAX;
        Tatoms[ii].id=index;
        atoms.push_back(Tatoms[ii]);
        Tatoms.erase(Tatoms.begin()+ii);
        break;
      }
    }
  }

  if(resname=="A") {
    for(auto str : init->RESIDUES_ATOMS["baseA"]) {
      for(size_t ii=0;ii<Tatoms.size();ii++) {
        if(Tatoms[ii].atomname==str) {
          nid++;
          int index=nid+index_one;
          if(index>ATOMID_MAX) index-=ATOMID_MAX;
          Tatoms[ii].id=index;
          atoms.push_back(Tatoms[ii]);
          Tatoms.erase(Tatoms.begin()+ii);
          break;
        }
      }
    }
    if(format=="std") resname_out="A";
    else if(format=="gmx") resname_out="ADE";
    else if(format=="2") resname_out="RA";
  }
  else if(resname=="G") {
    for(auto str : init->RESIDUES_ATOMS["baseG"]) {
      for(size_t ii=0;ii<Tatoms.size();ii++) {
        if(Tatoms[ii].atomname==str) {
          nid++;
          int index=nid+index_one;
          if(index>ATOMID_MAX) index-=ATOMID_MAX;
          Tatoms[ii].id=index;
          atoms.push_back(Tatoms[ii]);
          Tatoms.erase(Tatoms.begin()+ii);
          break;
        }
      }
    }
    if(format=="std") resname_out="G";
    else if(format=="gmx") resname_out="GUA";
    else if(format=="2") resname_out="RG";
  }
  else if(resname=="C") {
    for(auto str : init->RESIDUES_ATOMS["baseC"]) {
      for(size_t ii=0;ii<Tatoms.size();ii++) {
        if(Tatoms[ii].atomname==str) {
          nid++;
          int index=nid+index_one;
          if(index>ATOMID_MAX) index-=ATOMID_MAX;
          Tatoms[ii].id=index;
          atoms.push_back(Tatoms[ii]);
          Tatoms.erase(Tatoms.begin()+ii);
          break;
        }
      }
    }
    if(format=="std") resname_out="C";
    else if(format=="gmx") resname_out="CYT";
    else if(format=="2") resname_out="RC";
  }
  else if(resname=="U") {
    for(auto str : init->RESIDUES_ATOMS["baseU"]) {
      for(size_t ii=0;ii<Tatoms.size();ii++) {
        if(Tatoms[ii].atomname==str) {
          nid++;
          int index=nid+index_one;
          if(index>ATOMID_MAX) index-=ATOMID_MAX;
          Tatoms[ii].id=index;
          atoms.push_back(Tatoms[ii]);
          Tatoms.erase(Tatoms.begin()+ii);
          break;
        }
      }
    }
    if(format=="std") resname_out="U";
    else if(format=="gmx") resname_out="URA";
    else if(format=="2") resname_out="RU";
  }
  else {
    if(init->RESIDUES_ATOMS.count(resname)==0) {
      string str="residue.cpp: Undefined atoms list for residue "+resname;
      error->one(str);
    }

    for(auto str : init->RESIDUES_ATOMS[resname]) {
      for(size_t ii=0;ii<Tatoms.size();ii++) {
        if(Tatoms[ii].atomname==str) {
          nid++;
          int index=nid+index_one;
          if(index>ATOMID_MAX) index-=ATOMID_MAX;
          Tatoms[ii].id=index;
          atoms.push_back(Tatoms[ii]);
          Tatoms.erase(Tatoms.begin()+ii);
          break;
        }
      }
    }
    resname_out=resname;
  }

  for(size_t ii=0;ii<Tatoms.size();ii++) {
    nid++;
    int index=nid+index_one;
    if(index>ATOMID_MAX) index-=ATOMID_MAX;
    Tatoms[ii].id=index;
    atoms.push_back(Tatoms[ii]);
  }

  resname=resname_out;
  for(size_t ii=0;ii<atoms.size();ii++) {
    atoms[ii].resname=resname_out;
    if(atoms[ii].atomname=="OP1" && format=="gmx")
      atoms[ii].atomname="O1P";
    if(atoms[ii].atomname=="OP2" && format=="gmx")
      atoms[ii].atomname="O2P";
  }

  return index_after;
}

int Residue::set_atomid_empty(int index_one) {
  int index_after=index_one+atoms.size();
  if(index_after>ATOMID_MAX) index_after-=ATOMID_MAX;
  int nid=0;
  vector<Atom> Tatoms=atoms;
  atoms.clear();

  for(size_t ii=0;ii<Tatoms.size();ii++) {
    nid++;
    int index=nid+index_one;
    if(index>ATOMID_MAX) index-=ATOMID_MAX;
    Tatoms[ii].id=index;
    atoms.push_back(Tatoms[ii]);
  }

  return index_after;
}

void Residue::add_simplified_RNA_base_atoms() {
  if(resname=="A" || resname=="G" || resname=="R") {  //N9=B1,C8=B2,C4=B3
    vector<double> xyz;
    if(!Residue::check_atom("B1") && atomxyz.count("N9")) {
      xyz=atomxyz["N9"];
      Residue::add_atom("B1",xyz);
    }
    else {
      string str="residue.cpp: Failed to add simplified base atom B1";
      //error->warning(str);
    }
    if(!Residue::check_atom("B2") && atomxyz.count("C8")) {
      xyz=atomxyz["C8"];
      Residue::add_atom("B2",xyz);
    }
    else {
      string str="residue.cpp: Failed to add simplified base atom B2";
      //error->warning(str);
    }
    if(!Residue::check_atom("B3") && atomxyz.count("C4")) {
      xyz=atomxyz["C4"];
      Residue::add_atom("B3",xyz);
    }
    else {
      string str="residue.cpp: Failed to add simplified base atom B3";
      //error->warning(str);
    }
  }
  else if(resname=="C" || resname=="U" || resname=="Y") {//N1=B1,C6=B2,C2=B3
    vector<double> xyz;
    if(!Residue::check_atom("B1") && atomxyz.count("N1")) {
      xyz=atomxyz["N1"];
      Residue::add_atom("B1",xyz);
    }
    else {
      string str="residue.cpp: Failed to add simplified base atom B1";
      //error->warning(str);
    }
    if(!Residue::check_atom("B2") && atomxyz.count("C6")) {
      xyz=atomxyz["C6"];
      Residue::add_atom("B2",xyz);
    }
    else {
      string str="residue.cpp: Failed to add simplified base atom B2";
      //error->warning(str);
    }
    if(!Residue::check_atom("B3") && atomxyz.count("C2")) {
      xyz=atomxyz["C2"];
      Residue::add_atom("B3",xyz);
    }
    else {
      string str="residue.cpp: Failed to add simplified base atom B3";
      //error->warning(str);
    }
  }
  else {
    string str="residue.cpp: Failed to add simplified base atom B1/B2/B3";
    //error->warning(str);
  }
}

inline bool Residue::check_atoms(const vector<string> &Vname) {
  int check=0;
  for(size_t ii=0;ii<Vname.size();ii++) {
    string aname=Vname[ii];
    for(size_t jj=0;jj<atoms.size();jj++) {
      if(atoms[jj].atomname==aname) {
        check++;
        break;
      }
    }
  }

  if(check==((int)Vname.size())) return 1;
  else return 0;
}

inline bool Residue::add_atom(string name_one, const vector<double> &Vxyz) {
  if(Residue::check_atom(name_one)) {
    string str="residue.cpp: Failed to add atom --> atom "+name_one+" has been existed in residue "+get_residue_name();
    error->warning(str);
    return 0;
  }

  int nid=atoms.size()+1;
  Atom atom_one;
  double x=Vxyz[0];
  double y=Vxyz[1];
  double z=Vxyz[2];
  atom_one.set_atom(nid,chainid,resid,resname,name_one,x,y,z);
  atoms.push_back(atom_one);
  if(atomxyz.count(name_one)==0) {
    vector<double> xyz=Vxyz;
    atomxyz.insert(pair<string,vector<double> >(name_one,xyz));
  }

  return 1;
}

inline void Residue::remove_atom(string name_one) {
  vector<Atom> Tatoms=atoms;
  atoms.clear();
  for(size_t ii=0;ii<Tatoms.size();ii++) {
    if(Tatoms[ii].atomname==name_one) continue;
    atoms.push_back(Tatoms[ii]);
  }
}

ostream & operator <<(ostream &output, Residue &residue_one) {
  for(size_t ii=0;ii<residue_one.atoms.size();ii++) {
    output<<residue_one.atoms[ii];
  }

  return output;
}