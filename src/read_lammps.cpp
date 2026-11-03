#include "basic.h"
#include "constant.h"
#include "readfile.h"
#include "initialization.h"
#include "error.h"
#include "sort_vector.h"
#include "read_PDB.h"
#include "read_lammps.h"

using namespace CBS_NS;

ReadLAMMPS::ReadLAMMPS(Cobmms *cbs) : Pointers(cbs) {
  step_begin=0;
  step_end=-1;
}

ReadLAMMPS::ReadLAMMPS(const ReadLAMMPS &other) : Pointers(other.cbs) {
  step_begin=other.step_begin;
  step_end=other.step_end;
  VVatoms=other.VVatoms;
  VVsteps=other.VVsteps;
  VVenergies=other.VVenergies;
  atomlabel=other.atomlabel;
  Vbonds=other.Vbonds;
  Vangles=other.Vangles;
  Vdihedrals=other.Vdihedrals;
  Vimpropers=other.Vimpropers;
}

ReadLAMMPS& ReadLAMMPS::operator=(const ReadLAMMPS &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    step_begin=other.step_begin;
    step_end=other.step_end;
    VVatoms=other.VVatoms;
    VVsteps=other.VVsteps;
    VVenergies=other.VVenergies;
    atomlabel=atomlabel;
    Vbonds=other.Vbonds;
    Vangles=other.Vangles;
    Vdihedrals=other.Vdihedrals;
    Vimpropers=other.Vimpropers;
  }

  return *this;
}

ReadLAMMPS::~ReadLAMMPS() {
  step_begin=0;
  step_end=-1;
  VVatoms.clear();
  VVatoms.shrink_to_fit();
  VVsteps.clear();
  VVsteps.shrink_to_fit();
  VVenergies.clear();
  VVenergies.shrink_to_fit();
  atomlabel.clear();
  Vbonds.clear();
  Vbonds.shrink_to_fit();
  Vangles.clear();
  Vangles.shrink_to_fit();
  Vdihedrals.clear();
  Vdihedrals.shrink_to_fit();
  Vimpropers.clear();
  Vimpropers.shrink_to_fit();
}

void ReadLAMMPS::read_log(string filename, int step_begin_one, int step_end_one) {
  if(!is_file_exist(filename)) {
    string str="read_lammps.cpp: Cannot find file "+filename;
    error->one(str);
  }

  step_begin=step_begin_one;
  step_end=step_end_one;

  VVenergies.clear();
  if(VVsteps.size()==0 || VVsteps[0].size()==0) {
    string str="read_lammps.cpp: Please read .lammpstrj file(s) first!";
    error->one(str);
  }
  
  vector<string> Venergy;
  int flag=0;
  vector<int> Vstep;
  vector<string> Vlines=read_file(filename);
  for(size_t ii=0;ii<Vlines.size();ii++) {
    string sline=Vlines[ii];
    if(sline.find("Step")!=string::npos && sline.find("E_pair")!=string::npos && sline.find("E_mol")!=string::npos) {
      flag=1; 
      continue;
    }
    else if(sline.find("Loop")!=string::npos && sline.find("time")!=string::npos && sline.find("procs")!=string::npos) {
      flag=0; 
      continue;
    }
    if(flag==1) {
      vector<string> stoken=split_string(sline);
      if(stoken.size()<6 || stoken.size()>7) continue;
      if(stoken[0][0]<'0' || stoken[0][0]>'9') continue;
      int step=stoi(stoken[0]);
      if(count(Vstep.begin(),Vstep.end(),step)==0)
        Vstep.push_back(step);
      else
        continue;
      if(step<step_begin || (step_end>0 && step>step_end)) continue;
      if(count(VVsteps[0].begin(),VVsteps[0].end(),step)==0) continue;
      for(size_t mm=0;mm<stoken[2].size();mm++)
        stoken[2][mm]=tolower(stoken[2][mm]);
      if(stoken[2].find("nan")!=string::npos) continue;
      for(size_t mm=0;mm<stoken[3].size();mm++)
        stoken[3][mm]=tolower(stoken[3][mm]);
      if(stoken[3].find("nan")!=string::npos) continue;
      double E_pair=stod(stoken[2]);
      double E_mol=stod(stoken[3]);
      double E_pot=E_pair+E_mol;
      string str=to_string(step)+" "+to_string(E_pot);
      Venergy.push_back(str);
    }
  }
  VVenergies.push_back(Venergy);
}

void ReadLAMMPS::read_log(const vector<string> &Vfilename, int step_begin_one, int step_end_one) {
  step_begin=step_begin_one;
  step_end=step_end_one;

  VVenergies.clear();
  if(VVsteps.size()==0 || VVsteps[0].size()==0) {
    string str="read_lammps.cpp: Please read .lammpstrj file(s) first!";
    error->one(str);
  }

  for(size_t ii=0;ii<Vfilename.size();ii++) {
    if(Vfilename[ii].size()<2) continue;
    if(!is_file_exist(Vfilename[ii])) {
      string str="read_lammps.cpp: Cannot find file "+Vfilename[ii];
      error->one(str);
    }

    vector<string> Venergy;
    int flag=0;
    vector<int> Vstep;
    vector<string> Vlines=read_file(Vfilename[ii]);
    for(size_t jj=0;jj<Vlines.size();jj++) {
      string sline=Vlines[jj];
      if(sline.find("Step")!=string::npos && sline.find("E_pair")!=string::npos && sline.find("E_mol")!=string::npos) {
        flag=1; 
        continue;
      }
      else if(sline.find("Loop")!=string::npos && sline.find("time")!=string::npos && sline.find("procs")!=string::npos) {
        flag=0; 
        continue;
      }
      if(flag==1) {
        vector<string> stoken=split_string(sline);
        if(stoken.size()<6 || stoken.size()>7) continue;
        if(stoken[0][0]<'0' || stoken[0][0]>'9') continue;
        int step=stoi(stoken[0]);
        if(count(Vstep.begin(),Vstep.end(),step)==0)
          Vstep.push_back(step);
        else
          continue;
        if(step<step_begin || (step_end>0 && step>step_end)) continue;
        if(count(VVsteps[ii].begin(),VVsteps[ii].end(),step)==0) continue;
        for(size_t mm=0;mm<stoken[2].size();mm++)
          stoken[2][mm]=tolower(stoken[2][mm]);
        if(stoken[2].find("nan")!=string::npos) continue;
        for(size_t mm=0;mm<stoken[3].size();mm++)
          stoken[3][mm]=tolower(stoken[3][mm]);
        if(stoken[3].find("nan")!=string::npos) continue;
        double E_pair=stod(stoken[2]);
        double E_mol=stod(stoken[3]);
        double E_pot=E_pair+E_mol;
        string str=to_string(step)+" "+to_string(E_pot);
        Venergy.push_back(str);
      }
    }
    VVenergies.push_back(Venergy);
  }
}

void ReadLAMMPS::read_lammpstrj(string filename, int step_begin_one, int step_end_one) {
  if(!is_file_exist(filename)) {
    string str="read_lammps.cpp: Cannot find file "+filename;
    error->one(str);
  }
  
  step_begin=step_begin_one;
  step_end=step_end_one;

  VVatoms.clear();
  VVsteps.clear();

  vector<vector<Atom> > Vatoms;
  vector<int> Vstep;
  vector<string> Vlines=read_file(filename);
  size_t ii=0;
  while(ii<Vlines.size()) {
    if(Vlines[ii].size()<1) continue;
    if(Vlines[ii].find("TIMESTEP")!=string::npos) {
      ii++;
      vector<string> stoken=split_string(Vlines[ii]);
      int step=stoi(stoken[0]);
      ii=ii+2;
      stoken=split_string(Vlines[ii]);
      int natom=stoi(stoken[0]);
      if(step<step_begin || (step_end>0 && step>step_end)) {
        ii=ii+natom+6;
        continue;
      }
      ii=ii+5;
      vector<string> Vxyz;
      Vxyz.resize(natom);
      for(int jj=0;jj<natom;jj++) {
        ii++;
        stoken=split_string(Vlines[ii]);  //id mol type xu yu zu
        int id_one=stoi(stoken[0]);
        if(id_one>natom) {
          string str="read_lammps.cpp: Atom id should be renamed (id > natom) -->  "+to_string(id_one)+"/"+to_string(natom);
          error->one(str);
        }
        Vxyz[id_one-1]=Vlines[ii];
      }
      
      vector<Atom> atoms;
      for(int jj=0;jj<natom;jj++) {
        stoken=split_string(Vxyz[jj]); 
        int id_one=stoi(stoken[0]);
        string chainid_one=stoken[1];
        string atomname_one=stoken[2];
        for(size_t mm=0;mm<stoken[3].size();mm++)
          stoken[3][mm]=tolower(stoken[3][mm]);
        if(stoken[3].find("nan")!=string::npos) continue;
        for(size_t mm=0;mm<stoken[4].size();mm++)
          stoken[4][mm]=tolower(stoken[4][mm]);
        if(stoken[4].find("nan")!=string::npos) continue;
        for(size_t mm=0;mm<stoken[5].size();mm++)
          stoken[5][mm]=tolower(stoken[5][mm]);
        if(stoken[5].find("nan")!=string::npos) continue;
        double x_one=stod(stoken[3]);
        double y_one=stod(stoken[4]);
        double z_one=stod(stoken[5]);
        int resid_one=1;
        string resname_one=" ";
        Atom atom_one;
        atom_one.set_atom(id_one,chainid_one,resid_one,resname_one,atomname_one,x_one,y_one,z_one);
        atoms.push_back(atom_one);
      }
      if((int)atoms.size()==natom) {
        Vstep.push_back(step);
        Vatoms.push_back(atoms);
      }
      else {
        string str="read_lammps.cpp: Skip step="+to_string(step)+" in "+filename;
        error->warning(str);
      }
    }
    else {
      string str="read_lammps.cpp: Unknown line in "+filename+" --> "+Vlines[ii];
      error->warning(str);
    }
    ii++;
  }
  VVsteps.push_back(Vstep);  
  VVatoms.push_back(Vatoms);
}
     
void ReadLAMMPS::read_lammpstrj(const vector<string> &Vfilename, int step_begin_one, int step_end_one) {
  step_begin=step_begin_one;
  step_end=step_end_one;

  VVatoms.clear();
  VVsteps.clear();

  for(size_t ii=0;ii<Vfilename.size();ii++) {
    if(Vfilename[ii].size()<2) continue;
    if(!is_file_exist(Vfilename[ii])) {
      string str="read_lammps.cpp: Cannot find file "+Vfilename[ii];
      error->one(str);
    }

    vector<vector<Atom> > Vatoms;
    vector<int> Vstep;
    vector<string> Vlines=read_file(Vfilename[ii]);
    size_t jj=0;
    while(jj<Vlines.size()) {
      if(Vlines[jj].size()<1) continue;
      if(Vlines[jj].find("TIMESTEP")!=string::npos) {
        jj++;
        vector<string> stoken=split_string(Vlines[jj]);
        int step=stoi(stoken[0]);
        jj=jj+2;
        stoken=split_string(Vlines[jj]);
        int natom=stoi(stoken[0]);
        if(step<step_begin || (step_end>0 && step>step_end)) {
          jj=jj+natom+6;
          continue;
        }
        jj=jj+5;
        vector<string> Vxyz;
        Vxyz.resize(natom);
        for(int kk=0;kk<natom;kk++) {
          jj++;
          stoken=split_string(Vlines[jj]);  //id mol type xu yu zu
          int id_one=stoi(stoken[0]);
          if(id_one>natom) {
            string str="read_lammps.cpp: Atom id should be renamed (id > natom) -->  "+to_string(id_one)+"/"+to_string(natom);
            error->one(str);
          }
          Vxyz[id_one-1]=Vlines[jj];
        }

        vector<Atom> atoms;
        for(int kk=0;kk<natom;kk++) {
          stoken=split_string(Vxyz[kk]); 
          int id_one=stoi(stoken[0]);
          string chainid_one=stoken[1];
          string atomname_one=stoken[2];
          for(size_t mm=0;mm<stoken[3].size();mm++)
            stoken[3][mm]=tolower(stoken[3][mm]);
          if(stoken[3].find("nan")!=string::npos) continue;
          for(size_t mm=0;mm<stoken[4].size();mm++)
            stoken[4][mm]=tolower(stoken[4][mm]);
          if(stoken[4].find("nan")!=string::npos) continue;
          for(size_t mm=0;mm<stoken[5].size();mm++)
            stoken[5][mm]=tolower(stoken[5][mm]);
          if(stoken[5].find("nan")!=string::npos) continue;
          double x_one=stod(stoken[3]);
          double y_one=stod(stoken[4]);
          double z_one=stod(stoken[5]);
          int resid_one=1;
          string resname_one=" ";
          Atom atom_one;
          atom_one.set_atom(id_one,chainid_one,resid_one,resname_one,atomname_one,x_one,y_one,z_one);
          atoms.push_back(atom_one);
        }
        if((int)atoms.size()==natom) {
          Vstep.push_back(step);
          Vatoms.push_back(atoms);
        }
        else {
          string str="read_lammps.cpp: Skip step="+to_string(step)+" in "+Vfilename[ii];
          error->warning(str);
        }
      }
      else {
        string str="read_lammps.cpp: Unknown line in "+Vfilename[ii]+" --> "+Vlines[jj];
        error->warning(str);
      }
      jj++;
    }
    VVsteps.push_back(Vstep);  
    VVatoms.push_back(Vatoms);
  }
}

void ReadLAMMPS::recover_residue_info(string filename) {
  if(!is_file_exist(filename)) {
    string str="read_lammps.cpp: Cannot find file "+filename;
    error->one(str);
  }

  vector<Atom> Vatoms_init;
  ReadPDB readpdb_one(cbs);
  readpdb_one.read_PDB(filename);
  Vatoms_init=readpdb_one.Full_atoms[0];
  /*
  vector<string> Vlines=read_file(filename);
  for(size_t ii=0;ii<Vlines.size();ii++) {
    if(Vlines[ii].size()<2) continue;
    if(Vlines[ii].find("#")!=string::npos) continue;
    if(Vlines[ii].find("TER")!=string::npos) continue;
    if(Vlines[ii].find("END")!=string::npos) break;
    Atom atom_one;
    atom_one.set_atom(Vlines[ii]);
    Vatoms_init.push_back(atom_one);
  }
  */

  for(size_t ii=0;ii<VVatoms.size();ii++) {
    for(size_t jj=0;jj<VVatoms[ii].size();jj++) {
      if(VVatoms[ii][jj].size()!=Vatoms_init.size()) {
        string str="read_lammps.cpp: Numbers of atoms are not matched --> "+to_string(VVatoms[ii][jj].size())+"/"+to_string(Vatoms_init.size());
        error->one(str);
      }
      for(size_t kk=0;kk<Vatoms_init.size();kk++) {
        VVatoms[ii][jj][kk].chainid=Vatoms_init[kk].chainid;
        VVatoms[ii][jj][kk].resid=Vatoms_init[kk].resid;
        VVatoms[ii][jj][kk].resname=Vatoms_init[kk].resname;
        VVatoms[ii][jj][kk].atomname=Vatoms_init[kk].atomname;
      }
    }
  }
}

void ReadLAMMPS::recover_residue_info(string filename, vector<Atom> &Vatoms_one) {
  if(!is_file_exist(filename)) {
    string str="read_lammps.cpp: Cannot find file "+filename;
    error->one(str);
  }

  vector<Atom> Vatoms_init;
  ReadPDB readpdb_one(cbs);
  readpdb_one.read_PDB(filename);
  Vatoms_init=readpdb_one.Full_atoms[0];
  /*
  vector<string> Vlines=read_file(filename);
  for(size_t ii=0;ii<Vlines.size();ii++) {
    if(Vlines[ii].size()<2) continue;
    if(Vlines[ii].find("#")!=string::npos) continue;
    if(Vlines[ii].find("TER")!=string::npos) continue;
    if(Vlines[ii].find("END")!=string::npos) break;
    Atom atom_one;
    atom_one.set_atom(Vlines[ii]);
    Vatoms_init.push_back(atom_one);
  }
  */

  if(Vatoms_one.size()!=Vatoms_init.size()) {
    string str="read_lammps.cpp: Numbers of atoms are not matched --> "+to_string(Vatoms_one.size())+"/"+to_string(Vatoms_init.size());
    error->one(str);
  }
  for(size_t kk=0;kk<Vatoms_init.size();kk++) {
    Vatoms_one[kk].chainid=Vatoms_init[kk].chainid;
    Vatoms_one[kk].resid=Vatoms_init[kk].resid;
    Vatoms_one[kk].resname=Vatoms_init[kk].resname;
    Vatoms_one[kk].atomname=Vatoms_init[kk].atomname;
  }
}

vector<vector<Atom> > ReadLAMMPS::extract_top_structures(double Pe, double Ps, int nsp) const {
  if(VVenergies.size()==0) {
    string str="read_lammps.cpp: Please read .log file(s) first!";
    error->one(str);
  }
  
  vector<vector<Atom> > Vatoms_top;
  vector<double> Venergy_top;
  for(size_t ii=0;ii<VVsteps.size();ii++) {
    map<int,size_t> Mstep_index;
    for(size_t jj=0;jj<VVsteps[ii].size();jj++) {
      Mstep_index.emplace(VVsteps[ii][jj],jj);
    }

    vector<double> Venergy;
    int nsnap=static_cast<int>(VVenergies[ii].size());
    int index_begin=round((1.0-Ps)*nsnap);
    if(index_begin<0) index_begin=0;
    else if(index_begin>nsnap-1) index_begin=nsnap-1;
    for(int jj=index_begin;jj<nsnap;jj++) {
      vector<string> stoken=split_string(VVenergies[ii][jj]);
      double energy_i=stod(stoken[1]);
      Venergy.push_back(energy_i);
    }
    if(Pe<0.0) Pe=0.0;
    else if(Pe>1.0) Pe=1.0;
    vector<int> Vselected;
    int num_one=round(Venergy.size()/nsp);
    for(int kk=0;kk<nsp;kk++) {
      vector<double> Venergy_one;
      int index_b=kk*num_one;
      int index_e=(kk+1)*num_one;
      if(kk==nsp-1) 
        index_e=static_cast<int>(Venergy.size());
      for(int nn=index_b;nn<index_e;nn++) {
        Venergy_one.push_back(Venergy[nn]);
      }
      vector<size_t> Vsort=SortVector::sort_vector_ascending(Venergy_one);
      int nstru=ceil(Venergy_one.size()*Pe);
      if(nstru<1) nstru=1;
      for(int jj=0;jj<nstru;jj++) {
        int index_i=Vsort[jj]+index_begin+index_b;
        Vselected.push_back(index_i);
      }
    }
    for(size_t jj=0;jj<Vselected.size();jj++) {
      vector<string> stoken=split_string(VVenergies[ii][Vselected[jj]]);
      int step_one=stoi(stoken[0]);
      double energy_one=stod(stoken[1]);
      size_t index_one=Mstep_index[step_one];
      Vatoms_top.push_back(VVatoms[ii][index_one]);
      Venergy_top.push_back(energy_one);
    }
  }

  return Vatoms_top;
}

vector<Atom> ReadLAMMPS::extract_lowest_energy_structure() const {
  if(VVenergies.size()==0) {
    string str="read_lammps.cpp: Please read .log file(s) first!";
    error->one(str);
  }

  vector<Atom> Vatoms_lowest;
  double energy_lowest=MathConst::BIG;
  for(size_t ii=0;ii<VVenergies.size();ii++) {
    int index_lowest_i=-1;
    double energy_lowest_i=MathConst::BIG;
    for(size_t jj=0;jj<VVenergies[ii].size();jj++) {
      vector<string> stoken=split_string(VVenergies[ii][jj]);
      double energy_one=stod(stoken[1]);
      if(index_lowest_i<0) {
        index_lowest_i=jj;
        energy_lowest_i=energy_one;
      }
      if(energy_one<energy_lowest_i) {
        index_lowest_i=jj;
        energy_lowest_i=energy_one;
      }
    }
    if(ii==0) {
      Vatoms_lowest=VVatoms[ii][index_lowest_i];
      energy_lowest=energy_lowest_i;
    }
    else {
      if(energy_lowest_i<energy_lowest) {
        Vatoms_lowest=VVatoms[ii][index_lowest_i];
        energy_lowest=energy_lowest_i;
      }
    }
  }

  return Vatoms_lowest;
}

void ReadLAMMPS::read_data(string filename) {
  if(!is_file_exist(filename)) {
    string str="read_lammps.cpp: Cannot find file "+filename;
    error->one(str);
  }

  atomlabel.clear();
  Vbonds.clear();
  Vangles.clear();
  Vdihedrals.clear();
  Vimpropers.clear();

  size_t nbond=0;
  size_t nangle=0;
  size_t ndihedral=0;
  size_t nimproper=0;
  size_t natomtype=0;
  vector<string> Vlines=read_file(filename);
  for(size_t ii=0;ii<14;ii++) {
    if(Vlines[ii].size()<2) continue;
    if(Vlines[ii].find("#")!=string::npos) continue;
    vector<string> stoken=split_string(Vlines[ii]);
    if(stoken.size()<2) continue;
    if(stoken[1]=="bonds")
      nbond=stoi(stoken[0]);
    else if(stoken[1]=="angles")
      nangle=stoi(stoken[0]);
    else if(stoken[1]=="dihedrals")
      ndihedral=stoi(stoken[0]);
    else if(stoken[1]=="impropers")
      nimproper=stoi(stoken[0]);
    else if(stoken.size()==3 && stoken[1]=="atom" && stoken[2]=="types")
      natomtype=stoi(stoken[0]);
  }
  
  int flag=-1;
  enum {fBond,fAngle,fDihedral,fImproper,fAtomTypeLabel,fOther};
  for(size_t ii=12;ii<Vlines.size();ii++) {
    if(Vlines[ii].size()<2) continue;
    vector<string> stoken=split_string(Vlines[ii]);
    if(stoken[0].find("Bonds")!=string::npos) {
      flag=fBond;
      continue;
    }
    else if(stoken[0].find("Angles")!=string::npos) {
      flag=fAngle;
      continue;
    }
    else if(stoken[0].find("Dihedrals")!=string::npos) {
      flag=fDihedral;
      continue;
    }
    else if(stoken[0].find("Impropers")!=string::npos) {
      flag=fImproper;
      continue;
    }
    else if(stoken.size()==3 && stoken[0]=="Atom" && stoken[1]=="Type" && stoken[2]=="Labels") {
      flag=fAtomTypeLabel;
      continue;
    }
    else if(Vlines[ii].find("Coeffs")!=string::npos) {
      flag=fOther;
      continue;
    }
    if(flag==fBond && Vbonds.size()<nbond) {
      if(stoken.size()<4) continue;
      Bond bond_one;
      bond_one.type=stoken[1];
      bond_one.atom_i=stoi(stoken[2]);
      bond_one.atom_j=stod(stoken[3]);
      Vbonds.push_back(bond_one);
    }
    else if(flag==fAngle && Vangles.size()<nangle) {
      if(stoken.size()<5) continue;
      Angle angle_one;
      angle_one.type=stoken[1];
      angle_one.atom_i=stoi(stoken[2]);
      angle_one.atom_j=stoi(stoken[3]);
      angle_one.atom_k=stoi(stoken[4]);
      Vangles.push_back(angle_one);
    }
    else if(flag==fDihedral && Vdihedrals.size()<ndihedral) {
      if(stoken.size()<6) continue;
      Dihedral dihedral_one;
      dihedral_one.type=stoken[1];
      dihedral_one.atom_i=stoi(stoken[2]);
      dihedral_one.atom_j=stoi(stoken[3]);
      dihedral_one.atom_k=stoi(stoken[4]);
      dihedral_one.atom_l=stoi(stoken[5]);
      Vdihedrals.push_back(dihedral_one);
    }
    else if(flag==fImproper && Vimpropers.size()<nimproper) {
      if(stoken.size()<6) continue;
      Dihedral improper_one;
      improper_one.type=stoken[1];
      improper_one.atom_i=stoi(stoken[2]);
      improper_one.atom_j=stoi(stoken[3]);
      improper_one.atom_k=stoi(stoken[4]);
      improper_one.atom_l=stoi(stoken[5]);
      Vimpropers.push_back(improper_one);
    }
    else if(flag==fAtomTypeLabel && atomlabel.size()<natomtype) {
      if(stoken.size()!=2) continue;
      atomlabel.emplace(stoken[0],stoken[1]);
    }
  }
}

void ReadLAMMPS::clear() {
  step_begin=0;
  step_end=-1;
  VVatoms.clear();
  VVatoms.shrink_to_fit();
  VVsteps.clear();
  VVsteps.shrink_to_fit();
  VVenergies.clear();
  VVenergies.shrink_to_fit();
  atomlabel.clear();
  Vbonds.clear();
  Vbonds.shrink_to_fit();
  Vangles.clear();
  Vangles.shrink_to_fit();
  Vdihedrals.clear();
  Vdihedrals.shrink_to_fit();
  Vimpropers.clear();
  Vimpropers.shrink_to_fit();
}