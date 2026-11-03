#include "basic.h"
#include "constant.h"
#include "constant.h"
#include "initialization.h"
#include "error.h"
#include "RNA_coarse_grained.h"

using namespace CBS_NS;

RNACoarseGrained::RNACoarseGrained(Cobmms *cbs) : Pointers(cbs), rna_cg(cbs) {}

RNACoarseGrained::RNACoarseGrained(const RNACoarseGrained &other) : Pointers(other.cbs), rna_cg(other.cbs) {
  rna_cg=other.rna_cg;
  bonds=other.bonds;
  angles=other.angles;
  dihedrals=other.dihedrals;
  impropers=other.impropers;
  atomid_name=other.atomid_name;
}

RNACoarseGrained& RNACoarseGrained::operator =(const RNACoarseGrained &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    rna_cg=other.rna_cg;
    bonds=other.bonds;
    angles=other.angles;
    dihedrals=other.dihedrals;
    impropers=other.impropers;
    atomid_name=other.atomid_name;
  }

  return *this;
}

RNACoarseGrained::~RNACoarseGrained() {
  rna_cg.clear();
  bonds.clear();
  bonds.shrink_to_fit();
  angles.clear();
  angles.shrink_to_fit();
  dihedrals.clear();
  dihedrals.shrink_to_fit();
  impropers.clear();
  impropers.shrink_to_fit();
  atomid_name.clear();
}

// P: P,OP1,OP2,O5',O3'; S: C4',C5',C3',O4'; S2: C2',O2',C1';
// RC: N9,C8,N7,C5,C4,N3; AN: C6,N6; AC: N1,C2; GO: C6,O6,N1; GN: C2,N2;  
// YC: N1,C5,C6; YO: C2,O2; CN: N3,C4,N4; UO: N3,C4,O4;
void RNACoarseGrained::convert_RNA_to_CG(const Molecule &mole_one, size_t index_offset) {
  rna_cg.clear();
  atomid_name.clear();

  size_t index=0;
  for(size_t ii=0;ii<mole_one.chains.size();ii++) {
    if(mole_one.chains[ii].chain_label!=Chain::rna) {
      string str="RNA_coarse_grained.cpp: Not RNA Chains!  --> "+mole_one.chains[ii].chainid;
      error->warning(str);
      continue;
    }

    Chain chain_one(cbs);
    chain_one.chain_label=Chain::cg;
    for(size_t jj=0;jj<mole_one.chains[ii].residues.size();jj++) {
      Residue residue_one=mole_one.chains[ii].residues[jj];
      string chainid_one=residue_one.chainid;
      int resid_one=residue_one.resid;
      string resname_one="c"+residue_one.resname;
      vector<Atom> Vatoms;
      vector<double> xyz;

      if(residue_one.check_atom("P")) {
        Atom atom_P;
        xyz=residue_one["P"].get_atom_xyz();
        index++;
        atom_P.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"P",xyz);
        Vatoms.push_back(atom_P);
      }
      
      if(residue_one.check_atom("C4'")) {
        Atom atom_S;
        xyz=residue_one["C4'"].get_atom_xyz();
        index++;
        atom_S.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"S",xyz);
        Vatoms.push_back(atom_S);
      }

      if(residue_one.check_atom("C2'")) {
        Atom atom_S2;
        xyz=residue_one["C2'"].get_atom_xyz();
        index++;
        atom_S2.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"S2",xyz);
        Vatoms.push_back(atom_S2);
      }
      
      if(residue_one.resname=="A") {
        xyz.clear();
        xyz.resize(3);
        Atom atom_RC;
        xyz[0]=residue_one["N9"].x*MassConst::MASS_N+residue_one["C8"].x*MassConst::MASS_C+residue_one["N7"].x*MassConst::MASS_N
               +residue_one["C5"].x*MassConst::MASS_C+residue_one["C4"].x*MassConst::MASS_C+residue_one["N3"].x*MassConst::MASS_N;
        xyz[1]=residue_one["N9"].y*MassConst::MASS_N+residue_one["C8"].y*MassConst::MASS_C+residue_one["N7"].y*MassConst::MASS_N
               +residue_one["C5"].y*MassConst::MASS_C+residue_one["C4"].y*MassConst::MASS_C+residue_one["N3"].y*MassConst::MASS_N;  
        xyz[2]=residue_one["N9"].z*MassConst::MASS_N+residue_one["C8"].z*MassConst::MASS_C+residue_one["N7"].z*MassConst::MASS_N
               +residue_one["C5"].z*MassConst::MASS_C+residue_one["C4"].z*MassConst::MASS_C+residue_one["N3"].z*MassConst::MASS_N;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_N*3+MassConst::MASS_C*3);
        }
        index++;
        atom_RC.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"RC",xyz);
        Vatoms.push_back(atom_RC);

        Atom atom_AN;
        xyz[0]=residue_one["C6"].x*MassConst::MASS_C+residue_one["N6"].x*MassConst::MASS_N;
        xyz[1]=residue_one["C6"].y*MassConst::MASS_C+residue_one["N6"].y*MassConst::MASS_N;  
        xyz[2]=residue_one["C6"].z*MassConst::MASS_C+residue_one["N6"].z*MassConst::MASS_N;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_C+MassConst::MASS_N);
        }
        index++;
        atom_AN.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"AN",xyz);
        Vatoms.push_back(atom_AN);

        Atom atom_AC;
        xyz[0]=residue_one["N1"].x*MassConst::MASS_N+residue_one["C2"].x*MassConst::MASS_C;
        xyz[1]=residue_one["N1"].y*MassConst::MASS_N+residue_one["C2"].y*MassConst::MASS_C;  
        xyz[2]=residue_one["N1"].z*MassConst::MASS_N+residue_one["C2"].z*MassConst::MASS_C;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_N+MassConst::MASS_C);
        }
        index++;
        atom_AC.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"AC",xyz);
        Vatoms.push_back(atom_AC);
      }
      else if(residue_one.resname=="G") {
        xyz.clear();
        xyz.resize(3);
        Atom atom_RC;
        xyz[0]=residue_one["N9"].x*MassConst::MASS_N+residue_one["C8"].x*MassConst::MASS_C+residue_one["N7"].x*MassConst::MASS_N
               +residue_one["C5"].x*MassConst::MASS_C+residue_one["C4"].x*MassConst::MASS_C+residue_one["N3"].x*MassConst::MASS_N;
        xyz[1]=residue_one["N9"].y*MassConst::MASS_N+residue_one["C8"].y*MassConst::MASS_C+residue_one["N7"].y*MassConst::MASS_N
               +residue_one["C5"].y*MassConst::MASS_C+residue_one["C4"].y*MassConst::MASS_C+residue_one["N3"].y*MassConst::MASS_N;  
        xyz[2]=residue_one["N9"].z*MassConst::MASS_N+residue_one["C8"].z*MassConst::MASS_C+residue_one["N7"].z*MassConst::MASS_N
               +residue_one["C5"].z*MassConst::MASS_C+residue_one["C4"].z*MassConst::MASS_C+residue_one["N3"].z*MassConst::MASS_N;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_N*3+MassConst::MASS_C*3);
        }
        index++;
        atom_RC.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"RC",xyz);
        Vatoms.push_back(atom_RC);

        Atom atom_GO;
        xyz[0]=residue_one["C6"].x*MassConst::MASS_C+residue_one["O6"].x*MassConst::MASS_O+residue_one["N1"].x*MassConst::MASS_N;
        xyz[1]=residue_one["C6"].y*MassConst::MASS_C+residue_one["O6"].y*MassConst::MASS_O+residue_one["N1"].y*MassConst::MASS_N;  
        xyz[2]=residue_one["C6"].z*MassConst::MASS_C+residue_one["O6"].z*MassConst::MASS_O+residue_one["N1"].z*MassConst::MASS_N;    
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_C+MassConst::MASS_O+MassConst::MASS_N);
        }
        index++;
        atom_GO.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"GO",xyz);
        Vatoms.push_back(atom_GO);

        Atom atom_GN;
        xyz[0]=residue_one["C2"].x*MassConst::MASS_C+residue_one["N2"].x*MassConst::MASS_N;
        xyz[1]=residue_one["C2"].y*MassConst::MASS_C+residue_one["N2"].y*MassConst::MASS_N;  
        xyz[2]=residue_one["C2"].z*MassConst::MASS_C+residue_one["N2"].z*MassConst::MASS_N;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_C+MassConst::MASS_N);
        }
        index++;
        atom_GN.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"GN",xyz);
        Vatoms.push_back(atom_GN);
      }
      else if(residue_one.resname=="C") {
        xyz.clear();
        xyz.resize(3);
        Atom atom_YC;
        xyz[0]=residue_one["N1"].x*MassConst::MASS_N+residue_one["C5"].x*MassConst::MASS_C+residue_one["C6"].x*MassConst::MASS_C;
        xyz[1]=residue_one["N1"].y*MassConst::MASS_N+residue_one["C5"].y*MassConst::MASS_C+residue_one["C6"].y*MassConst::MASS_C;  
        xyz[2]=residue_one["N1"].z*MassConst::MASS_N+residue_one["C5"].z*MassConst::MASS_C+residue_one["C6"].z*MassConst::MASS_C;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_N+MassConst::MASS_C*2);
        }
        index++;
        atom_YC.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"YC",xyz);
        Vatoms.push_back(atom_YC);

        Atom atom_YO;
        xyz[0]=residue_one["C2"].x*MassConst::MASS_C+residue_one["O2"].x*MassConst::MASS_O;
        xyz[1]=residue_one["C2"].y*MassConst::MASS_C+residue_one["O2"].y*MassConst::MASS_O;  
        xyz[2]=residue_one["C2"].z*MassConst::MASS_C+residue_one["O2"].z*MassConst::MASS_O;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_C+MassConst::MASS_O);
        }
        index++;
        atom_YO.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"YO",xyz);
        Vatoms.push_back(atom_YO);

        Atom atom_CN;
        xyz[0]=residue_one["N3"].x*MassConst::MASS_N+residue_one["C4"].x*MassConst::MASS_C+residue_one["N4"].x*MassConst::MASS_N;
        xyz[1]=residue_one["N3"].y*MassConst::MASS_N+residue_one["C4"].y*MassConst::MASS_C+residue_one["N4"].y*MassConst::MASS_N;  
        xyz[2]=residue_one["N3"].z*MassConst::MASS_N+residue_one["C4"].z*MassConst::MASS_C+residue_one["N4"].z*MassConst::MASS_N;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_C+MassConst::MASS_N*2);
        }
        index++;
        atom_CN.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"CN",xyz);
        Vatoms.push_back(atom_CN);
      }
      else if(residue_one.resname=="U") {
        xyz.clear();
        xyz.resize(3);
        Atom atom_YC;
        xyz[0]=residue_one["N1"].x*MassConst::MASS_N+residue_one["C5"].x*MassConst::MASS_C+residue_one["C6"].x*MassConst::MASS_C;
        xyz[1]=residue_one["N1"].y*MassConst::MASS_N+residue_one["C5"].y*MassConst::MASS_C+residue_one["C6"].y*MassConst::MASS_C;  
        xyz[2]=residue_one["N1"].z*MassConst::MASS_N+residue_one["C5"].z*MassConst::MASS_C+residue_one["C6"].z*MassConst::MASS_C;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_N+MassConst::MASS_C*2);
        }
        index++;
        atom_YC.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"YC",xyz);
        Vatoms.push_back(atom_YC);

        Atom atom_YO;
        xyz[0]=residue_one["C2"].x*MassConst::MASS_C+residue_one["O2"].x*MassConst::MASS_O;
        xyz[1]=residue_one["C2"].y*MassConst::MASS_C+residue_one["O2"].y*MassConst::MASS_O;  
        xyz[2]=residue_one["C2"].z*MassConst::MASS_C+residue_one["O2"].z*MassConst::MASS_O;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_C+MassConst::MASS_O);
        }
        index++;
        atom_YO.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"YO",xyz);
        Vatoms.push_back(atom_YO);

        Atom atom_UO;
        xyz[0]=residue_one["N3"].x*MassConst::MASS_N+residue_one["C4"].x*MassConst::MASS_C+residue_one["O4"].x*MassConst::MASS_O;
        xyz[1]=residue_one["N3"].y*MassConst::MASS_N+residue_one["C4"].y*MassConst::MASS_C+residue_one["O4"].y*MassConst::MASS_O;  
        xyz[2]=residue_one["N3"].z*MassConst::MASS_N+residue_one["C4"].z*MassConst::MASS_C+residue_one["O4"].z*MassConst::MASS_O;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_N+MassConst::MASS_C+MassConst::MASS_O);
        }
        index++;
        atom_UO.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"UO",xyz);
        Vatoms.push_back(atom_UO);
      }
      else {
        string str="RNA_coarse_grained.cpp: Unknown RNA residue name --> "+residue_one.resname;
        error->warning(str);
      }

      Residue res_cg(cbs);
      res_cg.set_residue(Vatoms);
      chain_one.insert_residue(res_cg,Chain::end);

      for(size_t kk=0;kk<Vatoms.size();kk++) {
        atomid_name.emplace(Vatoms[kk].id,Vatoms[kk].atomname);
      }
    }
    rna_cg.insert_chain(chain_one);
  }
}

void RNACoarseGrained::convert_RNA_to_CG(const RNAMolecule &rna_one, size_t index_offset) {
  rna_cg.clear();
  atomid_name.clear();

  Molecule mole_one(cbs);
  mole_one.set_molecule(rna_one.chains);
  RNACoarseGrained::convert_RNA_to_CG(mole_one,index_offset);
}

// P: P,OP1,OP2,O5',O3'; S: C4',C5',C3',O4',C2',O2',C1';
// RC: N9,C8,N7,C5,C4,N3; AN: C6,N6; AC: N1,C2; GO: C6,O6,N1; GN: C2,N2;  
// YC: N1,C5,C6; YO: C2,O2; CN: N3,C4,N4; UO: N3,C4,O4;
void RNACoarseGrained::convert_RNA_to_CG_v2(const Molecule &mole_one, size_t index_offset) {
  rna_cg.clear();
  atomid_name.clear();

  size_t index=0;
  for(size_t ii=0;ii<mole_one.chains.size();ii++) {
    if(mole_one.chains[ii].chain_label!=Chain::rna) {
      string str="RNA_coarse_grained.cpp: Not RNA Chains!  --> "+mole_one.chains[ii].chainid;
      error->warning(str);
      continue;
    }

    Chain chain_one(cbs);
    chain_one.chain_label=Chain::cg;
    for(size_t jj=0;jj<mole_one.chains[ii].residues.size();jj++) {
      Residue residue_one=mole_one.chains[ii].residues[jj];
      string chainid_one=residue_one.chainid;
      int resid_one=residue_one.resid;
      string resname_one="c"+residue_one.resname;
      vector<Atom> Vatoms;
      vector<double> xyz;

      if(residue_one.check_atom("P")) {
        Atom atom_P;
        xyz=residue_one["P"].get_atom_xyz();
        index++;
        atom_P.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"P",xyz);
        Vatoms.push_back(atom_P);
      }
      
      if(residue_one.check_atom("C4'")) {
        Atom atom_S;
        xyz=residue_one["C4'"].get_atom_xyz();
        index++;
        atom_S.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"S",xyz);
        Vatoms.push_back(atom_S);
      }
      
      if(residue_one.resname=="A") {
        xyz.clear();
        xyz.resize(3);
        Atom atom_RC;
        xyz[0]=residue_one["N9"].x*MassConst::MASS_N+residue_one["C8"].x*MassConst::MASS_C+residue_one["N7"].x*MassConst::MASS_N
               +residue_one["C5"].x*MassConst::MASS_C+residue_one["C4"].x*MassConst::MASS_C+residue_one["N3"].x*MassConst::MASS_N;
        xyz[1]=residue_one["N9"].y*MassConst::MASS_N+residue_one["C8"].y*MassConst::MASS_C+residue_one["N7"].y*MassConst::MASS_N
               +residue_one["C5"].y*MassConst::MASS_C+residue_one["C4"].y*MassConst::MASS_C+residue_one["N3"].y*MassConst::MASS_N;  
        xyz[2]=residue_one["N9"].z*MassConst::MASS_N+residue_one["C8"].z*MassConst::MASS_C+residue_one["N7"].z*MassConst::MASS_N
               +residue_one["C5"].z*MassConst::MASS_C+residue_one["C4"].z*MassConst::MASS_C+residue_one["N3"].z*MassConst::MASS_N;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_N*3+MassConst::MASS_C*3);
        }
        index++;
        atom_RC.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"RC",xyz);
        Vatoms.push_back(atom_RC);

        Atom atom_AN;
        xyz[0]=residue_one["C6"].x*MassConst::MASS_C+residue_one["N6"].x*MassConst::MASS_N;
        xyz[1]=residue_one["C6"].y*MassConst::MASS_C+residue_one["N6"].y*MassConst::MASS_N;  
        xyz[2]=residue_one["C6"].z*MassConst::MASS_C+residue_one["N6"].z*MassConst::MASS_N;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_C+MassConst::MASS_N);
        }
        index++;
        atom_AN.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"AN",xyz);
        Vatoms.push_back(atom_AN);

        Atom atom_AC;
        xyz[0]=residue_one["N1"].x*MassConst::MASS_N+residue_one["C2"].x*MassConst::MASS_C;
        xyz[1]=residue_one["N1"].y*MassConst::MASS_N+residue_one["C2"].y*MassConst::MASS_C;  
        xyz[2]=residue_one["N1"].z*MassConst::MASS_N+residue_one["C2"].z*MassConst::MASS_C;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_N+MassConst::MASS_C);
        }
        index++;
        atom_AC.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"AC",xyz);
        Vatoms.push_back(atom_AC);
      }
      else if(residue_one.resname=="G") {
        xyz.clear();
        xyz.resize(3);
        Atom atom_RC;
        xyz[0]=residue_one["N9"].x*MassConst::MASS_N+residue_one["C8"].x*MassConst::MASS_C+residue_one["N7"].x*MassConst::MASS_N
               +residue_one["C5"].x*MassConst::MASS_C+residue_one["C4"].x*MassConst::MASS_C+residue_one["N3"].x*MassConst::MASS_N;
        xyz[1]=residue_one["N9"].y*MassConst::MASS_N+residue_one["C8"].y*MassConst::MASS_C+residue_one["N7"].y*MassConst::MASS_N
               +residue_one["C5"].y*MassConst::MASS_C+residue_one["C4"].y*MassConst::MASS_C+residue_one["N3"].y*MassConst::MASS_N;  
        xyz[2]=residue_one["N9"].z*MassConst::MASS_N+residue_one["C8"].z*MassConst::MASS_C+residue_one["N7"].z*MassConst::MASS_N
               +residue_one["C5"].z*MassConst::MASS_C+residue_one["C4"].z*MassConst::MASS_C+residue_one["N3"].z*MassConst::MASS_N;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_N*3+MassConst::MASS_C*3);
        }
        index++;
        atom_RC.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"RC",xyz);
        Vatoms.push_back(atom_RC);

        Atom atom_GO;
        xyz[0]=residue_one["C6"].x*MassConst::MASS_C+residue_one["O6"].x*MassConst::MASS_O+residue_one["N1"].x*MassConst::MASS_N;
        xyz[1]=residue_one["C6"].y*MassConst::MASS_C+residue_one["O6"].y*MassConst::MASS_O+residue_one["N1"].y*MassConst::MASS_N;  
        xyz[2]=residue_one["C6"].z*MassConst::MASS_C+residue_one["O6"].z*MassConst::MASS_O+residue_one["N1"].z*MassConst::MASS_N;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_C+MassConst::MASS_O+MassConst::MASS_N);
        }
        index++;
        atom_GO.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"GO",xyz);
        Vatoms.push_back(atom_GO);

        Atom atom_GN;
        xyz[0]=residue_one["C2"].x*MassConst::MASS_C+residue_one["N2"].x*MassConst::MASS_N;
        xyz[1]=residue_one["C2"].y*MassConst::MASS_C+residue_one["N2"].y*MassConst::MASS_N;  
        xyz[2]=residue_one["C2"].z*MassConst::MASS_C+residue_one["N2"].z*MassConst::MASS_N;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_C+MassConst::MASS_N);
        }
        index++;
        atom_GN.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"GN",xyz);
        Vatoms.push_back(atom_GN);
      }
      else if(residue_one.resname=="C") {
        xyz.clear();
        xyz.resize(3);
        Atom atom_YC;
        xyz[0]=residue_one["N1"].x*MassConst::MASS_N+residue_one["C5"].x*MassConst::MASS_C+residue_one["C6"].x*MassConst::MASS_C;
        xyz[1]=residue_one["N1"].y*MassConst::MASS_N+residue_one["C5"].y*MassConst::MASS_C+residue_one["C6"].y*MassConst::MASS_C;  
        xyz[2]=residue_one["N1"].z*MassConst::MASS_N+residue_one["C5"].z*MassConst::MASS_C+residue_one["C6"].z*MassConst::MASS_C;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_N+MassConst::MASS_C*2);
        }
        index++;
        atom_YC.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"YC",xyz);
        Vatoms.push_back(atom_YC);

        Atom atom_YO;
        xyz[0]=residue_one["C2"].x*MassConst::MASS_C+residue_one["O2"].x*MassConst::MASS_O;
        xyz[1]=residue_one["C2"].y*MassConst::MASS_C+residue_one["O2"].y*MassConst::MASS_O;  
        xyz[2]=residue_one["C2"].z*MassConst::MASS_C+residue_one["O2"].z*MassConst::MASS_O;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_C+MassConst::MASS_O);
        }
        index++;
        atom_YO.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"YO",xyz);
        Vatoms.push_back(atom_YO);

        Atom atom_CN;
        xyz[0]=residue_one["N3"].x*MassConst::MASS_N+residue_one["C4"].x*MassConst::MASS_C+residue_one["N4"].x*MassConst::MASS_N;
        xyz[1]=residue_one["N3"].y*MassConst::MASS_N+residue_one["C4"].y*MassConst::MASS_C+residue_one["N4"].y*MassConst::MASS_N;  
        xyz[2]=residue_one["N3"].z*MassConst::MASS_N+residue_one["C4"].z*MassConst::MASS_C+residue_one["N4"].z*MassConst::MASS_N;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_C+MassConst::MASS_N*2);
        }
        index++;
        atom_CN.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"CN",xyz);
        Vatoms.push_back(atom_CN);
      }
      else if(residue_one.resname=="U") {
        xyz.clear();
        xyz.resize(3);
        Atom atom_YC;
        xyz[0]=residue_one["N1"].x*MassConst::MASS_N+residue_one["C5"].x*MassConst::MASS_C+residue_one["C6"].x*MassConst::MASS_C;
        xyz[1]=residue_one["N1"].y*MassConst::MASS_N+residue_one["C5"].y*MassConst::MASS_C+residue_one["C6"].y*MassConst::MASS_C;  
        xyz[2]=residue_one["N1"].z*MassConst::MASS_N+residue_one["C5"].z*MassConst::MASS_C+residue_one["C6"].z*MassConst::MASS_C;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_N+MassConst::MASS_C*2);
        }
        index++;
        atom_YC.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"YC",xyz);
        Vatoms.push_back(atom_YC);

        Atom atom_YO;
        xyz[0]=residue_one["C2"].x*MassConst::MASS_C+residue_one["O2"].x*MassConst::MASS_O;
        xyz[1]=residue_one["C2"].y*MassConst::MASS_C+residue_one["O2"].y*MassConst::MASS_O;  
        xyz[2]=residue_one["C2"].z*MassConst::MASS_C+residue_one["O2"].z*MassConst::MASS_O;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_C+MassConst::MASS_O);
        }
        index++;
        atom_YO.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"YO",xyz);
        Vatoms.push_back(atom_YO);

        Atom atom_UO;
        xyz[0]=residue_one["N3"].x*MassConst::MASS_N+residue_one["C4"].x*MassConst::MASS_C+residue_one["O4"].x*MassConst::MASS_O;
        xyz[1]=residue_one["N3"].y*MassConst::MASS_N+residue_one["C4"].y*MassConst::MASS_C+residue_one["O4"].y*MassConst::MASS_O;  
        xyz[2]=residue_one["N3"].z*MassConst::MASS_N+residue_one["C4"].z*MassConst::MASS_C+residue_one["O4"].z*MassConst::MASS_O;      
        for(size_t kk=0;kk<3;kk++) {
          xyz[kk]=xyz[kk]/(MassConst::MASS_N+MassConst::MASS_C+MassConst::MASS_O);
        }
        index++;
        atom_UO.set_atom(index+index_offset,chainid_one,resid_one,resname_one,"UO",xyz);
        Vatoms.push_back(atom_UO);
      }
      else {
        string str="RNA_coarse_grained.cpp: Unknown RNA residue name --> "+residue_one.resname;
        error->warning(str);
      }

      Residue res_cg(cbs);
      res_cg.set_residue(Vatoms);
      chain_one.insert_residue(res_cg,Chain::end);

      for(size_t kk=0;kk<Vatoms.size();kk++) {
        atomid_name.emplace(Vatoms[kk].id,Vatoms[kk].atomname);
      }
    }
    rna_cg.insert_chain(chain_one);
  }
}

void RNACoarseGrained::convert_RNA_to_CG_v2(const RNAMolecule &rna_one, size_t index_offset) {
  rna_cg.clear();
  atomid_name.clear();

  Molecule mole_one(cbs);
  mole_one.set_molecule(rna_one.chains);
  RNACoarseGrained::convert_RNA_to_CG_v2(mole_one,index_offset);
}

void RNACoarseGrained::reset_RNA_CG(const vector<Atom> &Vatoms, size_t index_offset) {
  rna_cg.clear();
  atomid_name.clear();

  vector<Atom> Tatoms=Vatoms;
  for(size_t ii=0;ii<Tatoms.size();ii++) {
    Tatoms[ii].id=ii+1+index_offset;
    if(atomid_name.count(Tatoms[ii].id)==0)
      atomid_name.emplace(Tatoms[ii].id,Tatoms[ii].atomname);
  }

  rna_cg.set_molecule(Tatoms,Chain::cg);
}

void RNACoarseGrained::reset_RNA_CG(const vector<Atom> &Vatoms, const vector<string> &Vresnames, size_t index_offset) {
  rna_cg.clear();
  atomid_name.clear();

  vector<Atom> Tatoms;
  vector<string> Vatomnames; 
  size_t index=0;
  for(size_t ii=0;ii<Vresnames.size();ii++) {
    string resname_i=Vresnames[ii];
    int check=0;
    for(size_t jj=0;jj<Vatoms.size();jj++) {
      Atom atom_one=Vatoms[jj];
      if(atom_one.resname[0]=='c') atom_one.resname=atom_one.resname.substr(1);
      if(atom_one.get_residue_name()==resname_i) {
        string atomname_one=Vatoms[jj].get_atom_name();
        if(count(Vatomnames.begin(),Vatomnames.end(),atomname_one)==1) continue;
        Vatomnames.push_back(atomname_one);
        index++;
        atom_one=Vatoms[jj];
        atom_one.id=static_cast<int>(index+index_offset);
        Tatoms.push_back(atom_one);
        atomid_name.emplace(atom_one.id,atom_one.atomname);
        check=1;
      }
    }
    if(check==0) {
      string str="RNA_coarse_grained.cpp: Cannot find atoms for residue "+resname_i+"!";
      error->one(str);
    }
  } 

  rna_cg.set_molecule(Tatoms,Chain::cg);
  //for(size_t ii=0;ii<rna_cg.chains.size();ii++) {
    //index_offset=rna_cg.chains[ii].set_atomid(index_offset);
  //}
}

void RNACoarseGrained::reset_RNA_CG(const Molecule &mol_one) {
  rna_cg=mol_one;

  atomid_name.clear();
  for(size_t ii=0;ii<mol_one.chains.size();ii++) {
    for(size_t jj=0;jj<mol_one.chains[ii].residues.size();jj++) {
      for(size_t kk=0;kk<mol_one.chains[ii].residues[jj].atoms.size();kk++) {
        Atom atom_one=mol_one.chains[ii].residues[jj].atoms[kk];
        atomid_name.emplace(atom_one.id,atom_one.atomname);
      }
    }
  }
}

void RNACoarseGrained::build_RNA_topology() {
  string BONDTYPE[]={"P-S","S-P","S-S2","S2-RC","S2-YC","RC-AN","RC-AC","RC-GO","RC-GN","YC-YO","YC-CN","YC-UO","AN-AC","GO-GN","YO-CN","YO-UO"};
  string ANGLETYPE[]={"P-S-P","P-S-S2","S-P-S","S-S2-RC","S-S2-YC","S2-RC-AN","S2-RC-AC","S2-RC-GO","S2-RC-GN",
                      "S2-YC-YO","S2-YC-CN","S2-YC-UO","S2-S-P"}; //,"AN-RC-AC","GO-RC-GN","YO-YC-CN","YO-YC-UO"};
  string DIHEDRALTYPE[]={"P-S-P-S","P-S-S2-RC","P-S-S2-YC","S-S2-RC-AN","S-S2-RC-AC","S-S2-RC-GO","S-S2-RC-GN","S-S2-YC-YO",
                         "S-S2-YC-CN","S-S2-YC-UO","RC-S2-S-P","YC-S2-S-P","S2-S-P-S","S-P-S-P","S-P-S-S2",
                         "S2-RC-AN-AC","S2-RC-GO-GN","S2-YC-CN-YO","S2-YC-UO-YO"};
  string IMPROPERTYPE[]={"S-P-S2-P"};//,"RC-S2-AN-AC","RC-S2-GO-GN","YC-S2-YO-CN","YC-S2-YO-UO"};
  bonds.clear();
  angles.clear();
  dihedrals.clear();
  impropers.clear();

  //find bonds
  for(size_t ii=0;ii<rna_cg.chains.size();ii++) {
    for(size_t jj=0;jj<rna_cg.chains[ii].residues.size();jj++) {
      Residue &residue_one=rna_cg.chains[ii].residues[jj];
      Bond bond_one;
      if(jj>0) {
        bond_one.atom_i=rna_cg.chains[ii].residues[jj-1]["S"].id;
        bond_one.atom_j=residue_one.atoms[0].id;
        bond_one.type="S-"+residue_one.atoms[0].atomname;
        if(bond_one.type!="S-P") {
          string str="RNA_coarse_grained.cpp: Wrong bond type --> "+bond_one.type+" : \"S-P\"";
          error->one(str);
        }
        bonds.push_back(bond_one);
      }

      int index_P=-1;
      if(residue_one.check_atom("P")) {
        index_P=0;
        bond_one.atom_i=residue_one.atoms[0].id;
        bond_one.atom_j=residue_one.atoms[1].id;
        bond_one.type=residue_one.atoms[0].atomname+"-"+residue_one.atoms[1].atomname;
        if(bond_one.type!="P-S") {
          string str="RNA_coarse_grained.cpp: Wrong bond type --> "+bond_one.type+" : \"P-S\"";
          error->one(str);
        }
        bonds.push_back(bond_one);
      }

      bond_one.atom_i=residue_one.atoms[1+index_P].id;
      bond_one.atom_j=residue_one.atoms[2+index_P].id;
      bond_one.type=residue_one.atoms[1+index_P].atomname+"-"+residue_one.atoms[2+index_P].atomname;
      if(bond_one.type!="S-S2") {
        string str="RNA_coarse_grained.cpp: Wrong bond type --> "+bond_one.type+" : \"S-S2\"";
        error->one(str);
      }
      bonds.push_back(bond_one);

      if((int)residue_one.atoms.size()<=5+index_P) continue;

      bond_one.atom_i=residue_one.atoms[2+index_P].id;
      bond_one.atom_j=residue_one.atoms[3+index_P].id;
      bond_one.type=residue_one.atoms[2+index_P].atomname+"-"+residue_one.atoms[3+index_P].atomname;
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+bond_one.type;
        error->one(str);
      }
      bonds.push_back(bond_one);

      bond_one.atom_i=residue_one.atoms[3+index_P].id;
      bond_one.atom_j=residue_one.atoms[4+index_P].id;
      bond_one.type=residue_one.atoms[3+index_P].atomname+"-"+residue_one.atoms[4+index_P].atomname;
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+bond_one.type;
        error->one(str);
      }
      bonds.push_back(bond_one);

      bond_one.atom_i=residue_one.atoms[3+index_P].id;
      bond_one.atom_j=residue_one.atoms[5+index_P].id;
      bond_one.type=residue_one.atoms[3+index_P].atomname+"-"+residue_one.atoms[5+index_P].atomname;
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+bond_one.type;
        error->one(str);
      }
      bonds.push_back(bond_one);

      bond_one.atom_i=residue_one.atoms[4+index_P].id;
      bond_one.atom_j=residue_one.atoms[5+index_P].id;
      bond_one.type=residue_one.atoms[4+index_P].atomname+"-"+residue_one.atoms[5+index_P].atomname;
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+bond_one.type;
        error->one(str);
      }
      bonds.push_back(bond_one);
    }
  }

  map<int,vector<int> > mapBonds;
  for(size_t ii=0;ii<bonds.size();ii++) {
    int atom_i=bonds[ii].atom_i;
    int atom_j=bonds[ii].atom_j;
    if(mapBonds.count(atom_i)==0) {
      vector<int> Vatom_i;
      Vatom_i.push_back(atom_j);
      mapBonds.insert(pair<int,vector<int> >(atom_i,Vatom_i));
    }
    else {
      mapBonds[atom_i].push_back(atom_j);
    }
    if(mapBonds.count(atom_j)==0) {
      vector<int> Vatom_j;
      Vatom_j.push_back(atom_i);
      mapBonds.insert(pair<int,vector<int> >(atom_j,Vatom_j));
    }
    else {
      mapBonds[atom_j].push_back(atom_i);
    } 
  }
  
  //find angles
  for(map<int,vector<int> >::iterator it=mapBonds.begin(); it!=mapBonds.end(); it++) {
    int atom_j=it->first;
    vector<int> neigh_j=it->second;
    size_t nj=neigh_j.size();
    if(nj<2) continue;
    string name_j=atomid_name[atom_j];
    for(size_t ii=0;ii<nj-1;ii++) {
      int atom_i=neigh_j[ii];
      string name_i=atomid_name[atom_i];
      for(size_t kk=ii+1;kk<nj;kk++) {
        int atom_k=neigh_j[kk];
        string name_k=atomid_name[atom_k];
        Angle angle_one;
        if(atom_i<atom_k) {
          angle_one.atom_i=atom_i;
          angle_one.atom_j=atom_j;
          angle_one.atom_k=atom_k;
          angle_one.type=name_i+"-"+name_j+"-"+name_k;
        }
        else {
          angle_one.atom_i=atom_k;
          angle_one.atom_j=atom_j;
          angle_one.atom_k=atom_i;
          angle_one.type=name_k+"-"+name_j+"-"+name_i;
        }
        if(count(begin(ANGLETYPE),end(ANGLETYPE),angle_one.type))
          angles.push_back(angle_one);
      }
    }
  }
  
  //find dihedrals  
  for(map<int,vector<int> >::iterator it=mapBonds.begin(); it!=mapBonds.end(); it++) {
    int atom_j=it->first;
    vector<int> neigh_j=it->second;
    size_t nj=neigh_j.size();
    if(nj<2) continue;
    string name_j=atomid_name[atom_j];
    for(size_t ii=0;ii<nj;ii++) {
      int atom_i=neigh_j[ii];
      string name_i=atomid_name[atom_i];
      for(size_t kk=0;kk<nj;kk++) {
        if(ii==kk) continue;
        int atom_k=neigh_j[kk];
        vector<int> neigh_k=mapBonds[atom_k];
        size_t nk=neigh_k.size();
        if(nk<2) continue;
        string name_k=atomid_name[atom_k];
        for(size_t ll=0;ll<nk;ll++) {
          int atom_l=neigh_k[ll];
          if(atom_l==atom_i||atom_l==atom_j) continue;
          if(atom_l<atom_i) continue;
          string name_l=atomid_name[atom_l];
          Dihedral dihedral_one;
          dihedral_one.atom_i=atom_i;
          dihedral_one.atom_j=atom_j;
          dihedral_one.atom_k=atom_k;
          dihedral_one.atom_l=atom_l;
          dihedral_one.type=name_i+"-"+name_j+"-"+name_k+"-"+name_l;
          if(count(begin(DIHEDRALTYPE),end(DIHEDRALTYPE),dihedral_one.type))
            dihedrals.push_back(dihedral_one);
        }
      }
    }
  }

  //find impropers
  for(auto & it : mapBonds) {
    if(it.second.size()!=3) continue;
    int atom_i=it.first;
    string name_i=atomid_name[atom_i];
    vector<int> atom_jkl=it.second;
    sort(atom_jkl.begin(),atom_jkl.end());
    int atom_j=atom_jkl[0];
    string name_j=atomid_name[atom_j];
    int atom_k=atom_jkl[1];
    string name_k=atomid_name[atom_k];
    int atom_l=atom_jkl[2];
    string name_l=atomid_name[atom_l];
    Dihedral improper_one;
    improper_one.atom_i=atom_i;
    improper_one.atom_j=atom_j;
    improper_one.atom_k=atom_k;
    improper_one.atom_l=atom_l;
    improper_one.type=name_i+"-"+name_j+"-"+name_k+"-"+name_l;
    if(count(begin(IMPROPERTYPE),end(IMPROPERTYPE),improper_one.type))
      impropers.push_back(improper_one);
  }
}

void RNACoarseGrained::build_RNA_topology_v2() {
  string BONDTYPE[]={"P-S","S-P","S-RC","S-YC","RC-AN","RC-AC","RC-GO","RC-GN","YC-YO","YC-CN","YC-UO","AN-AC","GO-GN","YO-CN","YO-UO"};
  string ANGLETYPE[]={"P-S-P","S-P-S","P-S-RC","P-S-YC","S-RC-AN","S-RC-AC","S-RC-GO","S-RC-GN","S-YC-YO","S-YC-CN",
                      "S-YC-UO","RC-S-P","YC-S-P"};
  string DIHEDRALTYPE[]={"P-S-P-S","P-S-RC-AN","P-S-RC-AC","P-S-RC-GO","P-S-RC-GN","P-S-YC-YO","P-S-YC-CN","P-S-YC-UO",
                         "S-RC-AN-AC","S-RC-GO-GN","S-YC-CN-YO","S-YC-UO-YO","S-P-S-P","S-P-S-RC","S-P-S-YC","RC-S-P-S",
                         "AN-RC-S-P","AC-RC-S-P","GO-RC-S-P","GN-RC-S-P","YC-S-P-S","YO-YC-S-P","CN-YC-S-P","UO-YC-S-P"};
  bonds.clear();
  angles.clear();
  dihedrals.clear();
  impropers.clear();

  //find bonds
  for(size_t ii=0;ii<rna_cg.chains.size();ii++) {
    for(size_t jj=0;jj<rna_cg.chains[ii].residues.size();jj++) {
      Residue &residue_one=rna_cg.chains[ii].residues[jj];
      Bond bond_one;
      if(jj>0) {
        bond_one.atom_i=rna_cg.chains[ii].residues[jj-1]["S"].id;
        bond_one.atom_j=residue_one.atoms[0].id;
        bond_one.type="S-"+residue_one.atoms[0].atomname;
        if(bond_one.type!="S-P") {
          string str="RNA_coarse_grained.cpp: Wrong bond type --> "+bond_one.type+" : \"S-P\"";
          error->one(str);
        }
        bonds.push_back(bond_one);
      }

      int index_P=-1;
      if(residue_one.check_atom("P")) {
        index_P=0;
        bond_one.atom_i=residue_one.atoms[0].id;
        bond_one.atom_j=residue_one.atoms[1].id;
        bond_one.type=residue_one.atoms[0].atomname+"-"+residue_one.atoms[1].atomname;
        if(bond_one.type!="P-S") {
          string str="RNA_coarse_grained.cpp: Wrong bond type --> "+bond_one.type+" : \"P-S\"";
          error->one(str);
        }
        bonds.push_back(bond_one);
      }

      if((int)residue_one.atoms.size()<=4+index_P) continue;

      bond_one.atom_i=residue_one.atoms[1+index_P].id;
      bond_one.atom_j=residue_one.atoms[2+index_P].id;
      bond_one.type=residue_one.atoms[1+index_P].atomname+"-"+residue_one.atoms[2+index_P].atomname;
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+bond_one.type;
        error->one(str);
      }
      bonds.push_back(bond_one);

      bond_one.atom_i=residue_one.atoms[2+index_P].id;
      bond_one.atom_j=residue_one.atoms[3+index_P].id;
      bond_one.type=residue_one.atoms[2+index_P].atomname+"-"+residue_one.atoms[3+index_P].atomname;
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+bond_one.type;
        error->one(str);
      }
      bonds.push_back(bond_one);

      bond_one.atom_i=residue_one.atoms[2+index_P].id;
      bond_one.atom_j=residue_one.atoms[4+index_P].id;
      bond_one.type=residue_one.atoms[2+index_P].atomname+"-"+residue_one.atoms[4+index_P].atomname;
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+bond_one.type;
        error->one(str);
      }
      bonds.push_back(bond_one);

      bond_one.atom_i=residue_one.atoms[3+index_P].id;
      bond_one.atom_j=residue_one.atoms[4+index_P].id;
      bond_one.type=residue_one.atoms[3+index_P].atomname+"-"+residue_one.atoms[4+index_P].atomname;
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+bond_one.type;
        error->one(str);
      }
      bonds.push_back(bond_one);
    }
  }

  map<int,vector<int> > mapBonds;
  for(size_t ii=0;ii<bonds.size();ii++) {
    int atom_i=bonds[ii].atom_i;
    int atom_j=bonds[ii].atom_j;
    if(mapBonds.count(atom_i)==0) {
      vector<int> Vatom_i;
      Vatom_i.push_back(atom_j);
      mapBonds.insert(pair<int,vector<int> >(atom_i,Vatom_i));
    }
    else {
      mapBonds[atom_i].push_back(atom_j);
    }
    if(mapBonds.count(atom_j)==0) {
      vector<int> Vatom_j;
      Vatom_j.push_back(atom_i);
      mapBonds.insert(pair<int,vector<int> >(atom_j,Vatom_j));
    }
    else {
      mapBonds[atom_j].push_back(atom_i);
    } 
  }
  
  //find angles
  for(map<int,vector<int> >::iterator it=mapBonds.begin(); it!=mapBonds.end(); it++) {
    int atom_j=it->first;
    vector<int> neigh_j=it->second;
    size_t nj=neigh_j.size();
    if(nj<2) continue;
    string name_j=atomid_name[atom_j];
    for(size_t ii=0;ii<nj-1;ii++) {
      int atom_i=neigh_j[ii];
      string name_i=atomid_name[atom_i];
      for(size_t kk=ii+1;kk<nj;kk++) {
        int atom_k=neigh_j[kk];
        string name_k=atomid_name[atom_k];
        Angle angle_one;
        if(atom_i<atom_k) {
          angle_one.atom_i=atom_i;
          angle_one.atom_j=atom_j;
          angle_one.atom_k=atom_k;
          angle_one.type=name_i+"-"+name_j+"-"+name_k;
        }
        else {
          angle_one.atom_i=atom_k;
          angle_one.atom_j=atom_j;
          angle_one.atom_k=atom_i;
          angle_one.type=name_k+"-"+name_j+"-"+name_i;
        }
        if(count(begin(ANGLETYPE),end(ANGLETYPE),angle_one.type))
          angles.push_back(angle_one);
      }
    }
  }
  
  //find dihedrals  
  for(map<int,vector<int> >::iterator it=mapBonds.begin(); it!=mapBonds.end(); it++) {
    int atom_j=it->first;
    vector<int> neigh_j=it->second;
    size_t nj=neigh_j.size();
    if(nj<2) continue;
    string name_j=atomid_name[atom_j];
    for(size_t ii=0;ii<nj;ii++) {
      int atom_i=neigh_j[ii];
      string name_i=atomid_name[atom_i];
      for(size_t kk=0;kk<nj;kk++) {
        if(ii==kk) continue;
        int atom_k=neigh_j[kk];
        vector<int> neigh_k=mapBonds[atom_k];
        size_t nk=neigh_k.size();
        if(nk<2) continue;
        string name_k=atomid_name[atom_k];
        for(size_t ll=0;ll<nk;ll++) {
          int atom_l=neigh_k[ll];
          if(atom_l==atom_i||atom_l==atom_j) continue;
          if(atom_l<atom_i) continue;
          string name_l=atomid_name[atom_l];
          Dihedral dihedral_one;
          dihedral_one.atom_i=atom_i;
          dihedral_one.atom_j=atom_j;
          dihedral_one.atom_k=atom_k;
          dihedral_one.atom_l=atom_l;
          dihedral_one.type=name_i+"-"+name_j+"-"+name_k+"-"+name_l;
          if(count(begin(DIHEDRALTYPE),end(DIHEDRALTYPE),dihedral_one.type))
            dihedrals.push_back(dihedral_one);
        }
      }
    }
  }
}

void RNACoarseGrained::add_circular_RNA_topology() {
  for(size_t ii=0;ii<rna_cg.chains.size();ii++) {
    if(count(init->CircRNA_chains.begin(),init->CircRNA_chains.end(),rna_cg.chains[ii].chainid)==0) continue;
    Bond bond_one;
    Residue &residue_first=rna_cg.chains[ii].residues[0];
    Residue &residue_last=rna_cg.chains[ii].residues.back();
    bond_one.atom_i=residue_last["S"].id;
    if(!residue_first.check_atom("P")) {
      string str="RNA_coarse_grained.cpp: Cannot find atom P in the first residue!";
      error->one(str);
    }
    bond_one.atom_j=residue_first["P"].id;
    bond_one.type="S-P";
    bonds.push_back(bond_one);

    Angle angle_one;
    angle_one.atom_i=residue_last["P"].id;
    angle_one.atom_j=residue_last["S"].id;
    angle_one.atom_k=residue_first["P"].id;
    angle_one.type="P-S-P";
    angles.push_back(angle_one);
    angle_one.atom_i=residue_last["S2"].id;
    //angle_one.atom_j=residue_last["S"].id;
    //angle_one.atom_k=residue_first["P"].id;
    angle_one.type="S2-S-P";
    angles.push_back(angle_one);
    angle_one.atom_i=residue_last["S"].id;
    angle_one.atom_j=residue_first["P"].id;
    angle_one.atom_k=residue_first["S"].id;
    angle_one.type="S-P-S";
    angles.push_back(angle_one);

    Residue &residue_last_1p=rna_cg.chains[ii].residues[rna_cg.chains[ii].size()-2];
    Residue &residue_first_1n=rna_cg.chains[ii].residues[1];
    Dihedral dihedral_one;
    dihedral_one.atom_i=residue_last_1p["S"].id;
    dihedral_one.atom_j=residue_last["P"].id;
    dihedral_one.atom_k=residue_last["S"].id;
    dihedral_one.atom_l=residue_first["P"].id;
    dihedral_one.type="S-P-S-P";
    dihedrals.push_back(dihedral_one);
    dihedral_one.atom_i=residue_last["P"].id;
    dihedral_one.atom_j=residue_last["S"].id;
    dihedral_one.atom_k=residue_first["P"].id;
    dihedral_one.atom_l=residue_first["S"].id;
    dihedral_one.type="P-S-P-S";
    dihedrals.push_back(dihedral_one);
    dihedral_one.atom_i=residue_last["S2"].id;
    dihedral_one.type="S2-S-P-S";
    dihedrals.push_back(dihedral_one);
    dihedral_one.atom_i=residue_last["S"].id;
    dihedral_one.atom_j=residue_first["P"].id;
    dihedral_one.atom_k=residue_first["S"].id;
    dihedral_one.atom_l=residue_first["S2"].id;
    dihedral_one.type="S-P-S-S2";
    dihedrals.push_back(dihedral_one);
    dihedral_one.atom_l=residue_first_1n["P"].id;
    dihedral_one.type="S-P-S-P";
    dihedrals.push_back(dihedral_one);
  }
}

void RNACoarseGrained::add_circular_RNA_topology_v2() {
  for(size_t ii=0;ii<rna_cg.chains.size();ii++) {
    if(count(init->CircRNA_chains.begin(),init->CircRNA_chains.end(),rna_cg.chains[ii].chainid)==0) continue;
    Bond bond_one;
    Residue &residue_first=rna_cg.chains[ii].residues[0];
    Residue &residue_last=rna_cg.chains[ii].residues.back();
    bond_one.atom_i=residue_last["S"].id;
    if(!residue_first.check_atom("P")) {
      string str="RNA_coarse_grained.cpp: Cannot find atom P in the first residue!";
      error->one(str);
    }
    bond_one.atom_j=residue_first["P"].id;
    bond_one.type="S-P";
    bonds.push_back(bond_one);

    Angle angle_one;
    angle_one.atom_i=residue_last["P"].id;
    angle_one.atom_j=residue_last["S"].id;
    angle_one.atom_k=residue_first["P"].id;
    angle_one.type="P-S-P";
    angles.push_back(angle_one);
    if(residue_last.check_atom("RC")) {
      angle_one.atom_i=residue_last["RC"].id;
      angle_one.type="RC-S-P";
    }
    else if(residue_last.check_atom("YC")) {
      angle_one.atom_i=residue_last["YC"].id;
      angle_one.type="YC-S-P";
    }
    angles.push_back(angle_one);
    angle_one.atom_i=residue_last["S"].id;
    angle_one.atom_j=residue_first["P"].id;
    angle_one.atom_k=residue_first["S"].id;
    angle_one.type="S-P-S";
    angles.push_back(angle_one);

    Residue &residue_last_1p=rna_cg.chains[ii].residues[rna_cg.chains[ii].size()-2];
    Residue &residue_first_1n=rna_cg.chains[ii].residues[1];
    Dihedral dihedral_one;
    dihedral_one.atom_i=residue_last_1p["S"].id;
    dihedral_one.atom_j=residue_last["P"].id;
    dihedral_one.atom_k=residue_last["S"].id;
    dihedral_one.atom_l=residue_first["P"].id;
    dihedral_one.type="S-P-S-P";
    dihedrals.push_back(dihedral_one);
    dihedral_one.atom_i=residue_last["P"].id;
    dihedral_one.atom_j=residue_last["S"].id;
    dihedral_one.atom_k=residue_first["P"].id;
    dihedral_one.atom_l=residue_first["S"].id;
    dihedral_one.type="P-S-P-S";
    dihedrals.push_back(dihedral_one);
    if(residue_last.check_atom("RC")) {
      dihedral_one.atom_i=residue_last["RC"].id;
      dihedral_one.type="RC-S-P-S";
      dihedrals.push_back(dihedral_one);
    }
    else if(residue_last.check_atom("YC")) {
      dihedral_one.atom_i=residue_last["YC"].id;
      dihedral_one.type="YC-S-P-S";
      dihedrals.push_back(dihedral_one);
    }
    dihedral_one.atom_k=residue_last["S"].id;
    dihedral_one.atom_l=residue_first["P"].id;
    if(residue_last.check_atom("AN")) {
      dihedral_one.atom_j=residue_last["RC"].id;
      dihedral_one.atom_i=residue_last["AN"].id;
      dihedral_one.type="AN-RC-S-P";
      dihedrals.push_back(dihedral_one);
      dihedral_one.atom_i=residue_last["AC"].id;
      dihedral_one.type="AC-RC-S-P";
      dihedrals.push_back(dihedral_one);
    }
    else if(residue_last.check_atom("GO")) {
      dihedral_one.atom_j=residue_last["RC"].id;
      dihedral_one.atom_i=residue_last["GO"].id;
      dihedral_one.type="GO-RC-S-P";
      dihedrals.push_back(dihedral_one);
      dihedral_one.atom_i=residue_last["GN"].id;
      dihedral_one.type="GN-RC-S-P";
      dihedrals.push_back(dihedral_one);
    }
    else if(residue_last.check_atom("CN")) {
      dihedral_one.atom_j=residue_last["YC"].id;
      dihedral_one.atom_i=residue_last["YO"].id;
      dihedral_one.type="YO-YC-S-P";
      dihedrals.push_back(dihedral_one);
      dihedral_one.atom_i=residue_last["CN"].id;
      dihedral_one.type="CN-YC-S-P";
      dihedrals.push_back(dihedral_one);
    }
    else if(residue_last.check_atom("UO")) {
      dihedral_one.atom_j=residue_last["YC"].id;
      dihedral_one.atom_i=residue_last["YO"].id;
      dihedral_one.type="YO-YC-S-P";
      dihedrals.push_back(dihedral_one);
      dihedral_one.atom_i=residue_last["UO"].id;
      dihedral_one.type="UO-YC-S-P";
      dihedrals.push_back(dihedral_one);
    }
    dihedral_one.atom_i=residue_last["S"].id;
    dihedral_one.atom_j=residue_first["P"].id;
    dihedral_one.atom_k=residue_first["S"].id;
    if(residue_first.check_atom("RC")) {
      dihedral_one.atom_l=residue_first["RC"].id;
      dihedral_one.type="S-P-S-RC";
      dihedrals.push_back(dihedral_one);
    }
    else if(residue_first.check_atom("YC")) {
      dihedral_one.atom_l=residue_first["YC"].id;
      dihedral_one.type="S-P-S-YC";
      dihedrals.push_back(dihedral_one);
    }
    dihedral_one.atom_l=residue_first_1n["P"].id;
    dihedral_one.type="S-P-S-P";
    dihedrals.push_back(dihedral_one);
  }
}

void RNACoarseGrained::add_RNA_basepairs(const map<string,string> &mapBPs, string type_app) {
  string BONDTYPE[]={"GO-CN","GN-YO","AN-UO","AC-YO","GO-UO","GO-YO","RC-YC","S2-S2","S-S","P-P"};
  string BPTYPE[]={"GC","CG","AU","UA","GU","UG","cGcC","cCcG","cAcU","cUcA","cGcU","cUcG"};
  for(auto & it : mapBPs) {
    vector<string> stoken=split_string(it.first,"/");
    if(stoken.size()!=3) {
      string str="RNA_coarse_grained.cpp: Wrong residue name(chainid/resid/resname) --> "+it.first;
      error->one(str);
    }
    string resname_i=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    Residue residue_i=rna_cg.find_residue(resname_i);
    string bp_type=stoken[2];
    stoken=split_string(it.second,"/");
    if(stoken.size()!=3) {
      string str="RNA_coarse_grained.cpp: Wrong residue name(chainid/resid/resname) --> "+it.second;
      error->one(str);
    }
    string resname_j=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    Residue residue_j=rna_cg.find_residue(resname_j);
    bp_type+=stoken[2];
    if(count(begin(BPTYPE),end(BPTYPE),bp_type)==0) {
      string str="RNA_coarse_grained.cpp: Skip noncanonical base pair --> "+bp_type;
      error->warning(str);
      continue;
    }
    int bp_i_1=0,bp_i_2=0,bp_j_1=0,bp_j_2=0;
    int bp_i_3=0, bp_j_3=0;
    if(bp_type=="GC" || bp_type=="cGcC") {
      bp_i_1=residue_i["GO"].id; 
      bp_i_2=residue_i["GN"].id;
      bp_i_3=residue_i["RC"].id;
      bp_j_1=residue_j["CN"].id;
      bp_j_2=residue_j["YO"].id;
      bp_j_3=residue_j["YC"].id;
    }
    else if(bp_type=="CG" || bp_type=="cCcG") {
      bp_i_1=residue_j["GO"].id; 
      bp_i_2=residue_j["GN"].id;
      bp_i_3=residue_j["RC"].id;
      bp_j_1=residue_i["CN"].id;
      bp_j_2=residue_i["YO"].id;
      bp_j_3=residue_i["YC"].id;
    }
    else if(bp_type=="AU" || bp_type=="cAcU") {
      bp_i_1=residue_i["AN"].id; 
      bp_i_2=residue_i["AC"].id;
      bp_i_3=residue_i["RC"].id;
      bp_j_1=residue_j["UO"].id;
      bp_j_2=residue_j["YO"].id;
      bp_j_3=residue_j["YC"].id;
    }
    else if(bp_type=="UA" || bp_type=="cUcA") {
      bp_i_1=residue_j["AN"].id; 
      bp_i_2=residue_j["AC"].id;
      bp_i_3=residue_j["RC"].id;
      bp_j_1=residue_i["UO"].id;
      bp_j_2=residue_i["YO"].id;
      bp_j_3=residue_i["YC"].id;
    }
    else if(bp_type=="GU" || bp_type=="cGcU") {
      bp_i_1=residue_i["GO"].id; 
      bp_i_2=residue_i["GO"].id;
      bp_i_3=residue_i["RC"].id;
      bp_j_1=residue_j["UO"].id;
      bp_j_2=residue_j["YO"].id;
      bp_j_3=residue_j["YC"].id;
    }
    else if(bp_type=="UG" || bp_type=="cUcG") {
      bp_i_1=residue_j["GO"].id; 
      bp_i_2=residue_j["GO"].id;
      bp_i_3=residue_j["RC"].id;
      bp_j_1=residue_i["UO"].id;
      bp_j_2=residue_i["YO"].id;
      bp_j_3=residue_i["YC"].id;
    }
    
    Bond bond_one;
    bond_one.atom_i=bp_i_1;
    bond_one.atom_j=bp_j_1;
    bond_one.type=atomid_name[bp_i_1]+"-"+atomid_name[bp_j_1];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+it.first+" "+it.second;
      str+=" "+to_string(bp_i_1)+" "+to_string(bp_j_1)+" "+bond_one.type;
      error->one(str);
    }
    if(type_app.size()>0)
      bond_one.type+="_"+type_app;
    bonds.push_back(bond_one);
    
    bond_one.atom_i=bp_i_2;
    bond_one.atom_j=bp_j_2;
    bond_one.type=atomid_name[bp_i_2]+"-"+atomid_name[bp_j_2];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+it.first+" "+it.second;
      str+=" "+to_string(bp_i_2)+" "+to_string(bp_j_2)+" "+bond_one.type;
      error->one(str);
    }
    if(type_app.size()>0)
      bond_one.type+="_"+type_app;
    bonds.push_back(bond_one);

    bond_one.atom_i=bp_i_3;
    bond_one.atom_j=bp_j_3;
    bond_one.type=atomid_name[bp_i_3]+"-"+atomid_name[bp_j_3];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+it.first+" "+it.second;
      str+=" "+to_string(bp_i_3)+" "+to_string(bp_j_3)+" "+bond_one.type;
      error->one(str);
    }
    if(type_app.size()>0)
      bond_one.type+="_"+type_app;
    bonds.push_back(bond_one);

    if(residue_i.check_atom("S2") && residue_j.check_atom("S2")) {
      bond_one.atom_i=residue_i["S2"].id;
      bond_one.atom_j=residue_j["S2"].id;
      bond_one.type=atomid_name[bond_one.atom_i]+"-"+atomid_name[bond_one.atom_j];
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+it.first+" "+it.second;
        str+=" "+to_string(bond_one.atom_i)+" "+to_string(bond_one.atom_j)+" "+bond_one.type;
        error->one(str);
      }
      if(type_app.size()>0)
        bond_one.type+="_"+type_app;
      bonds.push_back(bond_one);
    }

    if(residue_i.check_atom("S") && residue_j.check_atom("S")) {
      bond_one.atom_i=residue_i["S"].id;
      bond_one.atom_j=residue_j["S"].id;
      bond_one.type=atomid_name[bond_one.atom_i]+"-"+atomid_name[bond_one.atom_j];
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+it.first+" "+it.second;
        str+=" "+to_string(bond_one.atom_i)+" "+to_string(bond_one.atom_j)+" "+bond_one.type;
        error->one(str);
      }
      if(type_app.size()>0)
        bond_one.type+="_"+type_app;
      bonds.push_back(bond_one);
    }

    if(residue_i.check_atom("P") && residue_j.check_atom("P")) {
      bond_one.atom_i=residue_i["P"].id;
      bond_one.atom_j=residue_j["P"].id;
      bond_one.type=atomid_name[bond_one.atom_i]+"-"+atomid_name[bond_one.atom_j];
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+it.first+" "+it.second;
        str+=" "+to_string(bond_one.atom_i)+" "+to_string(bond_one.atom_j)+" "+bond_one.type;
        error->one(str);
      }
      if(type_app.size()>0)
        bond_one.type+="_"+type_app;
      bonds.push_back(bond_one);
    }
  }//for-it
}

void RNACoarseGrained::add_RNA_basepairs(const vector<ResiduePair> &Vbps, string type_app) {
  string BONDTYPE[]={"GO-CN","GN-YO","AN-UO","AC-YO","GO-UO","GO-YO","RC-YC","S2-S2","S-S","P-P"};
  string BPTYPE[]={"GC","CG","AU","UA","GU","UG","cGcC","cCcG","cAcU","cUcA","cGcU","cUcG"};
  for(size_t ii=0;ii<Vbps.size();ii++) {
    vector<string> stoken=split_string(Vbps[ii].residue_i,"/");
    if(stoken.size()!=3) {
      string str="RNA_coarse_grained.cpp: Wrong residue name(chainid/resid/resname) --> "+Vbps[ii].residue_i;
      error->one(str);
    }
    string resname_i=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    Residue residue_i=rna_cg.find_residue(resname_i);
    string bp_type=stoken[2];
    stoken=split_string(Vbps[ii].residue_j,"/");
    if(stoken.size()!=3) {
      string str="RNA_coarse_grained.cpp: Wrong residue name(chainid/resid/resname) --> "+Vbps[ii].residue_j;
      error->one(str);
    }
    string resname_j=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    Residue residue_j=rna_cg.find_residue(resname_j);
    bp_type+=stoken[2];
    if(count(begin(BPTYPE),end(BPTYPE),bp_type)==0) {
      string str="RNA_coarse_grained.cpp: Skip noncanonical base pair --> "+bp_type;
      error->warning(str);
      continue;
    }
    int bp_i_1=0,bp_i_2=0,bp_j_1=0,bp_j_2=0;
    int bp_i_3=0, bp_j_3=0;
    if(bp_type=="GC" || bp_type=="cGcC") {
      bp_i_1=residue_i["GO"].id; 
      bp_i_2=residue_i["GN"].id;
      bp_i_3=residue_i["RC"].id;
      bp_j_1=residue_j["CN"].id;
      bp_j_2=residue_j["YO"].id;
      bp_j_3=residue_j["YC"].id;
    }
    else if(bp_type=="CG" || bp_type=="cCcG") {
      bp_i_1=residue_j["GO"].id; 
      bp_i_2=residue_j["GN"].id;
      bp_i_3=residue_j["RC"].id;
      bp_j_1=residue_i["CN"].id;
      bp_j_2=residue_i["YO"].id;
      bp_j_3=residue_i["YC"].id;
    }
    else if(bp_type=="AU" || bp_type=="cAcU") {
      bp_i_1=residue_i["AN"].id; 
      bp_i_2=residue_i["AC"].id;
      bp_i_3=residue_i["RC"].id;
      bp_j_1=residue_j["UO"].id;
      bp_j_2=residue_j["YO"].id;
      bp_j_3=residue_j["YC"].id;
    }
    else if(bp_type=="UA" || bp_type=="cUcA") {
      bp_i_1=residue_j["AN"].id; 
      bp_i_2=residue_j["AC"].id;
      bp_i_3=residue_j["RC"].id;
      bp_j_1=residue_i["UO"].id;
      bp_j_2=residue_i["YO"].id;
      bp_j_3=residue_i["YC"].id;
    }
    else if(bp_type=="GU" || bp_type=="cGcU") {
      bp_i_1=residue_i["GO"].id; 
      bp_i_2=residue_i["GO"].id;
      bp_i_3=residue_i["RC"].id;
      bp_j_1=residue_j["UO"].id;
      bp_j_2=residue_j["YO"].id;
      bp_j_3=residue_j["YC"].id;
    }
    else if(bp_type=="UG" || bp_type=="cUcG") {
      bp_i_1=residue_j["GO"].id; 
      bp_i_2=residue_j["GO"].id;
      bp_i_3=residue_j["RC"].id;
      bp_j_1=residue_i["UO"].id;
      bp_j_2=residue_i["YO"].id;
      bp_j_3=residue_i["YC"].id;
    }
    
    Bond bond_one;
    bond_one.atom_i=bp_i_1;
    bond_one.atom_j=bp_j_1;
    bond_one.type=atomid_name[bp_i_1]+"-"+atomid_name[bp_j_1];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+Vbps[ii].residue_i+" "+Vbps[ii].residue_j;
      str+=" "+to_string(bp_i_1)+" "+to_string(bp_j_1)+" "+bond_one.type;
      error->one(str);
    }
    if(type_app.size()>0)
      bond_one.type+="_"+type_app;
    bonds.push_back(bond_one);
    
    bond_one.atom_i=bp_i_2;
    bond_one.atom_j=bp_j_2;
    bond_one.type=atomid_name[bp_i_2]+"-"+atomid_name[bp_j_2];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+Vbps[ii].residue_i+" "+Vbps[ii].residue_j;
      str+=" "+to_string(bp_i_2)+" "+to_string(bp_j_2)+" "+bond_one.type;
      error->one(str);
    }
    if(type_app.size()>0)
      bond_one.type+="_"+type_app;
    bonds.push_back(bond_one);

    bond_one.atom_i=bp_i_3;
    bond_one.atom_j=bp_j_3;
    bond_one.type=atomid_name[bp_i_3]+"-"+atomid_name[bp_j_3];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+Vbps[ii].residue_i+" "+Vbps[ii].residue_j;
      str+=" "+to_string(bp_i_3)+" "+to_string(bp_j_3)+" "+bond_one.type;
      error->one(str);
    }
    if(type_app.size()>0)
      bond_one.type+="_"+type_app;
    bonds.push_back(bond_one);

    if(residue_i.check_atom("S2") && residue_j.check_atom("S2")) {
      bond_one.atom_i=residue_i["S2"].id;
      bond_one.atom_j=residue_j["S2"].id;
      bond_one.type=atomid_name[bond_one.atom_i]+"-"+atomid_name[bond_one.atom_j];
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+Vbps[ii].residue_i+" "+Vbps[ii].residue_j;
        str+=" "+to_string(bond_one.atom_i)+" "+to_string(bond_one.atom_j)+" "+bond_one.type;
        error->one(str);
      }
      if(type_app.size()>0)
        bond_one.type+="_"+type_app;
      bonds.push_back(bond_one);
    }

    if(residue_i.check_atom("S") && residue_j.check_atom("S")) {
      bond_one.atom_i=residue_i["S"].id;
      bond_one.atom_j=residue_j["S"].id;
      bond_one.type=atomid_name[bond_one.atom_i]+"-"+atomid_name[bond_one.atom_j];
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+Vbps[ii].residue_i+" "+Vbps[ii].residue_j;
        str+=" "+to_string(bond_one.atom_i)+" "+to_string(bond_one.atom_j)+" "+bond_one.type;
        error->one(str);
      }
      if(type_app.size()>0)
        bond_one.type+="_"+type_app;
      bonds.push_back(bond_one);
    }

    if(residue_i.check_atom("P") && residue_j.check_atom("P")) {
      bond_one.atom_i=residue_i["P"].id;
      bond_one.atom_j=residue_j["P"].id;
      bond_one.type=atomid_name[bond_one.atom_i]+"-"+atomid_name[bond_one.atom_j];
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+Vbps[ii].residue_i+" "+Vbps[ii].residue_j;
        str+=" "+to_string(bond_one.atom_i)+" "+to_string(bond_one.atom_j)+" "+bond_one.type;
        error->one(str);
      }
      if(type_app.size()>0)
        bond_one.type+="_"+type_app;
      bonds.push_back(bond_one);
    }
  }
}

void RNACoarseGrained::add_RNA_basepairs_v2(const map<string,string> &mapBPs, string type_app) {
  string BONDTYPE[]={"GO-CN","GN-YO","AN-UO","AC-YO","GO-UO","GO-YO","S-S","P-P"};
  string ANGLETYPE[]={"GN-YO-CN","AC-YO-UO","RC-GO-UO"};
  string DIHEDRALTYPE[]={"RC-GO-GN-YO","GO-GN-YO-CN","GN-YO-CN-YC","RC-GO-CN-YC","RC-GN-YO-YC",
                         "RC-AN-AC-YO","AN-AC-YO-UO","AC-YO-UO-YC","RC-AN-UO-YC","RC-AC-YO-YC",
                         "RC-GN-GO-YO","GN-GO-YO-UO","GO-YO-UO-YC","RC-GO-UO-YC"};
  string BPTYPE[]={"GC","CG","AU","UA","GU","UG","cGcC","cCcG","cAcU","cUcA","cGcU","cUcG"};
  for(auto & it : mapBPs) {
    vector<string> stoken=split_string(it.first,"/");
    if(stoken.size()!=3) {
      string str="RNA_coarse_grained.cpp: Wrong residue name(chainid/resid/resname) --> "+it.first;
      error->one(str);
    }
    string resname_i=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    Residue residue_i=rna_cg.find_residue(resname_i);
    string bp_type=stoken[2];
    stoken=split_string(it.second,"/");
    if(stoken.size()!=3) {
      string str="RNA_coarse_grained.cpp: Wrong residue name(chainid/resid/resname) --> "+it.second;
      error->one(str);
    }
    string resname_j=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    Residue residue_j=rna_cg.find_residue(resname_j);
    bp_type+=stoken[2];
    if(count(begin(BPTYPE),end(BPTYPE),bp_type)==0) {
      string str="RNA_coarse_grained.cpp: Skip noncanonical base pair --> "+bp_type;
      error->warning(str);
      continue;
    }
    int bp_i_1=0,bp_i_2=0,bp_i_3=0,bp_j_1=0,bp_j_2=0,bp_j_3=0;
    int flag_GU=0;
    if(bp_type=="GC" || bp_type=="cGcC") {
      bp_i_1=residue_i["GO"].id; 
      bp_i_2=residue_i["GN"].id;
      bp_i_3=residue_i["RC"].id;
      bp_j_1=residue_j["CN"].id;
      bp_j_2=residue_j["YO"].id;
      bp_j_3=residue_j["YC"].id;
    }
    else if(bp_type=="CG" || bp_type=="cCcG") {
      bp_i_1=residue_j["GO"].id; 
      bp_i_2=residue_j["GN"].id;
      bp_i_3=residue_j["RC"].id;
      bp_j_1=residue_i["CN"].id;
      bp_j_2=residue_i["YO"].id;
      bp_j_3=residue_i["YC"].id;
    }
    else if(bp_type=="AU" || bp_type=="cAcU") {
      bp_i_1=residue_i["AN"].id; 
      bp_i_2=residue_i["AC"].id;
      bp_i_3=residue_i["RC"].id;
      bp_j_1=residue_j["UO"].id;
      bp_j_2=residue_j["YO"].id;
      bp_j_3=residue_j["YC"].id;
    }
    else if(bp_type=="UA" || bp_type=="cUcA") {
      bp_i_1=residue_j["AN"].id; 
      bp_i_2=residue_j["AC"].id;
      bp_i_3=residue_j["RC"].id;
      bp_j_1=residue_i["UO"].id;
      bp_j_2=residue_i["YO"].id;
      bp_j_3=residue_i["YC"].id;
    }
    else if(bp_type=="GU" || bp_type=="cGcU") {
      bp_i_1=residue_i["GO"].id; 
      bp_i_2=residue_i["GN"].id;
      bp_i_3=residue_i["RC"].id;
      bp_j_1=residue_j["UO"].id;
      bp_j_2=residue_j["YO"].id;
      bp_j_3=residue_j["YC"].id;
      flag_GU=1;
    }
    else if(bp_type=="UG" || bp_type=="cUcG") {
      bp_i_1=residue_j["GO"].id; 
      bp_i_2=residue_j["GN"].id;
      bp_i_3=residue_j["RC"].id;
      bp_j_1=residue_i["UO"].id;
      bp_j_2=residue_i["YO"].id;
      bp_j_3=residue_i["YC"].id;
      flag_GU=1;
    }
    
    Bond bond_one;
    bond_one.atom_i=bp_i_1;
    bond_one.atom_j=bp_j_1;
    bond_one.type=atomid_name[bp_i_1]+"-"+atomid_name[bp_j_1];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+it.first+" "+it.second;
      str+=" "+to_string(bp_i_1)+" "+to_string(bp_j_1)+" "+bond_one.type;
      error->one(str);
    }
    if(type_app.size()>0)
      bond_one.type+="_"+type_app;
    bonds.push_back(bond_one);
    
    if(flag_GU==0)
      bond_one.atom_i=bp_i_2;
    else 
      bond_one.atom_i=bp_i_1;
    bond_one.atom_j=bp_j_2;
    bond_one.type=atomid_name[bond_one.atom_i]+"-"+atomid_name[bond_one.atom_j];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+it.first+" "+it.second;
      str+=" "+to_string(bond_one.atom_i)+" "+to_string(bond_one.atom_j)+" "+bond_one.type;
      error->one(str);
    }
    if(type_app.size()>0)
      bond_one.type+="_"+type_app;
    bonds.push_back(bond_one);

    if(residue_i.check_atom("S") && residue_j.check_atom("S")) {
      bond_one.atom_i=residue_i["S"].id;
      bond_one.atom_j=residue_j["S"].id;
      bond_one.type=atomid_name[bond_one.atom_i]+"-"+atomid_name[bond_one.atom_j];
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+it.first+" "+it.second;
        str+=" "+to_string(bond_one.atom_i)+" "+to_string(bond_one.atom_j)+" "+bond_one.type;
        error->one(str);
      }
      if(type_app.size()>0)
        bond_one.type+="_"+type_app;
      bonds.push_back(bond_one);
    }

    if(residue_i.check_atom("P") && residue_j.check_atom("P")) {
      bond_one.atom_i=residue_i["P"].id;
      bond_one.atom_j=residue_j["P"].id;
      bond_one.type=atomid_name[bond_one.atom_i]+"-"+atomid_name[bond_one.atom_j];
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+it.first+" "+it.second;
        str+=" "+to_string(bond_one.atom_i)+" "+to_string(bond_one.atom_j)+" "+bond_one.type;
        error->one(str);
      }
      if(type_app.size()>0)
        bond_one.type+="_"+type_app;
      bonds.push_back(bond_one);
    }
    
    Angle angle_one;
    if(flag_GU==0) {
      angle_one.atom_i=bp_i_2;
      angle_one.atom_j=bp_j_2;
      angle_one.atom_k=bp_j_1;
    }
    else {
      angle_one.atom_i=bp_i_3;
      angle_one.atom_j=bp_i_1;
      angle_one.atom_k=bp_j_1;
    }
    angle_one.type=atomid_name[angle_one.atom_i]+"-"+atomid_name[angle_one.atom_j]+"-"+atomid_name[angle_one.atom_k];
    if(count(begin(ANGLETYPE),end(ANGLETYPE),angle_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown angle type --> "+it.first+" "+it.second;
      str+=" "+to_string(angle_one.atom_i)+" "+to_string(angle_one.atom_j)+" "+to_string(angle_one.atom_k);
      str+=" "+angle_one.type;
      error->one(str);
    }
    angles.push_back(angle_one);

    vector<Dihedral> Tbpdihedral;
    Dihedral dihedral_one;
    dihedral_one.type="";
    if(flag_GU==0) {
      dihedral_one.atom_i=bp_i_3;
      dihedral_one.atom_j=bp_i_1;
      dihedral_one.atom_k=bp_i_2;
      dihedral_one.atom_l=bp_j_2;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_1;
      dihedral_one.atom_j=bp_i_2;
      dihedral_one.atom_k=bp_j_2;
      dihedral_one.atom_l=bp_j_1;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_2;
      dihedral_one.atom_j=bp_j_2;
      dihedral_one.atom_k=bp_j_1;
      dihedral_one.atom_l=bp_j_3;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_3;
      dihedral_one.atom_j=bp_i_1;
      dihedral_one.atom_k=bp_j_1;
      dihedral_one.atom_l=bp_j_3;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_3;
      dihedral_one.atom_j=bp_i_2;
      dihedral_one.atom_k=bp_j_2;
      dihedral_one.atom_l=bp_j_3;
      Tbpdihedral.push_back(dihedral_one);
    }
    else {
      dihedral_one.atom_i=bp_i_3;
      dihedral_one.atom_j=bp_i_2;
      dihedral_one.atom_k=bp_i_1;
      dihedral_one.atom_l=bp_j_2;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_2;
      dihedral_one.atom_j=bp_i_1;
      dihedral_one.atom_k=bp_j_2;
      dihedral_one.atom_l=bp_j_1;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_1;
      dihedral_one.atom_j=bp_j_2;
      dihedral_one.atom_k=bp_j_1;
      dihedral_one.atom_l=bp_j_3;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_3;
      dihedral_one.atom_j=bp_i_1;
      dihedral_one.atom_k=bp_j_1;
      dihedral_one.atom_l=bp_j_3;
      Tbpdihedral.push_back(dihedral_one);
    }

    for(size_t kk=0;kk<Tbpdihedral.size();kk++) {
      dihedral_one=Tbpdihedral[kk];
      dihedral_one.type=atomid_name[dihedral_one.atom_i]+"-"+atomid_name[dihedral_one.atom_j]+"-";
      dihedral_one.type+=atomid_name[dihedral_one.atom_k]+"-"+atomid_name[dihedral_one.atom_l];
      if(count(begin(DIHEDRALTYPE),end(DIHEDRALTYPE),dihedral_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown dihedral type --> "+it.first+" "+it.second;
        str+=" "+to_string(dihedral_one.atom_i)+" "+to_string(dihedral_one.atom_j)+" ";
        str+=to_string(dihedral_one.atom_k)+" "+to_string(dihedral_one.atom_l)+" "+dihedral_one.type;
        error->one(str);
      }
      dihedrals.push_back(dihedral_one);
    }
  }//for-it

}

void RNACoarseGrained::add_RNA_basepairs_v2(const vector<ResiduePair> &Vbps, string type_app) {
  string BONDTYPE[]={"GO-CN","GN-YO","AN-UO","AC-YO","GO-UO","GO-YO","S-S","P-P"};
  string ANGLETYPE[]={"GN-YO-CN","AC-YO-UO","RC-GO-UO"};
  string DIHEDRALTYPE[]={"RC-GO-GN-YO","GO-GN-YO-CN","GN-YO-CN-YC","RC-GO-CN-YC","RC-GN-YO-YC",
                         "RC-AN-AC-YO","AN-AC-YO-UO","AC-YO-UO-YC","RC-AN-UO-YC","RC-AC-YO-YC",
                         "RC-GN-GO-YO","GN-GO-YO-UO","GO-YO-UO-YC","RC-GO-UO-YC"};
  string BPTYPE[]={"GC","CG","AU","UA","GU","UG","cGcC","cCcG","cAcU","cUcA","cGcU","cUcG"};
  for(size_t ii=0;ii<Vbps.size();ii++) {
    vector<string> stoken=split_string(Vbps[ii].residue_i,"/");
    if(stoken.size()!=3) {
      string str="RNA_coarse_grained.cpp: Wrong residue name(chainid/resid/resname) --> "+Vbps[ii].residue_i;
      error->one(str);
    }
    string resname_i=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    Residue residue_i=rna_cg.find_residue(resname_i);
    string bp_type=stoken[2];
    stoken=split_string(Vbps[ii].residue_j,"/");
    if(stoken.size()!=3) {
      string str="RNA_coarse_grained.cpp: Wrong residue name(chainid/resid/resname) --> "+Vbps[ii].residue_j;
      error->one(str);
    }
    string resname_j=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    Residue residue_j=rna_cg.find_residue(resname_j);
    bp_type+=stoken[2];
    if(count(begin(BPTYPE),end(BPTYPE),bp_type)==0) {
      string str="RNA_coarse_grained.cpp: Skip noncanonical base pair --> "+bp_type;
      error->warning(str);
      continue;
    }
    int bp_i_1=0,bp_i_2=0,bp_i_3=0,bp_j_1=0,bp_j_2=0,bp_j_3=0;
    int flag_GU=0;
    if(bp_type=="GC" || bp_type=="cGcC") {
      bp_i_1=residue_i["GO"].id; 
      bp_i_2=residue_i["GN"].id;
      bp_i_3=residue_i["RC"].id;
      bp_j_1=residue_j["CN"].id;
      bp_j_2=residue_j["YO"].id;
      bp_j_3=residue_j["YC"].id;
    }
    else if(bp_type=="CG" || bp_type=="cCcG") {
      bp_i_1=residue_j["GO"].id; 
      bp_i_2=residue_j["GN"].id;
      bp_i_3=residue_j["RC"].id;
      bp_j_1=residue_i["CN"].id;
      bp_j_2=residue_i["YO"].id;
      bp_j_3=residue_i["YC"].id;
    }
    else if(bp_type=="AU" || bp_type=="cAcU") {
      bp_i_1=residue_i["AN"].id; 
      bp_i_2=residue_i["AC"].id;
      bp_i_3=residue_i["RC"].id;
      bp_j_1=residue_j["UO"].id;
      bp_j_2=residue_j["YO"].id;
      bp_j_3=residue_j["YC"].id;
    }
    else if(bp_type=="UA" || bp_type=="cUcA") {
      bp_i_1=residue_j["AN"].id; 
      bp_i_2=residue_j["AC"].id;
      bp_i_3=residue_j["RC"].id;
      bp_j_1=residue_i["UO"].id;
      bp_j_2=residue_i["YO"].id;
      bp_j_3=residue_i["YC"].id;
    }
    else if(bp_type=="GU" || bp_type=="cGcU") {
      bp_i_1=residue_i["GO"].id; 
      bp_i_2=residue_i["GN"].id;
      bp_i_3=residue_i["RC"].id;
      bp_j_1=residue_j["UO"].id;
      bp_j_2=residue_j["YO"].id;
      bp_j_3=residue_j["YC"].id;
    }
    else if(bp_type=="UG" || bp_type=="cUcG") {
      bp_i_1=residue_j["GO"].id; 
      bp_i_2=residue_j["GN"].id;
      bp_i_3=residue_j["RC"].id;
      bp_j_1=residue_i["UO"].id;
      bp_j_2=residue_i["YO"].id;
      bp_j_3=residue_i["YC"].id;
    }
    
    Bond bond_one;
    bond_one.atom_i=bp_i_1;
    bond_one.atom_j=bp_j_1;
    bond_one.type=atomid_name[bp_i_1]+"-"+atomid_name[bp_j_1];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+Vbps[ii].residue_i+" "+Vbps[ii].residue_j;
      str+=" "+to_string(bp_i_1)+" "+to_string(bp_j_1)+" "+bond_one.type;
      error->one(str);
    }
    if(type_app.size()>0)
      bond_one.type+="_"+type_app;
    bonds.push_back(bond_one);
    
    if(flag_GU==0)
      bond_one.atom_i=bp_i_2;
    else
      bond_one.atom_i=bp_i_1;
    bond_one.atom_j=bp_j_2;
    bond_one.type=atomid_name[bond_one.atom_i]+"-"+atomid_name[bond_one.atom_j];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+Vbps[ii].residue_i+" "+Vbps[ii].residue_j;
      str+=" "+to_string(bond_one.atom_i)+" "+to_string(bond_one.atom_j)+" "+bond_one.type;
      error->one(str);
    }
    if(type_app.size()>0)
      bond_one.type+="_"+type_app;
    bonds.push_back(bond_one);

    if(residue_i.check_atom("S") && residue_j.check_atom("S")) {
      bond_one.atom_i=residue_i["S"].id;
      bond_one.atom_j=residue_j["S"].id;
      bond_one.type=atomid_name[bond_one.atom_i]+"-"+atomid_name[bond_one.atom_j];
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+Vbps[ii].residue_i+" "+Vbps[ii].residue_j;
        str+=" "+to_string(bond_one.atom_i)+" "+to_string(bond_one.atom_j)+" "+bond_one.type;
        error->one(str);
      }
      if(type_app.size()>0)
        bond_one.type+="_"+type_app;
      bonds.push_back(bond_one);
    }

    if(residue_i.check_atom("P") && residue_j.check_atom("P")) {
      bond_one.atom_i=residue_i["P"].id;
      bond_one.atom_j=residue_j["P"].id;
      bond_one.type=atomid_name[bond_one.atom_i]+"-"+atomid_name[bond_one.atom_j];
      if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown bond type --> "+Vbps[ii].residue_i+" "+Vbps[ii].residue_j;
        str+=" "+to_string(bond_one.atom_i)+" "+to_string(bond_one.atom_j)+" "+bond_one.type;
        error->one(str);
      }
      if(type_app.size()>0)
        bond_one.type+="_"+type_app;
      bonds.push_back(bond_one);
    }
    
    Angle angle_one;
    if(flag_GU==0) {
      angle_one.atom_i=bp_i_2;
      angle_one.atom_j=bp_j_2;
      angle_one.atom_k=bp_j_1;
    }
    else {
      angle_one.atom_i=bp_i_3;
      angle_one.atom_j=bp_i_1;
      angle_one.atom_k=bp_j_1;
    }
    angle_one.type=atomid_name[angle_one.atom_i]+"-"+atomid_name[angle_one.atom_j]+"-"+atomid_name[angle_one.atom_k];
    if(count(begin(ANGLETYPE),end(ANGLETYPE),angle_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown angle type --> "+Vbps[ii].residue_i+" "+Vbps[ii].residue_j;
      str+=" "+to_string(angle_one.atom_i)+" "+to_string(angle_one.atom_j)+" "+to_string(angle_one.atom_k);
      str+=" "+angle_one.type;
      error->one(str);
    }
    angles.push_back(angle_one);

    vector<Dihedral> Tbpdihedral;
    Dihedral dihedral_one;
    dihedral_one.type="";
    if(flag_GU==0) {
      dihedral_one.atom_i=bp_i_3;
      dihedral_one.atom_j=bp_i_1;
      dihedral_one.atom_k=bp_i_2;
      dihedral_one.atom_l=bp_j_2;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_1;
      dihedral_one.atom_j=bp_i_2;
      dihedral_one.atom_k=bp_j_2;
      dihedral_one.atom_l=bp_j_1;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_2;
      dihedral_one.atom_j=bp_j_2;
      dihedral_one.atom_k=bp_j_1;
      dihedral_one.atom_l=bp_j_3;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_3;
      dihedral_one.atom_j=bp_i_1;
      dihedral_one.atom_k=bp_j_1;
      dihedral_one.atom_l=bp_j_3;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_3;
      dihedral_one.atom_j=bp_i_2;
      dihedral_one.atom_k=bp_j_2;
      dihedral_one.atom_l=bp_j_3;
      Tbpdihedral.push_back(dihedral_one);
    }
    else {
      dihedral_one.atom_i=bp_i_3;
      dihedral_one.atom_j=bp_i_2;
      dihedral_one.atom_k=bp_i_1;
      dihedral_one.atom_l=bp_j_2;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_2;
      dihedral_one.atom_j=bp_i_1;
      dihedral_one.atom_k=bp_j_2;
      dihedral_one.atom_l=bp_j_1;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_1;
      dihedral_one.atom_j=bp_j_2;
      dihedral_one.atom_k=bp_j_1;
      dihedral_one.atom_l=bp_j_3;
      Tbpdihedral.push_back(dihedral_one);

      dihedral_one.atom_i=bp_i_3;
      dihedral_one.atom_j=bp_i_1;
      dihedral_one.atom_k=bp_j_1;
      dihedral_one.atom_l=bp_j_3;
      Tbpdihedral.push_back(dihedral_one);
    }

    for(size_t kk=0;kk<Tbpdihedral.size();kk++) {
      dihedral_one=Tbpdihedral[kk];
      dihedral_one.type=atomid_name[dihedral_one.atom_i]+"-"+atomid_name[dihedral_one.atom_j]+"-";
      dihedral_one.type+=atomid_name[dihedral_one.atom_k]+"-"+atomid_name[dihedral_one.atom_l];
      if(count(begin(DIHEDRALTYPE),end(DIHEDRALTYPE),dihedral_one.type)==0) {
        string str="RNA_coarse_grained.cpp: Unknown dihedral type --> "+Vbps[ii].residue_i+" "+Vbps[ii].residue_j;
        str+=" "+to_string(dihedral_one.atom_i)+" "+to_string(dihedral_one.atom_j)+" ";
        str+=to_string(dihedral_one.atom_k)+" "+to_string(dihedral_one.atom_l)+" "+dihedral_one.type;
        error->one(str);
      }
      dihedrals.push_back(dihedral_one);
    }
  }
}

void RNACoarseGrained::add_RNA_basestacks(const vector<string> &Vbss) {
  string BONDTYPE[]={"AN-AN","AC-AC","GO-GO","GN-GN","AN-GO","AC-GN,","GO-CN","GN-YO","AN-UO","AC-YO",
                     "GO-UO","AN-CN","CN-CN","YO-YO","CN-UO","UO-UO"};
  string ATOMTYPE1[]={"AN","GO","CN","UO"};
  string ATOMTYPE2[]={"AC","GN","YO"};
  string BS_INV[]={"GA","cGcA","CA","cCcA","UA","cUcA","CG","cCcG","UG","cUcG","UC","cUcC"};
  for(size_t ii=0;ii<Vbss.size();ii++) {
    vector<string> stoken=split_string(Vbss[ii]);
    if(stoken.size()!=2) {
      string str="RNA_coarse_grained.cpp: Wrong base stacking format(residue_i residue_j) --> "+Vbss[ii];
      error->one(str);
    }
    string resname_i=stoken[0];
    string resname_j=stoken[1];
    stoken.clear();
    stoken=split_string(resname_i,"/");
    if(stoken.size()!=3) {
      string str="RNA_coarse_grained.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_i;
      error->one(str);
    }
    resname_i=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    Residue residue_i=rna_cg.find_residue(resname_i);
    string bs_type=stoken[2];
    stoken.clear();
    stoken=split_string(resname_j,"/");
    if(stoken.size()!=3) {
      string str="RNA_coarse_grained.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_j;
      error->one(str);
    }
    resname_j=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    Residue residue_j=rna_cg.find_residue(resname_j);
    bs_type+=stoken[2];

    int bs_i_1,bs_i_2,bs_j_1,bs_j_2;
    for(size_t kk=0;kk<residue_i.atoms.size();kk++) {
      string atomname_one=residue_i.atoms[kk].atomname;
      if(count(begin(ATOMTYPE1),end(ATOMTYPE1),atomname_one))
        bs_i_1=residue_i.atoms[kk].id;
      else if(count(begin(ATOMTYPE2),end(ATOMTYPE2),atomname_one))
        bs_i_2=residue_i.atoms[kk].id;
    }
    for(size_t kk=0;kk<residue_j.atoms.size();kk++) {
      string atomname_one=residue_j.atoms[kk].atomname;
      if(count(begin(ATOMTYPE1),end(ATOMTYPE1),atomname_one))
        bs_j_1=residue_j.atoms[kk].id;
      else if(count(begin(ATOMTYPE2),end(ATOMTYPE2),atomname_one))
        bs_j_2=residue_j.atoms[kk].id;
    }
    if(count(begin(BS_INV),end(BS_INV),bs_type)) {
      swap(bs_i_1,bs_j_1);
      swap(bs_i_2,bs_j_2);
    }
    
    Bond bond_one;
    bond_one.atom_i=bs_i_1;
    bond_one.atom_j=bs_j_1;
    bond_one.type=atomid_name[bs_i_1]+"-"+atomid_name[bs_j_1];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+bond_one.type;
      error->one(str);
    }
    bond_one.type+="_bs";
    bonds.push_back(bond_one);
    
    bond_one.atom_i=bs_i_2;
    bond_one.atom_j=bs_j_2;
    bond_one.type=atomid_name[bs_i_2]+"-"+atomid_name[bs_j_2];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+bond_one.type;
      error->one(str);
    }
    bond_one.type+="_bs";
    bonds.push_back(bond_one);
  }
}

void RNACoarseGrained::add_RNA_basestacks(const vector<ResiduePair> &Vbss) {
  string BONDTYPE[]={"AN-AN","AC-AC","GO-GO","GN-GN","AN-GO","AC-GN,","GO-CN","GN-YO","AN-UO","AC-YO",
                     "GO-UO","AN-CN","CN-CN","YO-YO","CN-UO","UO-UO"};
  string ATOMTYPE1[]={"AN","GO","CN","UO"};
  string ATOMTYPE2[]={"AC","GN","YO"};
  string BS_INV[]={"GA","cGcA","CA","cCcA","UA","cUcA","CG","cCcG","UG","cUcG","UC","cUcC"};
  for(size_t ii=0;ii<Vbss.size();ii++) {
    string resname_i=Vbss[ii].residue_i;
    string resname_j=Vbss[ii].residue_j;
    vector<string> stoken=split_string(resname_i,"/");
    if(stoken.size()!=3) {
      string str="RNA_coarse_grained.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_i;
      error->one(str);
    }
    resname_i=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    Residue residue_i=rna_cg.find_residue(resname_i);
    string bs_type=stoken[2];
    stoken.clear();
    stoken=split_string(resname_j,"/");
    if(stoken.size()!=3) {
      string str="RNA_coarse_grained.cpp: Wrong residue name(chainid/resid/resname) --> "+resname_j;
      error->one(str);
    }
    resname_j=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    Residue residue_j=rna_cg.find_residue(resname_j);
    bs_type+=stoken[2];

    int bs_i_1,bs_i_2,bs_j_1,bs_j_2;
    for(size_t kk=0;kk<residue_i.atoms.size();kk++) {
      string atomname_one=residue_i.atoms[kk].atomname;
      if(count(begin(ATOMTYPE1),end(ATOMTYPE1),atomname_one))
        bs_i_1=residue_i.atoms[kk].id;
      else if(count(begin(ATOMTYPE2),end(ATOMTYPE2),atomname_one))
        bs_i_2=residue_i.atoms[kk].id;
    }
    for(size_t kk=0;kk<residue_j.atoms.size();kk++) {
      string atomname_one=residue_j.atoms[kk].atomname;
      if(count(begin(ATOMTYPE1),end(ATOMTYPE1),atomname_one))
        bs_j_1=residue_j.atoms[kk].id;
      else if(count(begin(ATOMTYPE2),end(ATOMTYPE2),atomname_one))
        bs_j_2=residue_j.atoms[kk].id;
    }
    if(count(begin(BS_INV),end(BS_INV),bs_type)) {
      swap(bs_i_1,bs_j_1);
      swap(bs_i_2,bs_j_2);
    }
    
    Bond bond_one;
    bond_one.atom_i=bs_i_1;
    bond_one.atom_j=bs_j_1;
    bond_one.type=atomid_name[bs_i_1]+"-"+atomid_name[bs_j_1];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+bond_one.type;
      error->one(str);
    }
    bond_one.type+="_bs";
    bonds.push_back(bond_one);
    
    bond_one.atom_i=bs_i_2;
    bond_one.atom_j=bs_j_2;
    bond_one.type=atomid_name[bs_i_2]+"-"+atomid_name[bs_j_2];
    if(count(begin(BONDTYPE),end(BONDTYPE),bond_one.type)==0) {
      string str="RNA_coarse_grained.cpp: Unknown bond type --> "+bond_one.type;
      error->one(str);
    }
    bond_one.type+="_bs";
    bonds.push_back(bond_one);
  }
}

void RNACoarseGrained::add_distance_restraint(string atom_i, string atom_j, string type_one) {
  vector<string> stoken=split_string(atom_i,"/");
  if(stoken.size()!=4) {
    string str="RNA_coarse_grained.cpp: Wrong atom name(chainid/resid/resname/atomname) --> "+atom_i;
    error->one(str);
  }
  string resname_i=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
  string atomname_i=stoken[3];

  stoken=split_string(atom_j,"/");
  if(stoken.size()!=4) {
    string str="RNA_coarse_grained.cpp: Wrong atom name(chainid/resid/resname/atomname) --> "+atom_j;
    error->one(str);
  }
  string resname_j=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
  string atomname_j=stoken[3];

  Residue residue_i=rna_cg.find_residue(resname_i);
  Residue residue_j=rna_cg.find_residue(resname_j);
  Bond bond_one;
  bond_one.atom_i=residue_i[atomname_i].id;
  bond_one.atom_j=residue_j[atomname_j].id;
  bond_one.type=type_one;
  bonds.push_back(bond_one);
}

void RNACoarseGrained::add_distance_restraint(DistanceRestraint &disres_one) {
  vector<string> stoken=split_string(disres_one.atom_i,"/");
  if(stoken.size()!=4) {
    string str="RNA_coarse_grained.cpp: Wrong atom name(chainid/resid/resname/atomname) --> "+disres_one.atom_i;
    error->one(str);
  }
  string resname_i=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
  string atomname_i=stoken[3];

  stoken=split_string(disres_one.atom_j,"/");
  if(stoken.size()!=4) {
    string str="RNA_coarse_grained.cpp: Wrong atom name(chainid/resid/resname/atomname) --> "+disres_one.atom_j;
    error->one(str);
  }
  string resname_j=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
  string atomname_j=stoken[3];

  Residue residue_i=rna_cg.find_residue(resname_i);
  Residue residue_j=rna_cg.find_residue(resname_j);
  Bond bond_one;
  bond_one.atom_i=residue_i[atomname_i].id;
  bond_one.atom_j=residue_j[atomname_j].id;
  bond_one.type=disres_one.type+" "+disres_one.coeff;
  bonds.push_back(bond_one);
}

void RNACoarseGrained::add_distance_restraint(const vector<DistanceRestraint> &Vdisress) {
  for(size_t ii=0;ii<Vdisress.size();ii++) {
    DistanceRestraint disres_one=Vdisress[ii];
    vector<string> stoken=split_string(disres_one.atom_i,"/");
    if(stoken.size()!=4) {
      string str="RNA_coarse_grained.cpp: Wrong atom name(chainid/resid/resname/atomname) --> "+disres_one.atom_i;
      error->one(str);
    }
    string resname_i=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    string atomname_i=stoken[3];

    stoken=split_string(disres_one.atom_j,"/");
    if(stoken.size()!=4) {
      string str="RNA_coarse_grained.cpp: Wrong atom name(chainid/resid/resname/atomname) --> "+disres_one.atom_j;
      error->one(str);
    }
    string resname_j=stoken[0]+"/"+stoken[1]+"/c"+stoken[2];
    string atomname_j=stoken[3];

    Residue residue_i=rna_cg.find_residue(resname_i);
    Residue residue_j=rna_cg.find_residue(resname_j);
    Bond bond_one;
    bond_one.atom_i=residue_i[atomname_i].id;
    bond_one.atom_j=residue_j[atomname_j].id;
    bond_one.type=disres_one.type+" "+disres_one.coeff;
    bonds.push_back(bond_one);
  }
}

size_t RNACoarseGrained::number_of_CG_atoms() const {
  /*
  size_t num=0;
  for(size_t ii=0;ii<rna_cg.chains.size();ii++) {
    for(size_t jj=0;jj<rna_cg.chains[ii].size();jj++) {
      num+=rna_cg.chains[ii].residues[jj].atoms.size();
    }
  }

  return num;
  */
  return atomid_name.size();
}

void RNACoarseGrained::clear() {
  rna_cg.clear();
  bonds.clear();
  bonds.shrink_to_fit();
  angles.clear();
  angles.shrink_to_fit();
  dihedrals.clear();
  dihedrals.shrink_to_fit();
  impropers.clear();
  impropers.shrink_to_fit();
  atomid_name.clear();
}