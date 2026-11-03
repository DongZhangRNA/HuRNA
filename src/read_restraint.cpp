#include "basic.h"
#include "readfile.h"
#include "initialization.h"
#include "error.h"
#include "measure.h"
#include "read_PDB.h"
#include "RNA_molecule.h"
#include "read_restraint.h"

using namespace CBS_NS;

string BONDTYPE[]={"harmonic","harmonic/gauss","slope","morse"};

ReadRestraint::ReadRestraint(Cobmms *cbs) : Pointers(cbs) {}

ReadRestraint::ReadRestraint(const ReadRestraint &other) : Pointers(other.cbs) {
  RNA_templates=other.RNA_templates;
  RNA_basepairs=other.RNA_basepairs;
  RNA_basestacks=other.RNA_basestacks;
  rigid_bodies=other.rigid_bodies;
  distance_restraints=other.distance_restraints;
}

ReadRestraint& ReadRestraint::operator=(const ReadRestraint &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    RNA_templates=other.RNA_templates;
    RNA_basepairs=other.RNA_basepairs;
    RNA_basestacks=other.RNA_basestacks;
    rigid_bodies=other.rigid_bodies;
    distance_restraints=other.distance_restraints;
  }

  return *this;
}

ReadRestraint::~ReadRestraint() {
  RNA_templates.clear();
  RNA_templates.shrink_to_fit();
  RNA_basepairs.clear();
  RNA_basepairs.shrink_to_fit();
  RNA_basestacks.clear();
  RNA_basestacks.shrink_to_fit();
  rigid_bodies.clear();
  rigid_bodies.shrink_to_fit();
  distance_restraints.clear();
  distance_restraints.shrink_to_fit();
}

void ReadRestraint::read_restraint(string fstr) {
  vector<string> Vlines=read_file(fstr);
  vector<string> restraint_types;
  vector<vector<string> > restraint_infos;
  string restr_type="";
  vector<string> Vtmp;
  for(size_t ii=0;ii<Vlines.size();ii++) {
    string sline=Vlines[ii];
    size_t found=sline.find_first_of("#");
    if(found!=string::npos) {
      sline=sline.substr(0,found);
    }
    if(sline.size()<2) continue;
    if(sline.find("Restraint:")!=string::npos) {
      if(restr_type!="") {
        restraint_types.push_back(restr_type);
        restraint_infos.push_back(Vtmp);
      }
      vector<string> stoken=split_string(sline,":");
      restr_type=stoken[1];
      Vtmp.clear();
      continue;
    }
    Vtmp.push_back(sline);
  }
  if(restr_type!="") {
    restraint_types.push_back(restr_type);
    restraint_infos.push_back(Vtmp);
  }

  RNA_templates.clear();
  RNA_basepairs.clear();
  RNA_basestacks.clear();
  rigid_bodies.clear();
  distance_restraints.clear();

  for(size_t ii=0;ii<restraint_types.size();ii++) {
    if(restraint_types[ii]=="RNA_template")
      ReadRestraint::add_RNA_template(restraint_infos[ii]);
    else if(restraint_types[ii]=="RNA_basepair") {
      for(size_t jj=0;jj<restraint_infos[ii].size();jj++)
        ReadRestraint::add_RNA_basepair(restraint_infos[ii][jj]);
    }
    else if(restraint_types[ii]=="RNA_basestack") {
      for(size_t jj=0;jj<restraint_infos[ii].size();jj++)
        ReadRestraint::add_RNA_basestack(restraint_infos[ii][jj]);
    }
    else if(restraint_types[ii]=="Rigid_body") {
      for(size_t jj=0;jj<restraint_infos[ii].size();jj++)
        ReadRestraint::add_rigid_body(restraint_infos[ii][jj]);
    }
    else if(restraint_types[ii]=="Distance_restraint") {
      for(size_t jj=0;jj<restraint_infos[ii].size();jj++)
        ReadRestraint::add_distance_restraint(restraint_infos[ii][jj]);
    }
    else {
      string str="read_restraint.cpp: Unknown restraint type(RNA_template/Rigid_body/RNA_basepair/RNA_basestack/Distance_restraint) --> ";
      str+=restraint_types[ii];
      error->one(str);
    }
  }
}

void ReadRestraint::add_RNA_template(const vector<string> &Vstr) {
  ReadPDB readpdb_one(cbs);
  readpdb_one.parse_PDB(Vstr);
  RNAMolecule rna_one(cbs);
  rna_one.set_RNA_molecule(readpdb_one.RNA_atoms[0]);
  vector<Chain> Vchains;
  for(size_t ii=0;ii<rna_one.chains.size();ii++) {
    Chain chain_one(cbs);
    chain_one.chain_label=Chain::rna;
    for(size_t jj=0;jj<rna_one.chains[ii].residues.size();jj++) {
      Residue& residue_i=rna_one.chains[ii].residues[jj];
      if(jj==0) {
        chain_one.insert_residue(residue_i,Chain::end);
        continue;
      }
      if(residue_i.chainid!=chain_one.residues.back().chainid || residue_i.resid-chain_one.residues.back().resid!=1) {
        Vchains.push_back(chain_one);
        chain_one.clear();
        chain_one.chain_label=Chain::rna;
      }
      chain_one.insert_residue(residue_i,Chain::end);
    }
    if(chain_one.size()>0)
      Vchains.push_back(chain_one);
  }

  RNAMotif motif_one(cbs);
  motif_one.type="user_template";
  for(size_t ii=0;ii<Vchains.size();ii++) {
    motif_one.insert_RNA_segment(Vchains[ii],0);
    string secstr_one=".";
    motif_one.secstrs.push_back(secstr_one);
  }

  for(size_t ii=0;ii<RNA_templates.size();ii++) {
    int check=Measure::vector_relation(motif_one.Vresnames,RNA_templates[ii].Vresnames);
    if(check!=Measure::separated) {
      string str="read_restraint.cpp: There are overlapped residues between "+RNA_templates[ii].get_RNA_motif_info();
      str+=" and "+motif_one.get_RNA_motif_info();
      error->one(str);
    }
  }

  if(init->user_template_rigid)
    motif_one.flag_rigid=1;
  
  RNA_templates.push_back(motif_one);

  if(init->user_template_rigid) {
    vector<ResiduePair> nucleotide_pair;
    for(size_t ii=0;ii<motif_one.ends.size()/2;ii++) {
      ResiduePair nt_pair_one;
      nt_pair_one.residue_i=motif_one.ends[ii*2];
      nt_pair_one.residue_j=motif_one.ends[ii*2+1];
      nucleotide_pair.push_back(nt_pair_one);
    }
    rigid_bodies.push_back(nucleotide_pair);
  }
}

void ReadRestraint::add_rigid_body(string str_one) {
  vector<string> Vnt_pairs;
  for(size_t ii=0;ii<rigid_bodies.size();ii++) {
    for(size_t jj=0;jj<rigid_bodies[ii].size();jj++) {
      string spair=rigid_bodies[ii][jj].residue_i+"-"+rigid_bodies[ii][jj].residue_j;
      Vnt_pairs.push_back(spair);
    }
  }
  vector<ResiduePair> nucleotide_pair;
  vector<string> resids=split_string(str_one,",;");
  for(size_t ii=0;ii<resids.size();ii++) {
    vector<string> stoken=split_string(resids[ii],"-");
    string nt_i,nt_j;
    if(stoken.size()==1) {
      nt_i=stoken[0];
      nt_j=stoken[0];
    }
    else if(stoken.size()==2) {
      nt_i=stoken[0];
      nt_j=stoken[1];
    }
    else {
      string str="read_restraint.cpp: Wrong input format in add_rigid_body --> "+str_one;
      error->one(str);
    }
    stoken=split_string(nt_i,"/");
    if(stoken.size()!=3) {
      string str="read_restraint.cpp: Wrong residue name(chainid/resid/resname) --> "+nt_i;
      error->one(str);
    }
    stoken=split_string(nt_j,"/");
    if(stoken.size()!=3) {
      string str="read_restraint.cpp: Wrong residue name(chainid/resid/resname) --> "+nt_j;
      error->one(str);
    }

    string spair=nt_i+"-"+nt_j;
    if(count(Vnt_pairs.begin(),Vnt_pairs.end(),spair)==0) {
      ResiduePair nt_pair_one;
      nt_pair_one.residue_i=nt_i;
      nt_pair_one.residue_j=nt_j;
      nucleotide_pair.push_back(nt_pair_one);
    }
  }
  if(nucleotide_pair.size()>0)
    rigid_bodies.push_back(nucleotide_pair);
}

void ReadRestraint::add_RNA_basepair(string str_one) {
  vector<string> resids=split_string(str_one,"-");
  if(resids.size()!=2) {
    string str="read_restraint.cpp: Wrong input format in add_RNA_basepair --> "+str_one;
    error->one(str);
  }
  vector<string> stoken=split_string(resids[0],"/");
  if(stoken.size()!=3) {
    string str="read_restraint.cpp: Wrong residue name(chainid/resid/resname) --> "+resids[0];
    error->one(str);
  }
  stoken=split_string(resids[1],"/");
  if(stoken.size()!=3) {
    string str="read_restraint.cpp: Wrong residue name(chainid/resid/resname) --> "+resids[1];
    error->one(str);
  }
  
  ResiduePair nt_pair_one;
  nt_pair_one.residue_i=resids[0];
  nt_pair_one.residue_j=resids[1];
  RNA_basepairs.push_back(nt_pair_one);
}

void ReadRestraint::add_RNA_basestack(string str_one) {
  vector<string> resids=split_string(str_one,"-");
  if(resids.size()!=2) {
    string str="read_restraint.cpp: Wrong input format in add_RNA_basestack --> "+str_one;
    error->one(str);
  }
  vector<string> stoken=split_string(resids[0],"/");
  if(stoken.size()!=3) {
    string str="read_restraint.cpp: Wrong residue name(chainid/resid/resname) --> "+resids[0];
    error->one(str);
  }
  stoken=split_string(resids[1],"/");
  if(stoken.size()!=3) {
    string str="read_restraint.cpp: Wrong residue name(chainid/resid/resname) --> "+resids[1];
    error->one(str);
  }
  
  ResiduePair nt_pair_one;
  nt_pair_one.residue_i=resids[0];
  nt_pair_one.residue_j=resids[1];
  RNA_basestacks.push_back(nt_pair_one);
}

void ReadRestraint::add_distance_restraint(string str_one) {
  vector<string> restr_info=split_string(str_one);
  if(restr_info.size()<4) {
    string str="read_restraint.cpp: Wrong input format in add_distance_restraint --> "+str_one;
    error->one(str);
  }
  string atom_i=restr_info[0];
  string atom_j=restr_info[1];
  string bondtype=restr_info[2];
  string bondcoeff=restr_info[3];
  for(size_t ii=4;ii<restr_info.size();ii++) 
    bondcoeff+=" "+restr_info[ii];

  vector<string> stoken=split_string(atom_i,"/");
  if(stoken.size()!=4) {
    string str="read_restraint.cpp: Wrong atom name(chainid/resid/resname/atomname) --> "+atom_i;
    error->one(str);
  }
  stoken=split_string(atom_j,"/");
  if(stoken.size()!=4) {
    string str="read_restraint.cpp: Wrong atom name(chainid/resid/resname/atomname) --> "+atom_j;
    error->one(str);
  }
  if(count(begin(BONDTYPE),end(BONDTYPE),bondtype)==0) {
    string str="read_restraint.cpp: Unknown restraint type in add_distance_restraint --> "+str_one;
    error->one(str);
  }
  stoken=split_string(bondcoeff);
  if(bondtype=="harmonic" && stoken.size()!=2) {
    string str="read_restraint.cpp: Wrong input format in add_distance_restraint --> "+str_one;
    error->one(str);
  }
  else if(bondtype=="harmonic/gauss" && stoken.size()!=5) {
    string str="read_restraint.cpp: Wrong input format in add_distance_restraint --> "+str_one;
    error->one(str);
  }
  else if(bondtype=="slope" && stoken.size()!=4) {
    string str="read_restraint.cpp: Wrong input format in add_distance_restraint --> "+str_one;
    error->one(str);
  }
  else if(bondtype=="morse" && stoken.size()!=3) {
    string str="read_restraint.cpp: Wrong input format in add_distance_restraint --> "+str_one;
    error->one(str);
  }

  DistanceRestraint disrestr_one;
  disrestr_one.atom_i=atom_i;
  disrestr_one.atom_j=atom_j;
  disrestr_one.type=bondtype;
  disrestr_one.coeff=bondcoeff;
  distance_restraints.push_back(disrestr_one);
}

void ReadRestraint::clear() {
  RNA_templates.clear();
  RNA_templates.shrink_to_fit();
  RNA_basepairs.clear();
  RNA_basepairs.shrink_to_fit();
  RNA_basestacks.clear();
  RNA_basestacks.shrink_to_fit();
  rigid_bodies.clear();
  rigid_bodies.shrink_to_fit();
  distance_restraints.clear();
  distance_restraints.shrink_to_fit();
}

ostream& operator <<(ostream &output,ReadRestraint &restr_one) {
  if(restr_one.RNA_templates.size()>0) {
    output<<"#There are "<<restr_one.RNA_templates.size()<<" RNA_templates :"<<endl;
    for(size_t ii=0;ii<restr_one.RNA_templates.size();ii++) {
      output<<"Restraint: RNA_template ";
      vector<string> ends=restr_one.RNA_templates[ii].ends;
      string resids=ends[0]+"-"+ends[1];
      for(size_t jj=1;jj<ends.size()/2;jj++)
        resids+=", "+ends[jj*2]+"-"+ends[jj*2+1];
      output<<"  "<<resids<<endl;
      output<<restr_one.RNA_templates[ii];
      output<<endl;
    }
    output<<endl;
  }
  if(restr_one.rigid_bodies.size()>0) {
    output<<"#There are "<<restr_one.rigid_bodies.size()<<" Rigid_bodies :"<<endl;
    output<<"Restraint: Rigid_body "<<endl;
    for(size_t ii=0;ii<restr_one.rigid_bodies.size();ii++) {
      string resids=restr_one.rigid_bodies[ii][0].residue_i+"-"+restr_one.rigid_bodies[ii][0].residue_j;
      for(size_t jj=1;jj<restr_one.rigid_bodies[ii].size();jj++)
        resids+=", "+restr_one.rigid_bodies[ii][jj].residue_i+"-"+restr_one.rigid_bodies[ii][jj].residue_j;
      output<<"  "<<resids<<endl;
    }
    output<<endl;
  }
  if(restr_one.RNA_basepairs.size()>0) {
    output<<"#There are "<<restr_one.RNA_basepairs.size()<<" RNA_basepairs :"<<endl;
    output<<"Restraint: RNA_basepair "<<endl;
    for(size_t ii=0;ii<restr_one.RNA_basepairs.size();ii++)
      output<<"  "<<restr_one.RNA_basepairs[ii].residue_i<<" "<<restr_one.RNA_basepairs[ii].residue_j<<endl;
    output<<endl;
  }
  if(restr_one.RNA_basestacks.size()>0) {
    output<<"#There are "<<restr_one.RNA_basestacks.size()<<" RNA_basestacks :"<<endl;
    output<<"Restraint: RNA_basestack "<<endl;
    for(size_t ii=0;ii<restr_one.RNA_basestacks.size();ii++)
      output<<"  "<<restr_one.RNA_basestacks[ii].residue_i<<" "<<restr_one.RNA_basestacks[ii].residue_j<<endl;
    output<<endl;
  }
  if(restr_one.distance_restraints.size()>0) {
    output<<"#There are "<<restr_one.distance_restraints.size()<<" Distance_restraints :"<<endl;
    output<<"Restraint: Distance_restraint "<<endl;
    for(size_t ii=0;ii<restr_one.distance_restraints.size();ii++) {
      output<<"  "<<restr_one.distance_restraints[ii].atom_i<<" "<<restr_one.distance_restraints[ii].atom_j<<" ";
      output<<restr_one.distance_restraints[ii].type<<" "<<restr_one.distance_restraints[ii].coeff<<endl;
    }
    output<<endl;
  }
 
  return output;
}