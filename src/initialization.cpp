#include "basic.h"
#include "readfile.h"
#include "initialization.h"
#include "error.h"

using namespace CBS_NS;

Initialization::Initialization(Cobmms *ptr) : Pointers(ptr) {
  HOME_DIR="PATH_to_Cobmms/";
}

Initialization::~Initialization() {
  Excluded_PDBs.clear();
  Excluded_PDBs.shrink_to_fit();
  CircRNA_chains.clear();
  CircRNA_chains.shrink_to_fit();

  RESIDUES_XYZ.clear();
  RESIDUES_ATOMS.clear();
  NUCLEOTIDES_LIST.clear();
  RESIDUES_LIST.clear();
  RNA_STANDARD_NUCLEOTIDES.clear();
  RNA_STANDARD_NUCLEOTIDES.shrink_to_fit();
  PROTEIN_STANDARD_NUCLEOTIDES.clear();
  PROTEIN_STANDARD_NUCLEOTIDES.shrink_to_fit();

  SEQUENCE_MAP.clear();

  CHAIN_ID_LIST.clear();
  CHAIN_ID_LIST.shrink_to_fit();

  RNA_BP_NOTATION.clear();
  RNA_BP_NOTATION.shrink_to_fit();
}

void Initialization::set_parameters() {
  INPUT_DIR=HOME_DIR+"data/input/";
  RNA_TEMPLATE_DIR=HOME_DIR+"data/template/";
  DSSR=HOME_DIR+"external/DSSR/x3dna-dssr";
  MCAnnotate=HOME_DIR+"external/MC-Annotate/MC-Annotate";
  LMP=HOME_DIR+"external/lammps-29Aug2024/build/lmp";
  Contrafold=HOME_DIR+"external/contrafold/src/contrafold";
  USalign=HOME_DIR+"external/TMscore/USalign";
  MPIRUN_OPT="";
  RNA_DM="";

  WORK_DIR="./";

  noH=1;
  RNA_distance_O3p_P=1.602;
  RNA_distance_O3p_P_min=0.9;
  RNA_distance_O3p_P_max=3.1;
  RNA_check_chain_break=1;
  RNA_use_DM=0;
  RNA_BP_RMSD_cutoff=0.5;
  RNA_template_RMSD_cutoff=1.5;
  RNA_sequence_similarity_high=0.95;
  RNA_sequence_similarity_medium=0.5;
  RNA_template_nmax=5;
  Num_model_max=5;
  user_template_rigid=0;
  assembly_delta_distance=1.5;
  nthread=10;
  Num_remd_dup=10;
  output_format="std";
  RNA_version="latest";

  CHAIN_ID_LIST.clear();
  for(char ch='A';ch<='Z';ch++) {
    string cid=to_string(ch);
    CHAIN_ID_LIST.push_back(cid);
  }
  for(char ch='a';ch<='z';ch++) {
    string cid=to_string(ch);
    CHAIN_ID_LIST.push_back(cid);
  }
  for(char ch='0';ch<='9';ch++) {
    string cid=to_string(ch);
    CHAIN_ID_LIST.push_back(cid);
  }

  Initialization::set_residues_list();
  Initialization::set_standard_residues();
  Initialization::set_residue_name_map();
  Initialization::set_RNA_basepair_notation();
}

void Initialization::reset_parameters(string style, string arg) {
  if(style=="-format") output_format=arg;
  else if(style=="-noh") {
    noH=stoi(arg);
    Initialization::set_standard_residues();
  }
  else if(style=="-np") nthread=stoi(arg);
  else if(style=="-nout") Num_model_max=stoi(arg);
  else if(style=="-exclude") {
    vector<string> Vlines=read_file(arg);
    Excluded_PDBs.clear();
    for(size_t ii=0;ii<Vlines.size();ii++) {
      vector<string> stoken=split_string(Vlines[ii],",");
      for(size_t jj=0;jj<stoken.size();jj++) {
        string PDBid=stoken[jj];
        for(size_t kk=0;kk<PDBid.size();kk++)
          PDBid[kk]=tolower(PDBid[kk]);
        Excluded_PDBs.push_back(PDBid);
      }
    }
  }
  else if(style=="-circRNA") {
    CircRNA_chains.clear();
    vector<string> stoken=split_string(arg,",");
    for(size_t ii=0;ii<stoken.size();ii++)
      CircRNA_chains.push_back(stoken[ii]);
  }
  else if(style=="-workdir") {
    WORK_DIR=arg;
    if(WORK_DIR.back()!='/')
      WORK_DIR+="/";
    if(!is_directory_exist(WORK_DIR)) {
      string cmd="mkdir "+WORK_DIR;
      int sys=system(cmd.c_str());

      if(sys!=0) {
        string str="Fail to execute "+cmd;
        error->one(str);
      }
    }
  }
  else if(style=="-version") {
    RNA_version=arg;
    string VERSION[]={"v2","latest"};
    if(count(begin(VERSION),end(VERSION),RNA_version)==0) {
      string str="initialization.cpp: Wrong RNA version! --> "+RNA_version;
      error->one(str);
    }
  }
  else if(style=="-lmp") {
    LMP=arg;
    if(!is_file_exist(LMP)) {
      string str="initialization.cpp: Cannot find executable lmp! --> "+LMP;
      error->one(str);
    }
  }
  else if(style=="-mpi_opt") {
    MPIRUN_OPT=arg;
  }
}

vector<string> Initialization::read_water_ions_name() {
  vector<string> Vsol_names;
  string fstr=INPUT_DIR+"water_ions.dat";
  vector<string> Vlines=read_file(fstr);
  for(size_t ii=0;ii<Vlines.size();ii++) {
    vector<string> stoken=split_string(Vlines[ii],",");
    if(ii==0) {
      string snum=to_string(stoken.size());
      Vsol_names.push_back(snum);
    }
    for(size_t jj=0;jj<stoken.size();jj++) {
      string name_one=stoken[jj];
      for(size_t kk=0;kk<name_one.size();kk++)
        name_one[kk]=toupper(name_one[kk]);
      Vsol_names.push_back(name_one);
    }
  }

  return Vsol_names;
}

map<string,string> Initialization::read_mutation_atommap() {
  string fstr=INPUT_DIR+"mutation_atommap.dat";
  vector<string> Vlines=read_file(fstr);
  map<string,vector<string> > Amap;
  vector<string> Vtmp;
  string tname="";
  for(size_t ii=0;ii<Vlines.size();ii++) {
    if(Vlines[ii].size()<2) continue;
    if(Vlines[ii].find("#")!=string::npos) continue;
    if(Vlines[ii].find("[")!=string::npos && Vlines[ii].find("]")!=string::npos) {
      if(tname!="" && Vtmp.size()>0) {
        if(Amap.count(tname)==0) {
          Amap.emplace(tname,Vtmp);
          tname=Vlines[ii];
          Vtmp.clear();
        }
        else {
          string str="initialization.cpp: Mutant name has been used in mutation atommap --> "+Vlines[ii];
          error->one(str);
        }
      }
      else
        tname=Vlines[ii];
      continue;
    }
    else
      Vtmp.push_back(Vlines[ii]);
  }
  if(Amap.count(tname)==0)
    Amap.emplace(tname,Vtmp);
  else {
    string str="initialization.cpp: Mutant name has been used in mutation atommap --> "+tname;
    error->one(str);
  }
  map<string,string> Matommap;
  if(Amap.count("[CommonAtom]")==0) {
    string str="initialization.cpp: No common atoms in mutation atommap "+fstr;
    error->warning(str);
    string tmp="";
    Matommap.emplace("CommonAtom",tmp);
  }
  else {
    vector<string> Vtmp=Amap["[CommonAtom]"];
    string str=Vtmp[0];
    for(size_t jj=1;jj<Vtmp.size();jj++)
      str+=" "+Vtmp[jj];
    Matommap.emplace("CommonAtom",str);
    Amap.erase("[CommonAtom]");
  }
  for(auto & it : Amap) {
    string mname=it.first;
    vector<string> mlist=it.second;
    size_t found=mname.find("->");
    mname.replace(found,2," ");
    vector<string> stoken=split_string(mname,"[]");
    if(stoken.size()!=2) {
      string str="initialization.cpp: Wrong mutant name in mutation atommap --> "+mname;
      error->one(str);
    }
    string A2B=stoken[0]+"->"+stoken[1];
    string B2A=stoken[1]+"->"+stoken[0];
    for(auto line : mlist) {
      stoken=split_string(line,":");
      if(stoken.size()!=2) {
        cout<<"Wrong atom map --> "<<mname<<" : "<<line<<endl;
      }
      string atom1=A2B+"."+stoken[0];
      if(Matommap.count(atom1)) {
        string str="initialization.cpp: Atomname has been used in mutation atommap --> "+atom1;
        error->one(str);
      }
      Matommap.emplace(atom1,stoken[1]);
      atom1=B2A+"."+stoken[1];
      if(Matommap.count(atom1)) {
        string str="initialization.cpp: Atomname has been used in mutation atommap --> "+atom1;
        error->one(str);
      }
      Matommap.emplace(atom1,stoken[0]);
    }
  }

  return Matommap;
}

void Initialization::load_RNA_DM(int flag_use) {
  RNA_use_DM=flag_use;
  if(flag_use==0) return;

  //check RNA DM environment
  if(0) {
    string str="initialization.cpp: Cannot load RNA DM!";
    error->one(str);
  }

  //load RNA DM
}

void Initialization::set_residues_list() {
  NUCLEOTIDES_LIST.clear();
  string fstr=INPUT_DIR+"nucleotides_list.dat";
  vector<string> Vfile=read_file(fstr);
  for(size_t ii=0;ii<Vfile.size();ii++) {
    vector<string> stoken=split_string(Vfile[ii]);
    if(stoken.size()!=2) continue;
    string resname0=stoken[0];
    //for(size_t jj=0;jj<resname0.size();jj++)
      //resname0[jj]=toupper(resname0[jj]);
    string resname=stoken[1];
    for(size_t jj=0;jj<resname.size();jj++)
      resname[jj]=toupper(resname[jj]);
    if(resname=="T") resname="DT";
    NUCLEOTIDES_LIST.insert(pair<string,string>(resname0,resname));
  }

  RNA_STANDARD_NUCLEOTIDES.clear();
  string STD_NT[]={"A","ADE","RA","rA","RA3","RA5","C","CYT","RC","rC","RC3","RC5",
                   "G","GUA","RG","rG","RG3","RG5","U","URA","RU","rU","RU3","RU5",
                   "N"};
  for(auto nt : STD_NT)
    RNA_STANDARD_NUCLEOTIDES.push_back(nt);

  RESIDUES_LIST.clear();
  fstr=INPUT_DIR+"residues_list.dat";
  Vfile.clear();
  Vfile=read_file(fstr);
  for(size_t ii=0;ii<Vfile.size();ii++) {
    vector<string> stoken=split_string(Vfile[ii]);
    if(stoken.size()!=2) continue;
    string resname0=stoken[0];
    string resname=stoken[1];
    RESIDUES_LIST.emplace(resname0,resname);
  }

  PROTEIN_STANDARD_NUCLEOTIDES.clear();
  string STD_AA[]={"GLY","ALA","LEU","ILE","VAL","PRO","PHE","MET","TRP","SER","GLN","THR",
                   "CYS","ASN","TYR","ASP","GLU","LYS","ARG","HIS","HIP","HIE","HID"};
  for(auto aa : STD_AA)
    PROTEIN_STANDARD_NUCLEOTIDES.push_back(aa);
}

void Initialization::set_standard_residues() {
  RESIDUES_XYZ.clear();
  RESIDUES_ATOMS.clear();
  Atom atom;
  string fstr=INPUT_DIR+"RNA_nucleotide.pdb";
  vector<string> Vstr=read_file(fstr);
  vector<Atom> Vres;
  vector<string> Vname,Vrna,Vresname;
  string resname="";
  Vres.clear();
  Vname.clear();
  for(size_t ii=0;ii<Vstr.size();ii++) {
    string line=Vstr[ii];
    if(line.size()<=1) continue;
    if(line.find("#")!=string::npos) continue;
    if(line.find("[")!=string::npos && line.find("]")!=string::npos) {
      vector<string> stoken=split_string(line,"[]");
      if(Vres.size()!=0) {
        RESIDUES_XYZ.emplace(resname,Vres);
        RESIDUES_ATOMS.emplace(resname,Vname);
        Vresname.push_back(resname);
        Vres.clear();
        Vname.clear();
      }
      resname=stoken[0]; 
      continue;     
    }
    else if(line.find("ATOM")!=string::npos || line.find("HETATM")!=string::npos) {
      atom.set_atom(line);
      if(noH && atom.atomname[0]=='H') continue;
      atom.resname=resname;
      Vres.push_back(atom);
      Vname.push_back(atom.atomname);
      if(count(Vrna.begin(),Vrna.end(),atom.atomname)==0)
        Vrna.push_back(atom.atomname);  
      continue;
    }
  }
  if(Vres.size()>0) {
    RESIDUES_XYZ.emplace(resname,Vres);
    RESIDUES_ATOMS.emplace(resname,Vname);
  }
  RESIDUES_ATOMS.emplace("RNA",Vrna);
  Vresname.push_back(resname);
  Vres.clear();
  Vname.clear();

  if(RESIDUES_XYZ.count("rBB")==0) {
    string str="initialization.cpp: Error to set RNA nucleotide rBB!";
    error->one(str);
  }  
  if(RESIDUES_XYZ.count("A")==0) {
    string str="initialization.cpp: Error to set RNA nucleotide A!";
    error->one(str);
  }  
  if(RESIDUES_XYZ.count("C")==0) {
    string str="initialization.cpp: Error to set RNA nucleotide C!";
    error->one(str);
  } 
  if(RESIDUES_XYZ.count("G")==0) {
    string str="initialization.cpp: Error to set RNA nucleotide G!";
    error->one(str);
  } 
  if(RESIDUES_XYZ.count("U")==0) {
    string str="initialization.cpp: Error to set RNA nucleotide U!";
    error->one(str);
  } 
  if(RESIDUES_XYZ.count("BPi")==0) {
    string str="initialization.cpp: Error to set RNA nucleotide BPi!";
    error->one(str);
  } 
  if(RESIDUES_XYZ.count("BPj")==0) {
    string str="initialization.cpp: Error to set RNA nucleotide BPj!";
    error->one(str);
  }
  
  for(size_t ii=0;ii<Vresname.size();ii++) {
    string name=Vresname[ii];
    if(name=="rBB" || name=="BPi" || name=="BPj") continue;
    if(NUCLEOTIDES_LIST.count(name)==0)
      NUCLEOTIDES_LIST.emplace(name,name);
    if(count(RNA_STANDARD_NUCLEOTIDES.begin(),RNA_STANDARD_NUCLEOTIDES.end(),name)==0)
      RNA_STANDARD_NUCLEOTIDES.push_back(name);
  }

  /*
  for(auto& it : NUCLEOTIDES_LIST)
    cout<<it.first<<" : "<<it.second<<endl;

  for(auto nt : RNA_STANDARD_NUCLEOTIDES) 
    cout<<"Nt: "<<nt<<endl;
  
  for(auto aname : RESIDUES_ATOMS["RNA"]) 
    cout<<"Atom: "<<aname<<endl;
  */

  for(size_t ii=0;ii<RESIDUES_XYZ["A"].size();ii++) {
    atom=RESIDUES_XYZ["A"][ii];
    if(count(RESIDUES_ATOMS["rBB"].begin(),RESIDUES_ATOMS["rBB"].end(),atom.atomname)==0) {
      Vres.push_back(atom);
      Vname.push_back(atom.atomname);
    }
  }
  RESIDUES_XYZ.emplace("baseA",Vres);
  RESIDUES_ATOMS.emplace("baseA",Vname);
  Vres.clear();
  Vname.clear();

  for(size_t ii=0;ii<RESIDUES_XYZ["C"].size();ii++) {
    atom=RESIDUES_XYZ["C"][ii];
    if(count(RESIDUES_ATOMS["rBB"].begin(),RESIDUES_ATOMS["rBB"].end(),atom.atomname)==0) {
      Vres.push_back(atom);
      Vname.push_back(atom.atomname);
    }
  }
  RESIDUES_XYZ.emplace("baseC",Vres);
  RESIDUES_ATOMS.emplace("baseC",Vname);
  Vres.clear();
  Vname.clear();

  for(size_t ii=0;ii<RESIDUES_XYZ["G"].size();ii++) {
    atom=RESIDUES_XYZ["G"][ii];
    if(count(RESIDUES_ATOMS["rBB"].begin(),RESIDUES_ATOMS["rBB"].end(),atom.atomname)==0) {
      Vres.push_back(atom);
      Vname.push_back(atom.atomname);
    }
  }
  RESIDUES_XYZ.emplace("baseG",Vres);
  RESIDUES_ATOMS.emplace("baseG",Vname);
  Vres.clear();
  Vname.clear();

  for(size_t ii=0;ii<RESIDUES_XYZ["U"].size();ii++) {
    atom=RESIDUES_XYZ["U"][ii];
    if(count(RESIDUES_ATOMS["rBB"].begin(),RESIDUES_ATOMS["rBB"].end(),atom.atomname)==0) {
      Vres.push_back(atom);
      Vname.push_back(atom.atomname);
    }
  }
  RESIDUES_XYZ.emplace("baseU",Vres);
  RESIDUES_ATOMS.emplace("baseU",Vname);
  Vres.clear();
  Vname.clear();
}

void Initialization::set_residue_name_map() {
  SEQUENCE_MAP.clear();

  string fstr=INPUT_DIR+"sequence_map.dat";
  vector<string> Vlines=read_file(fstr);
  for(size_t ii=0;ii<Vlines.size();ii++) {
    string sline=Vlines[ii];
    if(sline.find("#")!=string::npos) {
      size_t pos=sline.find_first_of("#");
      sline=sline.substr(0,pos);
    }
    if(sline.size()<2) continue;
    vector<string> stoken=split_string(sline);
    if(stoken.size()!=3) {
      string fstr="initialization.cpp: Skip line in "+fstr+" --> "+sline;
      error->warning(fstr);
    }
    string type_flag=stoken[0];
    string full_resname=stoken[1];
    string one_letter=stoken[2];
    string sname=type_flag+"_"+full_resname;
    SEQUENCE_MAP.emplace(sname,one_letter);
    sname=type_flag+"_"+one_letter;
    SEQUENCE_MAP.emplace(sname,full_resname);
  }
}

void Initialization::set_RNA_basepair_notation() {
  RNA_BP_NOTATION.clear();
  char BPN[]={'.','(',')','[',']','{','}','<','>'};
  for(auto ch : BPN)
    RNA_BP_NOTATION.push_back(ch);
  char bpL0='A';
  char bpR0='a';
  for(int ii=0;ii<26;ii++) {
    char bpL=bpL0+ii;
    char bpR=bpR0+ii;
    RNA_BP_NOTATION.push_back(bpL);
    RNA_BP_NOTATION.push_back(bpR);
    //RNA_BP_NOTATION.push_back(bpR);
    //RNA_BP_NOTATION.push_back(bpL);
  }
}