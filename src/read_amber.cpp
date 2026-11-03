#include "basic.h"
#include "constant.h"
#include "readfile.h"
#include "initialization.h"
#include "error.h"
#include "read_amber.h"

using namespace CBS_NS;

const double LJ_CUT = 12.0;

ReadAmber::ReadAmber(Cobmms *cbs) : Pointers(cbs) {}

ReadAmber::ReadAmber(const ReadAmber &other) : Pointers(other.cbs) {
  Vatoms=other.Vatoms;
  Vbonds=other.Vbonds;
  Vangles=other.Vangles;
  Vdihedrals=other.Vdihedrals;
  Vimpropers=other.Vimpropers;
  Vmass=other.Vmass;
  Vbondcoeff=other.Vbondcoeff;
  Vanglecoeff=other.Vanglecoeff;
  Vdihedralcoeff=other.Vdihedralcoeff;
  Vimpropercoeff=other.Vimpropercoeff;
  Vpaircoeff=other.Vpaircoeff;
  Vboxsize=other.Vboxsize;
}

ReadAmber& ReadAmber::operator=(const ReadAmber &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    Vatoms=other.Vatoms;
    Vbonds=other.Vbonds;
    Vangles=other.Vangles;
    Vdihedrals=other.Vdihedrals;
    Vimpropers=other.Vimpropers;
    Vmass=other.Vmass;
    Vbondcoeff=other.Vbondcoeff;
    Vanglecoeff=other.Vanglecoeff;
    Vdihedralcoeff=other.Vdihedralcoeff;
    Vimpropercoeff=other.Vimpropercoeff;
    Vpaircoeff=other.Vpaircoeff;
    Vboxsize=other.Vboxsize;
  }

  return *this;
}

ReadAmber::~ReadAmber() {
  Vatoms.clear();
  Vatoms.shrink_to_fit();
  Vbonds.clear();
  Vbonds.shrink_to_fit();
  Vangles.clear();
  Vangles.shrink_to_fit();
  Vdihedrals.clear();
  Vdihedrals.shrink_to_fit();
  Vimpropers.clear();
  Vimpropers.shrink_to_fit();
  Vmass.clear();
  Vmass.shrink_to_fit();
  Vbondcoeff.clear();
  Vbondcoeff.shrink_to_fit();
  Vanglecoeff.clear();
  Vanglecoeff.shrink_to_fit();
  Vdihedralcoeff.clear();
  Vdihedralcoeff.shrink_to_fit();
  Vimpropercoeff.clear();
  Vimpropercoeff.shrink_to_fit();
  Vpaircoeff.clear();
  Vpaircoeff.shrink_to_fit();
  Vboxsize.clear();
  Vboxsize.shrink_to_fit();
}

bool ReadAmber::check_tleap() {
  string cmd="which tleap";
  string fout="./out.dat";
  cmd+=" > "+fout;

  int sys=system(cmd.c_str());
  if(sys!=0) {
    string str="Fail to execute "+cmd+" --> sys = "+to_string(sys);
    cout<<"Error: "<<str<<endl;
    return false;
  }

  vector<string> Vlines=read_file(fout);

  cmd="rm "+fout;
  sys=system(cmd.c_str());
  if(sys!=0) {
    string str="Fail to execute "+cmd;
    cout<<"Error: "<<str<<endl;
  }

  if(Vlines.size()==0) return false;
  else return true;
}

void ReadAmber::prepare_amber_files(const Molecule &mol_one, string fname, int fpause) {
  string fstr=fname+"_in.pdb";
  ofstream outF (fstr.c_str());
  init->output_format="amb";
  outF<<mol_one;
  outF.close(); outF.clear();
  init->output_format="std";
  
  ReadAmber::prepare_amber_files(fstr,fname,fpause);
}

void ReadAmber::prepare_amber_files(string infile, string fname, int fpause) {
  if(!ReadAmber::check_tleap()) {
    string str="read_amber.cpp: Cannot find tleap!";
    error->one(str);
  }

  string cmd="pdb4amber -i "+infile+" -o "+fname+"_amber.pdb -y --reduce --add-missing-atoms";
  cbs->execute_system_command(cmd);

  if(fpause) pause();
  
  string fleap=init->WORK_DIR+"leap.in";
  ofstream outF (fleap.c_str());
  outF<<"source leaprc.protein.ff14SB"<<endl;
  outF<<"source leaprc.DNA.OL15"<<endl;
  outF<<"source leaprc.RNA.OL3"<<endl;
  outF<<"source leaprc.water.tip3p"<<endl;
  outF<<"loadamberparams frcmod.ionsjc_tip3p"<<endl;
  outF<<endl;
  outF<<"mol = loadpdb "<<fname<<"_amber.pdb"<<endl;
  outF<<"saveamberparm mol "<<fname<<".prmtop "<<fname<<".inpcrd"<<endl;
  outF<<"quit"<<endl;
  outF.close(); outF.clear();

  cmd="tleap -s -f "+fleap+" > "+init->WORK_DIR+"leap.out";
  cbs->execute_system_command(cmd);
  cmd="rm reduce_info.log leap.log "+fleap+" "+init->WORK_DIR+"leap.out";
  cbs->execute_system_command(cmd);
}

void ReadAmber::prepare_gromacs_files(const Molecule &mol_one, string fname, const GmxBoxSet &boxset_one, int fpause) {
  string fstr=fname+"_in.pdb";
  ofstream outF (fstr.c_str());
  init->output_format="amb";
  outF<<mol_one;
  outF.close(); outF.clear();
  init->output_format="std";
  
  ReadAmber::prepare_gromacs_files(fstr,fname,boxset_one,fpause);
}

void ReadAmber::prepare_gromacs_files(string infile, string fname, const GmxBoxSet &boxset_one, int fpause) {
  if(!ReadAmber::check_tleap()) {
    string str="read_amber.cpp: Cannot find tleap!";
    error->one(str);
  }

  vector<string> Vlines=read_file(infile);
  string fstr=fname+"_in.pdb";
  ofstream outF (fstr.c_str());
  for(size_t jj=0;jj<Vlines.size();jj++) {
    string sline=Vlines[jj];
    string shead;
    if(sline.find("TER")!=string::npos) shead="TER";
    else {
      shead=sline.substr(0,6);
    }
    if(shead.find("ATOM")==string::npos && shead.find("HETATM")==string::npos && 
       shead.find("TER")==string::npos) continue;
    outF<<sline<<endl;
  }
  outF.close(); outF.clear();

  string cmd="pdb4amber -i "+fstr+" -o "+fname+"_amber.pdb -y --reduce --add-missing-atoms";
  cbs->execute_system_command(cmd);

  //specify for YTHDF1-m6A simulation
  string famb=fname+"_amber.pdb";
  Vlines.clear();
  Vlines=read_file(famb);
  outF.open(famb.c_str());
  for(size_t ll=0;ll<Vlines.size();ll++) {
    if(Vlines[ll].find("END")!=string::npos) continue;
    string sline=Vlines[ll];
    if(sline.find("TER")!=string::npos && (sline.find(" G ")!=string::npos || sline.find("M6A")!=string::npos)) 
      continue;
    outF<<Vlines[ll]<<endl;
  }
  outF.close(); outF.clear();


  if(fpause) pause();

  string fleap_in=init->WORK_DIR+"leap.in";
  string fleap_out=init->WORK_DIR+"leap.out";
  string ftail_out=init->WORK_DIR+"tail.dat";
  int number_ions=0;
  if(boxset_one.concentration>0.0) {
    outF.open(fleap_in.c_str());
    outF<<"source leaprc.protein.ff14SB"<<endl;
    outF<<"source leaprc.DNA.OL15"<<endl;
    outF<<"source leaprc.RNA.OL3"<<endl;
    outF<<"source leaprc.user"<<endl;
    outF<<"source leaprc.water.tip3p"<<endl;
    outF<<"source leaprc.gaff2"<<endl;
    outF<<"loadamberparams frcmod.ionsjc_tip3p"<<endl;
    outF<<endl;
    outF<<"mol = loadpdb "<<fname<<"_amber.pdb"<<endl;

    if(boxset_one.boxshape=="OCT") {
      outF<<"solvateOct mol TIP3PBOX "<<boxset_one.distance<<" iso"<<endl;
    }
    else if(boxset_one.boxshape=="CUB") {
      outF<<"solvateBox mol TIP3PBOX "<<boxset_one.distance<<" iso"<<endl;
    }

    outF<<"quit"<<endl;
    outF.close(); outF.clear();
    
    cmd="tleap -s -f "+fleap_in+" > "+fleap_out;
    cbs->execute_system_command(cmd);
    
    cmd="tail "+fleap_out+" -n 4 > "+ftail_out;
    cbs->execute_system_command(cmd);

    Vlines.clear();
    Vlines=read_file(ftail_out);
    vector<string> stoken=split_string(Vlines[0]);
    int nwater=stoi(stoken[1]);
    number_ions=round(nwater*0.0187*boxset_one.concentration);
  }

  outF.open(fleap_in.c_str());
  outF<<"source leaprc.protein.ff14SB"<<endl;
  outF<<"source leaprc.DNA.OL15"<<endl;
  outF<<"source leaprc.RNA.OL3"<<endl;
  outF<<"source leaprc.user"<<endl;
  outF<<"source leaprc.water.tip3p"<<endl;
  outF<<"source leaprc.gaff2"<<endl;
  outF<<"loadamberparams frcmod.ionsjc_tip3p"<<endl;
  outF<<endl;
  outF<<"mol = loadpdb "<<fname<<"_amber.pdb"<<endl;

  if(boxset_one.boxshape=="OCT") {
    outF<<"solvateOct mol TIP3PBOX "<<boxset_one.distance<<" iso"<<endl;
  }
  else if(boxset_one.boxshape=="CUB") {
    outF<<"solvateBox mol TIP3PBOX "<<boxset_one.distance<<" iso"<<endl;
  }

  if(boxset_one.neutralization==1) {
    outF<<"addions mol "<<boxset_one.positivename<<" 0"<<endl;
    outF<<"addions mol "<<boxset_one.negativename<<" 0"<<endl;
  }

  if(boxset_one.concentration>0.0) {
    outF<<"addions mol "<<boxset_one.positivename<<" "<<number_ions<<endl;
    outF<<"addions mol "<<boxset_one.negativename<<" 0"<<endl;
  }
  
  outF<<"saveamberparm mol "<<fname<<".prmtop "<<fname<<".inpcrd"<<endl;
  outF<<"quit"<<endl;
  outF.close(); outF.clear();

  cmd="tleap -s -f "+fleap_in+" > "+fleap_out;
  cbs->execute_system_command(cmd);

  cmd="amb2gro_top_gro.py -p "+fname+".prmtop -c "+fname+".inpcrd -t "+fname;
  cmd+=".top -g "+fname+".gro -b "+fname+"_water_ions.pdb";
  cbs->execute_system_command(cmd);

  cmd="rm reduce_info.log leap.log "+fleap_in+" "+fleap_out+" "+ftail_out;
  cbs->execute_system_command(cmd);
}

void ReadAmber::read_prmtop(string fname) {
  if(fname.find(".prmtop")==string::npos)
    fname+=".prmtop";

  if(!is_file_exist(fname)) {
    string str="read_amber.cpp: Cannot find file "+fname;
    error->one(str);
  }

  vector<string> Vlines=read_file(fname);
  map<string,vector<string> > Mblocks;
  size_t nline=0;
  string block_name="";
  vector<string> Vcontent;
  //int repeat_count;
  int field_width=4;
  while(nline<Vlines.size()) {
    string sline=Vlines[nline];
    nline++;
    if(sline[0]=='%' && sline.find("FLAG")==string::npos && sline.find("FORMAT")==string::npos) {
      continue;
    }
    else if(sline[0]=='%' && sline.find("FLAG")!=string::npos) {
      if(Vcontent.size()>0) {
        Mblocks.emplace(block_name,Vcontent);
        Vcontent.clear();
      }
      vector<string> stoken=split_string(sline.substr(6));
      block_name=stoken[0];
      continue;
    }
    else if(sline[0]=='%' && sline.find("FORMAT")!=string::npos) {
      pair<int,int> p_one=ReadAmber::parse_format(sline);
      //repeat_count=p_one.first;
      field_width=p_one.second;
      continue;
    }
    else {
      vector<string> Vtmp=ReadAmber::parse_line(sline,field_width);
      Vcontent.insert(Vcontent.end(),Vtmp.begin(),Vtmp.end());
      continue;
    }
  }
  if(Vcontent.size()>0) {
    Mblocks.emplace(block_name,Vcontent);
    Vcontent.clear();
  }
  
  if(Mblocks.count("POINTERS")==0) {
    string str="read_amber.cpp: Cannot find block POINTERS in "+fname;
    error->one(str);
  }
  Vcontent.clear();
  Vcontent=Mblocks["POINTERS"];
  if(Vcontent.size()<31) {
    string str="read_amber.cpp: block POINTERS in "+fname+" is too short!";
    error->one(str);
  }
  int Natom=stoi(Vcontent[0]);
  int Ntypes=stoi(Vcontent[1]);
  int Nbonh=stoi(Vcontent[2]);
  //int Mbona=stoi(Vcontent[3]);
  int Ntheth=stoi(Vcontent[4]);
  //int Mtheta=stoi(Vcontent[5]);
  int Nphih=stoi(Vcontent[6]);
  //int Mphia=stoi(Vcontent[7]);
  //int Nhparm=stoi(Vcontent[8]);
  //int Nparm=stoi(Vcontent[9]);
  //int Nnb=stoi(Vcontent[10]);
  int Nres=stoi(Vcontent[11]);
  int Nbona=stoi(Vcontent[12]);
  int Ntheta=stoi(Vcontent[13]);
  int Nphia=stoi(Vcontent[14]);
  int Numbnd=stoi(Vcontent[15]);
  int Numang=stoi(Vcontent[16]);
  int Nptra=stoi(Vcontent[17]);
  //int Natyp=stoi(Vcontent[18]);
  //int Nphb=stoi(Vcontent[19]);
  //int Ifpert=stoi(Vcontent[20]);
  //int Nbper=stoi(Vcontent[21]);
  //int Ngper=stoi(Vcontent[22]);
  //int Ndper=stoi(Vcontent[23]);
  //int Mbper=stoi(Vcontent[24]);
  //int Mgper=stoi(Vcontent[25]);
  //int Mdper=stoi(Vcontent[26]);
  int Ifbox=stoi(Vcontent[27]);
  //int Nmxrs=stoi(Vcontent[28]);
  //int Ifcap=stoi(Vcontent[29]);
  //int Numextra=stoi(Vcontent[30]);

  if(Mblocks.count("ATOM_NAME")==0) {
    string str="read_amber.cpp: Cannot find block ATOM_NAME in "+fname;
    error->one(str);
  }
  Vcontent.clear();
  Vcontent=Mblocks["ATOM_NAME"];
  if((int)Vcontent.size()!=Natom) {
    string str="read_amber.cpp: block ATOM_NAME in "+fname+" is wrong  -->  ";
    str+=to_string(Vcontent.size())+" / "+to_string(Natom);
    error->one(str);
  }
  vector<string> Vatom_name;
  for(size_t jj=0;jj<Vcontent.size();jj++) {
    vector<string> stoken=split_string(Vcontent[jj]);
    Vatom_name.push_back(stoken[0]);
  }

  if(Mblocks.count("CHARGE")==0) {
    string str="read_amber.cpp: Cannot find block CHARGE in "+fname;
    error->one(str);
  }
  Vcontent.clear();
  Vcontent=Mblocks["CHARGE"];
  if((int)Vcontent.size()!=Natom) {
    string str="read_amber.cpp: block CHARGE in "+fname+" is wrong  -->  ";
    str+=to_string(Vcontent.size())+" / "+to_string(Natom);
    error->one(str);
  }
  vector<double> Vcharge;
  for(size_t jj=0;jj<Vcontent.size();jj++) {
    double q_i=stod(Vcontent[jj]); ///18.2223;
    Vcharge.push_back(q_i);
  }

  if(Mblocks.count("MASS")==0) {
    string str="read_amber.cpp: Cannot find block MASS in "+fname;
    error->one(str);
  }
  Vcontent.clear();
  Vcontent=Mblocks["MASS"];
  if((int)Vcontent.size()!=Natom) {
    string str="read_amber.cpp: block MASS in "+fname+" is wrong  -->  ";
    str+=to_string(Vcontent.size())+" / "+to_string(Natom);
    error->one(str);
  }
  //vector<string> Vmass_all=Vcontent;
  vector<double> Vmass_all;
  for(size_t jj=0;jj<Vcontent.size();jj++) {
    double mass_i=stod(Vcontent[jj]);
    Vmass_all.push_back(mass_i);
  }

  if(Mblocks.count("ATOM_TYPE_INDEX")==0) {
    string str="read_amber.cpp: Cannot find block ATOM_TYPE_INDEX in "+fname;
    error->one(str);
  }
  Vcontent.clear();
  Vcontent=Mblocks["ATOM_TYPE_INDEX"];
  if((int)Vcontent.size()!=Natom) {
    string str="read_amber.cpp: block ATOM_TYPE_INDEX in "+fname+" is wrong  -->  ";
    str+=to_string(Vcontent.size())+" / "+to_string(Natom);
    error->one(str);
  }
  vector<int> Vatom_type_index;
  for(size_t jj=0;jj<Vcontent.size();jj++) {
    Vatom_type_index.push_back(stoi(Vcontent[jj]));
  }

  
  if(Mblocks.count("NONBONDED_PARM_INDEX")==0) {
    string str="read_amber.cpp: Cannot find block NONBONDED_PARM_INDEX in "+fname;
    error->one(str);
  }
  Vcontent.clear();
  Vcontent=Mblocks["NONBONDED_PARM_INDEX"];
  if((int)Vcontent.size()!=Ntypes*Ntypes) {
    string str="read_amber.cpp: block NONBONDED_PARM_INDEX in "+fname+" is wrong  -->  ";
    str+=to_string(Vcontent.size())+" / "+to_string(Ntypes*Ntypes);
    error->one(str);
  }
  vector<int> Vnonbonded_parm_index;
  for(size_t jj=0;jj<Vcontent.size();jj++) {
    Vnonbonded_parm_index.push_back(stoi(Vcontent[jj]));
  }

  vector<string> Vresidue_label;
  vector<int> Vresidue_pointer;
  if(Nres>0) {
    if(Mblocks.count("RESIDUE_LABEL")==0) {
      string str="read_amber.cpp: Cannot find block RESIDUE_LABEL in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["RESIDUE_LABEL"];
    if((int)Vcontent.size()!=Nres) {
      string str="read_amber.cpp: block RESIDUE_LABEL in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(Nres);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      vector<string> stoken=split_string(Vcontent[jj]);
      Vresidue_label.push_back(stoken[0]);
    }
  
    if(Mblocks.count("RESIDUE_POINTER")==0) {
      string str="read_amber.cpp: Cannot find block RESIDUE_POINTER in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["RESIDUE_POINTER"];
    if((int)Vcontent.size()!=Nres) {
      string str="read_amber.cpp: block RESIDUE_POINTER in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(Nres);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      Vresidue_pointer.push_back(stoi(Vcontent[jj]));
    }
  }

  vector<double> Vbond_force_constant, Vbond_equil_value;
  if(Numbnd>0) {
    if(Mblocks.count("BOND_FORCE_CONSTANT")==0) {
      string str="read_amber.cpp: Cannot find block BOND_FORCE_CONSTANT in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["BOND_FORCE_CONSTANT"];
    if((int)Vcontent.size()!=Numbnd) {
      string str="read_amber.cpp: block BOND_FORCE_CONSTANT in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(Numbnd);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      double k_i=stod(Vcontent[jj]); ///2.0;
      Vbond_force_constant.push_back(k_i);
    }

    if(Mblocks.count("BOND_EQUIL_VALUE")==0) {
      string str="read_amber.cpp: Cannot find block BOND_EQUIL_VALUE in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["BOND_EQUIL_VALUE"];
    if((int)Vcontent.size()!=Numbnd) {
      string str="read_amber.cpp: block BOND_EQUIL_VALUE in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(Numbnd);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      double b0_i=stod(Vcontent[jj]);
      Vbond_equil_value.push_back(b0_i);
    }
  }

  vector<double> Vangle_force_constant, Vangle_equil_value;
  if(Numang>0) {
    if(Mblocks.count("ANGLE_FORCE_CONSTANT")==0) {
      string str="read_amber.cpp: Cannot find block ANGLE_FORCE_CONSTANT in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["ANGLE_FORCE_CONSTANT"];
    if((int)Vcontent.size()!=Numang) {
      string str="read_amber.cpp: block ANGLE_FORCE_CONSTANT in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(Numang);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      double k_i=stod(Vcontent[jj]); ///2.0;
      Vangle_force_constant.push_back(k_i);
    }

    if(Mblocks.count("ANGLE_EQUIL_VALUE")==0) {
      string str="read_amber.cpp: Cannot find block ANGLE_EQUIL_VALUE in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["ANGLE_EQUIL_VALUE"];
    if((int)Vcontent.size()!=Numang) {
      string str="read_amber.cpp: block ANGLE_EQUIL_VALUE in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(Numang);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      double theta0_i=stod(Vcontent[jj]); //*180/PI;
      Vangle_equil_value.push_back(theta0_i);
    }
  }

  vector<double> Vdihedral_force_constant, Vdihedral_phase;
  vector<int> Vdihedral_periodicity;
  if(Nptra>0) {
    if(Mblocks.count("DIHEDRAL_FORCE_CONSTANT")==0) {
      string str="read_amber.cpp: Cannot find block DIHEDRAL_FORCE_CONSTANT in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["DIHEDRAL_FORCE_CONSTANT"];
    if((int)Vcontent.size()!=Nptra) {
      string str="read_amber.cpp: block DIHEDRAL_FORCE_CONSTANT in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(Nptra);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      double k_i=stod(Vcontent[jj]);
      Vdihedral_force_constant.push_back(k_i);
    }

    if(Mblocks.count("DIHEDRAL_PERIODICITY")==0) {
      string str="read_amber.cpp: Cannot find block DIHEDRAL_PERIODICITY in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["DIHEDRAL_PERIODICITY"];
    if((int)Vcontent.size()!=Nptra) {
      string str="read_amber.cpp: block DIHEDRAL_PERIODICITY in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(Nptra);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      double n_i=stod(Vcontent[jj]);
      Vdihedral_periodicity.push_back(round(n_i));
    }

    if(Mblocks.count("DIHEDRAL_PHASE")==0) {
      string str="read_amber.cpp: Cannot find block DIHEDRAL_PHASE in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["DIHEDRAL_PHASE"];
    if((int)Vcontent.size()!=Nptra) {
      string str="read_amber.cpp: block DIHEDRAL_PHASE in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(Nptra);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      double phi_i=stod(Vcontent[jj]);
      Vdihedral_phase.push_back(phi_i);
    }
  }

  if(Mblocks.count("LENNARD_JONES_ACOEF")==0) {
    string str="read_amber.cpp: Cannot find block LENNARD_JONES_ACOEF in "+fname;
    error->one(str);
  }
  Vcontent.clear();
  Vcontent=Mblocks["LENNARD_JONES_ACOEF"];
  if((int)Vcontent.size()!=(Ntypes*(Ntypes+1)/2)) {
    string str="read_amber.cpp: block LENNARD_JONES_ACOEF in "+fname+" is wrong  -->  ";
    str+=to_string(Vcontent.size())+" / "+to_string(Ntypes*(Ntypes+1)/2);
    error->one(str);
  }
  vector<double> Vlj_acoef;
  for(size_t jj=0;jj<Vcontent.size();jj++) {
    double a_i=stod(Vcontent[jj]);
    Vlj_acoef.push_back(a_i);
  }
 
  if(Mblocks.count("LENNARD_JONES_BCOEF")==0) {
    string str="read_amber.cpp: Cannot find block LENNARD_JONES_BCOEF in "+fname;
    error->one(str);
  }
  Vcontent.clear();
  Vcontent=Mblocks["LENNARD_JONES_BCOEF"];
  if((int)Vcontent.size()!=(Ntypes*(Ntypes+1)/2)) {
    string str="read_amber.cpp: block LENNARD_JONES_BCOEF in "+fname+" is wrong  -->  ";
    str+=to_string(Vcontent.size())+" / "+to_string(Ntypes*(Ntypes+1)/2);
    error->one(str);
  }
  vector<double> Vlj_bcoef;
  for(size_t jj=0;jj<Vcontent.size();jj++) {
    double b_i=stod(Vcontent[jj]);
    Vlj_bcoef.push_back(b_i);
  }
  
  vector<int> Vbonds_inc_hydrogen;
  if(Nbonh>0) {
    if(Mblocks.count("BONDS_INC_HYDROGEN")==0) {
      string str="read_amber.cpp: Cannot find block BONDS_INC_HYDROGEN in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["BONDS_INC_HYDROGEN"];
    if((int)Vcontent.size()!=3*Nbonh) {
      string str="read_amber.cpp: block BONDS_INC_HYDROGEN in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(3*Nbonh);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      Vbonds_inc_hydrogen.push_back(stoi(Vcontent[jj]));
    }
  }

  vector<int> Vbonds_without_hydrogen;
  if(Nbona>0) {
    if(Mblocks.count("BONDS_WITHOUT_HYDROGEN")==0) {
      string str="read_amber.cpp: Cannot find block BONDS_WITHOUT_HYDROGEN in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["BONDS_WITHOUT_HYDROGEN"];
    if((int)Vcontent.size()!=3*Nbona) {
      string str="read_amber.cpp: block BONDS_WITHOUT_HYDROGEN in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(3*Nbona);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      Vbonds_without_hydrogen.push_back(stoi(Vcontent[jj]));
    }
  }

  vector<int> Vangles_inc_hydrogen;
  if(Ntheth>0) {
    if(Mblocks.count("ANGLES_INC_HYDROGEN")==0) {
      string str="read_amber.cpp: Cannot find block ANGLES_INC_HYDROGEN in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["ANGLES_INC_HYDROGEN"];
    if((int)Vcontent.size()!=4*Ntheth) {
      string str="read_amber.cpp: block ANGLES_INC_HYDROGEN in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(4*Ntheth);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      Vangles_inc_hydrogen.push_back(stoi(Vcontent[jj]));
    }
  }

  vector<int> Vangles_without_hydrogen;
  if(Ntheta>0) {
    if(Mblocks.count("ANGLES_WITHOUT_HYDROGEN")==0) {
      string str="read_amber.cpp: Cannot find block ANGLES_WITHOUT_HYDROGEN in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["ANGLES_WITHOUT_HYDROGEN"];
    if((int)Vcontent.size()!=4*Ntheta) {
      string str="read_amber.cpp: block ANGLES_WITHOUT_HYDROGEN in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(4*Ntheta);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      Vangles_without_hydrogen.push_back(stoi(Vcontent[jj]));
    }
  }

  vector<int> Vdihedrals_inc_hydrogen;
  if(Nphih>0) {
    if(Mblocks.count("DIHEDRALS_INC_HYDROGEN")==0) {
      string str="read_amber.cpp: Cannot find block DIHEDRALS_INC_HYDROGEN in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["DIHEDRALS_INC_HYDROGEN"];
    if((int)Vcontent.size()!=5*Nphih) {
      string str="read_amber.cpp: block DIHEDRALS_INC_HYDROGEN in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(5*Nphih);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      Vdihedrals_inc_hydrogen.push_back(stoi(Vcontent[jj]));
    }
  }

  vector<int> Vdihedrals_without_hydrogen;
  if(Nphia>0) {
    if(Mblocks.count("DIHEDRALS_WITHOUT_HYDROGEN")==0) {
      string str="read_amber.cpp: Cannot find block DIHEDRALS_WITHOUT_HYDROGEN in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["DIHEDRALS_WITHOUT_HYDROGEN"];
    if((int)Vcontent.size()!=5*Nphia) {
      string str="read_amber.cpp: block DIHEDRALS_WITHOUT_HYDROGEN in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / "+to_string(5*Nphia);
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      Vdihedrals_without_hydrogen.push_back(stoi(Vcontent[jj]));
    }
  }

  if(Mblocks.count("AMBER_ATOM_TYPE")==0) {
    string str="read_amber.cpp: Cannot find block AMBER_ATOM_TYPE in "+fname;
    error->one(str);
  }
  Vcontent.clear();
  Vcontent=Mblocks["AMBER_ATOM_TYPE"];
  if((int)Vcontent.size()!=Natom) {
    string str="read_amber.cpp: block AMBER_ATOM_TYPE in "+fname+" is wrong  -->  ";
    str+=to_string(Vcontent.size())+" / "+to_string(Natom);
    error->one(str);
  }
  vector<string> Vamber_atom_type;
  for(size_t jj=0;jj<Vcontent.size();jj++) {
    vector<string> stoken=split_string(Vcontent[jj]);
    Vamber_atom_type.push_back(stoken[0]);
  }

  //int Ipters=0, Nspsol=0;
  int Nspm=0;
  vector<int> Vatoms_per_molecule;
  vector<double> Vbox_dimensions;
  if(Ifbox>0) {
    if(Mblocks.count("SOLVENT_POINTERS")==0) {
      string str="read_amber.cpp: Cannot find block SOLVENT_POINTERS in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["SOLVENT_POINTERS"];
    if((int)Vcontent.size()!=3) {
      string str="read_amber.cpp: block SOLVENT_POINTERS in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / 3";
      error->one(str);
    }
    //Ipters=stoi(Vcontent[0]);
    Nspm=stoi(Vcontent[1]);
    //Nspsol=stoi(Vcontent[2]);

    if(Nspm>0) {
      if(Mblocks.count("ATOMS_PER_MOLECULE")==0) {
        string str="read_amber.cpp: Cannot find block ATOMS_PER_MOLECULE in "+fname;
        error->one(str);
      }
      Vcontent.clear();
      Vcontent=Mblocks["ATOMS_PER_MOLECULE"];
      if((int)Vcontent.size()!=Nspm) {
        string str="read_amber.cpp: block ATOMS_PER_MOLECULE in "+fname+" is wrong  -->  ";
        str+=to_string(Vcontent.size())+" / "+to_string(Nspm);
        error->one(str);
      }
      for(size_t jj=0;jj<Vcontent.size();jj++) {
        Vatoms_per_molecule.push_back(stoi(Vcontent[jj]));
      }
    }

    if(Mblocks.count("BOX_DIMENSIONS")==0) {
      string str="read_amber.cpp: Cannot find block BOX_DIMENSIONS in "+fname;
      error->one(str);
    }
    Vcontent.clear();
    Vcontent=Mblocks["BOX_DIMENSIONS"];
    if((int)Vcontent.size()!=4) {
      string str="read_amber.cpp: block BOX_DIMENSIONS in "+fname+" is wrong  -->  ";
      str+=to_string(Vcontent.size())+" / 4";
      error->one(str);
    }
    for(size_t jj=0;jj<Vcontent.size();jj++) {
      double box_i=stod(Vcontent[jj]);
      Vbox_dimensions.push_back(box_i);
    }
  }

  Vatoms.clear();
  for(int jj=0;jj<Natom;jj++) {
    Atom atom_one;
    atom_one.id=jj+1;
    atom_one.resid=1;
    atom_one.chainid="0";
    atom_one.resname=" ";
    atom_one.atomname="t"+to_string(Vatom_type_index[jj]);
    atom_one.q=Vcharge[jj]/18.2223;
    Vatoms.push_back(atom_one);
  }
  for(int jj=0;jj<Nres;jj++) {
    int begin_atom_index=Vresidue_pointer[jj]-1;
    int end_atom_index;
    if(jj<Nres-1) {
      end_atom_index=Vresidue_pointer[jj+1]-1;
    }
    else 
      end_atom_index=Natom-1;
    for(int kk=begin_atom_index;kk<end_atom_index;kk++) {
      Vatoms[kk].resid=jj+1;
      Vatoms[kk].resname=Vresidue_label[jj];
    }
  }
  int natom_sum=0;
  for(size_t jj=0;jj<Vatoms_per_molecule.size();jj++) {
    for(int kk=0;kk<Vatoms_per_molecule[jj];kk++) {
      Vatoms[kk+natom_sum].chainid=to_string(jj);
    }
    natom_sum+=Vatoms_per_molecule[jj];
  }
  //cout<<"Vatoms : "<<Vatoms.size()<<endl;

  Vmass.clear();
  Vmass.resize(Ntypes);
  for(int jj=0;jj<Natom;jj++) {
    string smass="t"+to_string(Vatom_type_index[jj])+" "+to_string(Vmass_all[jj]);
    Vmass[Vatom_type_index[jj]-1]=smass;
  }
  //cout<<"Vmass : "<<Vmass.size()<<endl;

  Vbonds.clear();
  for(int jj=0;jj<Nbonh;jj++) {
    Bond bond_one;
    bond_one.atom_i=abs(Vbonds_inc_hydrogen[3*jj])/3+1;
    bond_one.atom_j=abs(Vbonds_inc_hydrogen[3*jj+1])/3+1;
    bond_one.type="b"+to_string(Vbonds_inc_hydrogen[3*jj+2]);
    Vbonds.push_back(bond_one);
  } 
  for(int jj=0;jj<Nbona;jj++) {
    Bond bond_one;
    bond_one.atom_i=abs(Vbonds_without_hydrogen[3*jj])/3+1;
    bond_one.atom_j=abs(Vbonds_without_hydrogen[3*jj+1])/3+1;
    bond_one.type="b"+to_string(Vbonds_without_hydrogen[3*jj+2]);
    Vbonds.push_back(bond_one);
  }
  //cout<<"Vbonds : "<<Vbonds.size()<<endl;

  Vangles.clear();
  for(int jj=0;jj<Ntheth;jj++) {
    Angle angle_one;
    angle_one.atom_i=abs(Vangles_inc_hydrogen[4*jj])/3+1;
    angle_one.atom_j=abs(Vangles_inc_hydrogen[4*jj+1])/3+1;
    angle_one.atom_k=abs(Vangles_inc_hydrogen[4*jj+2])/3+1;
    angle_one.type="a"+to_string(Vangles_inc_hydrogen[4*jj+3]);
    Vangles.push_back(angle_one);
  }
  for(int jj=0;jj<Ntheta;jj++) {
    Angle angle_one;
    angle_one.atom_i=abs(Vangles_without_hydrogen[4*jj])/3+1;
    angle_one.atom_j=abs(Vangles_without_hydrogen[4*jj+1])/3+1;
    angle_one.atom_k=abs(Vangles_without_hydrogen[4*jj+2])/3+1;
    angle_one.type="a"+to_string(Vangles_without_hydrogen[4*jj+3]);
    Vangles.push_back(angle_one);
  }
  //cout<<"Vangles : "<<Vangles.size()<<endl;

  Vdihedrals.clear();
  for(int jj=0;jj<Nphih;jj++) {
    Dihedral dihedral_one;
    dihedral_one.atom_i=abs(Vdihedrals_inc_hydrogen[5*jj])/3+1;
    dihedral_one.atom_j=abs(Vdihedrals_inc_hydrogen[5*jj+1])/3+1;
    dihedral_one.atom_k=abs(Vdihedrals_inc_hydrogen[5*jj+2])/3+1;
    dihedral_one.atom_l=abs(Vdihedrals_inc_hydrogen[5*jj+3])/3+1;
    dihedral_one.type="d"+to_string(Vdihedrals_inc_hydrogen[5*jj+4]);
    Vdihedrals.push_back(dihedral_one);
  }
  for(int jj=0;jj<Nphia;jj++) {
    Dihedral dihedral_one;
    dihedral_one.atom_i=abs(Vdihedrals_without_hydrogen[5*jj])/3+1;
    dihedral_one.atom_j=abs(Vdihedrals_without_hydrogen[5*jj+1])/3+1;
    dihedral_one.atom_k=abs(Vdihedrals_without_hydrogen[5*jj+2])/3+1;
    dihedral_one.atom_l=abs(Vdihedrals_without_hydrogen[5*jj+3])/3+1;
    dihedral_one.type="d"+to_string(Vdihedrals_without_hydrogen[5*jj+4]);
    Vdihedrals.push_back(dihedral_one);
  }
  //cout<<"Vimpropers : "<<Vimpropers.size()<<endl;

  Vimpropers.clear();

  Vbondcoeff.clear();
  for(int jj=0;jj<Numbnd;jj++) {
    string sbcoef="b"+to_string(jj+1)+" harmonic "+to_string(Vbond_force_constant[jj]/2.0);
    sbcoef+=" "+to_string(Vbond_equil_value[jj]);
    Vbondcoeff.push_back(sbcoef);
  }
  //cout<<"Vbondcoeff : "<<Vbondcoeff.size()<<endl;

  Vanglecoeff.clear();
  for(int jj=0;jj<Numang;jj++) {
    string sacoef="a"+to_string(jj+1)+" harmonic "+to_string(Vangle_force_constant[jj]/2.0);
    sacoef+=" "+to_string(Vangle_equil_value[jj]*180.0/MathConst::PI);
    Vanglecoeff.push_back(sacoef);
  }
  //cout<<"Vanglecoeff : "<<Vanglecoeff.size()<<endl;

  Vdihedralcoeff.clear();
  for(int jj=0;jj<Nptra;jj++) {
    string sdcoef="d"+to_string(jj+1)+" fourier 1 "+to_string(Vdihedral_force_constant[jj]);
    sdcoef+=" "+to_string(Vdihedral_periodicity[jj])+" "+to_string(-Vdihedral_phase[jj]*180/MathConst::PI);
    Vdihedralcoeff.push_back(sdcoef);
  }
  //cout<<"Vdihedralcoeff : "<<Vdihedralcoeff.size()<<endl;

  Vimpropercoeff.clear();

  Vpaircoeff.clear();
  for(int ii=0;ii<Ntypes;ii++) {
    for(int jj=ii;jj<Ntypes;jj++) {
      int index_ij=abs(Vnonbonded_parm_index[Ntypes*ii+jj])-1;
      double lja_ij=Vlj_acoef[index_ij];
      double ljb_ij=Vlj_bcoef[index_ij];
      double epsilon_ij,sigma_ij;
      if(lja_ij==0.0) 
        epsilon_ij=0.0;
      else 
        epsilon_ij=0.25*ljb_ij*ljb_ij/lja_ij;
      if(ljb_ij==0.0)
        sigma_ij=0.0;
      else
        sigma_ij=pow(lja_ij/ljb_ij,1.0/6); 
      string spcoef="t"+to_string(ii+1)+" t"+to_string(jj+1)+" lj/cut/coul/long ";
      spcoef+=to_string(epsilon_ij)+" "+to_string(sigma_ij)+" "+to_string(LJ_CUT,1);
      Vpaircoeff.push_back(spcoef);
    }
  }
  //cout<<"Vpaircoeff : "<<Vpaircoeff.size()<<endl;

  Vboxsize.clear();
  if(Vbox_dimensions.size()>0) {
    Vboxsize.push_back(0.0);
    Vboxsize.push_back(Vbox_dimensions[1]);
    Vboxsize.push_back(0.0);
    Vboxsize.push_back(Vbox_dimensions[2]);
    Vboxsize.push_back(0.0);
    Vboxsize.push_back(Vbox_dimensions[3]);
  }
  //cout<<"Vboxsize : "<<Vboxsize.size()<<endl;
}

void ReadAmber::read_inpcrd(string fname) {
  if(fname.find(".inpcrd")==string::npos)
    fname+=".inpcrd";
  
  if(!is_file_exist(fname)) {
    string str="read_amber.cpp: Cannot find file "+fname;
    error->one(str);
  }

  if(Vatoms.size()==0) {
    string str="read_amber.cpp: Please read .prmtop first!";
    error->one(str);
  }

  vector<string> Vlines=read_file(fname);

  if(Vlines.size()<2) {
    string str="read_amber.cpp: File "+fname+" is too short!";
    error->one(str);
  }
  
  int Natom=stoi(Vlines[1].substr(0,6));
  if(Natom!=(int)Vatoms.size()) {
    string str="read_amber.cpp: Number of atoms is not matched in .inpcrd and .prmtop files --> ";
    str+=to_string(Natom)+" / "+to_string(Vatoms.size());
    error->one(str);
  }

  vector<string> Vcontent;
  //int repeat_count=6;
  int field_width=12;
  for(size_t ii=2;ii<Vlines.size();ii++) {
    vector<string> Vtmp=ReadAmber::parse_line(Vlines[ii],field_width);
    Vcontent.insert(Vcontent.end(),Vtmp.begin(),Vtmp.end());
  }

  if((int)Vcontent.size()<Natom*3) {
    string str="read_amber.cpp: File "+fname+" is too short!";
    error->one(str);
  }

  for(int ii=0;ii<Natom;ii++) {
    Vatoms[ii].x=stod(Vcontent[3*ii]);
    Vatoms[ii].y=stod(Vcontent[3*ii+1]);
    Vatoms[ii].z=stod(Vcontent[3*ii+2]);
  }
}

void ReadAmber::clear() {
  Vatoms.clear();
  Vatoms.shrink_to_fit();
  Vbonds.clear();
  Vbonds.shrink_to_fit();
  Vangles.clear();
  Vangles.shrink_to_fit();
  Vdihedrals.clear();
  Vdihedrals.shrink_to_fit();
  Vimpropers.clear();
  Vimpropers.shrink_to_fit();
  Vmass.clear();
  Vmass.shrink_to_fit();
  Vbondcoeff.clear();
  Vbondcoeff.shrink_to_fit();
  Vanglecoeff.clear();
  Vanglecoeff.shrink_to_fit();
  Vdihedralcoeff.clear();
  Vdihedralcoeff.shrink_to_fit();
  Vimpropercoeff.clear();
  Vimpropercoeff.shrink_to_fit();
  Vpaircoeff.clear();
  Vpaircoeff.shrink_to_fit();
  Vboxsize.clear();
  Vboxsize.shrink_to_fit();
}

// Function to parse the %FORMAT line to get details
pair<int, int> ReadAmber::parse_format(string sline) {
  // format string looks like "%FORMAT(20a4)"
  sline=sline.substr(8); 
  size_t type_pos=sline.find_first_of("iIFEaADG"); // Find the data type character
  string srep=sline.substr(0,type_pos);
  string swidth=sline.substr(type_pos+1);
  vector<string> stoken=split_string(srep,"().");
  int repeat_count=stoi(stoken[0]);
  stoken.clear();
  stoken=split_string(swidth,"().");
  int field_width=stoi(stoken[0]);

  pair<int, int> p_one=make_pair(repeat_count, field_width);

  return p_one;
}

vector<string> ReadAmber::parse_line(string sline, int field_width) {
  vector<string> Vfields;
  for(size_t ii=0;ii<sline.length();ii+=field_width) {
    // Extract the fixed-width substring
    string field=sline.substr(ii,field_width);
    Vfields.push_back(field);
  }

  return Vfields;
}