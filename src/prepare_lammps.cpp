#include "basic.h"
#include "constant.h"
#include "readfile.h"
#include "initialization.h"
#include "error.h"
#include "read_amber.h"
#include "prepare_lammps.h"

using namespace CBS_NS;

enum{EM,MD,REMD};
const int PrepareLAMMPS::em=EM;
const int PrepareLAMMPS::md=MD;
const int PrepareLAMMPS::remd=REMD;

PrepareLAMMPS::PrepareLAMMPS(Cobmms *cbs) : Pointers(cbs), Rcoeff(cbs) {
  Voutput_logs.clear();
  Voutput_trjs.clear();
  Vbondstyle.clear();
  Vanglestyle.clear();
  Vdihedralstyle.clear();
  Vimproperstyle.clear();
  Vatomtype.clear();
  Vbondtype.clear();
  Vangletype.clear();
  Vdihedraltype.clear();
  Vimpropertype.clear();
  Vboxsize.clear();
}

PrepareLAMMPS::PrepareLAMMPS(const PrepareLAMMPS &other) : Pointers(other.cbs), Rcoeff(cbs) {
  md_label=other.md_label;
  Vtemperature=other.Vtemperature;
  Rcoeff=other.Rcoeff;
  Mrestraint_type=other.Mrestraint_type;
  Vhead=other.Vhead;
  Vio=other.Vio;
  Vgroup=other.Vgroup;
  Vem=other.Vem;
  Vfix=other.Vfix;
  Vrun=other.Vrun;
  Voutput_logs=other.Voutput_logs;
  Voutput_trjs=other.Voutput_trjs;
  Vboxsize=other.Vboxsize;
  Vbondstyle=other.Vbondstyle;
  Vanglestyle=other.Vanglestyle;
  Vdihedralstyle=other.Vdihedralstyle;
  Vimproperstyle=other.Vimproperstyle;
  Vatomtype=other.Vatomtype;
  Vbondtype=other.Vbondtype;
  Vangletype=other.Vangletype;
  Vdihedraltype=other.Vdihedraltype;
  Vimpropertype=other.Vimpropertype;
}

PrepareLAMMPS& PrepareLAMMPS::operator=(const PrepareLAMMPS &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    md_label=other.md_label;
    Vtemperature=other.Vtemperature;
    Rcoeff=other.Rcoeff;
    Mrestraint_type=other.Mrestraint_type;
    Vhead=other.Vhead;
    Vio=other.Vio;
    Vgroup=other.Vgroup;
    Vem=other.Vem;
    Vfix=other.Vfix;
    Vrun=other.Vrun;
    Voutput_logs=other.Voutput_logs;
    Voutput_trjs=other.Voutput_trjs;
    Vboxsize=other.Vboxsize;
    Vbondstyle=other.Vbondstyle;
    Vanglestyle=other.Vanglestyle;
    Vdihedralstyle=other.Vdihedralstyle;
    Vimproperstyle=other.Vimproperstyle;
    Vatomtype=other.Vatomtype;
    Vbondtype=other.Vbondtype;
    Vangletype=other.Vangletype;
    Vdihedraltype=other.Vdihedraltype;
    Vimpropertype=other.Vimpropertype;
  }

  return *this;
}

PrepareLAMMPS::~PrepareLAMMPS() {
  md_label=-1;
  Vtemperature.clear();
  Vtemperature.shrink_to_fit();
  Rcoeff.clear();
  Mrestraint_type.clear();
  Vhead.clear();
  Vhead.shrink_to_fit();
  Vio.clear();
  Vio.shrink_to_fit();
  Vgroup.clear();
  Vgroup.shrink_to_fit();
  Vem.clear();
  Vem.shrink_to_fit();
  Vfix.clear();
  Vfix.shrink_to_fit();
  Vrun.clear();
  Vrun.shrink_to_fit();
  Voutput_logs.clear();
  Voutput_logs.shrink_to_fit();
  Voutput_trjs.clear();
  Voutput_trjs.shrink_to_fit();
  Vboxsize.clear();
  Vboxsize.shrink_to_fit();
  Vbondstyle.clear();
  Vbondstyle.shrink_to_fit();
  Vanglestyle.clear();
  Vanglestyle.shrink_to_fit();
  Vdihedralstyle.clear();
  Vdihedralstyle.shrink_to_fit();
  Vimproperstyle.clear();
  Vimproperstyle.shrink_to_fit();
  Vatomtype.clear();
  Vatomtype.shrink_to_fit();
  Vbondtype.clear();
  Vbondtype.shrink_to_fit();
  Vangletype.clear();
  Vangletype.shrink_to_fit();
  Vdihedraltype.clear();
  Vdihedraltype.shrink_to_fit();
  Vimpropertype.clear();
  Vimpropertype.shrink_to_fit();
}

void PrepareLAMMPS::set_temperature(double T_start, double T_end) {
  Vtemperature.clear();
  Vtemperature.push_back(T_start);
  Vtemperature.push_back(T_end);
  md_label=PrepareLAMMPS::md;
}

void PrepareLAMMPS::set_remd_temperature() {
  Vtemperature.clear();

  //double TEMP[]={250,273,297,322,349,378,408,439};
  //for(int ii=0;ii<get_array_size(TEMP);ii++) {
    //Vtemperature.push_back(TEMP[ii]);
  //}

  double T0=275;
  double k=0.06;
  for(int ii=0;ii < init->Num_remd_dup;ii++) {
    double temp=round(T0*exp(k*ii))*1.0;
    Vtemperature.push_back(temp);
  }

  md_label=PrepareLAMMPS::remd;
}

void PrepareLAMMPS::set_coeff(const ReadCoeff &coeff_one) {
  if(coeff_one.Vmass.size()==0) {
    string str="prepare_lammps.cpp: Please read coeffs first!";
    error->one(str);
  }
  
  Rcoeff.clear();
  Rcoeff=coeff_one;
}

vector<string> PrepareLAMMPS::prepare_data(const Molecule &mol_one, const vector<Bond> &Vbonds, const vector<Angle> &Vangles,
 const vector<Dihedral> &Vdihedrals, const vector<Dihedral> &Vimpropers, double delxyz) {
  vector<Atom> Vatoms=mol_one.get_atom_list();
  vector<string> Vdata=PrepareLAMMPS::prepare_data(Vatoms,Vbonds,Vangles,Vdihedrals,Vimpropers,delxyz);

  return Vdata;
}

vector<string> PrepareLAMMPS::prepare_data(const vector<Atom> &Vatoms_one, const vector<Bond> &Vbonds, const vector<Angle> &Vangles,
 const vector<Dihedral> &Vdihedrals, const vector<Dihedral> &Vimpropers, double delxyz) {
  if(Rcoeff.Vmass.size()==0) {
    string str="prepare_lammps.cpp: Please set coeffs first!";
    error->one(str);
  }
  
  double xlo=MathConst::BIG;
  double xhi=-1.0*MathConst::BIG;
  double ylo=MathConst::BIG;
  double yhi=-1.0*MathConst::BIG;
  double zlo=MathConst::BIG;
  double zhi=-1.0*MathConst::BIG;
 
  int natom=0;
  Vatomtype.clear();
  map<string,int> Mchainid_out;
  for(size_t kk=0;kk<Vatoms_one.size();kk++) {
    natom++;
    Atom atom_one=Vatoms_one[kk];
    if(atom_one.x<xlo) xlo=atom_one.x;
    else if(atom_one.x>xhi) xhi=atom_one.x;
    if(atom_one.y<ylo) ylo=atom_one.y;
    else if(atom_one.y>yhi) yhi=atom_one.y;
    if(atom_one.z<zlo) zlo=atom_one.z;
    else if(atom_one.z>zhi) zhi=atom_one.z;
    string atomname_one=atom_one.atomname;
    if(count(Vatomtype.begin(),Vatomtype.end(),atomname_one)==0)
    Vatomtype.push_back(atomname_one);
    string chainid_one=atom_one.chainid;
    int flag_check=1;
    for(size_t nn=0;nn<chainid_one.length();nn++) {
      if(!isdigit(chainid_one[nn])) {
        flag_check=0;
        break;
      }
    }
    if(Mchainid_out.count(chainid_one)==0) {
      int molid_out;
      if(flag_check)
        molid_out=stoi(chainid_one);
      else 
        molid_out=static_cast<int>(Mchainid_out.size());
      Mchainid_out.emplace(chainid_one,molid_out);
    }
  }
 
  if(delxyz>=0.0 || Vboxsize.size()==0) {
    Vboxsize.clear();
    Vboxsize.push_back(xlo-delxyz);
    Vboxsize.push_back(xhi+delxyz);
    Vboxsize.push_back(ylo-delxyz);
    Vboxsize.push_back(yhi+delxyz);
    Vboxsize.push_back(zlo-delxyz);
    Vboxsize.push_back(zhi+delxyz);
  }
 
  Mrestraint_type.clear();
  map<string,int> Mlabel_type;
  Vbondtype.clear();
  Vbondstyle.clear();
  vector<string> Vbondcoeff_restr;
  int nbond=static_cast<int> (Vbonds.size());
  for(size_t ii=0;ii<Vbonds.size();ii++) {
    if(Vbonds[ii].atom_i<0 || Vbonds[ii].atom_j<0) nbond--;
    vector<string> stoken=split_string(Vbonds[ii].type);
    if(stoken.size()==1) {
      string bondtype_one=stoken[0];
      if(count(Vbondtype.begin(),Vbondtype.end(),bondtype_one)==0) {
        Vbondtype.push_back(bondtype_one);
      }
    }
    else if(stoken.size()>1) {
      string bondstyle=stoken[0];
      string bondcoeff=stoken[1];
      for(size_t jj=2;jj<stoken.size();jj++)
        bondcoeff+=" "+stoken[jj];
      if(count(Vbondstyle.begin(),Vbondstyle.end(),bondstyle)==0)
        Vbondstyle.push_back(bondstyle);
      if(Mrestraint_type.count(Vbonds[ii].type)==0) {
        string bondtype_one="UT"+to_string((int)Mrestraint_type.size()/2);
        Vbondtype.push_back(bondtype_one);
        string str=bondtype_one+" "+bondstyle+" "+bondcoeff;
        Vbondcoeff_restr.push_back(str);
        Mrestraint_type.emplace(Vbonds[ii].type,bondtype_one);
        Mrestraint_type.emplace(bondtype_one,Vbonds[ii].type);
      }
    } 
  }

  vector<string> Vbondcoeff_out;
  for(size_t kk=0;kk<Rcoeff.Vbondcoeff.size();kk++) {
    vector<string> stoken=split_string(Rcoeff.Vbondcoeff[kk]);
    string bondtype_one=stoken[0];
    if(count(Vbondtype.begin(),Vbondtype.end(),bondtype_one)) {
      Vbondcoeff_out.push_back(Rcoeff.Vbondcoeff[kk]);
      string bstyle=stoken[1];
      if(count(Vbondstyle.begin(),Vbondstyle.end(),bstyle)==0)
        Vbondstyle.push_back(bstyle);
    }
  }
  if(Vbondcoeff_restr.size()>0)
    Vbondcoeff_out.insert(Vbondcoeff_out.end(),Vbondcoeff_restr.begin(),Vbondcoeff_restr.end());
        
  if(Vbondtype.size()!=Vbondcoeff_out.size()) {
    string sbtype="";
    for(size_t jj=0;jj<Vbondtype.size();jj++) {
      int flag_check=0;
      for(size_t kk=0;kk<Vbondcoeff_out.size();kk++) {
        vector<string> stoken=split_string(Vbondcoeff_out[kk]);
        if(stoken[0]==Vbondtype[jj]) {
          flag_check=1;
          break;
        }
      }
      if(flag_check==0) sbtype+=" "+Vbondtype[jj];
    }
    if(sbtype.size()>0) {
      string str="prepare_lammps.cpp: Cannot find corresponding coeff for Bond type "+sbtype;
      error->one(str);
    }
  }

  if(Vbondstyle.size()==1) {
    for(size_t jj=0;jj<Vbondcoeff_out.size();jj++) {
      vector<string> stoken=split_string(Vbondcoeff_out[jj]);
      string bstr=stoken[0];
      for(size_t kk=2;kk<stoken.size();kk++)
        bstr+=" "+stoken[kk];
      Vbondcoeff_out[jj]=bstr;
    }
  }
  //cout<<"Restraints: "<<Mrestraint_type.size()/2<<endl;
  //cout<<"Bond type: "<<Vbondtype.size()<<endl;
 
  Vangletype.clear();
  Vanglestyle.clear();
  for(size_t ii=0;ii<Vangles.size();ii++) {
    vector<string> stoken=split_string(Vangles[ii].type);
    string angletype_one=stoken[0];
    if(count(Vangletype.begin(),Vangletype.end(),angletype_one)==0) {
      Vangletype.push_back(angletype_one);
    }
  }

  vector<string> Vanglecoeff_out;
  for(size_t kk=0;kk<Rcoeff.Vanglecoeff.size();kk++) {
    vector<string> stoken=split_string(Rcoeff.Vanglecoeff[kk]);
    string angletype_one=stoken[0];
    if(count(Vangletype.begin(),Vangletype.end(),angletype_one)) {
      Vanglecoeff_out.push_back(Rcoeff.Vanglecoeff[kk]);
      string astyle=stoken[1];
      if(count(Vanglestyle.begin(),Vanglestyle.end(),astyle)==0)
        Vanglestyle.push_back(astyle);
    }
  }

  if(Vangletype.size()!=Vanglecoeff_out.size()) {
    string satype="";
    for(size_t jj=0;jj<Vangletype.size();jj++) {
      int flag_check=0;
      for(size_t kk=0;kk<Vanglecoeff_out.size();kk++) {
        vector<string> stoken=split_string(Vanglecoeff_out[kk]);
        if(stoken[0]==Vangletype[jj]) {
          flag_check=1;
          break;
        }
      }
      if(flag_check==0) satype+=" "+Vangletype[jj];
    }
    if(satype.size()>0) {
      string str="prepare_lammps.cpp: Cannot find corresponding coeff for Angle type "+satype;
      error->one(str);
    }
  }
  
  if(Vanglestyle.size()==1) {
    for(size_t jj=0;jj<Vanglecoeff_out.size();jj++) {
      vector<string> stoken=split_string(Vanglecoeff_out[jj]);
      string astr=stoken[0];
      for(size_t kk=2;kk<stoken.size();kk++)
        astr+=" "+stoken[kk];
      Vanglecoeff_out[jj]=astr;
    }
  } 
  //cout<<"Angle type: "<<Vangletype.size()<<endl;
   
  Vdihedraltype.clear();
  Vdihedralstyle.clear();
  for(size_t ii=0;ii<Vdihedrals.size();ii++) {
    vector<string> stoken=split_string(Vdihedrals[ii].type);
    string dihedraltype_one=stoken[0];
    if(count(Vdihedraltype.begin(),Vdihedraltype.end(),dihedraltype_one)==0) {
      Vdihedraltype.push_back(dihedraltype_one);
    }  
  }

  vector<string> Vdihedralcoeff_out;
  for(size_t kk=0;kk<Rcoeff.Vdihedralcoeff.size();kk++) {
    vector<string> stoken=split_string(Rcoeff.Vdihedralcoeff[kk]);
    string dihedraltype_one=stoken[0];
    if(count(Vdihedraltype.begin(),Vdihedraltype.end(),dihedraltype_one)) {
      Vdihedralcoeff_out.push_back(Rcoeff.Vdihedralcoeff[kk]);
      string dstyle=stoken[1];
      if(count(Vdihedralstyle.begin(),Vdihedralstyle.end(),dstyle)==0) 
        Vdihedralstyle.push_back(dstyle);
    }
  }

  if(Vdihedraltype.size()!=Vdihedralcoeff_out.size()) {
    string sdtype="";
    for(size_t jj=0;jj<Vdihedraltype.size();jj++) {
      int flag_check=0;
      for(size_t kk=0;kk<Vdihedralcoeff_out.size();kk++) {
        vector<string> stoken=split_string(Vdihedralcoeff_out[kk]);
        if(stoken[0]==Vdihedraltype[jj]) {
          flag_check=1;
          break;
        }
      }
      if(flag_check==0) sdtype+=" "+Vdihedraltype[jj];
    }
    if(sdtype.size()>0) {
      string str="prepare_lammps.cpp: Cannot find corresponding coeff for Dihedral type "+sdtype;
      error->one(str);
    }
  }
  
  if(Vdihedralstyle.size()==1) {
    for(size_t jj=0;jj<Vdihedralcoeff_out.size();jj++) {
      vector<string> stoken=split_string(Vdihedralcoeff_out[jj]);
      string dstr=stoken[0];
      for(size_t kk=2;kk<stoken.size();kk++) 
        dstr+=" "+stoken[kk];
      Vdihedralcoeff_out[jj]=dstr;
    }
  }
  //cout<<"Dihedral type: "<<Vdihedraltype.size()<<endl;
 
  Vimpropertype.clear();
  Vimproperstyle.clear();
  for(size_t ii=0;ii<Vimpropers.size();ii++) {
    vector<string> stoken=split_string(Vimpropers[ii].type);
    string impropertype_one=stoken[0];
    if(count(Vimpropertype.begin(),Vimpropertype.end(),impropertype_one)==0) {
      Vimpropertype.push_back(impropertype_one);
    }
  }

  vector<string> Vimpropercoeff_out;
  for(size_t kk=0;kk<Rcoeff.Vimpropercoeff.size();kk++) {
    vector<string> stoken=split_string(Rcoeff.Vimpropercoeff[kk]);
    string impropertype_one=stoken[0];
    if(count(Vimpropertype.begin(),Vimpropertype.end(),impropertype_one)) {
      Vimpropercoeff_out.push_back(Rcoeff.Vimpropercoeff[kk]);
      string istyle=stoken[1];
      if(count(Vimproperstyle.begin(),Vimproperstyle.end(),istyle)==0)
        Vimproperstyle.push_back(istyle);
    }
  }

  if(Vimpropertype.size()!=Vimpropercoeff_out.size()) {
    string sitype="";
    for(size_t jj=0;jj<Vimpropertype.size();jj++) {
      int flag_check=0;
      for(size_t kk=0;kk<Vimpropercoeff_out.size();kk++) {
        vector<string> stoken=split_string(Vimpropercoeff_out[kk]);
        if(stoken[0]==Vimpropertype[jj]) {
          flag_check=1;
          break;
        }
      }
      if(flag_check==0) sitype+=" "+Vimpropertype[jj];
    }
    if(sitype.size()>0) {
      string str="prepare_lammps.cpp: Cannot find corresponding coeff for Improper type "+sitype;
      error->one(str);
    }
  }
  
  if(Vimproperstyle.size()==1) {
    for(size_t jj=0;jj<Vimpropercoeff_out.size();jj++) {
      vector<string> stoken=split_string(Vimpropercoeff_out[jj]);
      string istr=stoken[0];
      for(size_t kk=2;kk<stoken.size();kk++) 
        istr+=" "+stoken[kk];
      Vimpropercoeff_out[jj]=istr;
    }
  }
  //cout<<"Improper type: "<<Vimpropertype.size()<<endl;
  
  vector<string> Vdata;
  string sline;
  if(Vatomtype.size()!=Rcoeff.Vmass.size()) {
    sline="#Number of atom types is less than default --> "+to_string(Vatomtype.size())+"/"+to_string(Rcoeff.Vmass.size());
    Vdata.push_back(sline);
  }
  sline="#RNA CG Model";
  Vdata.push_back(sline);
  sline=to_string(natom)+" atoms";
  Vdata.push_back(sline);
  sline=to_string(nbond)+" bonds";
  Vdata.push_back(sline);
  sline=to_string(Vangles.size())+" angles";
  Vdata.push_back(sline);
  sline=to_string(Vdihedrals.size())+" dihedrals";
  Vdata.push_back(sline);
  sline=to_string(Vimpropers.size())+" impropers";
  Vdata.push_back(sline);
  sline=to_string(Rcoeff.Vmass.size())+" atom types";
  Vdata.push_back(sline);
  sline=to_string(Vbondtype.size())+" bond types";
  Vdata.push_back(sline);
  sline=to_string(Vangletype.size())+" angle types";
  Vdata.push_back(sline);
  sline=to_string(Vdihedraltype.size())+" dihedral types";
  Vdata.push_back(sline);
  sline=to_string(Vimpropertype.size())+" improper types";
  Vdata.push_back(sline);
  ostringstream oss;
  oss<<std::fixed<<std::setprecision(4)<<Vboxsize[0]<<"  "<<std::fixed<<std::setprecision(4)<<Vboxsize[1];
  sline=oss.str()+" xlo xhi";
  Vdata.push_back(sline);
  oss.str("");
  oss.clear();
  oss<<std::fixed<<std::setprecision(4)<<Vboxsize[2]<<"  "<<std::fixed<<std::setprecision(4)<<Vboxsize[3];
  sline=oss.str()+" ylo yhi";
  Vdata.push_back(sline);
  oss.str("");
  oss.clear();
  oss<<std::fixed<<std::setprecision(4)<<Vboxsize[4]<<"  "<<std::fixed<<std::setprecision(4)<<Vboxsize[5];
  sline=oss.str()+" zlo zhi";
  Vdata.push_back(sline);
  sline="";
  Vdata.push_back(sline);
  sline="Atom Type Labels";
  Vdata.push_back(sline);
  sline="";
  Vdata.push_back(sline);
  int type_one=0;
  for(size_t jj=0;jj<Rcoeff.Vmass.size();jj++) {
    vector<string> stoken=split_string(Rcoeff.Vmass[jj]);
    string atomtype=stoken[0];
    //if(count(Vatomtype.begin(),Vatomtype.end(),atomtype)) {
      type_one++;
      sline=to_string(type_one)+" "+atomtype;
      Vdata.push_back(sline);
      Mlabel_type.emplace(atomtype,type_one);
    //}
  }
  sline="";
  Vdata.push_back(sline);
  sline="Atoms";
  Vdata.push_back(sline);
  sline="";
  Vdata.push_back(sline);
  for(size_t ii=0;ii<Vatoms_one.size();ii++) {
    Atom atom_one=Vatoms_one[ii];
    oss.str("");
     oss.clear();
    //atomid chainid atomtype q x y z
    size_t molid_out=Mchainid_out[atom_one.chainid];
    oss<<atom_one.id<<" "<<molid_out<<" "<<atom_one.atomname<<" ";
    if(Rcoeff.Mtype_charge.count(atom_one.atomname))
      oss<<Rcoeff.Mtype_charge[atom_one.atomname];
    else {
      oss<<std::fixed<<std::setprecision(4)<<atom_one.q;
    }
    oss<<" "<<std::fixed<<std::setprecision(4)<<atom_one.x<<" "<<std::fixed<<std::setprecision(4)<<atom_one.y;
    oss<<" "<<std::fixed<<std::setprecision(4)<<atom_one.z;                      
    sline=oss.str(); 
    Vdata.push_back(sline); 
  }
  if(Vbondtype.size()>0) {
    sline="";
    Vdata.push_back(sline);
    sline="Bond Type Labels";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    type_one=0;
    for(size_t jj=0;jj<Vbondcoeff_out.size();jj++) {
      vector<string> stoken=split_string(Vbondcoeff_out[jj]);
      string bondtype=stoken[0];
      if(count(Vbondtype.begin(),Vbondtype.end(),bondtype)) {
        type_one++;
        sline=to_string(type_one)+" "+bondtype;
        Vdata.push_back(sline);
        Mlabel_type.emplace(bondtype,type_one);
      }
    }
    sline="";
    Vdata.push_back(sline);
    sline="Bonds";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    int nindex=0;
    for(size_t jj=0;jj<Vbonds.size();jj++) {
      if(Vbonds[jj].atom_i<0 || Vbonds[jj].atom_j<0) continue;
      nindex++;
      sline=to_string(nindex);
      string bondtype_one=Vbonds[jj].type;
      if(Mrestraint_type.count(bondtype_one))
       sline+=" "+Mrestraint_type[bondtype_one];
      else 
        sline+=" "+bondtype_one;
      sline+=" "+to_string(Vbonds[jj].atom_i)+" "+to_string(Vbonds[jj].atom_j);
      Vdata.push_back(sline);
    }
  }
  if(Vangletype.size()>0) {
    sline="";
    Vdata.push_back(sline);
    sline="Angle Type Labels";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    type_one=0;
    for(size_t jj=0;jj<Rcoeff.Vanglecoeff.size();jj++) {
      vector<string> stoken=split_string(Rcoeff.Vanglecoeff[jj]);
      string angletype=stoken[0];
      if(count(Vangletype.begin(),Vangletype.end(),angletype)) {
        type_one++;
        sline=to_string(type_one)+" "+angletype;
        Vdata.push_back(sline);
        Mlabel_type.emplace(angletype,type_one);
      }
    }
    sline="";
    Vdata.push_back(sline);
    sline="Angles";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    for(size_t jj=0;jj<Vangles.size();jj++) {
      sline=to_string(jj+1)+" "+Vangles[jj].type+" "+to_string(Vangles[jj].atom_i);
      sline+=" "+to_string(Vangles[jj].atom_j)+" "+to_string(Vangles[jj].atom_k);
      Vdata.push_back(sline);
    }
  }
  if(Vdihedraltype.size()>0) {
    sline="";
    Vdata.push_back(sline);
    sline="Dihedral Type Labels";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    type_one=0;
    for(size_t jj=0;jj<Rcoeff.Vdihedralcoeff.size();jj++) {
      vector<string> stoken=split_string(Rcoeff.Vdihedralcoeff[jj]);
      string dihedraltype=stoken[0];
      if(count(Vdihedraltype.begin(),Vdihedraltype.end(),dihedraltype)) {
        type_one++;
        sline=to_string(type_one)+" "+dihedraltype;
        Vdata.push_back(sline);
        Mlabel_type.emplace(dihedraltype,type_one);
      }
    }
    sline="";
    Vdata.push_back(sline);
    sline="Dihedrals";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    for(size_t jj=0;jj<Vdihedrals.size();jj++) {
      sline=to_string(jj+1)+" "+Vdihedrals[jj].type+" "+to_string(Vdihedrals[jj].atom_i);
      sline+=" "+to_string(Vdihedrals[jj].atom_j)+" "+to_string(Vdihedrals[jj].atom_k)+" "+to_string(Vdihedrals[jj].atom_l);
      Vdata.push_back(sline);
    }
  }
  if(Vimpropertype.size()>0) {
    sline="";
    Vdata.push_back(sline);
    sline="Improper Type Labels";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    type_one=0;
    for(size_t jj=0;jj<Rcoeff.Vimpropercoeff.size();jj++) {
      vector<string> stoken=split_string(Rcoeff.Vimpropercoeff[jj]);
      string impropertype=stoken[0];
      if(count(Vimpropertype.begin(),Vimpropertype.end(),impropertype)) {
        type_one++;
        sline=to_string(type_one)+" "+impropertype;
        Vdata.push_back(sline);
        Mlabel_type.emplace(impropertype,type_one);
      }
    }
    sline="";
    Vdata.push_back(sline);
    sline="Impropers";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    for(size_t jj=0;jj<Vimpropers.size();jj++) {
      sline=to_string(jj+1)+" "+Vimpropers[jj].type+" "+to_string(Vimpropers[jj].atom_i);
      sline+=" "+to_string(Vimpropers[jj].atom_j)+" "+to_string(Vimpropers[jj].atom_k)+" "+to_string(Vimpropers[jj].atom_l);
      Vdata.push_back(sline);
    }
  }
  
  if(Vbondtype.size()>0) {
    sline="";
    Vdata.push_back(sline);
    sline="Bond Coeffs";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    size_t btype=0;
    for(size_t jj=0;jj<Vbondcoeff_out.size();jj++) {
      vector<string> stoken=split_string(Vbondcoeff_out[jj]);
      if(count(Vbondtype.begin(),Vbondtype.end(),stoken[0])) {
        sline=to_string(Mlabel_type[stoken[0]]);
        for(size_t kk=1;kk<stoken.size();kk++)
          sline+=" "+stoken[kk];
        Vdata.push_back(sline);
        btype+=1;
      }
    }
    if(btype!=Vbondtype.size()) {
      string str="prepare_lammps.cpp: Please specify coeffs for all bond type --> "+to_string(Vbondtype.size()-btype)+" type(s) unsigned!";
      error->one(str);
    }
  }
  if(Vangletype.size()>0) {
    sline="";
    Vdata.push_back(sline);
    sline="Angle Coeffs";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    size_t atype=0;
    for(size_t jj=0;jj<Vanglecoeff_out.size();jj++) {
      vector<string> stoken=split_string(Vanglecoeff_out[jj]);
      if(count(Vangletype.begin(),Vangletype.end(),stoken[0])) {
        sline=to_string(Mlabel_type[stoken[0]]);
        for(size_t kk=1;kk<stoken.size();kk++)
          sline+=" "+stoken[kk];
        Vdata.push_back(sline);
        atype+=1;
      }
    }
    if(atype!=Vangletype.size()) {
      string str="prepare_lammps.cpp: Please specify coeffs for all angle type --> "+to_string(Vangletype.size()-atype)+" type(s) unsigned!";
      error->one(str);
    }
  }
  if(Vdihedraltype.size()>0) {
    sline="";
    Vdata.push_back(sline);
    sline="Dihedral Coeffs";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    size_t dtype=0;
    for(size_t jj=0;jj<Vdihedralcoeff_out.size();jj++) {
      vector<string> stoken=split_string(Vdihedralcoeff_out[jj]);
      if(count(Vdihedraltype.begin(),Vdihedraltype.end(),stoken[0])) {
        sline=to_string(Mlabel_type[stoken[0]]);
        for(size_t kk=1;kk<stoken.size();kk++)
          sline+=" "+stoken[kk];
        Vdata.push_back(sline);
        dtype+=1;
      }
    }
    if(dtype!=Vdihedraltype.size()) {
      string str="prepare_lammps.cpp: Please specify coeffs for all dihedral type --> "+to_string(Vdihedraltype.size()-dtype)+" type(s) unsigned!";
      error->one(str);
    }
  }
  if(Vimpropertype.size()>0) {
    sline="";
    Vdata.push_back(sline);
    sline="Improper Coeffs";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    size_t dtype=0;
    for(size_t jj=0;jj<Vimpropercoeff_out.size();jj++) {
      vector<string> stoken=split_string(Vimpropercoeff_out[jj]);
      if(count(Vimpropertype.begin(),Vimpropertype.end(),stoken[0])) {
        sline=to_string(Mlabel_type[stoken[0]]);
        for(size_t kk=1;kk<stoken.size();kk++)
          sline+=" "+stoken[kk];
        Vdata.push_back(sline);
        dtype+=1;
      }
    }
    if(dtype!=Vimpropertype.size()) {
      string str="prepare_lammps.cpp: Please specify coeffs for all improper type --> "+to_string(Vimpropertype.size()-dtype)+" type(s) unsigned!";
      error->one(str);
    }
  }
  
  if(Rcoeff.Vpaircoeff.size()==Vatomtype.size()) {
    sline="";
    Vdata.push_back(sline);
    sline="Pair Coeffs";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    size_t npair_style=Rcoeff.Mpair_style_args.size();
    if(Rcoeff.Mpair_style_args.count("neigh"))
      npair_style=npair_style-1;
    size_t index_begin=1;
    if(npair_style==1)
      index_begin=2;
    for(size_t jj=0;jj<Rcoeff.Vpaircoeff.size();jj++) {
      vector<string> stoken=split_string(Rcoeff.Vpaircoeff[jj]);
      if(Mlabel_type.count(stoken[0])==0) {
        string str="prepare_lammps.cpp: Unknown atom type for pair atoms --> "+stoken[0];
        error->warning(str);
        continue;
      }
      int atom_i_type=Mlabel_type[stoken[0]];
      sline=to_string(atom_i_type);
      for(size_t kk=index_begin;kk<stoken.size();kk++)
        sline+=" "+stoken[kk];
      Vdata.push_back(sline);
    }
  }
  else {
    sline="";
    Vdata.push_back(sline);
    sline="PairIJ Coeffs";
    Vdata.push_back(sline);
    sline="";
    Vdata.push_back(sline);
    size_t npair_style=Rcoeff.Mpair_style_args.size();
    if(Rcoeff.Mpair_style_args.count("neigh"))
      npair_style=npair_style-1;
    size_t index_begin=2;
    if(npair_style==1)
      index_begin=3;
    for(size_t jj=0;jj<Rcoeff.Vpaircoeff.size();jj++) {
      vector<string> stoken=split_string(Rcoeff.Vpaircoeff[jj]);
      if(Mlabel_type.count(stoken[0])==0 || Mlabel_type.count(stoken[1])==0) {
        string str="prepare_lammps.cpp: Unknown atom type for pair atoms --> "+stoken[0]+" "+stoken[1];
        error->warning(str);
        continue;
      }
      int atom_i_type=Mlabel_type[stoken[0]];
      int atom_j_type=Mlabel_type[stoken[1]];
      if(stoken[2].find(Rcoeff.default_pair_style)!=string::npos) {
        sline=to_string(atom_i_type)+" "+to_string(atom_j_type);
        for(size_t kk=index_begin;kk<stoken.size();kk++)
          sline+=" "+stoken[kk];
        Vdata.push_back(sline);
      }
    }
  }
  
  sline="";
  Vdata.push_back(sline);
  sline="Masses";
  Vdata.push_back(sline);
  sline="";
  Vdata.push_back(sline);
  for(size_t jj=0;jj<Rcoeff.Vmass.size();jj++) {
    vector<string> stoken=split_string(Rcoeff.Vmass[jj]);
    //if(count(Vatomtype.begin(),Vatomtype.end(),stoken[0])) {
      sline=Rcoeff.Vmass[jj];
      Vdata.push_back(sline);
    //}
  }
 
  return Vdata;
}

vector<string> PrepareLAMMPS::prepare_data_from_amber(const Molecule &mol_one, string fname) {
  ReadAmber readamb_one(cbs);
  readamb_one.prepare_amber_files(mol_one,fname);
  readamb_one.read_prmtop(fname);
  readamb_one.read_inpcrd(fname);
  Rcoeff.clear();
  Rcoeff.Vmass=readamb_one.Vmass;
  Rcoeff.Vbondcoeff=readamb_one.Vbondcoeff;
  Rcoeff.Vanglecoeff=readamb_one.Vanglecoeff;
  Rcoeff.Vdihedralcoeff=readamb_one.Vdihedralcoeff;
  Rcoeff.Vimpropercoeff=readamb_one.Vimpropercoeff;
  Rcoeff.Vpaircoeff=readamb_one.Vpaircoeff;
  Rcoeff.default_pair_style="lj/cut/coul/long";
  Vboxsize=readamb_one.Vboxsize;

  Rcoeff.get_interaction_style();

  vector<string> Vdata=PrepareLAMMPS::prepare_data(readamb_one.Vatoms,readamb_one.Vbonds,readamb_one.Vangles,
                                                   readamb_one.Vdihedrals,readamb_one.Vimpropers,-1.0);
  
  return Vdata;
}

void PrepareLAMMPS::update_data_xyz(vector<string> &Vdata, const vector<Atom> &Vatoms_xyz) {
  map<int,string> Matom_xyz;
  for(size_t ii=0;ii<Vatoms_xyz.size();ii++) {
    int atomid=Vatoms_xyz[ii].id;
    ostringstream oss;
    oss<<" "<<std::fixed<<std::setprecision(4)<<Vatoms_xyz[ii].x<<" "<<std::fixed<<std::setprecision(4)<<Vatoms_xyz[ii].y;
    oss<<" "<<std::fixed<<std::setprecision(4)<<Vatoms_xyz[ii].z;
    string sxyz=oss.str();
    oss.clear();
    oss.str("");
    if(Matom_xyz.count(atomid)==0)
      Matom_xyz.emplace(atomid,sxyz);
  }

  int natoms=0;
  for(size_t ii=0;ii<Vdata.size();ii++) {
    string sline=Vdata[ii];
    if(sline.size()<2) continue;
    vector<string> stoken=split_string(sline);
    if(stoken.size()==2 && stoken[1]=="atoms") {
      natoms=stoi(stoken[0]);
      break;
    }
  }

  if(natoms!=(int)Matom_xyz.size()) {
    string str="prepare_lammps.cpp: Number of atoms is not matched --> ";
    str+=to_string(natoms)+" / "+to_string(Matom_xyz.size());
    error->one(str);
  }

  int flag_check=0;
  int ncount=0;
  for(size_t ii=0;ii<Vdata.size();ii++) {
    if(Vdata[ii].size()<1) continue;
    if(Vdata[ii].find("Atoms")!=string::npos) {
      flag_check=1;
      continue;
    }
    if(ncount==natoms) {
      flag_check=0;
      break;
    }
    if(flag_check) {
      vector<string> stoken=split_string(Vdata[ii]);
      int atomid=stoi(stoken[0]);
      string sline=stoken[0];
      for(size_t kk=1;kk<stoken.size()-3;kk++)
        sline+=" "+stoken[kk];
      if(Matom_xyz.count(atomid)==0) {
        string str="prepare_lammps.cpp: Cannot find xyz for atom  --> "+Vdata[ii];
        error->one(str);
      }
      sline+=Matom_xyz[atomid];
      Vdata[ii]=sline;
      ncount++;
    } 
  }
}

void PrepareLAMMPS::prepare_in_script(string outDIR, string jobname, string nmodel, int Nstep, int Nevery){
  if(Rcoeff.Vmass.size()==0) {
    string str="prepare_lammps.cpp: Please set coeffs first!";
    error->one(str);
  }
  
  Vhead.clear();
  Vio.clear();
  Vgroup.clear();
  Vem.clear();
  Vfix.clear();
  Vrun.clear();

  string sline="#LAMMPS Script for CGMD simulation";
  Vhead.push_back(sline);
  sline="units           real";
  Vhead.push_back(sline);
  sline="dimension       3";
  Vhead.push_back(sline);
  sline="boundary        p p p";
  Vhead.push_back(sline);
  sline="atom_style      full";
  Vhead.push_back(sline);
  if(init->RNA_version=="latest") {
    sline="special_bonds   lj/coul 0.0 0.0 0.0 angle yes dihedral yes";
    Vhead.push_back(sline);
  }
  else if(init->RNA_version=="v2") {
    sline="special_bonds   lj/coul 0.0 0.0 0.0 angle no dihedral no";
    Vhead.push_back(sline);
  }
  sline="";
  Vhead.push_back(sline);
  sline="neighbor        2.0 bin"; //2.0 multi
  Vhead.push_back(sline);
  sline="neigh_modify    every 10 check yes"; //collection/interval 3 8.0 12.0 "+Mpair_style_args["neigh"];
  Vhead.push_back(sline);
  sline="";
  Vhead.push_back(sline);
  if(Vbondstyle.size()==1) {
    sline="bond_style      "+Vbondstyle[0];
    Vhead.push_back(sline);
  } 
  else if(Vbondstyle.size()>1) {
    sline="bond_style      hybrid";
    for(size_t kk=0;kk<Vbondstyle.size();kk++)
      sline+=" "+Vbondstyle[kk];
    Vhead.push_back(sline);
  }
  if(Vanglestyle.size()==1) {
    sline="angle_style     "+Vanglestyle[0];
    Vhead.push_back(sline);
  }
  else if(Vanglestyle.size()>1) {
    sline="angle_style     hybrid";
    for(size_t kk=0;kk<Vanglestyle.size();kk++) {
      sline+=" "+Vanglestyle[kk];
      if(Vanglestyle[kk]=="zero")
        sline+=" nocoeff";
    }
    Vhead.push_back(sline);
  }
  if(Vdihedralstyle.size()==1) {
    sline="dihedral_style  "+Vdihedralstyle[0];
    Vhead.push_back(sline);
  }
  else if(Vdihedralstyle.size()>1) {
    sline="dihedral_style  hybrid";
    for(size_t kk=0;kk<Vdihedralstyle.size();kk++) {
      sline+=" "+Vdihedralstyle[kk];
      if(Vdihedralstyle[kk]=="zero")
        sline+=" nocoeff";
    }
    Vhead.push_back(sline);
  }
  if(Vimproperstyle.size()==1) {
    sline="improper_style  "+Vimproperstyle[0];
    Vhead.push_back(sline);
  }
  else if(Vimproperstyle.size()>1) {
    sline="improper_style  hybrid";
    for(size_t kk=0;kk<Vimproperstyle.size();kk++) {
      sline+=" "+Vimproperstyle[kk];
      if(Vimproperstyle[kk]=="zero")
        sline+=" nocoeff";
    }
    Vhead.push_back(sline);
  }
  size_t npair_style=Rcoeff.Mpair_style_args.size();
  if(Rcoeff.Mpair_style_args.count("neigh"))
    npair_style=npair_style-1;
  if(npair_style==1) {
    sline="pair_style     ";
  }
  else {
    sline="pair_style      hybrid/overlay";
  }
  for(auto & it : Rcoeff.Mpair_style_args) {
    if(it.first=="neigh") continue;
    sline+=" "+it.first+" "+it.second;
  }
  Vhead.push_back(sline);
  sline="pair_modify     shift yes";
  Vhead.push_back(sline);
  sline="";
  Vhead.push_back(sline);

  sline="variable        t world";
  if(md_label==PrepareLAMMPS::remd) {
    ostringstream oss;
    for(size_t kk=0;kk<Vtemperature.size();kk++)
      oss<<" "<<std::fixed<<std::setprecision(0)<<Vtemperature[kk];
    sline+=oss.str();
    oss.str("");
    oss.clear();
  }
  else {
    sline+=" "+nmodel;
  }
  Vio.push_back(sline);
  sline="read_data       "+outDIR+"data."+jobname+".";
  if(md_label==PrepareLAMMPS::remd) {
    sline+=nmodel;
  }
  else
    sline+="$t";
  Vio.push_back(sline);
  sline="";
  Vio.push_back(sline);
  if(Rcoeff.Vpaircoeff.size()>Rcoeff.Vmass.size()) {
    size_t index_begin=2;
    if(npair_style==1)
      index_begin=3;
    for(size_t jj=0;jj<Rcoeff.Vpaircoeff.size();jj++) {
      vector<string> stoken=split_string(Rcoeff.Vpaircoeff[jj]);
      if(stoken[2].find(Rcoeff.default_pair_style)!=string::npos) continue;
      string sout=stoken[0]+" "+stoken[1];
      for(size_t kk=index_begin;kk<stoken.size();kk++)
        sout+=" "+stoken[kk];
      sline="pair_coeff      "+sout;
      Vio.push_back(sline);
    }
  }
  sline="dielectric      78.285";
  Vio.push_back(sline);
  sline="";
  Vio.push_back(sline);
  sline="thermo          "+to_string(Nevery);
  Vio.push_back(sline);
  sline="log             "+outDIR+"log."+jobname+".";
  if(md_label==PrepareLAMMPS::remd) {
    sline+=nmodel+".$t";
    for(size_t kk=0;kk<Vtemperature.size();kk++) {
      string ostr=outDIR+"log."+jobname+"."+nmodel+".";
      ostringstream oss;
      oss<<std::fixed<<std::setprecision(0)<<Vtemperature[kk];
      ostr+=oss.str();
      Voutput_logs.push_back(ostr);
      oss.str("");
      oss.clear();
    }
  }
  else {
    sline+="$t";
    vector<string> stoken=split_string(nmodel);
    for(size_t kk=0;kk<stoken.size();kk++) {
      string ostr=outDIR+"log."+jobname+".";
      ostr+=stoken[kk];
      Voutput_logs.push_back(ostr);
    }
  }
  Vio.push_back(sline);
  sline="dump            1 all custom "+to_string(Nevery)+" "+outDIR+"dump."+jobname+".";
  if(md_label==PrepareLAMMPS::remd) {
    sline+=nmodel+".";
    for(size_t kk=0;kk<Vtemperature.size();kk++) {
      string ostr=outDIR+"dump."+jobname+"."+nmodel+".";
      ostringstream oss;
      oss<<std::fixed<<std::setprecision(0)<<Vtemperature[kk];
      ostr+=oss.str();
      ostr+=".lammpstrj";
      Voutput_trjs.push_back(ostr);
      oss.str("");
      oss.clear();
    }
  }
  else {
    vector<string> stoken=split_string(nmodel);
    for(size_t kk=0;kk<stoken.size();kk++) {
      string ostr=outDIR+"dump."+jobname+".";
      ostr+=stoken[kk];
      ostr+=".lammpstrj";
      Voutput_trjs.push_back(ostr);
    }
  }
  sline+="$t.lammpstrj id mol type xu yu zu";
  Vio.push_back(sline);
  sline="";
  Vio.push_back(sline);
  sline="restart         "+to_string(Nevery)+" "+outDIR+"restart."+jobname+".";
  if(md_label==PrepareLAMMPS::remd) {
    sline+=nmodel+".$t.1";
    sline+=" "+outDIR+"restart."+jobname+"."+nmodel+".$t.2";
  }
  else {
    sline+="$t.1";
    sline+=" "+outDIR+"restart."+jobname+".$t.2";
  }
  //Vio.push_back(sline);

  sline="group           stat id 1";
  Vgroup.push_back(sline);
  sline="group           other subtract all stat";
  Vgroup.push_back(sline);

  
  sline="";
  Vem.push_back(sline);
  sline="min_style       cg";
  Vem.push_back(sline);
  sline="minimize        1.0e-4 1.0e-6 20000 100000";
  Vem.push_back(sline);
  sline="";
  Vem.push_back(sline);
  
  if(md_label==PrepareLAMMPS::remd) {
    ostringstream oss;
    oss<<"fix             mylang other langevin $t $t 100.0 "<<57412+(int)(743891*rand_one->uniform());
    sline=oss.str();
    Vfix.push_back(sline);
  }
  else if(md_label==PrepareLAMMPS::md) {
    ostringstream oss;
    oss<<"fix             mylang other langevin "<<std::fixed<<std::setprecision(2)<<Vtemperature[0]<<" ";
    oss<<std::fixed<<std::setprecision(2)<<Vtemperature[1]<<" 100.0 "<<57412+(int)(743891*rand_one->uniform());
    sline=oss.str();
    Vfix.push_back(sline);
  }
  sline="fix             mynve other nve/limit 0.2";
  Vfix.push_back(sline);
  sline="";
  Vfix.push_back(sline);
  sline="thermo_style    custom step temp epair emol etotal press f_mynve";
  Vfix.push_back(sline);

  if(md_label==PrepareLAMMPS::remd) {
    ostringstream oss;
    oss<<"velocity        other create $t "<<37624+(int)(913240*rand_one->uniform())<<" rot yes dist gaussian";
    sline=oss.str();
    Vrun.push_back(sline);
  }
  else if(md_label==PrepareLAMMPS::md) {
    ostringstream oss;
    oss<<"velocity        other create "<<std::fixed<<std::setprecision(2)<<Vtemperature[0]<<" ";
    oss<<37624+(int)(913240*rand_one->uniform())<<" rot yes dist gaussian";
    sline=oss.str();
    Vrun.push_back(sline);
  }
  sline="timestep        0.1";
  Vrun.push_back(sline);
  sline="run             50000 upto";
  Vrun.push_back(sline);

  sline="";
  Vrun.push_back(sline);
  /*
  string SHAKETYPE[]={"S-S2","RC-AN","RC-AC","RC-GO","RC-GN","YC-YO","YC-CN","YC-UO","AN-RC-AC","GO-RC-GN","YO-YC-CN","YO-YC-UO"};
  string btype="";
  for(size_t jj=0;jj<Vbondtype.size();jj++) {
    if(count(begin(SHAKETYPE),end(SHAKETYPE),Vbondtype[jj])) {
      btype+=" "+Vbondtype[jj];
    }
  }
  string atype="";
  for(size_t jj=0;jj<Vangletype.size();jj++) {
    if(count(begin(SHAKETYPE),end(SHAKETYPE),Vangletype[jj])) {
      atype+=" "+Vangletype[jj];
    }
  }
  sline="fix             myshake other shake 0.01 100 0 b"+btype+" &";
  Vrun.push_back(sline);
  sline="                a "+atype;
  Vrun.push_back(sline);
  */

  sline="timestep        10.0";
  Vrun.push_back(sline);
  sline="run_style       respa 3 2 4 bond 1 angle 1 dihedral 2";
  if(npair_style>1) {
    string shybrid="hybrid";
    for(auto & it : Rcoeff.Mpair_style_args) {
      if(it.first=="neigh") continue;
      else if(it.first=="hbond/alternative")
        shybrid+=" 1";
      else {
        shybrid+=" 3";
      }
    }
    sline+=" "+shybrid;
  }
  Vrun.push_back(sline);
  //sline="fix_modify      myshake respa 1"; //for shake, the default respa level is 1
  //Vrun.push_back(sline);
  /*
  sline="unfix          mylang";
  Vrun.push_back(sline);
  if(md_label==PrepareLAMMPS::remd) {
    ostringstream oss;
    oss<<"fix             mylang other langevin $t $t 100.0 "<<83546+(int)(627043*rand_one->uniform());
    sline=oss.str();
    Vrun.push_back(sline);
  }
  else if(md_label==PrepareLAMMPS::md) {
    ostringstream oss;
    oss<<"fix             mylang other langevin "<<std::fixed<<std::setprecision(2)<<Vtemperature[0]<<" ";
    oss<<std::fixed<<std::setprecision(2)<<Vtemperature[1]<<" 100.0 "<<83546+(int)(627043*rand_one->uniform());
    sline=oss.str();
    Vrun.push_back(sline);
  }
  */
  if(md_label==PrepareLAMMPS::remd) {
    ostringstream oss;
    oss<<"temper          "<<Nstep<<" 5000 $t mylang "<<35824+(int)(165243*rand_one->uniform())<<" "<<86529+(int)(43280*rand_one->uniform());
    sline=oss.str();
    Vrun.push_back(sline);
    oss.str("");
    oss.clear();
  }
  else if(md_label==PrepareLAMMPS::md) {
    sline="run             "+to_string(Nstep);
    Vrun.push_back(sline);
  }
}

void PrepareLAMMPS::prepare_in_script_allatom(string outDIR, string jobname, string nmodel, int Nstep, int Nevery){
  Vhead.clear();
  Vio.clear();
  Vgroup.clear();
  Vem.clear();
  Vfix.clear();
  Vrun.clear();

  string sline="#LAMMPS Script for all atom energy minimization";
  Vhead.push_back(sline);
  sline="units           real";
  Vhead.push_back(sline);
  sline="dimension       3";
  Vhead.push_back(sline);
  sline="boundary        p p p";
  Vhead.push_back(sline);
  sline="atom_style      full";
  Vhead.push_back(sline);
  sline="special_bonds   amber";
  Vhead.push_back(sline);
  sline="";
  Vhead.push_back(sline);
  sline="neighbor        1.0 bin";
  Vhead.push_back(sline);
  sline="neigh_modify    every 2 check yes";
  Vhead.push_back(sline);
  sline="";
  Vhead.push_back(sline);
  if(Vbondstyle.size()==1) {
    sline="bond_style      "+Vbondstyle[0];
    Vhead.push_back(sline);
  } 
  else if(Vbondstyle.size()>1) {
    sline="bond_style      hybrid";
    for(size_t kk=0;kk<Vbondstyle.size();kk++)
      sline+=" "+Vbondstyle[kk];
    Vhead.push_back(sline);
  }
  if(Vanglestyle.size()==1) {
    sline="angle_style     "+Vanglestyle[0];
    Vhead.push_back(sline);
  }
  else if(Vanglestyle.size()>1) {
    sline="angle_style     hybrid";
    for(size_t kk=0;kk<Vanglestyle.size();kk++) {
      sline+=" "+Vanglestyle[kk];
      if(Vanglestyle[kk]=="zero")
        sline+=" nocoeff";
    }
    Vhead.push_back(sline);
  }
  if(Vdihedralstyle.size()==1) {
    sline="dihedral_style  "+Vdihedralstyle[0];
    Vhead.push_back(sline);
  }
  else if(Vdihedralstyle.size()>1) {
    sline="dihedral_style  hybrid";
    for(size_t kk=0;kk<Vdihedralstyle.size();kk++) {
      sline+=" "+Vdihedralstyle[kk];
      if(Vdihedralstyle[kk]=="zero")
        sline+=" nocoeff";
    }
    Vhead.push_back(sline);
  }
  if(Vimproperstyle.size()==1) {
    sline="improper_style  "+Vimproperstyle[0];
    Vhead.push_back(sline);
  }
  else if(Vimproperstyle.size()>1) {
    sline="improper_style  hybrid";
    for(size_t kk=0;kk<Vimproperstyle.size();kk++) {
      sline+=" "+Vimproperstyle[kk];
      if(Vimproperstyle[kk]=="zero")
        sline+=" nocoeff";
    }
    Vhead.push_back(sline);
  }
  size_t npair_style=Rcoeff.Mpair_style_args.size();
  if(Rcoeff.Mpair_style_args.count("neigh"))
    npair_style=npair_style-1;
  if(npair_style==1) {
    sline="pair_style      ";
  }
  else {
    sline="pair_style      hybrid/overlay";
  }
  for(auto & it : Rcoeff.Mpair_style_args) {
    if(it.first=="neigh") continue;
    sline+=" "+it.first+" "+it.second;
  }
  Vhead.push_back(sline);
  sline="";
  Vhead.push_back(sline);

  sline="variable        t world";
  if(md_label==PrepareLAMMPS::remd) {
    ostringstream oss;
    for(size_t kk=0;kk<Vtemperature.size();kk++)
      oss<<" "<<std::fixed<<std::setprecision(0)<<Vtemperature[kk];
    sline+=oss.str();
    oss.str("");
    oss.clear();
  }
  else {
    sline+=" "+nmodel;
  }
  Vio.push_back(sline);
  sline="read_data       "+outDIR+"data."+jobname+".";
  if(md_label==PrepareLAMMPS::remd) {
    sline+=nmodel;
  }
  else
    sline+="$t";
  Vio.push_back(sline);
  sline="";
  Vio.push_back(sline);
  if(Rcoeff.Vpaircoeff.size()>Rcoeff.Vmass.size()) {
    size_t index_begin=2;
    if(npair_style==1)
      index_begin=3;
    for(size_t jj=0;jj<Rcoeff.Vpaircoeff.size();jj++) {
      vector<string> stoken=split_string(Rcoeff.Vpaircoeff[jj]);
      if(stoken[2].find(Rcoeff.default_pair_style)!=string::npos) continue;
      string sout=stoken[0]+" "+stoken[1];
      for(size_t kk=index_begin;kk<stoken.size();kk++)
        sout+=" "+stoken[kk];
      sline="pair_coeff      "+sout;
      Vio.push_back(sline);
    }
  }
  sline="";
  Vio.push_back(sline);
  sline="kspace_style    pppm 1e-4";
  Vio.push_back(sline);
  sline="dielectric      78.285";
  Vio.push_back(sline);
  sline="";
  Vio.push_back(sline);
  sline="thermo          "+to_string(Nevery);
  Vio.push_back(sline);
  sline="log             "+outDIR+"log."+jobname+".$t";
  vector<string> stoken=split_string(nmodel);
  for(size_t kk=0;kk<stoken.size();kk++) {
    string ostr=outDIR+"log."+jobname+".";
    ostr+=stoken[kk];
    Voutput_logs.push_back(ostr);
  }
  Vio.push_back(sline);
  sline="dump            1 all custom "+to_string(Nevery)+" "+outDIR+"dump."+jobname+".";
  stoken.clear();
  stoken=split_string(nmodel);
  for(size_t kk=0;kk<stoken.size();kk++) {
    string ostr=outDIR+"dump."+jobname+".";
    ostr+=stoken[kk];
    ostr+=".lammpstrj";
    Voutput_trjs.push_back(ostr);
  }
  sline+="$t.lammpstrj id mol type xu yu zu";
  Vio.push_back(sline);
  sline="";
  Vio.push_back(sline);

  sline="group           stat id 1";
  Vgroup.push_back(sline);
  sline="group           other subtract all stat";
  Vgroup.push_back(sline);

  sline="";
  Vem.push_back(sline);
  sline="min_style       cg";
  Vem.push_back(sline);
  sline="minimize        1.0e-4 1.0e-6 "+to_string(Nstep)+" 100000";
  Vem.push_back(sline);
}

void PrepareLAMMPS::write_in_script(string filename, int reset_random) { //reset_random=1: reset random number in velocity and temper
  if(reset_random) {
    for(size_t ii=0;ii<Vrun.size();ii++) {
      string shead=Vrun[ii].substr(0,10);
      if(shead.find("velocity")!=string::npos) {
        int rand_num1=37624+(int)(913240*rand_one->uniform());
        vector<string> stoken=split_string(Vrun[ii]);
        if(stoken[2].find("create")!=string::npos) {
          string sline="velocity        "+stoken[1]+" "+stoken[2]+" "+stoken[3]+" "+to_string(rand_num1);
          for(size_t jj=5;jj<stoken.size();jj++)
            sline+=" "+stoken[jj];
          Vrun[ii]=sline;
        }
      }
      else if(shead.find("temper")!=string::npos) {
        int rand_num1=35824+(int)(165243*rand_one->uniform());
        int rand_num2=86529+(int)(43280*rand_one->uniform());
        vector<string> stoken=split_string(Vrun[ii]);
        string sline="temper          "+stoken[1]+" "+stoken[2]+" "+stoken[3]+" "+stoken[4]+" "+to_string(rand_num1)+" "+to_string(rand_num2);
        Vrun[ii]=sline;
      }
    }
  }
  
  ofstream outF (filename.c_str());
  for(auto sline : Vhead)
    outF<<sline<<endl;
  for(auto sline : Vio)
    outF<<sline<<endl;
  for(auto sline : Vgroup)
    outF<<sline<<endl;
  for(auto sline : Vem)
    outF<<sline<<endl;
  if(md_label!=PrepareLAMMPS::em) {
    for(auto sline : Vfix)
      outF<<sline<<endl;
    outF<<endl;
    for(auto sline : Vrun)
      outF<<sline<<endl;
  }
  outF.close(); outF.clear();
}

vector<string> PrepareLAMMPS::write_in_script(string outDIR, string jobname, size_t nmodel, int Nstep, int Nevery) {
  vector<string> Vrun_files;
  int nthread=init->nthread;
  int nloop=nmodel/nthread;
  int nlast=nmodel%nthread;
  for(int jj=0;jj<nloop;jj++) {
    string str=to_string(jj*nthread+1);
    for(int kk=1;kk<nthread;kk++)
      str+=" "+to_string(jj*nthread+1+kk);
    PrepareLAMMPS::prepare_in_script(outDIR,jobname,str,Nstep,Nevery);
    string fstr=outDIR+"in."+jobname+"."+to_string(jj+1);
    PrepareLAMMPS::write_in_script(fstr);
    string srun=fstr+" "+to_string(nthread)+" 1";
    Vrun_files.push_back(srun);
  }
  if(nlast>0) {
    string str=to_string(nloop*nthread+1);
    for(int kk=1;kk<nlast;kk++)
      str+=" "+to_string(nloop*nthread+1+kk);
    PrepareLAMMPS::prepare_in_script(outDIR,jobname,str,Nstep,Nevery);
    string fstr=outDIR+"in."+jobname+"."+to_string(nloop+1);
    PrepareLAMMPS::write_in_script(fstr);
    string srun=fstr+" "+to_string(nlast)+" 1";
    Vrun_files.push_back(srun);
  }

  return Vrun_files;
}

vector<string> PrepareLAMMPS::write_run_script(const vector<string> &Vruns, string filename, int gpu_label) {
  vector<string> Vcmd;
  string lmp=init->LMP;
  string lmp_gpu=lmp+"_gpu";
  if(!is_file_exist(lmp) && !is_file_exist(lmp_gpu)) {
    string str="prepare_lammps.cpp: Cannot find executable lmp! Please specify it by -lmp!";
    error->one(str);
  }
  if(!is_file_exist(lmp) && is_file_exist(lmp_gpu)) {
    lmp=lmp_gpu;
  }
  if(gpu_label && is_file_exist(lmp_gpu)) {
    lmp=lmp_gpu;
  }
  ofstream outF (filename.c_str());
  outF<<"#!/bin/bash"<<endl;
  for(size_t ii=0;ii<Vruns.size();ii++) {
    string str=" Running "+to_string(ii+1)+"/"+to_string(Vruns.size())+" job";
    outF<<"echo "<<"\""<<str<<"\""<<endl;
    vector<string> stoken=split_string(Vruns[ii]);
    if(stoken.size()!=3) {
      string str="prepare_lammps.cpp: Wrong line in Vruns --> "+Vruns[ii];
      error->one(str);
    }
    string infile=stoken[0];
    int nn=stoi(stoken[1]);
    int nm=stoi(stoken[2]);
    string cmd="mpirun ";
    if(init->MPIRUN_OPT.size()>0)
      cmd+=init->MPIRUN_OPT+" ";
    cmd+="-np "+to_string(nn*nm)+" "+lmp+" -p "+to_string(nn)+"x"+to_string(nm)+" -in "+infile;
    if(gpu_label) 
      cmd+=" -sf gpu -pk gpu 0 neigh no";
    outF<<cmd<<endl;
    Vcmd.push_back(cmd);
  }
  outF.close(); outF.clear();

  return Vcmd;
}

void PrepareLAMMPS::clear() {
  md_label=-1;
  Vtemperature.clear();
  Vtemperature.shrink_to_fit();
  Rcoeff.clear();
  Mrestraint_type.clear();
  Vhead.clear();
  Vhead.shrink_to_fit();
  Vio.clear();
  Vio.shrink_to_fit();
  Vgroup.clear();
  Vgroup.shrink_to_fit();
  Vem.clear();
  Vem.shrink_to_fit();
  Vfix.clear();
  Vfix.shrink_to_fit();
  Vrun.clear();
  Vrun.shrink_to_fit();
  Voutput_logs.clear();
  Voutput_logs.shrink_to_fit();
  Voutput_trjs.clear();
  Voutput_trjs.shrink_to_fit();
  Vbondstyle.clear();
  Vbondstyle.shrink_to_fit();
  Vanglestyle.clear();
  Vanglestyle.shrink_to_fit();
  Vdihedralstyle.clear();
  Vdihedralstyle.shrink_to_fit();
  Vimproperstyle.clear();
  Vimproperstyle.shrink_to_fit();
  Vatomtype.clear();
  Vatomtype.shrink_to_fit();
  Vbondtype.clear();
  Vbondtype.shrink_to_fit();
  Vangletype.clear();
  Vangletype.shrink_to_fit();
  Vdihedraltype.clear();
  Vdihedraltype.shrink_to_fit();
  Vimpropertype.clear();
  Vimpropertype.shrink_to_fit();
}