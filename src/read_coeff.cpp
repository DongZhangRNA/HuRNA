#include "basic.h"
#include "constant.h"
#include "readfile.h"
#include "initialization.h"
#include "error.h"
#include "read_coeff.h"

using namespace CBS_NS;

ReadCoeff::ReadCoeff(Cobmms *cbs) : Pointers(cbs) {
  Vmass.clear();
  Vbondcoeff.clear();
  Vanglecoeff.clear();
  Vdihedralcoeff.clear();
  Vimpropercoeff.clear();
  Vpaircoeff.clear();
  Mtype_charge.clear();
  Mpair_style_args.clear();

  default_pair_style="morse/gauss2/lj";
}

ReadCoeff::ReadCoeff(const ReadCoeff &other) : Pointers(other.cbs) {
  default_pair_style=other.default_pair_style;
  Vmass=other.Vmass;
  Vbondcoeff=other.Vbondcoeff;
  Vanglecoeff=other.Vanglecoeff;
  Vdihedralcoeff=other.Vdihedralcoeff;
  Vimpropercoeff=other.Vimpropercoeff;
  Vpaircoeff=other.Vpaircoeff;
  Mtype_charge=other.Mtype_charge;
  Mpair_style_args=other.Mpair_style_args;
}

ReadCoeff& ReadCoeff::operator=(const ReadCoeff &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    default_pair_style=other.default_pair_style;
    Vmass=other.Vmass;
    Vbondcoeff=other.Vbondcoeff;
    Vanglecoeff=other.Vanglecoeff;
    Vdihedralcoeff=other.Vdihedralcoeff;
    Vimpropercoeff=other.Vimpropercoeff;
    Vpaircoeff=other.Vpaircoeff;
    Mtype_charge=other.Mtype_charge;
    Mpair_style_args=other.Mpair_style_args;
  }

  return *this;
}

ReadCoeff::~ReadCoeff() {
  default_pair_style="";
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
  Mtype_charge.clear();
  Mpair_style_args.clear();
}

void ReadCoeff::read_RNA_coeff() {
  string fstr;
  if(init->RNA_version=="v2") {
    fstr=init->INPUT_DIR+"RNA_coeff_v2.dat";
    default_pair_style="morse/gauss2/lj9";
  }
  else if(init->RNA_version=="latest") {
    fstr=init->INPUT_DIR+"RNA_coeff.dat";
    default_pair_style="morse/gauss2/lj";
  }
  
  if(!is_file_exist(fstr)) {
    string str="read_coeff.cpp: Cannot find file "+fstr;
    error->one(str);
  }

  vector<string> Vlines=read_file(fstr);
  map<string,vector<string> > Mcoeffs;
  string name="";
  vector<string> Vcoeff;
  for(size_t ii=0;ii<Vlines.size();ii++) {
    vector<string> stoken=split_string(Vlines[ii]);
    if(stoken.size()<1) continue;
    string sline=stoken[0];
    for(size_t kk=1;kk<stoken.size();kk++)
      sline+=" "+stoken[kk];
    if(sline.find("#")!=string::npos) {
      size_t found=sline.find_first_of("#");
      sline=sline.substr(0,found);
    }
    if(sline.size()<2) continue;
    if(sline.find("[")!=string::npos && sline.find("]")!=string::npos) {
      if(name!="" && Vcoeff.size()>0)
        Mcoeffs.emplace(name,Vcoeff);
      vector<string> stoken=split_string(sline,"[]");
      name=stoken[0];
      Vcoeff.clear();
      continue;
    }
    stoken.clear();
    stoken=split_string(sline);
    string str="";
    for(size_t kk=0;kk<stoken.size();kk++) {
      if(stoken[kk][0]=='/' || stoken[kk][0]=='#') break;
      if(kk==0) str+=stoken[0];
      else str+=" "+stoken[kk];
    }
    Vcoeff.push_back(str);
  }
  Mcoeffs.emplace(name,Vcoeff);

  for(auto& it : Mcoeffs) {
    string name=it.first;
    for(size_t kk=0;kk<name.size();kk++)
      name[kk]=tolower(name[kk]);
    if(name=="mass") {
      for(auto sline : it.second) {
        if(sline.find("/")!=string::npos) {
          size_t found=sline.find_first_of("/");
          sline=sline.substr(0,found);
        }
        vector<string> stoken=split_string(sline);
        sline=stoken[0]+" "+stoken[1];
        if(count(Vmass.begin(),Vmass.end(),sline)==0)
          Vmass.push_back(sline);
        if(stoken.size()>2) {
          if(Mtype_charge.count(stoken[0])==0)
            Mtype_charge.emplace(stoken[0],stoken[2]);
        }
      }
    }
    else if(name=="bondcoeff") {
      for(auto sline : it.second) {
        if(count(Vbondcoeff.begin(),Vbondcoeff.end(),sline)==0)
          Vbondcoeff.push_back(sline);
      }
    }
    else if(name=="anglecoeff") {
      for(auto sline : it.second) {
        if(count(Vanglecoeff.begin(),Vanglecoeff.end(),sline)==0)
          Vanglecoeff.push_back(sline);
      }
    }
    else if(name=="dihedralcoeff") {
      for(auto sline : it.second) {
        if(count(Vdihedralcoeff.begin(),Vdihedralcoeff.end(),sline)==0)
          Vdihedralcoeff.push_back(sline);
      }
    }
    else if(name=="impropercoeff") {
      for(auto sline : it.second) {
        if(count(Vimpropercoeff.begin(),Vimpropercoeff.end(),sline)==0)
          Vimpropercoeff.push_back(sline);
      }
    }
    else if(name=="paircoeff") {
      for(auto sline : it.second) {
        if(count(Vpaircoeff.begin(),Vpaircoeff.end(),sline)==0)
          Vpaircoeff.push_back(sline);
      }
    }
  }
}

void ReadCoeff::add_coeff(string sline) {
  vector<string> stoken=split_string(sline);
  if(stoken.size()<3) {
    string str="read_coeff.cpp: Wrong coeff line--> "+sline;
    error->one(str);
  }
  string stype=stoken[0];
  string scoeff=stoken[1];
  for(size_t jj=2;jj<stoken.size();jj++)
    scoeff+=" "+stoken[jj];
  if(stype=="bond_coeff")
    Vbondcoeff.push_back(scoeff);
  else if(stype=="angle_coeff")
    Vanglecoeff.push_back(scoeff);
  else if(stype=="dihedral_coeff")
    Vdihedralcoeff.push_back(scoeff);
  else if(stype=="improper_coeff")
    Vimpropercoeff.push_back(scoeff);
  else if(stype=="pair_coeff")
    Vpaircoeff.push_back(scoeff);
  else {
    string str="read_coeff.cpp: Wrong coeff line (bond/angle/dihedral/improper/pair_coeff label style coeff) --> "+sline;
    error->one(str);
  }
}

void ReadCoeff::reset_paircoeff_lj() {
  vector<string> Tpaircoeff=Vpaircoeff;
  Vpaircoeff.clear();

  string pairtype_one="lj/cut";
  for(size_t ii=0;ii<Tpaircoeff.size();ii++) {
    string sline=Tpaircoeff[ii];
    if(sline.find(default_pair_style)!=string::npos) {
      vector<string> stoken=split_string(sline);
      if(stoken.size()<14) {
        string str="read_coeff.cpp: Wrong pair coeff --> "+sline;
        error->one(str);
      }
      string sigma=stoken[stoken.size()-2];
      string scoeff=stoken[0]+" "+stoken[1]+" "+pairtype_one+" 1.0 "+sigma+" "+sigma;
      Vpaircoeff.push_back(scoeff);
    }
    else if(sline.find("coul/debye")!=string::npos) {
      Vpaircoeff.push_back(sline);
    }
    else
      continue;
  }
  default_pair_style=pairtype_one;
}

void ReadCoeff::reset_paircoeff_soft(double A_one, double sigma_offset) {
  vector<string> Tpaircoeff=Vpaircoeff;
  Vpaircoeff.clear();

  string pairtype_one="soft";
  for(size_t ii=0;ii<Tpaircoeff.size();ii++) {
    string sline=Tpaircoeff[ii];
    if(sline.find(default_pair_style)!=string::npos) {
      vector<string> stoken=split_string(sline);
      if(stoken.size()<5) {
        string str="read_coeff.cpp: Wrong pair coeff --> "+sline;
        error->one(str);
      }
      if(stoken[2]==pairtype_one) {
        double sigma=stod(stoken[4]);
        sigma+=sigma_offset;
        ostringstream oss;
        oss<<std::fixed<<std::setprecision(1)<<A_one;
        oss<<" "<<std::fixed<<std::setprecision(1)<<sigma;
        string scoeff=stoken[0]+" "+stoken[1]+" "+pairtype_one+" "+oss.str();
        Vpaircoeff.push_back(scoeff);
        oss.str("");
        oss.clear();
      }
      else {
        if(stoken.size()<14) {
          string str="read_coeff.cpp: Wrong pair coeff --> "+sline;
          error->one(str);
        }
        double sigma=stod(stoken[stoken.size()-2]);
        sigma+=sigma_offset;
        ostringstream oss;
        oss<<std::fixed<<std::setprecision(1)<<A_one;
        oss<<" "<<std::fixed<<std::setprecision(1)<<sigma;
        string scoeff=stoken[0]+" "+stoken[1]+" "+pairtype_one+" "+oss.str();
        Vpaircoeff.push_back(scoeff);
        oss.str("");
        oss.clear();
      }
    }
    //else if(sline.find("coul/debye")!=string::npos) {
      //Vpaircoeff.push_back(sline);
    //}
    else
      continue;
  }
  default_pair_style=pairtype_one;
}

void ReadCoeff::get_interaction_style() {
  Mpair_style_args.clear();
  struct PCUT {
    double cut1,cut2,cut3;
  };
  map<string,PCUT> Mtmp;
  string kappa;
  int flag_pairij=1;
  if(Vpaircoeff.size()==Vmass.size()) flag_pairij=0;
  for(size_t ii=0;ii<Vpaircoeff.size();ii++) {
    string sline=Vpaircoeff[ii];
    vector<string> stoken=split_string(Vpaircoeff[ii]);
    string style_one;
    if(flag_pairij==0)
      style_one=stoken[1];
    else
      style_one=stoken[2];
    PCUT cut_one;
    if(style_one=="hbond/alternative") {
      if(stoken.size()!=13) {
        string str="read_coeff.cpp: Wrong pair coeff --> "+Vpaircoeff[ii];
        error->one(str);
      }
      size_t len=stoken.size();
      cut_one.cut1=stod(stoken[len-3]);
      cut_one.cut2=stod(stoken[len-2]);
      cut_one.cut3=stod(stoken[len-1]);
    }
    else {
      cut_one.cut1=stod(stoken.back());
      cut_one.cut2=0.0;
      cut_one.cut3=0.0;
    }
    if(Mtmp.count(style_one)==0) {
      Mtmp.emplace(style_one,cut_one);
    }
    else {
      PCUT cut_tmp=Mtmp[style_one];
      if(cut_tmp.cut1<cut_one.cut1)
        cut_tmp.cut1=cut_one.cut1;
      if(cut_tmp.cut2<cut_one.cut2)
        cut_tmp.cut2=cut_one.cut2;
      if(cut_tmp.cut3<cut_one.cut3)
        cut_tmp.cut3=cut_one.cut3;
      Mtmp[style_one]=cut_tmp;
    }
    if(style_one.find("coul/debye")!=string::npos) {
      if((flag_pairij==0 && stoken.size()!=4) || (flag_pairij==1 && stoken.size()!=5)) {
        string str="read_coeff.cpp: Wrong pair coeff --> "+Vpaircoeff[ii];
        error->one(str);
      }
      if(flag_pairij==0) {
        kappa=stoken[2];
        sline=stoken[0]+" "+stoken[1]+" "+stoken[3];
      }
      else {
        kappa=stoken[3];
        sline=stoken[0]+" "+stoken[1]+" "+stoken[2]+" "+stoken[4];
      }
      Vpaircoeff[ii]=sline;
    }
  }
  double cut_max=0.0;
  for(auto& it :Mtmp) {
    string style_one=it.first;
    ostringstream oss;
    if(style_one=="hbond/alternative") {
      oss<<std::fixed<<std::setprecision(1)<<it.second.cut1<<" ";
      oss<<std::fixed<<std::setprecision(1)<<it.second.cut2<<" ";
      oss<<std::fixed<<std::setprecision(1)<<it.second.cut3;
    }
    else {
      oss<<std::fixed<<std::setprecision(1)<<it.second.cut1;
      if(cut_max<it.second.cut1)
        cut_max=it.second.cut1;
    }
    string sargs=oss.str();
    oss.str("");
    oss.clear();
    if(style_one.find("coul/debye")!=string::npos)
      sargs=kappa+" "+sargs;
    if(Mpair_style_args.count(style_one)==0) {
      Mpair_style_args.emplace(style_one,sargs);
    }
    else {
      Mpair_style_args[style_one]=sargs;
    }

    string str="Pair style args: "+style_one+" "+sargs;
    error->message(str);
  }

  if(Mpair_style_args.size()==0) {
    string str="read_coeff.cpp: No pair coeff!";
    error->one(str);
  }
  else if(Mpair_style_args.size()==1) {
    auto it=Mpair_style_args.begin();
    default_pair_style=it->first;
  }

  if(Mpair_style_args.size()>0) {
    if(Mpair_style_args.count(default_pair_style)==0) {
      string str="read_coeff.cpp: Cannot find default_pair_style \""+default_pair_style+"\" in paircoeff!";
      error->one(str);
    }
  }

  cut_max+=2.0;
  ostringstream oss;
  oss<<std::fixed<<std::setprecision(1)<<cut_max;
  string sarg=oss.str();
  oss.str("");
  oss.clear();
  Mpair_style_args.emplace("neigh",sarg);
}

void ReadCoeff::add_extra_RNA_basepair_coeff(string bp_app) {
  string BONDTYPE[]={"GO-CN","GN-YO","AN-UO","AC-YO","GO-UO","GO-YO","S-S","P-P"};
  string BONDTYPE_EXT[]={"RC-YC","S2-S2"};
  if(Vbondcoeff.size()==0) {
    string str="read_coeff.cpp: Please read RNA coeffs first!";
    error->one(str);
  }

  map<string,string> Mbondtype_coeff;
  for(size_t ii=0;ii<Vbondcoeff.size();ii++) {
    vector<string> stoken=split_string(Vbondcoeff[ii]);
    string type_one=stoken[0];
    string coeff_one=stoken[1];
    for(size_t kk=2;kk<stoken.size();kk++)
      coeff_one+=" "+stoken[kk];
    if(Mbondtype_coeff.count(type_one)==0)
      Mbondtype_coeff.emplace(type_one,coeff_one);
  }
  
  vector<string> Vbpcoeff;
  for(int ii=0;ii<get_array_size(BONDTYPE);ii++) {
    string type_one=BONDTYPE[ii];
    string type_one_app=type_one+"_"+bp_app;
    if(Mbondtype_coeff.count(type_one)==0) {
      string str="read_coeff.cpp: Cannot find bond coeff for base pairing interaction "+type_one;
      error->one(str);
    }
    if(Mbondtype_coeff.count(type_one_app)==0) {
      string coeff_one=type_one_app+" "+Mbondtype_coeff[type_one];
      Vbpcoeff.push_back(coeff_one);
    }
  }
  
  if(init->RNA_version=="latest") {
    for(int ii=0;ii<get_array_size(BONDTYPE_EXT);ii++) {
      string type_one=BONDTYPE_EXT[ii];
      string type_one_app=type_one+"_"+bp_app;
      if(Mbondtype_coeff.count(type_one)==0) {
        string str="read_coeff.cpp: Cannot find bond coeff for base pairing interaction "+type_one;
        error->one(str);
      }
      if(Mbondtype_coeff.count(type_one_app)==0) {
        string coeff_one=type_one_app+" "+Mbondtype_coeff[type_one];
        Vbpcoeff.push_back(coeff_one);
      }
    }
  }
  
  for(size_t ii=0;ii<Vbpcoeff.size();ii++)
    Vbondcoeff.push_back(Vbpcoeff[ii]);
}

vector<string> ReadCoeff::get_weighted_RNA_basepair_coeff(string bp_app, double wbp) const{
  string BONDTYPE[]={"GO-CN","GN-YO","AN-UO","AC-YO","GO-UO","GO-YO","RC-YC","S2-S2","S-S","P-P"};
  if(Vbondcoeff.size()==0) {
    string str="read_coeff.cpp: Please read RNA coeffs first!";
    error->one(str);
  }

  vector<string> Vbpcoeff;
  for(size_t ii=0;ii<Vbondcoeff.size();ii++) {
    vector<string> stoken=split_string(Vbondcoeff[ii]);
    string type_one=stoken[0];
    if(count(begin(BONDTYPE),end(BONDTYPE),type_one)) {
      if(stoken.size()<7) {
        string str="read_coeff.cpp: Wrong bond coeffs for RNA basepair type "+type_one;
        error->one(str);
      }
      type_one+="_"+bp_app;
      string coeff_one=type_one+" "+stoken[1];
      double k=wbp*stod(stoken[2]);
      if(k<0.001) k=0.001;
      string r0=stoken[3];
      double H=wbp*stod(stoken[4]);
      string r1=stoken[5];
      string sigma=stoken[6];
      ostringstream oss;
      oss<<std::fixed<<std::setprecision(3)<<k<<" "<<r0<<" ";
      oss<<std::fixed<<std::setprecision(3)<<H<<" "<<r1<<" "<<sigma;
      coeff_one+=" "+oss.str();
      Vbpcoeff.push_back(coeff_one);
      oss.str("");
      oss.clear();
    }
  }

  if(Vbpcoeff.size()==0) {
    string str="read_coeff.cpp: Please read RNA coeffs first!";
    error->one(str);
  }

  return Vbpcoeff;
}

vector<string> ReadCoeff::get_pair_coeff() const {
  vector<string> Vpairstyle;
  for(size_t jj=0;jj<Vpaircoeff.size();jj++) {
    vector<string> stoken=split_string(Vpaircoeff[jj]);
    if(stoken.size()<3) continue;
    string style_one=stoken[2];
    if(count(Vpairstyle.begin(),Vpairstyle.end(),style_one)==0)
      Vpairstyle.push_back(style_one);
  }
  vector<string> Vcoeff;
  size_t index_begin=2;
  if(Vpairstyle.size()==1) index_begin=3;
  for(size_t jj=0;jj<Vpaircoeff.size();jj++) {
    vector<string> stoken=split_string(Vpaircoeff[jj]);
    string sline="pair_coeff      "+stoken[0]+" "+stoken[1];
    if(stoken[2].find(default_pair_style)!=string::npos) {
      for(size_t kk=index_begin;kk<stoken.size();kk++)
        sline+=" "+stoken[kk];
      Vcoeff.push_back(sline);
    }
  }
  
  return Vcoeff;
}

void ReadCoeff::clear() {
  default_pair_style="";
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
  Mtype_charge.clear();
  Mpair_style_args.clear();
}