#include "basic.h"
#include "constant.h"
#include "atom.h"

using namespace CBS_NS;

Atom::Atom() {
  id=-1;
  resid=-1;
  chainid="";
  resname="";
  atomname="";
  atomalt="";
  iCode=' ';
  q=0.0;
  x=0.0;
  y=0.0;
  z=0.0;
}

Atom::Atom(string sline) {
  Atom::set_atom(sline);
}

Atom::Atom(const Atom &other) {
  Atom::clear();
  id=other.id;
  chainid=other.chainid;
  resid=other.resid;
  resname=other.resname;
  atomname=other.atomname;
  q=other.q;
  x=other.x;
  y=other.y;
  z=other.z;
  atomalt=other.atomalt;
  iCode=other.iCode;
}

Atom::~Atom() {
  id=-1;
  resid=-1;
  chainid="";
  resname="";
  atomname="";
  atomalt="";
  iCode=' ';
  q=0.0;
  x=0.0;
  y=0.0;
  z=0.0;
}

void Atom::set_atom(int id_one, string cid_one, int resid_one, string resname_one, string atomname_one, double x_one, double y_one, double z_one) {
  id=id_one;
  chainid=cid_one;
  resid=resid_one;
  resname=resname_one;
  atomname=atomname_one;
  x=x_one;
  y=y_one;
  z=z_one;
  atomalt="";
  iCode=' ';
  q=0.0;
}

void Atom::set_atom(string atomname_one, double x_one, double y_one, double z_one) {
  id=0;
  chainid="";
  resid=0;
  resname="";
  atomname=atomname_one;
  x=x_one;
  y=y_one;
  z=z_one;
  atomalt="";
  iCode=' ';
  q=0.0;
}

void Atom::set_atom(size_t id_one, string cid_one, int resid_one, string resname_one, string atomname_one, vector<double> &Vxyz) {
  id=static_cast<int>(id_one);
  chainid=cid_one;
  resid=resid_one;
  resname=resname_one;
  atomname=atomname_one;
  x=Vxyz[0];
  y=Vxyz[1];
  z=Vxyz[2];
  atomalt="";
  iCode=' ';
  q=0.0;
}

void Atom::set_atom(string sline) {  //sline in pdb format
  if(sline.size()<54) {
    string str="atom.cpp: Not standard PDB format --> "+sline;
    cout<<"Error: "<<str<<endl;
    exit(0);
  }

  string serial,aname,altLoc,rName,chainID,resSeq,sx,sy,sz,head,occupy,tempfactor;

  string str0=sline.substr(0,6);
  for(size_t i=0;i<str0.size();i++)
    if(!isspace(str0[i])) head+=str0[i];

  string str1=sline.substr(6,5);
  for(size_t i=0;i<str1.size();i++)
    if(!isspace(str1[i])) {serial+=str1[i];}

  string str2=sline.substr(12,4);
  for(size_t i=0;i<str2.size();i++)
    if(!isspace(str2[i])) {aname+=str2[i];}

  altLoc=sline[16];

  string str3=sline.substr(17,3);
  for(size_t i=0;i<str3.size();i++)
    if(!isspace(str3[i])) {rName+=str3[i];}

  chainID=sline[21];

  string str4=sline.substr(22,4);
  for(size_t i=0;i<str4.size();i++)
    if(!isspace(str4[i])) {resSeq+=str4[i];}

  iCode=sline[26];
  //resSeq+=iCode;

  string str5=sline.substr(30,8);
  for(size_t i=0;i<str5.size();i++)
    if(!isspace(str5[i])) {sx+=str5[i];}

  string str6=sline.substr(38,8);
  for(size_t i=0;i<str6.size();i++)
    if(!isspace(str6[i])) {sy+=str6[i];}

  string str7=sline.substr(46,8);
  for(size_t i=0;i<str7.size();i++)
    if(!isspace(str7[i])) {sz+=str7[i];}

  if(sline.size()>=54) {
    string str8=sline.substr(54,6);
    for(size_t i=0;i<str8.size();i++)
      if(!isspace(str8[i])) {occupy+=str8[i];}
  }

  if(sline.size()>=60) {
    string str9=sline.substr(60,6);
    for(size_t i=0;i<str9.size();i++)
      if(!isspace(str9[i])) {tempfactor+=str9[i];}
  }

  if(aname[0]>='0' && aname[0]<='9' && aname.length()>1 && aname[1]=='H') {
    string tmp=aname;
    aname=tmp.substr(1);
    aname+=to_string(tmp[0]);
  }

  id=stoi(serial);
  atomname=aname;
  resname=rName;
  chainid=chainID;
  resid=stoi(resSeq);
  x=stod(sx);
  y=stod(sy);
  z=stod(sz);
  atomalt=altLoc;
  q=0.0;
}

void Atom::set_atom(const Atom &atom_one) {
  id=atom_one.id;
  chainid=atom_one.chainid;
  resid=atom_one.resid;
  resname=atom_one.resname;
  atomname=atom_one.atomname;
  x=atom_one.x;
  y=atom_one.y;
  z=atom_one.z;
  atomalt=atom_one.atomalt;
  iCode=atom_one.iCode;
  q=atom_one.q;
}

string Atom::get_residue_name() const{
  string name=chainid+"/"+to_string(resid);
  if(!isspace(iCode)) name+=to_string((char)iCode);
  name+="/"+resname;
  if(chainid.size()==0 && resname.size()==0) name=ResidueConst::RESIDUE_NULL;

  return name;
}

string Atom::get_atom_name() const {
  string name=chainid+"/"+to_string(resid);
  if(!isspace(iCode)) name+=to_string((char)iCode);
  name+="/"+resname;
  if(chainid.size()==0 && resname.size()==0) name=ResidueConst::RESIDUE_NULL;
  name+="/"+atomname;

  return name;
}

void Atom::reset_residue_name(string name_one) {
  vector<string> stoken=split_string(name_one,"/");
  if(stoken.size()!=3) {
    string str="atom.cpp: Wrong residue name(chainid/resid/resname) --> "+name_one;
    cout<<"Error: "<<str<<endl;
    exit(0);
  }

  chainid=stoken[0];
  resid=stoi(stoken[1]);
  resname=stoken[2];
}

void Atom::reset_residue_name(string chainid_one, int resid_one, string resname_one) {
  chainid=chainid_one;
  resid=resid_one;
  resname=resname_one;
}

void Atom::reset_atom_xyz(const vector<double> &Vxyz) {
  x=Vxyz[0];
  y=Vxyz[1];
  z=Vxyz[2];
}

void Atom::reset_atom_xyz(const Atom &atom_one) {
  x=atom_one.x;
  y=atom_one.y;
  z=atom_one.z;
}

void Atom::set_atom_xyz(const Atom &atom_one) {
  x=atom_one.x;
  y=atom_one.y;
  z=atom_one.z;
}

void Atom::move_atom(double dx, double dy, double dz) {
  x+=dx;
  y+=dy;
  z+=dz;
}

void Atom::move_atom(const vector<double> & Vxyz) {
  if(Vxyz.size()!=3) {
    string str="atom.cpp: Vxyz.size()!=3 --> "+to_string(Vxyz.size());
    cout<<"Error: "<<str<<endl;
    exit(0);
  }

  x+=Vxyz[0];
  y+=Vxyz[1];
  z+=Vxyz[2];
}

vector<double> Atom::get_atom_xyz() const {
  vector<double> Vxyz;
  Vxyz.push_back(x);
  Vxyz.push_back(y);
  Vxyz.push_back(z);

  return Vxyz;
}

double Atom::distance(const Atom &atom_one) {
  double dis=(x-atom_one.x)*(x-atom_one.x)+(y-atom_one.y)*(y-atom_one.y)+(z-atom_one.z)*(z-atom_one.z);
  dis=sqrt(dis);

  return dis;
}

void Atom::clear() {
  id=-1;
  resid=-1;
  chainid="";
  resname="";
  atomname="";
  atomalt="";
  iCode=' ';
  q=0.0;
  x=0.0;
  y=0.0;
  z=0.0;
}

ostream& operator <<(ostream &output, Atom &atom_one) {
  char aname[5]="    ";
  size_t len=atom_one.atomname.size();
  if(len==4) {
    for(size_t ll=0;ll<len;ll++)
      aname[ll]=atom_one.atomname[ll];
  }
  else {
    for(size_t ll=0;ll<len;ll++)
      aname[ll+1]=atom_one.atomname[ll];
  }
  char resseq[4]="   ";
  string resname=atom_one.resname;
  if(resname.size()>3) 
    resname=resname.substr(resname.size()-3);
  len=resname.size();
  for(int ii=len-1;ii>=0;ii--)
    resseq[3-(len-ii)]=resname[ii];
  char cid;
  if(atom_one.chainid.size()==0) 
    cid=' ';
  else
    cid=atom_one.chainid[0]; 
  char str[128];
  sprintf(str,"ATOM  %5d %s %s %c%4d%c   %8.3f%8.3f%8.3f  1.00  0.00           %c",
              atom_one.id,aname,resseq,cid,atom_one.resid,atom_one.iCode,atom_one.x,atom_one.y,atom_one.z,atom_one.atomname[0]);
  output<<str<<endl;

  return output;
}