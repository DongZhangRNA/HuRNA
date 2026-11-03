#include "basic.h"
#include "constant.h"
#include "measure.h"

using namespace CBS_NS;

enum{EQUAL,CONTAINED,INTERSECTED,SEPARATED};
const int Measure::equal=EQUAL;
const int Measure::contained=CONTAINED;
const int Measure::intersected=INTERSECTED;
const int Measure::separated=SEPARATED;

Measure::Measure() {}

Measure::~Measure() {}

double Measure::distance(double x1[3], double x2[3]) {
  double dx=x1[0]-x2[0];
  double dy=x1[1]-x2[1];
  double dz=x1[2]-x2[2];
  double rsq=dx*dx+dy*dy+dz*dz;
  double r=sqrt(rsq);

  return r;
}

double Measure::distance(const vector<double> &x1, const vector<double> &x2) {
  double dx=x1[0]-x2[0];
  double dy=x1[1]-x2[1];
  double dz=x1[2]-x2[2];
  double rsq=dx*dx+dy*dy+dz*dz;
  double r=sqrt(rsq);

  return r;
}

double Measure::distance(const Atom &atom_i, const Atom &atom_j) {
  vector<double> x_i=atom_i.get_atom_xyz();
  vector<double> x_j=atom_j.get_atom_xyz();

  double dx=x_i[0]-x_j[0];
  double dy=x_i[1]-x_j[1];
  double dz=x_i[2]-x_j[2];
  double rsq=dx*dx+dy*dy+dz*dz;
  double r=sqrt(rsq);

  return r;
}

double Measure::angle(double x1[3], double x2[3], double x3[3]) {
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double rsq1,rsq2,r1,r2,c,theta;

  // 1st bond
  delx1=x1[0]-x2[0];
  dely1=x1[1]-x2[1];
  delz1=x1[2]-x2[2];
  rsq1=delx1*delx1+dely1*dely1+delz1*delz1;
  r1=sqrt(rsq1);
  // 2nd bond
  delx2=x3[0]-x2[0];
  dely2=x3[1]-x2[1];
  delz2=x3[2]-x2[2];
  rsq2=delx2*delx2+dely2*dely2+delz2*delz2;
  r2=sqrt(rsq2);
  // c = cosine of angle
  c=delx1*delx2+dely1*dely2+delz1*delz2;
  c/=r1*r2;
  if(c>1.0) c=1.0;
  if(c<-1.0) c=-1.0;
  theta=acos(c)*180/MathConst::PI;

  if(theta<=0) theta=0;
  else if(theta>=180) theta=180;

  return theta;
}

double Measure::angle(const vector<double> &x1, const vector<double> &x2, const vector<double> &x3) {
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double rsq1,rsq2,r1,r2,c,theta;

  // 1st bond
  delx1=x1[0]-x2[0];
  dely1=x1[1]-x2[1];
  delz1=x1[2]-x2[2];
  rsq1=delx1*delx1+dely1*dely1+delz1*delz1;
  r1=sqrt(rsq1);
  // 2nd bond
  delx2=x3[0]-x2[0];
  dely2=x3[1]-x2[1];
  delz2=x3[2]-x2[2];
  rsq2=delx2*delx2+dely2*dely2+delz2*delz2;
  r2=sqrt(rsq2);
  // c = cosine of angle
  c=delx1*delx2+dely1*dely2+delz1*delz2;
  c/=r1*r2;
  if(c>1.0) c=1.0;
  if(c<-1.0) c=-1.0;
  theta=acos(c)*180/MathConst::PI;

  if(theta<=0) theta=0;
  else if(theta>=180) theta=180;

  return theta;
}

double Measure::angle(const Atom &atom_i, const Atom &atom_j, const Atom &atom_k) {
  vector<double> x_i=atom_i.get_atom_xyz();
  vector<double> x_j=atom_j.get_atom_xyz();
  vector<double> x_k=atom_k.get_atom_xyz();

  // 1st bond
  double delx1=x_i[0]-x_j[0];
  double dely1=x_i[1]-x_j[1];
  double delz1=x_i[2]-x_j[2];
  double rsq1=delx1*delx1+dely1*dely1+delz1*delz1;
  double r1=sqrt(rsq1);
  // 2nd bond
  double delx2=x_k[0]-x_j[0];
  double dely2=x_k[1]-x_j[1];
  double delz2=x_k[2]-x_j[2];
  double rsq2=delx2*delx2+dely2*dely2+delz2*delz2;
  double r2=sqrt(rsq2);
  // c = cosine of angle
  double c=delx1*delx2+dely1*dely2+delz1*delz2;
  c/=r1*r2;
  if(c>1.0) c=1.0;
  if(c<-1.0) c=-1.0;
  double theta=acos(c)*180/MathConst::PI;

  if(theta<=0) theta=0;
  else if(theta>=180) theta=180;

  return theta;
}

double Measure::dihedral(double x1[3], double x2[3], double x3[3], double x4[3]) {
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double sb1,sb3,rb1,rb3,c0,b1mag2,b1mag,b2mag2;
  double b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2;
  double c2mag,sin2,sc1,sc2,s12,c,cx,cy,cz,cmag,dx,phi;

  // 1st bond
  vb1x=x1[0]-x2[0];
  vb1y=x1[1]-x2[1];
  vb1z=x1[2]-x2[2];
  // 2nd bond
  vb2x=x3[0]-x2[0];
  vb2y=x3[1]-x2[1];
  vb2z=x3[2]-x2[2];
  vb2xm=-vb2x;
  vb2ym=-vb2y;
  vb2zm=-vb2z;
  // 3rd bond
  vb3x=x4[0]-x3[0];
  vb3y=x4[1]-x3[1];
  vb3z=x4[2]-x3[2];
  // c0 calculation
  sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
  sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);
  rb1 = sqrt(sb1);
  rb3 = sqrt(sb3);
  c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;
  // 1st and 2nd angle
  b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
  b1mag = sqrt(b1mag2);
  b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
  b2mag = sqrt(b2mag2);
  b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
  b3mag = sqrt(b3mag2);
  ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z;
  r12c1 = 1.0 / (b1mag*b2mag);
  c1mag = ctmp * r12c1;
  ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z;
  r12c2 = 1.0 / (b2mag*b3mag);
  c2mag = ctmp * r12c2;
  // cos and sin of 2 angles and final c
  sin2 = max(1.0 - c1mag*c1mag,0.0);
  sc1 = sqrt(sin2);
  if (sc1 < MathConst::EPS) sc1 = MathConst::EPS;
  sc1 = 1.0/sc1;
  sin2 = max(1.0 - c2mag*c2mag,0.0);
  sc2 = sqrt(sin2);
  if (sc2 < MathConst::EPS) sc2 = MathConst::EPS;
  sc2 = 1.0/sc2;
  s12 = sc1 * sc2;
  c = (c0 + c1mag*c2mag) * s12;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  cx = vb1y*vb2z - vb1z*vb2y;
  cy = vb1z*vb2x - vb1x*vb2z;
  cz = vb1x*vb2y - vb1y*vb2x;
  cmag = sqrt(cx*cx + cy*cy + cz*cz);
  dx = (cx*vb3x + cy*vb3y + cz*vb3z)/cmag/b3mag;

  //phi
  phi=acos(c)*180/MathConst::PI;

  if(phi>=180) phi=180;
  else if(phi<=0) phi=0;

  if(dx>0.0) phi=-phi;

  return phi;
}

double Measure::dihedral(const vector<double> &x1, const vector<double> &x2, const vector<double> &x3, const vector<double> &x4) {
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double sb1,sb3,rb1,rb3,c0,b1mag2,b1mag,b2mag2;
  double b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2;
  double c2mag,sin2,sc1,sc2,s12,c,cx,cy,cz,cmag,dx,phi;

  // 1st bond
  vb1x=x1[0]-x2[0];
  vb1y=x1[1]-x2[1];
  vb1z=x1[2]-x2[2];
  // 2nd bond
  vb2x=x3[0]-x2[0];
  vb2y=x3[1]-x2[1];
  vb2z=x3[2]-x2[2];
  vb2xm=-vb2x;
  vb2ym=-vb2y;
  vb2zm=-vb2z;
  // 3rd bond
  vb3x=x4[0]-x3[0];
  vb3y=x4[1]-x3[1];
  vb3z=x4[2]-x3[2];
  // c0 calculation
  sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
  sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);
  rb1 = sqrt(sb1);
  rb3 = sqrt(sb3);
  c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;
  // 1st and 2nd angle
  b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
  b1mag = sqrt(b1mag2);
  b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
  b2mag = sqrt(b2mag2);
  b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
  b3mag = sqrt(b3mag2);
  ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z;
  r12c1 = 1.0 / (b1mag*b2mag);
  c1mag = ctmp * r12c1;
  ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z;
  r12c2 = 1.0 / (b2mag*b3mag);
  c2mag = ctmp * r12c2;
  // cos and sin of 2 angles and final c
  sin2 = max(1.0 - c1mag*c1mag,0.0);
  sc1 = sqrt(sin2);
  if (sc1 < MathConst::EPS) sc1 = MathConst::EPS;
  sc1 = 1.0/sc1;
  sin2 = max(1.0 - c2mag*c2mag,0.0);
  sc2 = sqrt(sin2);
  if (sc2 < MathConst::EPS) sc2 = MathConst::EPS;
  sc2 = 1.0/sc2;
  s12 = sc1 * sc2;
  c = (c0 + c1mag*c2mag) * s12;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  cx = vb1y*vb2z - vb1z*vb2y;
  cy = vb1z*vb2x - vb1x*vb2z;
  cz = vb1x*vb2y - vb1y*vb2x;
  cmag = sqrt(cx*cx + cy*cy + cz*cz);
  dx = (cx*vb3x + cy*vb3y + cz*vb3z)/cmag/b3mag;

  //phi
  phi=acos(c)*180/MathConst::PI;

  if(phi>=180) phi=180;
  else if(phi<=0) phi=0;

  if(dx>0.0) phi=-phi;

  return phi;
}

double Measure::dihedral(const Atom &atom_i, const Atom &atom_j, const Atom &atom_k, const Atom &atom_l) {
  vector<double> x_i=atom_i.get_atom_xyz();
  vector<double> x_j=atom_j.get_atom_xyz();
  vector<double> x_k=atom_k.get_atom_xyz();
  vector<double> x_l=atom_l.get_atom_xyz();

  // 1st bond
  double vb1x=x_i[0]-x_j[0];
  double vb1y=x_i[1]-x_j[1];
  double vb1z=x_i[2]-x_j[2];
  // 2nd bond
  double vb2x=x_k[0]-x_j[0];
  double vb2y=x_k[1]-x_j[1];
  double vb2z=x_k[2]-x_j[2];
  double vb2xm=-vb2x;
  double vb2ym=-vb2y;
  double vb2zm=-vb2z;
  // 3rd bond
  double vb3x=x_l[0]-x_k[0];
  double vb3y=x_l[1]-x_k[1];
  double vb3z=x_l[2]-x_k[2];
  // c0 calculation
  double sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
  double sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);
  double rb1 = sqrt(sb1);
  double rb3 = sqrt(sb3);
  double c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;
  // 1st and 2nd angle
  double b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
  double b1mag = sqrt(b1mag2);
  double b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
  double b2mag = sqrt(b2mag2);
  double b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
  double b3mag = sqrt(b3mag2);
  double ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z;
  double r12c1 = 1.0 / (b1mag*b2mag);
  double c1mag = ctmp * r12c1;
  ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z;
  double r12c2 = 1.0 / (b2mag*b3mag);
  double c2mag = ctmp * r12c2;
  // cos and sin of 2 angles and final c
  double sin2 = max(1.0 - c1mag*c1mag,0.0);
  double sc1 = sqrt(sin2);
  if (sc1 < MathConst::EPS) sc1 = MathConst::EPS;
  sc1 = 1.0/sc1;
  sin2 = max(1.0 - c2mag*c2mag,0.0);
  double sc2 = sqrt(sin2);
  if (sc2 < MathConst::EPS) sc2 = MathConst::EPS;
  sc2 = 1.0/sc2;
  double s12 = sc1 * sc2;
  double c = (c0 + c1mag*c2mag) * s12;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  double cx = vb1y*vb2z - vb1z*vb2y;
  double cy = vb1z*vb2x - vb1x*vb2z;
  double cz = vb1x*vb2y - vb1y*vb2x;
  double cmag = sqrt(cx*cx + cy*cy + cz*cz);
  double dx = (cx*vb3x + cy*vb3y + cz*vb3z)/cmag/b3mag;

  //phi
  double phi=acos(c)*180/MathConst::PI;

  if(phi>=180) phi=180;
  else if(phi<=0) phi=0;

  if(dx>0.0) phi=-phi;

  return phi;
}

double Measure::improper(double x1[3], double x2[3], double x3[3], double x4[3]) {
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2;
  double s12,c,chi;
  
  vb1x = x1[0] - x2[0];
  vb1y = x1[1] - x2[1];
  vb1z = x1[2] - x2[2];

  vb2x = x3[0] - x2[0];
  vb2y = x3[1] - x2[1];
  vb2z = x3[2] - x2[2];

  vb3x = x4[0] - x3[0];
  vb3y = x4[1] - x3[1];
  vb3z = x4[2] - x3[2];

  ss1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
  ss2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
  ss3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);

  r1 = sqrt(ss1);
  r2 = sqrt(ss2);
  r3 = sqrt(ss3);

  // sin and cos of angle

  c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3;
  c1 = (vb1x * vb2x + vb1y * vb2y + vb1z * vb2z) * r1 * r2;
  c2 = -(vb3x * vb2x + vb3y * vb2y + vb3z * vb2z) * r3 * r2;

  s1 = 1.0 - c1*c1;
  if (s1 < MathConst::EPS) s1 = MathConst::EPS;
  s1 = 1.0 / s1;

  s2 = 1.0 - c2*c2;
  if (s2 < MathConst::EPS) s2 = MathConst::EPS;
  s2 = 1.0 / s2;

  s12 = sqrt(s1*s2);
  c = (c1*c2 + c0) * s12;

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;  
  chi=acos(c)*180/MathConst::PI;

  return chi;
}

double Measure::improper(const vector<double> &x1, const vector<double> &x2, const vector<double> &x3, const vector<double> &x4) {
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2;
  double s12,c,chi;
  
  vb1x = x1[0] - x2[0];
  vb1y = x1[1] - x2[1];
  vb1z = x1[2] - x2[2];

  vb2x = x3[0] - x2[0];
  vb2y = x3[1] - x2[1];
  vb2z = x3[2] - x2[2];

  vb3x = x4[0] - x3[0];
  vb3y = x4[1] - x3[1];
  vb3z = x4[2] - x3[2];

  ss1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
  ss2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
  ss3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);

  r1 = sqrt(ss1);
  r2 = sqrt(ss2);
  r3 = sqrt(ss3);

  // sin and cos of angle

  c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3;
  c1 = (vb1x * vb2x + vb1y * vb2y + vb1z * vb2z) * r1 * r2;
  c2 = -(vb3x * vb2x + vb3y * vb2y + vb3z * vb2z) * r3 * r2;

  s1 = 1.0 - c1*c1;
  if (s1 < MathConst::EPS) s1 = MathConst::EPS;
  s1 = 1.0 / s1;

  s2 = 1.0 - c2*c2;
  if (s2 < MathConst::EPS) s2 = MathConst::EPS;
  s2 = 1.0 / s2;

  s12 = sqrt(s1*s2);
  c = (c1*c2 + c0) * s12;

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;  
  chi=acos(c)*180/MathConst::PI;

  return chi;
}

double Measure::improper(const Atom &atom_i, const Atom &atom_j, const Atom &atom_k, const Atom &atom_l) {
  vector<double> x_i=atom_i.get_atom_xyz();
  vector<double> x_j=atom_j.get_atom_xyz();
  vector<double> x_k=atom_k.get_atom_xyz();
  vector<double> x_l=atom_l.get_atom_xyz();
  
  double vb1x = x_i[0] - x_j[0];
  double vb1y = x_i[1] - x_j[1];
  double vb1z = x_i[2] - x_j[2];

  double vb2x = x_k[0] - x_j[0];
  double vb2y = x_k[1] - x_j[1];
  double vb2z = x_k[2] - x_j[2];

  double vb3x = x_l[0] - x_k[0];
  double vb3y = x_l[1] - x_k[1];
  double vb3z = x_l[2] - x_k[2];

  double ss1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
  double ss2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
  double ss3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);

  double r1 = sqrt(ss1);
  double r2 = sqrt(ss2);
  double r3 = sqrt(ss3);

  // sin and cos of angle

  double c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3;
  double c1 = (vb1x * vb2x + vb1y * vb2y + vb1z * vb2z) * r1 * r2;
  double c2 = -(vb3x * vb2x + vb3y * vb2y + vb3z * vb2z) * r3 * r2;

  double s1 = 1.0 - c1*c1;
  if (s1 < MathConst::EPS) s1 = MathConst::EPS;
  s1 = 1.0 / s1;

  double s2 = 1.0 - c2*c2;
  if (s2 < MathConst::EPS) s2 = MathConst::EPS;
  s2 = 1.0 / s2;

  double s12 = sqrt(s1*s2);
  double c = (c1*c2 + c0) * s12;

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;  
  double chi=acos(c)*180/MathConst::PI;

  return chi;
}

double Measure::point_to_line_distance(double xp[3], double xl1[3], double xl2[3]) {
  double theta=Measure::angle(xp,xl1,xl2);
  double rpl1=Measure::distance(xp,xl1);
  double pdis=rpl1*sin(theta*MathConst::PI/180);

  return pdis;
}

double Measure::point_to_line_distance(const vector<double> &xp, const vector<double> &xl1, const vector<double> &xl2) {
  double theta=Measure::angle(xp,xl1,xl2);
  double rpl1=Measure::distance(xp,xl1);
  double pdis=rpl1*sin(theta*MathConst::PI/180);

  return pdis;
}

pair<double,double> Measure::statistics(const vector<double> &Vx) {
  double aver=0.0;
  int num=(int)Vx.size();
  for(int ii=0;ii<num;ii++)
    aver+=Vx[ii];
  if(num>0)
    aver/=1.0*num;

  double std=0.0;
  for(int ii=0;ii<num;ii++) 
    std+=(Vx[ii]-aver)*(Vx[ii]-aver);
  if(num>1)  
    std=sqrt(std*1.0/(num-1));

  pair<double,double> p_one=make_pair(aver,std);
  
  return p_one;
}

pair<double,double> Measure::statistics(const vector<int> &Vx) {
  double aver=0.0;
  int num=(int)Vx.size();
  for(int ii=0;ii<num;ii++)
    aver+=Vx[ii];
  if(num>0)
    aver/=1.0*num;

  double std=0.0;
  for(int ii=0;ii<num;ii++) 
    std+=(Vx[ii]-aver)*(Vx[ii]-aver);
  if(num>1)  
    std=sqrt(std*1.0/(num-1));

  pair<double,double> p_one=make_pair(aver,std);
  
  return p_one;
}

double Measure::average(const vector<double> &Vx) {
  double aver=0.0;
  int num=(int)Vx.size();
  for(int ii=0;ii<num;ii++)
    aver+=Vx[ii];
  aver/=1.0*num;

  return aver;
}

double Measure::average(const vector<int> &Vx) {
  double aver=0.0;
  int num=(int)Vx.size();
  for(int ii=0;ii<num;ii++)
    aver+=Vx[ii];
  aver=aver*1.0/num;

  return aver;
}

double Measure::standard_deviation(const vector<double> &Vx) {
  double aver=0.0;
  int num=(int)Vx.size();
  for(int ii=0;ii<num;ii++)
    aver+=Vx[ii];
  aver=aver*1.0/num;

  double std=0.0;
  for(int ii=0;ii<num;ii++) 
    std+=(Vx[ii]-aver)*(Vx[ii]-aver);
  
  std=sqrt(std*1.0/(num-1));

  return std;
}

double Measure::standard_deviation(const vector<int> &Vx) {
  double aver=0.0;
  int num=(int)Vx.size();
  for(int ii=0;ii<num;ii++)
    aver+=Vx[ii];
  aver=aver*1.0/num;

  double std=0.0;
  for(int ii=0;ii<num;ii++) 
    std+=(Vx[ii]-aver)*(Vx[ii]-aver);
  
  std=sqrt(std*1.0/(num-1));

  return std;
}

double Measure::standard_deviation(const vector<double> &Vx, double aver) {
  int num=(int)Vx.size();
  double std=0.0;
  for(int ii=0;ii<num;ii++) 
    std+=(Vx[ii]-aver)*(Vx[ii]-aver);
  
  std=sqrt(std*1.0/(num-1));

  return std;
}

double Measure::standard_deviation(const vector<int> &Vx, double aver) {
  int num=(int)Vx.size();
  double std=0.0;
  for(int ii=0;ii<num;ii++) 
    std+=(Vx[ii]-aver)*(Vx[ii]-aver);
  
  std=sqrt(std*1.0/(num-1));

  return std;
}

vector<double> Measure::probability_distribution(const vector<double> &Vx, double delx) {
  vector<double> Vdist;
  size_t num=Vx.size();
  if(num==0) {
    string str="measure.cpp (probability_distribution): Vx.size()==0";
    cout<<"Warning: "<<str<<endl;
    return Vdist;
  }

  double xmax=Vx[0], xmin=Vx[0];
  for(size_t ii=1;ii<num;ii++) {
    if(Vx[ii]>xmax) xmax=Vx[ii];
    if(Vx[ii]<xmin) xmin=Vx[ii];
  }
  double xhi=ceil(xmax/delx)*delx;
  double xlo=floor(xmin/delx)*delx;
  int nbin=round((xhi-xlo)/delx)+1;
  Vdist.resize(nbin+1);
  Vdist[0]=xlo+delx/2.0; //x_start
  for(int ii=1;ii<=nbin;ii++)
    Vdist[ii]=0.0;

  for(size_t ii=0;ii<num;ii++) {
    for(int jj=1;jj<=nbin;jj++) {
      if(Vx[ii]>=xlo+(jj-1)*delx && Vx[ii]<xlo+jj*delx)
        Vdist[jj]=Vdist[jj]+1.0;
    }
  }
  for(int jj=1;jj<=nbin;jj++)
    Vdist[jj]=Vdist[jj]*1.0/num;
  
  return Vdist;
}

vector<double> Measure::probability_distribution(const vector<double> &Vx, double delx, double xlo, double xhi) {
  vector<double> Vdist;
  size_t num=Vx.size();
  if(num==0) {
    string str="measure.cpp (probability_distribution): Vx.size()==0";
    cout<<"Warning: "<<str<<endl;
    return Vdist;
  }

  int nbin=ceil((xhi-xlo)*1.0/delx)+2; //one bin for <xlo & one bin for >xhi;
  Vdist.resize(nbin+1);
  Vdist[0]=xlo-delx/2.0; //x_start
  for(int ii=1;ii<=nbin;ii++)
    Vdist[ii]=0.0;
  
  for(size_t ii=0;ii<num;ii++) {
    if(Vx[ii]<xlo) Vdist[1]=Vdist[1]+1.0;
    else if(Vx[ii]>=xhi) Vdist[nbin]=Vdist[nbin]+1.0;
    else {
      for(int jj=2;jj<nbin;jj++) {
        if(Vx[ii]>=xlo+(jj-2)*delx && Vx[ii]<xlo+(jj-1)*delx)
          Vdist[jj]=Vdist[jj]+1.0;
      }
    }
  }
  for(int jj=1;jj<=nbin;jj++)
    Vdist[jj]=Vdist[jj]*1.0/num;
 
  return Vdist;
}

vector<double> Measure::probability_distribution_limited(const vector<double> &Vx, double delx, double xlo, double xhi) {
  vector<double> Vdist;
  size_t num=Vx.size();
  if(num==0) {
    string str="measure.cpp (probability_distribution_limited): Vx.size()==0";
    cout<<"Warning: "<<str<<endl;
    return Vdist;
  }

  int nbin=ceil((xhi-xlo)*1.0/delx); //Limited range [xlo,xhi) 
  Vdist.resize(nbin+1);
  Vdist[0]=xlo+delx/2.0; //x_start
  for(int ii=1;ii<=nbin;ii++)
    Vdist[ii]=0.0;
  
  int ncount=0;
  for(size_t ii=0;ii<num;ii++) {
    for(int jj=1;jj<=nbin;jj++) {
      if(Vx[ii]>=xlo+(jj-1)*delx && Vx[ii]<xlo+jj*delx) {
        Vdist[jj]=Vdist[jj]+1.0;
        ncount++;
      }
    }
  }
  for(int jj=1;jj<=nbin;jj++)
    Vdist[jj]=Vdist[jj]*1.0/ncount;
 
  return Vdist;
}

vector<double> Measure::probability_distribution(const vector<int> &Vx, double delx) {
  vector<double> Vdist;
  size_t num=Vx.size();
  if(num==0) {
    string str="measure.cpp (probability_distribution): Vx.size()==0";
    cout<<"Warning: "<<str<<endl;
    return Vdist;
  }

  double xmax=Vx[0], xmin=Vx[0];
  for(size_t ii=1;ii<num;ii++) {
    if(Vx[ii]>xmax) xmax=Vx[ii];
    if(Vx[ii]<xmin) xmin=Vx[ii];
  }
  double xhi=xmax+delx/2.0;
  double xlo=xmin-delx/2.0;
  int nbin=round((xhi-xlo)/delx);
  Vdist.resize(nbin+1);
  Vdist[0]=xlo+delx/2.0; //x_start
  for(int ii=1;ii<=nbin;ii++)
    Vdist[ii]=0.0;

  for(size_t ii=0;ii<num;ii++) {
    for(int jj=1;jj<=nbin;jj++) {
      if(Vx[ii]>=xlo+(jj-1)*delx && Vx[ii]<xlo+jj*delx)
        Vdist[jj]=Vdist[jj]+1.0;
    }
  }
  for(int jj=1;jj<=nbin;jj++)
    Vdist[jj]=Vdist[jj]*1.0/num;
  
  return Vdist;
}

vector<double> Measure::probability_distribution(const vector<int> &Vx, double delx, double xlo, double xhi) {
  vector<double> Vdist;
  size_t num=Vx.size();
  if(num==0) {
    string str="measure.cpp (probability_distribution): Vx.size()==0";
    cout<<"Warning: "<<str<<endl;
    return Vdist;
  }

  int nbin=ceil((xhi-xlo)*1.0/delx)+2; //one bin for <xlo & one bin for >xhi;
  Vdist.resize(nbin+1);
  Vdist[0]=xlo-delx/2.0; //x_start
  for(int ii=1;ii<=nbin;ii++)
    Vdist[ii]=0.0;
  
  for(size_t ii=0;ii<num;ii++) {
    if(Vx[ii]<xlo) Vdist[1]=Vdist[1]+1.0;
    else if(Vx[ii]>=xhi) Vdist[nbin]=Vdist[nbin]+1.0;
    else {
      for(int jj=2;jj<nbin;jj++) {
        if(Vx[ii]>=xlo+(jj-2)*delx && Vx[ii]<xlo+(jj-1)*delx)
          Vdist[jj]=Vdist[jj]+1.0;
      }
    }
  }
  for(int jj=1;jj<=nbin;jj++)
    Vdist[jj]=Vdist[jj]*1.0/num;
 
  return Vdist;
}

vector<double> Measure::probability_distribution_limited(const vector<int> &Vx, double delx, double xlo, double xhi) {
  vector<double> Vdist;
  size_t num=Vx.size();
  if(num==0) {
    string str="measure.cpp (probability_distribution_limited): Vx.size()==0";
    cout<<"Warning: "<<str<<endl;
    return Vdist;
  }

  int nbin=ceil((xhi-xlo)*1.0/delx); //Limited range [xlo,xhi) 
  Vdist.resize(nbin+1);
  Vdist[0]=xlo+delx/2.0; //x_start
  for(int ii=1;ii<=nbin;ii++)
    Vdist[ii]=0.0;
  
  int ncount=0;
  for(size_t ii=0;ii<num;ii++) {
    for(int jj=1;jj<=nbin;jj++) {
      if(Vx[ii]>=xlo+(jj-1)*delx && Vx[ii]<xlo+jj*delx) {
        Vdist[jj]=Vdist[jj]+1.0;
        ncount++;
      }
    }
  }
  for(int jj=1;jj<=nbin;jj++)
    Vdist[jj]=Vdist[jj]*1.0/ncount;
 
  return Vdist;
}

double Measure::Pearson_correlation_coefficient(const vector<double> &Vx, const vector<double> &Vy) {
  int num=(int)Vx.size();
  if(num != (int)Vy.size()) {
    string str="measure.cpp (Pearson_correlation_coefficient): Numbers of data do not matched --> "+to_string(num)+" : "+to_string(Vy.size());
    cout<<"Error: "<<str<<endl;
    exit(0);
  }

  double Ave1=0.0, Ave2=0.0;
  for(int jj=0;jj<num;jj++) {
    Ave1+=Vx[jj];
    Ave2+=Vy[jj];
  }
  Ave1=Ave1/(1.0*num);
  Ave2=Ave2/(1.0*num);
  double SUMxx=0.0, SUMyy=0.0, SUMxy=0.0;
  for(int jj=0;jj<num;jj++){
    SUMxx+=(Vx[jj]-Ave1)*(Vx[jj]-Ave1);
    SUMyy+=(Vy[jj]-Ave2)*(Vy[jj]-Ave2);
    SUMxy+=(Vx[jj]-Ave1)*(Vy[jj]-Ave2);
  }
  double cc=SUMxy/sqrt(SUMxx*SUMyy);

  return cc;
}

double Measure::Pearson_correlation_coefficient(const vector<int> &Vx, const vector<int> &Vy) {
  int num=(int)Vx.size();
  if(num != (int)Vy.size()) {
    string str="measure.cpp (Pearson_correlation_coefficient): Numbers of data do not matched --> "+to_string(num)+" : "+to_string(Vy.size());
    cout<<"Error: "<<str<<endl;
    exit(0);
  }

  double Ave1=0.0, Ave2=0.0;
  for(int jj=0;jj<num;jj++) {
    Ave1+=Vx[jj];
    Ave2+=Vy[jj];
  }
  Ave1=Ave1*1.0/num;
  Ave2=Ave2*1.0/num;
  double SUMxx=0.0, SUMyy=0.0, SUMxy=0.0;
  for(int jj=0;jj<num;jj++){
    SUMxx+=(Vx[jj]-Ave1)*(Vx[jj]-Ave1);
    SUMyy+=(Vy[jj]-Ave2)*(Vy[jj]-Ave2);
    SUMxy+=(Vx[jj]-Ave1)*(Vy[jj]-Ave2);
  }
  double cc=SUMxy*1.0/sqrt(SUMxx*SUMyy);

  return cc;
}

double Measure::Spearman_rank_correlation_coefficient(const vector<double> &Vx, const vector<double> &Vy) {
  int num=(int)Vx.size();
  if(num != (int)Vy.size()) {
    string str="measure.cpp (Spearman_rank_correlation_coefficient): Numbers of data do not matched --> "+to_string(num)+" : "+to_string(Vy.size());
    cout<<"Error: "<<str<<endl;
    exit(0);
  }

  vector<double> Vrank_x,Vrank_y;
  for(int ii=0;ii<num;ii++) {
    int r=1;
    int s=1;
    for(int jj=0;jj<ii;jj++) {
      if(Vx[jj]<Vx[ii]) r++;
      if(fabs(Vx[jj]-Vx[ii])<=MathConst::EPS) s++;
    }
    for(int jj=ii+1;jj<num;jj++) {
      if(Vx[jj]<Vx[ii]) r++;
      if(fabs(Vx[jj]-Vx[ii])<=MathConst::EPS) s++;
    }
    double ranki=r+(s-1)*0.5;
    Vrank_x.push_back(ranki);
  }
  for(int ii=0;ii<num;ii++) {
    int r=1;
    int s=1;
    for(int jj=0;jj<ii;jj++) {
      if(Vy[jj]<Vy[ii]) r++;
      if(fabs(Vy[jj]-Vy[ii])<=MathConst::EPS) s++;
    }
    for(int jj=ii+1;jj<num;jj++) {
      if(Vy[jj]<Vy[ii]) r++;
      if(fabs(Vy[jj]-Vy[ii])<=MathConst::EPS) s++;
    }
    double ranki=r+(s-1)*0.5;
    Vrank_y.push_back(ranki);
  }
  
  double Ave1=0.0, Ave2=0.0;
  for(int ii=0;ii<num;ii++) {
    Ave1+=Vrank_x[ii];
    Ave2+=Vrank_y[ii];
  }
  Ave1=Ave1*1.0/num;
  Ave2=Ave2*1.0/num;
  double SUMxx=0.0, SUMyy=0.0, SUMxy=0.0;
  for(int ii=0;ii<num;ii++){
    SUMxx+=(Vrank_x[ii]-Ave1)*(Vrank_x[ii]-Ave1);
    SUMyy+=(Vrank_y[ii]-Ave2)*(Vrank_y[ii]-Ave2);
    SUMxy+=(Vrank_x[ii]-Ave1)*(Vrank_y[ii]-Ave2);
  }
  double spcc=SUMxy*1.0/sqrt(SUMxx*SUMyy);

  return spcc;
}

double Measure::Spearman_rank_correlation_coefficient(const vector<int> &Vx, const vector<int> &Vy) {
  int num=(int)Vx.size();
  if(num != (int)Vy.size()) {
    string str="measure.cpp (Spearman_rank_correlation_coefficient): Numbers of data do not matched --> "+to_string(num)+" : "+to_string(Vy.size());
    cout<<"Error: "<<str<<endl;
    exit(0);
  }

  vector<double> Vrank_x,Vrank_y;
  for(int ii=0;ii<num;ii++) {
    int r=1;
    int s=1;
    for(int jj=0;jj<ii;jj++) {
      if(Vx[jj]<Vx[ii]) r++;
      if(Vx[jj]==Vx[ii]) s++;
    }
    for(int jj=ii+1;jj<num;jj++) {
      if(Vx[jj]<Vx[ii]) r++;
      if(Vx[jj]==Vx[ii]) s++;
    }
    double ranki=r+(s-1)*0.5;
    Vrank_x.push_back(ranki);
  }
  for(int ii=0;ii<num;ii++) {
    int r=1;
    int s=1;
    for(int jj=0;jj<ii;jj++) {
      if(Vy[jj]<Vy[ii]) r++;
      if(Vy[jj]==Vy[ii]) s++;
    }
    for(int jj=ii+1;jj<num;jj++) {
      if(Vy[jj]<Vy[ii]) r++;
      if(Vy[jj]==Vy[ii]) s++;
    }
    double ranki=r+(s-1)*0.5;
    Vrank_y.push_back(ranki);
  }
  
  double Ave1=0.0, Ave2=0.0;
  for(int ii=0;ii<num;ii++) {
    Ave1+=Vrank_x[ii];
    Ave2+=Vrank_y[ii];
  }
  Ave1=Ave1*1.0/num;
  Ave2=Ave2*1.0/num;
  double SUMxx=0.0, SUMyy=0.0, SUMxy=0.0;
  for(int ii=0;ii<num;ii++){
    SUMxx+=(Vrank_x[ii]-Ave1)*(Vrank_x[ii]-Ave1);
    SUMyy+=(Vrank_y[ii]-Ave2)*(Vrank_y[ii]-Ave2);
    SUMxy+=(Vrank_x[ii]-Ave1)*(Vrank_y[ii]-Ave2);
  }
  double spcc=SUMxy*1.0/sqrt(SUMxx*SUMyy);

  return spcc;
}

pair<double,double> Measure::linear_fitting(const vector<double> &Vx, const vector<double> &Vy) {
  int num=(int)Vx.size();
  if(num != (int)Vy.size()) {
    string str="measure.cpp (linear_fitting): Numbers of data do not matched --> "+to_string(num)+" : "+to_string(Vy.size());
    cout<<"Error: "<<str<<endl;
    exit(0);
  }

  double xcm=0.0, ycm=0.0;
  for(int ii=0;ii<num;ii++) {
    xcm+=Vx[ii];
    ycm+=Vy[ii];
  }
  xcm=xcm*1.0/num;
  ycm=ycm*1.0/num;
  double SUMxy=0.0, SUMxsq=0.0;
  for(int ii=0;ii<num;ii++) {
    SUMxy+=(Vx[ii]-xcm)*(Vy[ii]-ycm);
    SUMxsq+=(Vx[ii]-xcm)*(Vx[ii]-xcm);
  }
  double k=SUMxy*1.0/SUMxsq;
  double b=ycm-k*xcm;

  pair<double,double> p_one=make_pair(k,b);

  return p_one;
}

pair<double,double> Measure::linear_fitting(const vector<int> &Vx, const vector<int> &Vy) {
  int num=(int)Vx.size();
  if(num != (int)Vy.size()) {
    string str="measure.cpp (linear_fitting): Numbers of data do not matched --> "+to_string(num)+" : "+to_string(Vy.size());
    cout<<"Error: "<<str<<endl;
    exit(0);
  }

  double xcm=0.0, ycm=0.0;
  for(int ii=0;ii<num;ii++) {
    xcm+=Vx[ii];
    ycm+=Vy[ii];
  }
  xcm=xcm*1.0/num;
  ycm=ycm*1.0/num;
  double SUMxy=0.0, SUMxsq=0.0;
  for(int ii=0;ii<num;ii++) {
    SUMxy+=(Vx[ii]-xcm)*(Vy[ii]-ycm);
    SUMxsq+=(Vx[ii]-xcm)*(Vx[ii]-xcm);
  }
  double k=SUMxy*1.0/SUMxsq;
  double b=ycm-k*xcm;
  
  pair<double,double> p_one=make_pair(k,b);

  return p_one;
}

double Measure::root_mean_square_error(const vector<double> &Vx, const vector<double> &Vy) {
  int num=(int)Vx.size();
  if(num != (int)Vy.size()) {
    string str="measure.cpp (root_mean_square_error): Numbers of data do not matched --> "+to_string(num)+" : "+to_string(Vy.size());
    cout<<"Error: "<<str<<endl;
    exit(0);
  }

  double mse=0.0;
  for(int ii=0;ii<num;ii++)
    mse+=(Vx[ii]-Vy[ii])*(Vx[ii]-Vy[ii]);
  double rmse=sqrt(mse*1.0/num);

  return rmse;
}

double Measure::root_mean_square_error(const vector<int> &Vx, const vector<int> &Vy) {
  int num=(int)Vx.size();
  if(num != (int)Vy.size()) {
    string str="measure.cpp (root_mean_square_error): Numbers of data do not matched --> "+to_string(num)+" : "+to_string(Vy.size());
    cout<<"Error: "<<str<<endl;
    exit(0);
  }

  double mse=0.0;
  for(int ii=0;ii<num;ii++)
    mse+=(Vx[ii]-Vy[ii])*(Vx[ii]-Vy[ii]);
  double rmse=sqrt(mse*1.0/num);

  return rmse;
}

double Measure::cross_entropy(const vector<double> &Vp, const vector<double> &Vq) {  //H(p,q)=Sum{-p(x)*log[q(x)]}
  size_t num=Vp.size();
  if(num != Vq.size()) {
    string str="measure.cpp (cross_entropy): Numbers of data do not matched --> "+to_string(num)+" : "+to_string(Vp.size());
    cout<<"Error: "<<str<<endl;
    exit(0);
  }

  double loss=0.0;
  for(size_t ii=0;ii<num;ii++) {
    double q_i=Vq[ii];
    if(q_i<MathConst::EPS) q_i=MathConst::EPS;
    loss-=Vp[ii]*log(q_i);
  }
    
  return loss;
}

double Measure::relative_entropy(const vector<double> &Vp, const vector<double> &Vq) { //D(p||q)=Sum{p(x)*[log(p(x)/q(x))]}
  size_t num=Vp.size();
  if(num != Vq.size()) {
    string str="measure.cpp (relative_entropy): Numbers of data do not matched --> "+to_string(num)+" : "+to_string(Vq.size());
    cout<<"Error: "<<str<<endl;
    exit(0);
  }

  double loss=0.0;
  for(size_t ii=0;ii<num;ii++) {
    double q_i=Vq[ii];
    if(q_i<MathConst::EPS) q_i=MathConst::EPS;
    double ratio_i=Vp[ii]/q_i;
    if(ratio_i<MathConst::EPS) ratio_i=MathConst::EPS;
    loss+=Vp[ii]*log(ratio_i);
  }
    
  return loss;
}

int Measure::vector_relation(vector<string> Vec1, vector<string> Vec2) {
  if(Vec1.size()<Vec2.size())
    swap(Vec1,Vec2);
    
  size_t num=0;
  for(size_t ii=0;ii<Vec2.size();ii++) {
    if(count(Vec1.begin(),Vec1.end(),Vec2[ii]))
      num++;
  }
  
  if(num==0) return Measure::separated;
  else if(num>0 && num<Vec2.size()) return Measure::intersected;
  else if(num==Vec2.size() && Vec2.size()==Vec1.size()) return Measure::equal;
  else return Measure::contained;
}

int Measure::vector_relation(vector<int> Vec1, vector<int> Vec2) {
  if(Vec1.size()<Vec2.size())
    swap(Vec1,Vec2);
    
  size_t num=0;
  for(size_t ii=0;ii<Vec2.size();ii++) {
    if(count(Vec1.begin(),Vec1.end(),Vec2[ii]))
      num++;
  }
  
  if(num==0) return Measure::separated;
  else if(num>0 && num<Vec2.size()) return Measure::intersected;
  else if(num==Vec2.size() && Vec2.size()==Vec1.size()) return Measure::equal;
  else return Measure::contained;
}

ConfusionMatrixScore Measure::confusion_matrix_score(const vector<string> &Vpred, const vector<string> &Vtrue) {
  ConfusionMatrixScore score;
  for(size_t ii=0;ii<Vpred.size();ii++) {
    if(count(Vtrue.begin(),Vtrue.end(),Vpred[ii]))
      score.TP+=1;
    else 
      score.FP+=1;
  }
  for(size_t ii=0;ii<Vtrue.size();ii++) {
    if(count(Vpred.begin(),Vpred.end(),Vtrue[ii])==0)
      score.FN+=1;
  }

  if(score.TP+score.FP>0)
    score.PPV=score.TP*1.0/(score.TP+score.FP);
  
  if(score.TP+score.FN>0)
    score.STY=score.TP*1.0/(score.TP+score.FN);

  // RNA-Puzzles INF definition
  if(score.TP==0 && (score.FP==0 || score.FN==0))
    score.INF=-1.0;
  else
    score.INF=sqrt(score.PPV*score.STY);

  if(score.PPV+score.STY>0)
    score.F1=2*score.PPV*score.STY/(score.PPV+score.STY);

  return score;
}
