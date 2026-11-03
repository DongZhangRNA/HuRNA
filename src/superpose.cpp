#include "basic.h"
#include <Eigen/Dense>
#include "superpose.h"

using namespace CBS_NS;
using namespace Eigen;

//const double RMSD_phais=0.00001;

Superpose::Superpose() {
  x1.clear();
  x2.clear();
  lenx=0;
  for(int ii=0;ii<3;ii++) {
    x1cen.push_back(0.0);
    x2cen.push_back(0.0);
    vector<double> vm;
    for(int jj=0;jj<3;jj++) {
      if(ii==jj) vm.push_back(0.0);
      else vm.push_back(0.0);
    }
    M.push_back(vm);
  }
}

Superpose::Superpose(const vector<double> &Vx1, const vector<double> &Vx2) {
  Superpose::set_x1x2(Vx1,Vx2);
}

Superpose::~Superpose() {
  lenx=-1;
  M.clear();
  M.shrink_to_fit();
  x1cen.clear();
  x1cen.shrink_to_fit();
  x2cen.clear();
  x2cen.shrink_to_fit();
  x1.clear();
  x1.shrink_to_fit();
  x2.clear();
  x2.shrink_to_fit();
}

void Superpose::set_x1x2(const vector<double> &Vx1, const vector<double> &Vx2) {
  if(Vx1.size()!=Vx2.size()) {
    string str="superposition.cpp: Numbers of data are not matched --> "+to_string(Vx1.size()/3)+"/"+to_string(Vx2.size()/3);
    cout<<"Error in "<<str<<endl;
    exit(0);
  }
  if(Vx1.size()%3!=0) {
    string str="superposition.cpp: Number of data is wrong in x1 --> "+to_string(Vx1.size());
    cout<<"Error in "<<str<<endl;
    exit(0);
  }
  if(Vx2.size()%3!=0) {
    string str="superposition.cpp: Number of data is wrong in x2 --> "+to_string(Vx2.size());
    cout<<"Error in "<<str<<endl;
    exit(0);
  }

  lenx=Vx1.size()/3;
  if(lenx<2) {
    string str="superposition.cpp: Number of data < 2 --> "+to_string(lenx);
    cout<<"Error in "<<str<<endl;
    exit(0);
  }

  for(int ii=0;ii<3;ii++) {
    x1cen.push_back(0.0);
    x2cen.push_back(0.0);
    vector<double> vm;
    for(int jj=0;jj<3;jj++) {
      if(ii==jj) vm.push_back(0.0);
      else vm.push_back(0.0);
    }
    M.push_back(vm);
  }

  x1.clear();
  x2.clear();
  for(int ii=0;ii<lenx;ii++) {
    vector<double> v1,v2;
    for(int jj=0;jj<3;jj++) {
      v1.push_back(Vx1[ii*3+jj]);
      v2.push_back(Vx2[ii*3+jj]);
      x1cen[jj]+=Vx1[ii*3+jj];
      x2cen[jj]+=Vx2[ii*3+jj];
    }
    x1.push_back(v1);
    x2.push_back(v2);
  }
  for(int jj=0;jj<3;jj++) {
    x1cen[jj]/=lenx*1.0;
    x2cen[jj]/=lenx*1.0;
  }
}

void Superpose::set_x1(const vector<double> &Vx1) {
  if(Vx1.size()%3!=0) {
    string str="superposition.cpp: Number of data is wrong in x1 --> "+to_string(Vx1.size());
    cout<<"Error in "<<str<<endl;
    exit(0);
  }

  lenx=Vx1.size()/3;
  if(lenx<2) {
    string str="superposition.cpp: Number of data < 2 --> "+to_string(lenx);
    cout<<"Error in "<<str<<endl;
    exit(0);
  }

  for(int ii=0;ii<3;ii++) {
    x1cen[ii]=0.0;
  }

  x1.clear();
  for(int ii=0;ii<lenx;ii++) {
    vector<double> v1;
    for(int jj=0;jj<3;jj++) {
      v1.push_back(Vx1[ii*3+jj]);
      x1cen[jj]+=Vx1[ii*3+jj];
    }
    x1.push_back(v1);
  }
  for(int jj=0;jj<3;jj++) {
    x1cen[jj]/=lenx*1.0;
  }
}

void Superpose::set_x2(const vector<double> &Vx2) {
  if(Vx2.size()%3!=0) {
    string str="superposition.cpp: Number of data is wrong in x2 --> "+to_string(Vx2.size());
    cout<<"Error in "<<str<<endl;
    exit(0);
  }

  if(lenx!=(int)Vx2.size()/3) {
    string str="superposition.cpp: Numbers of data are not matched --> "+to_string(lenx)+"/"+to_string(Vx2.size()/3);
    cout<<"Error in "<<str<<endl;
    exit(0);
  }
  
  for(int ii=0;ii<3;ii++) {
    x2cen[ii]=0.0;
    for(int jj=0;jj<3;jj++)
      M[ii][jj]=0.0;
  }

  x2.clear();
  for(int ii=0;ii<lenx;ii++) {
    vector<double> v2;
    for(int jj=0;jj<3;jj++) {
      v2.push_back(Vx2[ii*3+jj]);
      x2cen[jj]+=Vx2[ii*3+jj];
    }
    x2.push_back(v2);
  }
  for(int jj=0;jj<3;jj++) {
    x2cen[jj]/=lenx*1.0;
  }
}

void Superpose::calculate_rotation_matrix() {  //Using Kabsch algorithm to calculate the optimal rotation matrix
  Eigen::Matrix3d H;
  for(int ii=0;ii<3;ii++) {
    for(int jj=0;jj<3;jj++) {
      H(ii,jj)=0.0;
      for(int kk=0;kk<lenx;kk++) {
        H(ii,jj)+=(x2[kk][ii]-x2cen[ii])*(x1[kk][jj]-x1cen[jj]);
      }
    }
  }
  Eigen::JacobiSVD<Matrix3d> svd(H, ComputeFullU | ComputeFullV);
  Eigen::Matrix3d U = svd.matrixU();
  Eigen::Matrix3d V = svd.matrixV();

  Eigen::Matrix3d R = V * U.transpose();
  if (R.determinant() < 0) {
    V.col(2) *= -1;
    R = V * U.transpose();
  }

  for(int ii=0;ii<3;ii++)
    for(int jj=0;jj<3;jj++)
      M[ii][jj]=R(ii,jj);
}

double Superpose::get_RMSD() {
  Superpose::calculate_rotation_matrix();
  double rms=0.0;
  for(int kk=0;kk<lenx;kk++) {
    for(int ii=0;ii<3;ii++) {
      double yplus=0;
      for(int jj=0;jj<3;jj++)
        yplus+=M[ii][jj]*(x2[kk][jj]-x2cen[jj]);
      rms+=((x1[kk][ii]-x1cen[ii])-yplus)*((x1[kk][ii]-x1cen[ii])-yplus);
    }
  }
  rms=sqrt(rms/lenx);

  return rms;
}

/*
double Superpose::get_RMSD() {
  double phai[3] = {1.0}, seta[3] = {1.0}, delta[3] = {1.0};
  double setap = 0.0, deltap = 0.0, mx[3] = {0.0}, my[3] = {0.0}, mz[3] = {0.0};
  double C[3][3] = {0.0}, C2[3][3] = {0.0};

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++) {
      C[i][j]=0;
      for(int k=0;k<lenx;k++)
        C[i][j]+=(x1[k][j]-x1cen[j])*(x2[k][i]-x2cen[i]);
    }

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++) {
      C2[i][j]=0;
      for(int k=0;k<3;k++)
        C2[i][j]+=M[i][k]*C[k][j];
    }

  double rms=RMSD_rmsdx();

  int p = 0, q = 0, r = 0;
  while(fabs(phai[0])>RMSD_phais || fabs(phai[1])>RMSD_phais || fabs(phai[2])>RMSD_phais) {
    if(p==0) {q=1; r=2;}
    else if(p==1) {q=2; r=0;}
    else if(p==2) {q=0; r=1;}

    setap=C2[r][q]-C2[q][r];
    deltap=C2[q][q]+C2[r][r];

    for(int k=0;k<3;k++) {
      mx[k]=M[p][k];
      my[k]=(deltap*M[q][k]+setap*M[r][k])/sqrt(deltap*deltap+setap*setap);
      mz[k]=(-setap*M[q][k]+deltap*M[r][k])/sqrt(deltap*deltap+setap*setap);
    }

    for(int k=0;k<3;k++) {
      M[p][k]=mx[k];
      M[q][k]=my[k];
      M[r][k]=mz[k];
    }

    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++) {
        C2[i][j]=0;
        for(int k=0;k<3;k++)
          C2[i][j]+=M[i][k]*C[k][j];
      }

    rms=RMSD_rmsdx();

    seta[0]=C2[2][1]-C2[1][2];
    delta[0]=C2[1][1]+C2[2][2];
    seta[1]=C2[0][2]-C2[2][0];
    delta[1]=C2[2][2]+C2[0][0];
    seta[2]=C2[1][0]-C2[0][1];
    delta[2]=C2[0][0]+C2[1][1];

    for(int k=0;k<3;k++)
      phai[k]=atan(seta[k]/delta[k]);

    if(p==2) p=0;
    else p++;
  }

  return rms;
}
*/

void Superpose::reset_x2(const vector<double> &Vx2) {
  if(Vx2.size()%3!=0) {
    string str="superposition.cpp: Number of data is wrong in x2 --> "+to_string(Vx2.size());
    cout<<"Error in "<<str<<endl;
    exit(0);
  }

  lenx=Vx2.size()/3;
  x2.clear();
  for(int ii=0;ii<lenx;ii++) {
    vector<double> v2;
    for(int jj=0;jj<3;jj++)
      v2.push_back(Vx2[ii*3+jj]);
    x2.push_back(v2);
  }
}

vector<double> Superpose::rotate_x2() {
  vector<double> Vx3;
  for(int kk=0;kk<lenx;kk++){
    for(int ii=0;ii<3;ii++){
      double yyplus=0.0;
      for(int jj=0;jj<3;jj++)
        yyplus+=M[ii][jj]*(x2[kk][jj]-x2cen[jj]);
      yyplus+=x1cen[ii];
      Vx3.push_back(yyplus);
    }
  }

  return Vx3;
}

vector<double> Superpose::rotate_x2(const vector<double> &Vx2) {
  if(Vx2.size()%3!=0) {
    string str="superposition.cpp: Number of data is wrong in x2 --> "+to_string(Vx2.size());
    cout<<"Error in "<<str<<endl;
    exit(0);
  }

  lenx=Vx2.size()/3;
  x2.clear();
  for(int ii=0;ii<lenx;ii++) {
    vector<double> v2;
    for(int jj=0;jj<3;jj++)
      v2.push_back(Vx2[ii*3+jj]);
    x2.push_back(v2);
  }
  
  vector<double> Vx3;
  for(int kk=0;kk<lenx;kk++){
    for(int ii=0;ii<3;ii++){
      double yyplus=0.0;
      for(int jj=0;jj<3;jj++)
        yyplus+=M[ii][jj]*(x2[kk][jj]-x2cen[jj]);
      yyplus+=x1cen[ii];
      Vx3.push_back(yyplus);
    }
  }

  return Vx3;
}

void Superpose::clear() {
  lenx=-1;
  M.clear();
  M.shrink_to_fit();
  x1cen.clear();
  x1cen.shrink_to_fit();
  x2cen.clear();
  x2cen.shrink_to_fit();
  x1.clear();
  x1.shrink_to_fit();
  x2.clear();
  x2.shrink_to_fit();
}

/*
inline double Superpose::RMSD_rmsdx() {
  double rms=0;
  for(int kk=0;kk<lenx;kk++) {
    for(int ii=0;ii<3;ii++) {
      double yplus=0;
      for(int jj=0;jj<3;jj++)
        yplus+=M[ii][jj]*(x2[kk][jj]-x2cen[jj]);
      rms+=((x1[kk][ii]-x1cen[ii])-yplus)*((x1[kk][ii]-x1cen[ii])-yplus);
    }
  }
  rms=sqrt(rms/lenx);

  return rms;
}
*/