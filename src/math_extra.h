#ifndef _MATH_EXTRA_H
#define _MATH_EXTRA_H

#include "basic.h"

namespace MathExtra {
inline void zero3(double *);
inline void zero3(int *);
inline void zero3(vector<double> &);
inline void zero3(vector<int> &);
inline void norm3(double *);
inline void norm3(vector<double> &);
inline vector<double> norm3(const vector<int> &);
inline double length3(double [3]);
inline double length3(const vector<double> &);
inline double length3(const vector<int> &);
inline void scale3(const double, double *);
inline void scale3(const double, vector<double> &);
inline void negate3(double *);
inline void negate3(vector<double> &);
inline void negate3(vector<int> &);
inline vector<double> add3(const vector<double> &, const vector<double> &);
inline vector<int> add3(const vector<int> &, const vector<int> &);
inline vector<double> sub3(const vector<double> &, const vector<double> &);
inline vector<int> sub3(const vector<int> &, const vector<int> &);
inline double dot3(double [3], double [3]);
inline double dot3(const vector<double> &, const vector<double> &);
inline int dot3(const vector<int> &, const vector<int> &);
inline vector<double> cross3(double [3], double [3]);
inline vector<double> cross3(const vector<double> &, const vector<double> &);
inline vector<int> cross3(const vector<int> &, const vector<int> &);
inline vector<double> to_double3(const vector<int> &);

inline double erf(double x);
}

inline void MathExtra::zero3(double *v) {
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 0.0;
}

inline void MathExtra::zero3(int *v) {
  v[0] = 0;
  v[1] = 0;
  v[2] = 0;  
}

inline void MathExtra::zero3(vector<double> &vec) {
  vec.clear();
  vec.resize(3);

  vec[0] = 0.0;
  vec[1] = 0.0;
  vec[2] = 0.0;
}

inline void MathExtra::zero3(vector<int> &vec) {
  vec.clear();
  vec.resize(3);
  
  vec[0] = 0;
  vec[1] = 0;
  vec[2] = 0;
}

inline void MathExtra::norm3(double *v) {
  double val = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; 
  if(val > 0.0) {
    double scale = 1.0 / sqrt(val);
    v[0] *= scale;
    v[1] *= scale;
    v[2] *= scale;
  } 
}

inline void MathExtra::norm3(vector<double> &vec) {
  double val = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]; 
  if(val > 0.0) {
    double scale = 1.0 / sqrt(val);
    vec[0] *= scale;
    vec[1] *= scale;
    vec[2] *= scale;
  }   
}

inline vector<double> MathExtra::norm3(const vector<int> &vec) {
  double val = vec[0]*vec[0]*1.0 + vec[1]*vec[1]*1.0 + vec[2]*vec[2]*1.0; 
  vector<double> Vnorm;
  if(val > 0.0) {
    double scale = 1.0 / sqrt(val);
    Vnorm.push_back(vec[0]*scale);
    Vnorm.push_back(vec[1]*scale);
    Vnorm.push_back(vec[2]*scale);
  }
  
  return Vnorm;
}

inline double MathExtra::length3(double x1[3]) {
  double rsq = x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2];
  double r = sqrt(rsq);
  
  return r;
}
  
inline double MathExtra::length3(const vector<double> &Vx1) {
  double rsq = 0.0;
  for(size_t ii=0;ii<Vx1.size();ii++)
    rsq += Vx1[ii]*Vx1[ii];
  double r = sqrt(rsq);
  
  return r;
}
  
inline double MathExtra::length3(const vector<int> &Vx1) {
  double rsq = 0.0;
  for(size_t ii=0;ii<Vx1.size();ii++)
    rsq += Vx1[ii]*Vx1[ii];
  double r = sqrt(rsq);
  
  return r;
}

inline void MathExtra::scale3(const double s, double *v) {
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

inline void MathExtra::scale3(const double s, vector<double> &vec) {
  vec[0] *= s;
  vec[1] *= s;
  vec[2] *= s;
}

inline void MathExtra::negate3(double *v) {
  v[0] = -v[0];
  v[1] = -v[1];
  v[2] = -v[2];
}

inline void MathExtra::negate3(vector<double> &vec) {
  for(size_t ii=0;ii<vec.size();ii++)
    vec[ii] = -vec[ii];
}

inline void MathExtra::negate3(vector<int> &vec) {
  for(size_t ii=0;ii<vec.size();ii++)
    vec[ii] = -vec[ii];
}

inline vector<double> MathExtra::add3(const vector<double> &vx1, const vector<double> &vx2) {
  vector<double> Vadd;
  for(size_t ii=0;ii<vx1.size();ii++)
    Vadd.push_back(vx1[ii] + vx2[ii]);
  
  return Vadd;
}

inline vector<int> MathExtra::add3(const vector<int> &vx1, const vector<int> &vx2) {
  vector<int> Vadd;
  for(size_t ii=0;ii<vx1.size();ii++)
    Vadd.push_back(vx1[ii] + vx2[ii]);
  
  return Vadd;
}

inline vector<double> MathExtra::sub3(const vector<double> &vx1, const vector<double> &vx2) {
  vector<double> Vsub;
  for(size_t ii=0;ii<vx1.size();ii++)
    Vsub.push_back(vx1[ii] - vx2[ii]);
  
  return Vsub;
}

inline vector<int> MathExtra::sub3(const vector<int> &vx1, const vector<int> &vx2) {
  vector<int> Vsub;
  for(size_t ii=0;ii<vx1.size();ii++)
    Vsub.push_back(vx1[ii] - vx2[ii]);
  
  return Vsub;
}

inline double MathExtra::dot3(double x1[3], double x2[3]) {
  double cs = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
    
  return cs; 
}
  
inline double MathExtra::dot3(const vector<double> &x1, const vector<double> &x2) {
  double cs = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
  
  return cs; 
}
  
inline int MathExtra::dot3(const vector<int> &x1, const vector<int> &x2) {
  int cs = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
  
  return cs; 
}
  
inline vector<double> MathExtra::cross3(double x1[3], double x2[3]) {
  double x3_0 = x1[1]*x2[2] - x1[2]*x2[1];
  double x3_1 = x1[2]*x2[0] - x1[0]*x2[2];
  double x3_2 = x1[0]*x2[1] - x1[1]*x2[0];
    
  vector<double> Vx3;
  Vx3.push_back(x3_0);
  Vx3.push_back(x3_1);
  Vx3.push_back(x3_2);
    
  return Vx3; 
}
  
inline vector<double> MathExtra::cross3(const vector<double> &x1, const vector<double> &x2) {
  double x3_0 = x1[1]*x2[2] - x1[2]*x2[1];
  double x3_1 = x1[2]*x2[0] - x1[0]*x2[2];
  double x3_2 = x1[0]*x2[1] - x1[1]*x2[0];
    
  vector<double> Vx3;
  Vx3.push_back(x3_0);
  Vx3.push_back(x3_1);
  Vx3.push_back(x3_2);
    
  return Vx3; 
}
  
inline vector<int> MathExtra::cross3(const vector<int> &x1, const vector<int> &x2) {
  int x3_0 = x1[1]*x2[2] - x1[2]*x2[1];
  int x3_1 = x1[2]*x2[0] - x1[0]*x2[2];
  int x3_2 = x1[0]*x2[1] - x1[1]*x2[0];
    
  vector<int> Vx3;
  Vx3.push_back(x3_0);
  Vx3.push_back(x3_1);
  Vx3.push_back(x3_2);
    
  return Vx3; 
}

inline vector<double> MathExtra::to_double3(const vector<int> &vx) {
  vector<double> vx1;
  for(size_t ii=0;ii<vx.size();ii++) 
    vx1.push_back(vx[ii]*1.0);
  
  return vx1;
}

inline double MathExtra::erf(double x) {
  double t = 1.0 / (1.0 + 0.5 * fabs(x));
  // use Horner's method
  double ans = 1 - t * exp( -x*x -  1.26551223 +
                t * ( 1.00002368 +
                t * ( 0.37409196 + 
                t * ( 0.09678418 + 
                t * (-0.18628806 + 
                t * ( 0.27886807 + 
                t * (-1.13520398 + 
                t * ( 1.48851587 + 
                t * (-0.82215223 + 
                t * ( 0.17087277))))))))));
  if(x >= 0.0)
      return ans;
  else
      return -ans;
}

#endif