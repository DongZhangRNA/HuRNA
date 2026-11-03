#include "basic.h"
#include "sort_vector.h"

using namespace CBS_NS;

SortVector::SortVector() {}

SortVector::~SortVector() {}

vector<size_t> SortVector::sort_vector_ascending(vector<double> &vec_one) {
  vector<size_t> Vsort;
  for(size_t ii=0;ii<vec_one.size();ii++)
    Vsort.push_back(ii);

  if(vec_one.size()<2) return Vsort;

  for(int ii=0;ii<(int)vec_one.size()-1;ii++) {
    for(int jj=(int)vec_one.size()-1;jj>ii;jj--) {
      if(vec_one[jj]<vec_one[jj-1]) {
        double tmp=vec_one[jj];
        vec_one[jj]=vec_one[jj-1];
        vec_one[jj-1]=tmp;
        size_t tmp1=Vsort[jj];
        Vsort[jj]=Vsort[jj-1];
        Vsort[jj-1]=tmp1;
      }
    }
  }

  return Vsort;
}

vector<size_t> SortVector::sort_vector_ascending(vector<int> &vec_one) {
  vector<size_t> Vsort;
  for(size_t ii=0;ii<vec_one.size();ii++)
    Vsort.push_back(ii);
  
  if(vec_one.size()<2) return Vsort;

  for(int ii=0;ii<(int)vec_one.size()-1;ii++) {
    for(int jj=(int)vec_one.size()-1;jj>ii;jj--) {
      if(vec_one[jj]<vec_one[jj-1]) {
        int tmp=vec_one[jj];
        vec_one[jj]=vec_one[jj-1];
        vec_one[jj-1]=tmp;
        size_t tmp1=Vsort[jj];
        Vsort[jj]=Vsort[jj-1];
        Vsort[jj-1]=tmp1;
      }
    }
  }

  return Vsort;
}

vector<size_t> SortVector::sort_vector_ascending(vector<size_t> &vec_one) {
  vector<size_t> Vsort;
  for(size_t ii=0;ii<vec_one.size();ii++)
    Vsort.push_back(ii);
  
  if(vec_one.size()<2) return Vsort;

  for(int ii=0;ii<(int)vec_one.size()-1;ii++) {
    for(int jj=(int)vec_one.size()-1;jj>ii;jj--) {
      if(vec_one[jj]<vec_one[jj-1]) {
        size_t tmp=vec_one[jj];
        vec_one[jj]=vec_one[jj-1];
        vec_one[jj-1]=tmp;
        size_t tmp1=Vsort[jj];
        Vsort[jj]=Vsort[jj-1];
        Vsort[jj-1]=tmp1;
      }
    }
  }

  return Vsort;
}

void SortVector::sort_vector_ascending(vector<string> &vec_one, size_t nth) {
  if(vec_one.size()<2) return;

  for(int ii=0;ii<(int)vec_one.size()-1;ii++) {
    for(int jj=(int)vec_one.size()-1;jj>ii;jj--) {
      vector<string> stoken=split_string(vec_one[jj]);
      if(nth>=stoken.size()) {
        string str="sort_vector.cpp: Wrong index nth="+to_string(nth)+" !  --> "+vec_one[jj];
        cout<<"Error: "<<str<<endl;
        exit(0);
      }
      double val_j=stod(stoken[nth]);
      stoken=split_string(vec_one[jj-1]);
      if(nth>=stoken.size()) {
        string str="sort_vector.cpp: Wrong index nth="+to_string(nth)+" !  --> "+vec_one[jj-1];
        cout<<"Error: "<<str<<endl;
        exit(0);
      }
      double val_j1=stod(stoken[nth]);
      if(val_j<val_j1) {
        string str=vec_one[jj];
        vec_one[jj]=vec_one[jj-1];
        vec_one[jj-1]=str;
      }
    }
  }

}

vector<size_t> SortVector::sort_vector_descending(vector<double> &vec_one) {
  vector<size_t> Vsort;
  for(size_t ii=0;ii<vec_one.size();ii++)
    Vsort.push_back(ii);
  
  if(vec_one.size()<2) return Vsort;

  for(int ii=0;ii<(int)vec_one.size()-1;ii++) {
    for(int jj=(int)vec_one.size()-1;jj>ii;jj--) {
      if(vec_one[jj]>vec_one[jj-1]) {
        double tmp=vec_one[jj];
        vec_one[jj]=vec_one[jj-1];
        vec_one[jj-1]=tmp;
        size_t tmp1=Vsort[jj];
        Vsort[jj]=Vsort[jj-1];
        Vsort[jj-1]=tmp1;
      }
    }
  }

  return Vsort;
}

vector<size_t> SortVector::sort_vector_descending(vector<int> &vec_one) {
  vector<size_t> Vsort;
  for(size_t ii=0;ii<vec_one.size();ii++)
    Vsort.push_back(ii);
  
  if(vec_one.size()<2) return Vsort;

  for(int ii=0;ii<(int)vec_one.size()-1;ii++) {
    for(int jj=(int)vec_one.size()-1;jj>ii;jj--) {
      if(vec_one[jj]>vec_one[jj-1]) {
        int tmp=vec_one[jj];
        vec_one[jj]=vec_one[jj-1];
        vec_one[jj-1]=tmp;
        size_t tmp1=Vsort[jj];
        Vsort[jj]=Vsort[jj-1];
        Vsort[jj-1]=tmp1;
      }
    }
  }

  return Vsort;
}

vector<size_t> SortVector::sort_vector_descending(vector<size_t> &vec_one) {
  vector<size_t> Vsort;
  for(size_t ii=0;ii<vec_one.size();ii++)
    Vsort.push_back(ii);
  
  if(vec_one.size()<2) return Vsort;

  for(int ii=0;ii<(int)vec_one.size()-1;ii++) {
    for(int jj=(int)vec_one.size()-1;jj>ii;jj--) {
      if(vec_one[jj]>vec_one[jj-1]) {
        size_t tmp=vec_one[jj];
        vec_one[jj]=vec_one[jj-1];
        vec_one[jj-1]=tmp;
        size_t tmp1=Vsort[jj];
        Vsort[jj]=Vsort[jj-1];
        Vsort[jj-1]=tmp1;
      }
    }
  }

  return Vsort;
}

void SortVector::sort_vector_descending(vector<string> &vec_one, size_t nth) {
  if(vec_one.size()<2) return;

  for(int ii=0;ii<(int)vec_one.size()-1;ii++) {
    for(int jj=(int)vec_one.size()-1;jj>ii;jj--) {
      vector<string> stoken=split_string(vec_one[jj]);
      if(nth>=stoken.size()) {
        string str="sort_vector.cpp: Wrong index nth="+to_string(nth)+" !  --> "+vec_one[jj];
        cout<<"Error: "<<str<<endl;
        exit(0);
      }
      double val_j=stod(stoken[nth]);
      stoken=split_string(vec_one[jj-1]);
      if(nth>=stoken.size()) {
        string str="sort_vector.cpp: Wrong index nth="+to_string(nth)+" !  --> "+vec_one[jj-1];
        cout<<"Error: "<<str<<endl;
        exit(0);
      }
      double val_j1=stod(stoken[nth]);
      if(val_j>val_j1) {
        string str=vec_one[jj];
        vec_one[jj]=vec_one[jj-1];
        vec_one[jj-1]=str;
      }
    }
  }
}