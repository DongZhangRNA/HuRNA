#include "basic.h"
#include "error.h"

using namespace CBS_NS;

Error::Error(Cobmms *ptr) : Pointers(ptr) {}

Error::~Error() {
  //outERR.close();
  //outERR.clear();
}

/*
void Error::output_error(string fstr) {
  outERR.open(fstr.c_str());
  if(!outERR.is_open()) {
    cout<<"Can not open file "<<fstr<<" !"<<endl;
    exit(0);
  }
}
*/

void Error::one(string str) {
  cout<<"Error: "<<str<<endl;
  //outERR<<"Error: "<<str<<endl;
  (*logfile)<<"Error: "<<str<<endl;
  exit(0);
}

void Error::warning(string str) {
  //cout<<"WARNING: "<<str<<endl;
  //outERR<<"Warning: "<<str<<endl;
  (*logfile)<<"Warning: "<<str<<endl;
}

void Error::message(string str) {
  cout<<"Messsage: "<<str<<endl;
}

void Error::log(string str) {
  (*logfile)<<str<<endl;
}