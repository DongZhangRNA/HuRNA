#ifndef CBS_POINTERS_H
#define CBS_POINTERS_H

#include "cobmms.h"

namespace CBS_NS {

class Pointers {   //share log and initialization
 public:
  Pointers(Cobmms *ptr) :
    cbs(ptr),
    error(ptr->error),
    init(ptr->init),

    rand_one(ptr->rand_one),
    
    logfile(ptr->logfile) {}
  Pointers& operator=(const Pointers &other) {
    if(this != &other) {
      cbs = other.cbs;
      error = other.cbs->error;
      init = other.cbs->init;
      
      rand_one = other.cbs->rand_one;
      
      logfile = other.cbs->logfile;
    }

    return *this;
  }
  virtual ~Pointers() {};

 protected:
  Cobmms *cbs;
  Error *&error;
  Initialization *&init;

  Random *&rand_one;

  ofstream *&logfile;
};

}

#endif