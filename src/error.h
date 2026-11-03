#ifndef CBS_ERROR_H
#define CBS_ERROR_H

#include "pointers.h"

namespace CBS_NS {

class Error : protected Pointers {
 public:
  Error(Cobmms *);

  ~Error();

  //ofstream outERR;
  //void output_error(string);
  void one(string);
  void warning(string);
  void message(string);
  void log(string);
};

}

#endif