#ifndef CBS_COBMMS_H
#define CBS_COBMMS_H

#include "random.h"

namespace CBS_NS {

class Cobmms {
 public:
  class Error *error;
  class Initialization *init;
  
  Random *rand_one;

  ofstream *logfile;

  Cobmms(string);
  ~Cobmms();

  void parse_path_to_Cobmms(char *);
  void set_path_to_Cobmms(string);
  void initialization();
  void parse_command(map<string,vector<string> > &);
  void execute_system_command(string);
  void run_lammps(string);
  vector<string> execute_command_to_string(string);
  string calculate_runtime(double);

 private:
  Cobmms() {};  // prohibit using the default constructor
  Cobmms(const Cobmms &) {};  // prohibit using the copy constructor
};

}

#endif