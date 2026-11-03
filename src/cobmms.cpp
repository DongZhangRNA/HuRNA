#include "basic.h"
#include "readfile.h"
#include "initialization.h"
#include "error.h"
#include "cobmms.h"

using namespace CBS_NS;

Cobmms::Cobmms(string flog) 
 : error(nullptr), init(nullptr), rand_one(nullptr), logfile(nullptr) 
{
  error = new Error(this);
  init = new Initialization(this);

  rand_one = new Random();
  
  if(flog.find("Cobmms.log")==string::npos && is_file_exist(flog))
    logfile = new ofstream(flog.c_str(),ios::app);
  else
    logfile = new ofstream(flog.c_str());
}

Cobmms::~Cobmms() {
  delete error;
  delete init;

  delete rand_one;

  logfile->close(); logfile->clear();
  
}

void Cobmms::parse_path_to_Cobmms(char *arg) {
  string cmd=arg;
  if(cmd.find("/bin/")==string::npos) {
    string str="Fail to detect Path_to_Cobmms! Please give Path_to_Cobmms using -homedir!";
    error->one(str);
  }
  size_t found=cmd.find_last_of("/");
  string homedir=cmd.substr(0,found);
  found=homedir.find_last_of("/");
  homedir=homedir.substr(0,found+1);

  if(!is_directory_exist(homedir)) {
    string str="Fail to detect Path_to_Cobmms! Please give Path_to_Cobmms using -homedir!";
    error->one(str);
  }

  init->HOME_DIR=homedir;
}

void Cobmms::set_path_to_Cobmms(string homedir) {
  if(!is_directory_exist(homedir)) {
    string str="Cannot find Path_to_Cobmms: "+homedir;
    error->one(str);
  }
  if(homedir.back()!='/')
    homedir+="/";
  init->HOME_DIR=homedir;
}

void Cobmms::initialization() {
  init->set_parameters();
}

void Cobmms::parse_command(map<string,vector<string> > &Mcmd) {
  if(Mcmd.size()==0) return;
  string KEYs[]={"-format","-noh","-np","-nout","-exclude","-circRNA","-workdir","-version","-lmp","-mpi_opt"};
  for(int ii=0;ii<get_array_size(KEYs);ii++) {
    string key=KEYs[ii];
    string keyvalue;
    if(Mcmd.count(key)) {
      vector<string> Vstr=Mcmd[key];
      keyvalue=Vstr[0];
      if(key=="-circRNA") {
        for(size_t kk=1;kk<Vstr.size();kk++)
          keyvalue+=" "+Vstr[kk];
      }
      else if(key=="-mpi_opt") {
        keyvalue=keyvalue.substr(1);
        for(size_t kk=1;kk<Vstr.size();kk++)
          keyvalue+=" "+Vstr[kk];
      }
      else {
        if(Vstr.size()!=1) {
          string str="Unknown command option -->  "+key;
          for(auto ss : Vstr)
            str+=" "+ss;
          cout<<str<<endl;
          exit(0);
        }
      }
      Mcmd.erase(key);
      init->reset_parameters(key,keyvalue);
    }
  }

  if(Mcmd.size()!=0) {
    string str="Unknown command option -->  ";
    for(auto & it :Mcmd) {
      str+=" "+it.first;
      for(auto ss : it.second)
        str+=" "+ss;
    }
    cout<<str<<endl;
    exit(0);
  }
}

void Cobmms::execute_system_command(string cmd) {
  string str="Command: "+cmd;
  cout<<str<<endl;
  error->log(str);
  
  int sys=system(cmd.c_str());

  if(sys!=0) {
    string str="Fail to execute "+cmd;
    error->one(str);
  }
}

void Cobmms::run_lammps(string cmd) {
  string str="Command: "+cmd;
  cout<<str<<endl;
  error->log(str);
  
  int sys=system(cmd.c_str());

  if(sys!=0) {
    string str="Fail to execute "+cmd;
    error->warning(str);
    str="Please reduce the timestep in lammps configure file!";
    error->warning(str);
  }
}

vector<string> Cobmms::execute_command_to_string(string cmd) {
  string fout=init->WORK_DIR+"out.dat";
  cmd+=" > "+fout;

  int sys=system(cmd.c_str());
  if(sys!=0) {
    string str="Fail to execute "+cmd;
    error->one(str);
  }

  vector<string> Vlines=read_file(fout);

  string cmd1="rm "+fout;
  sys=system(cmd1.c_str());
  if(sys!=0) {
    string str="Fail to execute "+cmd1;
    cout<<str<<endl;
  }
  
  return Vlines;
}

string Cobmms::calculate_runtime(double dtime) {
  int dH=dtime/3600;
  int dM=(dtime-dH*3600)/60;
  int dS=dtime-dH*3600-dM*60;
  string stime="Runtime: "+to_string(dH)+" Hours "+to_string(dM)+" Minutes "+to_string(dS)+" Seconds.";

  return stime;
}