#ifndef _BASIC_H_
#define _BASIC_H_

#include <cstddef>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <regex>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <ctype.h>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <iterator>
#include <time.h>
#include <assert.h>
#include <limits>
#include <thread>
#include <memory>
using namespace std;

inline string to_string (const char aa) {
  string bb;
  stringstream ss;
  ss << aa;
  ss >> bb;

  return bb;
}

inline string to_string(const double x, const int len) {
  string str;
  ostringstream oss;
  oss<<std::fixed<<std::setprecision(len)<<x;
  str=oss.str();
  oss.str("");
  oss.clear();

  return str;
}

inline vector<string> split_string(string sline, string sp=" ") {
  size_t found=sline.find_first_of(sp);
  while(found!=string::npos) {
    sline[found]=' ';
    found=sline.find_first_of(sp,found+1);
  }

  istringstream ss (sline);
  string buf;
  vector<string> token;
  while(ss >> buf) token.push_back(buf);

  return token;
}

template <class T>
inline int get_array_size(const T &arr) {
  return (sizeof(arr)/sizeof(arr[0]));
}

inline void pause() {
  using std::numeric_limits;
  using std::streamsize;
  cout <<"Pause: Press ENTER to continue..." << flush;
  cin.clear();  // use only if a human is involved
  cin.ignore(numeric_limits<streamsize>::max(), '\n');
}

/*
inline bool is_valid_integer(const string& str) {
  if(str.empty()) {
    return false;
  }
    
  try {
    size_t pos=0;
    int value=stoi(str, &pos);
    return pos == str.length();
  }
  catch (const std::invalid_argument& e) {
    return false;
  }
  catch (const std::out_of_range& e) {
    return false;
  }
}
*/

#endif
