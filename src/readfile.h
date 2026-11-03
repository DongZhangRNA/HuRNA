#ifndef _READFILE_H
#define _READFILE_H

inline vector<string> read_file(string fstr) {
  vector<string> Vstr;
  string sline;
  ifstream inF (fstr.c_str());
  if(!inF.is_open()) {
    cout<<"Failed to read file: "<<fstr<<endl;
    exit(0);
  }
  while(getline(inF,sline)) {
    if(sline.size()<1) continue;
    Vstr.push_back(sline);
  }
  inF.close(); inF.clear();

  return Vstr;
}

inline vector<string> get_dir_files (string dir) {
  vector<string> files;
  DIR *dp, *dp_tmp;
  struct dirent *dirp;
  if((dp  = opendir(dir.c_str())) == NULL) {
    cout << "Error(" << errno << ") opening " << dir << endl; //return errno;
  }

  while ((dirp = readdir(dp)) != NULL) {
    string dir_tmp = string(dirp->d_name);
    if(dir_tmp != "." && dir_tmp != "..") {
      if((dp_tmp  = opendir(dir_tmp.c_str())) == NULL) files.push_back(dir_tmp);
      else closedir(dp_tmp);
    }
  }
  closedir(dp);

  return files;
}

inline bool is_file_exist(const std::string& filename) {
  struct stat buf;
  if (stat(filename.c_str(), &buf) != -1)
    return true;
  else
    return false;
}

inline bool is_directory_exist(const std::string& directoryPath) {
  struct stat info;
  if (stat(directoryPath.c_str(), &info) != 0)
    return false;
  else if (info.st_mode & S_IFDIR)
    return true;
  else
    return false;
}

inline string trim(const string& str) {
  const char* ws = " \t\r\n\v\f";
  size_t start = str.find_first_not_of(ws);
  if (start == std::string::npos) return "";
  size_t end = str.find_last_not_of(ws);
  return str.substr(start, end - start + 1);
}

inline bool line_starts_with(const string& str, const string& prefix) {
  return str.size() >= prefix.size() && str.compare(0, prefix.size(), prefix) == 0;
}

inline bool line_ends_with(const string& str, const string& suffix) {
  return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

#endif