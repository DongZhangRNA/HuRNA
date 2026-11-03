#include "basic.h"
#include "readfile.h"

int main(int narg, char **arg) {
  map<string,vector<string> > Mcmd;
  string style="";
  vector<string> Vargs;
  for(int ii=1;ii<narg;ii++) {
    if(arg[ii][0]=='-') {
      if(style!="") {
        if(Mcmd.count(style)) {
          cout<<"Wrong commands!"<<endl;
          exit(0);
        }
        Mcmd.emplace(style,Vargs);
      }
      style=arg[ii];
      Vargs.clear();
      continue;
    }
    Vargs.push_back(arg[ii]);
  }
  if(Mcmd.count(style)) {
    cout<<"Wrong commands!"<<endl;
    exit(0);
  }
  Mcmd.emplace(style,Vargs);

  if(Mcmd.count("-i")==0 && Mcmd.count("-list")==0) {
    cout<<"Wrong commands: No input (PDBid or PDBlist) !"<<endl;
    exit(0);
  }
  if(Mcmd.count("-o")==0) {
    cout<<"Wrong commands: No output directory!"<<endl;
    exit(0);
  }
  
  string outDIR=Mcmd["-o"][0];
  if(outDIR.back()!='/') 
    outDIR+="/";
  string type="cif";
  if(Mcmd.count("-type")) 
    type=Mcmd["-type"][0];
  if(type!="cif" && type!="pdb") {
    cout<<"Error: Cannot fecth non pdb/cif file --> ."+type;
    return 1;
  }
  
  string Tstr="./FP_out.dat";
  string cmd="which curl > "+Tstr;
  int sys=system(cmd.c_str());

  if(sys!=0) {
    string str="Error: Fail to execute "+cmd;
    cout<<str<<endl;
    return 1;
  }

  vector<string> lines=read_file(Tstr.c_str());
  cmd="rm "+Tstr;
  sys=system(cmd.c_str());

  if(sys!=0) {
    string str="Error: Fail to execute "+cmd;
    cout<<str<<endl;
    return 1;
  }

  if(lines.size()==0) {
    cout<<"You need to install 'curl' for this program to work"<<endl;
    return 1;
  }


  vector<string> Vpdbs;
  if(Mcmd.count("-i")==1) {
    for(size_t ii=0;ii<Mcmd["-i"].size();ii++)
      Vpdbs.push_back(Mcmd["-i"][ii]);
  }
  if(Mcmd.count("-list")==1) {
    string fstr=Mcmd["-list"][0];
    vector<string> Vlines=read_file(fstr);
    for(size_t ii=0;ii<Vlines.size();ii++) {
      vector<string> stoken=split_string(Vlines[ii],",");
      for(size_t jj=0;jj<stoken.size();jj++)
        Vpdbs.push_back(stoken[jj]);
    }
  }
  vector<string> PDBs;
  for(size_t ii=0;ii<Vpdbs.size();ii++) {
    string PDBid=Vpdbs[ii];
    int check=0;
    for(size_t jj=0;jj<PDBid.size();jj++) {
      char ch=PDBid[jj];
      if(ch>='0' && ch<='9') continue;
      else if(ch>='a' && ch<='z') continue;
      else if (ch>='A' && ch<='Z') 
        PDBid[jj]=tolower(PDBid[jj]);
      else 
        check=1;
    }
     
    if(PDBid.size()!=4 && check) {
      cout<<"Skip invalid PDBid --> "<<PDBid<<endl;
    }
    PDBs.push_back(PDBid);
  }

  int N_exi=0;
  int N_dow=0;
  for(size_t ii=0;ii<PDBs.size();ii++) {
    if(ii>0 && ii%100==0) {
      cout<<"Downloaded "<<ii<<"/"<<PDBs.size()<<" PDBs."<<endl;
    }
    string PDBid=PDBs[ii]+"."+type;
    string fstr=outDIR+PDBid;
    if(is_file_exist(fstr)) {
      N_exi++;
      continue;
    }
    cmd="curl -s -f https://files.rcsb.org/download/"+PDBid+" -o "+fstr;
    cmd+=" || echo "+PDBid+ " >> "+Tstr;
    sys=system(cmd.c_str());

    if(sys!=0) {
      string str="Error: Fail to execute "+cmd;
      cout<<str<<endl;
      return 1;
    }
  }
  
  lines.clear();
  int N_fail=0;
  if(is_file_exist(Tstr)) {
    lines=read_file(Tstr);
    cmd="rm "+Tstr;
    sys=system(cmd.c_str());

    if(sys!=0) {
      string str="Error: Fail to execute "+cmd;
      cout<<str<<endl;
      return 1;
    }

    N_fail=(int)lines.size();
  }
  
  N_dow=(int)PDBs.size()-N_exi-N_fail;
  cout<<"Total "<<PDBs.size()<<" PDBs: "<<N_dow<<" downloaded, "<<N_exi<<" already exist, ";
  cout<<N_fail<<" download failed."<<endl;

  if(lines.size()>0) {
    cout<<endl<<"####################"<<endl;
    cout<<"Fail to download PDBs: "<<endl;
    for(size_t ii=0;ii<lines.size();ii++) {
      cout<<lines[ii]<<endl;
      string fstr=outDIR+lines[ii];
      if(is_file_exist(fstr)) {
        cmd="rm "+fstr;
        sys=system(cmd.c_str());

        if(sys!=0) {
          string str="Error: Fail to execute "+cmd;
          cout<<str<<endl;
          return 1;
        }
      }
    }
  }
  
  return 0;
}