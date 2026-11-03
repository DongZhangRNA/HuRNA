#include "basic.h"
#include "readfile.h"
#include "error.h"
#include "initialization.h"

using namespace CBS_NS;

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

  string fin="";
  string seq="",secstr="";
  if(Mcmd.count("-i")) {
    fin=Mcmd["-i"][0];
    Mcmd.erase("-i");
  }
  vector<string> Vlines=read_file(fin);
  int ncount=0;
  for(size_t ii=0;ii<Vlines.size();ii++) {
    string sline=Vlines[ii];
    if(sline.find("#")!=string::npos) continue;
    ncount++;
    if(ncount==1) seq=sline;
    else if(ncount==2) secstr=sline;
    else {
      cout<<"Wrong input file --> "<<fin<<endl;
      exit(0);
    }
  }
  if(seq.size()==0) {
    cout<<"No input sequence (-i) !"<<endl;
    exit(0);
  }
  
  string jobname="rna";
  if(Mcmd.count("-name")) {
    jobname=Mcmd["-name"][0];
    Mcmd.erase("-name");
  }

  string DIR="./"+jobname+"_prediction";
  if(Mcmd.count("-DIR")) {
    DIR=Mcmd["-DIR"][0];
    Mcmd.erase("-DIR");
  }
  
  string fexclude="";
  if(Mcmd.count("-exclude")) {
    fexclude=Mcmd["-exclude"][0];
    Mcmd.erase("-exclude");
  }
  
  string frestr="";
  if(Mcmd.count("-restr")) {
    frestr=Mcmd["-restr"][0];
    Mcmd.erase("-restr");
  }
 
  int Nstep=-1;
  if(Mcmd.count("-nstep")) {
    Nstep=stoi(Mcmd["-nstep"][0]);
    Mcmd.erase("-nstep");
  }

  string _version="v2";
  if(Mcmd.count("-version")) {
    _version=Mcmd["-version"][0];
    Mcmd.erase("-version");
  }

  int flag_gpu=0;
  if(Mcmd.count("-gpu")) {
    flag_gpu=1;
    Mcmd.erase("-gpu");
  }

  int flag_RNA_use_DM=0;
  if(Mcmd.count("-use_DM")) {
    flag_RNA_use_DM=1;
    Mcmd.erase("-use_DM");
  }

  int flag_pe=0;
  string pe="";
  if(Mcmd.count("-pe")) {
    flag_pe=1;
    pe=Mcmd["-pe"][0];
    Mcmd.erase("-pe");
  }
  int flag_cutoff=0;
  string cutoff="";
  if(Mcmd.count("-cutoff")) {
    flag_cutoff=1;
    cutoff=Mcmd["-cutoff"][0];
    Mcmd.erase("-cutoff");
  }
  int flag_pc=0;
  string pc="";
  if(Mcmd.count("-pc")) {
    flag_pc=1;
    pc=Mcmd["-pc"][0];
    Mcmd.erase("-pc");
  }

  string secstr_pred_method="contrafold";
  if(Mcmd.count("-2d_pred")) {
    secstr_pred_method = Mcmd["-2d_pred"][0];
    Mcmd.erase("-2d_pred");
  }

  int flag_exe=1;
  if(Mcmd.count("-exe")) {
    flag_exe=stoi(Mcmd["-exe"][0]);
    Mcmd.erase("-exe");
  }

  int flag_template_rigid=0;
  if(Mcmd.count("-template_rigid")) {
    flag_template_rigid=1;
    Mcmd.erase("-template_rigid");
  }

  string flog="./"+jobname+".log";
  if(Mcmd.count("-log")) {
    vector<string> Vlog=Mcmd["-log"];
    if(Vlog.size()!=1) {
      string str="Unknown command option -->  -log";
      for(auto ss : Vlog)
        str+=" "+ss;
      cout<<str<<endl;
      exit(0);
    }
    flog=Vlog[0];
    Mcmd.erase("-log");
  } 
  if(flog.find(".log")==string::npos)
    flog+=".log";

  Cobmms *cbs= new Cobmms(flog);
  if(Mcmd.count("-homedir")==0) {
    cbs->parse_path_to_Cobmms(arg[0]);
  }
  else {
    cbs->set_path_to_Cobmms(Mcmd["-homedir"][0]);
  }
  cbs->initialization();

  string BIN=cbs->init->HOME_DIR+"bin";
  if(Mcmd.count("-bin")) {
    BIN=Mcmd["-bin"][0];
    Mcmd.erase("-bin");
  }

  string cmd_add="";
  for(auto & it : Mcmd) {
    cmd_add+=" "+it.first;
    for(auto str : it.second) {
      if(str[0]=='-' && str[1]=='-') 
        str="\""+str+"\"";
      cmd_add+=" "+str;
    }
  }

  if(secstr.size()==0) {
    string slog="No secondary structure input, predicted by "+secstr_pred_method+"!";
    cbs->error->log(slog);
    if(!is_directory_exist(DIR)) {
      string cmd="mkdir "+DIR;
      cbs->execute_system_command(cmd);
    }
    vector<string> stoken=split_string(seq);
    seq="";
    for(size_t jj=0;jj<stoken.size();jj++)
      seq+=stoken[jj];
    if(secstr_pred_method=="contrafold") {
      string contrafold=cbs->init->Contrafold;
      if(!is_file_exist(contrafold)) {
        string str="Cannot find contrafold for 2D structure prediction! --> "+contrafold;
        cbs->error->one(str);
      }
      string fseq=DIR+"/"+jobname+".seq";
      ofstream outS (fseq.c_str());
      outS<<seq<<endl;
      outS.close(); outS.clear();
      string fpred_2d=DIR+"/"+jobname+"_contrafold.2d";
      string contra_cmd=contrafold+" predict "+fseq+" | sed -n '2,3p' > "+fpred_2d;
      cbs->execute_system_command(contra_cmd);
      fin=fpred_2d;
    }
    else if(secstr_pred_method=="mxfold2") {
      string mxfold2_cmd="conda run -n mxfold2_cobmms mxfold2 predict ";
      string fseq=DIR+"/"+jobname+".seq";
      ofstream outS (fseq.c_str());
      outS<<">"<<jobname<<endl;
      outS<<seq<<endl;
      outS.close(); outS.clear();
      string fpred_2d=DIR+"/"+jobname+"_mxfold2.2d";
      string mxfold2_run=mxfold2_cmd+fseq+" | awk 'NR==2; NR==3 {print $1}' > "+fpred_2d;
      cbs->execute_system_command(mxfold2_run);
      fin=fpred_2d;
    } 
    else {
      string str="Only 'contrafold' or 'mxfold2' was supported for RNA 2D structure prediction!";
      cbs->error->one(str);
    }
  }

  string fstr="./predict_"+jobname+".sh";
  ofstream outF (fstr.c_str());
  outF<<"#!/bin/bash"<<endl;
  outF<<"BIN="<<BIN<<endl;
  outF<<"jobname="<<jobname<<endl;
  outF<<"DIR="<<DIR<<endl;
  string sline="${BIN}/Construct_RNA_initial_structure -i "+fin+" -log ${jobname}.log ";
  sline+="-name ${jobname} -DIR ${DIR} -version "+_version;
  if(fexclude.size()>0)
    sline+=" -exclude "+fexclude;
  if(frestr.size()>0)
    sline+=" -restr "+frestr;
  if(flag_gpu)
    sline+=" -gpu";
  if(flag_RNA_use_DM)
    sline+=" -use_DM";
  if(cmd_add.size()>0)
    sline+=cmd_add;
  outF<<sline<<endl;
  sline="${BIN}/Perform_CGMD_simulation -i "+fin+" -if ${DIR}/${jobname}_build_cg_init.pdb";
  sline+=" -iformat cg -restr ${DIR}/${jobname}_restraint.dat -DIR ${DIR} -name ${jobname} ";
  sline+=" -log ${jobname}.log -ti_fold -ss_restr -version "+_version;
  if(Nstep>0)
    sline+=" -nstep "+to_string(Nstep);
  if(flag_gpu)
    sline+=" -gpu";
  if(flag_template_rigid)
    sline+=" -template_rigid";
  if(cmd_add.size()>0)
    sline+=cmd_add;
  outF<<sline<<endl;
  sline="${BIN}/Obtain_top_structure -lt ${DIR}/${jobname}_prod_lammpstrj_list.dat -ll ";
  sline+="${DIR}/${jobname}_prod_lammpslog_list.dat -ref ${DIR}/${jobname}_cg.pdb -log ${jobname}.log";
  sline+=" -o ${DIR}/${jobname}_predicted_cg.pdb -b 0";
  if(flag_pe)
    sline+=" -pe "+pe;
  if(flag_cutoff)
    sline+=" -cutoff "+cutoff;
  if(flag_pc)
    sline+=" -pc "+pc;
  if(cmd_add.size()>0)
    sline+=cmd_add;
  outF<<sline<<endl;
  sline="${BIN}/Reconstruct_allatom_structure -i ${DIR}/${jobname}_predicted_cg.pdb -ss "+fin+" -DIR ${DIR} ";
  sline+="-name ${jobname} -o ${DIR}/${jobname}_predicted.pdb -log ${jobname}.log -min -version "+_version;
  if(cmd_add.size()>0)
    sline+=cmd_add;
  outF<<sline<<endl;
  outF<<"if [ -f ${DIR}/${jobname}_predicted_min.pdb ]; then"<<endl;
  outF<<"  cp ${DIR}/${jobname}_predicted_min.pdb ."<<endl;
  outF<<"elif [ -f ${DIR}/${jobname}_predicted.pdb ]; then"<<endl;
  outF<<"  cp ${DIR}/${jobname}_predicted.pdb ."<<endl;
  outF<<"fi"<<endl;
  outF.close(); outF.clear();

  if(flag_exe) {
    string cmd="bash "+fstr;
    cbs->execute_system_command(cmd);
  }

  delete cbs;

  return 0;
}