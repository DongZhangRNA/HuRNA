#include "basic.h"
#include "error.h"
#include "readfile.h"
#include "initialization.h"
#include "molecule.h"
#include "read_PDB.h"
#include "read_lammps.h"
#include "cluster.h"

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

  string flist_trj="";
  if(Mcmd.count("-lt")) {
    flist_trj=Mcmd["-lt"][0];
    Mcmd.erase("-lt");
  }
  string fintrj="";
  if(Mcmd.count("-it")) {
    fintrj=Mcmd["-it"][0];
    Mcmd.erase("-it");
  }
  if(flist_trj.length()==0 && fintrj.length()==0) {
    cout<<"Wrong commands: No input trj files (-it or -lt)!"<<endl;
    exit(0);
  }
  string flist_log="";
  if(Mcmd.count("-ll")) {
    flist_log=Mcmd["-ll"][0];
    Mcmd.erase("-ll");
  }
  string finlog="";
  if(Mcmd.count("-il")) {
    finlog=Mcmd["-il"][0];
    Mcmd.erase("-il");
  }
  if(flist_log.length()==0 && finlog.length()==0) {
    cout<<"Wrong commands: No input log files (-il or -ll)!"<<endl;
    exit(0);
  }
  string inDIR="";
  if(Mcmd.count("-inDIR")) {
    inDIR=Mcmd["-inDIR"][0];
    Mcmd.erase("-inDIR");
  }
  if(inDIR.size()!=0 && inDIR.back()!='/')
    inDIR+="/";
  string fout="";
  if(Mcmd.count("-o")) {
    fout=Mcmd["-o"][0];
    Mcmd.erase("-o");
  }
  else {
    cout<<"Wrong commands: No output (-o)!"<<endl;
    exit(0);
  }
  if(fout.find(".pdb")==string::npos) {
    fout+=".pdb";
  }
  string fref;
  if(Mcmd.count("-ref")) {
    fref=Mcmd["-ref"][0];
    Mcmd.erase("-ref");
  }
  else {
    cout<<"Wrong commands: No input for residue info (-ref)!"<<endl;
    exit(0);
  }
  int step_begin=0;
  int flag_check=0;
  if(Mcmd.count("-b")) {
    flag_check=1;
    step_begin=stoi(Mcmd["-b"][0]);
    Mcmd.erase("-b");
  }
  int step_end=-1;
  if(Mcmd.count("-e")) {
    flag_check=1;
    step_end=stoi(Mcmd["-e"][0]);
    Mcmd.erase("-e");
  }
  double Pstru=0.5; // last percentage of snapshots as structure ensemble
  if(flag_check) Pstru=1.0;
  if(Mcmd.count("-ps")) {
    Pstru=stod(Mcmd["-ps"][0]);
    Mcmd.erase("-ps");
  }
  double Peng=0.1; // top potential energy structrues
  if(Mcmd.count("-pe")) {
    Peng=stod(Mcmd["-pe"][0]);
    Mcmd.erase("-pe");
  }
  double rmsdcutoff=-1.0;
  if(Mcmd.count("-cutoff")) {
    rmsdcutoff=stod(Mcmd["-cutoff"][0]);
    Mcmd.erase("-cutoff");
  }
  double Pclu=0.5; // percentage of structures in largest cluster
  if(Mcmd.count("-pc")) {
    Pclu=stod(Mcmd["-pc"][0]);
    Mcmd.erase("-pc");
  } 
  double Pclu_min=0.15; // min percentage of structures in largest cluster
  if(Mcmd.count("-pcmin")) {
    Pclu_min=stod(Mcmd["-pcmin"][0]);
    Mcmd.erase("-pcmin");
  }
  string frmsdmatrix="";
  int flag_matrix=0;
  if(Mcmd.count("-mat")) {
    flag_matrix=1;
    frmsdmatrix=Mcmd["-mat"][0];
    Mcmd.erase("-mat");
  }
  if(frmsdmatrix.length()==0) {
    size_t found=fout.find(".pdb");
    frmsdmatrix=fout.substr(0,found)+"_rmsd_matrix.dat";
  }
  
  string flog="./Cobmms.log";
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
    Mcmd.erase("-homedir");
  }
  cbs->initialization();
  cbs->parse_command(Mcmd);

  time_t tstart,tend;
  time (&tstart);

  string slog="##### Obtain_top_structure #####";
  cbs->error->log(slog);

  vector<string> Vtrjs;
  if(fintrj.length()>0) Vtrjs.push_back(fintrj);
  if(flist_trj.length()>0) {
    vector<string> Vlines=read_file(flist_trj);
    for(size_t ii=0;ii<Vlines.size();ii++) {
      if(Vlines[ii].size()<1) continue;
      vector<string> stoken=split_string(Vlines[ii]);
      string fname=stoken[0];
      if(inDIR.size()>0)
        fname=inDIR+fname;
      if(is_file_exist(fname))
        Vtrjs.push_back(fname);
      else {
        string str="Cannot find file "+fname+" !";
        cbs->error->warning(str);
      }
    }
  }
  vector<string> Vlogs;
  if(finlog.length()>0) Vlogs.push_back(finlog);
  if(flist_log.length()>0) {
    vector<string> Vlines=read_file(flist_log);
    for(size_t ii=0;ii<Vlines.size();ii++) {
      if(Vlines[ii].size()<1) continue;
      vector<string> stoken=split_string(Vlines[ii]);
      string fname=stoken[0];
      if(inDIR.size()>0)
        fname=inDIR+fname;
      if(is_file_exist(fname))
        Vlogs.push_back(fname);
      else {
        string str="Cannot find file "+fname+" !";
        cbs->error->warning(str);
      }
    }
  }

  ReadPDB readpdb_one(cbs);
  readpdb_one.read_PDB(fref);
  vector<Atom> Vatoms=readpdb_one.Full_atoms[0];

  Molecule mol_ref(cbs);
  mol_ref.set_molecule(Vatoms,Chain::cg);
  mol_ref.set_atomid(0);
  size_t len=mol_ref.number_of_residues();

  ReadLAMMPS readlmp_one(cbs);
  readlmp_one.read_lammpstrj(Vtrjs,step_begin,step_end);
  readlmp_one.read_log(Vlogs,step_begin,step_end);
  vector<vector<Atom> > Vstrus=readlmp_one.extract_top_structures(Peng,Pstru,1);
  
  Cluster clu_one(cbs);
  clu_one.input_data(Vstrus);
  if(flag_matrix==0 && is_file_exist(frmsdmatrix)) {
    string cmd="rm "+frmsdmatrix;
    cbs->execute_system_command(cmd);
  }
  clu_one.calculate_RMSD_matrix(frmsdmatrix);
  if(rmsdcutoff>0.0) {
    clu_one.get_clusters(rmsdcutoff);
  }
  else {
    rmsdcutoff=0.1*len;
    clu_one.get_clusters(rmsdcutoff,Pclu,Pclu_min);
  }
  
  (*(cbs->logfile))<<clu_one;

  int nout=cbs->init->Num_model_max;
  if((int)clu_one.size()<nout)
    nout=static_cast<int>(clu_one.size());
  
  ofstream outR (fout.c_str());
  for(int ii=0;ii<nout;ii++) {
    size_t index_one=clu_one.clusters[ii][0];
    vector<Atom> Vstru_center=Vstrus[index_one];
    vector<double> Vxyz;
    for(size_t jj=0;jj<Vstru_center.size();jj++) {
      vector<double> xyz=Vstru_center[jj].get_atom_xyz();
      Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
    }
    mol_ref.reset_molecule_xyz(Vxyz);
    vector<double> Vxyz_cog=mol_ref.center_of_geometry();
    mol_ref.move_molecule(-Vxyz_cog[0],-Vxyz_cog[1],-Vxyz_cog[2]);
    outR<<"MODEL     "<<ii+1<<endl;
    outR<<mol_ref;
    outR<<"ENDMDL"<<endl;
  }
  outR<<"END"<<endl;
  outR.close(); outR.clear();  

  time (&tend);

  double dtime=difftime(tend,tstart);
  string stime=cbs->calculate_runtime(dtime);
  cbs->error->log(stime);
  stime="";
  cbs->error->log(stime);

  delete cbs;
  
  return 0;
}