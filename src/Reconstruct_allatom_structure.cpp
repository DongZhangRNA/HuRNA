#include "basic.h"
#include "error.h"
#include "readfile.h"
#include "initialization.h"
#include "molecule.h"
#include "read_PDB.h"
#include "RNA_coarse_grained.h"
#include "RNA_recovery.h"
#include "RNA_molecule.h"
#include "RNA_basepair.h"
#include "read_amber.h"
#include "read_coeff.h"
#include "prepare_lammps.h"
#include "read_lammps.h"

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

  string fin;
  if(Mcmd.count("-i")) {
    fin=Mcmd["-i"][0];
    Mcmd.erase("-i");
  }
  else {
    cout<<"Wrong commands: No input (-i)!"<<endl;
    exit(0);
  }
  string DIR="./prediction/";
  if(Mcmd.count("-DIR")) {
    DIR=Mcmd["-DIR"][0];
    Mcmd.erase("-DIR");
  }
  if(DIR.back()!='/')
    DIR+="/";
  string fout="";
  if(Mcmd.count("-o")) {
    fout=Mcmd["-o"][0];
    Mcmd.erase("-o");
  }
  else {
    cout<<"Wrong commands: No output (-o)!"<<endl;
    exit(0);
  }
  size_t pos_fmt=max(0,(int)fout.size()-4);
  string fout_fmt=fout.substr(pos_fmt);
  if(fout_fmt!=".pdb") 
    fout+=".pdb";

  string fss="";
  if(Mcmd.count("-ss")) {
    fss=Mcmd["-ss"][0];
    Mcmd.erase("-ss");
  }

  int nmodel=-1;
  if(Mcmd.count("-nmodel")) {
    nmodel=stoi(Mcmd["-nmodel"][0]);
    Mcmd.erase("-nmodel");
  }
  
  int flag_osep=0; //output separated model 
  if(Mcmd.count("-os")) {
    flag_osep=1;
    Mcmd.erase("-os");
  }
  int flag_min=0;  //energy minimization for allatom structure 
  if(Mcmd.count("-min")) {
    flag_min=1;
    Mcmd.erase("-min");
  }
  string jobname="min";
  if(Mcmd.count("-name")) {
    jobname=Mcmd["-name"][0];
    Mcmd.erase("-name");
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
    Mcmd.erase("-homedir");
  }
  cbs->initialization();
  cbs->parse_command(Mcmd);

  time_t tstart,tend;
  time (&tstart);

  string slog="##### Reconstruct_allatom_structure #####";
  cbs->error->log(slog);

  string cmd;
  if(!is_directory_exist(DIR)) {
    cmd="mkdir "+DIR;
    cbs->execute_system_command(cmd);
  }

  cbs->init->WORK_DIR=DIR;

  ReadPDB readpdb_one(cbs);
  readpdb_one.read_PDB(fin);
  if(readpdb_one.RNA_atoms.size()==0 || readpdb_one.RNA_atoms[0].size()==0) {
    cout<<"No RNA atoms in "<<fin<<endl;
    exit(0);
  }
  vector<Molecule> Vmol_cg;
  if(nmodel<0) 
    nmodel=static_cast<int>(readpdb_one.RNA_atoms.size());
  if((int)readpdb_one.RNA_atoms.size()<nmodel)
    nmodel=static_cast<int>(readpdb_one.RNA_atoms.size());
  for(int ii=0;ii<nmodel;ii++) {
    Molecule mol_one(cbs);
    mol_one.set_molecule(readpdb_one.RNA_atoms[ii],Chain::cg);
    mol_one.set_atomid(0);
    Vmol_cg.push_back(mol_one);
  }

  map<string,string> mapBP_one;
  if(fss.size()>0) {
    vector<string> Vlines=read_file(fss);
    int ncount=0;
    string seq,secstr;
    for(size_t ii=0;ii<Vlines.size();ii++) {
      if(Vlines[ii].find("#")!=string::npos) continue;
      vector<string> stoken=split_string(Vlines[ii]);
      string sline=stoken[0];
      for(size_t jj=1;jj<stoken.size();jj++)
        sline+=" "+stoken[jj];
      if(sline.size()<1) continue;
      ncount++;
      if(ncount==1) seq=sline;
      else if(ncount==2) secstr=sline;
    }
    if(ncount!=2) {
      string str="Wrong input file for sequence & secondary structure --> "+fss;
      cbs->error->warning(str);
    }
    else {
      RNAMolecule rna_one(cbs);
      rna_one.set_RNA_molecule(seq);
      vector<vector<string> > Vres=rna_one.get_residue_name();
      RNABasePair rna_bp_one(cbs);
      rna_bp_one.parse_RNA_basepairs(Vres,secstr);
      mapBP_one=rna_bp_one.mapBP;
    }
  }

  ReadCoeff readcoeff_one(cbs);
  readcoeff_one.read_RNA_coeff();
  readcoeff_one.get_interaction_style();
  PrepareLAMMPS prelmp_one(cbs);
  prelmp_one.set_coeff(readcoeff_one);
  string fname=jobname+"_min_cg";
  for(size_t ii=0;ii<Vmol_cg.size();ii++) {
    RNACoarseGrained rnacg_one(cbs);
    rnacg_one.reset_RNA_CG(Vmol_cg[ii]);
    if(cbs->init->RNA_version=="latest")
      rnacg_one.build_RNA_topology();
    else if(cbs->init->RNA_version=="v2")
      rnacg_one.build_RNA_topology_v2();
    if(mapBP_one.size()>0) {
      if(cbs->init->RNA_version=="latest")
        rnacg_one.add_RNA_basepairs(mapBP_one);
      else if(cbs->init->RNA_version=="v2")
        rnacg_one.add_RNA_basepairs_v2(mapBP_one);
    }
    if(cbs->init->CircRNA_chains.size()>0) {
      if(cbs->init->RNA_version=="latest")
        rnacg_one.add_circular_RNA_topology();
      else if(cbs->init->RNA_version=="v2")
        rnacg_one.add_circular_RNA_topology_v2();
    }
    vector<string> Vdata=prelmp_one.prepare_data(rnacg_one.rna_cg,rnacg_one.bonds,rnacg_one.angles,rnacg_one.dihedrals,rnacg_one.impropers);
    string fstr=DIR+"data."+fname+"."+to_string(ii+1);
    ofstream outF (fstr.c_str());
    for(auto sline : Vdata)
      outF<<sline<<endl;
    outF.close(); outF.clear();
  }
  prelmp_one.md_label=PrepareLAMMPS::em;
  vector<string> Vrun_files=prelmp_one.write_in_script(DIR,fname,Vmol_cg.size());
  vector<string> Vtrjs_list=prelmp_one.Voutput_trjs;
  string fsh=cbs->init->WORK_DIR+"run_lammps_"+fname+".sh";
  vector<string> Vcmd=prelmp_one.write_run_script(Vrun_files,fsh);
  string flammpslog=DIR+"run_lammps_"+jobname+".log";
  for(size_t jj=0;jj<Vcmd.size();jj++) {
    string str=" Running "+to_string(jj+1)+"/"+to_string(Vcmd.size())+" job...";
    cmd="echo \""+str+"\" >> "+flammpslog;
    cbs->execute_system_command(cmd);
    cmd=Vcmd[jj]+" -sc none >> "+flammpslog;
    cbs->run_lammps(cmd);
  }

  ReadLAMMPS readlmp_one(cbs);
  readlmp_one.read_lammpstrj(Vtrjs_list);
  for(size_t ii=0;ii<Vmol_cg.size();ii++) {
    vector<Atom> Vatoms_xyz=readlmp_one.VVatoms[ii].back();
    if(Vatoms_xyz.size()==0) continue;
    vector<double> Vxyz;
    for(size_t jj=0;jj<Vatoms_xyz.size();jj++) {
      vector<double> xyz=Vatoms_xyz[jj].get_atom_xyz();
      Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
    }
    Vmol_cg[ii].reset_molecule_xyz(Vxyz);
  }

  vector<Molecule> Vmol_aa;
  for(size_t ii=0;ii<Vmol_cg.size();ii++) {
    RNARecovery rnarec_one(cbs);
    if(cbs->init->RNA_version=="latest")
      rnarec_one.recover_RNA_to_AA(Vmol_cg[ii]);
    else if(cbs->init->RNA_version=="v2")
      rnarec_one.recover_RNA_to_AA_v2(Vmol_cg[ii]);
    Vmol_aa.push_back(rnarec_one.rna_aa);
  }

  ofstream outP (fout.c_str());
  for(size_t ii=0;ii<Vmol_aa.size();ii++) {
    outP<<"MODEL     "<<ii+1<<endl;
    outP<<Vmol_aa[ii];
    outP<<"ENDMDL"<<endl;
  }
  outP<<"END"<<endl;
  outP.close(); outP.clear();

  if(flag_osep) {
    for(size_t ii=0;ii<Vmol_aa.size();ii++) {
      string fstr=fout.substr(0,fout.size()-4);
      fstr+="_"+to_string(ii+1)+".pdb";
      outP.open(fstr.c_str());
      outP<<Vmol_aa[ii];
      outP.close(); outP.clear();
    }
  }

  if(flag_min==0) {
    time (&tend);

    double dtime=difftime(tend,tstart);
    int dH=dtime/3600;
    int dM=(dtime-dH*3600)/60;
    int dS=dtime-dH*3600-dM*60;
    string stime="Runtime: "+to_string(dH)+" Hours "+to_string(dM)+" Minutes "+to_string(dS)+" Seconds.";
    cbs->error->log(stime);
    stime="";
    cbs->error->log(stime);

    delete cbs;
    
    return 0;
  }
    
  if(!ReadAmber::check_tleap()) {
    string str="Cannot find tleap, skip energy minimization for allatom structure!";
    cbs->error->warning(str);

    time (&tend);

    double dtime=difftime(tend,tstart);
    int dH=dtime/3600;
    int dM=(dtime-dH*3600)/60;
    int dS=dtime-dH*3600-dM*60;
    string stime="Runtime: "+to_string(dH)+" Hours "+to_string(dM)+" Minutes "+to_string(dS)+" Seconds.";
    cbs->error->log(stime);
    stime="";
    cbs->error->log(stime);

    delete cbs;
    
    return 1;
  }

  prelmp_one.clear();
  prelmp_one.md_label=PrepareLAMMPS::em;
  vector<vector<Atom> > Vmol_atoms;
  fname=jobname+"_min_aa";
  for(size_t ii=0;ii<Vmol_aa.size();ii++) {
    string fname_amber=DIR+fname+to_string(ii+1);
    vector<string> Vdata=prelmp_one.prepare_data_from_amber(Vmol_aa[ii],fname_amber);
    string fstr=DIR+"data."+fname+"."+to_string(ii+1);
    ofstream outF (fstr.c_str());
    for(auto sline : Vdata)
      outF<<sline<<endl;
    outF.close(); outF.clear();

    fstr=fname_amber+"_amber.pdb"; //all atom structure with H atoms, consistent with data.min.$i
    vector<string> Vlines=read_file(fstr);
    vector<Atom> Vatoms;
    for(auto sline : Vlines) {
      if(sline.find("ATOM")==string::npos && sline.find("HETATM")==string::npos) continue;
      Atom atom_one;
      atom_one.set_atom(sline);
      Vatoms.push_back(atom_one);
    }
    Vmol_atoms.push_back(Vatoms);
  }

  vector<string> Vruns;
  int nthread=cbs->init->nthread;
  int nloop=(int)Vmol_aa.size()/nthread;
  int nlast=(int)Vmol_aa.size()%nthread;
  for(int jj=0;jj<nloop;jj++) {
    string str=to_string(jj*nthread+1);
    for(int kk=1;kk<nthread;kk++)
      str+=" "+to_string(jj*nthread+1+kk);
    prelmp_one.prepare_in_script_allatom(DIR,fname,str);
    string fstr=DIR+"in."+fname+"."+to_string(jj+1);
    prelmp_one.write_in_script(fstr);
    string srun=fstr+" "+to_string(nthread)+" 1";
    Vruns.push_back(srun);
  }
  if(nlast>0) {
    string str=to_string(nloop*nthread+1);
    for(int kk=1;kk<nlast;kk++)
      str+=" "+to_string(nloop*nthread+1+kk);
    prelmp_one.prepare_in_script_allatom(DIR,fname,str);
    string fstr=DIR+"in."+fname+"."+to_string(nloop+1);
    prelmp_one.write_in_script(fstr);
    string srun=fstr+" "+to_string(nlast)+" 1";
    Vruns.push_back(srun);
  }

  fsh=cbs->init->WORK_DIR+"run_lammps_"+fname+".sh";
  Vcmd=prelmp_one.write_run_script(Vruns,fsh);
  for(size_t jj=0;jj<Vcmd.size();jj++) {
    string str=" Running "+to_string(jj+1)+"/"+to_string(Vcmd.size())+" job...";
    cmd="echo \""+str+"\" >> "+flammpslog;
    cbs->execute_system_command(cmd);
    cmd=Vcmd[jj]+" -sc none >> "+flammpslog;
    cbs->run_lammps(cmd);
  }

  readlmp_one.clear();
  readlmp_one.read_lammpstrj(prelmp_one.Voutput_trjs);

  if(readlmp_one.VVatoms.size()!=Vmol_atoms.size()) {
    string str="Number of models is not matched during energy minimization --> ";
    str+=to_string(readlmp_one.VVatoms.size())+" / "+to_string(Vmol_atoms.size());
    cbs->error->one(str);
  }
  for(size_t ii=0;ii<Vmol_atoms.size();ii++) {
    vector<Atom> &Vatoms=readlmp_one.VVatoms[ii].back();
    for(size_t kk=0;kk<Vmol_atoms[ii].size();kk++) {
      Vmol_atoms[ii][kk].reset_atom_xyz(Vatoms[kk]);
    }
  }

  for(size_t ii=0;ii<Vmol_aa.size();ii++) {
    Vmol_aa[ii].reset_molecule_xyz(Vmol_atoms[ii]);
  }

  string fout_min=fout.substr(0,fout.size()-4);
  fout_min+="_min.pdb";
  outP.open(fout_min.c_str());
  for(size_t ii=0;ii<Vmol_aa.size();ii++) {
    outP<<"MODEL     "<<ii+1<<endl;
    outP<<Vmol_aa[ii];
    outP<<"ENDMDL"<<endl;
  }
  outP<<"END"<<endl;
  outP.close(); outP.clear();

  if(flag_osep) {
    for(size_t ii=0;ii<Vmol_aa.size();ii++) {
      string fstr=fout_min.substr(0,fout.size()-4);
      fstr+="_"+to_string(ii+1)+".pdb";
      outP.open(fstr.c_str());
      outP<<Vmol_aa[ii];
      outP.close(); outP.clear();
    }
  }

  time (&tend);

  double dtime=difftime(tend,tstart);
  string stime=cbs->calculate_runtime(dtime);
  cbs->error->log(stime);
  stime="";
  cbs->error->log(stime);

  delete cbs;
  
  return 0;
}