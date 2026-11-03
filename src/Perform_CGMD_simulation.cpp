#include "basic.h"
#include "error.h"
#include "readfile.h"
#include "initialization.h"
#include "read_PDB.h"
#include "sort_vector.h"
#include "measure.h"
#include "read_restraint.h"
#include "RNA_molecule.h"
#include "RNA_coarse_grained.h"
#include "RNA_basepair.h"
#include "RNA_motif.h"
#include "parse_RNA_motif.h"
#include "prepare_lammps.h"
#include "read_amber.h"
#include "read_lammps.h"

using namespace CBS_NS;

const size_t NMAX_RIGID = 29;
//const int NSTEP_MAX = 30000000;

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
  string seq="", secstr="";
  if(Mcmd.count("-i")) {
    fin=Mcmd["-i"][0];
    Mcmd.erase("-i");
  }
  if(Mcmd.count("-seq")) {
    seq=Mcmd["-seq"][0];
    Mcmd.erase("-seq");
  }
  if(Mcmd.count("-secstr")) {
    secstr=Mcmd["-secstr"][0];
    Mcmd.erase("-secstr");
  }
  if(fin.length()>1) {
    vector<string> Vlines=read_file(fin);
    int ncount=0;
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
      else {
        cout<<"Wrong input file --> "<<fin<<endl;
        exit(0);
      }
    }
  }

  if(seq.length()==0) {
    cout<<"Wrong commands: No input sequence (-i or -seq)!"<<endl;
    exit(0);
  }

  string inPDB;
  if(Mcmd.count("-if")) {
    inPDB=Mcmd["-if"][0];
    Mcmd.erase("-if");
  }
  else {
    cout<<"Wrong commands: No input initial structure (-if)!"<<endl;
    exit(0);
  }

  string input_format="aa";
  if(Mcmd.count("-iformat")) {
    string str=Mcmd["-iformat"][0];
    Mcmd.erase("-iformat");
    if(str=="aa") input_format="aa";
    else if(str=="cg") input_format="cg";
    else {
      cout<<"Wrong commands: Unknown input initial structure format (only aa/cg)!"<<endl;
      exit(0);
    }
  }

  int nmodel=-1;
  if(Mcmd.count("-nmodel")) {
    nmodel=stoi(Mcmd["-nmodel"][0]);
    Mcmd.erase("-nmodel");
  }

  int flag_break_check=0;
  if(Mcmd.count("-break_check")) {
    flag_break_check=1;
    Mcmd.erase("-break_check");
  }

  string jobname="rna";
  if(Mcmd.count("-name")) {
    jobname=Mcmd["-name"][0];
    Mcmd.erase("-name");
  }
  string outDIR="./simulation/";
  if(Mcmd.count("-DIR")) {
    outDIR=Mcmd["-DIR"][0];
    Mcmd.erase("-DIR");
  }
  if(outDIR.size()!=0 && outDIR.back()!='/')
    outDIR+="/";
  
  string frestr="";
  if(Mcmd.count("-restr")) {
    frestr=Mcmd["-restr"][0];
    Mcmd.erase("-restr");
  }

  int flag_template_rigid=0;
  if(Mcmd.count("-template_rigid")) {
    flag_template_rigid=1;
    Mcmd.erase("-template_rigid");
  }
 
  int Nstep=0;
  int flag_step=0;
  if(Mcmd.count("-nstep")) {
    Nstep=stoi(Mcmd["-nstep"][0]);
    flag_step=1;
    Mcmd.erase("-nstep");
  }

  int Nstep_restr=-1;
  if(Mcmd.count("-nstep_restr")) {
    Nstep_restr=stoi(Mcmd["-nstep_restr"][0]);
    Mcmd.erase("-nstep_restr");
  }

  double Tstart=300.0;
  if(Mcmd.count("-temp")) {
    Tstart=stod(Mcmd["-temp"][0]);
    Mcmd.erase("-temp");
  }

  int md_label=PrepareLAMMPS::remd;
  if(Mcmd.count("-m")) {
    string slabel=Mcmd["-m"][0];
    if(slabel=="md")
      md_label=PrepareLAMMPS::md;
    else if(slabel=="remd")
      md_label=PrepareLAMMPS::remd;
    else if(slabel=="em")
      md_label=PrepareLAMMPS::em;
    else {
      cout<<"Wrong commands: only md/remd/em allowed! --> "<<slabel<<endl;
      exit(0);
    }
    Mcmd.erase("-m");
  }

  int flag_ti=0;  //folding teritary interactions
  if(Mcmd.count("-ti_fold")) {
    flag_ti=1;
    Mcmd.erase("-ti_fold");
  }

  if(secstr.size()==0 && flag_ti) {
    cout<<"Wrong commands: No input 2D structure (-secstr or -i) for teritary interactions folding!"<<endl;
    exit(0);
  }

  int flag_ssrestr=0; //2D structure as restraints
  if(Mcmd.count("-ss_restr")) {
    flag_ssrestr=1;
    Mcmd.erase("-ss_restr");
  }

  if(secstr.size()==0 && flag_ssrestr) {
    cout<<"Wrong commands: No input 2D structure (-secstr or -i) as restraints!"<<endl;
    exit(0);
  }

  int Nout=1000;
  if(Mcmd.count("-nout")) {
    Nout=stoi(Mcmd["-nout"][0]);
    Mcmd.erase("-nout");
  }

  int Nout_restr=50;
  if(Mcmd.count("-nout_restr")) {
    Nout_restr=stoi(Mcmd["-nout_restr"][0]);
    Mcmd.erase("-nout_restr");
  }

  int flag_gpu=0;
  if(Mcmd.count("-gpu")) {
    flag_gpu=1;
    Mcmd.erase("-gpu");
  }

  int flag_exe=1;
  if(Mcmd.count("-noexe")) {
    flag_exe=0;
    Mcmd.erase("-noexe");
  }
  int flag_exe_restr=1;
  if(Mcmd.count("-noexe_restr")) {
    flag_exe_restr=0;
    Mcmd.erase("-noexe_restr");
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

  string slog="##### Perform_CGMD_simulation #####";
  cbs->error->log(slog);

  if(!is_directory_exist(outDIR)) {
    string cmd="mkdir "+outDIR;
    cbs->execute_system_command(cmd);
  }

  cbs->init->WORK_DIR=outDIR;
  cbs->init->RNA_check_chain_break=flag_break_check;
  cbs->init->user_template_rigid=flag_template_rigid;

  if(cbs->init->RNA_version=="v2") {
    flag_ssrestr=1;
    if(secstr.length()==0) {
      string str="No input 2D structure as restraints for IsRNA2!";
      cbs->error->warning(str);
    }
  }

  RNAMolecule rna_one(cbs);
  rna_one.set_RNA_molecule(seq);
  size_t nres_total=rna_one.number_of_residues();
  string seq_in=rna_one.get_sequences();
  vector<vector<string> > Vres=rna_one.get_residue_name();
  RNABasePair rna_bp_one(cbs);
  rna_bp_one.parse_RNA_basepairs(Vres,secstr);
  map<string,string> MBPs_ti;
  for(auto & it : rna_bp_one.mapBP) {
    int bp_level=rna_bp_one.get_basepair_level(it.first,it.second);
    if(bp_level>1) {
      MBPs_ti.emplace(it.first,it.second);
    }
  }

  ParseRNAMotif PRmotif_one(cbs);
  PRmotif_one.parse_RNA_motifs(rna_bp_one);
  int topo=PRmotif_one.get_RNA_topology();
  vector<string> Vresnames_helix;
  for(size_t jj=0;jj<PRmotif_one.motifs.size();jj++) {
    string type=PRmotif_one.motifs[jj].type;
    if(type=="Helix") {
      vector<string> &Vresnames_one=PRmotif_one.motifs[jj].Vresnames;
      Vresnames_helix.insert(Vresnames_helix.end(),Vresnames_one.begin(),Vresnames_one.end());
    }
  }

  ReadPDB readpdb_one(cbs);
  readpdb_one.read_PDB(inPDB);
  if(readpdb_one.RNA_atoms.size()==0 || readpdb_one.RNA_atoms[0].size()==0) {
    string str="No RNA atoms in "+inPDB;
    cbs->error->one(str);
  }

  if(nmodel<0) 
    nmodel=static_cast<int>(readpdb_one.RNA_atoms.size());
  if((int)readpdb_one.RNA_atoms.size()<nmodel)
    nmodel=static_cast<int>(readpdb_one.RNA_atoms.size());
  
  vector<RNACoarseGrained> Vrnacg;
  for(int ii=0;ii<nmodel;ii++) {
    RNACoarseGrained rnacg_one(cbs);
    if(input_format=="aa") {
      RNAMolecule rna_aa(cbs);
      rna_aa.set_RNA_molecule(readpdb_one.RNA_atoms[ii]);
      if(cbs->init->RNA_version=="latest")
        rnacg_one.convert_RNA_to_CG(rna_aa);
      else if(cbs->init->RNA_version=="v2")
        rnacg_one.convert_RNA_to_CG_v2(rna_aa);
    }
    else if(input_format=="cg") {
      Molecule rna_cg(cbs);
      for(size_t jj=0;jj<readpdb_one.RNA_atoms[ii].size();jj++) {
        string resname_one=readpdb_one.RNA_atoms[ii][jj].resname;
        if(resname_one.size()==1) resname_one="c"+resname_one;
        readpdb_one.RNA_atoms[ii][jj].resname=resname_one;
      }
      rnacg_one.reset_RNA_CG(readpdb_one.RNA_atoms[ii]);
    }
    Vrnacg.push_back(rnacg_one);
  }  

  string foutname=outDIR+jobname;
  string fstr=foutname+"_cg.pdb";
  ofstream outF (fstr.c_str());
  for(int ii=0;ii<nmodel;ii++) {
    outF<<"MODEL   "<<ii+1<<endl;
    outF<<Vrnacg[ii].rna_cg;
    outF<<"ENDMDL"<<endl;
  }
  outF<<"END"<<endl;
  outF.close(); outF.clear();

  for(int ii=0;ii<nmodel;ii++) {
    string seq_stru=Vrnacg[ii].rna_cg.get_sequence();
    if(seq_stru!=seq_in) {
      string str="Sequences are not matched --> "+seq_in+"(-seq) / "+seq_stru;
      str+="(Model #"+to_string(ii+1)+" in "+inPDB+")";
      cbs->error->one(str);
    }
    if(cbs->init->RNA_version=="latest")
      Vrnacg[ii].build_RNA_topology();
    else if(cbs->init->RNA_version=="v2")
      Vrnacg[ii].build_RNA_topology_v2();
    if(cbs->init->CircRNA_chains.size()>0) {
      if(cbs->init->RNA_version=="latest")
        Vrnacg[ii].add_circular_RNA_topology();
      else if(cbs->init->RNA_version=="v2")
      Vrnacg[ii].add_circular_RNA_topology_v2();
    }
  }
 
  ReadCoeff readcoeff_one(cbs);
  readcoeff_one.read_RNA_coeff();
  readcoeff_one.get_interaction_style();
  PrepareLAMMPS prelmp_one(cbs);
  prelmp_one.set_coeff(readcoeff_one);
  if(md_label==PrepareLAMMPS::md) {
    prelmp_one.set_temperature(Tstart,Tstart);
  }
  else if(md_label==PrepareLAMMPS::remd) {
    prelmp_one.set_remd_temperature();
  }
  else 
    prelmp_one.md_label=PrepareLAMMPS::em;

  vector<string> Vgroup,Vfix,Vrun,Vrun_timestep;
  ReadRestraint readrestr_one(cbs);
  if(frestr.size()>0) {
    cbs->init->RNA_check_chain_break=1;
    readrestr_one.read_restraint(frestr);
    cbs->init->RNA_check_chain_break=flag_break_check;
  }

  if(readrestr_one.RNA_templates.size()>0) {
    for(size_t jj=0;jj<readrestr_one.RNA_templates.size();jj++) {
      for(size_t kk=0;kk<readrestr_one.RNA_templates[jj].chains.size();kk++) {
        string res_begin=readrestr_one.RNA_templates[jj].chains[kk].first_of_residue();
        string res_end=readrestr_one.RNA_templates[jj].chains[kk].last_of_residue();
        string secstr_one=rna_bp_one.get_segment_secondary_structure(res_begin,res_end);
        readrestr_one.RNA_templates[jj].secstrs[kk]=secstr_one;
      }
    }
  }
  
  size_t nrigidbody=readrestr_one.rigid_bodies.size();
  
  //for(size_t ll=0;ll<nrigidbody;ll++) {
    //for(size_t mm=0;mm<readrestr_one.rigid_bodies[ll].size();mm++) 
      //cout<<readrestr_one.rigid_bodies[ll][mm].residue_i<<":"<<readrestr_one.rigid_bodies[ll][mm].residue_j<<"; ";
    //cout<<endl;
  //}

  if(nrigidbody>NMAX_RIGID) {
    nrigidbody=NMAX_RIGID;
    vector<vector<ResiduePair> > Trigid_bodies=readrestr_one.rigid_bodies;
    readrestr_one.rigid_bodies.clear();
    vector<double> Vnres_total;
    for(size_t jj=0;jj<Trigid_bodies.size();jj++) {
      double nres_total=0.0;
      for(size_t kk=0;kk<Trigid_bodies[jj].size();kk++) {
        string &resname_begin=Trigid_bodies[jj][kk].residue_i;
        string &resname_end=Trigid_bodies[jj][kk].residue_j;
        nres_total+=rna_one.number_of_residues(resname_begin,resname_end);
      }
      nres_total*=(1.0+Trigid_bodies[jj].size())/2.0;
      Vnres_total.push_back(nres_total);
    }
    vector<size_t> Vsort=SortVector::sort_vector_descending(Vnres_total);
    for(size_t jj=0;jj<nrigidbody;jj++) {
      readrestr_one.rigid_bodies.push_back(Trigid_bodies[Vsort[jj]]);
    }

    string str="Reserve top "+to_string(NMAX_RIGID)+" rigidbodies with most residues:";
    cbs->error->warning(str);
    for(size_t jj=0;jj<nrigidbody;jj++) {
      str="   "+readrestr_one.rigid_bodies[jj][0].residue_i+"-"+readrestr_one.rigid_bodies[jj][0].residue_j;
      for(size_t kk=1;kk<readrestr_one.rigid_bodies[jj].size();kk++)
        str+=", "+readrestr_one.rigid_bodies[jj][kk].residue_i+"-"+readrestr_one.rigid_bodies[jj][kk].residue_j;
      cbs->error->warning(str);
    }
  }

  size_t nres_rigid=0;
  vector<string> Vresnames_rigid;
  vector<vector<int> > Vrigid_atomids;
  if(nrigidbody==1) {
    string atomid_list="";
    vector<int> Vatomids_rigid_i;
    for(size_t jj=0;jj<readrestr_one.rigid_bodies[0].size();jj++) {
      string resname_begin=readrestr_one.rigid_bodies[0][jj].residue_i;
      string resname_end=readrestr_one.rigid_bodies[0][jj].residue_j;
      //nres_rigid+=rna_one.number_of_residues(resname_begin,resname_end);
      vector<string> Vresnames_one=rna_one.get_residue_list(resname_begin,resname_end);
      Vresnames_rigid.insert(Vresnames_rigid.end(),Vresnames_one.begin(),Vresnames_one.end());
      string atomid_list_j=Vrnacg[0].rna_cg.get_atomid_list(resname_begin,resname_end);
      atomid_list+=" "+atomid_list_j;
      vector<string> stoken=split_string(atomid_list_j,":");
      if(stoken.size()%2!=0) {
        string str="Wrong atomid list (id_begin:id_end) --> "+atomid_list_j;
        cbs->error->one(str);
      }
      for(size_t kk=0;kk<stoken.size()/2;kk++) {
        int atomid_begin=stoi(stoken[2*kk]);
        int atomid_end=stoi(stoken[2*kk+1]);
        for(int ll=atomid_begin;ll<=atomid_end;ll++)
          Vatomids_rigid_i.push_back(ll);
      }
    }
    Vrigid_atomids.push_back(Vatomids_rigid_i);
    string sline="group           rb id"+atomid_list;
    Vgroup.push_back(sline);
    
    ostringstream oss;
    oss<<"fix             myrigid rb rigid single langevin "<<std::fixed<<std::setprecision(2)<<Tstart<<" ";
    oss<<std::fixed<<std::setprecision(2)<<Tstart<<" 100.0 "<<57412+(int)(743891*cbs->rand_one->uniform());
    sline=oss.str();
    Vfix.push_back(sline);
    oss.str("");
    oss.clear();    
  }
  else if(nrigidbody>1) {
    vector<int> Vatomids_rigid_all;
    for(size_t ii=0;ii<readrestr_one.rigid_bodies.size();ii++) {
      string atomid_list="";
      vector<int> Vatomids_rigid_i;
      for(size_t jj=0;jj<readrestr_one.rigid_bodies[ii].size();jj++) {
        string resname_begin=readrestr_one.rigid_bodies[ii][jj].residue_i;
        string resname_end=readrestr_one.rigid_bodies[ii][jj].residue_j;
        //nres_rigid+=rna_one.number_of_residues(resname_begin,resname_end);
        vector<string> Vresnames_one=rna_one.get_residue_list(resname_begin,resname_end);
        Vresnames_rigid.insert(Vresnames_rigid.end(),Vresnames_one.begin(),Vresnames_one.end());
        string atomid_list_j=Vrnacg[0].rna_cg.get_atomid_list(resname_begin,resname_end);
        atomid_list+=" "+atomid_list_j;
        vector<string> stoken=split_string(atomid_list_j,":");
        if(stoken.size()%2!=0) {
          string str="Wrong atomid list (id_begin:id_end) --> "+atomid_list_j;
          cbs->error->one(str);
        }
        for(size_t kk=0;kk<stoken.size()/2;kk++) {
          int atomid_begin=stoi(stoken[2*kk]);
          int atomid_end=stoi(stoken[2*kk+1]);
          for(int ll=atomid_begin;ll<=atomid_end;ll++)
            Vatomids_rigid_i.push_back(ll);
        }
      }
      Vrigid_atomids.push_back(Vatomids_rigid_i);
      if(ii==0) {
        Vatomids_rigid_all=Vatomids_rigid_i;
      }
      else {
        int relation=Measure::vector_relation(Vatomids_rigid_all,Vatomids_rigid_i);
        if(relation!=Measure::separated) {
          string str="There are overlapped atoms between rigid bodies --> ";
          str+=atomid_list+" / "+to_string(Vatomids_rigid_i[0])+":"+to_string(Vatomids_rigid_i.back());
          cbs->error->one(str);
        }
      }
      string sline="group           rb"+to_string(ii+1)+" id"+atomid_list;
      Vgroup.push_back(sline);
    }
    string sline="group           rb union";
    string srb="";
    for(size_t ii=0;ii<nrigidbody;ii++)
      srb+=" rb"+to_string(ii+1);
    sline+=srb;
    Vgroup.push_back(sline);

    ostringstream oss;
    oss<<"fix             myrigid rb rigid group "<<nrigidbody<<srb;
    oss<<" langevin "<<std::fixed<<std::setprecision(2)<<Tstart<<" ";
    oss<<std::fixed<<std::setprecision(2)<<Tstart<<" 100.0 "<<57412+(int)(743891*cbs->rand_one->uniform());
    sline=oss.str();
    Vfix.push_back(sline);
    oss.str("");
    oss.clear();
  }

  nres_rigid=Vresnames_rigid.size();
  int flag_group_other=1;
  if(nres_rigid==nres_total) {
    flag_group_other=0;
    md_label=PrepareLAMMPS::md;
    prelmp_one.set_temperature(Tstart,Tstart);
    for(size_t jj=0;jj<Vgroup.size()-1;jj++) {
      vector<string> stoken=split_string(Vgroup[jj]);
      int flag_check=0;
      for(size_t kk=3;kk<stoken.size();kk++){
        vector<string> Vatomid_k=split_string(stoken[kk],":");
        if(Vatomid_k[0]=="1") {
          flag_check=1;
          stoken[kk]="2:"+Vatomid_k[1];
          break;
        }
      }
      if(flag_check) {
        string sline="group           "+stoken[1]+" "+stoken[2];
        for(size_t kk=3;kk<stoken.size();kk++)
          sline+=" "+stoken[kk];
        Vgroup[jj]=sline;
        break;
      }
    }
  }
  
  if(flag_group_other && Vgroup.size()>0) {
    string sline="group           other subtract all rb";
    Vgroup.push_back(sline);
  }  
  
  if(nrigidbody>0) {
    ostringstream oss;
    string sline;
    if(md_label==PrepareLAMMPS::remd && flag_group_other) {
      oss<<"fix             mylang other langevin $t $t 100.0 "<<57412+(int)(743891*cbs->rand_one->uniform());
      sline=oss.str();
      Vfix.push_back(sline);
      oss.str("");
      oss.clear();
    }
    else if(md_label==PrepareLAMMPS::md && flag_group_other) {
      oss<<"fix             mylang other langevin "<<std::fixed<<std::setprecision(2)<<Tstart<<" ";
      oss<<std::fixed<<std::setprecision(2)<<Tstart<<" 100.0 "<<57412+(int)(743891*cbs->rand_one->uniform());
      sline=oss.str();
      Vfix.push_back(sline);
      oss.str("");
      oss.clear();
    }
    if(flag_group_other) {
      sline="fix             mynve other nve/limit 0.2";
      Vfix.push_back(sline);
      sline="";
      Vfix.push_back(sline);
      sline="thermo_style    custom step temp epair emol etotal press f_mynve";
      Vfix.push_back(sline);
    }
    else {
      sline="thermo_style    custom step temp epair emol etotal press";
      Vfix.push_back(sline);
    }

    oss<<"velocity        rb create "<<std::fixed<<std::setprecision(2)<<Tstart<<" ";
    oss<<37624+(int)(913240*cbs->rand_one->uniform())<<" rot yes dist gaussian";
    sline=oss.str();
    Vrun.push_back(sline);
    oss.str("");
    oss.clear();

    if(md_label==PrepareLAMMPS::remd && flag_group_other) {
      oss<<"velocity        other create $t "<<37624+(int)(913240*cbs->rand_one->uniform())<<" rot yes dist gaussian";
      sline=oss.str();
      Vrun.push_back(sline);
      oss.str("");
      oss.clear();
    }
    else if(md_label==PrepareLAMMPS::md && flag_group_other) {
      oss<<"velocity        other create "<<std::fixed<<std::setprecision(2)<<Tstart<<" ";
      oss<<37624+(int)(913240*cbs->rand_one->uniform())<<" rot yes dist gaussian";
      sline=oss.str();
      Vrun.push_back(sline);
      oss.str("");
      oss.clear();
    }
  }

  size_t nres_helix=0;
  for(size_t jj=0;jj<Vresnames_helix.size();jj++) {
    if(count(Vresnames_rigid.begin(),Vresnames_rigid.end(),Vresnames_helix[jj])==0)
      nres_helix++;
  }

  size_t nres_eff=nres_total-nres_rigid-nres_helix/2;
  if(flag_step==0) {
    Nstep=10000000;
    if(topo==ParseRNAMotif::sl) {
      if(nres_eff<=5) Nstep=500000;
      else if(nres_eff>5 && nres_eff<=10) Nstep=1000000;
      else if(nres_eff>10 && nres_eff<=20) Nstep=2000000;
      else if(nres_eff>20 && nres_eff<=50) Nstep=3000000;
      else if(nres_eff>50 && nres_eff<=100) Nstep=4500000;
      else if(nres_eff>100 && nres_eff<=150) Nstep=6000000;
      else if(nres_eff>150 && nres_eff<=200) Nstep=7500000;
      else if(nres_eff>200 && nres_eff<=250) Nstep=9000000;
      else if(nres_eff>250 && nres_eff<=300) Nstep=10500000;
      else Nstep=12000000;
      //Nstep=round(8+0.35*nres_eff)*100000;
    }
    else if(topo==ParseRNAMotif::slfl) {
      if(nres_eff<=5) Nstep=1000000;
      else if(nres_eff>5 && nres_eff<=10) Nstep=2000000;
      else if(nres_eff>10 && nres_eff<=20) Nstep=4000000;
      else if(nres_eff>20 && nres_eff<=50) Nstep=6000000;
      else if(nres_eff>50 && nres_eff<=100) Nstep=8000000;
      else if(nres_eff>100 && nres_eff<=150) Nstep=10000000;
      else if(nres_eff>150 && nres_eff<=200) Nstep=12000000;
      else if(nres_eff>200 && nres_eff<=250) Nstep=14000000;
      else if(nres_eff>250 && nres_eff<=300) Nstep=16000000;
      else Nstep=18000000;
      //Nstep=round(30+0.55*nres_eff)*100000;
    }
    else if(topo==ParseRNAMotif::mw) {
      if(nres_eff<=5) Nstep=500000;
      else if(nres_eff>5 && nres_eff<=10) Nstep=1000000;
      else if(nres_eff>10 && nres_eff<=20) Nstep=3500000;
      else if(nres_eff>20 && nres_eff<=50) Nstep=6000000;
      else if(nres_eff>50 && nres_eff<=100) Nstep=8500000;
      else if(nres_eff>100 && nres_eff<=150) Nstep=11000000;
      else if(nres_eff>150 && nres_eff<=200) Nstep=13500000;
      else if(nres_eff>200 && nres_eff<=250) Nstep=16000000;
      else if(nres_eff>250 && nres_eff<=300) Nstep=18500000;
      else Nstep=21000000;
      //Nstep=round(30+0.55*nres_eff)*100000;
    }
    else if(topo==ParseRNAMotif::mwfl) {
      if(nres_eff<=5) Nstep=1000000;
      else if(nres_eff>5 && nres_eff<=10) Nstep=2000000;
      else if(nres_eff>10 && nres_eff<=20) Nstep=4500000;
      else if(nres_eff>20 && nres_eff<=50) Nstep=7000000;
      else if(nres_eff>50 && nres_eff<=100) Nstep=9500000;
      else if(nres_eff>100 && nres_eff<=150) Nstep=12000000;
      else if(nres_eff>150 && nres_eff<=200) Nstep=14500000;
      else if(nres_eff>200 && nres_eff<=250) Nstep=17000000;
      else if(nres_eff>250 && nres_eff<=300) Nstep=19500000;
      else Nstep=22000000;
      //Nstep=round(40+0.75*nres_eff)*100000;
    }
    else if(topo==ParseRNAMotif::ti) {
      if(nres_eff<=5) Nstep=2000000;
      else if(nres_eff>5 && nres_eff<=10) Nstep=4000000;
      else if(nres_eff>10 && nres_eff<=20) Nstep=6000000;
      else if(nres_eff>20 && nres_eff<=50) Nstep=10000000;
      else if(nres_eff>50 && nres_eff<=100) Nstep=14000000;
      else if(nres_eff>100 && nres_eff<=150) Nstep=18000000;
      else if(nres_eff>150 && nres_eff<=200) Nstep=22000000;
      else if(nres_eff>200 && nres_eff<=250) Nstep=26000000;
      else if(nres_eff>250 && nres_eff<=300) Nstep=30000000;
      else Nstep=34000000;
      //Nstep=round(50+1.00*nres_eff)*100000;
    }
    else if(topo==ParseRNAMotif::tifl) {
      if(nres_eff<=5) Nstep=4000000;
      else if(nres_eff>5 && nres_eff<=10) Nstep=6000000;
      else if(nres_eff>10 && nres_eff<=20) Nstep=8000000;
      else if(nres_eff>20 && nres_eff<=50) Nstep=10000000;
      else if(nres_eff>50 && nres_eff<=100) Nstep=14000000;
      else if(nres_eff>100 && nres_eff<=150) Nstep=18000000;
      else if(nres_eff>150 && nres_eff<=200) Nstep=22000000;
      else if(nres_eff>200 && nres_eff<=250) Nstep=26000000;
      else if(nres_eff>250 && nres_eff<=300) Nstep=30000000;
      else Nstep=34000000;
      //Nstep=round(70+1.20*nres_eff)*100000;
    }
 
    //if(Nstep > NSTEP_MAX) Nstep=NSTEP_MAX;
  }

  Vrun_timestep.clear();
  if(Nstep>=15000000) {
    int Nstep_one=round(Nstep/5);
    int Nstep_two=round(Nstep/5);
    int Nstep_three=Nstep-Nstep_one-Nstep_two;
    ostringstream oss;
    string sline;
    sline="timestep        4.0";
    Vrun_timestep.push_back(sline);
    if(md_label==PrepareLAMMPS::remd) {
      oss<<"temper          "<<Nstep_one<<" 5000 $t mylang "<<35824+(int)(165243*cbs->rand_one->uniform())<<" "<<86529+(int)(43280*cbs->rand_one->uniform());
      sline=oss.str();
      Vrun_timestep.push_back(sline);
      oss.str("");
      oss.clear();
    }
    else if(md_label==PrepareLAMMPS::md) {
      sline="run             "+to_string(Nstep_one);
      Vrun_timestep.push_back(sline);
    }
    sline="timestep        8.0";
    Vrun_timestep.push_back(sline);
    if(md_label==PrepareLAMMPS::remd) {
      oss<<"temper          "<<Nstep_two<<" 5000 $t mylang "<<35824+(int)(165243*cbs->rand_one->uniform())<<" "<<86529+(int)(43280*cbs->rand_one->uniform());
      sline=oss.str();
      Vrun_timestep.push_back(sline);
      oss.str("");
      oss.clear();
    }
    else if(md_label==PrepareLAMMPS::md) {
      sline="run             "+to_string(Nstep_two);
      Vrun_timestep.push_back(sline);
    }
    sline="timestep        10.0";
    Vrun_timestep.push_back(sline);
    if(md_label==PrepareLAMMPS::remd) {
      oss<<"temper          "<<Nstep_three<<" 5000 $t mylang "<<43817+(int)(297034*cbs->rand_one->uniform())<<" "<<91478+(int)(69035*cbs->rand_one->uniform());
      sline=oss.str();
      Vrun_timestep.push_back(sline);
      oss.str("");
      oss.clear();
    }
    else if(md_label==PrepareLAMMPS::md) {
      sline="run             "+to_string(Nstep_three);
      Vrun_timestep.push_back(sline);
    }
  }
  else if(Nstep<15000000 && Nstep>=8000000) {
    int Nstep_one=round(Nstep/5);
    int Nstep_two=Nstep-Nstep_one;
    ostringstream oss;
    string sline;
    sline="timestep        6.0";
    Vrun_timestep.push_back(sline);
    if(md_label==PrepareLAMMPS::remd) {
      oss<<"temper          "<<Nstep_one<<" 5000 $t mylang "<<35824+(int)(165243*cbs->rand_one->uniform())<<" "<<86529+(int)(43280*cbs->rand_one->uniform());
      sline=oss.str();
      Vrun_timestep.push_back(sline);
      oss.str("");
      oss.clear();
    }
    else if(md_label==PrepareLAMMPS::md) {
      sline="run             "+to_string(Nstep_one);
      Vrun_timestep.push_back(sline);
    }
    sline="timestep        10.0";
    Vrun_timestep.push_back(sline);
    if(md_label==PrepareLAMMPS::remd) {
      oss<<"temper          "<<Nstep_two<<" 5000 $t mylang "<<43817+(int)(297034*cbs->rand_one->uniform())<<" "<<91478+(int)(69035*cbs->rand_one->uniform());
      sline=oss.str();
      Vrun_timestep.push_back(sline);
      oss.str("");
      oss.clear();
    }
    else if(md_label==PrepareLAMMPS::md) {
      sline="run             "+to_string(Nstep_two);
      Vrun_timestep.push_back(sline);
    }
  }
  else {
    ostringstream oss;
    string sline;
    if(md_label==PrepareLAMMPS::remd) {
      sline="timestep        10.0";
      Vrun_timestep.push_back(sline);
      oss<<"temper          "<<Nstep<<" 5000 $t mylang "<<35824+(int)(165243*cbs->rand_one->uniform())<<" "<<86529+(int)(43280*cbs->rand_one->uniform());
      sline=oss.str();
      Vrun_timestep.push_back(sline);
      oss.str("");
      oss.clear();
    }
    else if(md_label==PrepareLAMMPS::md) {
      sline="timestep        8.0";
      Vrun_timestep.push_back(sline);
      sline="run             "+to_string(Nstep);
      Vrun_timestep.push_back(sline);
    }
  }

  vector<string> Vem_skip;  //skip energy minimization
  Vem_skip.push_back("");

  size_t ndistancerestraint=readrestr_one.distance_restraints.size();
  size_t nbp_restr=readrestr_one.RNA_basepairs.size();
  size_t nbs_restr=readrestr_one.RNA_basestacks.size();

  map<string,string> TBPs_ti=MBPs_ti;
  MBPs_ti.clear();
  for(auto & it : TBPs_ti) {
    string atomid_list_i=Vrnacg[0].rna_cg.get_atomid_list(it.first,it.first);
    string atomid_list_j=Vrnacg[0].rna_cg.get_atomid_list(it.second,it.second);
    vector<int> Vatomids_BP_one;
    vector<string> stoken=split_string(atomid_list_i,":");
    if(stoken.size()!=2) {
      string str="Wrong atomid list (id_begin:id_end) --> "+atomid_list_j;
      cbs->error->one(str);
    }
    int atomid_begin=stoi(stoken[0]);
    int atomid_end=stoi(stoken[1]);
    for(int kk=atomid_begin;kk<=atomid_end;kk++)
       Vatomids_BP_one.push_back(kk);
    stoken.clear();
    stoken=split_string(atomid_list_j,":");
    if(stoken.size()!=2) {
      string str="Wrong atomid list (id_begin:id_end) --> "+atomid_list_j;
      cbs->error->one(str);
    }
    atomid_begin=stoi(stoken[0]);
    atomid_end=stoi(stoken[1]);
    for(int kk=atomid_begin;kk<=atomid_end;kk++)
       Vatomids_BP_one.push_back(kk);
    
    int flag_check=0; 
    for(size_t jj=0;jj<Vrigid_atomids.size();jj++) {
      int relation=Measure::vector_relation(Vrigid_atomids[jj],Vatomids_BP_one);
      if(relation==Measure::equal || relation==Measure::contained) {
        flag_check=1;
        break;
      }
    }
    if(flag_check==0) {
      MBPs_ti.emplace(it.first,it.second);
    }
  }
  TBPs_ti.clear();
  
  if(ndistancerestraint>0 || nbp_restr>0 || nbs_restr>0 || (MBPs_ti.size()>0 && flag_ti)) {  //folding restraint interactions
    double distance_max=0.0;
    vector<RNACoarseGrained> Vrnacg_restr=Vrnacg;
    string sjobname=jobname+"_restr";
    if(nbp_restr>0 || (MBPs_ti.size()>0 && flag_ti))
      readcoeff_one.add_extra_RNA_basepair_coeff("res");
    prelmp_one.set_coeff(readcoeff_one);
    for(int ii=0;ii<nmodel;ii++) {
      if(ndistancerestraint>0) {
        Vrnacg_restr[ii].add_distance_restraint(readrestr_one.distance_restraints);
      }
      if(nbp_restr>0) {
        if(cbs->init->RNA_version=="latest")
          Vrnacg_restr[ii].add_RNA_basepairs(readrestr_one.RNA_basepairs,"res");
        else if(cbs->init->RNA_version=="v2")
          Vrnacg_restr[ii].add_RNA_basepairs_v2(readrestr_one.RNA_basepairs,"res");
      }
      if(nbs_restr>0) {
        Vrnacg_restr[ii].add_RNA_basestacks(readrestr_one.RNA_basestacks);
      }
      if(MBPs_ti.size()>0 && flag_ti) {
        if(cbs->init->RNA_version=="latest")
          Vrnacg_restr[ii].add_RNA_basepairs(MBPs_ti,"res");
        else if(cbs->init->RNA_version=="v2")
          Vrnacg_restr[ii].add_RNA_basepairs_v2(MBPs_ti,"res");
      }

      if(cbs->init->RNA_version=="v2") { //add base pairing interactions for IsRNA2
        map<string,string> MBPs_one;
        for(auto & it : rna_bp_one.mapBP) {
          if(MBPs_ti.count(it.first)) continue;
          MBPs_one.emplace(it.first,it.second);
        }
        Vrnacg_restr[ii].add_RNA_basepairs_v2(MBPs_one);
      }

      vector<int> Vatomids_restr;
      for(size_t jj=Vrnacg[ii].bonds.size();jj<Vrnacg_restr[ii].bonds.size();jj++) {
        int atom_i=Vrnacg_restr[ii].bonds[jj].atom_i;
        int atom_j=Vrnacg_restr[ii].bonds[jj].atom_j;
        Vatomids_restr.push_back(atom_i);
        Vatomids_restr.push_back(atom_j);
      }

      vector<Atom> Vatoms_list=Vrnacg_restr[ii].rna_cg.get_atom_list();
      map<int,vector<double> > Matomid_xyz;
      for(size_t kk=0;kk<Vatoms_list.size();kk++) {
        int atomid_one=Vatoms_list[kk].id;
        if(count(Vatomids_restr.begin(),Vatomids_restr.end(),atomid_one)) {
          vector<double> xyz=Vatoms_list[kk].get_atom_xyz();
          Matomid_xyz.emplace(atomid_one,xyz);
        }
      }
      for(size_t kk=0;kk<Vatomids_restr.size()/2;kk++) {
        int atomid_i=Vatomids_restr[2*kk];
        int atomid_j=Vatomids_restr[2*kk+1];
        if(Matomid_xyz.count(atomid_i)==0 || Matomid_xyz.count(atomid_j)==0) {
          string str="Cannot find xyz for restraint atom "+to_string(atomid_i)+"/"+to_string(atomid_j);
          cbs->error->one(str);
        }
        double dist_one=Measure::distance(Matomid_xyz[atomid_i],Matomid_xyz[atomid_j]);
        if(dist_one > distance_max)
          distance_max=dist_one;
      }

      vector<string> Vdata_one=prelmp_one.prepare_data(Vatoms_list,Vrnacg_restr[ii].bonds,Vrnacg_restr[ii].angles,Vrnacg_restr[ii].dihedrals,
                                                       Vrnacg_restr[ii].impropers);
      if(Vdata_one[0].find("Number of atom types")!=string::npos) {
        cout<<"Warning: Model #"<<ii+1<<" "<<Vdata_one[0]<<endl;
        Vdata_one.erase(Vdata_one.begin());
      }
      if(flag_exe_restr) {
        fstr=outDIR+"data."+sjobname+"."+to_string(ii+1);
        outF.open(fstr.c_str());
        for(size_t kk=0;kk<Vdata_one.size();kk++)
          outF<<Vdata_one[kk]<<endl;
        outF.close(); outF.clear();
      }
    }

    vector<string> Vrestr_bond_lable;
    vector<string> Vrestr_bond_style;
    vector<double> Vrestr_bond_k;
    vector<string> Vrestr_bond_b0;
    for(size_t jj=0;jj<ndistancerestraint;jj++) {
      string stype=readrestr_one.distance_restraints[jj].type+" "+readrestr_one.distance_restraints[jj].coeff;
      if(prelmp_one.Mrestraint_type.count(stype)==0) {
        string str="Cannot find bond type label for distance restraint "+stype;
        cbs->error->one(str);
      }
      string slabel=prelmp_one.Mrestraint_type[stype];
      double k_one=0.0;
      vector<string> stoken=split_string(readrestr_one.distance_restraints[jj].coeff);
      if(stoken.size()>0)
        k_one=stod(stoken[0]);
      string b0_one="";
      if(stoken.size()>1) {
        b0_one=stoken[1];
        for(size_t kk=2;kk<stoken.size();kk++)
          b0_one+=" "+stoken[kk];
      }
      if(count(Vrestr_bond_lable.begin(),Vrestr_bond_lable.end(),slabel)==0) {
        Vrestr_bond_lable.push_back(slabel);
        Vrestr_bond_style.push_back(readrestr_one.distance_restraints[jj].type);
        Vrestr_bond_k.push_back(k_one);
        Vrestr_bond_b0.push_back(b0_one);
      }
    }

    double WEIGHT[]={0.001,0.01,0.1,1.0};
    double WEIFHT_NSTEP[]={0.2,0.3,0.3,0.2};
    string TIMESTEP[]={"0.2","0.5","1.0","1.0"};
    int n_splited=get_array_size(WEIGHT);

    vector<vector<string> > Vrestr_bond_coeff;
    for(int jj=0;jj<n_splited;jj++) {
      vector<string> restr_bond_coeff;
      if(ndistancerestraint>0) {
        for(size_t kk=0;kk<Vrestr_bond_lable.size();kk++) {
          string scoeff=Vrestr_bond_lable[kk];
          if(prelmp_one.Vbondstyle.size()>1)
            scoeff+=" "+Vrestr_bond_style[kk];
          scoeff+=" "+to_string(Vrestr_bond_k[kk]*WEIGHT[jj])+" "+Vrestr_bond_b0[kk];
          restr_bond_coeff.push_back(scoeff);
        }
      }
      if(nbp_restr>0 || (MBPs_ti.size()>0 && flag_ti)) {
        vector<string> Vbp_coeff_weighted=readcoeff_one.get_weighted_RNA_basepair_coeff("res",WEIGHT[jj]);
        for(size_t kk=0;kk<Vbp_coeff_weighted.size();kk++) {
          vector<string> stoken=split_string(Vbp_coeff_weighted[kk]);
          string type_one=stoken[0];
          string style_one=stoken[1];
          string coeff_one=stoken[2];
          for(size_t mm=3;mm<stoken.size();mm++)
            coeff_one+=" "+stoken[mm];
          if(count(prelmp_one.Vbondtype.begin(),prelmp_one.Vbondtype.end(),type_one)) {
            string scoeff=type_one;
            if(prelmp_one.Vbondstyle.size()>1)
              scoeff+=" "+style_one;
            scoeff+=" "+coeff_one;
            restr_bond_coeff.push_back(scoeff);
          }
        }
      }
      
      Vrestr_bond_coeff.push_back(restr_bond_coeff);
    }

    if(Nstep_restr<0) {
      distance_max=distance_max-2.0;
      if(distance_max<0.0) distance_max=0.0;
      Nstep_restr=1000000+round(distance_max*n_splited)*150000;
      if(Nstep_restr>n_splited*5000000) Nstep_restr=n_splited*5000000;
    }
    int Nevery_restr=Nstep_restr/Nout_restr;
    vector<int> Nstep_restr_splited;
    for(int jj=0;jj<n_splited;jj++) {
      int Nstep_one=Nstep_restr*WEIFHT_NSTEP[jj];
      Nstep_restr_splited.push_back(Nstep_one);
    }

    vector<string> Vrun_restr=Vrun;
    for(int jj=0;jj<n_splited;jj++) {
      string sline="";
      Vrun_restr.push_back(sline);
      for(size_t kk=0;kk<Vrestr_bond_coeff[jj].size();kk++) {
        sline="bond_coeff      "+Vrestr_bond_coeff[jj][kk];
        Vrun_restr.push_back(sline);
      }
      sline="";
      Vrun_restr.push_back(sline);
      sline="timestep        "+TIMESTEP[jj];
      Vrun_restr.push_back(sline);
      if(jj==0) {
        sline="run             100000 upto";
        Vrun_restr.push_back(sline);
      }
      if(md_label==PrepareLAMMPS::remd) {
        ostringstream oss;
        oss<<"temper          "<<Nstep_restr_splited[jj]<<" 5000 $t mylang "<<35824+(int)(165243*cbs->rand_one->uniform())<<" "<<86529+(int)(43280*cbs->rand_one->uniform());
        sline=oss.str();
        Vrun_restr.push_back(sline);
        oss.str("");
        oss.clear();
      }
      else if(md_label==PrepareLAMMPS::md) {
        sline="run             "+to_string(Nstep_restr_splited[jj]);
        Vrun_restr.push_back(sline);
      }
    }
    
    vector<string> Vin_scripts;
    if(md_label==PrepareLAMMPS::md) {
      int nthread=cbs->init->nthread;
      int nloop=nmodel/nthread;
      int nlast=nmodel%nthread;
      for(int jj=0;jj<nloop;jj++) {
        string str=to_string(jj*nthread+1);
        for(int kk=1;kk<nthread;kk++)
          str+=" "+to_string(jj*nthread+1+kk);
        prelmp_one.prepare_in_script(outDIR,sjobname,str,Nstep_restr,Nevery_restr);
        if(Vgroup.size()>0)
          prelmp_one.Vgroup=Vgroup;
        if(Vfix.size()>0)
          prelmp_one.Vfix=Vfix;
        prelmp_one.Vrun=Vrun_restr;
        prelmp_one.Vem=Vem_skip;
        fstr=outDIR+"in."+sjobname+"."+to_string(jj+1);
        if(flag_exe_restr)
          prelmp_one.write_in_script(fstr,1);
        string srun=fstr+" "+to_string(nthread)+" 1";
        Vin_scripts.push_back(srun);
      }
      if(nlast>0) {
        string str=to_string(nloop*nthread+1);
        for(int kk=1;kk<nlast;kk++)
          str+=" "+to_string(nloop*nthread+1+kk);
        prelmp_one.prepare_in_script(outDIR,sjobname,str,Nstep_restr,Nevery_restr);
        if(Vgroup.size()>0)
          prelmp_one.Vgroup=Vgroup;
        if(Vfix.size()>0)
          prelmp_one.Vfix=Vfix;
        prelmp_one.Vrun=Vrun_restr;
        prelmp_one.Vem=Vem_skip;
        fstr=outDIR+"in."+sjobname+"."+to_string(nloop+1);
        if(flag_exe_restr)
          prelmp_one.write_in_script(fstr,1);
        string srun=fstr+" "+to_string(nlast)+" 1";
        Vin_scripts.push_back(srun);
      }
    }
    else if(md_label==PrepareLAMMPS::remd) {
      for(int jj=1;jj<=nmodel;jj++) {
        string str=to_string(jj);
        prelmp_one.prepare_in_script(outDIR,sjobname,str,Nstep_restr,Nevery_restr);
        fstr=outDIR+"in."+sjobname+"."+to_string(jj);
        if(Vgroup.size()>0)
          prelmp_one.Vgroup=Vgroup;
        if(Vfix.size()>0)
          prelmp_one.Vfix=Vfix;
        prelmp_one.Vrun=Vrun_restr;
        prelmp_one.Vem=Vem_skip;
        if(flag_exe_restr)
          prelmp_one.write_in_script(fstr,1);
        string srun=fstr+" "+to_string(prelmp_one.Vtemperature.size())+" 1";
        Vin_scripts.push_back(srun);
      }
    }

    slog="Nstep_restr: "+to_string(Nstep_restr+100000);
    cbs->error->log(slog);

    string fsh=cbs->init->WORK_DIR+"run_lammps_"+sjobname+".sh";
    vector<string> Vcmd;
    if(flag_exe_restr)
      Vcmd=prelmp_one.write_run_script(Vin_scripts,fsh,flag_gpu);
    string flammpslog=outDIR+"run_lammps_"+jobname+".log";
    if(flag_exe_restr) {
      for(size_t jj=0;jj<Vcmd.size();jj++) {
        string str=" Running "+to_string(jj+1)+"/"+to_string(Vcmd.size())+" job...";
        string cmd="echo \""+str+"\" >> "+flammpslog;
        cbs->execute_system_command(cmd);
        cmd=Vcmd[jj]+" >> "+flammpslog;
        cbs->run_lammps(cmd);
      }
    }

    int step_begin=0;
    for(int jj=0;jj<n_splited-1;jj++)
      step_begin+=Nstep_restr_splited[jj];
    int step_end=-1;
    vector<string> Vtrjs_list=prelmp_one.Voutput_trjs;
    vector<string> Vlogs_list=prelmp_one.Voutput_logs;
    if(md_label==PrepareLAMMPS::md) {
      for(int ii=0;ii<nmodel;ii++) {
        ReadLAMMPS readlmp_one(cbs);
        readlmp_one.read_lammpstrj(Vtrjs_list[ii],step_begin,step_end);
        if(readlmp_one.VVsteps.size()==0 || readlmp_one.VVsteps[0].size()==0) {
          readlmp_one.clear();
          readlmp_one.read_lammpstrj(Vtrjs_list[ii],0,step_end);
          if(readlmp_one.VVsteps.size()==0 || readlmp_one.VVsteps[0].size()==0) continue;
          readlmp_one.read_log(Vlogs_list[ii],0,step_end);
        }
        else {
          readlmp_one.read_log(Vlogs_list[ii],step_begin,step_end);
        }
        vector<Atom> Vatoms_updated;
        if(readlmp_one.VVenergies.size()==0 || readlmp_one.VVenergies[0].size()==0) {
          if(readlmp_one.VVatoms.size()>0 && readlmp_one.VVatoms[0].size()>0) 
            Vatoms_updated=readlmp_one.VVatoms[0].back();
        }
        else {
          Vatoms_updated=readlmp_one.extract_lowest_energy_structure();
        }
        if(Vatoms_updated.size()>0) {
          vector<double> Vxyz;
          for(size_t jj=0;jj<Vatoms_updated.size();jj++) {
            vector<double> xyz=Vatoms_updated[jj].get_atom_xyz();
            Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
          }
          Vrnacg[ii].rna_cg.reset_molecule_xyz(Vxyz);
        }
      }
    }
    else if(md_label==PrepareLAMMPS::remd) {
      for(int ii=0;ii<nmodel;ii++) {
        vector<string> Vtrjs_list_i,Vlogs_list_i;
        for(size_t jj=0;jj<prelmp_one.Vtemperature.size();jj++) {
          size_t index_j=ii*prelmp_one.Vtemperature.size()+jj;
          Vtrjs_list_i.push_back(Vtrjs_list[index_j]);
          Vlogs_list_i.push_back(Vlogs_list[index_j]);
        }
        ReadLAMMPS readlmp_one(cbs);
        readlmp_one.read_lammpstrj(Vtrjs_list_i,step_begin,step_end);
        if(readlmp_one.VVsteps.size()==0 || readlmp_one.VVsteps[0].size()==0) {
          readlmp_one.clear();
          readlmp_one.read_lammpstrj(Vtrjs_list_i,0,step_end);
          if(readlmp_one.VVsteps.size()==0 || readlmp_one.VVsteps[0].size()==0) continue;
          readlmp_one.read_log(Vlogs_list_i,0,step_end);
        }
        else {
          readlmp_one.read_log(Vlogs_list_i,step_begin,step_end);
        }
        vector<Atom> Vatoms_updated;
        if(readlmp_one.VVenergies.size()==0 || readlmp_one.VVenergies[0].size()==0) {
          if(readlmp_one.VVatoms.size()>0 && readlmp_one.VVatoms[0].size()>0) 
            Vatoms_updated=readlmp_one.VVatoms[0].back();
        }
        else {
          Vatoms_updated=readlmp_one.extract_lowest_energy_structure();
        }
        if(Vatoms_updated.size()>0) {
          vector<double> Vxyz;
          for(size_t jj=0;jj<Vatoms_updated.size();jj++) {
            vector<double> xyz=Vatoms_updated[jj].get_atom_xyz();
            Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
          }
          Vrnacg[ii].rna_cg.reset_molecule_xyz(Vxyz);
        }
      }
    }
  }

  prelmp_one.Voutput_trjs.clear();
  prelmp_one.Voutput_logs.clear();
  string fname=jobname+"_prod";
  for(int ii=0;ii<nmodel;ii++) {
    if(ndistancerestraint>0)
      Vrnacg[ii].add_distance_restraint(readrestr_one.distance_restraints);
    if(nbp_restr>0) {
      if(cbs->init->RNA_version=="latest")
        Vrnacg[ii].add_RNA_basepairs(readrestr_one.RNA_basepairs);
      else if(cbs->init->RNA_version=="v2")
        Vrnacg[ii].add_RNA_basepairs_v2(readrestr_one.RNA_basepairs);
    }
    if(nbs_restr>0)
      Vrnacg[ii].add_RNA_basestacks(readrestr_one.RNA_basestacks);
    if(flag_ssrestr) {
      if(cbs->init->RNA_version=="latest")
        Vrnacg[ii].add_RNA_basepairs(rna_bp_one.mapBP);
      else if(cbs->init->RNA_version=="v2")
        Vrnacg[ii].add_RNA_basepairs_v2(rna_bp_one.mapBP);
    }
    vector<string> Vdata_one;
    Vdata_one=prelmp_one.prepare_data(Vrnacg[ii].rna_cg,Vrnacg[ii].bonds,Vrnacg[ii].angles,Vrnacg[ii].dihedrals,Vrnacg[ii].impropers);
    if(Vdata_one[0].find("Number of atom types")!=string::npos) {
      cout<<"Warning: Model #"<<ii+1<<" "<<Vdata_one[0]<<endl;
      Vdata_one.erase(Vdata_one.begin());
    }
    Vdata_one[0]+=inPDB+" Model #"+to_string(ii+1);
    fstr=outDIR+"data."+fname+"."+to_string(ii+1);
    outF.open(fstr.c_str());
    for(auto sline : Vdata_one)
      outF<<sline<<endl;
    outF.close(); outF.clear();
  }

  int Nevery=Nstep/Nout;
  vector<string> Vin_scripts;
  if(md_label==PrepareLAMMPS::md) {
    int nthread=cbs->init->nthread;
    int nloop=nmodel/nthread;
    int nlast=nmodel%nthread;
    for(int jj=0;jj<nloop;jj++) {
      string str=to_string(jj*nthread+1);
      for(int kk=1;kk<nthread;kk++)
        str+=" "+to_string(jj*nthread+1+kk);
      prelmp_one.prepare_in_script(outDIR,fname,str,Nstep,Nevery);
      if(Vgroup.size()>0)
        prelmp_one.Vgroup=Vgroup;
      if(Vfix.size()>0)
        prelmp_one.Vfix=Vfix;
      if(Vrun.size()>0) {
        vector<string> Vrun_full=Vrun;
        for(size_t kk=0;kk<prelmp_one.Vrun.size();kk++) {
          string sline=prelmp_one.Vrun[kk];
          if(sline.find("velocity")!=string::npos && sline.find("create")!=string::npos) continue;
          Vrun_full.push_back(sline);
        }
        prelmp_one.Vrun=Vrun_full;
      }
      if(Vrun_timestep.size()>0) {
        vector<string> Vrun_extend;
        size_t nline=prelmp_one.Vrun.size();
        for(size_t kk=0;kk<nline;kk++) {
          string sline=prelmp_one.Vrun[kk];
          if(kk<nline-3)
            Vrun_extend.push_back(sline);
          else {
            if(sline.find("run_style")!=string::npos) 
              Vrun_extend.push_back(sline);
          }
        }
        for(size_t kk=0;kk<Vrun_timestep.size();kk++)
          Vrun_extend.push_back(Vrun_timestep[kk]);
        prelmp_one.Vrun=Vrun_extend;
      }
      prelmp_one.Vem=Vem_skip;
      fstr=outDIR+"in."+fname+"."+to_string(jj+1);
      prelmp_one.write_in_script(fstr,1);
      string srun=fstr+" "+to_string(nthread)+" 1";
      Vin_scripts.push_back(srun);
    }
    if(nlast>0) {
      string str=to_string(nloop*nthread+1);
      for(int kk=1;kk<nlast;kk++)
        str+=" "+to_string(nloop*nthread+1+kk);
      prelmp_one.prepare_in_script(outDIR,fname,str,Nstep,Nevery);
      if(Vgroup.size()>0)
        prelmp_one.Vgroup=Vgroup;
      if(Vfix.size()>0)
        prelmp_one.Vfix=Vfix;
      if(Vrun.size()>0) {
        vector<string> Vrun_full=Vrun;
        for(size_t kk=0;kk<prelmp_one.Vrun.size();kk++) {
          string sline=prelmp_one.Vrun[kk];
          if(sline.find("velocity")!=string::npos && sline.find("create")!=string::npos) continue;
          Vrun_full.push_back(sline);
        }
        prelmp_one.Vrun=Vrun_full;
      }
      if(Vrun_timestep.size()>0) {
        vector<string> Vrun_extend;
        size_t nline=prelmp_one.Vrun.size();
        for(size_t kk=0;kk<nline;kk++) {
          string sline=prelmp_one.Vrun[kk];
          if(kk<nline-3)
            Vrun_extend.push_back(sline);
          else {
            if(sline.find("run_style")!=string::npos) 
              Vrun_extend.push_back(sline);
          }
        }
        for(size_t kk=0;kk<Vrun_timestep.size();kk++)
          Vrun_extend.push_back(Vrun_timestep[kk]);
        prelmp_one.Vrun=Vrun_extend;
      }
      prelmp_one.Vem=Vem_skip;
      fstr=outDIR+"in."+fname+"."+to_string(nloop+1);
      prelmp_one.write_in_script(fstr,1);
      string srun=fstr+" "+to_string(nlast)+" 1";
      Vin_scripts.push_back(srun);
    }
  }
  else if(md_label==PrepareLAMMPS::remd) {
    int nm=cbs->init->nthread / prelmp_one.Vtemperature.size();
    for(int jj=1;jj<=nmodel;jj++) {
      string str=to_string(jj);
      prelmp_one.prepare_in_script(outDIR,fname,str,Nstep,Nevery);
      if(Vgroup.size()>0)
        prelmp_one.Vgroup=Vgroup;
      if(Vfix.size()>0)
        prelmp_one.Vfix=Vfix;
      if(Vrun.size()>0) {
        vector<string> Vrun_full=Vrun;
        for(size_t kk=0;kk<prelmp_one.Vrun.size();kk++) {
          string sline=prelmp_one.Vrun[kk];
          if(sline.find("velocity")!=string::npos && sline.find("create")!=string::npos) continue;
          Vrun_full.push_back(sline);
        }
        prelmp_one.Vrun=Vrun_full;
      }
      if(Vrun_timestep.size()>0) {
        vector<string> Vrun_extend;
        size_t nline=prelmp_one.Vrun.size();
        for(size_t kk=0;kk<nline;kk++) {
          string sline=prelmp_one.Vrun[kk];
          if(kk<nline-3)
            Vrun_extend.push_back(sline);
          else {
            if(sline.find("run_style")!=string::npos) 
              Vrun_extend.push_back(sline);
          }
        }
        for(size_t kk=0;kk<Vrun_timestep.size();kk++)
          Vrun_extend.push_back(Vrun_timestep[kk]);
        prelmp_one.Vrun=Vrun_extend;
      }
      prelmp_one.Vem=Vem_skip;
      string fstr=outDIR+"in."+fname+"."+to_string(jj);
      prelmp_one.write_in_script(fstr,1);
      string srun=fstr+" "+to_string(prelmp_one.Vtemperature.size())+" "+to_string(nm);
      Vin_scripts.push_back(srun);
    }
  }

  fstr=outDIR+fname+"_lammpstrj_list.dat";
  outF.open(fstr.c_str());
  for(auto sline : prelmp_one.Voutput_trjs)
    outF<<sline<<endl;
  outF.close();outF.clear();

  fstr=outDIR+fname+"_lammpslog_list.dat";
  outF.open(fstr.c_str());
  for(auto sline : prelmp_one.Voutput_logs)
    outF<<sline<<endl;
  outF.close(); outF.clear();
  
  slog="Nstep_prod: "+to_string(Nstep+50000);
  cbs->error->log(slog);

  string fsh=cbs->init->WORK_DIR+"run_lammps_"+fname+".sh";
  vector<string> Vcmd=prelmp_one.write_run_script(Vin_scripts,fsh,flag_gpu);
  string flammpslog=outDIR+"run_lammps_"+jobname+".log";
  if(flag_exe) {
    for(size_t jj=0;jj<Vcmd.size();jj++) {
      string str=" Running "+to_string(jj+1)+"/"+to_string(Vcmd.size())+" job...";
      string cmd="echo \""+str+"\" >> "+flammpslog;
      cbs->execute_system_command(cmd);
      cmd=Vcmd[jj]+" >> "+flammpslog;
      cbs->run_lammps(cmd);
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