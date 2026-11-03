#include "basic.h"
#include "error.h"
#include "readfile.h"
#include "initialization.h"
#include "struct_definition.h"
#include "cluster.h"
#include "measure.h"
#include "superpose.h"
#include "RNA_molecule.h"
#include "RNA_motif.h"
#include "RNA_basepair.h"
#include "RNA_coarse_grained.h"
#include "RNA_recovery.h"
#include "read_RNA_template.h"
#include "read_restraint.h"
#include "parse_RNA_motif.h"
#include "match_RNA_motif.h"
#include "find_RNA_template.h"
#include "assemble_RNA_motif.h"
#include "build_RNA_model.h"
#include "read_amber.h"
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

  string smotif="";
  if(Mcmd.count("-im")) {
    string fim=Mcmd["-im"][0];
    Mcmd.erase("-im");
    vector<string> Vlines=read_file(fim);
    smotif=Vlines[0];
  }

  if(seq.length()==0 && smotif.length()==0) {
    cout<<"Wrong commands: No input (-i or -seq or -im)!"<<endl;
    exit(0);
  }
  
  string jobname="rnamotif";
  if(Mcmd.count("-name")) {
    jobname=Mcmd["-name"][0];
    Mcmd.erase("-name");
  }
  string outDIR="./prediction/";
  if(Mcmd.count("-DIR")) {
    outDIR=Mcmd["-DIR"][0];
    Mcmd.erase("-DIR");
  }
  if(outDIR.size()!=0 && outDIR.back()!='/')
    outDIR+="/";
 
  int Nstep=-1;
  if(Mcmd.count("-nstep")) {
    Nstep=stoi(Mcmd["-nstep"][0]);
    Mcmd.erase("-nstep");
  }

  int nmax_template=-1;
  if(Mcmd.count("-nmax_template")) {
    nmax_template=stoi(Mcmd["-nmax_template"][0]);
    Mcmd.erase("-nmax_template");
  }
  int flag_gpu=0;
  if(Mcmd.count("-gpu")) {
    flag_gpu=1;
    Mcmd.erase("-gpu");
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

  string slog="##### Build_RNA_motif_structure #####";
  cbs->error->log(slog);

  if(!is_directory_exist(outDIR)) {
    string cmd="mkdir "+outDIR;
    cbs->execute_system_command(cmd);
  }

  cbs->init->WORK_DIR=outDIR;

  if(nmax_template>0)
    cbs->init->RNA_template_nmax=nmax_template;

  vector<RNAMotif> Vmotifs;
  if(seq.length()>0) {
    vector<string> Vseq=split_string(seq,"-&");
    vector<string> Vsecstr=split_string(secstr,"-&");
    if(Vseq.size()!=Vsecstr.size()) {
      string str="Input sequence and 2D structure are not match! --> "+seq+" "+secstr;
      cbs->error->one(str);
    }
    size_t napp=0;
    for(size_t ii=0;ii<Vseq.size();ii++) {
      if(Vseq[ii].size()!=Vsecstr[ii].size()) {
        string str="Input sequence and 2D structure are not match! --> "+seq+" "+secstr;
        cbs->error->one(str);
      }
      Vseq[ii]=to_string((char)Vseq[ii][0])+Vseq[ii]+to_string((char)Vseq[ii].back());
      Vsecstr[ii]=to_string((char)Vsecstr[ii][0])+Vsecstr[ii]+to_string((char)Vsecstr[ii].back());
      napp+=2;
    }

    string seq_app=Vseq[0];
    for(size_t jj=1;jj<Vseq.size();jj++)
      seq_app+="&"+Vseq[jj];
    string secstr_app=Vsecstr[0];
    for(size_t jj=1;jj<Vsecstr.size();jj++)
      secstr_app+="&"+Vsecstr[jj];

    RNAMolecule rna_one(cbs);
    rna_one.set_RNA_molecule(seq_app);
    vector<vector<string> > Vres=rna_one.get_residue_name();
    RNABasePair rna_bp_one(cbs);
    rna_bp_one.parse_RNA_basepairs(Vres,secstr_app);
    ParseRNAMotif PRmotif_one(cbs);
    PRmotif_one.parse_RNA_motifs(rna_bp_one);
    string topo=PRmotif_one.get_RNA_topology_string();
    for(size_t jj=0;jj<PRmotif_one.motifs.size();jj++) {
      RNAMotif motif_one=PRmotif_one.motifs[jj];
      if(motif_one.type.find("Helix")!=string::npos || motif_one.type.find("Tertiary_helix")!=string::npos) continue;
      Vmotifs.push_back(motif_one);
    }
  }
  else if(seq.length()==0 && smotif.length()>0) {
    RNAMotif motif_one(cbs);
    motif_one.set_RNA_motif(smotif);
    if(motif_one.type.find("way_junction")!=string::npos) {
      RNAMotif Tmotif_one(cbs);
      Tmotif_one=motif_one;
      motif_one.clear();
      motif_one.type=Tmotif_one.type;
      int len_min=Tmotif_one.lengths[0];
      size_t index_min=0;
      for(size_t kk=1;kk<Tmotif_one.lengths.size();kk++) {
        if(Tmotif_one.lengths[kk]<len_min) {
          len_min=Tmotif_one.lengths[kk];
          index_min=kk;
        }
      }
      for(size_t kk=index_min;kk<Tmotif_one.lengths.size();kk++) {
        motif_one.insert_RNA_segment(Tmotif_one.chains[kk],-2);
        motif_one.secstrs.push_back(Tmotif_one.secstrs[kk]);
      }
      for(size_t kk=0;kk<index_min;kk++) {
        motif_one.insert_RNA_segment(Tmotif_one.chains[kk],-2);
        motif_one.secstrs.push_back(Tmotif_one.secstrs[kk]);
      }
    }
    Vmotifs.push_back(motif_one);
  }

  slog=Vmotifs[0].get_RNA_motif_info();
  cbs->error->log(slog);

  ReadRestraint readrestr_one(cbs);
  string frestr=outDIR+jobname+"_restraint.dat";
  
  if(cbs->init->Excluded_PDBs.size()>0) {
    (*cbs->logfile)<<"Excluded PDBs:";
    for(size_t ii=0;ii < cbs->init->Excluded_PDBs.size();ii++)
      (*cbs->logfile)<<" "<<cbs->init->Excluded_PDBs[ii];
    (*cbs->logfile)<<endl<<endl;
  }

  if(Vmotifs[0].type.find("-way_junction")!=string::npos) {  //add Lone_BP to close loops when failed to find template
    for(size_t jj=0;jj<Vmotifs[0].chains.size();jj++) {
      Residue residue_last=Vmotifs[0].chains[jj].residues.back();
      string secstr_last=to_string((char)Vmotifs[0].secstrs[jj].back());
      Residue residue_first_next(cbs);
      string secstr_first_next;
      if(jj<Vmotifs[0].chains.size()-1) {
        residue_first_next=Vmotifs[0].chains[jj+1].residues[0];
        secstr_first_next=to_string((char)Vmotifs[0].secstrs[jj+1][0]);
      }
      else {
        residue_first_next=Vmotifs[0].chains[0].residues[0];
        secstr_first_next=to_string((char)Vmotifs[0].secstrs[0][0]);
      }
      Chain Lbp_1(cbs),Lbp_2(cbs);
      Lbp_1.chain_label=Chain::rna;
      Lbp_2.chain_label=Chain::rna;
      Lbp_1.insert_residue(residue_last,Chain::begin);
      Lbp_2.insert_residue(residue_first_next,Chain::begin);
      RNAMotif Lbp(cbs);
      Lbp.type="Lone_BP";
      Lbp.insert_RNA_segment(Lbp_1,0);
      Lbp.insert_RNA_segment(Lbp_2,0);
      Lbp.secstrs.push_back(secstr_last);
      Lbp.secstrs.push_back(secstr_first_next);
      Vmotifs.push_back(Lbp);
    }
  }
    
  BuildRNAModel BRmodel_one(cbs);
  double moved_distance_max=BRmodel_one.build_RNA_model(Vmotifs,readrestr_one);
  (*cbs->logfile)<<endl;
  (*cbs->logfile)<<BRmodel_one;
  (*cbs->logfile)<<"Build "<<BRmodel_one.models.size()<<" initial_structure."<<endl;
  (*cbs->logfile)<<endl;
  string fstr=outDIR+jobname+"_build.pdb";
  ofstream outF (fstr.c_str());
  for(size_t ii=0;ii<BRmodel_one.models.size();ii++) {
    outF<<"MODEL     "<<ii+1<<endl;
    for(size_t jj=0;jj<BRmodel_one.models[ii].size();jj++)
      outF<<BRmodel_one.models[ii][jj];
    outF<<"ENDMDL"<<endl;
  }
  outF<<"END"<<endl;
  outF.close(); outF.clear();

  outF.open(frestr.c_str());
  outF<<readrestr_one;
  outF.close();outF.clear();

  ReadCoeff readcoeff_one(cbs);
  readcoeff_one.read_RNA_coeff();
  readcoeff_one.reset_paircoeff_lj();
  readcoeff_one.get_interaction_style();
  PrepareLAMMPS prelmp_one(cbs);
  double Temp=300.0;
  prelmp_one.set_temperature(Temp,Temp);
  prelmp_one.set_coeff(readcoeff_one);
  
  double moved_distance_max_limited=moved_distance_max;
  if(Nstep<0) {
    Nstep=1000000+round(10*moved_distance_max)*100000;
    if(Nstep>50000000) Nstep=50000000;
    if(moved_distance_max_limited > 49.0)
      moved_distance_max_limited=49.0;
  }

  string restr_type="harmonic 0.001 "+to_string(moved_distance_max/4);
  string FLEXIBLE_LOOP[]={"Loop","Loop_kiss_loop","Unpaired_chain","Long_loop"};
  int flag_loop=0;
  int flag_motif=0;
  int nmodel_total=0;
  int nbond_restr=0;
  vector<vector<Atom> > VVatoms_cg;
  for(size_t ii=0;ii<BRmodel_one.models.size();ii++) {
    vector<Atom> Vatoms_cg, Vatoms_cg_modified;
    vector<Bond> Vbonds;
    vector<Angle> Vangles;
    vector<Dihedral> Vdihedrals, Vimpropers;
    int nmotif=1;
    size_t index_offset=0;
    for(size_t jj=0;jj<BRmodel_one.models[ii].size();jj++) {
      RNAMotif &mtf_one=BRmodel_one.models[ii][jj];
      Molecule mol_one=mtf_one;
      RNACoarseGrained rnacg_one(cbs);
      if(cbs->init->RNA_version=="latest")
        rnacg_one.convert_RNA_to_CG(mol_one,index_offset);
      else if(cbs->init->RNA_version=="v2")
        rnacg_one.convert_RNA_to_CG_v2(mol_one,index_offset);
      index_offset+=rnacg_one.number_of_CG_atoms();
      vector<Atom> Vatoms_list=rnacg_one.rna_cg.get_atom_list();
      Vatoms_cg.insert(Vatoms_cg.end(),Vatoms_list.begin(),Vatoms_list.end());
      if(count(begin(FLEXIBLE_LOOP),end(FLEXIBLE_LOOP),mtf_one.type)) {
        if(cbs->init->RNA_version=="latest")
          rnacg_one.build_RNA_topology();
        else if(cbs->init->RNA_version=="v2")
          rnacg_one.build_RNA_topology_v2();
        Vbonds.insert(Vbonds.end(),rnacg_one.bonds.begin(),rnacg_one.bonds.end());
        Vangles.insert(Vangles.end(),rnacg_one.angles.begin(),rnacg_one.angles.end());
        Vdihedrals.insert(Vdihedrals.end(),rnacg_one.dihedrals.begin(),rnacg_one.dihedrals.end());
        Vimpropers.insert(Vimpropers.end(),rnacg_one.impropers.begin(),rnacg_one.impropers.end());
        flag_loop=1;
        for(size_t kk=0;kk<Vatoms_list.size();kk++) {
          Atom atom_one=Vatoms_list[kk];
          atom_one.chainid="0";
          Vatoms_cg_modified.push_back(atom_one);
        }
      }
      else {
        for(size_t kk=0;kk<Vatoms_list.size();kk++) {
          Atom atom_one=Vatoms_list[kk];
          atom_one.chainid=to_string(nmotif);
          Vatoms_cg_modified.push_back(atom_one);
        }
        nmotif++;
        flag_motif=1;
      }
    }

    VVatoms_cg.push_back(Vatoms_cg);
    nmodel_total++;

    vector<pair<int,string> > Vatomid_name;
    for(size_t jj=0;jj<Vatoms_cg.size();jj++) {
      int atomid_one=Vatoms_cg[jj].id;
      string atomname_full=Vatoms_cg[jj].get_atom_name();
      pair<int,string> p_one=make_pair(atomid_one,atomname_full);
      Vatomid_name.push_back(p_one);
    }
    for(size_t jj=0;jj<Vatomid_name.size()-1;jj++) {
      int atomid_i=Vatomid_name[jj].first;
      string atomname_i=Vatomid_name[jj].second;
      for(size_t kk=jj+1;kk<Vatomid_name.size();kk++) {
        if(atomname_i==Vatomid_name[kk].second) {
          Bond bond_one;
          bond_one.atom_i=atomid_i;
          bond_one.atom_j=Vatomid_name[kk].first;
          bond_one.type=restr_type;
          Vbonds.push_back(bond_one);
          nbond_restr++;
        }
      }
    }

    vector<string> Vdata=prelmp_one.prepare_data(Vatoms_cg_modified,Vbonds,Vangles,Vdihedrals,Vimpropers);
    if(Vdata[0].find("Number of atom types")!=string::npos) {
      cbs->error->warning(Vdata[0]);
      Vdata.erase(Vdata.begin());
    }
    fstr=outDIR+"data."+jobname+"_build."+to_string(nmodel_total);
    outF.open(fstr.c_str());
    for(size_t kk=0;kk<Vdata.size();kk++)
      outF<<Vdata[kk]<<endl;
    outF.close(); outF.clear();

    fstr=outDIR+jobname+"_build_cg."+to_string(nmodel_total)+".pdb";
    outF.open(fstr.c_str());
    for(size_t kk=0;kk<Vatoms_cg_modified.size();kk++) {
      Atom atom_one=Vatoms_cg_modified[kk];
      int chainid_index=stoi(atom_one.chainid);
      chainid_index=chainid_index%cbs->init->CHAIN_ID_LIST.size();
      atom_one.chainid=cbs->init->CHAIN_ID_LIST[chainid_index];
      outF<<atom_one;
    }
    outF.close(); outF.clear();
  }
  
  if(nbond_restr==0) Nstep=0;
  int Nevery=Nstep/100;
  if(Nevery<1) Nevery=1;
  double WEIGHT_NSTEP[]={0.1,0.1,0.3,0.3,0.2};
  int n_splited=get_array_size(WEIGHT_NSTEP);
  vector<int> Nstep_splited;
  for(int jj=0;jj<n_splited;jj++) {
    int Nstep_one=Nstep*WEIGHT_NSTEP[jj];
    Nstep_splited.push_back(Nstep_one);
  }
  string KBOND[]={"0.001","0.01","0.1","1.0","10.0"};
  double B0[]={moved_distance_max_limited/4,moved_distance_max_limited/8,moved_distance_max_limited/16,0.0,0.0};
  string TIMESTEP[]={"0.2","0.5","0.75","1.0","1.0"};
  string restr_bond_label="";
  if(nbond_restr>0 && prelmp_one.Mrestraint_type.count(restr_type)) {
    restr_bond_label=prelmp_one.Mrestraint_type[restr_type];
  } 

  vector<string> Vgroup;
  string sline;
  if(flag_loop) {
    string sline="group           loop molecule 0";
    Vgroup.push_back(sline);
  }
  if(flag_motif) {
    string sline="group           motif molecule > 0";
    Vgroup.push_back(sline);
  }

  vector<string> Vfix;
  ostringstream oss;
  if(flag_motif) {
    oss<<"fix             1 motif rigid molecule langevin "<<std::fixed<<std::setprecision(2)<<Temp<<" ";
    oss<<std::fixed<<std::setprecision(2)<<Temp<<" 100.0 "<<83546+(int)(627043*cbs->rand_one->uniform());
    sline=oss.str();
    Vfix.push_back(sline);
  }
  if(flag_loop) {
    oss.clear();
    oss.str("");
    oss<<"fix             2 loop langevin "<<std::fixed<<std::setprecision(2)<<Temp<<" ";
    oss<<std::fixed<<std::setprecision(2)<<Temp<<" 100.0 "<<49073+(int)(937625*cbs->rand_one->uniform());
    sline=oss.str();
    Vfix.push_back(sline);
    sline="fix             3 loop nve/limit 0.5";
    Vfix.push_back(sline);
  }

  vector<string> Vrun;
  if(flag_motif) {
    oss.clear();
    oss.str("");
    oss<<"velocity        motif create "<<std::fixed<<std::setprecision(2)<<Temp<<" ";
    oss<<54627+(int)(913240*cbs->rand_one->uniform())<<" rot yes dist gaussian";
    sline=oss.str();
    Vrun.push_back(sline);
  }
  if(flag_loop) {
    oss.clear();
    oss.str("");
    oss<<"velocity        loop create "<<std::fixed<<std::setprecision(2)<<Temp<<" ";
    oss<<29748+(int)(821507*cbs->rand_one->uniform())<<" rot yes dist gaussian";
    sline=oss.str();
    Vrun.push_back(sline);
  }
  sline="timestep        "+TIMESTEP[0];
  Vrun.push_back(sline);
  if(Nstep_splited[0]>0) {
    sline="run             "+to_string(Nstep_splited[0])+" upto";
    Vrun.push_back(sline);
  }
  if(nbond_restr>0) {
    for(int kk=1;kk<n_splited;kk++) {
      sline="";
      Vrun.push_back(sline);
      if(restr_bond_label.length()>0) {
        if(prelmp_one.Vbondstyle.size()==1) {
          sline="bond_coeff      "+restr_bond_label+" "+KBOND[kk]+" "+to_string(B0[kk]);
        }
        else {
          sline="bond_coeff      "+restr_bond_label+" harmonic "+KBOND[kk]+" "+to_string(B0[kk]);
        }
        Vrun.push_back(sline);
      }
      sline="";
      Vrun.push_back(sline);
      sline="timestep        "+TIMESTEP[kk];
      Vrun.push_back(sline);
      sline="run             "+to_string(Nstep_splited[kk]);
      Vrun.push_back(sline);
    }
  }

  string fname=jobname+"_build";
  vector<string> Vrun_files;
  int nthread=cbs->init->nthread;
  int nloop=nmodel_total/nthread;
  int nlast=nmodel_total%nthread;
  for(int jj=0;jj<nloop;jj++) {
    string str=to_string(jj*nthread+1);
    for(int kk=1;kk<nthread;kk++)
      str+=" "+to_string(jj*nthread+1+kk);
    prelmp_one.prepare_in_script(outDIR,fname,str,Nstep_splited[0],Nevery);
    prelmp_one.Vgroup=Vgroup;
    prelmp_one.Vfix=Vfix;
    prelmp_one.Vrun=Vrun;
    string fstr=outDIR+"in."+fname+"."+to_string(jj+1);
    prelmp_one.write_in_script(fstr);
    string srun=fstr+" "+to_string(nthread)+" 1";
    Vrun_files.push_back(srun);
  }
  if(nlast>0) {
    string str=to_string(nloop*nthread+1);
    for(int kk=1;kk<nlast;kk++)
      str+=" "+to_string(nloop*nthread+1+kk);
    prelmp_one.prepare_in_script(outDIR,fname,str,Nstep_splited[0],Nevery);
    prelmp_one.Vgroup=Vgroup;
    prelmp_one.Vfix=Vfix;
    prelmp_one.Vrun=Vrun;
    string fstr=outDIR+"in."+fname+"."+to_string(nloop+1);
    prelmp_one.write_in_script(fstr);
    string srun=fstr+" "+to_string(nlast)+" 1";
    Vrun_files.push_back(srun);
  }

  slog="Nstep_build: "+to_string(Nstep);
  cbs->error->log(slog);

  string fsh=cbs->init->WORK_DIR+"run_lammps_"+fname+".sh";
  vector<string> Vcmd=prelmp_one.write_run_script(Vrun_files,fsh,flag_gpu);
  string flammpslog=outDIR+"run_lammps_"+jobname+".log";
  for(size_t jj=0;jj<Vcmd.size();jj++) {
    string str=" Running "+to_string(jj+1)+"/"+to_string(Vcmd.size())+" job...";
    string cmd="echo \""+str+"\" >> "+flammpslog;
    cbs->execute_system_command(cmd);
    cmd=Vcmd[jj]+" -sc none >> "+flammpslog;
    cbs->run_lammps(cmd);
  }

  int step_begin=Nstep*(n_splited-1)/n_splited;
  int step_end=-1;
  vector<string> Vtrjs_list=prelmp_one.Voutput_trjs;
  ReadLAMMPS readlmp_one(cbs);
  readlmp_one.read_lammpstrj(Vtrjs_list,step_begin,step_end);
  for(size_t ii=0;ii<VVatoms_cg.size();ii++) {
    vector<Atom> Vatoms_xyz=readlmp_one.VVatoms[ii].back();
    for(size_t jj=0;jj<VVatoms_cg[ii].size();jj++) {
      VVatoms_cg[ii][jj].set_atom_xyz(Vatoms_xyz[jj]);
    }
  }

  vector<vector<Atom> > VVatoms_cg_sorted;
  vector<string> Vresnames=Vmotifs[0].Vresnames;
  for(size_t ii=0;ii<VVatoms_cg.size();ii++) {
    vector<Atom> Vatoms_sorted;
    vector<string> Vatomnames_full;
    for(size_t jj=0;jj<Vresnames.size();jj++) {
      for(size_t kk=0;kk<VVatoms_cg[ii].size();kk++) {
        Atom atom_one=VVatoms_cg[ii][kk];
        atom_one.resname=atom_one.resname.substr(1);
        if(atom_one.get_residue_name()==Vresnames[jj]) {
          string atomname_full=atom_one.get_atom_name();
          if(count(Vatomnames_full.begin(),Vatomnames_full.end(),atomname_full)==0) {
            Vatoms_sorted.push_back(VVatoms_cg[ii][kk]);
            Vatomnames_full.push_back(atomname_full);
          }
        }
      }
    }
    for(size_t kk=0;kk<Vatoms_sorted.size();kk++) {
      Vatoms_sorted[kk].id=static_cast<int>(kk+1);
    }
    VVatoms_cg_sorted.push_back(Vatoms_sorted);
  }

  size_t nres=Vmotifs[0].number_of_residues();
  if((int)VVatoms_cg_sorted.size()> cbs->init->Num_model_max) {
    double rmsd_cutoff=0.1*nres;
    double Pclu=0.5;
    Cluster clu_one(cbs);
    clu_one.input_data(VVatoms_cg_sorted);
    clu_one.calculate_RMSD_matrix();
    clu_one.get_clusters(rmsd_cutoff,Pclu);
    
    (*cbs->logfile)<<endl;
    (*cbs->logfile)<<clu_one;
    (*cbs->logfile)<<endl;

    vector<vector<Atom> > VVtmp=VVatoms_cg_sorted;
    VVatoms_cg_sorted.clear();
    int nclu=(int)clu_one.clusters.size();
    if(nclu > cbs->init->Num_model_max)
      nclu=cbs->init->Num_model_max;
    for(int jj=0;jj<nclu;jj++) {
      size_t index_clu=clu_one.clusters[jj][0];
      VVatoms_cg_sorted.push_back(VVtmp[index_clu]);
    }
  }

  fstr=outDIR+fname+"_cg_init.pdb";
  outF.open(fstr.c_str());
  for(size_t ii=0;ii<VVatoms_cg_sorted.size();ii++) {
    outF<<"MODEL       "<<(ii+1)<<endl;
    for(size_t jj=0;jj<VVatoms_cg_sorted[ii].size();jj++)
      outF<<VVatoms_cg_sorted[ii][jj];
    outF<<"ENDMDL"<<endl;
  }
  outF<<"END"<<endl;
  outF.close(); outF.clear();

  vector<Molecule> Vmol_cg;
  for(size_t ii=0;ii<VVatoms_cg_sorted.size();ii++) {
    Molecule mol_one(cbs);
    mol_one.set_molecule(VVatoms_cg_sorted[ii],Chain::cg);
    mol_one.set_atomid(0);
    Vmol_cg.push_back(mol_one);
  }

  prelmp_one.clear();
  readcoeff_one.clear();
  readcoeff_one.read_RNA_coeff();
  readcoeff_one.get_interaction_style();
  prelmp_one.set_coeff(readcoeff_one);
  string fmin=fname+"_min_cg";
  for(size_t ii=0;ii<Vmol_cg.size();ii++) {
    RNACoarseGrained rnacg_one(cbs);
    rnacg_one.rna_cg=Vmol_cg[ii];
    if(cbs->init->RNA_version=="latest")
      rnacg_one.build_RNA_topology();
    else if(cbs->init->RNA_version=="v2")
      rnacg_one.build_RNA_topology_v2();
    vector<string> Vdata=prelmp_one.prepare_data(rnacg_one.rna_cg,rnacg_one.bonds,rnacg_one.angles,rnacg_one.dihedrals,rnacg_one.impropers);
    fstr=outDIR+"data."+fmin+"."+to_string(ii+1);
    outF.open(fstr.c_str());
    for(auto sline : Vdata)
      outF<<sline<<endl;
    outF.close(); outF.clear();
  }
  prelmp_one.md_label=PrepareLAMMPS::em;
  Vrun_files.clear();
  Vrun_files=prelmp_one.write_in_script(outDIR,fmin,Vmol_cg.size());
  Vtrjs_list.clear();
  Vtrjs_list=prelmp_one.Voutput_trjs;

  fsh=cbs->init->WORK_DIR+"run_lammps_"+fmin+".sh";
  Vcmd.clear();
  Vcmd=prelmp_one.write_run_script(Vrun_files,fsh);
  for(size_t jj=0;jj<Vcmd.size();jj++) {
    string str=" Running "+to_string(jj+1)+"/"+to_string(Vcmd.size())+" job...";
    string cmd="echo \""+str+"\" >> "+flammpslog;
    cbs->execute_system_command(cmd);
    cmd=Vcmd[jj]+" -sc none >> "+flammpslog;
    cbs->run_lammps(cmd);
  }

  readlmp_one.clear();
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

  fstr=outDIR+fname+"_aa_init.pdb";
  outF.open(fstr.c_str());
  for(size_t ii=0;ii<Vmol_aa.size();ii++) {
    outF<<"MODEL     "<<ii+1<<endl;
    outF<<Vmol_aa[ii];
    outF<<"ENDMDL"<<endl;
  }
  outF<<"END"<<endl;
  outF.close(); outF.clear();

  time (&tend);

  double dtime=difftime(tend,tstart);
  string stime=cbs->calculate_runtime(dtime);
  cbs->error->log(stime);
  stime="";
  cbs->error->log(stime);

  delete cbs;
  
  return 0;
}