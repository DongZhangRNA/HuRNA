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
  if(seq.length()==0) {
    cout<<"Wrong commands: No input sequence (-i or -seq)!"<<endl;
    exit(0);
  }
  
  string jobname="rna";
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
  
  string frestr="";
  if(Mcmd.count("-restr")) {
    frestr=Mcmd["-restr"][0];
    Mcmd.erase("-restr");
  }
 
  int Nstep=-1;
  if(Mcmd.count("-nstep")) {
    Nstep=stoi(Mcmd["-step"][0]);
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
  int flag_RNA_use_DM=0;
  if(Mcmd.count("-use_DM")) {
    flag_RNA_use_DM=1;
    Mcmd.erase("-use_DM");
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

  string slog="##### Construct_RNA_initial_structure #####";
  cbs->error->log(slog);

  if(!is_directory_exist(outDIR)) {
    string cmd="mkdir "+outDIR;
    cbs->execute_system_command(cmd);
  }

  cbs->init->WORK_DIR=outDIR;
  cbs->init->load_RNA_DM(flag_RNA_use_DM);

  if(nmax_template>0)
    cbs->init->RNA_template_nmax=nmax_template;

  RNAMolecule rna_one(cbs);
  rna_one.set_RNA_molecule(seq);
  size_t nres=rna_one.number_of_residues();
  if(secstr.length()==0) {
    string str="No input 2D structure (-i or -secstr)!";
    cbs->error->warning(str);
    for(size_t jj=0;jj<nres;jj++)
      secstr+=".";
  }
  vector<vector<string> > Vres=rna_one.get_residue_name();
  RNABasePair rna_bp_one(cbs);
  rna_bp_one.parse_RNA_basepairs(Vres,secstr);
  ParseRNAMotif PRmotif_one(cbs);
  PRmotif_one.parse_RNA_motifs(rna_bp_one);
  string topo=PRmotif_one.get_RNA_topology_string();

  (*cbs->logfile)<<"Seq:"<<rna_one.get_sequences()<<endl;
  (*cbs->logfile)<<" 2D:"<<rna_bp_one.secstr_out<<endl;
  (*cbs->logfile)<<"Number_of_nucleotides: "<<nres<<endl;
  (*cbs->logfile)<<"Number_of_basepairs: "<<rna_bp_one.size()<<endl;
  (*cbs->logfile)<<"Topology: "<<topo<<endl;
  if(cbs->init->CircRNA_chains.size()>0) {
    string str="Circular_RNA_chain:";
    for(size_t jj=0;jj < cbs->init->CircRNA_chains.size();jj++)
      str+=" "+cbs->init->CircRNA_chains[jj];
    (*cbs->logfile)<<str<<endl;
  }
  (*cbs->logfile)<<endl;
  (*cbs->logfile)<<PRmotif_one;
  (*cbs->logfile)<<endl;

  ReadRestraint readrestr_one(cbs);
  if(frestr.length()>1) {
    readrestr_one.read_restraint(frestr);

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

    (*cbs->logfile)<<readrestr_one;
  }
  frestr=outDIR+jobname+"_restraint.dat";

  if(readrestr_one.RNA_templates.size()>0) {
    vector<string> Vresnames_UT;
    for(size_t ii=0;ii<readrestr_one.RNA_templates.size();ii++) {
      vector<string> Vres=readrestr_one.RNA_templates[ii].Vresnames;
      Vresnames_UT.insert(Vresnames_UT.end(),Vres.begin(),Vres.end());
    }

    vector<string> Vresnames_model=rna_one.Vresnames;
    string sname="";
    int flag_check=0;
    for(size_t ii=0;ii<Vresnames_UT.size();ii++) {
      if(count(Vresnames_model.begin(),Vresnames_model.end(),Vresnames_UT[ii])) continue;
      else {
        flag_check=1;
        sname+=" "+Vresnames_UT[ii];
      }
    }
    if(flag_check) {
      string str="RNA_template residues not existed! --> "+sname;
      cbs->error->one(str);
    }
  }
  

  if(cbs->init->Excluded_PDBs.size()>0) {
    (*cbs->logfile)<<"Excluded PDBs:";
    for(size_t ii=0;ii < cbs->init->Excluded_PDBs.size();ii++)
      (*cbs->logfile)<<" "<<cbs->init->Excluded_PDBs[ii];
    (*cbs->logfile)<<endl;
  }
    
  BuildRNAModel BRmodel_one(cbs);
  double moved_distance_max=BRmodel_one.build_RNA_model(PRmotif_one.motifs,readrestr_one);
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
  double A_soft=100.0;
  readcoeff_one.reset_paircoeff_soft(A_soft,0.0);
  vector<string> Vpaircoeff_strong=readcoeff_one.get_pair_coeff();
  A_soft=50.0;
  readcoeff_one.reset_paircoeff_soft(A_soft,0.0);
  vector<string> Vpaircoeff_mediun=readcoeff_one.get_pair_coeff();
  A_soft=10.0;
  double sigma_offset=1.5;
  readcoeff_one.reset_paircoeff_soft(A_soft,sigma_offset);
  readcoeff_one.get_interaction_style();
  PrepareLAMMPS prelmp_one(cbs);
  double Temp=300.0;
  prelmp_one.set_temperature(Temp,Temp);

  double moved_distance_max_limited=moved_distance_max;
  if(Nstep<0) {
    Nstep=1000000+round(10*moved_distance_max)*100000;
    if(Nstep>60000000) Nstep=60000000;
    if(moved_distance_max_limited > 60.0)
      moved_distance_max_limited=60.0;
  }

  string restr_TI_type="harmonic 0.0001 "+to_string(moved_distance_max_limited/4);
  string FLEXIBLE_LOOP[]={"Loop","Loop_kiss_loop","Unpaired_chain","Long_loop"};
  string HELIX[]={"Helix","Tertiary_helix","THelix"};
  //double bcut0=3*cbs->init->assembly_delta_distance;
  //if(moved_distance_max_limited/2 < bcut0)
    //bcut0=moved_distance_max_limited/2;
  double bcut0=moved_distance_max_limited/4;
  string restr_bond0_label="DisRes0";
  string restr_bond0_coeff="bond_coeff "+restr_bond0_label+" harmonic 0.01 "+to_string(bcut0/4);
  readcoeff_one.add_coeff(restr_bond0_coeff);
  string restr_bond_label="DisRes";
  string restr_bond_coeff="bond_coeff "+restr_bond_label+" harmonic 0.001 "+to_string(moved_distance_max_limited/4);
  readcoeff_one.add_coeff(restr_bond_coeff);
  prelmp_one.set_coeff(readcoeff_one);
  int flag_loop=0;
  int flag_motif=0;
  int nmodel_total=0;
  int nbond_restr=0;
  vector<vector<Atom> > VVatoms_cg;
  vector<int> Vreserved_atomids;
  vector<string> Vreserved_atomnames;
  vector<string> VTHelix_resnames;
  for(size_t ii=0;ii<BRmodel_one.models.size();ii++) {
    vector<Atom> Vatoms_cg, Vatoms_cg_modified;
    vector<Bond> Vbonds;
    vector<Angle> Vangles;
    vector<Dihedral> Vdihedrals, Vimpropers;
    int nmotif=1;
    size_t index_offset=0;
    vector<string> Vresnames_user_template;
    for(size_t jj=0;jj<BRmodel_one.models[ii].size();jj++) {
      RNAMotif &mtf_one=BRmodel_one.models[ii][jj];
      if(mtf_one.type!="user_template") continue;
      Molecule mol_one=mtf_one;
      vector<string>& Vresnames_one=mol_one.Vresnames;
      Vresnames_user_template.insert(Vresnames_user_template.end(),Vresnames_one.begin(),Vresnames_one.end());
      RNACoarseGrained rnacg_one(cbs);
      if(cbs->init->RNA_version=="latest")
        rnacg_one.convert_RNA_to_CG(mol_one,index_offset);
      else if(cbs->init->RNA_version=="v2")
        rnacg_one.convert_RNA_to_CG_v2(mol_one,index_offset);
      index_offset+=rnacg_one.number_of_CG_atoms();
      vector<Atom> Vatoms_list=rnacg_one.rna_cg.get_atom_list();
      if(ii==0) {
        for(size_t kk=0;kk<Vatoms_list.size();kk++) {
          string atomname_one=Vatoms_list[kk].get_atom_name();
          Vreserved_atomnames.push_back(atomname_one);
          Vreserved_atomids.push_back(Vatoms_list[kk].id);
        }
      }
      Vatoms_cg.insert(Vatoms_cg.end(),Vatoms_list.begin(),Vatoms_list.end());

      for(size_t kk=0;kk<Vatoms_list.size();kk++) {
        Atom atom_one=Vatoms_list[kk];
        atom_one.chainid=to_string(nmotif);
        Vatoms_cg_modified.push_back(atom_one);
      }
      nmotif++;
      flag_motif=1;
    }

    for(size_t jj=0;jj<BRmodel_one.models[ii].size();jj++) {
      RNAMotif &mtf_one=BRmodel_one.models[ii][jj];
      if(mtf_one.type=="user_template" || count(begin(HELIX),end(HELIX),mtf_one.type)) continue;
      Molecule mol_one=mtf_one;
      vector<string> Vresnames_j=mol_one.Vresnames;
      RNACoarseGrained rnacg_one(cbs);
      if(cbs->init->RNA_version=="latest")
        rnacg_one.convert_RNA_to_CG(mol_one,index_offset);
      else if(cbs->init->RNA_version=="v2")
        rnacg_one.convert_RNA_to_CG_v2(mol_one,index_offset);
      index_offset+=rnacg_one.number_of_CG_atoms();
      vector<Atom> Vatoms_list=rnacg_one.rna_cg.get_atom_list();
      if(ii==0) {
        for(size_t kk=0;kk<Vatoms_list.size();kk++) {
          string atomname_one=Vatoms_list[kk].get_atom_name();
          if(count(Vreserved_atomnames.begin(),Vreserved_atomnames.end(),atomname_one)==0) {
            Vreserved_atomnames.push_back(atomname_one);
            Vreserved_atomids.push_back(Vatoms_list[kk].id);
          }
        }
      }
      Vatoms_cg.insert(Vatoms_cg.end(),Vatoms_list.begin(),Vatoms_list.end());

      int relation=Measure::vector_relation(Vresnames_user_template,Vresnames_j);

      if(count(begin(FLEXIBLE_LOOP),end(FLEXIBLE_LOOP),mtf_one.type) || 
         (Vresnames_user_template.size()>0 && (relation==Measure::contained || relation==Measure::equal))) {
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

    for(size_t jj=0;jj<BRmodel_one.models[ii].size();jj++) {
      RNAMotif &mtf_one=BRmodel_one.models[ii][jj];
      if(mtf_one.type=="user_template" || count(begin(HELIX),end(HELIX),mtf_one.type)==0) continue;
      Molecule mol_one=mtf_one;
      RNACoarseGrained rnacg_one(cbs);
      if(cbs->init->RNA_version=="latest")
        rnacg_one.convert_RNA_to_CG(mol_one,index_offset);
      else if(cbs->init->RNA_version=="v2")
        rnacg_one.convert_RNA_to_CG_v2(mol_one,index_offset);
      index_offset+=rnacg_one.number_of_CG_atoms();
      vector<Atom> Vatoms_list=rnacg_one.rna_cg.get_atom_list();
      if(ii==0 && mtf_one.type!="THelix") {
        for(size_t kk=0;kk<Vatoms_list.size();kk++) {
          string atomname_one=Vatoms_list[kk].get_atom_name();
          if(count(Vreserved_atomnames.begin(),Vreserved_atomnames.end(),atomname_one)==0) {
            Vreserved_atomnames.push_back(atomname_one);
            Vreserved_atomids.push_back(Vatoms_list[kk].id);
          }
        }
      }
      Vatoms_cg.insert(Vatoms_cg.end(),Vatoms_list.begin(),Vatoms_list.end());

      if(ii==0 && mtf_one.type=="THelix") {
        if(mtf_one.chains[0].chainid==mtf_one.chains[1].chainid) {
          vector<string> Vresnames_one=rnacg_one.rna_cg.Vresnames;
          VTHelix_resnames.insert(VTHelix_resnames.end(),Vresnames_one.begin(),Vresnames_one.end());
        }
      }

      for(size_t kk=0;kk<Vatoms_list.size();kk++) {
        Atom atom_one=Vatoms_list[kk];
        atom_one.chainid=to_string(nmotif);
        Vatoms_cg_modified.push_back(atom_one);
      }
      nmotif++;
      flag_motif=1;
    }

    VVatoms_cg.push_back(Vatoms_cg);
    nmodel_total++;

    vector<pair<int,string> > Vatomid_name;
    map<int,string> Matom_resname;
    vector<string> Vresname_bond0;
    map<int,vector<double> > Matomid_xyz;
    for(size_t jj=0;jj<Vatoms_cg.size();jj++) {
      int atomid_one=Vatoms_cg[jj].id;
      string atomname_full=Vatoms_cg[jj].get_atom_name();
      pair<int,string> p_one=make_pair(atomid_one,atomname_full);
      Vatomid_name.push_back(p_one);
      string resname_one=Vatoms_cg[jj].get_residue_name();
      Matom_resname.emplace(atomid_one,resname_one);
      vector<double> xyz=Vatoms_cg[jj].get_atom_xyz();
      Matomid_xyz.emplace(atomid_one,xyz);
    }
    for(size_t jj=0;jj<Vatomid_name.size()-1;jj++) {
      int atomid_i=Vatomid_name[jj].first;
      string atomname_i=Vatomid_name[jj].second;
      string resname_i=Matom_resname[atomid_i];
      for(size_t kk=jj+1;kk<Vatomid_name.size();kk++) {
        if(atomname_i==Vatomid_name[kk].second) {
          Bond bond_one;
          bond_one.atom_i=atomid_i;
          bond_one.atom_j=Vatomid_name[kk].first;
          if(count(VTHelix_resnames.begin(),VTHelix_resnames.end(),resname_i))
            bond_one.type=restr_TI_type;
          else {
            //bond_one.type=restr_type;
            vector<double> xyz_i=Matomid_xyz[bond_one.atom_i];
            vector<double> xyz_j=Matomid_xyz[bond_one.atom_j];
            double dis_one=Measure::distance(xyz_i,xyz_j);
            if(dis_one <= bcut0) {
              bond_one.type=restr_bond0_label;
              if(count(Vresname_bond0.begin(),Vresname_bond0.end(),resname_i)==0)
                Vresname_bond0.push_back(resname_i);
            }
            else {
              bond_one.type=restr_bond_label;
            }
          }
          Vbonds.push_back(bond_one);
          if(ii==0) nbond_restr++;
        }
      }
    }

    for(size_t kk=0;kk<Vbonds.size();kk++) {
      if(Vbonds[kk].type!=restr_bond_label) continue;
      if(Matom_resname.count(Vbonds[kk].atom_i)==0) continue;
      string resname_i=Matom_resname[Vbonds[kk].atom_i];
      if(count(Vresname_bond0.begin(),Vresname_bond0.end(),resname_i)) {
        Vbonds[kk].type=restr_bond0_label;
      }
    }

    if(Vbonds.size()>0) {
      Bond bond_one;
      bond_one.atom_i=-1;
      bond_one.atom_j=-1;
      bond_one.type=restr_bond0_label;
      Vbonds.push_back(bond_one);
      bond_one.type=restr_bond_label;
      Vbonds.push_back(bond_one);
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

  //for(size_t kk=0;kk<Vreserved_atomids.size();kk++) {
    //cout<<Vreserved_atomids[kk]<<" "<<Vreserved_atomnames[kk]<<endl;
  //}
  
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
  string KBOND0[]={"0.01","0.1","1.0","5.0","10.0"};
  double DR0[]={bcut0/4,bcut0/8,bcut0/16,0.0,0.0};
  string KBOND[]={"0.001","0.01","0.1","1.0","10.0"};
  double B0[]={moved_distance_max_limited/4,moved_distance_max_limited/8,moved_distance_max_limited/16,0.0,0.0};
  string KBOND_TI[]={"0.0001","0.001","0.01","0.1","1.0"};
  double B0_TI[]={moved_distance_max_limited/4,moved_distance_max_limited/8,moved_distance_max_limited/16,moved_distance_max_limited/32,0.0};
  string TIMESTEP[]={"0.2","0.5","0.75","1.0","1.0"};
  string restr_TI_bond_label="";
  if(nbond_restr>0 && prelmp_one.Mrestraint_type.count(restr_TI_type)) {
    restr_TI_bond_label=prelmp_one.Mrestraint_type[restr_TI_type];
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
      if(restr_bond0_label.length()>0) {
        if(prelmp_one.Vbondstyle.size()==1) {
          sline="bond_coeff      "+restr_bond0_label+" "+KBOND0[kk]+" "+to_string(DR0[kk]);
        }
        else {
          sline="bond_coeff      "+restr_bond0_label+" harmonic "+KBOND0[kk]+" "+to_string(DR0[kk]);
        }
        Vrun.push_back(sline);
      }
      if(restr_bond_label.length()>0) {
        if(prelmp_one.Vbondstyle.size()==1) {
          sline="bond_coeff      "+restr_bond_label+" "+KBOND[kk]+" "+to_string(B0[kk]);
        }
        else {
          sline="bond_coeff      "+restr_bond_label+" harmonic "+KBOND[kk]+" "+to_string(B0[kk]);
        }
        Vrun.push_back(sline);
      }
      if(restr_TI_bond_label.length()>0) {
        if(prelmp_one.Vbondstyle.size()==1) {
          sline="bond_coeff      "+restr_TI_bond_label+" "+KBOND_TI[kk]+" "+to_string(B0_TI[kk]);
        }
        else {
          sline="bond_coeff      "+restr_TI_bond_label+" harmonic "+KBOND_TI[kk]+" "+to_string(B0_TI[kk]);
        }
        Vrun.push_back(sline);
      }
      
      if(kk==n_splited-2) {
        sline="";
        Vrun.push_back(sline);
        for(size_t nn=0;nn<Vpaircoeff_mediun.size();nn++) {
          Vrun.push_back(Vpaircoeff_mediun[nn]);
        }
      }
      else if(kk==n_splited-1) {
        sline="";
        Vrun.push_back(sline);
        for(size_t nn=0;nn<Vpaircoeff_strong.size();nn++) {
          Vrun.push_back(Vpaircoeff_strong[nn]);
        }
      }

      sline="";
      Vrun.push_back(sline);
      sline="timestep        "+TIMESTEP[kk];
      Vrun.push_back(sline);
      sline="run             "+to_string(Nstep_splited[kk]);
      Vrun.push_back(sline);
    }
  }
  else {
    sline="run             0";
    Vrun.push_back(sline);
  }

  vector<string> Vem_skip;  //skip energy minimization
  sline="";
  Vem_skip.push_back(sline);

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
    prelmp_one.Vem=Vem_skip;
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
    prelmp_one.Vem=Vem_skip;
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
  vector<string> Vresnames=rna_one.Vresnames;
  for(size_t ii=0;ii<VVatoms_cg.size();ii++) {
    vector<Atom> Vreserved_atoms;
    for(size_t jj=0;jj<VVatoms_cg[ii].size();jj++) {
      int id_one=VVatoms_cg[ii][jj].id;
      if(count(Vreserved_atomids.begin(),Vreserved_atomids.end(),id_one))
        Vreserved_atoms.push_back(VVatoms_cg[ii][jj]);
    }
    vector<Atom> Vatoms_sorted;
    for(size_t jj=0;jj<Vresnames.size();jj++) {
      for(size_t kk=0;kk<Vreserved_atoms.size();kk++) {
        Atom atom_one=Vreserved_atoms[kk];
        atom_one.resname=atom_one.resname.substr(1);
        if(atom_one.get_residue_name()==Vresnames[jj]) {
          Vatoms_sorted.push_back(Vreserved_atoms[kk]);
        }
      }
    }
    for(size_t kk=0;kk<Vatoms_sorted.size();kk++) {
      Vatoms_sorted[kk].id=static_cast<int>(kk+1);
    }
    VVatoms_cg_sorted.push_back(Vatoms_sorted);
  }

  if((int)VVatoms_cg_sorted.size()> cbs->init->Num_model_max) {
    double rmsd_cutoff=0.1*nres;
    double Pclu=0.5;
    double Pclu_min=0.2;
    Cluster clu_one(cbs);
    clu_one.input_data(VVatoms_cg_sorted);
    clu_one.calculate_RMSD_matrix();
    clu_one.get_clusters(rmsd_cutoff,Pclu,Pclu_min);
    
    (*cbs->logfile)<<endl;
    (*cbs->logfile)<<clu_one;
    (*cbs->logfile)<<endl;

    vector<vector<Atom> > VVtmp=VVatoms_cg_sorted;
    VVatoms_cg_sorted.clear();
    int nclu=static_cast<int>(clu_one.clusters.size());
    if(nclu > cbs->init->Num_model_max)
      nclu=cbs->init->Num_model_max;
    for(int jj=0;jj<nclu;jj++) {
      size_t index_clu=clu_one.clusters[jj][0];
      VVatoms_cg_sorted.push_back(VVtmp[index_clu]);
    }
  }

  for(size_t ii=0;ii<VVatoms_cg_sorted.size();ii++) {
    double xcm=0.0, ycm=0.0, zcm=0.0;
    size_t natoms=VVatoms_cg_sorted[ii].size();
    for(size_t jj=0;jj<natoms;jj++) {
      xcm+=VVatoms_cg_sorted[ii][jj].x;
      ycm+=VVatoms_cg_sorted[ii][jj].y;
      zcm+=VVatoms_cg_sorted[ii][jj].z;
    }
    if(natoms>0) {
      xcm/=natoms;
      ycm/=natoms;
      zcm/=natoms;
    }
    for(size_t jj=0;jj<natoms;jj++)
      VVatoms_cg_sorted[ii][jj].move_atom(-xcm,-ycm,-zcm);
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

  time (&tend);

  double dtime=difftime(tend,tstart);
  string stime=cbs->calculate_runtime(dtime);
  cbs->error->log(stime);
  stime="";
  cbs->error->log(stime);

  delete cbs;
  
  return 0;
}