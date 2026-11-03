#include "basic.h"
#include "error.h"
#include "readfile.h"
#include "initialization.h"
#include "superpose.h"
#include "cluster.h"
#include "sort_vector.h"
#include "Update_RNA_template.h"

using namespace CBS_NS;

const int RESID_MAX = 9999;
const size_t LENGTH_MAX = 10000;
const size_t N_CLUSTER_MAX = 30000;

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

  if(Mcmd.count("-inDIR")==0) {
    cout<<"Wrong commands: Please give the directory for pdb/cif files (-inDIR)!"<<endl;
    exit(0);
  }
  string inDIR=Mcmd["-inDIR"][0];
  Mcmd.erase("-inDIR");
  if(inDIR.back()!='/')
    inDIR+="/";
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
  int check_template=0;
  if(Mcmd.count("-check")) {
    check_template=1;
    Mcmd.erase("-check");
  }

  cbs= new Cobmms(flog);
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

  time_t now=time(0);
  tm *ltm=localtime(&now);
  int year=1900+ltm->tm_year;
  int month=1+ltm->tm_mon;
  int day=ltm->tm_mday;
  _version_ = "Version: "+to_string(year)+"-"+to_string(month)+"-"+to_string(day);

  vector<string> Tfile=get_dir_files(inDIR);
  if(Tfile.size()==0) return 0;
  
  string tempDIR=cbs->init->RNA_TEMPLATE_DIR;

  if(!is_directory_exist(tempDIR)) {
    string cmd="mkdir "+tempDIR;
    cbs->execute_system_command(cmd);
  }

  string pdbDIR=tempDIR+"pdb";
  if(!is_directory_exist(pdbDIR)) {
    string cmd="mkdir "+pdbDIR;
    cbs->execute_system_command(cmd);
  }
  
  if(check_template)
    CBS_NS::check_template_list(tempDIR);
  
  PDB_list.clear();
  string fpdblist=tempDIR+"pdb_list.dat";
  if(is_file_exist(fpdblist)) {
    vector<string> Vlines=read_file(fpdblist);
    for(size_t ii=0;ii<Vlines.size();ii++) {
      vector<string> stoken=split_string(Vlines[ii]);
      PDB_list.push_back(stoken[0]);
    }
  }

  vector<string> Vfile;
  int n_skip=0;
  for(size_t ii=0;ii<Tfile.size();ii++) {
    vector<string> stoken=split_string(Tfile[ii]);
    if(stoken[0].find(".pdb")==string::npos && stoken[0].find(".cif")==string::npos) {
      string str="Skip file "+stoken[0];
      cbs->error->log(str);
      cout<<str<<endl;
      continue;
    }
    string PDBid=stoken[0].substr(0,stoken[0].size()-4);
    for(size_t ii=0;ii<PDBid.size();ii++)
      PDBid[ii]=tolower(PDBid[ii]);
    if(count(PDB_list.begin(),PDB_list.end(),PDBid)) {
      n_skip++;
      continue;
    }
    Vfile.push_back(stoken[0]);
  }
  if(n_skip>0) {
    string str="Skip "+to_string(n_skip)+" existed PDBs!";
    cbs->error->log(str);
    cout<<str<<endl;
  }

  string str="Step 1: Extracting RNA templates from "+to_string(Vfile.size())+" PDBs.";
  cbs->error->log(str);
  cout<<str<<endl;

  int ncount=0;
  int nthread=cbs->init->nthread;
  size_t total_length=0;
  vector<string> pdbs;
  vector<string> seqs;
  vector<vector<string> > Vresnames;
  vector<RNAMolecule> rnas;
  vector<RNABasePair> basepairs;
  vector<ExtractRNATemplate> temps;
  sort(Vfile.begin(),Vfile.end());
  for(int ii=(int)Vfile.size()-1;ii>=0;ii--) {
    cbs->error->log(Vfile[ii]);
    cout<<Vfile[ii]<<endl;
    vector<string> stoken=split_string(Vfile[ii]);
    string fstr=inDIR+stoken[0];
    
    ReadPDB readpdb_one(cbs,fstr);
    string PDBid=readpdb_one.PDBid;
    PDB_list.push_back(PDBid);
   
    if(readpdb_one.RNA_atoms.size()==0 || readpdb_one.RNA_atoms[0].size()==0) {
      cbs->error->log("**********************************************");
      cbs->error->log(PDBid);
      cbs->error->log("Skipped!");
      continue;
    }

    RNAMolecule rna(cbs);
    rna.set_RNA_molecule(readpdb_one.RNA_atoms[0]);
    if(rna.number_of_residues()<5) {
      cbs->error->log("**********************************************");
      cbs->error->log(PDBid);
      cbs->error->log("Skipped!");
      continue;
    }

    if(rna.chains.size()>cbs->init->CHAIN_ID_LIST.size()) {
      RNAMolecule rna_tmp=rna;
      rna.clear();
      for(size_t ii=0;ii<rna_tmp.chains.size();ii++) {
        if(rna_tmp.chains[ii].residues.size()<=3) continue;
        rna.insert_chain(rna_tmp.chains[ii]);
      }
    }

    vector<string> Vresnames_old=rna.Vresnames;
    rna.standardize_molecule();

    RNABasePair basepair_one(cbs);
    string secstr_one=rna.extract_secondary_structure();
    vector<vector<string> > Vres=rna.get_residue_name();
    basepair_one.parse_RNA_basepairs(Vres,secstr_one);
    //cout<<basepair_one;
    if(basepair_one.size()<2) {
      cbs->error->log("**********************************************");
      cbs->error->log(PDBid);
      cbs->error->log("Skipped!");
      continue;
    }

    //for(size_t kk=0;kk<basepair_one.Vresnames.size();kk++)
      //cout<<kk<<" :" <<basepair_one.Vresnames[kk]<<endl;
    
    rna.simplify_sidechain();
    
    pdbs.push_back(PDBid);
    seqs.push_back(rna.get_sequences());
    Vresnames.push_back(Vresnames_old);
    rnas.push_back(rna);
    basepairs.push_back(basepair_one);
    ExtractRNATemplate ERtemp_one(cbs);
    temps.push_back(ERtemp_one);
    ncount++;
    total_length+=rna.number_of_residues();
    if(ncount%nthread==0 || total_length > LENGTH_MAX) {
      vector<std::thread> threads;
      for(int ni=0;ni<ncount;ni++)
        threads.emplace_back(std::thread(CBS_NS::extract_template,std::ref(rnas[ni]),std::ref(basepairs[ni]),std::ref(temps),ni));
      for(int ni=0;ni<ncount;ni++)
        threads[ni].join();
      
      for(int kk=0;kk<ncount;kk++) {
        cbs->error->log("**********************************************");
        cbs->error->log(pdbs[kk]);
        cbs->error->log(seqs[kk]);
        cbs->error->log(basepairs[kk].secstr_out);
        if(temps[kk].Vinfo.size()!=0) {
          for(auto sline : temps[kk].Vinfo)
            (*cbs->logfile)<<sline<<endl;
        }
        else
          (*cbs->logfile)<<temps[kk];
      }
      
      Vpdbs.insert(Vpdbs.end(),pdbs.begin(),pdbs.end());
      VVresnames_old.insert(VVresnames_old.end(),Vresnames.begin(),Vresnames.end());
      Vrnas.insert(Vrnas.end(),rnas.begin(),rnas.end());
      Vtemplates.insert(Vtemplates.end(),temps.begin(),temps.end());

      ncount=0;
      total_length=0;
      pdbs.clear();
      seqs.clear();
      Vresnames.clear();
      rnas.clear();
      basepairs.clear();
      temps.clear();

      CBS_NS::write_info(tempDIR);
    }
  }
  if(ncount>0) {
    vector<std::thread> threads;
    for(int ni=0;ni<ncount;ni++)
      threads.emplace_back(std::thread(CBS_NS::extract_template,std::ref(rnas[ni]),std::ref(basepairs[ni]),std::ref(temps),ni));
    for(int ni=0;ni<ncount;ni++)
      threads[ni].join();
    
    for(int kk=0;kk<ncount;kk++) {
      cbs->error->log("**********************************************");
      cbs->error->log(pdbs[kk]);
      cbs->error->log(seqs[kk]);
      cbs->error->log(basepairs[kk].secstr_out);
      if(temps[kk].Vinfo.size()!=0) {
        for(auto sline : temps[kk].Vinfo)
          (*cbs->logfile)<<sline<<endl;
      }
      else
        (*cbs->logfile)<<temps[kk];
    }
    
    Vpdbs.insert(Vpdbs.end(),pdbs.begin(),pdbs.end());
    VVresnames_old.insert(VVresnames_old.end(),Vresnames.begin(),Vresnames.end());
    Vrnas.insert(Vrnas.end(),rnas.begin(),rnas.end());
    Vtemplates.insert(Vtemplates.end(),temps.begin(),temps.end());

    ncount=0;
    total_length=0;
    pdbs.clear();
    seqs.clear();
    Vresnames.clear();
    rnas.clear();
    basepairs.clear();
    temps.clear();

    CBS_NS::write_info(tempDIR);
  }

  ostringstream oss;
  oss<<std::fixed<<std::setprecision(2)<<cbs->init->RNA_template_RMSD_cutoff;
  str="Step 2: Clustering templates through RMSD cutoff "+oss.str()+" Angstrom.";
  cbs->error->log(str);
  cout<<str<<endl;
  map<string,int> template_count;
  string fversion=tempDIR+"version.txt";
  if(is_file_exist(fversion)) {
    vector<string> Vlines=read_file(fversion);
    int check=0;
    for(size_t ii=0;ii<Vlines.size();ii++) {
      if(Vlines[ii].find("#List of templates")!=string::npos)
        {check=1; continue;}
      if(check==0) continue;
      vector<string> stoken=split_string(Vlines[ii],":");
      template_count.insert(pair<string,int>(stoken[0],stoi(stoken[1])));
    }
  }
  
  map<string,vector<string> > Mtemplates;
  string flist=tempDIR+"template_list.dat";
  if(!is_file_exist(flist)) return 0;
  vector<string> Vlines=read_file(flist);
  for(auto sline : Vlines) {
    vector<string> stoken=split_string(sline);
    string type=stoken[1];
    if(type=="Helix") continue;
    string name=type+"#"+stoken[2];
    if(Mtemplates.count(name)==0) {
      vector<string> Vtmp;
      Vtmp.push_back(sline);
      Mtemplates.emplace(name,Vtmp);
    }
    else {
      Mtemplates[name].push_back(sline);
    }
  }

  for(auto& it : Mtemplates) {
    string name=it.first;
    vector<string> Vlist=it.second;
    ReadRNATemplate readRT_one(cbs,name);

    if(template_count.count(name)==0) {
      readRT_one.read_templates();
      int ncount=static_cast<int>(readRT_one.template_list.size());
      template_count.emplace(name,ncount);
    }

    readRT_one.read_full_template_list();
    size_t count_old=readRT_one.full_template_list.size();
    for(auto sline : Vlist) {
      if(count(readRT_one.full_template_list.begin(),readRT_one.full_template_list.end(),sline)==0)
        readRT_one.full_template_list.push_back(sline);
    }
    size_t count_new=readRT_one.full_template_list.size();
    if(count_new == count_old) continue;
    readRT_one.write_full_template_list();

    string sout="  "+name+" : "+to_string(count_old)+" --> "+to_string(count_new);
    cbs->error->log(sout);
    cout<<sout<<endl;

    vector<string> Vcluster_list,Vappend_list;
    if(count_new<=N_CLUSTER_MAX) {
      Vcluster_list=readRT_one.full_template_list;
    }
    else if(count_new>N_CLUSTER_MAX && count_old<N_CLUSTER_MAX) {
      for(size_t kk=0;kk<N_CLUSTER_MAX;kk++)
        Vcluster_list.push_back(readRT_one.full_template_list[kk]);
      for(size_t kk=N_CLUSTER_MAX;kk<count_new;kk++)
        Vappend_list.push_back(readRT_one.full_template_list[kk]);
    }
    else if(count_old>=N_CLUSTER_MAX) {
      for(size_t kk=count_old;kk<count_new;kk++)
        Vappend_list.push_back(readRT_one.full_template_list[kk]);
    }

    readRT_one.template_list.clear();
    vector<vector<double> > VVxyz;

    int ntemplates=(int)Vcluster_list.size();  //Clustering templates by RMSD cutoff
    string rmsdDIR=tempDIR+"rmsd_matrix/";
    if(!is_directory_exist(rmsdDIR)) {
      string cmd="mkdir "+rmsdDIR;
      cbs->execute_system_command(cmd);
    }
    if(ntemplates>0) {
      //if(ntemplates > 2000) {
        //sout="    Reading template xyz...";
        //cbs->error->log(sout);
        //cout<<sout<<endl;
      //}
      Vpdbs.clear();
      ncount=0;
      vector<string> Vtemp_info;
      map<string,int> Mpdbs;
      size_t Nnts=0;
      for(int jj=0;jj<ntemplates;jj++) {
        vector<string> stoken=split_string(Vcluster_list[jj]);
        string PDBid=stoken[0];
        if(jj==0) {
          string str_seq=stoken[4];
          vector<string> stoken_s=split_string(str_seq,"-&");
          for(size_t kk=0;kk<stoken_s.size();kk++)
            Nnts+=stoken_s[kk].size();
        }
        if(count(Vpdbs.begin(),Vpdbs.end(),PDBid)==0) {
          Vpdbs.push_back(PDBid);
          Mpdbs.emplace(PDBid,ncount);
          ncount++;
        }
        Vtemp_info.push_back(Vcluster_list[jj]);
        if(ncount%nthread==0 || jj==ntemplates-1) {
          vector<string> Vfstr;
          vector<ReadPDB> Vreadpdbs;
          vector<Molecule> Vmoles;
          for(int kk=0;kk<ncount;kk++) {
            string fstr=tempDIR+"pdb/"+Vpdbs[kk]+".pdb";
            Vfstr.push_back(fstr);
            ReadPDB readpdb_one(cbs);
            Vreadpdbs.push_back(readpdb_one);
            Molecule mole_one(cbs);
            Vmoles.push_back(mole_one);
          }
          vector<std::thread> threads;
          for(int ni=0;ni<ncount;ni++)
            threads.emplace_back(std::thread(CBS_NS::read_pdb_xyz,Vfstr[ni],std::ref(Vreadpdbs[ni]),std::ref(Vmoles),ni));
          for(int ni=0;ni<ncount;ni++)
            threads[ni].join();
          Vfstr.clear();
          Vreadpdbs.clear();
        
          int nloop=(int)Vtemp_info.size()/nthread;
          int nlast=(int)Vtemp_info.size()%nthread;
          for(int nl=0;nl<nloop;nl++) {
            int nl0=nl*nthread;
            vector<string> templates_info;
            templates_info.resize(nthread);
            pdbs.clear();
            for(int kk=0;kk<nthread;kk++) {
              templates_info[kk]=Vtemp_info[nl0+kk];
              stoken=split_string(templates_info[kk]);
              pdbs.push_back(stoken[0]);
            }
            vector<vector<double> > Vxyz;
            Vxyz.resize(nthread);
            vector<std::thread> threads;
            for(int ni=0;ni<nthread;ni++) {
              Molecule &mole_ni=Vmoles[Mpdbs[pdbs[ni]]];
              threads.emplace_back(std::thread(CBS_NS::read_template_xyz,std::ref(mole_ni),templates_info[ni],std::ref(Vxyz),ni));
            }
            for(int ni=0;ni<nthread;ni++)
              threads[ni].join();
            for(int kk=0;kk<nthread;kk++) {
              VVxyz.push_back(Vxyz[kk]);
            }  
          }
          if(nlast!=0) {
            int nl0=nloop*nthread;
            vector<string> templates_info;
            templates_info.resize(nlast);
            pdbs.clear();
            for(int kk=0;kk<nlast;kk++) {
              templates_info[kk]=Vtemp_info[nl0+kk];
              stoken=split_string(templates_info[kk]);
              pdbs.push_back(stoken[0]);
            }
            vector<vector<double> > Vxyz_last;
            Vxyz_last.resize(nlast);
            vector<std::thread> threads;
            for(int ni=0;ni<nlast;ni++) {
              Molecule &mole_ni=Vmoles[Mpdbs[pdbs[ni]]];
              threads.emplace_back(std::thread(CBS_NS::read_template_xyz,std::ref(mole_ni),templates_info[ni],std::ref(Vxyz_last),ni));
            }
            for(int ni=0;ni<nlast;ni++)
              threads[ni].join();
            for(int kk=0;kk<nlast;kk++) {
              VVxyz.push_back(Vxyz_last[kk]);
            }
          }

          ncount=0;
          Vpdbs.clear();
          Vtemp_info.clear();
          Mpdbs.clear();
        }
      }  
      /*
      for(auto sline : readRT_one.full_template_list) {
        vector<double> Vxyz=readRT_one.read_template_xyz(sline).get_molecule_xyz();
        VVxyz.push_back(Vxyz);
      }
      */

      if((int)VVxyz.size()!=ntemplates) {
        string str="Update_RNA_template.cpp: VVxyz.size()!=ntemplates for clustering templates --> "+to_string(VVxyz.size())+" : "+to_string(ntemplates);
        cbs->error->one(str);
      }

      if(VVxyz[0].size()!=45*Nnts) {
        string str="Update_RNA_template.cpp: VVxyz[0].size()!=15*Nnts --> Nnts = "+to_string(Nnts)+" : Ndata = "+to_string(VVxyz[0].size());
        cbs->error->one(str);
      }

      if(ntemplates >= 2000) {
        sout="    Clustering templates...";
        cbs->error->log(sout);
        cout<<sout<<endl;
      }
      Cluster cluster_one(cbs,VVxyz);
      string frmsd=rmsdDIR+readRT_one.template_name+"_rmsd.dat";
      cluster_one.calculate_RMSD_matrix(frmsd);
      cluster_one.get_clusters(cbs->init->RNA_template_RMSD_cutoff);
      for(size_t ii=0;ii<cluster_one.clusters.size();ii++) {
        vector<string> Vtemps;
        for(size_t jj=0;jj<cluster_one.clusters[ii].size();jj++) {
          size_t index=cluster_one.clusters[ii][jj];
          Vtemps.push_back(readRT_one.full_template_list[index]);
        }
        readRT_one.template_list.push_back(Vtemps);
      }
    
      readRT_one.write_templates();

      int ncount=static_cast<int>(readRT_one.template_list.size());
      if(template_count.count(name)) {
        template_count[name]=ncount;
      }
      else {
        template_count.emplace(name,ncount);
      }
    }

    ntemplates=(int)Vappend_list.size();  //Appending new templates to existed templates
    VVxyz.clear();
    if(ntemplates>0) {
      //if(ntemplates > 2000) {
        //sout="    Reading template xyz...";
        //cbs->error->log(sout);
        //cout<<sout<<endl;
      //}
      Vpdbs.clear();
      ncount=0;
      vector<string> Vtemp_info;
      map<string,int> Mpdbs;
      size_t Nnts=0;
      for(int jj=0;jj<ntemplates;jj++) {
        vector<string> stoken=split_string(Vappend_list[jj]);
        string PDBid=stoken[0];
        if(jj==0) {
          string str_seq=stoken[4];
          vector<string> stoken_s=split_string(str_seq,"-&");
          for(size_t kk=0;kk<stoken_s.size();kk++)
            Nnts+=stoken_s[kk].size();
        }
        if(count(Vpdbs.begin(),Vpdbs.end(),PDBid)==0) {
          Vpdbs.push_back(PDBid);
          Mpdbs.emplace(PDBid,ncount);
          ncount++;
        }
        Vtemp_info.push_back(Vappend_list[jj]);
        if(ncount%nthread==0 || jj==ntemplates-1) {
          vector<string> Vfstr;
          vector<ReadPDB> Vreadpdbs;
          vector<Molecule> Vmoles;
          for(int kk=0;kk<ncount;kk++) {
            string fstr=tempDIR+"pdb/"+Vpdbs[kk]+".pdb";
            Vfstr.push_back(fstr);
            ReadPDB readpdb_one(cbs);
            Vreadpdbs.push_back(readpdb_one);
            Molecule mole_one(cbs);
            Vmoles.push_back(mole_one);
          }
          vector<std::thread> threads;
          for(int ni=0;ni<ncount;ni++)
            threads.emplace_back(std::thread(CBS_NS::read_pdb_xyz,Vfstr[ni],std::ref(Vreadpdbs[ni]),std::ref(Vmoles),ni));
          for(int ni=0;ni<ncount;ni++)
            threads[ni].join();
          Vfstr.clear();
          Vreadpdbs.clear();
        
          int nloop=(int)Vtemp_info.size()/nthread;
          int nlast=(int)Vtemp_info.size()%nthread;
          for(int nl=0;nl<nloop;nl++) {
            int nl0=nl*nthread;
            vector<string> templates_info;
            templates_info.resize(nthread);
            pdbs.clear();
            for(int kk=0;kk<nthread;kk++) {
              templates_info[kk]=Vtemp_info[nl0+kk];
              stoken=split_string(templates_info[kk]);
              pdbs.push_back(stoken[0]);
            }
            vector<vector<double> > Vxyz;
            Vxyz.resize(nthread);
            vector<std::thread> threads;
            for(int ni=0;ni<nthread;ni++) {
              Molecule &mole_ni=Vmoles[Mpdbs[pdbs[ni]]];
              threads.emplace_back(std::thread(CBS_NS::read_template_xyz,std::ref(mole_ni),templates_info[ni],std::ref(Vxyz),ni));
            }
            for(int ni=0;ni<nthread;ni++)
              threads[ni].join();
            for(int kk=0;kk<nthread;kk++) {
              VVxyz.push_back(Vxyz[kk]);
            }  
          }
          if(nlast!=0) {
            int nl0=nloop*nthread;
            vector<string> templates_info;
            templates_info.resize(nlast);
            pdbs.clear();
            for(int kk=0;kk<nlast;kk++) {
              templates_info[kk]=Vtemp_info[nl0+kk];
              stoken=split_string(templates_info[kk]);
              pdbs.push_back(stoken[0]);
            }
            vector<vector<double> > Vxyz_last;
            Vxyz_last.resize(nlast);
            vector<std::thread> threads;
            for(int ni=0;ni<nlast;ni++) {
              Molecule &mole_ni=Vmoles[Mpdbs[pdbs[ni]]];
              threads.emplace_back(std::thread(CBS_NS::read_template_xyz,std::ref(mole_ni),templates_info[ni],std::ref(Vxyz_last),ni));
            }
            for(int ni=0;ni<nlast;ni++)
              threads[ni].join();
            for(int kk=0;kk<nlast;kk++) {
              VVxyz.push_back(Vxyz_last[kk]);
            }
          }

          ncount=0;
          Vpdbs.clear();
          Vtemp_info.clear();
          Mpdbs.clear();
        }
      }  
      /*
      for(auto sline : readRT_one.full_template_list) {
        vector<double> Vxyz=readRT_one.read_template_xyz(sline).get_molecule_xyz();
        VVxyz.push_back(Vxyz);
      }
      */

      if((int)VVxyz.size()!=ntemplates) {
        string str="Update_RNA_template.cpp: VVxyz.size()!=ntemplates for appending templates --> "+to_string(VVxyz.size())+" : "+to_string(ntemplates);
        cbs->error->one(str);
      }

      if(ntemplates >= 2000) {
        sout="    Appending templates...";
        cbs->error->log(sout);
        cout<<sout<<endl;
      }

      vector<vector<double> > Vxyz_clusters;
      readRT_one.read_templates();
      for(int jj=0;jj<ntemplates;jj++) {
        Superpose sp_one;
        sp_one.set_x1(VVxyz[jj]);
        size_t index_min=-1;
        double rmsd_min=1000.0;
        for(size_t kk=0;kk<readRT_one.template_list.size();kk++) {
          vector<double> xyz_cluster;
          if(kk<Vxyz_clusters.size())
            xyz_cluster=Vxyz_clusters[kk];
          else {
            xyz_cluster=readRT_one.read_template_xyz(readRT_one.template_list[kk][0]);
            Vxyz_clusters.push_back(xyz_cluster);
          }
          sp_one.set_x2(xyz_cluster);
          double rmsd=sp_one.get_RMSD();
          if(rmsd<rmsd_min) {
            rmsd_min=rmsd;
            index_min=kk;
          }
          if(rmsd_min <= cbs->init->RNA_template_RMSD_cutoff/4) {
            break;
          }
        }
       
        if(rmsd_min <= cbs->init->RNA_template_RMSD_cutoff) {
          readRT_one.template_list[index_min].push_back(Vappend_list[jj]);
        }
        else {
          vector<string> Vtmp;
          Vtmp.push_back(Vappend_list[jj]);
          readRT_one.template_list.push_back(Vtmp);
        }
      } 

      vector<vector<string> > Tlist=readRT_one.template_list;
      readRT_one.template_list.clear();
      vector<size_t> Vsize;
      for(size_t jj=0;jj<Tlist.size();jj++) {
        size_t num=Tlist[jj].size();
        Vsize.push_back(num);
      }
      vector<size_t> Vsort=SortVector::sort_vector_descending(Vsize);
      for(size_t jj=0;jj<Vsort.size();jj++) {
        size_t index_i=Vsort[jj];
        readRT_one.template_list.push_back(Tlist[index_i]);
      }
      /*
      vector<string> cluster_size;
      for(size_t jj=0;jj<Tlist.size();jj++) {
        string str=to_string(jj)+" "+to_string(Tlist[jj].size());
        cluster_size.push_back(str);
      }
      SortVector::sort_vector_descending(cluster_size,1);
      for(size_t jj=0;jj<cluster_size.size();jj++) {
        vector<string> stoken=split_string(cluster_size[jj]);
        size_t index=stoi(stoken[0]);
        readRT_one.template_list.push_back(Tlist[index]);
      }
      */
      readRT_one.write_templates();

      template_count[name]=static_cast<int>(readRT_one.template_list.size());
    }  
  }

  map<string,int> type_count;
  for(auto& it : template_count) {
    vector<string> stoken=split_string(it.first,"#");
    string type=stoken[0];
    if(type_count.count(type))
      type_count[type]=type_count[type]+it.second;
    else
      type_count.emplace(type,it.second);
  }

  ofstream outF (fversion.c_str());
  outF<<_version_<<endl;
  for(auto& it : type_count)
    outF<<it.first<<" : "<<it.second<<endl;
  outF<<"#List of templates (name : number of clusters)"<<endl;
  for(auto& it : template_count)
    outF<<it.first<<" : "<<it.second<<endl;
  outF.close(); outF.clear();

  time (&tend);

  double dtime=difftime(tend,tstart);
  string stime=cbs->calculate_runtime(dtime);
  cbs->error->log(stime);

  delete cbs;

  return 0;
}

int CBS_NS::extract_template(const RNAMolecule &rna_one, const RNABasePair &basepair_one, vector<ExtractRNATemplate> &ERtemp, int nth) {
  ERtemp[nth].extract_RNA_templates(rna_one,basepair_one,1);
  
  return 0;
}

void CBS_NS::write_info(string &outDIR) {
  string fstr=outDIR+"template_list.dat";
  ofstream outT(fstr.c_str(),ios::app);
  for(size_t ii=0;ii<VVresnames_old.size();ii++) {
    string PDBid=Vpdbs[ii];
    int cid=0;
    int resid=1;
    vector<string> Vresnames_new;
    vector<string> stoken=split_string(VVresnames_old[ii][0],"/");
    string chainid_p=stoken[0];
    string resid_p=stoken[1];
    string res_one=cbs->init->CHAIN_ID_LIST[cid]+"/"+to_string(resid)+"/"+stoken[2];
    Vresnames_new.push_back(res_one);
    for(size_t jj=1;jj<VVresnames_old[ii].size();jj++) {
      resid++;
      stoken=split_string(VVresnames_old[ii][jj],"/");
      if(stoken[0]!=chainid_p) {
        cid++;
        chainid_p=stoken[0];
        resid=1;
      }
      if(cid >= (int)cbs->init->CHAIN_ID_LIST.size()) {
        cout<<"Warning: Number of chains is "<<cid<<" in "<<PDBid<<endl;
        cid=cid%((int)cbs->init->CHAIN_ID_LIST.size());
      }
      if(resid > RESID_MAX) {
        cout<<"Warning: Number of residues is "<<resid<<" in Chain "<<stoken[0]<<" of "<<PDBid<<endl;
        resid=(resid-1)%RESID_MAX+1;
      }
      res_one=cbs->init->CHAIN_ID_LIST[cid]+"/"+to_string(resid)+"/"+stoken[2];
      Vresnames_new.push_back(res_one);
    }

    map<string,string> Mresnames;
    int nres=0;
    for(size_t jj=0;jj<Vrnas[ii].chains.size();jj++) {
      for(size_t kk=0;kk<Vrnas[ii].chains[jj].residues.size();kk++) {
        string res_one=Vrnas[ii].chains[jj].residues[kk].get_residue_name();
        Mresnames.emplace(res_one,Vresnames_new[nres]);
        Vrnas[ii].chains[jj].residues[kk].reset_residue_name(Vresnames_new[nres]);
        if(kk==0) {
          vector<string> stoken1=split_string(Vresnames_new[nres],"/");
          Vrnas[ii].chains[jj].chainid=stoken1[0];
        }
        nres++;
      }
    }
    Vrnas[ii].update_residue_name();

    vector<string> Vresname_temp;
    for(size_t jj=0;jj<Vtemplates[ii].templates.size();jj++) {
      for(size_t kk=0;kk<Vtemplates[ii].templates[jj].ends.size();kk++) {
        string end_one=Vtemplates[ii].templates[jj].ends[kk];
        Vtemplates[ii].templates[jj].ends[kk]=Mresnames[end_one];
      }
      string temp_info=Vtemplates[ii].templates[jj].get_RNA_motif_info();
      string type=Vtemplates[ii].templates[jj].type;
      if(type=="Loop" || type=="Loop_kiss_loop") {
        string nt_b_name=Vtemplates[ii].templates[jj].ends[0];
        string nt_e_name=Vtemplates[ii].templates[jj].ends[1];
        Residue nt_b=Vrnas[ii].find_residue(nt_b_name);
        Residue nt_e=Vrnas[ii].find_residue(nt_e_name);
        vector<double> xyz1=nt_b["C4'"].get_atom_xyz();
        vector<double> xyz2=nt_e["C4'"].get_atom_xyz();
        double e2e_dis=(xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0])+(xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1])+(xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2]);
        e2e_dis=sqrt(e2e_dis);
        temp_info+=" "+to_string(e2e_dis);
      }
      outT<<PDBid<<" "<<temp_info<<endl;

      if(type!="Helix") {
        vector<string> Vres=Vtemplates[ii].templates[jj].Vresnames;
        for(auto res : Vres) {
          string res_new=Mresnames[res];
          if(count(Vresname_temp.begin(),Vresname_temp.end(),res_new)==0)
            Vresname_temp.push_back(res_new);
        }
      }
    }
      
    for(size_t jj=0;jj<Vtemplates[ii].extended_templates.size();jj++) {
      for(size_t kk=0;kk<Vtemplates[ii].extended_templates[jj].ends.size();kk++){
        string end_one=Vtemplates[ii].extended_templates[jj].ends[kk];
        Vtemplates[ii].extended_templates[jj].ends[kk]=Mresnames[end_one];
      }
      string temp_info=Vtemplates[ii].extended_templates[jj].get_RNA_motif_info();
      string type=Vtemplates[ii].extended_templates[jj].type;
      if(type=="Loop" || type=="Loop_kiss_loop") {
        string nt_b_name=Vtemplates[ii].extended_templates[jj].ends[0];
        string nt_e_name=Vtemplates[ii].extended_templates[jj].ends[1];
        Residue nt_b=Vrnas[ii].find_residue(nt_b_name);
        Residue nt_e=Vrnas[ii].find_residue(nt_e_name);
        vector<double> xyz1=nt_b["C4'"].get_atom_xyz();
        vector<double> xyz2=nt_e["C4'"].get_atom_xyz();
        double e2e_dis=(xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0])+(xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1])+(xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2]);
        e2e_dis=sqrt(e2e_dis);
        temp_info+=" "+to_string(e2e_dis);
      }
      outT<<PDBid<<" "<<temp_info<<endl;

      if(type!="Helix") {
        vector<string> Vres=Vtemplates[ii].extended_templates[jj].Vresnames;
        for(auto res : Vres) {
          string res_new=Mresnames[res];
          if(count(Vresname_temp.begin(),Vresname_temp.end(),res_new)==0)
            Vresname_temp.push_back(res_new);
        }
      }
    } 
  
    fstr=outDIR+"pdb/"+PDBid+".pdb";
    ofstream outP(fstr.c_str());
    for(size_t jj=0;jj<Vrnas[ii].chains.size();jj++) {
      int check=0;
      for(size_t kk=0;kk<Vrnas[ii].chains[jj].size();kk++) {
        string res_one=Vrnas[ii].chains[jj].residues[kk].get_residue_name();
        if(count(Vresname_temp.begin(),Vresname_temp.end(),res_one)) {
          check=1;

          outP<<Vrnas[ii].chains[jj].residues[kk];
        }
      }
      if(check) outP<<"TER  "<<endl;
    }
    outP.close(); outP.clear();
  }
  outT.close();outT.clear();

  fstr=outDIR+"pdb_list.dat";
  ofstream outF (fstr.c_str());
  for(auto pdb : PDB_list)
    outF<<pdb<<endl;
  outF.close(); outF.clear();
  
  Vpdbs.clear();
  VVresnames_old.clear();
  Vrnas.clear();
  Vtemplates.clear();
}

void CBS_NS::check_template_list (string &outDIR) {
  string fstr=outDIR+"pdb_list.dat";
  if(!is_file_exist(fstr)) return;
  vector<string> Vlines=read_file(fstr);
  vector<string> pdbids;
  for(size_t ii=0;ii<Vlines.size();ii++) {
    vector<string> stoken=split_string(Vlines[ii]);
    string pdbid=stoken[0];
    if(count(pdbids.begin(),pdbids.end(),pdbid)==0)
    pdbids.push_back(pdbid);
  }

  fstr=outDIR+"template_list.dat";
  if(!is_file_exist(fstr)) return;
  Vlines.clear();
  Vlines=read_file(fstr);
  vector<string> template_list;
  for(size_t ii=0;ii<Vlines.size();ii++) {
    vector<string> stoken=split_string(Vlines[ii]);
    string pdbid=stoken[0];
    string template_name=stoken[0];
    for(size_t jj=1;jj<stoken.size();jj++)
      template_name+=" "+stoken[jj];
    if(count(pdbids.begin(),pdbids.end(),pdbid)==0)
      pdbids.push_back(pdbid);
    if(count(template_list.begin(),template_list.end(),template_name)==0)
      template_list.push_back(template_name);
  }

  ofstream outF (fstr.c_str());
  for(auto sline : template_list)
    outF<<sline<<endl;
  outF.close(); outF.clear();

  fstr=outDIR+"pdb_list.dat";
  outF.open(fstr.c_str());
  for(auto sline : pdbids)
    outF<<sline<<endl;
  outF.close(); outF.clear();
}

int CBS_NS::read_pdb_xyz(string fstr, ReadPDB &readpdb_one, vector<Molecule> &Vmoles, int nth) {
  readpdb_one.read_PDB(fstr);
  Vmoles[nth].set_molecule(readpdb_one.RNA_atoms[0],Chain::rna);
  
  return 0;
}

int CBS_NS::read_template_xyz(const Molecule &mole_one, string temp_info, vector<vector<double> > &VVxyz, int nth) {
  vector<double> Vxyz;
  vector<string> stoken=split_string(temp_info);
  string PDBid=stoken[0];
  string type=stoken[1];
  string str_size=stoken[2];
  string str_end=stoken[3];
  vector<string> Vends;
  stoken=split_string(str_end,"-&");
  for(size_t ii=0;ii<stoken.size();ii++)
    Vends.push_back(stoken[ii]);
  vector<int> lengths;
  stoken=split_string(str_size,"-&");
  for(size_t ii=0;ii<stoken.size();ii++)
    lengths.push_back(stoi(stoken[ii]));
  if(type=="Hairpin_kiss" || type=="Loop_kiss_loop")
    lengths[0]=lengths[0]+lengths[1]+lengths[2];
  else if(type=="Tail5" || type=="Tail3")
    lengths[0]=lengths[0]-1;
  
  size_t nchain=Vends.size()/2;
  for(size_t ii=0;ii<nchain;ii++) {
    string residue_begin=Vends[ii*2];
    string residue_end=Vends[ii*2+1];
    vector<string> stoken=split_string(residue_begin,"/");
    string chainid_one=stoken[0];
    int check=0;
    for(size_t jj=0;jj<mole_one.size();jj++) {
      for(int kk=0;kk<(int)mole_one.chains[jj].size();kk++) {
        if(mole_one.chains[jj].chainid!=chainid_one) continue;
        if(mole_one.chains[jj].residues[kk].get_residue_name()==residue_begin) {
          if(kk+lengths[ii]+1>=(int)mole_one.chains[jj].size()) break;
          else if(mole_one.chains[jj].residues[kk+lengths[ii]+1].get_residue_name()!=residue_end) continue;
          else {
            vector<double> xyz=mole_one.chains[jj].residues[kk].get_residue_xyz();
            Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
            check=1;
            continue;
          } 
        }
        if(check==1 && mole_one.chains[jj].residues[kk].get_residue_name()==residue_end) {
          vector<double> xyz=mole_one.chains[jj].residues[kk].get_residue_xyz();
          Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
          check=2;
          break;
        }
        if(check==1) {
          vector<double> xyz=mole_one.chains[jj].residues[kk].get_residue_xyz();
          Vxyz.insert(Vxyz.end(),xyz.begin(),xyz.end());
        }
      }
      if(check==2) break;
    }
    if(check==0) {
      string str="read_template_xyz: Cannot find residues "+residue_begin+"-"+residue_end+" in PDB "+PDBid;
      cout<<str<<endl;
      exit(0);
    }
  }
  
  VVxyz[nth]=Vxyz;

  return 0;
}

/*
int CBS_NS::read_template_xyz(ReadRNATemplate &readRT_one, string template_info, vector<vector<double> > &VVxyz, int nth) {
  VVxyz[nth]=readRT_one.read_template_xyz(template_info);
  return 0;
}
*/