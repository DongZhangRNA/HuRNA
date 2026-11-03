#include "basic.h"
#include "constant.h" 
#include "initialization.h"
#include "error.h"
#include "readfile.h"
#include "constant.h"
#include "sort_vector.h"
#include "cluster.h"

using namespace CBS_NS;

const int NOUT = 2000;
const float RMSD_SMALL = 0.01;
const double SCALE_R2 = 0.5;

Cluster::Cluster(Cobmms *cbs) : Pointers(cbs) {}

Cluster::Cluster(Cobmms *cbs, const vector<vector<double> > & Vxyz_one) : Pointers(cbs) {
  Cluster::input_data(Vxyz_one); 
}

Cluster::Cluster(const Cluster &other) : Pointers(other.cbs) {
  Cluster::clear();
  Vxyz=other.Vxyz;
  clusters=other.clusters;
  RMSD_matrix=other.RMSD_matrix;
  RMSD_cutoff=other.RMSD_cutoff;
  Nmember=other.Nmember;
}

Cluster& Cluster::operator =(const Cluster &other) {
  if(this != &other) {
    Pointers::operator=(other.cbs);
    Vxyz=other.Vxyz;
    clusters=other.clusters;
    RMSD_matrix=other.RMSD_matrix;
    RMSD_cutoff=other.RMSD_cutoff;
    Nmember=other.Nmember;
  }

  return *this;
}

Cluster::~Cluster() {
  Vxyz.clear();
  Vxyz.shrink_to_fit();
  clusters.clear();
  clusters.shrink_to_fit();
  RMSD_matrix.clear();
  RMSD_matrix.shrink_to_fit();
  
  RMSD_cutoff=0.0;
  Nmember=0;
}

void Cluster::input_data(const vector<vector<double> > & Vxyz_one) {
  Vxyz=Vxyz_one;

  Nmember=Vxyz.size();
  size_t ndata=Vxyz[0].size();
  for(size_t ii=1;ii<Nmember;ii++) {
    if(Vxyz[ii].size()!=ndata) {
      string str="cluster.cpp: Numbers of data are not match  --> "+to_string(ndata);
      str+=" : "+to_string(Vxyz[ii].size())+" in member #"+to_string(ii);
      error->one(str);
    }
  }

  RMSD_matrix.clear();
  for(size_t ii=0;ii<Nmember-1;ii++) {
    vector<float> rmsdij;
    for(size_t jj=ii+1;jj<Nmember;jj++)
      rmsdij.push_back(0.0);
    RMSD_matrix.push_back(rmsdij);
  }
}

void Cluster::input_data(const vector<vector<Atom> > & Vatom_one) {
  Vxyz.clear();
  Nmember=Vatom_one.size();
  for(size_t ii=0;ii<Nmember;ii++) {
    vector<double> xyz;
    for(size_t jj=0;jj<Vatom_one[ii].size();jj++) {
      vector<double> xyz_one=Vatom_one[ii][jj].get_atom_xyz();
      xyz.insert(xyz.end(),xyz_one.begin(),xyz_one.end());
    }
    Vxyz.push_back(xyz);
  }

  size_t ndata=Vxyz[0].size();
  for(size_t ii=1;ii<Nmember;ii++) {
    if(Vxyz[ii].size()!=ndata) {
      string str="cluster.cpp: Numbers of data are not match  --> "+to_string(ndata);
      str+=" : "+to_string(Vxyz[ii].size())+" in member #"+to_string(ii);
      error->one(str);
    }
  }

  RMSD_matrix.clear();
  for(size_t ii=0;ii<Nmember-1;ii++) {
    vector<float> rmsdij;
    for(size_t jj=ii+1;jj<Nmember;jj++)
      rmsdij.push_back(0.0);
    RMSD_matrix.push_back(rmsdij);
  }
}

void Cluster::calculate_RMSD_matrix() {
  int nthread=init->nthread;
  map<size_t,size_t> Midentical;
  for(size_t ii=0;ii<Nmember;ii++) {
    if(ii>0 && ii%NOUT==0) {
      string str="      Calculating pairwise RMSDs: "+to_string(ii)+"/"+to_string(Nmember);
      error->log(str);
    }

    if(Midentical.count(ii)==1) {
      for(size_t jj=ii+1;jj<Nmember;jj++)
        RMSD_matrix[ii][jj-ii-1]=RMSD_matrix[Midentical[ii]][jj-Midentical[ii]-1];
      continue;
    }

    //for(size_t jj=0;jj<ii;jj++)
      //RMSD_matrix[ii][jj]=RMSD_matrix[jj][ii];
    int ni_member=Nmember-ii-1;
    vector<Superpose> Vsp;
    int nsp=nthread;
    if(nsp>ni_member)
      nsp=ni_member;
    for(int jj=0;jj<nsp;jj++) {
      Superpose sp_one;
      sp_one.set_x1(Vxyz[ii]);
      Vsp.push_back(sp_one);
    }
    int nloop=ni_member/nthread;
    int nlast=ni_member%nthread;
    for(int nl=0;nl<nloop;nl++) {
      int nl0=ii+1+nl*nthread;
      vector<vector<double> > Vxyz_loop;
      Vxyz_loop.resize(nthread);
      for(int kk=0;kk<nthread;kk++)
        Vxyz_loop[kk]=Vxyz[nl0+kk];
      vector<double> rmsd_loop;
      rmsd_loop.resize(nthread);
      vector<std::thread> threads;
      for(int ni=0;ni<nthread;ni++)
        threads.emplace_back(std::thread(Cluster::calculate_pairwise_RMSD,std::ref(Vsp[ni]),std::ref(Vxyz_loop[ni]),std::ref(rmsd_loop),ni));
      for(int ni=0;ni<nthread;ni++)
        threads[ni].join();
      for(int kk=0;kk<nthread;kk++)
        RMSD_matrix[ii][nl0+kk-ii-1]=static_cast<float>(rmsd_loop[kk]);
    }
    if(nlast!=0) {
      int nl0=ii+1+nloop*nthread;
      vector<vector<double> > Vxyz_last;
      Vxyz_last.resize(nlast);
      for(int kk=0;kk<nlast;kk++)
        Vxyz_last[kk]=Vxyz[nl0+kk];
      vector<double> rmsd_last;
      rmsd_last.resize(nlast);
      vector<std::thread> threads;
      for(int ni=0;ni<nlast;ni++)
        threads.emplace_back(std::thread(Cluster::calculate_pairwise_RMSD,std::ref(Vsp[ni]),std::ref(Vxyz_last[ni]),std::ref(rmsd_last),ni));
      for(int ni=0;ni<nlast;ni++)
        threads[ni].join();
      for(int kk=0;kk<nlast;kk++)
        RMSD_matrix[ii][nl0+kk-ii-1]=static_cast<float>(rmsd_last[kk]);
    }

    for(size_t jj=ii+1;jj<Nmember;jj++) {  //if rmsdij<RMSD_SMALL, treat i&j as identical structures, set rmsd_ij(j,n)=rmsd_ij(i,n)
      if(RMSD_matrix[ii][jj-ii-1]<RMSD_SMALL)
        Midentical.emplace(jj,ii);
    }
  }
}

void Cluster::calculate_RMSD_matrix(string &fstr) {
  size_t n_skip=0;
  if(is_file_exist(fstr)) {
    vector<string> Vlines=read_file(fstr);
    for(size_t ii=0;ii<Vlines.size();ii++) {
      vector<string> stoken=split_string(Vlines[ii]);
      if(ii==0 && stoken.size()!=Vlines.size()) {
        string str="cluster.cpp: Wrong data for RMSD_matrix in "+fstr;
        error->message(str);
        break;
      }
      else
        n_skip=Vlines.size();
      for(size_t jj=0;jj<stoken.size();jj++) {
        float rmsd_one=stof(stoken[jj]);
        RMSD_matrix[ii][jj]=rmsd_one;
        //RMSD_matrix[ii][ii+jj+1]=rmsd_one;
        //RMSD_matrix[ii+jj+1][ii]=rmsd_one;
      } 
    }

    //string cmd="mv "+fstr+" "+fstr+".bak";
    //cbs->execute_system_command(cmd);
  }
  
  int nthread=init->nthread;
  map<size_t,size_t> Midentical;
  //cout<<"    "<<"Number of thread: "<<nthread<<endl;
  for(size_t ii=0;ii<Nmember;ii++) {
    if(ii>0 && ii%NOUT==0) {
      string str="      Calculating pairwise RMSDs: "+to_string(ii)+"/"+to_string(Nmember);
      error->log(str);
    }

    if(Midentical.count(ii)==1) {
      for(size_t jj=ii+1;jj<Nmember;jj++)
        RMSD_matrix[ii][jj-ii-1]=RMSD_matrix[Midentical[ii]][jj-Midentical[ii]-1];
      continue;
    }

    int ni_member=Nmember-ii-1;
    int ni_start=ii+1;
    if(ii<=n_skip) {
      ni_member=Nmember-1-n_skip;
      ni_start=n_skip+1;
    }
    else {
      //for(size_t jj=0;jj<ii;jj++)
        //RMSD_matrix[ii][jj]=RMSD_matrix[jj][ii]; 
    }
    vector<Superpose> Vsp;
    int nsp=nthread;
    if(nsp>ni_member)
      nsp=ni_member;
    for(int jj=0;jj<nsp;jj++) {
      Superpose sp_one;
      sp_one.set_x1(Vxyz[ii]);
      Vsp.push_back(sp_one);
    }
    int nloop=ni_member/nthread;
    int nlast=ni_member%nthread;
    for(int nl=0;nl<nloop;nl++) {
      int nl0=ni_start+nl*nthread;
      vector<vector<double> > Vxyz_loop;
      Vxyz_loop.resize(nthread);
      for(int kk=0;kk<nthread;kk++)
        Vxyz_loop[kk]=Vxyz[nl0+kk];
      vector<double> rmsd_loop;
      rmsd_loop.resize(nthread);
      vector<std::thread> threads;
      for(int ni=0;ni<nthread;ni++)
        threads.emplace_back(std::thread(Cluster::calculate_pairwise_RMSD,std::ref(Vsp[ni]),std::ref(Vxyz_loop[ni]),std::ref(rmsd_loop),ni));
      for(int ni=0;ni<nthread;ni++)
        threads[ni].join();
      for(int kk=0;kk<nthread;kk++)
        RMSD_matrix[ii][nl0+kk-ii-1]=static_cast<float>(rmsd_loop[kk]);
    }
    if(nlast!=0) {
      int nl0=ni_start+nloop*nthread;
      vector<vector<double> > Vxyz_last;
      Vxyz_last.resize(nlast);
      for(int kk=0;kk<nlast;kk++)
        Vxyz_last[kk]=Vxyz[nl0+kk];
      vector<double> rmsd_last;
      rmsd_last.resize(nlast);
      vector<std::thread> threads;
      for(int ni=0;ni<nlast;ni++)
        threads.emplace_back(std::thread(Cluster::calculate_pairwise_RMSD,std::ref(Vsp[ni]),std::ref(Vxyz_last[ni]),std::ref(rmsd_last),ni));
      for(int ni=0;ni<nlast;ni++)
        threads[ni].join();
      for(int kk=0;kk<nlast;kk++)
        RMSD_matrix[ii][nl0+kk-ii-1]=static_cast<float>(rmsd_last[kk]);
    }
    
    for(size_t jj=ii+1;jj<Nmember;jj++) {  //if rmsdij<RMSD_SMALL, treat i&j as identical structures, set rmsd_ij(j,n)=rmsd_ij(i,n)
      if(RMSD_matrix[ii][jj-ii-1]<RMSD_SMALL)
        Midentical.emplace(jj,ii);
    }
  }

  ofstream outF (fstr.c_str());
  for(size_t ii=0;ii<Nmember-1;ii++) {
    for(size_t jj=ii+1;jj<Nmember;jj++)
      outF<<std::fixed<<std::setprecision(2)<<RMSD_matrix[ii][jj-ii-1]<<" ";
    outF<<endl;
  } 
  outF.close(); outF.clear(); 
}

void Cluster::get_RMSD_cutoff(double percentage_one) {
  int ncut_one=round(Nmember*percentage_one);
  double rmsd_cutoff_one=MathConst::BIG;
  for(size_t ii=0;ii<Nmember;ii++) {
    vector<float> rmsdi;
    for(size_t jj=0;jj<Nmember;jj++)
      rmsdi.push_back(rmsd_ij(ii,jj));
    sort(rmsdi.begin(),rmsdi.end());
    double rmsd=rmsdi[ncut_one-1]*1.0;
    if(rmsd<rmsd_cutoff_one)
      rmsd_cutoff_one=rmsd;
  }

  RMSD_cutoff=rmsd_cutoff_one;
}

void Cluster::get_clusters(double rmsd_cutoff_one) {
  clusters.clear();
  RMSD_cutoff=rmsd_cutoff_one;

  size_t cluster_center_index;
  int cluster_size=0;
  for(size_t ii=0;ii<Nmember;ii++) {
    int num=0;
    for(size_t jj=0;jj<Nmember;jj++) {
      if(Cluster::rmsd_ij(ii,jj)<=rmsd_cutoff_one)
        num++;
    }
    if(num>cluster_size) {
      cluster_size=num;
      cluster_center_index=ii;
    }
  }
  vector<size_t> Vids;
  Vids.push_back(cluster_center_index);
  for(size_t ii=0;ii<Nmember;ii++) {
    if(ii==cluster_center_index) continue;
    if(Cluster::rmsd_ij(cluster_center_index,ii)<=rmsd_cutoff_one)
      Vids.push_back(ii);
  }
  clusters.push_back(Vids);   

  map<size_t,vector<size_t> > Mneighbors;
  cluster_center_index=0;
  cluster_size=0;
  vector<size_t> Vsingle;
  for(size_t ii=0;ii<Nmember;ii++) {
    if(count(Vids.begin(),Vids.end(),ii)) continue;
    vector<size_t> Vneighbor;
    for(size_t jj=0;jj<Nmember;jj++) {
      if(jj==ii) continue;
      if(count(Vids.begin(),Vids.end(),jj)) continue;
      if(Cluster::rmsd_ij(ii,jj)<=rmsd_cutoff_one)
        Vneighbor.push_back(jj);
    }
    int num=(int)Vneighbor.size();
    if(num==0) {
      Vsingle.push_back(ii);
    }
    else {
      Mneighbors.emplace(ii,Vneighbor);
      if(num>cluster_size) {
        cluster_center_index=ii;
        cluster_size=num;
      }
    } 
  }

  while(Mneighbors.size()>0) {
    Vids=Mneighbors[cluster_center_index];
    Vids.insert(Vids.begin(),cluster_center_index);
    clusters.push_back(Vids);
    for(auto id : Vids) {
      Mneighbors.erase(id);
    }

    cluster_center_index=0;
    cluster_size=0;
    vector<size_t> Vsingle_new;
    for(auto& it : Mneighbors) {
      vector<size_t> Vids_new;
      for(auto id : it.second) {
        if(count(Vids.begin(),Vids.end(),id)==0)
          Vids_new.push_back(id);
      }
      if(Vids_new.size()==0) {
        Vsingle.push_back(it.first);
        Vsingle_new.push_back(it.first);
        continue;
      }
      if((int)Vids_new.size()>cluster_size) {
        cluster_size=Vids_new.size();
        cluster_center_index=it.first;
      }
      Mneighbors[it.first]=Vids_new;
    }
    if(Vsingle_new.size()>0) {
      for(auto id : Vsingle_new)
        Mneighbors.erase(id);
    }
  }
  for(auto id : Vsingle) {
    Vids.clear();
    Vids.push_back(id);
    clusters.push_back(Vids);
  }
}

void Cluster::get_clusters(double rmsd_cutoff_one, double percentage_one, double percentage_min) {
  int ncut_one=round(Nmember*percentage_one);
  int ncut_min=round(Nmember*percentage_min);
  double rmsd_one_min=MathConst::BIG, rmsd_min_min=MathConst::BIG;
  for(size_t ii=0;ii<Nmember;ii++) {
    vector<float> rmsdi;
    for(size_t jj=0;jj<Nmember;jj++)
      rmsdi.push_back(rmsd_ij(ii,jj));
    std::sort(rmsdi.begin(),rmsdi.end());
    double rmsd_one=rmsdi[ncut_one-1]*1.0;
    double rmsd_min=rmsdi[ncut_min-1]*1.0;
    if(rmsd_one<rmsd_one_min)
      rmsd_one_min=rmsd_one;
    if(rmsd_min<rmsd_min_min)
      rmsd_min_min=rmsd_min;
  }

  if(rmsd_one_min<rmsd_cutoff_one)
    rmsd_cutoff_one=rmsd_one_min;
  if(rmsd_min_min>rmsd_cutoff_one)
    rmsd_cutoff_one=rmsd_min_min;
  
  Cluster::get_clusters(rmsd_cutoff_one);
}

//Cluster big dataset through two rounds
void Cluster::cluster_big_dataset(const vector<vector<double> > &Vxyz_full, double rmsd_cutoff_one, int Nevery) {
  //Round #1 clustering
  vector<vector<double> > Vxyz_one;
  int nloop=Vxyz_full.size()/Nevery;
  int nlast=Vxyz_full.size()%Nevery;
  int Ncal=nloop;
  if(nlast>0) Ncal+=1;
  vector<vector<size_t> >Vclusters_R1;
  for(int ii=0;ii<nloop;ii++) {
    string str="cluster.cpp (Round #1): "+to_string(ii+1)+"/"+to_string(Ncal);
    error->log(str);
    size_t index_offset=ii*Nevery;
    Vxyz_one.clear();
    for(int jj=0;jj<Nevery;jj++)
      Vxyz_one.push_back(Vxyz_full[jj+index_offset]);
    Cluster::input_data(Vxyz_one);
    Cluster::calculate_RMSD_matrix();
    Cluster::get_clusters(rmsd_cutoff_one);
    for(size_t jj=0;jj<clusters.size();jj++) {
      vector<size_t> cluster_one;
      for(size_t kk=0;kk<clusters[jj].size();kk++)
        cluster_one.push_back(clusters[jj][kk]+index_offset);
      Vclusters_R1.push_back(cluster_one);
    }
  }
  if(nlast>0) {
    string str="cluster.cpp (Round #1): "+to_string(Ncal)+"/"+to_string(Ncal);
    error->log(str);
    size_t index_offset=nloop*Nevery;
    Vxyz_one.clear();
    for(int jj=0;jj<nlast;jj++)
      Vxyz_one.push_back(Vxyz_full[jj+index_offset]);
    Cluster::input_data(Vxyz_one);
    Cluster::calculate_RMSD_matrix();
    Cluster::get_clusters(rmsd_cutoff_one);
    for(size_t jj=0;jj<clusters.size();jj++) {
      vector<size_t> cluster_one;
      for(size_t kk=0;kk<clusters[jj].size();kk++)
        cluster_one.push_back(clusters[jj][kk]+index_offset);
      Vclusters_R1.push_back(cluster_one);
    }
  }

  //Round #2 clustering
  string str="cluster.cpp (Round #2): Clustering "+to_string(Vclusters_R1.size())+" center structures...";
  error->log(str);
  
  Vxyz_one.clear();
  for(size_t ii=0;ii<Vclusters_R1.size();ii++) {
    Vxyz_one.push_back(Vxyz_full[Vclusters_R1[ii][0]]);
  }
  Cluster::input_data(Vxyz_one);
  Cluster::calculate_RMSD_matrix();
  Cluster::get_clusters(rmsd_cutoff_one*SCALE_R2);

  //Merge Round #1 and #2 clustering
  vector<vector<size_t> > Tclusters;
  for(size_t ii=0;ii<clusters.size();ii++) {
    vector<size_t> cluster_i;
    for(size_t jj=0;jj<clusters[ii].size();jj++) {
      size_t index_j=clusters[ii][jj];
      cluster_i.insert(cluster_i.end(),Vclusters_R1[index_j].begin(),Vclusters_R1[index_j].end());
    }
    Tclusters.push_back(cluster_i);
  }
  vector<size_t> Vsize;
  for(size_t ii=0;ii<Tclusters.size();ii++) {
    size_t num=Tclusters[ii].size();
    Vsize.push_back(num);
  }
  vector<size_t> Vsort=SortVector::sort_vector_descending(Vsize);
  clusters.clear();
  for(size_t ii=0;ii<Vsort.size();ii++) {
    clusters.push_back(Tclusters[Vsort[ii]]);
  }

  /*
  vector<string> Vsizes;
  for(size_t ii=0;ii<Tclusters.size();ii++) {
    string str=to_string(ii)+" "+to_string(Tclusters[ii].size());
    Vsizes.push_back(str);
  }
  SortVector::sort_vector_descending(Vsizes,1);
  clusters.clear();
  for(size_t ii=0;ii<Vsizes.size();ii++) {
    vector<string> stoken=split_string(Vsizes[ii]);
    int index_i=stoi(stoken[0]);
    clusters.push_back(Tclusters[index_i]);
  }
  */
}

void Cluster::clear() {
  Vxyz.clear();
  clusters.clear();
  RMSD_matrix.clear();
  
  RMSD_cutoff=0.0;
  Nmember=0;
}

size_t Cluster::size() const {
  return clusters.size();
}

int Cluster::calculate_pairwise_RMSD(Superpose &sp_one, const vector<double> &Vtar, vector<double> &Vrmsd, int nth) {
  //cout<<"Thread #"<<nth<<" starts!"<<endl;
  sp_one.set_x2(Vtar);
  Vrmsd[nth]=sp_one.get_RMSD();
  //cout<<"Thread #"<<nth<<" ends!"<<endl;

  return 0;
}

inline float Cluster::rmsd_ij(size_t row_i, size_t col_j) {
  if(row_i==col_j) return 0.0;
  else if(row_i<col_j) return RMSD_matrix[row_i][col_j-row_i-1];
  else return RMSD_matrix[col_j][row_i-col_j-1];
}

ostream & operator <<(ostream &output, Cluster &cluster_one) {
  string str="Number_of_structures : "+to_string(cluster_one.Vxyz.size());
  output<<str<<endl;
  str="Number_of_clusters : "+to_string(cluster_one.size());
  output<<str<<endl;
  str="RMSD_cutoff_for_clustering : "+to_string(cluster_one.RMSD_cutoff);
  output<<str<<endl;
  output<<"----------"<<endl;
  size_t Nclu=cluster_one.size();
  if(Nclu>50) {
    Nclu=50;
    str="#Cluster(Top50) Number_of_members Percentage";
    output<<str<<endl;
  }
  else {
    str="#Cluster Number_of_members Percentage";
    output<<str<<endl;
  }
  for(size_t ii=0;ii<Nclu;ii++) {
    double per=(1.0*cluster_one.clusters[ii].size())/cluster_one.Vxyz.size();
    str="#"+to_string(ii+1)+" "+to_string(cluster_one.clusters[ii].size())+" "+to_string(per);
    output<<str<<endl;
  }

  return output;
}