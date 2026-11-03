#ifndef CBS_CLUSTER_H
#define CBS_CLUSTER_H

#include "pointers.h"
#include "atom.h"
#include "superpose.h"

namespace CBS_NS {class Cluster;}
ostream& operator <<(ostream &, CBS_NS::Cluster &);

namespace CBS_NS {

class Cluster : protected Pointers{
 public:
  vector<vector<double> > Vxyz;
  vector<vector<size_t> > clusters;
  vector<vector<float> > RMSD_matrix; //Using float type and upper triangle matrix to reduce memory cost
  double RMSD_cutoff;

  Cluster(Cobmms *);
  Cluster(Cobmms *, const vector<vector<double> > &);
  Cluster(const Cluster &);
  Cluster& operator=(const Cluster &);
  ~Cluster();

  void input_data(const vector<vector<double> > &);
  void input_data(const vector<vector<Atom> > &);
  void calculate_RMSD_matrix();
  void calculate_RMSD_matrix(string &);
  void get_RMSD_cutoff(double);
  void get_clusters(double);
  void get_clusters(double, double, double =0.01);
  void cluster_big_dataset(const vector<vector<double> > &, double, int = 10000);
  void clear();
  size_t size() const;
  friend ostream& ::operator <<(ostream &, Cluster &);

 private:
  size_t Nmember;
  static int calculate_pairwise_RMSD(Superpose &, const vector<double> &, vector<double> &, int); 
  inline float rmsd_ij(size_t, size_t);
};

}

#endif