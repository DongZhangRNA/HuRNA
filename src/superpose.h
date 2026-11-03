#ifndef CBS_SUPERPOSE_H
#define CBS_SUPERPOSE_H

namespace CBS_NS {

class Superpose {
 public:
  vector<vector<double> > M;
  vector<double> x1cen,x2cen;

  Superpose();
  Superpose(const vector<double> &, const vector<double> &);
  ~Superpose();

  void set_x1x2(const vector<double> &, const vector<double> &);
  void set_x1(const vector<double> &);
  void set_x2(const vector<double> &);
  void calculate_rotation_matrix();
  double get_RMSD();
  void reset_x2(const vector<double> &);
  vector<double> rotate_x2();
  vector<double> rotate_x2(const vector<double> &);
  void clear();

 private:
  int lenx;
  vector<vector<double> > x1,x2;

  //inline double RMSD_rmsdx();
};

}

#endif