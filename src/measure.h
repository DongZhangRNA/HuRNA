#ifndef CBS_MEASURE_H
#define CBS_MEASURE_H

#include "struct_definition.h"
#include "atom.h"

namespace CBS_NS {

class Measure {
 public:
  static const int equal,contained,intersected,separated;

  Measure();
  ~Measure();

  static double distance(double [3], double [3]);
  static double distance(const vector<double> &, const vector<double> &);
  static double distance(const Atom &, const Atom &);
  static double angle(double [3], double [3], double [3]);
  static double angle(const vector<double> &, const vector<double> &, const vector<double> &);
  static double angle(const Atom &, const Atom &, const Atom &);
  static double dihedral(double [3], double [3], double [3], double [3]);
  static double dihedral(const vector<double> &, const vector<double> &, const vector<double> &, const vector<double> &);
  static double dihedral(const Atom &, const Atom &, const Atom &, const Atom &);
  static double improper(double [3], double [3], double [3], double [3]);
  static double improper(const vector<double> &, const vector<double> &, const vector<double> &, const vector<double> &);
  static double improper(const Atom &, const Atom &, const Atom &, const Atom &);
  static double point_to_line_distance(double [3], double [3], double [3]);
  static double point_to_line_distance(const vector<double> &, const vector<double> &, const vector<double> &);
  static pair<double,double> statistics(const vector<double> &);
  static pair<double,double> statistics(const vector<int> &);
  static double average(const vector<double> &);
  static double average(const vector<int> &);
  static double standard_deviation(const vector<double> &);
  static double standard_deviation(const vector<int> &);
  static double standard_deviation(const vector<double> &, double);
  static double standard_deviation(const vector<int> &, double); 
  static vector<double> probability_distribution(const vector<double> &, double);
  static vector<double> probability_distribution(const vector<double> &, double, double, double);
  static vector<double> probability_distribution_limited(const vector<double> &, double, double, double);
  static vector<double> probability_distribution(const vector<int> &, double);
  static vector<double> probability_distribution(const vector<int> &, double, double, double);
  static vector<double> probability_distribution_limited(const vector<int> &, double, double, double);
  static double Pearson_correlation_coefficient(const vector<double> &, const vector<double> &);
  static double Pearson_correlation_coefficient(const vector<int> &, const vector<int> &);
  static double Spearman_rank_correlation_coefficient(const vector<double> &, const vector<double> &);
  static double Spearman_rank_correlation_coefficient(const vector<int> &, const vector<int> &);
  static pair<double,double> linear_fitting(const vector<double> &, const vector<double> &);
  static pair<double,double> linear_fitting(const vector<int> &, const vector<int> &);
  static double root_mean_square_error(const vector<double> &, const vector<double> &);
  static double root_mean_square_error(const vector<int> &, const vector<int> &);
  static double cross_entropy(const vector<double> &, const vector<double> &);
  static double relative_entropy(const vector<double> &, const vector<double> &);
  static int vector_relation(vector<string>, vector<string>);
  static int vector_relation(vector<int>, vector<int>);
  static ConfusionMatrixScore confusion_matrix_score(const vector<string> &, const vector<string> &);
};

}

#endif