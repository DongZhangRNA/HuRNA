#ifndef CBS_SORT_VECTOR_H
#define CBS_SORT_VECTOR_H

namespace CBS_NS {

class SortVector {
 public:

  SortVector();
  ~SortVector();
  
  static vector<size_t> sort_vector_ascending(vector<double> &);
  static vector<size_t> sort_vector_ascending(vector<int> &);
  static vector<size_t> sort_vector_ascending(vector<size_t> &);
  static void sort_vector_ascending(vector<string> &, size_t =0);
  static vector<size_t> sort_vector_descending(vector<double> &);
  static vector<size_t> sort_vector_descending(vector<int> &);
  static vector<size_t> sort_vector_descending(vector<size_t> &);
  static void sort_vector_descending(vector<string> &, size_t =0);
};

}

#endif