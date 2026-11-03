#ifndef CBS_RANDOM_H
#define CBS_RANDOM_H

namespace CBS_NS {

class Random {
 public:
  Random();

  double uniform();
  double gaussian();
  void reset(int);

 private:
  int seed,save;
  double second;
};

}

#endif