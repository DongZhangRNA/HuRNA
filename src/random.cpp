#include "basic.h"
#include "random.h"

using namespace CBS_NS;

const int IA = 16807;
const int IM = 2147483647;
const double AM = (1.0/IM);
const int IQ = 127773;
const int IR = 2836;

Random::Random() {
  srand((unsigned)time(NULL));
  seed=(int) rand();
  save=0;
}

double Random::uniform() {
  int k=seed/IQ;
  seed=IA*(seed-k*IQ)-IR*k;
  if(seed<0) seed+=IM;
  double ans=AM*seed;

  return ans;
}

double Random::gaussian() {
  double first,v1,v2,rsq,fac;

  if(!save) {
    int again=1;
    while(again) {
      v1=2.0*Random::uniform()-1.0;
      v2=2.0*Random::uniform()-1.0;
      rsq=v1*v1+v2*v2;
      if(rsq<1.0 && rsq!=0.0) again=0;
    }
    fac=sqrt(-2.0*log(rsq)/rsq);
    second=v1*fac;
    first=v2*fac;
    save=1;
  }
  else {
    first=second;
    save=0;
  }

  return first;
}

void Random::reset(int seed_init) {
  if(seed_init<=0) {
    cout<<"Error in random.cpp: Invalid seed for random # generator"<<endl;
    exit(0);
  }
  seed=seed_init;
  save=0;
}