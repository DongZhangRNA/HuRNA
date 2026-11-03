#ifndef CBS_CONSTANT_H
#define CBS_CONSTANT_H

namespace CBS_NS {

namespace MathConst {
  const double PI  = 3.14159265358979323846;
  const double SMALL = 0.001;
  const double EPS = 1E-6;
  const double BIG = 1E20;

  const double KBT = 0.596; //T=300 K (kcal/mol)
  const double KBT_310 = 0.616; //T=310K (kcal/mol)
}

namespace ResidueConst {
  const string RESIDUE_NULL = "RESNULL";
}

namespace MassConst {
  const double MASS_C = 12.011;
  const double MASS_N = 14.007;
  const double MASS_O = 15.999;
  const double MASS_S = 32.065;
}

}

#endif