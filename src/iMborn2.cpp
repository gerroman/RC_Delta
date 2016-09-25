#include "iMborn2.h"
#include "constants.h"
#include "proton_formfactors.h"
#include <cmath>

double iMborn2(double nu, double q2) {
  return Z2 * e4 / pow(q2, 2) * (
       4 * Mp2 * pow(GE(q2), 2) * (pow(nu, 2) + q2 * (4*Mp2 - q2))
        - q2 * pow(GM(q2), 2) * (pow(nu, 2) - (4*me2 + q2) * (4*Mp2 - q2))
      ) / (4 * Mp2 - q2);
}
