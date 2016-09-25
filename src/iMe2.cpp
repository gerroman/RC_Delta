#include "constants.h"
#include "proton_formfactors.h"
#include <cmath>

double iMe2(double nu, double, double qp2, double kK, double kP, double kq) {
return (e6*Z2*pow(kK - kq,-2)*pow(kK + kq,-2)*pow(qp2,-2)*pow(-4*Mp2 + qp2,-1)*
  (8*kq*(4*Mp2*(4*kK*kP*me2*nu + (4*Mp2 - qp2)*qp2*(pow(kK,2) - pow(kq,2)))*
        pow(GE(qp2),2) + qp2*(-4*kK*kP*me2*nu + 
          (4*Mp2 - qp2)*(4*me2 + qp2)*(pow(kK,2) - pow(kq,2)))*pow(GM(qp2),2)\
) + 4*(4*Mp2*(4*me2*pow(kK,2)*pow(kP,2) + 
          qp2*pow(kP,2)*(pow(kK,2) - pow(kq,2)) + 
          (4*Mp2 - qp2)*(pow(kK,4) - pow(kq,4)))*pow(GE(qp2),2) + 
       qp2*(-4*me2*pow(kK,2)*pow(kP,2) - 
          qp2*pow(kP,2)*(pow(kK,2) - pow(kq,2)) + 
          (4*Mp2 - qp2)*(pow(kK,4) - pow(kq,4)))*pow(GM(qp2),2)) - 
    4*(qp2*(pow(kK,2) - pow(kq,2)) + 4*me2*pow(kq,2))*
     (-4*Mp2*((4*Mp2 - qp2)*qp2 + pow(nu,2))*pow(GE(qp2),2) + 
       qp2*(-((4*Mp2 - qp2)*(4*me2 + qp2)) + pow(nu,2))*pow(GM(qp2),2)))
);
}

