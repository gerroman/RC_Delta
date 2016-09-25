#include "constants.h"
#include "proton_formfactors.h"
#include <cmath>

double iMp2(double nu, double qe2, double, double kK, double kP, double kq) {
return ((e6*Z4*pow(kP - kq,-2)*pow(kP + kq,-2)*pow(Mp,-4)*pow(qe2,-2)*
    (-8*Mp2*F1(qe2)*(F1(0) + F2(0))*F2(qe2)*(pow(kP,2) - pow(kq,2))*
       (4*Mp2*(4*me2*(F1(0) - F2(0))*pow(kq,2) + 
            qe2*((F1(0) - F2(0))*pow(kK,2) + F1(0)*pow(kq,2))) + 
         F2(0)*(qe2*(3*kK*kP*nu + 8*me2*(-pow(kP,2) + pow(kq,2))) + 
            pow(kq,2)*pow(nu,2) + 
            (2*pow(kK,2) - pow(kP,2) + pow(kq,2))*pow(qe2,2))) - 
      8*Mp2*pow(F1(qe2),2)*(8*Mp4*pow(kK,2)*pow(kP,2)*pow(F1(0),2) + 
         (pow(kP,2) - pow(kq,2))*
          (2*qe2*(kK*kP*nu + 2*me2*(-pow(kP,2) + pow(kq,2))) + 
            pow(kq,2)*pow(nu,2) + pow(kK,2)*pow(qe2,2))*pow(F2(0),2) + 
         2*Mp2*(pow(kP,2) - pow(kq,2))*
          (4*me2*(pow(kP,2)*pow(F1(0),2) + 
               pow(kq,2)*(2*F1(0)*F2(0) + pow(F1(0),2) - pow(F2(0),2))) + 
            qe2*(pow(kP,2)*pow(F1(0),2) + 
               pow(kK,2)*(2*F1(0)*F2(0) + pow(F1(0),2) - pow(F2(0),2)) + 
               pow(kq,2)*pow(F1(0) + F2(0),2)))) - 
      F2(0)*(pow(kP,2) - pow(kq,2))*
       (8*Mp2*F1(qe2)*(F1(0) + F2(0))*F2(qe2)*pow(kK,2)*pow(kP,2) + 
         8*Mp2*F2(0)*pow(kK,2)*pow(kP,2)*pow(F1(qe2),2) + 
         (4*Mp2*(2*F1(0) + F2(0))*pow(kK,2)*pow(kP,2) + 
            F2(0)*(pow(kP,2) - pow(kq,2))*
             (4*me2*(pow(kP,2) - pow(kq,2)) - 
               qe2*(pow(kK,2) - pow(kP,2) + pow(kq,2))))*pow(F2(qe2),2)) - 
      4*Mp2*(qe2*(pow(kP,2) - pow(kq,2)) + 4*Mp2*pow(kq,2))*pow(F1(0),2)*
       (16*Mp2*qe2*(2*me2 + qe2)*F1(qe2)*F2(qe2) + 
         4*Mp2*(qe2*(4*me2 + 4*Mp2 + qe2) + pow(nu,2))*pow(F1(qe2),2) + 
         qe2*(16*me2*Mp2 + 4*Mp2*qe2 - pow(nu,2) + pow(qe2,2))*
          pow(F2(qe2),2)) + 8*kq*Mp2*pow(F1(0),2)*
       (16*Mp2*qe2*(2*me2 + qe2)*F1(qe2)*F2(qe2)*(pow(kP,2) - pow(kq,2)) + 
         4*Mp2*(4*Mp2*(kK*kP*nu + qe2*(pow(kP,2) - pow(kq,2))) + 
            qe2*(4*me2 + qe2)*(pow(kP,2) - pow(kq,2)))*pow(F1(qe2),2) + 
         qe2*(-4*Mp2*(kK*kP*nu - (4*me2 + qe2)*(pow(kP,2) - pow(kq,2))) + 
            (pow(kP,2) - pow(kq,2))*pow(qe2,2))*pow(F2(qe2),2)) + 
      2*kq*(pow(kP,2) - pow(kq,2))*
       (4*Mp2*F1(qe2)*F2(0)*(F1(0) + F2(0))*F2(qe2)*
          (qe2*(3*pow(kK,2) - pow(kP,2) + pow(kq,2)) + 
            2*(kK*kP*nu + 6*me2*(-pow(kP,2) + pow(kq,2)))) + 
         8*Mp2*(kK*kP*nu + qe2*pow(kK,2) + 4*me2*(-pow(kP,2) + pow(kq,2)))*
          pow(F1(qe2),2)*pow(F2(0),2) + 
         ((-pow(kP,2) + pow(kq,2))*pow(qe2,2)*pow(F2(0),2) + 
            4*Mp2*(kK*kP*nu*F2(0)*(2*F1(0) + F2(0)) - 
               4*me2*(pow(kP,2) - pow(kq,2))*
                (3*F1(0)*F2(0) + pow(F1(0),2) + pow(F2(0),2)) + 
               qe2*(-(F1(0)*(F1(0) + F2(0))*(pow(kP,2) - pow(kq,2))) + 
                  pow(kK,2)*(3*F1(0)*F2(0) + pow(F1(0),2) + pow(F2(0),2))))\
)*pow(F2(qe2),2)) + pow(F2(qe2),2)*
       (-4*Mp2*(pow(kP,2) - pow(kq,2))*
          (F2(0)*(2*F1(0) + F2(0))*pow(kq,2)*pow(nu,2) + 
            pow(qe2,2)*(F1(0)*(F1(0) + 2*F2(0))*pow(kq,2) + 
               pow(kP,2)*(-2*F1(0)*F2(0) + pow(F1(0),2) - pow(F2(0),2)) + 
               pow(kK,2)*(4*F1(0)*F2(0) + pow(F1(0),2) + 2*pow(F2(0),2))) + 
            2*qe2*(kK*kP*nu*(3*F1(0)*F2(0) + pow(F1(0),2) + pow(F2(0),2)) - 
               2*me2*(pow(kP,2) - pow(kq,2))*
                (4*F1(0)*F2(0) + pow(F1(0),2) + 2*pow(F2(0),2)))) + 
         16*Mp4*(4*me2*F2(0)*(2*F1(0) + F2(0))*(pow(kP,2) - pow(kq,2))*
             pow(kq,2) + qe2*pow(kK,2)*
             (-(F2(0)*(2*F1(0) + F2(0))*pow(kq,2)) + 
               pow(kP,2)*pow(F1(0) + F2(0),2))) + 
         qe2*(-pow(nu,2) + pow(qe2,2))*pow(F2(0),2)*
          pow(pow(kP,2) - pow(kq,2),2))))/4.
);
}

