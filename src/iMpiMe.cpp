#include "constants.h"
#include "proton_formfactors.h"
#include <cmath>

double iMpiMe(double nu, double qe2, double qp2, double kK, double kP, double kq) {
return (2*e6*Z3*pow(kK - kq,-1)*pow(kK + kq,-1)*pow(-kP + kq,-1)*pow(kP + kq,-1)*
  pow(Mp,-2)*pow(qe2,-1)*(kK*kP*
     (F1(qp2)*F2(0)*F2(qe2)*pow(kq,2)*(pow(kK,2) - pow(kP,2) + pow(kq,2)) + 
       F2(qp2)*(F2(qe2)*((2*F1(0) + 3*F2(0))*pow(kq,2)*
              (-pow(kP,2) + pow(kq,2)) + 
             pow(kK,2)*((F1(0) + F2(0))*pow(kP,2) + F2(0)*pow(kq,2))) + 
          F1(qe2)*F2(0)*(pow(kK,2)*pow(kP,2) - 2*pow(kP,2)*pow(kq,2) + 
             2*pow(kq,4)))) - kq*
     (F2(qe2)*(F2(qp2)*(4*kK*kP*
              (-2*me2*(F1(0) + 2*F2(0))*(pow(kP,2) - pow(kq,2)) + 
                Mp2*F1(0)*pow(kq,2)) + 
             kK*kP*qe2*(2*(F1(0) + 2*F2(0))*pow(kK,2) - 
                (F1(0) + 3*F2(0))*pow(kP,2) + 3*(F1(0) + F2(0))*pow(kq,2)\
) + nu*((2*F1(0) + 3*F2(0))*pow(kq,2)*(-pow(kP,2) + pow(kq,2)) + 
                pow(kK,2)*((F1(0) + 3*F2(0))*pow(kP,2) + F2(0)*pow(kq,2)))\
) + F1(qp2)*(4*kK*kP*Mp2*(2*F1(0) + F2(0))*pow(kq,2) + 
             F2(0)*(8*kK*kP*me2*(-pow(kP,2) + pow(kq,2)) + 
                2*kK*kP*qe2*(pow(kK,2) - pow(kP,2) + pow(kq,2)) + 
                nu*(-(pow(kP,2)*pow(kq,2)) + 
                   pow(kK,2)*(2*pow(kP,2) + pow(kq,2)) + pow(kq,4))))) + 
       F1(qe2)*(4*kK*kP*Mp2*F1(qp2)*
           ((F1(0) + F2(0))*pow(kK,2) + F1(0)*pow(kP,2) + 
             (F1(0) + F2(0))*pow(kq,2)) + 
          F2(qp2)*(4*kK*kP*Mp2*
              ((3*F1(0) + F2(0))*pow(kK,2) + F1(0)*pow(kq,2)) + 
             F2(0)*(8*kK*kP*me2*(-pow(kP,2) + pow(kq,2)) + 
                kK*kP*qe2*(2*pow(kK,2) - pow(kP,2) + pow(kq,2)) + 
                nu*(pow(kK,2)*pow(kP,2) - 2*pow(kP,2)*pow(kq,2) + 
                   2*pow(kq,4)))))) + 
    F1(0)*(kK*kP*qe2 + nu*pow(kq,2))*
     (4*Mp2*F1(qe2)*(2*qe2*(2*me2 + qe2)*F2(qp2) + 
          F1(qp2)*(qe2*(4*me2 + 4*Mp2 + qe2) + pow(nu,2))) + 
       qe2*F2(qe2)*(8*Mp2*(2*me2 + qe2)*F1(qp2) + 
          F2(qp2)*(16*me2*Mp2 + 4*Mp2*qe2 - pow(nu,2) + pow(qe2,2)))) + 
    F1(qe2)*(4*Mp2*F1(qp2)*(4*kP*Mp2*F1(0)*pow(kK,3) + 
          4*kK*me2*F1(0)*pow(kP,3) + nu*F1(0)*pow(kK,2)*pow(kq,2) + 
          nu*F2(0)*pow(kK,2)*pow(kq,2) + nu*F1(0)*pow(kP,2)*pow(kq,2) + 
          kK*kP*qe2*((F1(0) + F2(0))*pow(kK,2) + F1(0)*pow(kP,2) + 
             (3*F1(0) + F2(0))*pow(kq,2)) + nu*F1(0)*pow(kq,4) + 
          nu*F2(0)*pow(kq,4)) + 
       F2(qp2)*(4*Mp2*(kK*kP*qe2*
              (2*F1(0)*pow(kK,2) + (7*F1(0) + F2(0))*pow(kq,2)) + 
             pow(kq,2)*(-4*kK*kP*me2*F2(0) + 
                nu*((3*F1(0) + F2(0))*pow(kK,2) + F1(0)*pow(kq,2)))) + 
          F2(0)*(qe2*(pow(kP,2) - pow(kq,2))*
              (-4*kK*kP*me2 + 2*nu*pow(kK,2) - 3*nu*pow(kq,2)) - 
             nu*pow(kq,2)*(kK*kP*nu + 8*me2*(-pow(kP,2) + pow(kq,2))) + 
             kP*pow(kK,3)*pow(qe2,2)))) + 
    F2(qe2)*(F2(qp2)*(-(pow(kq,2)*
             (nu*(kK*kP*nu*(F1(0) - 2*F2(0)) - 4*Mp2*F1(0)*pow(kq,2)) + 
               8*me2*(2*kK*kP*Mp2*(F1(0) + 2*F2(0)) + 
                  nu*F1(0)*(-pow(kP,2) + pow(kq,2))))) + 
          qe2*(4*kK*kP*(-(me2*(F1(0) + 2*F2(0))*
                   (pow(kP,2) - pow(kq,2))) + 
                Mp2*(-(F2(0)*pow(kK,2)) + (3*F1(0) + F2(0))*pow(kq,2))) + 
             nu*(pow(kK,2)*((F1(0) + 3*F2(0))*pow(kP,2) + 
                   F2(0)*pow(kq,2)) + 
                pow(kq,2)*(-((F1(0) + 2*F2(0))*pow(kP,2)) + 
                   (3*F1(0) + 2*F2(0))*pow(kq,2)))) + 
          kK*kP*((F1(0) + 2*F2(0))*pow(kK,2) - F2(0)*pow(kP,2) + 
             (5*F1(0) + F2(0))*pow(kq,2))*pow(qe2,2)) + 
       F1(qp2)*(4*Mp2*(-4*kK*kP*me2*F2(0)*pow(kq,2) + 
             kK*kP*qe2*(F1(0)*pow(kK,2) + (4*F1(0) + F2(0))*pow(kq,2)) + 
             nu*(2*F1(0) + F2(0))*pow(kq,4)) + 
          F2(0)*(nu*pow(kq,2)*
              (3*kK*kP*nu + 8*me2*(-pow(kP,2) + pow(kq,2))) + 
             qe2*(4*kK*kP*me2*(-pow(kP,2) + pow(kq,2)) + 
                nu*((pow(kP,2) - pow(kq,2))*pow(kq,2) + 
                   pow(kK,2)*(pow(kP,2) + 3*pow(kq,2)))) + 
             kK*kP*(pow(kK,2) - pow(kP,2) + pow(kq,2))*pow(qe2,2)))) - 
    kq*(F1(qe2)*(4*Mp2*F1(0)*F1(qp2)*
           (kK*kP*qe2*(8*me2 + 3*qe2) + 4*kK*Mp2*(kK*nu + 2*kP*qe2) + 
             nu*(-4*me2*(pow(kP,2) - 2*pow(kq,2)) + 2*qe2*pow(kq,2)) + 
             kK*kP*pow(nu,2)) + 
          F2(qp2)*(4*Mp2*(kK*kP*qe2*(8*me2 + 7*qe2)*F1(0) + 
                nu*(4*me2*(2*F1(0) + F2(0))*pow(kq,2) + 
                   qe2*((F1(0) + F2(0))*pow(kK,2) + 
                      (4*F1(0) - F2(0))*pow(kq,2)))) - 
             nu*F2(0)*(2*qe2*(kK*kP*nu + 2*me2*(-pow(kP,2) + pow(kq,2))) + 
                pow(kq,2)*pow(nu,2) + 
                (pow(kK,2) + pow(kP,2) - pow(kq,2))*pow(qe2,2)))) + 
       F2(qe2)*(F1(qp2)*(-4*Mp2*
              (-(kK*kP*qe2*(8*me2 + 5*qe2)*F1(0)) + 
                nu*(4*me2*F2(0)*pow(kq,2) + 
                   qe2*((F1(0) + F2(0))*pow(kK,2) - 
                      (4*F1(0) + F2(0))*pow(kq,2)))) + 
             nu*F2(0)*(2*qe2*(kK*kP*nu + 2*me2*(-pow(kP,2) + pow(kq,2))) + 
                pow(kq,2)*pow(nu,2) + 
                (pow(kK,2) + pow(kP,2) - pow(kq,2))*pow(qe2,2))) + 
          F1(0)*F2(qp2)*(nu*qe2*
              (4*me2*(pow(kP,2) - pow(kq,2)) + 3*qe2*pow(kq,2)) + 
             4*Mp2*(kK*kP*qe2*(8*me2 + 3*qe2) + 
                nu*(-(qe2*pow(kK,2)) + 4*me2*pow(kq,2) + 3*qe2*pow(kq,2))) - 
             2*kK*kP*qe2*pow(nu,2) - pow(kq,2)*pow(nu,3) + 4*kK*kP*pow(qe2,3)\
))))*pow(qp2,-1)
);
}

