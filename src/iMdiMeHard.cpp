#include "constants.h"
#include "proton_formfactors.h"
#include "delta_formfactors.h"
#include <cmath>

double iMdiMeHard1(double nu, double qe2, double qp2, double kK, double kP, double kq) {
return ((e6*(kP + kq - Md2 + Mp2)*Z3*pow(kK - kq,-1)*pow(kK + kq,-1)*pow(Md,-5)*
    pow(Mp,-1)*pow(qe2,-1)*(2*Mp*F1(qp2)*
       (kK*(kP + kq)*qe2*G3(qe2)*
          (Mp*((-12*kP*G1(0) - 24*kq*G1(0) + 15*kP*G2(0) + 26*kq*G2(0))*
                pow(kK,2) + kq*
                (4*kP*kq*G1(0) - 9*kP*kq*G2(0) + G2(0)*pow(kP,2) + 
                  24*G1(0)*pow(kq,2) - 26*G2(0)*pow(kq,2))) + 
            Md*((6*kP*G1(0) + 16*kq*G1(0) - 9*kP*G2(0) - 2*kq*G2(0))*
                pow(kK,2) + kq*
                (kP*kq*(-2*G1(0) + 7*G2(0)) + 
                  (2*G1(0) - G2(0))*pow(kP,2) + 
                  2*(-8*G1(0) + G2(0))*pow(kq,2)))) + 
         Md*G1(qe2)*(2*kK*kP*Md*Mp*
             (2*kq*(kP + kq)*(-4*kq*G1(0) + kP*G2(0)) + 
               (8*kP*G1(0) + 8*kq*G1(0) - 5*kP*G2(0) - 8*kq*G2(0))*
                pow(kK,2)) - 2*kK*kP*Mp2*
             ((4*kP*G1(0) + 10*kq*G1(0) - 9*kP*G2(0) - 6*kq*G2(0))*
                pow(kK,2) + 2*kq*
                (kP*kq*(3*G1(0) + 2*G2(0)) + 
                  (-2*G1(0) + G2(0))*pow(kP,2) + 
                  (5*G1(0) - 2*G2(0))*pow(kq,2))) + 
            (kP + kq)*(-4*kK*me2*
                (kP*(2*G1(0) - G2(0)) - kq*(10*G1(0) + G2(0)))*
                (pow(kP,2) - pow(kq,2)) + 
               nu*((kP - kq)*
                   (kP*(2*G1(0) - G2(0)) - kq*(10*G1(0) + G2(0)))*
                   pow(kK,2) + 
                  pow(kq,2)*(2*kP*kq*(4*G1(0) + G2(0)) + 
                     (-2*G1(0) + G2(0))*pow(kP,2) + 
                     (10*G1(0) - 11*G2(0))*pow(kq,2))) + 
               kK*qe2*((2*kP*G1(0) - 22*kq*G1(0) + 5*kP*G2(0) + 
                     5*kq*G2(0))*pow(kK,2) + 
                  kq*(4*G1(0) - 2*G2(0))*pow(kP,2) + 
                  (-2*G1(0) + G2(0))*pow(kP,3) + 
                  4*kP*(7*G1(0) - 5*G2(0))*pow(kq,2) + 
                  (22*G1(0) - 5*G2(0))*pow(kq,3)))) + 
         Md*G2(qe2)*(-2*kK*kP*Md*Mp*
             (kq*(kP + kq)*(-2*kP*G1(0) + 2*kP*G2(0) + kq*G2(0)) + 
               (5*kP*G1(0) + 8*kq*G1(0) - 6*kP*G2(0) - 9*kq*G2(0))*
                pow(kK,2)) + 2*kK*kP*Mp2*
             (-(kq*(kP + kq)*(2*kP*G1(0) + 2*kq*G1(0) - 2*kP*G2(0) - 
                    kq*G2(0))) + 
               (9*kP*G1(0) + 12*kq*G1(0) - 10*kP*G2(0) - 7*kq*G2(0))*
                pow(kK,2)) + (kP + kq)*
             (nu*((kP + kq)*(kP*(G1(0) - 2*G2(0)) + kq*(G1(0) + G2(0)))*
                   pow(kq,2) + 
                  pow(kK,2)*(3*kP*kq*G2(0) - 
                     (G1(0) - 2*G2(0))*pow(kP,2) + 
                     (-11*G1(0) + G2(0))*pow(kq,2))) + 
               kK*qe2*(-3*(-3*kq*G1(0) + kP*(G1(0) + G2(0)))*pow(kK,2) + 
                  kq*(-2*G1(0) + 7*G2(0))*pow(kP,2) + 
                  (G1(0) - 2*G2(0))*pow(kP,3) + 
                  6*kP*(-2*G1(0) + G2(0))*pow(kq,2) - 9*G1(0)*pow(kq,3)) + 
               4*kK*(kP - kq)*me2*(G1(0) - 2*G2(0))*pow(kP + kq,2)))) + 
      kK*(2*kP*(kP + kq)*Md*Mp*F1(qp2)*
          (G2(qe2)*(-(kq*(kP + kq)*
                  (kP*(G1(0) - 2*G2(0)) + kq*(G1(0) + G2(0)))) + 
               (5*kP*G1(0) + 11*kq*G1(0) - 7*kP*G2(0) - kq*G2(0))*
                pow(kK,2)) + G1(qe2)*
             ((-(kq*(10*G1(0) + G2(0))) + kP*(2*G1(0) + 5*G2(0)))*
                pow(kK,2) - kq*
                (2*kP*kq*(4*G1(0) + G2(0)) + 
                  (-2*G1(0) + G2(0))*pow(kP,2) + 
                  (10*G1(0) - 11*G2(0))*pow(kq,2)))) + 
         F2(qp2)*(2*kP*Md*G1(qe2)*
             ((kP + kq)*Mp*((2*kP*G1(0) - 8*kq*G1(0) + kP*G2(0) - 
                     6*kq*G2(0))*pow(kK,2) + 
                  2*(-4*kP*G1(0) - 4*kq*G1(0) + kP*G2(0) + 7*kq*G2(0))*
                   pow(kq,2)) + 
               Md*(-((kP + kq)*
                     (kq*(-8*G1(0) + G2(0)) + kP*(4*G1(0) + G2(0)))*
                     pow(kq,2)) + 
                  pow(kK,2)*(2*kP*kq*(-2*G1(0) + G2(0)) + 
                     (4*G1(0) - 2*G2(0))*pow(kP,2) + 
                     (-8*G1(0) + 7*G2(0))*pow(kq,2)))) - 
            (kP + kq)*qe2*G3(qe2)*
             (pow(kK,2)*(-4*kP*kq*(G1(0) - 2*G2(0)) + 
                  (6*G1(0) - 3*G2(0))*pow(kP,2) + 
                  2*(-10*G1(0) + 11*G2(0))*pow(kq,2)) + 
               kq*(3*kq*(-2*G1(0) + G2(0))*pow(kP,2) + 
                  (2*G1(0) - G2(0))*pow(kP,3) - 6*kP*G2(0)*pow(kq,2) + 
                  2*(10*G1(0) - 11*G2(0))*pow(kq,3))) - 
            2*kP*Md*G2(qe2)*(-((kP + kq)*Mp*
                  ((kP*G1(0) + 6*kq*G1(0) - 3*kP*G2(0) + 5*kq*G2(0))*
                     pow(kK,2) + (kP + kq)*(2*G1(0) - 3*G2(0))*pow(kq,2))) \
+ Md*(pow(kK,2)*(-2*kP*kq*(G1(0) - 2*G2(0)) + 
                     (2*G1(0) - G2(0))*pow(kP,2) + 
                     (-7*G1(0) + 8*G2(0))*pow(kq,2)) + 
                  (G1(0) - 2*G2(0))*pow(kq,2)*pow(kP + kq,2))))) + 
      F2(qp2)*(Md*G1(qe2)*(-(Mp*
               (4*kK*kP*Mp2*((2*kP*G1(0) + 4*kq*G1(0) - 2*kP*G2(0) + 
                       3*kq*G2(0))*pow(kK,2) + 
                    (4*kP*G1(0) + 4*kq*G1(0) - kP*G2(0) - 7*kq*G2(0))*
                     pow(kq,2)) - 
                 kK*qe2*(pow(kK,2)*
                     (13*kP*kq*(-2*G1(0) + G2(0)) + 
                       4*(2*G1(0) + 3*G2(0))*pow(kP,2) - 
                       2*(14*G1(0) + G2(0))*pow(kq,2)) + 
                    kq*(kP + kq)*
                     (kP*kq*(34*G1(0) - 33*G2(0)) + 
                       7*(2*G1(0) - G2(0))*pow(kP,2) + 
                       2*(14*G1(0) + G2(0))*pow(kq,2))) - 
                 4*kq*(kP + kq)*
                  (2*kK*(kP + kq)*me2*
                     (kP*(10*G1(0) - 3*G2(0)) - kq*(10*G1(0) + G2(0))) + 
                    nu*((-6*kP*G1(0) + 4*kq*G1(0) + 2*kP*G2(0) + 
                        3*kq*G2(0))*pow(kK,2) + 
                       (4*kP*G1(0) + 4*kq*G1(0) - kP*G2(0) - 7*kq*G2(0))*
                        pow(kq,2))))) + 
            Md*(4*kK*kP*Mp2*((2*kP*G1(0) - 4*kq*G1(0) - 2*kP*G2(0) + 
                     3*kq*G2(0))*pow(kK,2) + 
                  (-2*kP*G1(0) + 4*kq*G1(0) - kP*G2(0) - kq*G2(0))*
                   pow(kq,2)) + 
               kK*qe2*(pow(kK,2)*
                   (-3*kP*kq*G2(0) + (8*G1(0) - 4*G2(0))*pow(kP,2) - 
                     8*(G1(0) - 2*G2(0))*pow(kq,2)) + 
                  kq*(kP + kq)*
                   (-(kP*kq*(8*G1(0) + 5*G2(0))) + 9*G2(0)*pow(kP,2) + 
                     8*(G1(0) - 2*G2(0))*pow(kq,2))) + 
               2*kq*(4*kK*kP*(kP + kq)*me2*
                   (kq*(-8*G1(0) + G2(0)) + kP*(4*G1(0) + G2(0))) + 
                  nu*((kP + kq)*
                      (kq*(-8*G1(0) + G2(0)) + kP*(4*G1(0) + G2(0)))*
                      pow(kq,2) - 
                     pow(kK,2)*
                      (2*kP*kq*(-2*G1(0) + G2(0)) + 
                        (4*G1(0) + G2(0))*pow(kP,2) + 
                        (-8*G1(0) + 7*G2(0))*pow(kq,2)))))) - 
         qe2*G3(qe2)*(-2*kK*Md*Mp*
             (kq*(kP + kq)*(kP*kq*(4*G1(0) + 3*G2(0)) + 
                  (G1(0) - G2(0))*pow(kP,2) - 
                  2*(5*G1(0) + 3*G2(0))*pow(kq,2)) + 
               pow(kK,2)*(kP*kq*(12*G1(0) + G2(0)) + 
                  (G1(0) - 5*G2(0))*pow(kP,2) + 
                  2*(5*G1(0) + 3*G2(0))*pow(kq,2))) + 
            2*kK*Mp2*(pow(kK,2)*
                (kP*kq*(-6*G1(0) + 9*G2(0)) + 
                  (11*G1(0) - 9*G2(0))*pow(kP,2) + 
                  2*(-11*G1(0) + 12*G2(0))*pow(kq,2)) + 
               kq*(kq*(-7*G1(0) + 4*G2(0))*pow(kP,2) + 
                  (G1(0) - G2(0))*pow(kP,3) + 
                  kP*(8*G1(0) - 13*G2(0))*pow(kq,2) + 
                  2*(11*G1(0) - 12*G2(0))*pow(kq,3))) + 
            (kP + kq)*(-16*kK*kP*kq*(kP + kq)*me2*(G1(0) - 2*G2(0)) + 
               nu*(pow(kK,2)*
                   (kP*kq*(14*G1(0) - 19*G2(0)) + 
                     (-2*G1(0) + G2(0))*pow(kP,2) + 
                     2*(4*G1(0) - 5*G2(0))*pow(kq,2)) + 
                  pow(kq,2)*(kP*kq*(-14*G1(0) + 19*G2(0)) + 
                     (2*G1(0) - G2(0))*pow(kP,2) - 
                     4*(G1(0) - 2*G2(0))*pow(kq,2))) + 
               kK*qe2*((2*kP*G1(0) + 14*kq*G1(0) - 4*kP*G2(0) - 
                     19*kq*G2(0))*pow(kK,2) + 
                  (-2*G1(0) + G2(0))*pow(kP,3) + 
                  kP*(4*G1(0) + G2(0))*pow(kq,2) + 
                  (-14*G1(0) + 19*G2(0))*pow(kq,3)))) - 
         Md*G2(qe2)*(Mp*(-4*kK*kP*Mp2*
                ((2*kP*(G1(0) - G2(0)) + 3*kq*(G1(0) + G2(0)))*
                   pow(kK,2) + (kP + kq)*(G1(0) - G2(0))*pow(kq,2)) + 
               2*kq*(kP + kq)*
                (4*kK*(kP + kq)*me2*
                   (3*kP*G1(0) + kq*G1(0) - 2*kP*G2(0) - 2*kq*G2(0)) + 
                  nu*((-4*kP*G1(0) + 6*kq*G1(0) + kP*G2(0) + 
                        5*kq*G2(0))*pow(kK,2) + 
                     (kP + kq)*(2*G1(0) - 3*G2(0))*pow(kq,2))) + 
               kK*qe2*(pow(kK,2)*
                   (kP*kq*(-13*G1(0) + 4*G2(0)) + 
                     2*(2*G1(0) + 5*G2(0))*pow(kP,2) - 
                     2*(7*G1(0) + 4*G2(0))*pow(kq,2)) + 
                  kq*(kP + kq)*
                   (kP*kq*(17*G1(0) - 2*G2(0)) + 
                     (7*G1(0) - 16*G2(0))*pow(kP,2) + 
                     2*(7*G1(0) + 4*G2(0))*pow(kq,2)))) + 
            Md*(4*kK*kP*Mp2*(G1(0) - G2(0))*
                ((2*kP - 3*kq)*pow(kK,2) + (kP + kq)*pow(kq,2)) + 
               kK*qe2*(kq*(kP + kq)*
                   (kP*kq*(G1(0) - 2*G2(0)) + 
                     (-7*G1(0) + 8*G2(0))*pow(kP,2) + 
                     12*(G1(0) - G2(0))*pow(kq,2)) + 
                  pow(kK,2)*(kP*kq*(7*G1(0) - 8*G2(0)) + 
                     6*(G1(0) - G2(0))*pow(kP,2) + 
                     12*(-G1(0) + G2(0))*pow(kq,2))) - 
               2*kq*(4*kK*kP*me2*(G1(0) - 2*G2(0))*pow(kP + kq,2) + 
                  nu*(pow(kK,2)*
                      (-2*kP*kq*(G1(0) - 2*G2(0)) - 
                        (G1(0) - 2*G2(0))*pow(kP,2) + 
                        (-7*G1(0) + 8*G2(0))*pow(kq,2)) + 
                     (G1(0) - 2*G2(0))*pow(kq,2)*pow(kP + kq,2)))))) + 
      2*Mp*F1(qp2)*(qe2*G3(qe2)*
          (Md*(4*kK*Mp2*((7*kP*G1(0) + 10*kq*G1(0) - 8*kP*G2(0) - 
                     5*kq*G2(0))*pow(kK,2) + 
                  (-7*kP*(G1(0) - G2(0)) + 5*kq*(-2*G1(0) + G2(0)))*
                   pow(kq,2)) - 
               (kP + kq)*(kK*qe2*(2*G1(0) - G2(0))*
                   (3*kP*kq + 5*pow(kK,2) + pow(kP,2) - 5*pow(kq,2)) - 
                  8*kK*me2*(G1(0) - 2*G2(0))*(pow(kP,2) - pow(kq,2)) + 
                  nu*((6*kP*G1(0) - 2*kq*G1(0) - 9*kP*G2(0) + 
                        7*kq*G2(0))*pow(kK,2) - 
                     3*(kP - kq)*(2*G1(0) - 3*G2(0))*pow(kq,2)))) + 
            Mp*(4*kK*Mp2*((-11*kP*G1(0) - 14*kq*G1(0) + 12*kP*G2(0) + 
                     15*kq*G2(0))*pow(kK,2) + 
                  (9*kP*G1(0) + 14*kq*G1(0) - 11*kP*G2(0) - 15*kq*G2(0))*
                   pow(kq,2)) + 
               (kP + kq)*(-16*kK*(kP + kq)*me2*
                   (-(kq*G1(0)) + 2*kP*(G1(0) - G2(0))) + 
                  kK*qe2*(kP*kq*(8*G1(0) - 7*G2(0)) + 
                     (4*G1(0) - 7*G2(0))*pow(kK,2) - G2(0)*pow(kP,2) + 
                     (-4*G1(0) + 7*G2(0))*pow(kq,2)) + 
                  nu*((16*kP*G1(0) + 4*kq*G1(0) - 17*kP*G2(0) - 
                        7*kq*G2(0))*pow(kK,2) + 
                     (-16*kP*G1(0) + 4*kq*G1(0) + 17*kP*G2(0) + kq*G2(0))*
                      pow(kq,2))))) + 
         Md*G1(qe2)*(2*Md*Mp*(8*kK*kP*Mp2*
                ((G1(0) - G2(0))*pow(kK,2) - G1(0)*pow(kq,2)) + 
               kK*qe2*((4*kP*G1(0) + 4*kq*G1(0) - 3*kP*G2(0) - 
                     9*kq*G2(0))*pow(kK,2) - 
                  (kP + kq)*(-6*kP*kq*G2(0) + 2*G2(0)*pow(kP,2) + 
                     (4*G1(0) - 9*G2(0))*pow(kq,2))) + 
               2*(-2*kK*me2*(kP*G2(0) + kq*(-4*G1(0) + G2(0)))*
                   (kP*kq + 2*pow(kP,2) - pow(kq,2)) + 
                  nu*((kP + kq)*(4*kq*G1(0) - kP*G2(0))*pow(kq,2) + 
                     pow(kK,2)*
                      (2*kP*kq*(-2*G1(0) + G2(0)) + G2(0)*pow(kP,2) + 
                        4*(-G1(0) + G2(0))*pow(kq,2))))) + 
            2*Mp2*(2*(-2*kK*(kP + kq)*me2*
                   (-(kP*kq*(8*G1(0) + 3*G2(0))) + 
                     (4*G1(0) - 2*G2(0))*pow(kP,2) + 
                     (12*G1(0) - G2(0))*pow(kq,2)) + 
                  nu*(pow(kK,2)*
                      (kP*kq*(G1(0) - 4*G2(0)) + 
                        (2*G1(0) - G2(0))*pow(kP,2) + 
                        (5*G1(0) - 3*G2(0))*pow(kq,2)) + 
                     pow(kq,2)*
                      (kP*kq*(3*G1(0) + 2*G2(0)) + 
                        (-2*G1(0) + G2(0))*pow(kP,2) + 
                        (5*G1(0) - 2*G2(0))*pow(kq,2)))) + 
               kK*qe2*((-16*kP*G1(0) - 28*kq*G1(0) + 13*kP*G2(0) + 
                     13*kq*G2(0))*pow(kK,2) + 
                  6*kq*(2*G1(0) - G2(0))*pow(kP,2) + 
                  (-4*G1(0) + 2*G2(0))*pow(kP,3) + 
                  kP*(44*G1(0) - 27*G2(0))*pow(kq,2) + 
                  (28*G1(0) - 13*G2(0))*pow(kq,3))) - 
            16*kK*kP*((G1(0) - G2(0))*pow(kK,2) + G1(0)*pow(kq,2))*
             pow(Mp,4) - (kP + kq)*(2*G1(0) - G2(0))*
             (qe2*(nu*(kP*kq*(kP + kq) + (-7*kP + 5*kq)*pow(kK,2)) + 
                  12*kK*me2*(pow(kP,2) - pow(kq,2))) + 
               kK*(kP - 5*kq)*kq*pow(nu,2) - 
               kK*(-5*kP*kq + 6*pow(kK,2) + pow(kP,2) - 6*pow(kq,2))*
                pow(qe2,2))) + 
         Md*G2(qe2)*(-2*Md*Mp*
             (-4*kP*kq*nu*G1(0)*pow(kK,2) + 5*kP*kq*nu*G2(0)*pow(kK,2) + 
               8*kP*Mp2*(G1(0) - G2(0))*pow(kK,3) + 
               12*kK*kq*me2*G1(0)*pow(kP,2) - 
               16*kK*kq*me2*G2(0)*pow(kP,2) - 
               2*nu*G1(0)*pow(kK,2)*pow(kP,2) + 
               2*nu*G2(0)*pow(kK,2)*pow(kP,2) + 
               8*kK*me2*G1(0)*pow(kP,3) - 8*kK*me2*G2(0)*pow(kP,3) - 
               8*kK*kP*me2*G2(0)*pow(kq,2) - 
               8*nu*G1(0)*pow(kK,2)*pow(kq,2) + 
               9*nu*G2(0)*pow(kK,2)*pow(kq,2) + 
               2*nu*G1(0)*pow(kP,2)*pow(kq,2) - 
               2*nu*G2(0)*pow(kP,2)*pow(kq,2) - 
               4*kK*me2*G1(0)*pow(kq,3) + 2*kP*nu*G1(0)*pow(kq,3) - 
               3*kP*nu*G2(0)*pow(kq,3) + 
               kK*qe2*((3*kP*G1(0) + 7*kq*G1(0) - 3*kP*G2(0) - 
                     8*kq*G2(0))*pow(kK,2) + 
                  4*kq*(-G1(0) + G2(0))*pow(kP,2) + 
                  2*(G1(0) - G2(0))*pow(kP,3) + 
                  kP*(-11*G1(0) + 13*G2(0))*pow(kq,2) + 
                  (-7*G1(0) + 8*G2(0))*pow(kq,3)) - nu*G2(0)*pow(kq,4)) + 
            2*Mp2*(nu*((kP + kq)*
                   (2*kP*G1(0) + 2*kq*G1(0) - 2*kP*G2(0) - kq*G2(0))*
                   pow(kq,2) + 
                  pow(kK,2)*(kP*kq*(-8*G1(0) + 9*G2(0)) - 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     (-12*G1(0) + 7*G2(0))*pow(kq,2))) + 
               kK*qe2*((5*kP*G1(0) + 11*kq*G1(0) - 7*kP*G2(0) - 
                     6*kq*G2(0))*pow(kK,2) + 
                  kq*(-6*G1(0) + 8*G2(0))*pow(kP,2) + 
                  2*(G1(0) - G2(0))*pow(kP,3) + 
                  kP*(-19*G1(0) + 15*G2(0))*pow(kq,2) + 
                  (-11*G1(0) + 6*G2(0))*pow(kq,3)) + 
               4*kK*me2*(kq*G1(0) + 2*kP*(G1(0) - G2(0)))*pow(kP + kq,2)) \
+ 16*kP*(G1(0) - G2(0))*pow(kK,3)*pow(Mp,4) - 
            (kP + kq)*(qe2*(-4*kK*me2*(G1(0) - 2*G2(0))*
                   (pow(kP,2) - pow(kq,2)) + 
                  nu*(3*(-3*kq*G1(0) + kP*(G1(0) - G2(0)) + 2*kq*G2(0))*
                      pow(kK,2) + 
                     kq*(3*kP*kq*G1(0) - (G1(0) - 2*G2(0))*pow(kP,2) + 
                        (4*G1(0) - 5*G2(0))*pow(kq,2)))) + 
               kK*kq*(-(kP*G1(0)) + 5*kq*G1(0) + 2*kP*G2(0) - 4*kq*G2(0))*
                pow(nu,2) + kK*
                (kP*kq*(-5*G1(0) + G2(0)) + 
                  (6*G1(0) - 3*G2(0))*pow(kK,2) + 
                  (G1(0) + G2(0))*pow(kP,2) + 
                  3*(-2*G1(0) + G2(0))*pow(kq,2))*pow(qe2,2)))) + 
      F2(qp2)*(qe2*G3(qe2)*(-2*Mp2*
             (-8*kK*kq*(kP + kq)*me2*
                (2*kq*G1(0) + kP*(G1(0) - 2*G2(0))) + 
               nu*(pow(kK,2)*
                   (kP*kq*(11*G1(0) - 12*G2(0)) + 
                     (-G1(0) + G2(0))*pow(kP,2) + 
                     (8*G1(0) - 7*G2(0))*pow(kq,2)) + 
                  pow(kq,2)*(kP*kq*(-15*G1(0) + 16*G2(0)) + 
                     (G1(0) - G2(0))*pow(kP,2) + 
                     (-10*G1(0) + 11*G2(0))*pow(kq,2))) + 
               kK*qe2*((9*kP*G1(0) + 15*kq*G1(0) - 11*kP*G2(0) - 
                     19*kq*G2(0))*pow(kK,2) + 
                  kq*(-4*G1(0) + 6*G2(0))*pow(kP,2) + 
                  (-G1(0) + G2(0))*pow(kP,3) - 
                  2*kP*(5*G1(0) - 8*G2(0))*pow(kq,2) + 
                  (-15*G1(0) + 19*G2(0))*pow(kq,3))) + 
            2*Md*Mp*(7*kP*kq*nu*G1(0)*pow(kK,2) - 
               10*kP*kq*nu*G2(0)*pow(kK,2) - 
               24*kK*kq*me2*G1(0)*pow(kP,2) + 
               16*kK*kq*me2*G2(0)*pow(kP,2) - 
               nu*G1(0)*pow(kK,2)*pow(kP,2) + 
               nu*G2(0)*pow(kK,2)*pow(kP,2) - 
               32*kK*kP*me2*G1(0)*pow(kq,2) + 
               32*kK*kP*me2*G2(0)*pow(kq,2) + 
               6*nu*G1(0)*pow(kK,2)*pow(kq,2) - 
               11*nu*G2(0)*pow(kK,2)*pow(kq,2) + 
               nu*G1(0)*pow(kP,2)*pow(kq,2) - 
               nu*G2(0)*pow(kP,2)*pow(kq,2) + 
               4*kK*Mp2*((2*kP*(G1(0) - G2(0)) + 3*kq*(G1(0) + G2(0)))*
                   pow(kK,2) - 
                  (kP*(G1(0) - G2(0)) + 3*kq*(G1(0) + G2(0)))*pow(kq,2)) \
- 8*kK*me2*G1(0)*pow(kq,3) - 11*kP*nu*G1(0)*pow(kq,3) + 
               16*kK*me2*G2(0)*pow(kq,3) + 12*kP*nu*G2(0)*pow(kq,3) - 
               kK*qe2*((kP*(5*G1(0) + G2(0)) + kq*(G1(0) + 3*G2(0)))*
                   pow(kK,2) + kq*(8*G1(0) - 4*G2(0))*pow(kP,2) + 
                  (G1(0) - G2(0))*pow(kP,3) + 
                  4*kP*(G1(0) - G2(0))*pow(kq,2) - 
                  (G1(0) + 3*G2(0))*pow(kq,3)) - 12*nu*G1(0)*pow(kq,4) + 
               13*nu*G2(0)*pow(kq,4)) - 
            8*kK*(G1(0) - G2(0))*
             ((2*kP - 3*kq)*pow(kK,2) + (-kP + 3*kq)*pow(kq,2))*pow(Mp,4) \
+ (kP + kq)*(qe2*(-8*kK*kP*(kP + kq)*me2*(G1(0) - 2*G2(0)) + 
                  nu*((4*kP*G1(0) + 4*kq*G1(0) - 8*kP*G2(0) - 
                        5*kq*G2(0))*pow(kK,2) + 
                     kq*(-4*kP*kq*(G1(0) - 2*G2(0)) + 
                        (2*G1(0) - G2(0))*pow(kP,2) - 
                        2*(G1(0) - 2*G2(0))*pow(kq,2)))) + 
               kK*kP*kq*(2*G1(0) - G2(0))*pow(nu,2) + 
               kK*(kq*(2*kP*G1(0) - 2*kq*G1(0) - kP*G2(0) + 
                     4*kq*G2(0)) + 2*(G1(0) - 2*G2(0))*pow(kK,2))*
                pow(qe2,2))) + 
         Md*G1(qe2)*(Mp*(qe2*
                (-4*kK*me2*(10*kP*G1(0) + 6*kq*G1(0) - 3*kP*G2(0) - 
                     3*kq*G2(0))*(-(kP*kq) + pow(kP,2) - 2*pow(kq,2)) + 
                  nu*((kP + kq)*
                      (6*kP*G1(0) + 14*kq*G1(0) - 7*kP*G2(0) - 
                        3*kq*G2(0))*pow(kq,2) + 
                     pow(kK,2)*
                      (6*kP*kq*(-2*G1(0) + G2(0)) + 
                        (18*G1(0) - 5*G2(0))*pow(kP,2) + 
                        (-18*G1(0) + 5*G2(0))*pow(kq,2)))) - 
               4*Mp2*(kK*qe2*
                   ((10*kP*G1(0) + 10*kq*G1(0) - 9*kP*G2(0) - 
                        kq*G2(0))*pow(kK,2) + 
                     kq*(-4*kP*kq*(4*G1(0) - 3*G2(0)) + 
                        (-8*G1(0) + 4*G2(0))*pow(kP,2) + 
                        (-10*G1(0) + G2(0))*pow(kq,2))) + 
                  kq*(-4*kK*(kP + kq)*me2*
                      (2*kP*G1(0) + kq*(-12*G1(0) + G2(0))) + 
                     nu*((-4*kq*G1(0) + kP*G2(0) - 3*kq*G2(0))*
                        pow(kK,2) + 
                        (-4*kP*G1(0) - 4*kq*G1(0) + kP*G2(0) + 
                        7*kq*G2(0))*pow(kq,2)))) + 
               10*kK*(kP + kq)*(2*G1(0) - G2(0))*pow(kq,2)*pow(nu,2) + 
               kK*(2*(8*kP*G1(0) + 2*kq*G1(0) - 3*kP*G2(0))*pow(kK,2) - 
                  (kP + kq)*(5*kP*kq*(2*G1(0) - G2(0)) + 
                     (2*G1(0) - G2(0))*pow(kP,2) + 4*G1(0)*pow(kq,2)))*
                pow(qe2,2)) + 
            Md*(qe2*(-4*kK*kP*(kP + kq)*me2*
                   (4*kP*G1(0) - 8*kq*G1(0) + kP*G2(0) - 5*kq*G2(0)) + 
                  nu*((kP + kq)*
                      (-8*kP*G1(0) + 8*kq*G1(0) + 7*kP*G2(0) + 
                        5*kq*G2(0))*pow(kq,2) + 
                     pow(kK,2)*
                      (-12*kP*kq*G2(0) + (8*G1(0) - G2(0))*pow(kP,2) + 
                        (-8*G1(0) + 7*G2(0))*pow(kq,2)))) + 
               4*Mp2*(kK*qe2*
                   ((2*kP*G1(0) - 2*kq*G1(0) - 3*kP*G2(0) + 
                        5*kq*G2(0))*pow(kK,2) + 
                     kq*(-2*kP*kq*G1(0) + 4*G2(0)*pow(kP,2) + 
                        (2*G1(0) - 5*G2(0))*pow(kq,2))) + 
                  kq*(4*kK*me2*
                      (kP*kq*(-8*G1(0) + G2(0)) + 2*G1(0)*pow(kP,2) + 
                        (-4*G1(0) + G2(0))*pow(kq,2)) + 
                     nu*((-2*kP*G1(0) + 4*kq*G1(0) + kP*G2(0) - 
                        3*kq*G2(0))*pow(kK,2) + 
                        (2*kP*G1(0) - 4*kq*G1(0) + kP*G2(0) + kq*G2(0))*
                         pow(kq,2)))) + 
               6*kK*kP*G2(0)*pow(kq,2)*pow(nu,2) + 
               kK*(2*(2*kP*G1(0) + 2*kq*G1(0) - kP*G2(0) - 7*kq*G2(0))*
                   pow(kK,2) - 
                  (kP + kq)*(-3*kP*kq*G2(0) + 3*G2(0)*pow(kP,2) + 
                     2*(2*G1(0) - 7*G2(0))*pow(kq,2)))*pow(qe2,2))) + 
         Md*G2(qe2)*(Md*(nu*qe2*
                ((kP + kq)*(5*kP*G1(0) + kq*G1(0) - 4*kP*G2(0) - 
                     2*kq*G2(0))*pow(kq,2) + 
                  pow(kK,2)*(6*kP*kq*(-G1(0) + G2(0)) + 
                     (G1(0) - 2*G2(0))*pow(kP,2) + 
                     (7*G1(0) - 8*G2(0))*pow(kq,2))) + 
               kK*qe2*(-4*kP*me2*(G1(0) - 2*G2(0))*
                   (pow(kP,2) - pow(kq,2)) + 
                  qe2*(-((G1(0) - 2*G2(0))*pow(kP,3)) + 
                     3*kP*(3*G1(0) - 4*G2(0))*pow(kq,2) + 
                     2*kq*(4*G1(0) - 5*G2(0))*(-pow(kK,2) + pow(kq,2)))) \
- 4*Mp2*(-(kq*(4*kK*kq*(kP + kq)*me2*G1(0) + 
                       nu*(G1(0) - G2(0))*
                        ((kP - 3*kq)*pow(kK,2) + (kP + kq)*pow(kq,2)))) \
+ kK*qe2*((-3*kq*G1(0) + 5*kP*(G1(0) - G2(0)) + 2*kq*G2(0))*pow(kK,2) + 
                     kq*(kP*kq*G2(0) - 4*(G1(0) - G2(0))*pow(kP,2) + 
                        (3*G1(0) - 2*G2(0))*pow(kq,2)))) + 
               6*kK*kP*(G1(0) - G2(0))*pow(kq,2)*pow(nu,2)) + 
            Mp*(qe2*(4*kK*(kP + kq)*me2*
                   (kP*kq*(-3*G1(0) + 2*G2(0)) + 
                     (3*G1(0) - 2*G2(0))*pow(kP,2) - 
                     2*(G1(0) - 2*G2(0))*pow(kq,2)) + 
                  nu*(-((kP + kq)*
                        (7*kP*G1(0) + 11*kq*G1(0) - 6*kP*G2(0) - 
                         12*kq*G2(0))*pow(kq,2)) + 
                     pow(kK,2)*
                      (2*kP*kq*(7*G1(0) - 5*G2(0)) - 5*G1(0)*pow(kP,2) + 
                        (13*G1(0) - 6*G2(0))*pow(kq,2)))) + 
               4*Mp2*(kq*(4*kK*kq*(kP + kq)*me2*G1(0) - 
                     nu*((kP*(G1(0) - G2(0)) + 3*kq*(G1(0) + G2(0)))*
                         pow(kK,2) + (kP + kq)*(G1(0) - G2(0))*pow(kq,2))\
) + kK*qe2*((5*kP*G1(0) + 5*kq*G1(0) - 5*kP*G2(0) + 2*kq*G2(0))*
                      pow(kK,2) + 
                     kq*(kP*kq*(-8*G1(0) + G2(0)) - 
                        4*(G1(0) - G2(0))*pow(kP,2) - 
                        (5*G1(0) + 2*G2(0))*pow(kq,2)))) - 
               2*kK*(kP + kq)*(5*G1(0) - 4*G2(0))*pow(kq,2)*pow(nu,2) + 
               kK*(-2*(5*kP*G1(0) + 2*kq*G1(0) - kP*G2(0) + kq*G2(0))*
                   pow(kK,2) + 
                  (kP + kq)*(kP*kq*(5*G1(0) + 4*G2(0)) + 
                     (G1(0) - 4*G2(0))*pow(kP,2) + 
                     2*(2*G1(0) + G2(0))*pow(kq,2)))*pow(qe2,2)))) + 
      2*Mp*F1(qp2)*(qe2*G3(qe2)*
          (Mp*(-((kP + kq)*(kq*nu - kK*qe2)*(kK*nu + kP*qe2)*G2(0)) + 
               4*Mp2*(-4*kK*(kP + kq)*me2*
                   (-(kq*G1(0)) + 2*kP*(G1(0) - G2(0))) + 
                  kK*qe2*(-(kq*
                        (-2*kP*G1(0) + kq*G1(0) + kP*G2(0) - 
                        2*kq*G2(0))) + (G1(0) - 2*G2(0))*pow(kK,2)) + 
                  nu*((4*kP + kq)*(G1(0) - G2(0))*pow(kK,2) + 
                     (-4*kP*G1(0) + kq*G1(0) + 4*kP*G2(0))*pow(kq,2))) - 
               32*kK*(G1(0) - G2(0))*(pow(kK,2) - pow(kq,2))*pow(Mp,4)) + 
            Md*(-((kP + kq)*(kq*nu - kK*qe2)*(kK*nu + kP*qe2)*
                  (2*G1(0) - G2(0))) - 
               4*Mp2*(-4*kK*(kP + kq)*me2*
                   (kq*G1(0) + 2*kP*(G1(0) - G2(0))) + 
                  kK*qe2*(-(kq*(3*kq*G1(0) + kP*G2(0) - 2*kq*G2(0))) + 
                     (3*G1(0) - 2*G2(0))*pow(kK,2)) + 
                  nu*((4*kP + kq)*(G1(0) - G2(0))*pow(kK,2) - 
                     (4*kP*G1(0) + kq*G1(0) - 4*kP*G2(0))*pow(kq,2))) + 
               32*kK*(G1(0) - G2(0))*(pow(kK,2) - pow(kq,2))*pow(Mp,4))) - 
         2*Md*Mp*G2(qe2)*(Md*(-8*kK*Mp2*(G1(0) - G2(0))*
                (kK*kq*nu + qe2*(kq*(kP + kq) - pow(kK,2))) - 
               qe2*(-4*kK*(kP + kq)*me2*
                   (-(kq*G1(0)) + 2*kP*(G1(0) - G2(0))) + 
                  nu*((3*kP*G1(0) - kq*G1(0) - 3*kP*G2(0) + 
                        2*kq*G2(0))*pow(kK,2) + 
                     kq*(3*kP*kq*(-G1(0) + G2(0)) - 
                        2*(G1(0) - G2(0))*pow(kP,2) + G1(0)*pow(kq,2)))) \
+ kK*(2*kP - kq)*kq*(G1(0) - G2(0))*pow(nu,2) - 
               kK*(2*kP*kq*G2(0) + (G1(0) - 2*G2(0))*pow(kK,2) + 
                  (-G1(0) + G2(0))*pow(kP,2) - 
                  (G1(0) - 2*G2(0))*pow(kq,2))*pow(qe2,2)) - 
            Mp*(-8*kK*Mp2*(G1(0) - G2(0))*
                (kK*kq*nu + qe2*(kq*(kP + kq) - pow(kK,2))) + 
               qe2*(4*kK*(kP + kq)*me2*
                   (kq*G1(0) + 2*kP*(G1(0) - G2(0))) + 
                  nu*((-3*kP*G1(0) + 3*kq*G1(0) + 3*kP*G2(0) - 
                        2*kq*G2(0))*pow(kK,2) + 
                     kq*(3*kP*kq*(G1(0) - G2(0)) + 
                        2*(G1(0) - G2(0))*pow(kP,2) + G1(0)*pow(kq,2)))) \
+ kK*(2*kP - kq)*kq*(G1(0) - G2(0))*pow(nu,2) + 
               kK*(kP*kq*(4*G1(0) - 2*G2(0)) + 
                  (-3*G1(0) + 2*G2(0))*pow(kK,2) + 
                  (G1(0) - G2(0))*pow(kP,2) + 
                  (3*G1(0) - 2*G2(0))*pow(kq,2))*pow(qe2,2))) + 
         2*Md*Mp*G1(qe2)*(Mp*(-4*Mp2*
                (kK*qe2*(kq*(-4*kP*G1(0) - 5*kq*G1(0) + 2*kP*G2(0) + 
                        3*kq*G2(0)) + (5*G1(0) - 3*G2(0))*pow(kK,2)) - 
                  2*kq*(2*kK*(2*kP - kq)*me2*G1(0) + 
                     nu*((G1(0) - G2(0))*pow(kK,2) + G1(0)*pow(kq,2)))) + 
               (2*G1(0) - G2(0))*
                (qe2*(-12*kK*kP*(kP + kq)*me2 + 
                     nu*((5*kP - kq)*pow(kK,2) - 
                        kq*(5*kP*kq + 2*pow(kP,2) + 3*pow(kq,2)))) + 
                  kK*kq*(-2*kP + kq)*pow(nu,2) + 
                  kK*(-4*kP*kq + 3*pow(kK,2) - pow(kP,2) - 3*pow(kq,2))*
                   pow(qe2,2))) + 
            Md*(4*Mp2*(kK*qe2*
                   (kq*(-(kq*G1(0)) + 2*kP*G2(0) + 3*kq*G2(0)) + 
                     (G1(0) - 3*G2(0))*pow(kK,2)) + 
                  2*kq*(2*kK*(2*kP - kq)*me2*G1(0) + 
                     nu*((-G1(0) + G2(0))*pow(kK,2) + G1(0)*pow(kq,2)))) + 
               G2(0)*(qe2*(-12*kK*kP*(kP + kq)*me2 + 
                     nu*((5*kP - kq)*pow(kK,2) - 
                        kq*(5*kP*kq + 2*pow(kP,2) + 3*pow(kq,2)))) + 
                  kK*kq*(-2*kP + kq)*pow(nu,2) + 
                  kK*(-4*kP*kq + 3*pow(kK,2) - pow(kP,2) - 3*pow(kq,2))*
                   pow(qe2,2))))) + 
      F2(qp2)*(2*Mp*qe2*G3(qe2)*
          (Mp*(-4*Mp2*(kK*qe2*(G1(0) - G2(0))*
                   (-(kq*(kP + 2*kq)) + 2*pow(kK,2)) + 
                  kq*(-4*kK*kq*me2*G1(0) + 
                     nu*(G1(0) - G2(0))*(pow(kK,2) - 2*pow(kq,2)))) + 
               qe2*(-4*kK*kP*(kP + kq)*me2*(G1(0) - 2*G2(0)) + 
                  nu*((2*kP*G1(0) + 2*kq*G1(0) - 3*kP*G2(0) - kq*G2(0))*
                      pow(kK,2) + 
                     kq*(kP*kq*(-2*G1(0) + 3*G2(0)) + 
                        (G1(0) - G2(0))*pow(kP,2) - 
                        (G1(0) - 2*G2(0))*pow(kq,2)))) - 
               kK*kq*(kq*G2(0) + kP*(-G1(0) + G2(0)))*pow(nu,2) + 
               kK*(kP*kq*(G1(0) + G2(0)) + (G1(0) - 2*G2(0))*pow(kK,2) + 
                  G2(0)*pow(kP,2) - (G1(0) - 2*G2(0))*pow(kq,2))*
                pow(qe2,2)) + Md*
             (4*Mp2*(kK*qe2*(G1(0) - G2(0))*
                   (-(kq*(kP + 2*kq)) + 2*pow(kK,2)) + 
                  kq*(4*kK*kq*me2*G1(0) + 
                     nu*(G1(0) - G2(0))*(pow(kK,2) - 2*pow(kq,2)))) + 
               qe2*(4*kK*kP*(kP + kq)*me2*(3*G1(0) - 2*G2(0)) + 
                  nu*((-4*kP*G1(0) + 3*kP*G2(0) + kq*G2(0))*pow(kK,2) + 
                     kq*(kP*kq*(4*G1(0) - 3*G2(0)) + 
                        (-G1(0) + G2(0))*pow(kP,2) + 
                        (3*G1(0) - 2*G2(0))*pow(kq,2)))) + 
               kK*kq*(kq*(-2*G1(0) + G2(0)) + kP*(-G1(0) + G2(0)))*
                pow(nu,2) - kK*
                (kP*kq*(-3*G1(0) + G2(0)) + 
                  (3*G1(0) - 2*G2(0))*pow(kK,2) + 
                  (-2*G1(0) + G2(0))*pow(kP,2) + 
                  (-3*G1(0) + 2*G2(0))*pow(kq,2))*pow(qe2,2))) + 
         Md*G2(qe2)*(Mp*(-(qe2*(3*G1(0) - 2*G2(0))*
                  (qe2*(-4*kK*kP*(kP + kq)*me2 + 
                       nu*(2*kP*pow(kK,2) - kq*pow(kP + kq,2))) - 
                    kK*kP*kq*pow(nu,2) + 
                    kK*(-(kq*(kP + kq)) + pow(kK,2))*pow(qe2,2))) + 
               4*Mp2*(qe2*(4*kK*me2*G1(0)*pow(kq,2) + 
                     nu*((kP*G1(0) + kq*G1(0) - kP*G2(0))*pow(kK,2) - 
                        (kP + 2*kq)*(G1(0) - G2(0))*pow(kq,2))) + 
                  kK*(-G1(0) + G2(0))*pow(kq,2)*pow(nu,2) + 
                  kK*(-(kP*kq*(G1(0) - 2*G2(0))) + 
                     2*(G1(0) - G2(0))*pow(kK,2) + 
                     (G1(0) - G2(0))*pow(kP,2) + 
                     2*(-G1(0) + G2(0))*pow(kq,2))*pow(qe2,2))) - 
            Md*(-(qe2*(G1(0) - 2*G2(0))*
                  (qe2*(-4*kK*kP*(kP + kq)*me2 + 
                       nu*(2*kP*pow(kK,2) - kq*pow(kP + kq,2))) - 
                    kK*kP*kq*pow(nu,2) + 
                    kK*(-(kq*(kP + kq)) + pow(kK,2))*pow(qe2,2))) + 
               4*Mp2*(qe2*(-4*kK*me2*G1(0)*pow(kq,2) + 
                     nu*((kP*G1(0) - kq*G1(0) - kP*G2(0))*pow(kK,2) - 
                        (kP + 2*kq)*(G1(0) - G2(0))*pow(kq,2))) + 
                  kK*(-G1(0) + G2(0))*pow(kq,2)*pow(nu,2) + 
                  kK*(kP*kq*(-3*G1(0) + 2*G2(0)) + 
                     2*(G1(0) - G2(0))*pow(kK,2) + 
                     (G1(0) - G2(0))*pow(kP,2) + 
                     2*(-G1(0) + G2(0))*pow(kq,2))*pow(qe2,2)))) - 
         Md*G1(qe2)*(Md*(64*kK*me2*G1(0)*pow(kq,2)*pow(Mp,4) + 
               3*qe2*G2(0)*(qe2*
                   (4*kK*kP*(kP + kq)*me2 + 
                     nu*(-2*kP*pow(kK,2) + kq*pow(kP + kq,2))) + 
                  kK*kP*kq*pow(nu,2) + 
                  kK*(kq*(kP + kq) - pow(kK,2))*pow(qe2,2)) - 
               4*Mp2*(qe2*(-4*kK*kP*(kP - 2*kq)*me2*G1(0) + 
                     nu*((kP - kq)*(2*G1(0) - G2(0))*pow(kK,2) + 
                        (-2*kP*G1(0) + 2*kq*G1(0) + kP*G2(0))*pow(kq,2))) \
+ kK*G2(0)*pow(kq,2)*pow(nu,2) + 
                  kK*(kP*kq*G2(0) + (G1(0) - 3*G2(0))*pow(kK,2) - 
                     G2(0)*pow(kP,2) - (G1(0) - 3*G2(0))*pow(kq,2))*
                   pow(qe2,2))) + 
            Mp*(64*kK*me2*G1(0)*pow(kq,2)*pow(Mp,4) - 
               3*qe2*(2*G1(0) - G2(0))*
                (qe2*(-4*kK*kP*(kP + kq)*me2 + 
                     nu*(2*kP*pow(kK,2) - kq*pow(kP + kq,2))) - 
                  kK*kP*kq*pow(nu,2) + 
                  kK*(-(kq*(kP + kq)) + pow(kK,2))*pow(qe2,2)) + 
               4*Mp2*(qe2*(4*kK*kP*(kP - 2*kq)*me2*G1(0) + 
                     nu*(-(kP*G2(0)*pow(kK,2)) + kq*G2(0)*pow(kK,2) + 
                        kP*G2(0)*pow(kq,2) - 2*G1(0)*pow(kq,3))) + 
                  kK*(-2*G1(0) + G2(0))*pow(kq,2)*pow(nu,2) + 
                  kK*(kP*kq*(-2*G1(0) + G2(0)) + 
                     (5*G1(0) - 3*G2(0))*pow(kK,2) + 
                     (2*G1(0) - G2(0))*pow(kP,2) + 
                     (-5*G1(0) + 3*G2(0))*pow(kq,2))*pow(qe2,2))))))*
    pow(qp2,-1)*pow(Gd2*Md2 + pow(kP + kq - Md2 + Mp2,2),-1))/36.
);
}
double iMdiMeHard2(double nu, double qe2, double qp2, double kK, double kP, double kq) {
return ((e6*Z3*pow(kK - kq,-1)*pow(kK + kq,-1)*pow(Md,-5)*pow(Mp,-1)*
    pow(-kP + kq - Md2 + Mp2,-1)*pow(qe2,-1)*
    (F2(qp2)*(Md*G1(qe2)*(Mp*(4*kK*kP*Mp2*
                ((2*kP*G1(0) - 4*kq*G1(0) - 2*kP*G2(0) - 3*kq*G2(0))*
                   pow(kK,2) + 
                  (4*kP*G1(0) - 4*kq*G1(0) - kP*G2(0) + 7*kq*G2(0))*
                   pow(kq,2)) - 
               kK*qe2*(pow(kK,2)*
                   (13*kP*kq*(2*G1(0) - G2(0)) + 
                     4*(2*G1(0) + 3*G2(0))*pow(kP,2) - 
                     2*(14*G1(0) + G2(0))*pow(kq,2)) - 
                  (kP - kq)*kq*
                   (kP*kq*(-34*G1(0) + 33*G2(0)) + 
                     7*(2*G1(0) - G2(0))*pow(kP,2) + 
                     2*(14*G1(0) + G2(0))*pow(kq,2))) - 
               4*(kP - kq)*kq*
                (-2*kK*(kP - kq)*me2*
                   (10*kP*G1(0) + 10*kq*G1(0) - 3*kP*G2(0) + kq*G2(0)) + 
                  nu*((6*kP*G1(0) + 4*kq*G1(0) - 2*kP*G2(0) + 
                        3*kq*G2(0))*pow(kK,2) + 
                     (-4*kP*G1(0) + 4*kq*G1(0) + kP*G2(0) - 7*kq*G2(0))*
                      pow(kq,2)))) + 
            Md*(-4*kK*kP*Mp2*((2*kP*G1(0) + 4*kq*G1(0) - 2*kP*G2(0) - 
                     3*kq*G2(0))*pow(kK,2) - 
                  (2*kP*G1(0) + 4*kq*G1(0) + kP*G2(0) - kq*G2(0))*
                   pow(kq,2)) + 
               kK*qe2*((kP - kq)*kq*
                   (kP*kq*(8*G1(0) + 5*G2(0)) + 9*G2(0)*pow(kP,2) + 
                     8*(G1(0) - 2*G2(0))*pow(kq,2)) + 
                  pow(kK,2)*(-3*kP*kq*G2(0) + 
                     (-8*G1(0) + 4*G2(0))*pow(kP,2) + 
                     8*(G1(0) - 2*G2(0))*pow(kq,2))) + 
               2*kq*(4*kK*kP*(kP - kq)*me2*
                   (4*kP*G1(0) + 8*kq*G1(0) + kP*G2(0) - kq*G2(0)) + 
                  nu*((kP - kq)*
                      (4*kP*G1(0) + 8*kq*G1(0) + kP*G2(0) - kq*G2(0))*
                      pow(kq,2) - 
                     pow(kK,2)*
                      (kP*kq*(4*G1(0) - 2*G2(0)) + 
                        (4*G1(0) + G2(0))*pow(kP,2) + 
                        (-8*G1(0) + 7*G2(0))*pow(kq,2)))))) + 
         Md*G2(qe2)*(Mp*(-4*kK*kP*Mp2*
                ((2*kP*(G1(0) - G2(0)) - 3*kq*(G1(0) + G2(0)))*
                   pow(kK,2) + (kP - kq)*(G1(0) - G2(0))*pow(kq,2)) - 
               2*(kP - kq)*kq*
                (4*kK*(kP - kq)*me2*
                   (3*kP*G1(0) - kq*G1(0) - 2*kP*G2(0) + 2*kq*G2(0)) + 
                  nu*((-4*kP*G1(0) - 6*kq*G1(0) + kP*G2(0) - 
                        5*kq*G2(0))*pow(kK,2) + 
                     (kP - kq)*(2*G1(0) - 3*G2(0))*pow(kq,2))) + 
               kK*qe2*(pow(kK,2)*
                   (kP*kq*(13*G1(0) - 4*G2(0)) + 
                     2*(2*G1(0) + 5*G2(0))*pow(kP,2) - 
                     2*(7*G1(0) + 4*G2(0))*pow(kq,2)) - 
                  (kP - kq)*kq*
                   (kP*kq*(-17*G1(0) + 2*G2(0)) + 
                     (7*G1(0) - 16*G2(0))*pow(kP,2) + 
                     2*(7*G1(0) + 4*G2(0))*pow(kq,2)))) + 
            Md*(4*kK*kP*Mp2*(G1(0) - G2(0))*
                ((2*kP + 3*kq)*pow(kK,2) + (kP - kq)*pow(kq,2)) + 
               kK*qe2*((kP - kq)*kq*
                   (kP*kq*(G1(0) - 2*G2(0)) + 
                     (7*G1(0) - 8*G2(0))*pow(kP,2) + 
                     12*(-G1(0) + G2(0))*pow(kq,2)) + 
                  pow(kK,2)*(kP*kq*(-7*G1(0) + 8*G2(0)) + 
                     6*(G1(0) - G2(0))*pow(kP,2) + 
                     12*(-G1(0) + G2(0))*pow(kq,2))) + 
               2*kq*(4*kK*kP*me2*(G1(0) - 2*G2(0))*pow(kP - kq,2) + 
                  nu*((G1(0) - 2*G2(0))*pow(kP - kq,2)*pow(kq,2) + 
                     pow(kK,2)*
                      (2*kP*kq*(G1(0) - 2*G2(0)) - 
                        (G1(0) - 2*G2(0))*pow(kP,2) + 
                        (-7*G1(0) + 8*G2(0))*pow(kq,2)))))) + 
         qe2*G3(qe2)*(-2*kK*Md*Mp*
             (-((kP - kq)*kq*(-(kP*kq*(4*G1(0) + 3*G2(0))) + 
                    (G1(0) - G2(0))*pow(kP,2) - 
                    2*(5*G1(0) + 3*G2(0))*pow(kq,2))) + 
               pow(kK,2)*(-(kP*kq*(12*G1(0) + G2(0))) + 
                  (G1(0) - 5*G2(0))*pow(kP,2) + 
                  2*(5*G1(0) + 3*G2(0))*pow(kq,2))) + 
            (kP - kq)*(16*kK*kP*(kP - kq)*kq*me2*(G1(0) - 2*G2(0)) + 
               nu*(pow(kK,2)*
                   (kP*kq*(-14*G1(0) + 19*G2(0)) + 
                     (-2*G1(0) + G2(0))*pow(kP,2) + 
                     2*(4*G1(0) - 5*G2(0))*pow(kq,2)) + 
                  pow(kq,2)*(kP*kq*(14*G1(0) - 19*G2(0)) + 
                     (2*G1(0) - G2(0))*pow(kP,2) - 
                     4*(G1(0) - 2*G2(0))*pow(kq,2))) + 
               kK*qe2*((-14*kq*G1(0) + 2*kP*(G1(0) - 2*G2(0)) + 
                     19*kq*G2(0))*pow(kK,2) + 
                  (-2*G1(0) + G2(0))*pow(kP,3) + 
                  kP*(4*G1(0) + G2(0))*pow(kq,2) + 
                  (14*G1(0) - 19*G2(0))*pow(kq,3))) + 
            2*kK*Mp2*(pow(kK,2)*
                (kP*kq*(6*G1(0) - 9*G2(0)) + 
                  (11*G1(0) - 9*G2(0))*pow(kP,2) + 
                  2*(-11*G1(0) + 12*G2(0))*pow(kq,2)) + 
               kq*(kq*(-7*G1(0) + 4*G2(0))*pow(kP,2) + 
                  (-G1(0) + G2(0))*pow(kP,3) + 
                  kP*(-8*G1(0) + 13*G2(0))*pow(kq,2) + 
                  2*(11*G1(0) - 12*G2(0))*pow(kq,3))))) + 
      kK*(2*kP*(kP - kq)*Md*Mp*F1(qp2)*
          (G2(qe2)*((kP - kq)*kq*
                (kP*(G1(0) - 2*G2(0)) - kq*(G1(0) + G2(0))) + 
               (5*kP*G1(0) - 11*kq*G1(0) - 7*kP*G2(0) + kq*G2(0))*
                pow(kK,2)) + G1(qe2)*
             ((2*kP*G1(0) + 10*kq*G1(0) + 5*kP*G2(0) + kq*G2(0))*
                pow(kK,2) + kq*
                (-2*kP*kq*(4*G1(0) + G2(0)) + 
                  (-2*G1(0) + G2(0))*pow(kP,2) + 
                  (10*G1(0) - 11*G2(0))*pow(kq,2)))) + 
         F2(qp2)*(2*kP*Md*G1(qe2)*
             ((kP - kq)*Mp*((2*kP*G1(0) + 8*kq*G1(0) + kP*G2(0) + 
                     6*kq*G2(0))*pow(kK,2) + 
                  2*(-4*kP*G1(0) + 4*kq*G1(0) + kP*G2(0) - 7*kq*G2(0))*
                   pow(kq,2)) + 
               Md*(-((kP - kq)*
                     (4*kP*G1(0) + 8*kq*G1(0) + kP*G2(0) - kq*G2(0))*
                     pow(kq,2)) + 
                  pow(kK,2)*(kP*kq*(4*G1(0) - 2*G2(0)) + 
                     (4*G1(0) - 2*G2(0))*pow(kP,2) + 
                     (-8*G1(0) + 7*G2(0))*pow(kq,2)))) - 
            2*kP*Md*G2(qe2)*(-((kP - kq)*Mp*
                  ((kP*G1(0) - 6*kq*G1(0) - 3*kP*G2(0) - 5*kq*G2(0))*
                     pow(kK,2) + (kP - kq)*(2*G1(0) - 3*G2(0))*pow(kq,2))\
) + Md*((G1(0) - 2*G2(0))*pow(kP - kq,2)*pow(kq,2) + 
                  pow(kK,2)*(2*kP*kq*(G1(0) - 2*G2(0)) + 
                     (2*G1(0) - G2(0))*pow(kP,2) + 
                     (-7*G1(0) + 8*G2(0))*pow(kq,2)))) + 
            (kP - kq)*qe2*G3(qe2)*
             (pow(kK,2)*(-4*kP*kq*(G1(0) - 2*G2(0)) + 
                  (-6*G1(0) + 3*G2(0))*pow(kP,2) + 
                  2*(10*G1(0) - 11*G2(0))*pow(kq,2)) - 
               kq*(3*kq*(-2*G1(0) + G2(0))*pow(kP,2) + 
                  (-2*G1(0) + G2(0))*pow(kP,3) + 6*kP*G2(0)*pow(kq,2) + 
                  2*(10*G1(0) - 11*G2(0))*pow(kq,3))))) + 
      2*Mp*F1(qp2)*(kK*(kP - kq)*qe2*G3(qe2)*
          (Mp*((12*kP*G1(0) - 24*kq*G1(0) - 15*kP*G2(0) + 26*kq*G2(0))*
                pow(kK,2) + kq*
                (-4*kP*kq*G1(0) + 9*kP*kq*G2(0) + G2(0)*pow(kP,2) + 
                  24*G1(0)*pow(kq,2) - 26*G2(0)*pow(kq,2))) + 
            Md*((-6*kP*G1(0) + 16*kq*G1(0) + 9*kP*G2(0) - 2*kq*G2(0))*
                pow(kK,2) + kq*
                (kP*kq*(2*G1(0) - 7*G2(0)) + 
                  (2*G1(0) - G2(0))*pow(kP,2) + 
                  2*(-8*G1(0) + G2(0))*pow(kq,2)))) + 
         Md*G2(qe2)*(-2*kK*kP*Mp2*
             ((kP - kq)*kq*(2*kP*(G1(0) - G2(0)) + 
                  kq*(-2*G1(0) + G2(0))) + 
               (9*kP*G1(0) - 12*kq*G1(0) - 10*kP*G2(0) + 7*kq*G2(0))*
                pow(kK,2)) + 2*kK*kP*Md*Mp*
             ((kP - kq)*kq*(2*kP*(G1(0) - G2(0)) + kq*G2(0)) + 
               (5*kP*G1(0) - 8*kq*G1(0) - 6*kP*G2(0) + 9*kq*G2(0))*
                pow(kK,2)) - (kP - kq)*
             (4*kK*(kP + kq)*me2*(G1(0) - 2*G2(0))*pow(kP - kq,2) + 
               nu*((kP - kq)*
                   (kP*(G1(0) - 2*G2(0)) - kq*(G1(0) + G2(0)))*pow(kq,2) \
+ pow(kK,2)*(-3*kP*kq*G2(0) - (G1(0) - 2*G2(0))*pow(kP,2) + 
                     (-11*G1(0) + G2(0))*pow(kq,2))) + 
               kK*qe2*(-3*(3*kq*G1(0) + kP*(G1(0) + G2(0)))*pow(kK,2) + 
                  kq*(2*G1(0) - 7*G2(0))*pow(kP,2) + 
                  (G1(0) - 2*G2(0))*pow(kP,3) + 
                  6*kP*(-2*G1(0) + G2(0))*pow(kq,2) + 9*G1(0)*pow(kq,3)))) \
+ Md*G1(qe2)*(-2*kK*kP*Md*Mp*(-2*(kP - kq)*kq*(4*kq*G1(0) + kP*G2(0)) + 
               (8*kP*G1(0) - 8*kq*G1(0) - 5*kP*G2(0) + 8*kq*G2(0))*
                pow(kK,2)) + 2*kK*kP*Mp2*
             ((4*kP*G1(0) - 10*kq*G1(0) - 9*kP*G2(0) + 6*kq*G2(0))*
                pow(kK,2) + 2*kq*
                (kP*kq*(3*G1(0) + 2*G2(0)) + 
                  (2*G1(0) - G2(0))*pow(kP,2) + 
                  (-5*G1(0) + 2*G2(0))*pow(kq,2))) - 
            (kP - kq)*(-4*kK*me2*
                (2*kP*G1(0) + 10*kq*G1(0) - kP*G2(0) + kq*G2(0))*
                (pow(kP,2) - pow(kq,2)) + 
               nu*((kP + kq)*(2*kP*G1(0) + 10*kq*G1(0) - kP*G2(0) + 
                     kq*G2(0))*pow(kK,2) + 
                  pow(kq,2)*(-2*kP*kq*(4*G1(0) + G2(0)) + 
                     (-2*G1(0) + G2(0))*pow(kP,2) + 
                     (10*G1(0) - 11*G2(0))*pow(kq,2))) + 
               kK*qe2*((2*kP*G1(0) + 22*kq*G1(0) + 5*kP*G2(0) - 
                     5*kq*G2(0))*pow(kK,2) + 
                  2*kq*(-2*G1(0) + G2(0))*pow(kP,2) + 
                  (-2*G1(0) + G2(0))*pow(kP,3) + 
                  4*kP*(7*G1(0) - 5*G2(0))*pow(kq,2) + 
                  (-22*G1(0) + 5*G2(0))*pow(kq,3))))) - 
      2*Mp*F1(qp2)*(qe2*G3(qe2)*
          (-(Md*(4*kK*Mp2*((kP*(7*G1(0) - 8*G2(0)) + 
                       5*kq*(-2*G1(0) + G2(0)))*pow(kK,2) + 
                    (-7*kP*G1(0) + 10*kq*G1(0) + 7*kP*G2(0) - 
                       5*kq*G2(0))*pow(kq,2)) - 
                 (kP - kq)*(kK*qe2*(2*G1(0) - G2(0))*
                     (-3*kP*kq + 5*pow(kK,2) + pow(kP,2) - 5*pow(kq,2)) \
- 8*kK*me2*(G1(0) - 2*G2(0))*(pow(kP,2) - pow(kq,2)) + 
                    nu*((6*kP*G1(0) + 2*kq*G1(0) - 9*kP*G2(0) - 
                        7*kq*G2(0))*pow(kK,2) - 
                       3*(kP + kq)*(2*G1(0) - 3*G2(0))*pow(kq,2))))) - 
            Mp*(4*kK*Mp2*((-11*kP*G1(0) + 14*kq*G1(0) + 12*kP*G2(0) - 
                     15*kq*G2(0))*pow(kK,2) + 
                  (9*kP*G1(0) - 14*kq*G1(0) - 11*kP*G2(0) + 15*kq*G2(0))*
                   pow(kq,2)) - 
               (kP - kq)*(16*kK*(kP - kq)*me2*
                   (kq*G1(0) + 2*kP*(G1(0) - G2(0))) + 
                  kK*qe2*(kP*kq*(8*G1(0) - 7*G2(0)) + 
                     (-4*G1(0) + 7*G2(0))*pow(kK,2) + G2(0)*pow(kP,2) + 
                     (4*G1(0) - 7*G2(0))*pow(kq,2)) + 
                  nu*((-16*kP*G1(0) + 4*kq*G1(0) + 17*kP*G2(0) - 
                        7*kq*G2(0))*pow(kK,2) + 
                     (16*kP*G1(0) + 4*kq*G1(0) - 17*kP*G2(0) + kq*G2(0))*
                      pow(kq,2))))) + 
         Md*G1(qe2)*(-2*Md*Mp*
             (8*kK*kP*Mp2*((G1(0) - G2(0))*pow(kK,2) - 
                  G1(0)*pow(kq,2)) + 
               kK*qe2*((4*kP*G1(0) - 4*kq*G1(0) - 3*kP*G2(0) + 
                     9*kq*G2(0))*pow(kK,2) - 
                  (kP - kq)*(6*kP*kq*G2(0) + 2*G2(0)*pow(kP,2) + 
                     (4*G1(0) - 9*G2(0))*pow(kq,2))) + 
               2*(-2*kK*me2*(4*kq*G1(0) + kP*G2(0) - kq*G2(0))*
                   (-(kP*kq) + 2*pow(kP,2) - pow(kq,2)) + 
                  nu*(-((kP - kq)*(4*kq*G1(0) + kP*G2(0))*pow(kq,2)) + 
                     pow(kK,2)*
                      (kP*kq*(4*G1(0) - 2*G2(0)) + G2(0)*pow(kP,2) + 
                        4*(-G1(0) + G2(0))*pow(kq,2))))) + 
            2*Mp2*(2*(2*kK*(kP - kq)*me2*
                   (kP*kq*(8*G1(0) + 3*G2(0)) + 
                     (4*G1(0) - 2*G2(0))*pow(kP,2) + 
                     (12*G1(0) - G2(0))*pow(kq,2)) + 
                  nu*(pow(kq,2)*
                      (kP*kq*(3*G1(0) + 2*G2(0)) + 
                        (2*G1(0) - G2(0))*pow(kP,2) + 
                        (-5*G1(0) + 2*G2(0))*pow(kq,2)) + 
                     pow(kK,2)*
                      (kP*kq*(G1(0) - 4*G2(0)) + 
                        (-2*G1(0) + G2(0))*pow(kP,2) + 
                        (-5*G1(0) + 3*G2(0))*pow(kq,2)))) + 
               kK*qe2*((16*kP*G1(0) - 28*kq*G1(0) - 13*kP*G2(0) + 
                     13*kq*G2(0))*pow(kK,2) + 
                  6*kq*(2*G1(0) - G2(0))*pow(kP,2) + 
                  (4*G1(0) - 2*G2(0))*pow(kP,3) + 
                  kP*(-44*G1(0) + 27*G2(0))*pow(kq,2) + 
                  (28*G1(0) - 13*G2(0))*pow(kq,3))) + 
            16*kK*kP*((G1(0) - G2(0))*pow(kK,2) + G1(0)*pow(kq,2))*
             pow(Mp,4) - (kP - kq)*(2*G1(0) - G2(0))*
             (qe2*(nu*(kP*(kP - kq)*kq + (7*kP + 5*kq)*pow(kK,2)) + 
                  12*kK*me2*(-pow(kP,2) + pow(kq,2))) + 
               kK*kq*(kP + 5*kq)*pow(nu,2) + 
               kK*(5*kP*kq + 6*pow(kK,2) + pow(kP,2) - 6*pow(kq,2))*
                pow(qe2,2))) + 
         Md*G2(qe2)*(-2*Mp2*(4*kK*me2*
                (-(kq*G1(0)) + 2*kP*(G1(0) - G2(0)))*pow(kP - kq,2) + 
               nu*((kP - kq)*
                   (2*kP*(G1(0) - G2(0)) + kq*(-2*G1(0) + G2(0)))*
                   pow(kq,2) + 
                  pow(kK,2)*(kP*kq*(8*G1(0) - 9*G2(0)) - 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     (-12*G1(0) + 7*G2(0))*pow(kq,2))) + 
               kK*qe2*((5*kP*G1(0) - 11*kq*G1(0) - 7*kP*G2(0) + 
                     6*kq*G2(0))*pow(kK,2) + 
                  kq*(6*G1(0) - 8*G2(0))*pow(kP,2) + 
                  2*(G1(0) - G2(0))*pow(kP,3) + 
                  kP*(-19*G1(0) + 15*G2(0))*pow(kq,2) + 
                  (11*G1(0) - 6*G2(0))*pow(kq,3))) + 
            2*Md*Mp*(4*kP*kq*nu*G1(0)*pow(kK,2) - 
               5*kP*kq*nu*G2(0)*pow(kK,2) + 
               8*kP*Mp2*(G1(0) - G2(0))*pow(kK,3) - 
               12*kK*kq*me2*G1(0)*pow(kP,2) + 
               16*kK*kq*me2*G2(0)*pow(kP,2) - 
               2*nu*G1(0)*pow(kK,2)*pow(kP,2) + 
               2*nu*G2(0)*pow(kK,2)*pow(kP,2) + 
               8*kK*me2*G1(0)*pow(kP,3) - 8*kK*me2*G2(0)*pow(kP,3) - 
               8*kK*kP*me2*G2(0)*pow(kq,2) - 
               8*nu*G1(0)*pow(kK,2)*pow(kq,2) + 
               9*nu*G2(0)*pow(kK,2)*pow(kq,2) + 
               2*nu*G1(0)*pow(kP,2)*pow(kq,2) - 
               2*nu*G2(0)*pow(kP,2)*pow(kq,2) + 
               4*kK*me2*G1(0)*pow(kq,3) - 2*kP*nu*G1(0)*pow(kq,3) + 
               3*kP*nu*G2(0)*pow(kq,3) + 
               kK*qe2*((-7*kq*G1(0) + 3*kP*(G1(0) - G2(0)) + 
                     8*kq*G2(0))*pow(kK,2) + 
                  4*kq*(G1(0) - G2(0))*pow(kP,2) + 
                  2*(G1(0) - G2(0))*pow(kP,3) + 
                  kP*(-11*G1(0) + 13*G2(0))*pow(kq,2) + 
                  (7*G1(0) - 8*G2(0))*pow(kq,3)) - nu*G2(0)*pow(kq,4)) - 
            16*kP*(G1(0) - G2(0))*pow(kK,3)*pow(Mp,4) + 
            (kP - kq)*(qe2*(-4*kK*me2*(G1(0) - 2*G2(0))*
                   (pow(kP,2) - pow(kq,2)) + 
                  nu*(3*(kP*G1(0) + 3*kq*G1(0) - kP*G2(0) - 2*kq*G2(0))*
                      pow(kK,2) + 
                     kq*(3*kP*kq*G1(0) + (G1(0) - 2*G2(0))*pow(kP,2) + 
                        (-4*G1(0) + 5*G2(0))*pow(kq,2)))) + 
               kK*kq*(kP*G1(0) + 5*kq*G1(0) - 2*kP*G2(0) - 4*kq*G2(0))*
                pow(nu,2) + kK*
                (kP*kq*(5*G1(0) - G2(0)) + 
                  (6*G1(0) - 3*G2(0))*pow(kK,2) + 
                  (G1(0) + G2(0))*pow(kP,2) + 
                  3*(-2*G1(0) + G2(0))*pow(kq,2))*pow(qe2,2)))) + 
      F2(qp2)*(qe2*G3(qe2)*(-2*Mp2*
             (8*kK*(kP - kq)*kq*me2*
                (-2*kq*G1(0) + kP*(G1(0) - 2*G2(0))) + 
               nu*(pow(kK,2)*
                   (kP*kq*(-11*G1(0) + 12*G2(0)) + 
                     (-G1(0) + G2(0))*pow(kP,2) + 
                     (8*G1(0) - 7*G2(0))*pow(kq,2)) + 
                  pow(kq,2)*(kP*kq*(15*G1(0) - 16*G2(0)) + 
                     (G1(0) - G2(0))*pow(kP,2) + 
                     (-10*G1(0) + 11*G2(0))*pow(kq,2))) + 
               kK*qe2*((9*kP*G1(0) - 15*kq*G1(0) - 11*kP*G2(0) + 
                     19*kq*G2(0))*pow(kK,2) + 
                  kq*(4*G1(0) - 6*G2(0))*pow(kP,2) + 
                  (-G1(0) + G2(0))*pow(kP,3) - 
                  2*kP*(5*G1(0) - 8*G2(0))*pow(kq,2) + 
                  (15*G1(0) - 19*G2(0))*pow(kq,3))) + 
            2*Md*Mp*(-7*kP*kq*nu*G1(0)*pow(kK,2) + 
               10*kP*kq*nu*G2(0)*pow(kK,2) + 
               24*kK*kq*me2*G1(0)*pow(kP,2) - 
               16*kK*kq*me2*G2(0)*pow(kP,2) - 
               nu*G1(0)*pow(kK,2)*pow(kP,2) + 
               nu*G2(0)*pow(kK,2)*pow(kP,2) - 
               32*kK*kP*me2*G1(0)*pow(kq,2) + 
               32*kK*kP*me2*G2(0)*pow(kq,2) + 
               6*nu*G1(0)*pow(kK,2)*pow(kq,2) - 
               11*nu*G2(0)*pow(kK,2)*pow(kq,2) + 
               nu*G1(0)*pow(kP,2)*pow(kq,2) - 
               nu*G2(0)*pow(kP,2)*pow(kq,2) + 
               4*kK*Mp2*((2*kP*(G1(0) - G2(0)) - 3*kq*(G1(0) + G2(0)))*
                   pow(kK,2) + 
                  (kP*(-G1(0) + G2(0)) + 3*kq*(G1(0) + G2(0)))*pow(kq,2)\
) + 8*kK*me2*G1(0)*pow(kq,3) + 11*kP*nu*G1(0)*pow(kq,3) - 
               16*kK*me2*G2(0)*pow(kq,3) - 12*kP*nu*G2(0)*pow(kq,3) - 
               kK*qe2*((kP*(5*G1(0) + G2(0)) - kq*(G1(0) + 3*G2(0)))*
                   pow(kK,2) + 4*kq*(-2*G1(0) + G2(0))*pow(kP,2) + 
                  (G1(0) - G2(0))*pow(kP,3) + 
                  4*kP*(G1(0) - G2(0))*pow(kq,2) + 
                  (G1(0) + 3*G2(0))*pow(kq,3)) - 12*nu*G1(0)*pow(kq,4) + 
               13*nu*G2(0)*pow(kq,4)) - 
            8*kK*(G1(0) - G2(0))*
             ((2*kP + 3*kq)*pow(kK,2) - (kP + 3*kq)*pow(kq,2))*pow(Mp,4) \
+ (kP - kq)*(qe2*(-8*kK*kP*(kP - kq)*me2*(G1(0) - 2*G2(0)) + 
                  nu*((-4*kq*G1(0) + 4*kP*(G1(0) - 2*G2(0)) + 
                        5*kq*G2(0))*pow(kK,2) + 
                     kq*(-4*kP*kq*(G1(0) - 2*G2(0)) + 
                        (-2*G1(0) + G2(0))*pow(kP,2) + 
                        2*(G1(0) - 2*G2(0))*pow(kq,2)))) + 
               kK*kP*kq*(-2*G1(0) + G2(0))*pow(nu,2) + 
               kK*(kq*(-2*kP*G1(0) - 2*kq*G1(0) + kP*G2(0) + 
                     4*kq*G2(0)) + 2*(G1(0) - 2*G2(0))*pow(kK,2))*
                pow(qe2,2))) + 
         Md*G1(qe2)*(Mp*(qe2*
                (-4*kK*me2*(10*kP*G1(0) - 6*kq*G1(0) - 3*kP*G2(0) + 
                     3*kq*G2(0))*(kP*kq + pow(kP,2) - 2*pow(kq,2)) + 
                  nu*((kP - kq)*
                      (6*kP*G1(0) - 14*kq*G1(0) - 7*kP*G2(0) + 
                        3*kq*G2(0))*pow(kq,2) + 
                     pow(kK,2)*
                      (6*kP*kq*(2*G1(0) - G2(0)) + 
                        (18*G1(0) - 5*G2(0))*pow(kP,2) + 
                        (-18*G1(0) + 5*G2(0))*pow(kq,2)))) - 
               4*Mp2*(kK*qe2*
                   ((10*kP*G1(0) - 10*kq*G1(0) - 9*kP*G2(0) + 
                        kq*G2(0))*pow(kK,2) + 
                     kq*(-4*kP*kq*(4*G1(0) - 3*G2(0)) + 
                        (8*G1(0) - 4*G2(0))*pow(kP,2) + 
                        (10*G1(0) - G2(0))*pow(kq,2))) - 
                  kq*(-4*kK*(kP - kq)*me2*
                      (2*kP*G1(0) + 12*kq*G1(0) - kq*G2(0)) + 
                     nu*((4*kq*G1(0) + kP*G2(0) + 3*kq*G2(0))*
                        pow(kK,2) + 
                        (-4*kP*G1(0) + 4*kq*G1(0) + kP*G2(0) - 
                        7*kq*G2(0))*pow(kq,2)))) + 
               10*kK*(kP - kq)*(2*G1(0) - G2(0))*pow(kq,2)*pow(nu,2) + 
               kK*(2*(8*kP*G1(0) - 2*kq*G1(0) - 3*kP*G2(0))*pow(kK,2) - 
                  (kP - kq)*(5*kP*kq*(-2*G1(0) + G2(0)) + 
                     (2*G1(0) - G2(0))*pow(kP,2) + 4*G1(0)*pow(kq,2)))*
                pow(qe2,2)) + 
            Md*(qe2*(-4*kK*kP*(kP - kq)*me2*
                   (4*kP*G1(0) + 8*kq*G1(0) + kP*G2(0) + 5*kq*G2(0)) + 
                  nu*(-((kP - kq)*
                        (8*kP*G1(0) + 8*kq*G1(0) - 7*kP*G2(0) + 
                        5*kq*G2(0))*pow(kq,2)) + 
                     pow(kK,2)*
                      (12*kP*kq*G2(0) + (8*G1(0) - G2(0))*pow(kP,2) + 
                        (-8*G1(0) + 7*G2(0))*pow(kq,2)))) + 
               4*Mp2*(kK*qe2*
                   ((2*kP*G1(0) + 2*kq*G1(0) - 3*kP*G2(0) - 
                        5*kq*G2(0))*pow(kK,2) - 
                     kq*(2*kP*kq*G1(0) + 4*G2(0)*pow(kP,2) + 
                        (2*G1(0) - 5*G2(0))*pow(kq,2))) - 
                  kq*(4*kK*me2*
                      (kP*kq*(8*G1(0) - G2(0)) + 2*G1(0)*pow(kP,2) + 
                        (-4*G1(0) + G2(0))*pow(kq,2)) + 
                     nu*((-2*kP*G1(0) - 4*kq*G1(0) + kP*G2(0) + 
                        3*kq*G2(0))*pow(kK,2) + 
                        (2*kP*G1(0) + 4*kq*G1(0) + kP*G2(0) - kq*G2(0))*
                         pow(kq,2)))) + 
               6*kK*kP*G2(0)*pow(kq,2)*pow(nu,2) + 
               kK*(2*(2*kP*G1(0) - 2*kq*G1(0) - kP*G2(0) + 7*kq*G2(0))*
                   pow(kK,2) - 
                  (kP - kq)*(3*kP*kq*G2(0) + 3*G2(0)*pow(kP,2) + 
                     2*(2*G1(0) - 7*G2(0))*pow(kq,2)))*pow(qe2,2))) + 
         Md*G2(qe2)*(Md*(nu*qe2*
                ((kP - kq)*(5*kP*G1(0) - kq*G1(0) - 4*kP*G2(0) + 
                     2*kq*G2(0))*pow(kq,2) + 
                  pow(kK,2)*(6*kP*kq*(G1(0) - G2(0)) + 
                     (G1(0) - 2*G2(0))*pow(kP,2) + 
                     (7*G1(0) - 8*G2(0))*pow(kq,2))) + 
               kK*qe2*(-4*kP*me2*(G1(0) - 2*G2(0))*
                   (pow(kP,2) - pow(kq,2)) + 
                  qe2*(-((G1(0) - 2*G2(0))*pow(kP,3)) + 
                     2*kq*(4*G1(0) - 5*G2(0))*(pow(kK,2) - pow(kq,2)) + 
                     3*kP*(3*G1(0) - 4*G2(0))*pow(kq,2))) - 
               4*Mp2*(kq*(4*kK*kq*(-kP + kq)*me2*G1(0) + 
                     nu*(G1(0) - G2(0))*
                      ((kP + 3*kq)*pow(kK,2) + (kP - kq)*pow(kq,2))) + 
                  kK*qe2*((5*kP*G1(0) + 3*kq*G1(0) - 5*kP*G2(0) - 
                        2*kq*G2(0))*pow(kK,2) + 
                     kq*(kP*kq*G2(0) + 4*(G1(0) - G2(0))*pow(kP,2) + 
                        (-3*G1(0) + 2*G2(0))*pow(kq,2)))) + 
               6*kK*kP*(G1(0) - G2(0))*pow(kq,2)*pow(nu,2)) + 
            Mp*(4*Mp2*(-(kq*(4*kK*kq*(-kP + kq)*me2*G1(0) + 
                       nu*((kP*(-G1(0) + G2(0)) + 3*kq*(G1(0) + G2(0)))*
                         pow(kK,2) - (kP - kq)*(G1(0) - G2(0))*pow(kq,2))\
)) + kK*qe2*((5*kP*G1(0) - 5*kq*G1(0) - 5*kP*G2(0) - 2*kq*G2(0))*
                      pow(kK,2) + 
                     kq*(kP*kq*(-8*G1(0) + G2(0)) + 
                        4*(G1(0) - G2(0))*pow(kP,2) + 
                        (5*G1(0) + 2*G2(0))*pow(kq,2)))) + 
               qe2*(nu*(-((kP - kq)*
                        (7*kP*G1(0) - 11*kq*G1(0) - 6*kP*G2(0) + 
                        12*kq*G2(0))*pow(kq,2)) + 
                     pow(kK,2)*
                      (-2*kP*kq*(7*G1(0) - 5*G2(0)) - 
                        5*G1(0)*pow(kP,2) + 
                        (13*G1(0) - 6*G2(0))*pow(kq,2))) + 
                  4*kK*me2*((3*G1(0) - 2*G2(0))*pow(kP,3) + 
                     kP*(-5*G1(0) + 6*G2(0))*pow(kq,2) + 
                     2*(G1(0) - 2*G2(0))*pow(kq,3))) - 
               2*kK*(kP - kq)*(5*G1(0) - 4*G2(0))*pow(kq,2)*pow(nu,2) + 
               kK*(2*(kP*(-5*G1(0) + G2(0)) + kq*(2*G1(0) + G2(0)))*
                   pow(kK,2) + 
                  (kP - kq)*(-(kP*kq*(5*G1(0) + 4*G2(0))) + 
                     (G1(0) - 4*G2(0))*pow(kP,2) + 
                     2*(2*G1(0) + G2(0))*pow(kq,2)))*pow(qe2,2)))) + 
      2*Mp*F1(qp2)*(qe2*G3(qe2)*
          (Md*(-((kP - kq)*(kq*nu + kK*qe2)*(kK*nu + kP*qe2)*
                  (2*G1(0) - G2(0))) + 
               4*Mp2*(-4*kK*(kP - kq)*me2*
                   (-(kq*G1(0)) + 2*kP*(G1(0) - G2(0))) + 
                  kK*qe2*(kq*(-3*kq*G1(0) + kP*G2(0) + 2*kq*G2(0)) + 
                     (3*G1(0) - 2*G2(0))*pow(kK,2)) + 
                  nu*((4*kP - kq)*(G1(0) - G2(0))*pow(kK,2) + 
                     (-4*kP*G1(0) + kq*G1(0) + 4*kP*G2(0))*pow(kq,2))) - 
               32*kK*(G1(0) - G2(0))*(pow(kK,2) - pow(kq,2))*pow(Mp,4)) + 
            Mp*(-((kP - kq)*(kq*nu + kK*qe2)*(kK*nu + kP*qe2)*G2(0)) - 
               4*Mp2*(-4*kK*(kP - kq)*me2*
                   (kq*G1(0) + 2*kP*(G1(0) - G2(0))) + 
                  kK*qe2*(kq*
                      (-2*kP*G1(0) - kq*G1(0) + kP*G2(0) + 2*kq*G2(0)) \
+ (G1(0) - 2*G2(0))*pow(kK,2)) + 
                  nu*((4*kP - kq)*(G1(0) - G2(0))*pow(kK,2) - 
                     (4*kP*G1(0) + kq*G1(0) - 4*kP*G2(0))*pow(kq,2))) + 
               32*kK*(G1(0) - G2(0))*(pow(kK,2) - pow(kq,2))*pow(Mp,4))) + 
         2*Md*Mp*G2(qe2)*(Md*(8*kK*Mp2*(G1(0) - G2(0))*
                (kK*kq*nu + qe2*((kP - kq)*kq + pow(kK,2))) + 
               qe2*(4*kK*(kP - kq)*me2*
                   (kq*G1(0) + 2*kP*(G1(0) - G2(0))) + 
                  nu*((-3*kP*G1(0) - kq*G1(0) + 3*kP*G2(0) + 
                        2*kq*G2(0))*pow(kK,2) + 
                     kq*(3*kP*kq*(G1(0) - G2(0)) - 
                        2*(G1(0) - G2(0))*pow(kP,2) + G1(0)*pow(kq,2)))) \
- kK*kq*(2*kP + kq)*(G1(0) - G2(0))*pow(nu,2) + 
               kK*(2*kP*kq*G2(0) - (G1(0) - 2*G2(0))*pow(kK,2) + 
                  (G1(0) - G2(0))*pow(kP,2) + 
                  (G1(0) - 2*G2(0))*pow(kq,2))*pow(qe2,2)) + 
            Mp*(-8*kK*Mp2*(G1(0) - G2(0))*
                (kK*kq*nu + qe2*((kP - kq)*kq + pow(kK,2))) + 
               qe2*(-4*kK*(kP - kq)*me2*
                   (-(kq*G1(0)) + 2*kP*(G1(0) - G2(0))) + 
                  nu*((3*kP*G1(0) + 3*kq*G1(0) - 3*kP*G2(0) - 
                        2*kq*G2(0))*pow(kK,2) + 
                     kq*(3*kP*kq*(-G1(0) + G2(0)) + 
                        2*(G1(0) - G2(0))*pow(kP,2) + G1(0)*pow(kq,2)))) \
+ kK*kq*(2*kP + kq)*(G1(0) - G2(0))*pow(nu,2) + 
               kK*(kP*kq*(4*G1(0) - 2*G2(0)) + 
                  (3*G1(0) - 2*G2(0))*pow(kK,2) + 
                  (-G1(0) + G2(0))*pow(kP,2) + 
                  (-3*G1(0) + 2*G2(0))*pow(kq,2))*pow(qe2,2))) - 
         2*Md*Mp*G1(qe2)*(-(Mp*
               (4*Mp2*(kK*qe2*
                     (kq*(4*kP*G1(0) - 5*kq*G1(0) - 2*kP*G2(0) + 
                        3*kq*G2(0)) + (5*G1(0) - 3*G2(0))*pow(kK,2)) + 
                    2*kq*(2*kK*(2*kP + kq)*me2*G1(0) + 
                       nu*((G1(0) - G2(0))*pow(kK,2) + G1(0)*pow(kq,2)))) \
- (2*G1(0) - G2(0))*(qe2*(12*kK*kP*(-kP + kq)*me2 + 
                       nu*((5*kP + kq)*pow(kK,2) + 
                         kq*(-5*kP*kq + 2*pow(kP,2) + 3*pow(kq,2)))) + 
                    kK*kq*(2*kP + kq)*pow(nu,2) + 
                    kK*(4*kP*kq + 3*pow(kK,2) - pow(kP,2) - 3*pow(kq,2))*
                     pow(qe2,2)))) + 
            Md*(4*Mp2*(kK*qe2*
                   (-(kq*(kq*G1(0) + 2*kP*G2(0) - 3*kq*G2(0))) + 
                     (G1(0) - 3*G2(0))*pow(kK,2)) - 
                  2*kq*(2*kK*(2*kP + kq)*me2*G1(0) + 
                     nu*((-G1(0) + G2(0))*pow(kK,2) + G1(0)*pow(kq,2)))) + 
               G2(0)*(qe2*(12*kK*kP*(-kP + kq)*me2 + 
                     nu*((5*kP + kq)*pow(kK,2) + 
                        kq*(-5*kP*kq + 2*pow(kP,2) + 3*pow(kq,2)))) + 
                  kK*kq*(2*kP + kq)*pow(nu,2) + 
                  kK*(4*kP*kq + 3*pow(kK,2) - pow(kP,2) - 3*pow(kq,2))*
                   pow(qe2,2))))) + 
      F2(qp2)*(2*Mp*qe2*G3(qe2)*
          (Mp*(4*Mp2*(kK*qe2*(G1(0) - G2(0))*
                   ((kP - 2*kq)*kq + 2*pow(kK,2)) + 
                  kq*(-4*kK*kq*me2*G1(0) - 
                     nu*(G1(0) - G2(0))*(pow(kK,2) - 2*pow(kq,2)))) - 
               qe2*(-4*kK*kP*(kP - kq)*me2*(G1(0) - 2*G2(0)) + 
                  nu*((2*kP*G1(0) - 2*kq*G1(0) - 3*kP*G2(0) + kq*G2(0))*
                      pow(kK,2) + 
                     kq*(kP*kq*(-2*G1(0) + 3*G2(0)) + 
                        (-G1(0) + G2(0))*pow(kP,2) + 
                        (G1(0) - 2*G2(0))*pow(kq,2)))) + 
               kK*kq*(kP*(G1(0) - G2(0)) + kq*G2(0))*pow(nu,2) + 
               kK*(kP*kq*(G1(0) + G2(0)) - (G1(0) - 2*G2(0))*pow(kK,2) - 
                  G2(0)*pow(kP,2) + (G1(0) - 2*G2(0))*pow(kq,2))*
                pow(qe2,2)) + Md*
             (-4*Mp2*(kK*qe2*(G1(0) - G2(0))*
                   ((kP - 2*kq)*kq + 2*pow(kK,2)) + 
                  kq*(4*kK*kq*me2*G1(0) - 
                     nu*(G1(0) - G2(0))*(pow(kK,2) - 2*pow(kq,2)))) + 
               qe2*(-4*kK*kP*(kP - kq)*me2*(3*G1(0) - 2*G2(0)) + 
                  nu*((4*kP*G1(0) - 3*kP*G2(0) + kq*G2(0))*pow(kK,2) + 
                     kq*(kP*kq*(-4*G1(0) + 3*G2(0)) + 
                        (-G1(0) + G2(0))*pow(kP,2) + 
                        (3*G1(0) - 2*G2(0))*pow(kq,2)))) + 
               kK*kq*(-(kP*G1(0)) + 2*kq*G1(0) + kP*G2(0) - kq*G2(0))*
                pow(nu,2) + kK*
                (kP*kq*(3*G1(0) - G2(0)) + 
                  (3*G1(0) - 2*G2(0))*pow(kK,2) + 
                  (-2*G1(0) + G2(0))*pow(kP,2) + 
                  (-3*G1(0) + 2*G2(0))*pow(kq,2))*pow(qe2,2))) + 
         Md*G2(qe2)*(Mp*(qe2*(3*G1(0) - 2*G2(0))*
                (qe2*(4*kK*kP*(-kP + kq)*me2 + 
                     nu*(2*kP*pow(kK,2) + kq*pow(kP - kq,2))) + 
                  kK*kP*kq*pow(nu,2) + 
                  kK*((kP - kq)*kq + pow(kK,2))*pow(qe2,2)) - 
               4*Mp2*(qe2*(4*kK*me2*G1(0)*pow(kq,2) + 
                     nu*((kP*G1(0) - kq*G1(0) - kP*G2(0))*pow(kK,2) - 
                        (kP - 2*kq)*(G1(0) - G2(0))*pow(kq,2))) + 
                  kK*(-G1(0) + G2(0))*pow(kq,2)*pow(nu,2) + 
                  kK*(kP*kq*(G1(0) - 2*G2(0)) + 
                     2*(G1(0) - G2(0))*pow(kK,2) + 
                     (G1(0) - G2(0))*pow(kP,2) + 
                     2*(-G1(0) + G2(0))*pow(kq,2))*pow(qe2,2))) + 
            Md*(-(qe2*(G1(0) - 2*G2(0))*
                  (qe2*(4*kK*kP*(-kP + kq)*me2 + 
                       nu*(2*kP*pow(kK,2) + kq*pow(kP - kq,2))) + 
                    kK*kP*kq*pow(nu,2) + 
                    kK*((kP - kq)*kq + pow(kK,2))*pow(qe2,2))) + 
               4*Mp2*(qe2*(-4*kK*me2*G1(0)*pow(kq,2) + 
                     nu*((kP*G1(0) + kq*G1(0) - kP*G2(0))*pow(kK,2) - 
                        (kP - 2*kq)*(G1(0) - G2(0))*pow(kq,2))) + 
                  kK*(-G1(0) + G2(0))*pow(kq,2)*pow(nu,2) + 
                  kK*(kP*kq*(3*G1(0) - 2*G2(0)) + 
                     2*(G1(0) - G2(0))*pow(kK,2) + 
                     (G1(0) - G2(0))*pow(kP,2) + 
                     2*(-G1(0) + G2(0))*pow(kq,2))*pow(qe2,2)))) + 
         Md*G1(qe2)*(Md*(64*kK*me2*G1(0)*pow(kq,2)*pow(Mp,4) - 
               3*qe2*G2(0)*(qe2*
                   (4*kK*kP*(-kP + kq)*me2 + 
                     nu*(2*kP*pow(kK,2) + kq*pow(kP - kq,2))) + 
                  kK*kP*kq*pow(nu,2) + 
                  kK*((kP - kq)*kq + pow(kK,2))*pow(qe2,2)) - 
               4*Mp2*(qe2*(-4*kK*kP*(kP + 2*kq)*me2*G1(0) + 
                     nu*((kP + kq)*(2*G1(0) - G2(0))*pow(kK,2) + 
                        (-2*kP*G1(0) - 2*kq*G1(0) + kP*G2(0))*pow(kq,2))) \
+ kK*G2(0)*pow(kq,2)*pow(nu,2) + 
                  kK*(-(kP*kq*G2(0)) + (G1(0) - 3*G2(0))*pow(kK,2) - 
                     G2(0)*pow(kP,2) - G1(0)*pow(kq,2) + 3*G2(0)*pow(kq,2)\
)*pow(qe2,2))) + Mp*(64*kK*me2*G1(0)*pow(kq,2)*pow(Mp,4) - 
               3*qe2*(2*G1(0) - G2(0))*
                (qe2*(4*kK*kP*(-kP + kq)*me2 + 
                     nu*(2*kP*pow(kK,2) + kq*pow(kP - kq,2))) + 
                  kK*kP*kq*pow(nu,2) + 
                  kK*((kP - kq)*kq + pow(kK,2))*pow(qe2,2)) + 
               4*Mp2*(qe2*(4*kK*kP*(kP + 2*kq)*me2*G1(0) + 
                     nu*(-(kP*G2(0)*pow(kK,2)) - kq*G2(0)*pow(kK,2) + 
                        kP*G2(0)*pow(kq,2) + 2*G1(0)*pow(kq,3))) + 
                  kK*(-2*G1(0) + G2(0))*pow(kq,2)*pow(nu,2) + 
                  kK*(kP*kq*(2*G1(0) - G2(0)) + 
                     (5*G1(0) - 3*G2(0))*pow(kK,2) + 
                     (2*G1(0) - G2(0))*pow(kP,2) + 
                     (-5*G1(0) + 3*G2(0))*pow(kq,2))*pow(qe2,2))))))*
    pow(qp2,-1))/36.
);
}

