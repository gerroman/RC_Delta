#include "constants.h"
#include "proton_formfactors.h"
#include "delta_formfactors.h"
#include <cmath>

double iMdiMe1(double nu, double qe2, double qp2, double kK, double kP, double kq) {
return ((e6*(kP + kq - Md2 + Mp2)*Z3*pow(kK - kq,-1)*pow(kK + kq,-1)*pow(Md,-5)*
    pow(Mp,-1)*pow(qe2,-1)*(2*Mp*F1(qp2)*
       (kK*qe2*G3(qe2)*(Md*(-2*kq*(kP + kq)*
                (-3*kP*kq*G1(0) + (G1(0) - 2*G2(0))*pow(kP,2) + 
                  (8*G1(0) - G2(0))*pow(kq,2)) + 
               pow(kK,2)*(2*kP*kq*(7*G1(0) - 2*G2(0)) + 
                  (8*G1(0) - 7*G2(0))*pow(kP,2) + 
                  2*(8*G1(0) - G2(0))*pow(kq,2))) + 
            Mp*(2*kq*(kP + kq)*
                (-3*kP*kq*G2(0) + 4*(G1(0) - G2(0))*pow(kP,2) + 
                  (12*G1(0) - 13*G2(0))*pow(kq,2)) + 
               pow(kK,2)*(kP*kq*(-32*G1(0) + 38*G2(0)) + 
                  (-8*G1(0) + 11*G2(0))*pow(kP,2) + 
                  2*(-12*G1(0) + 13*G2(0))*pow(kq,2)))) + 
         Md*G1(qe2)*(4*kK*kP*Md*Mp*
             (kq*(kP + kq)*(-4*kq*G1(0) + kP*G2(0)) + 
               (4*kP*G1(0) + 4*kq*G1(0) - 2*kP*G2(0) - 4*kq*G2(0))*
                pow(kK,2)) - 2*kP*kq*qe2*G1(0)*pow(kK,3) - 
            5*kP*kq*qe2*G2(0)*pow(kK,3) + 
            4*kq*nu*G1(0)*pow(kK,2)*pow(kP,2) - 
            8*kq*nu*G2(0)*pow(kK,2)*pow(kP,2) - 
            6*qe2*G1(0)*pow(kK,3)*pow(kP,2) + 
            9*qe2*G2(0)*pow(kK,3)*pow(kP,2) + 
            16*kK*kq*me2*G1(0)*pow(kP,3) + 10*kK*kq*qe2*G1(0)*pow(kP,3) + 
            16*kK*kq*me2*G2(0)*pow(kP,3) - 11*kK*kq*qe2*G2(0)*pow(kP,3) - 
            8*kK*me2*G1(0)*pow(kP,4) - 2*kK*qe2*G1(0)*pow(kP,4) + 
            4*kK*me2*G2(0)*pow(kP,4) + kK*qe2*G2(0)*pow(kP,4) - 
            22*kP*nu*G1(0)*pow(kK,2)*pow(kq,2) + 
            11*kP*nu*G2(0)*pow(kK,2)*pow(kq,2) - 
            32*qe2*G1(0)*pow(kK,3)*pow(kq,2) + 
            4*qe2*G2(0)*pow(kK,3)*pow(kq,2) + 
            120*kK*me2*G1(0)*pow(kP,2)*pow(kq,2) + 
            22*kK*qe2*G1(0)*pow(kP,2)*pow(kq,2) - 
            36*kK*me2*G2(0)*pow(kP,2)*pow(kq,2) - 
            29*kK*qe2*G2(0)*pow(kP,2)*pow(kq,2) - 
            2*nu*G1(0)*pow(kP,3)*pow(kq,2) + 
            nu*G2(0)*pow(kP,3)*pow(kq,2) - 
            4*kK*kP*Mp2*((kP*G1(0) + 5*kq*G1(0) - 4*kP*G2(0) - 
                  3*kq*G2(0))*pow(kK,2) + 
               kq*(kP*kq*(3*G1(0) + 2*G2(0)) + 
                  (-2*G1(0) + G2(0))*pow(kP,2) + 
                  (5*G1(0) - 2*G2(0))*pow(kq,2))) + 
            16*kK*kP*me2*G1(0)*pow(kq,3) + 42*kK*kP*qe2*G1(0)*pow(kq,3) - 
            56*kK*kP*me2*G2(0)*pow(kq,3) - 21*kK*kP*qe2*G2(0)*pow(kq,3) + 
            10*nu*G1(0)*pow(kK,2)*pow(kq,3) + 
            nu*G2(0)*pow(kK,2)*pow(kq,3) + 
            6*nu*G1(0)*pow(kP,2)*pow(kq,3) + 
            3*nu*G2(0)*pow(kP,2)*pow(kq,3) - 80*kK*me2*G1(0)*pow(kq,4) + 
            18*kP*nu*G1(0)*pow(kq,4) + 32*kK*qe2*G1(0)*pow(kq,4) - 
            8*kK*me2*G2(0)*pow(kq,4) - 9*kP*nu*G2(0)*pow(kq,4) - 
            4*kK*qe2*G2(0)*pow(kq,4) + 10*nu*G1(0)*pow(kq,5) - 
            11*nu*G2(0)*pow(kq,5)) + 
         Md*G2(qe2)*(-2*kK*kP*Md*Mp*
             (kq*(kP + kq)*(-2*kP*G1(0) + 2*kP*G2(0) + kq*G2(0)) + 
               (4*kP*G1(0) + 8*kq*G1(0) - 5*kP*G2(0) - 9*kq*G2(0))*
                pow(kK,2)) + 2*kK*kP*Mp2*
             (-(kq*(kP + kq)*(2*kP*G1(0) + 2*kq*G1(0) - 2*kP*G2(0) - 
                    kq*G2(0))) + 
               (8*kP*G1(0) + 12*kq*G1(0) - 9*kP*G2(0) - 7*kq*G2(0))*
                pow(kK,2)) - 11*kP*kq*qe2*G1(0)*pow(kK,3) + 
            7*kP*kq*qe2*G2(0)*pow(kK,3) - 
            8*kq*nu*G1(0)*pow(kK,2)*pow(kP,2) + 
            4*kq*nu*G2(0)*pow(kK,2)*pow(kP,2) - 
            qe2*G1(0)*pow(kK,3)*pow(kP,2) - 
            4*qe2*G2(0)*pow(kK,3)*pow(kP,2) + 
            16*kK*kq*me2*G1(0)*pow(kP,3) - 9*kK*kq*qe2*G1(0)*pow(kP,3) - 
            8*kK*kq*me2*G2(0)*pow(kP,3) + 12*kK*kq*qe2*G2(0)*pow(kP,3) + 
            4*kK*me2*G1(0)*pow(kP,4) + kK*qe2*G1(0)*pow(kP,4) - 
            8*kK*me2*G2(0)*pow(kP,4) - 2*kK*qe2*G2(0)*pow(kP,4) - 
            kP*nu*G1(0)*pow(kK,2)*pow(kq,2) - 
            4*kP*nu*G2(0)*pow(kK,2)*pow(kq,2) + 
            8*qe2*G1(0)*pow(kK,3)*pow(kq,2) + 
            2*qe2*G2(0)*pow(kK,3)*pow(kq,2) - 
            36*kK*me2*G1(0)*pow(kP,2)*pow(kq,2) - 
            5*kK*qe2*G1(0)*pow(kP,2)*pow(kq,2) + 
            24*kK*me2*G2(0)*pow(kP,2)*pow(kq,2) + 
            13*kK*qe2*G2(0)*pow(kP,2)*pow(kq,2) + 
            nu*G1(0)*pow(kP,3)*pow(kq,2) - 
            2*nu*G2(0)*pow(kP,3)*pow(kq,2) - 
            56*kK*kP*me2*G1(0)*pow(kq,3) - 3*kK*kP*qe2*G1(0)*pow(kq,3) + 
            40*kK*kP*me2*G2(0)*pow(kq,3) - 3*kK*kP*qe2*G2(0)*pow(kq,3) - 
            11*nu*G1(0)*pow(kK,2)*pow(kq,3) + 
            nu*G2(0)*pow(kK,2)*pow(kq,3) + 
            3*nu*G1(0)*pow(kP,2)*pow(kq,3) - 
            3*nu*G2(0)*pow(kP,2)*pow(kq,3) - 8*kK*me2*G1(0)*pow(kq,4) + 
            3*kP*nu*G1(0)*pow(kq,4) - 8*kK*qe2*G1(0)*pow(kq,4) + 
            16*kK*me2*G2(0)*pow(kq,4) - 2*kK*qe2*G2(0)*pow(kq,4) + 
            nu*G1(0)*pow(kq,5) + nu*G2(0)*pow(kq,5))) + 
      kK*(2*kP*Md*Mp*F1(qp2)*(G1(qe2)*
             (-(kq*(kP + kq)*
                  (2*kP*kq*(4*G1(0) + G2(0)) + 
                    (-2*G1(0) + G2(0))*pow(kP,2) + 
                    (10*G1(0) - 11*G2(0))*pow(kq,2))) + 
               pow(kK,2)*(kP*kq*(2*G1(0) - G2(0)) + 6*G2(0)*pow(kP,2) - 
                  (10*G1(0) + G2(0))*pow(kq,2))) + 
            G2(qe2)*(pow(kK,2)*
                (kP*kq*(11*G1(0) - 4*G2(0)) + 
                  6*(G1(0) - G2(0))*pow(kP,2) + 
                  (11*G1(0) - G2(0))*pow(kq,2)) - 
               kq*(kP*(G1(0) - 2*G2(0)) + kq*(G1(0) + G2(0)))*
                pow(kP + kq,2))) + 
         F2(qp2)*(-((kP + kq)*qe2*G3(qe2)*
               (2*pow(kq,2)*(-3*kP*kq*G2(0) + 
                    2*(G1(0) - 2*G2(0))*pow(kP,2) + 
                    (10*G1(0) - 11*G2(0))*pow(kq,2)) + 
                 pow(kK,2)*(-4*kP*kq*(G1(0) - 2*G2(0)) + 
                    (6*G1(0) - 3*G2(0))*pow(kP,2) + 
                    2*(-10*G1(0) + 11*G2(0))*pow(kq,2)))) + 
            2*kP*Md*G1(qe2)*(2*Mp*
                (-((kP + kq)*
                     (4*kP*G1(0) + 4*kq*G1(0) - kP*G2(0) - 7*kq*G2(0))*
                     pow(kq,2)) + 
                  pow(kK,2)*(kP*kq*(2*G1(0) - 5*G2(0)) + 
                     G2(0)*pow(kP,2) - (4*G1(0) + 3*G2(0))*pow(kq,2))) + 
               Md*(-((kP + kq)*
                     (kq*(-8*G1(0) + G2(0)) + kP*(4*G1(0) + G2(0)))*
                     pow(kq,2)) + 
                  pow(kK,2)*(2*kP*kq*(-2*G1(0) + G2(0)) + 
                     (4*G1(0) - 2*G2(0))*pow(kP,2) + 
                     (-8*G1(0) + 7*G2(0))*pow(kq,2)))) - 
            2*kP*Md*G2(qe2)*(-(Mp*
                  (pow(kK,2)*(2*kP*kq*(G1(0) + 3*G2(0)) + 
                       2*(G1(0) - G2(0))*pow(kP,2) + 
                       (6*G1(0) + 5*G2(0))*pow(kq,2)) + 
                    (2*G1(0) - 3*G2(0))*pow(kq,2)*pow(kP + kq,2))) + 
               Md*(pow(kK,2)*(-2*kP*kq*(G1(0) - 2*G2(0)) + 
                     (2*G1(0) - G2(0))*pow(kP,2) + 
                     (-7*G1(0) + 8*G2(0))*pow(kq,2)) + 
                  (G1(0) - 2*G2(0))*pow(kq,2)*pow(kP + kq,2))))) + 
      F2(qp2)*(qe2*G3(qe2)*(2*kK*Md*Mp*
             (-2*(kP + kq)*(-6*kP*G1(0) + 5*kq*G1(0) + 2*kP*G2(0) + 
                  3*kq*G2(0))*pow(kq,2) + 
               pow(kK,2)*(4*kP*kq*(G1(0) + 2*G2(0)) + 
                  3*(G1(0) - G2(0))*pow(kP,2) + 
                  2*(5*G1(0) + 3*G2(0))*pow(kq,2))) - 
            2*kK*Mp2*(2*pow(kq,2)*
                (kP*kq*(6*G1(0) - 8*G2(0)) + 
                  (G1(0) - 2*G2(0))*pow(kP,2) + 
                  (11*G1(0) - 12*G2(0))*pow(kq,2)) + 
               pow(kK,2)*(-2*kP*kq*(5*G1(0) - 6*G2(0)) + 
                  (7*G1(0) - 5*G2(0))*pow(kP,2) + 
                  2*(-11*G1(0) + 12*G2(0))*pow(kq,2))) - 
            (kP + kq)*(-16*kK*kP*kq*(kP + kq)*me2*(G1(0) - 2*G2(0)) + 
               nu*(pow(kK,2)*
                   (4*kP*kq*(G1(0) - 2*G2(0)) + 
                     (-2*G1(0) + G2(0))*pow(kP,2) + 
                     2*(4*G1(0) - 5*G2(0))*pow(kq,2)) - 
                  2*pow(kq,2)*
                   (3*kP*kq*(4*G1(0) - 5*G2(0)) + 
                     (-2*G1(0) + G2(0))*pow(kP,2) + 
                     2*(G1(0) - 2*G2(0))*pow(kq,2))) + 
               kK*qe2*((2*kP*G1(0) + 14*kq*G1(0) - 4*kP*G2(0) - 
                     19*kq*G2(0))*pow(kK,2) + 
                  kq*(3*kP*kq*(4*G1(0) - 3*G2(0)) - 
                     18*(G1(0) - G2(0))*pow(kP,2) + 
                     (-14*G1(0) + 19*G2(0))*pow(kq,2))))) + 
         Md*G1(qe2)*(Md*(4*kK*kP*Mp2*
                ((2*kP*G1(0) - 4*kq*G1(0) - kP*G2(0) + 3*kq*G2(0))*
                   pow(kK,2) + 
                  (-2*kP*G1(0) + 4*kq*G1(0) - kP*G2(0) - kq*G2(0))*
                   pow(kq,2)) + 
               kK*qe2*(pow(kK,2)*
                   (-3*kP*kq*G2(0) + (8*G1(0) - 4*G2(0))*pow(kP,2) - 
                     8*(G1(0) - 2*G2(0))*pow(kq,2)) + 
                  kq*(kP + kq)*
                   (kP*kq*(8*G1(0) - 19*G2(0)) + 
                     (-4*G1(0) + 11*G2(0))*pow(kP,2) + 
                     8*(G1(0) - 2*G2(0))*pow(kq,2))) + 
               2*kq*(4*kK*kP*(kP + kq)*me2*
                   (kq*(-8*G1(0) + G2(0)) + kP*(4*G1(0) + G2(0))) + 
                  nu*((kP + kq)*
                      (kq*(-8*G1(0) + G2(0)) + kP*(4*G1(0) + G2(0)))*
                      pow(kq,2) + 
                     pow(kK,2)*
                      (kP*kq*(4*G1(0) - 2*G2(0)) + 
                        (-4*G1(0) + 2*G2(0))*pow(kP,2) + 
                        (8*G1(0) - 7*G2(0))*pow(kq,2))))) + 
            Mp*(kK*qe2*(kq*(kP + kq)*
                   (-(kP*kq*(6*G1(0) + 17*G2(0))) + 
                     (22*G1(0) - 13*G2(0))*pow(kP,2) + 
                     4*(12*G1(0) + G2(0))*pow(kq,2)) - 
                  pow(kK,2)*(kP*kq*(-10*G1(0) + 17*G2(0)) + 
                     4*(2*G1(0) - 5*G2(0))*pow(kP,2) + 
                     4*(12*G1(0) + G2(0))*pow(kq,2))) + 
               4*(kK*kP*Mp2*((-4*kq*G1(0) + kP*G2(0) - 3*kq*G2(0))*
                      pow(kK,2) + 
                     (-4*kP*G1(0) - 4*kq*G1(0) + kP*G2(0) + 7*kq*G2(0))*
                      pow(kq,2)) + 
                  kq*(2*kK*(kP + kq)*me2*
                      (kP*kq*(22*G1(0) - 15*G2(0)) + 
                        (6*G1(0) - G2(0))*pow(kP,2) - 
                        2*(10*G1(0) + G2(0))*pow(kq,2)) + 
                     nu*((kP + kq)*
                         (4*kP*G1(0) + 4*kq*G1(0) - kP*G2(0) - 
                         7*kq*G2(0))*pow(kq,2) + 
                        pow(kK,2)*
                         (-2*kP*kq*(6*G1(0) - 5*G2(0)) + 
                         2*(G1(0) - G2(0))*pow(kP,2) + 
                         (4*G1(0) + 3*G2(0))*pow(kq,2))))))) - 
         Md*G2(qe2)*(Mp*(-4*kK*kP*Mp2*
                ((kP*(G1(0) - G2(0)) + 3*kq*(G1(0) + G2(0)))*pow(kK,2) + 
                  (kP + kq)*(G1(0) - G2(0))*pow(kq,2)) + 
               kK*qe2*(pow(kK,2)*
                   (kP*kq*(21*G1(0) - 16*G2(0)) + 12*G2(0)*pow(kP,2) - 
                     12*(G1(0) + G2(0))*pow(kq,2)) + 
                  kq*(kP + kq)*
                   (-3*kP*kq*(9*G1(0) - 8*G2(0)) + 
                     (13*G1(0) - 18*G2(0))*pow(kP,2) + 
                     12*(G1(0) + G2(0))*pow(kq,2))) + 
               2*kq*(4*kK*(kP + kq)*me2*
                   (kP*kq*(15*G1(0) - 8*G2(0)) + 
                     (G1(0) - 4*G2(0))*pow(kP,2) + 
                     2*(G1(0) - 2*G2(0))*pow(kq,2)) + 
                  nu*(pow(kK,2)*
                      (kP*kq*(-8*G1(0) + 14*G2(0)) + 4*G1(0)*pow(kP,2) + 
                        (6*G1(0) + 5*G2(0))*pow(kq,2)) + 
                     (2*G1(0) - 3*G2(0))*pow(kq,2)*pow(kP + kq,2)))) + 
            Md*(4*kK*kP*Mp2*(G1(0) - G2(0))*
                ((kP - 3*kq)*pow(kK,2) + (kP + kq)*pow(kq,2)) + 
               kK*qe2*(kq*(kP + kq)*
                   (3*kP*kq*(5*G1(0) - 6*G2(0)) + 
                     (-11*G1(0) + 10*G2(0))*pow(kP,2) + 
                     12*(G1(0) - G2(0))*pow(kq,2)) + 
                  pow(kK,2)*(kP*kq*(7*G1(0) - 8*G2(0)) + 
                     6*(G1(0) - G2(0))*pow(kP,2) + 
                     12*(-G1(0) + G2(0))*pow(kq,2))) - 
               2*kq*(4*kK*kP*me2*(G1(0) - 2*G2(0))*pow(kP + kq,2) + 
                  nu*(pow(kK,2)*
                      (-2*kP*kq*(G1(0) - 2*G2(0)) + 
                        (2*G1(0) - G2(0))*pow(kP,2) + 
                        (-7*G1(0) + 8*G2(0))*pow(kq,2)) + 
                     (G1(0) - 2*G2(0))*pow(kq,2)*pow(kP + kq,2)))))) + 
      2*Mp*F1(qp2)*(qe2*G3(qe2)*
          (Mp*(8*kP*kq*nu*G1(0)*pow(kK,2) - 
               10*kP*kq*nu*G2(0)*pow(kK,2) - 
               16*kK*kq*me2*G1(0)*pow(kP,2) + 
               32*kK*kq*me2*G2(0)*pow(kP,2) + 
               8*nu*G1(0)*pow(kK,2)*pow(kP,2) - 
               9*nu*G2(0)*pow(kK,2)*pow(kP,2) - 
               32*kK*me2*G1(0)*pow(kP,3) + 32*kK*me2*G2(0)*pow(kP,3) + 
               kK*qe2*((12*kP*G1(0) + 12*kq*G1(0) - 13*kP*G2(0) - 
                     11*kq*G2(0))*pow(kK,2) - 
                  (8*kP*G1(0) - 12*kq*G1(0) - 8*kP*G2(0) + 11*kq*G2(0))*
                   (pow(kP,2) - pow(kq,2))) + 
               64*kK*kP*me2*G1(0)*pow(kq,2) - 
               16*kK*kP*me2*G2(0)*pow(kq,2) - 
               4*nu*G2(0)*pow(kK,2)*pow(kq,2) - 
               24*nu*G1(0)*pow(kP,2)*pow(kq,2) + 
               26*nu*G2(0)*pow(kP,2)*pow(kq,2) + 
               4*kK*Mp2*((-10*kP*G1(0) - 14*kq*G1(0) + 11*kP*G2(0) + 
                     15*kq*G2(0))*pow(kK,2) + 
                  kq*(2*kP*kq*(4*G1(0) - 5*G2(0)) + 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     (14*G1(0) - 15*G2(0))*pow(kq,2))) + 
               48*kK*me2*G1(0)*pow(kq,3) - 16*kP*nu*G1(0)*pow(kq,3) - 
               16*kK*me2*G2(0)*pow(kq,3) + 24*kP*nu*G2(0)*pow(kq,3) + 
               8*nu*G1(0)*pow(kq,4) - 2*nu*G2(0)*pow(kq,4)) + 
            Md*(-12*kP*kq*nu*G1(0)*pow(kK,2) + 
               24*kK*kq*me2*G1(0)*pow(kP,2) - 
               4*nu*G1(0)*pow(kK,2)*pow(kP,2) + 
               5*nu*G2(0)*pow(kK,2)*pow(kP,2) + 
               8*kK*me2*G1(0)*pow(kP,3) - 16*kK*me2*G2(0)*pow(kP,3) - 
               64*kK*kP*me2*G1(0)*pow(kq,2) + 
               56*kK*kP*me2*G2(0)*pow(kq,2) + 
               10*nu*G1(0)*pow(kK,2)*pow(kq,2) - 
               14*nu*G2(0)*pow(kK,2)*pow(kq,2) + 
               10*nu*G1(0)*pow(kP,2)*pow(kq,2) - 
               14*nu*G2(0)*pow(kP,2)*pow(kq,2) + 
               kK*qe2*((-28*kq*G1(0) + 11*kq*G2(0) + 
                     kP*(-8*G1(0) + G2(0)))*pow(kK,2) + 
                  (kP + kq)*(-3*kP*kq*(8*G1(0) - 5*G2(0)) + 
                     2*(G1(0) - 2*G2(0))*pow(kP,2) + 
                     (28*G1(0) - 11*G2(0))*pow(kq,2))) + 
               4*kK*Mp2*((6*kP*G1(0) + 10*kq*G1(0) - 7*kP*G2(0) - 
                     5*kq*G2(0))*pow(kK,2) + 
                  kq*(6*kP*kq*(-G1(0) + G2(0)) - 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     5*(-2*G1(0) + G2(0))*pow(kq,2))) - 
               80*kK*me2*G1(0)*pow(kq,3) - 4*kP*nu*G1(0)*pow(kq,3) + 
               40*kK*me2*G2(0)*pow(kq,3) + 2*kP*nu*G2(0)*pow(kq,3) - 
               14*nu*G1(0)*pow(kq,4) + 16*nu*G2(0)*pow(kq,4))) + 
         Md*G2(qe2)*(-6*kP*kq*nu*qe2*G1(0)*pow(kK,2) + 
            3*kP*kq*nu*qe2*G2(0)*pow(kK,2) + 
            32*kK*kq*me2*qe2*G1(0)*pow(kP,2) - 
            16*kK*kq*me2*qe2*G2(0)*pow(kP,2) + 
            6*nu*qe2*G1(0)*pow(kK,2)*pow(kP,2) - 
            3*nu*qe2*G2(0)*pow(kK,2)*pow(kP,2) + 
            4*kq*me2*nu*G1(0)*pow(kP,3) - 8*kq*me2*nu*G2(0)*pow(kP,3) - 
            12*kK*me2*qe2*G2(0)*pow(kP,3) - 
            12*kK*kP*me2*qe2*G1(0)*pow(kq,2) + 
            12*kK*kP*me2*qe2*G2(0)*pow(kq,2) + 
            24*nu*qe2*G1(0)*pow(kK,2)*pow(kq,2) - 
            12*nu*qe2*G2(0)*pow(kK,2)*pow(kq,2) + 
            4*nu*qe2*G1(0)*pow(kP,2)*pow(kq,2) - 
            24*me2*nu*G2(0)*pow(kP,2)*pow(kq,2) - 
            5*nu*qe2*G2(0)*pow(kP,2)*pow(kq,2) + 
            36*kP*me2*nu*G1(0)*pow(kq,3) - 
            44*kK*me2*qe2*G1(0)*pow(kq,3) - 
            16*kP*nu*qe2*G1(0)*pow(kq,3) - 24*kP*me2*nu*G2(0)*pow(kq,3) + 
            16*kK*me2*qe2*G2(0)*pow(kq,3) + 5*kP*nu*qe2*G2(0)*pow(kq,3) - 
            2*Mp2*(kq*nu*(-(kq*(kP + kq)*
                     (2*kP*G1(0) + 2*kq*G1(0) - 2*kP*G2(0) - kq*G2(0))) \
+ (6*kP*G1(0) + 12*kq*G1(0) - 7*kP*G2(0) - 7*kq*G2(0))*pow(kK,2)) - 
               4*kK*me2*(kq*(3*G1(0) - 2*G2(0))*pow(kP,2) + 
                  2*(G1(0) - G2(0))*pow(kP,3) - 3*kP*G1(0)*pow(kq,2) + 
                  2*G1(0)*pow(kq,3)) + 
               kK*qe2*((6*kq*(-2*G1(0) + G2(0)) + 
                     kP*(3*G1(0) + 4*G2(0)))*pow(kK,2) + 
                  3*kq*(5*G1(0) - 6*G2(0))*pow(kP,2) - 
                  2*(G1(0) - G2(0))*pow(kP,3) + 
                  kP*(17*G1(0) - 19*G2(0))*pow(kq,2) + 
                  6*(2*G1(0) - G2(0))*pow(kq,3))) + 
            40*me2*nu*G1(0)*pow(kq,4) - 20*nu*qe2*G1(0)*pow(kq,4) - 
            8*me2*nu*G2(0)*pow(kq,4) + 10*nu*qe2*G2(0)*pow(kq,4) - 
            2*Md*Mp*(-2*kP*kq*nu*G1(0)*pow(kK,2) + 
               3*kP*kq*nu*G2(0)*pow(kK,2) + 
               8*kP*Mp2*(G1(0) - G2(0))*pow(kK,3) + 
               4*kK*kq*me2*G1(0)*pow(kP,2) - 
               8*kK*kq*me2*G2(0)*pow(kP,2) + 8*kK*me2*G1(0)*pow(kP,3) - 
               8*kK*me2*G2(0)*pow(kP,3) - 12*kK*kP*me2*G1(0)*pow(kq,2) - 
               8*nu*G1(0)*pow(kK,2)*pow(kq,2) + 
               9*nu*G2(0)*pow(kK,2)*pow(kq,2) + 
               2*nu*G1(0)*pow(kP,2)*pow(kq,2) - 
               2*nu*G2(0)*pow(kP,2)*pow(kq,2) - 
               8*kK*me2*G1(0)*pow(kq,3) + 2*kP*nu*G1(0)*pow(kq,3) - 
               3*kP*nu*G2(0)*pow(kq,3) - 
               kK*qe2*((kP*G1(0) - 6*kq*G1(0) + 8*kq*G2(0))*pow(kK,2) + 
                  kq*(9*G1(0) - 10*G2(0))*pow(kP,2) - 
                  2*(G1(0) - G2(0))*pow(kP,3) + 
                  kP*(15*G1(0) - 19*G2(0))*pow(kq,2) + 
                  2*(3*G1(0) - 4*G2(0))*pow(kq,3)) - nu*G2(0)*pow(kq,4)) \
+ 16*kP*(G1(0) - G2(0))*pow(kK,3)*pow(Mp,4) - 
            kK*kq*G1(0)*pow(kP,2)*pow(nu,2) + 
            2*kK*kq*G2(0)*pow(kP,2)*pow(nu,2) + 
            7*kK*kP*G1(0)*pow(kq,2)*pow(nu,2) + 
            kK*kP*G2(0)*pow(kq,2)*pow(nu,2) - 
            10*kK*G1(0)*pow(kq,3)*pow(nu,2) + 
            8*kK*G2(0)*pow(kq,3)*pow(nu,2) + 
            2*kP*G1(0)*pow(kK,3)*pow(qe2,2) - 
            16*kq*G1(0)*pow(kK,3)*pow(qe2,2) - 
            4*kP*G2(0)*pow(kK,3)*pow(qe2,2) + 
            5*kq*G2(0)*pow(kK,3)*pow(qe2,2) - 
            7*kK*kq*G1(0)*pow(kP,2)*pow(qe2,2) - 
            kK*kq*G2(0)*pow(kP,2)*pow(qe2,2) + 
            6*kK*G1(0)*pow(kP,3)*pow(qe2,2) - 
            9*kK*G2(0)*pow(kP,3)*pow(qe2,2) + 
            3*kK*kP*G1(0)*pow(kq,2)*pow(qe2,2) + 
            3*kK*kP*G2(0)*pow(kq,2)*pow(qe2,2) + 
            16*kK*G1(0)*pow(kq,3)*pow(qe2,2) - 
            5*kK*G2(0)*pow(kq,3)*pow(qe2,2)) + 
         Md*G1(qe2)*(40*kP*kq*nu*qe2*G1(0)*pow(kK,2) - 
            14*kP*kq*nu*qe2*G2(0)*pow(kK,2) - 
            112*kK*kq*me2*qe2*G1(0)*pow(kP,2) + 
            56*kK*kq*me2*qe2*G2(0)*pow(kP,2) + 
            4*nu*qe2*G1(0)*pow(kK,2)*pow(kP,2) + 
            4*nu*qe2*G2(0)*pow(kK,2)*pow(kP,2) - 
            8*kq*me2*nu*G1(0)*pow(kP,3) - 16*kK*me2*qe2*G1(0)*pow(kP,3) + 
            4*kq*me2*nu*G2(0)*pow(kP,3) + 8*kK*me2*qe2*G2(0)*pow(kP,3) + 
            40*kK*kP*me2*qe2*G1(0)*pow(kq,2) + 
            4*kK*kP*me2*qe2*G2(0)*pow(kq,2) - 
            36*nu*qe2*G1(0)*pow(kK,2)*pow(kq,2) + 
            18*nu*qe2*G2(0)*pow(kK,2)*pow(kq,2) + 
            48*me2*nu*G1(0)*pow(kP,2)*pow(kq,2) - 
            20*nu*qe2*G1(0)*pow(kP,2)*pow(kq,2) + 
            10*nu*qe2*G2(0)*pow(kP,2)*pow(kq,2) + 
            4*Md*Mp*(kq*nu*(kq*(kP + kq)*(4*kq*G1(0) - kP*G2(0)) + 
                  (-4*kP*G1(0) - 4*kq*G1(0) + kP*G2(0) + 4*kq*G2(0))*
                   pow(kK,2)) + 
               4*kK*kP*Mp2*((G1(0) - G2(0))*pow(kK,2) - 
                  G1(0)*pow(kq,2)) - 
               kK*(qe2*(kq*(4*G1(0) - 5*G2(0))*pow(kP,2) + 
                     G2(0)*pow(kP,3) + 
                     4*kq*G2(0)*(pow(kK,2) - pow(kq,2)) + 
                     2*kP*(2*G1(0) - 5*G2(0))*pow(kq,2)) + 
                  2*(kP + kq)*me2*
                   (-(kP*kq*(8*G1(0) + G2(0))) + 2*G2(0)*pow(kP,2) + 
                     2*(4*G1(0) - G2(0))*pow(kq,2)))) - 
            24*kP*me2*nu*G1(0)*pow(kq,3) + 
            136*kK*me2*qe2*G1(0)*pow(kq,3) - 4*kP*nu*qe2*G1(0)*pow(kq,3) + 
            36*kP*me2*nu*G2(0)*pow(kq,3) - 44*kK*me2*qe2*G2(0)*pow(kq,3) + 
            8*kP*nu*qe2*G2(0)*pow(kq,3) - 
            4*Mp2*(kq*nu*((-5*kq*G1(0) + 3*kq*G2(0) + 
                     kP*(G1(0) + 3*G2(0)))*pow(kK,2) - 
                  kq*(kP*kq*(3*G1(0) + 2*G2(0)) + 
                     (-2*G1(0) + G2(0))*pow(kP,2) + 
                     (5*G1(0) - 2*G2(0))*pow(kq,2))) + 
               2*kK*me2*(-(kq*(8*G1(0) + 3*G2(0))*pow(kP,2)) + 
                  (4*G1(0) - 2*G2(0))*pow(kP,3) + 
                  3*kP*(-2*G1(0) + G2(0))*pow(kq,2) + 
                  2*(12*G1(0) - G2(0))*pow(kq,3)) + 
               kK*qe2*((3*kP*G1(0) + 20*kq*G1(0) - 3*kP*G2(0) - 
                     7*kq*G2(0))*pow(kK,2) + 
                  kq*(-9*G1(0) + 8*G2(0))*pow(kP,2) + 
                  (2*G1(0) - G2(0))*pow(kP,3) + 
                  2*kP*(-11*G1(0) + 8*G2(0))*pow(kq,2) + 
                  (-20*G1(0) + 7*G2(0))*pow(kq,3))) - 
            80*me2*nu*G1(0)*pow(kq,4) + 16*nu*qe2*G1(0)*pow(kq,4) + 
            40*me2*nu*G2(0)*pow(kq,4) - 2*nu*qe2*G2(0)*pow(kq,4) - 
            16*kK*kP*((G1(0) - G2(0))*pow(kK,2) + G1(0)*pow(kq,2))*
             pow(Mp,4) + 2*kK*kq*G1(0)*pow(kP,2)*pow(nu,2) - 
            kK*kq*G2(0)*pow(kP,2)*pow(nu,2) - 
            14*kK*kP*G1(0)*pow(kq,2)*pow(nu,2) + 
            7*kK*kP*G2(0)*pow(kq,2)*pow(nu,2) + 
            20*kK*G1(0)*pow(kq,3)*pow(nu,2) - 
            10*kK*G2(0)*pow(kq,3)*pow(nu,2) + 
            4*kP*G1(0)*pow(kK,3)*pow(qe2,2) + 
            40*kq*G1(0)*pow(kK,3)*pow(qe2,2) + 
            4*kP*G2(0)*pow(kK,3)*pow(qe2,2) - 
            14*kq*G2(0)*pow(kK,3)*pow(qe2,2) + 
            6*kK*kq*G1(0)*pow(kP,2)*pow(qe2,2) + 
            9*kK*kq*G2(0)*pow(kP,2)*pow(qe2,2) - 
            4*kK*G1(0)*pow(kP,3)*pow(qe2,2) + 
            8*kK*G2(0)*pow(kP,3)*pow(qe2,2) - 
            30*kK*kP*G1(0)*pow(kq,2)*pow(qe2,2) + 
            15*kK*kP*G2(0)*pow(kq,2)*pow(qe2,2) - 
            40*kK*G1(0)*pow(kq,3)*pow(qe2,2) + 
            14*kK*G2(0)*pow(kq,3)*pow(qe2,2))) + 
      F2(qp2)*(qe2*G3(qe2)*(-2*Mp2*
             (-8*kK*kq*(kP + kq)*me2*
                (5*kP*G1(0) + 6*kq*G1(0) - 6*kP*G2(0) - 2*kq*G2(0)) + 
               nu*(2*(G1(0) - G2(0))*
                   (-12*kP*kq + pow(kP,2) - 7*pow(kq,2))*pow(kq,2) + 
                  pow(kK,2)*(2*kP*kq*(5*G1(0) - 6*G2(0)) + 
                     (-G1(0) + G2(0))*pow(kP,2) + 
                     2*(6*G1(0) - 5*G2(0))*pow(kq,2))) + 
               kK*qe2*((kP*G1(0) + 7*kq*G1(0) - 5*kP*G2(0) - 
                     15*kq*G2(0))*pow(kK,2) + 
                  kq*(-2*kP*kq*(G1(0) - 5*G2(0)) - 
                     17*(G1(0) - G2(0))*pow(kP,2) + 
                     (-7*G1(0) + 15*G2(0))*pow(kq,2)))) + 
            2*Md*Mp*(-4*kP*kq*nu*G1(0)*pow(kK,2) - 
               8*kP*kq*nu*G2(0)*pow(kK,2) - 
               8*kK*kq*me2*G1(0)*pow(kP,2) + 
               32*kK*kq*me2*G2(0)*pow(kP,2) - 
               nu*G1(0)*pow(kK,2)*pow(kP,2) + 
               nu*G2(0)*pow(kK,2)*pow(kP,2) - 
               88*kK*kP*me2*G1(0)*pow(kq,2) + 
               72*kK*kP*me2*G2(0)*pow(kq,2) + 
               14*nu*G1(0)*pow(kK,2)*pow(kq,2) - 
               18*nu*G2(0)*pow(kK,2)*pow(kq,2) + 
               2*nu*G1(0)*pow(kP,2)*pow(kq,2) - 
               2*nu*G2(0)*pow(kP,2)*pow(kq,2) + 
               kK*qe2*(-((3*kP*G1(0) + 19*kq*G1(0) + 5*kP*G2(0) - 
                       3*kq*G2(0))*pow(kK,2)) + 
                  kq*(6*kP*kq*G2(0) + 
                     (-23*G1(0) + 13*G2(0))*pow(kP,2) + 
                     (19*G1(0) - 3*G2(0))*pow(kq,2))) - 
               80*kK*me2*G1(0)*pow(kq,3) - 18*kP*nu*G1(0)*pow(kq,3) + 
               40*kK*me2*G2(0)*pow(kq,3) + 18*kP*nu*G2(0)*pow(kq,3) + 
               4*kK*Mp2*((kP*(G1(0) - G2(0)) + 3*kq*(G1(0) + G2(0)))*
                   pow(kK,2) - 3*(G1(0) + G2(0))*pow(kq,3)) - 
               20*nu*G1(0)*pow(kq,4) + 20*nu*G2(0)*pow(kq,4)) - 
            8*kK*(G1(0) - G2(0))*((kP - 3*kq)*pow(kK,2) + 3*pow(kq,3))*
             pow(Mp,4) + (kP + kq)*
             (nu*(16*(kP + kq)*me2*(G1(0) - 2*G2(0)) + 
                  kK*nu*(-10*G1(0) + 11*G2(0)))*pow(kq,2) - 
               qe2*(8*kK*kP*(kP + kq)*me2*(G1(0) - 2*G2(0)) + 
                  nu*((-8*kq*G1(0) + 13*kq*G2(0) + 
                        kP*(4*G1(0) + G2(0)))*pow(kK,2) + 
                     kq*(kP*kq*(16*G1(0) - 23*G2(0)) + 
                        (-2*G1(0) + G2(0))*pow(kP,2) - 
                        2*(G1(0) + G2(0))*pow(kq,2)))) + 
               kK*(2*kP*kq*(7*G1(0) - 8*G2(0)) + 
                  2*(G1(0) - 2*G2(0))*pow(kK,2) + 
                  (-8*G1(0) + 7*G2(0))*pow(kP,2) - 
                  2*(G1(0) - 2*G2(0))*pow(kq,2))*pow(qe2,2))) + 
         Md*G2(qe2)*(Md*(kK*qe2*
                (-4*kP*me2*(G1(0) - 2*G2(0))*(pow(kP,2) - pow(kq,2)) + 
                  qe2*(kq*(17*G1(0) - 22*G2(0))*pow(kP,2) + 
                     (-5*G1(0) + 4*G2(0))*pow(kP,3) + 
                     6*kP*(5*G1(0) - 6*G2(0))*pow(kq,2) + 
                     2*kq*(4*G1(0) - 5*G2(0))*(-pow(kK,2) + pow(kq,2)))) \
- 4*Mp2*(kq*(-8*kK*(kP + kq)*me2*(kq*G1(0) + kP*(G1(0) - G2(0))) + 
                     nu*(G1(0) - G2(0))*
                      ((kP + 3*kq)*pow(kK,2) - (kP + kq)*pow(kq,2))) + 
                  kK*qe2*((kP*(G1(0) - 2*G2(0)) + 
                        2*kq*(-2*G1(0) + G2(0)))*pow(kK,2) + 
                     kq*(kP*kq*(7*G1(0) - 5*G2(0)) - 
                        4*(G1(0) - G2(0))*pow(kP,2) + 
                        2*(2*G1(0) - G2(0))*pow(kq,2)))) + 
               nu*(qe2*((kP + kq)*
                      (-(kP*G1(0)) + 13*kq*G1(0) + 2*kP*G2(0) - 
                        14*kq*G2(0))*pow(kq,2) + 
                     pow(kK,2)*
                      (4*kP*kq*(G1(0) - 2*G2(0)) + 
                        (-4*G1(0) + 2*G2(0))*pow(kP,2) + 
                        3*(3*G1(0) - 4*G2(0))*pow(kq,2))) + 
                  8*me2*(G1(0) - 2*G2(0))*pow(kq,2)*pow(kP + kq,2)) + 
               6*kK*kP*(-G1(0) + G2(0))*pow(kq,2)*pow(nu,2)) + 
            Mp*(2*nu*pow(kq,2)*
                (kK*nu*(8*kP*G1(0) - 10*kq*G1(0) - kP*G2(0) + 
                     8*kq*G2(0)) + 
                  4*me2*G1(0)*(2*kP*kq - 5*pow(kP,2) + 7*pow(kq,2))) + 
               qe2*(4*kK*(kP + kq)*me2*
                   (13*kP*kq*G1(0) + (G1(0) - 4*G2(0))*pow(kP,2) + 
                     2*(-11*G1(0) + 4*G2(0))*pow(kq,2)) + 
                  nu*((kP + kq)*
                      (5*kP*G1(0) - 47*kq*G1(0) + 30*kq*G2(0))*pow(kq,2) \
+ pow(kK,2)*(kP*kq*(-26*G1(0) + 18*G2(0)) + 
                        (6*G1(0) - 2*G2(0))*pow(kP,2) + 
                        (37*G1(0) - 14*G2(0))*pow(kq,2)))) - 
               4*Mp2*(kq*(-8*kK*me2*
                      (kP*kq*(-3*G1(0) + G2(0)) + 
                        (-G1(0) + G2(0))*pow(kP,2) + G1(0)*pow(kq,2)) + 
                     nu*((kP*(-G1(0) + G2(0)) + 3*kq*(G1(0) + G2(0)))*
                        pow(kK,2) + (kP + kq)*(G1(0) - G2(0))*pow(kq,2))\
) + kK*qe2*((-2*kq*(3*G1(0) + G2(0)) + kP*(3*G1(0) + 2*G2(0)))*
                      pow(kK,2) + 
                     kq*(kP*kq*(-3*G1(0) + 5*G2(0)) + 
                        4*(G1(0) - G2(0))*pow(kP,2) + 
                        2*(3*G1(0) + G2(0))*pow(kq,2)))) + 
               kK*(2*(3*kP*(G1(0) - 2*G2(0)) + kq*(-12*G1(0) + G2(0)))*
                   pow(kK,2) + 
                  (kP + kq)*(kP*kq*(-46*G1(0) + 40*G2(0)) + 
                     (5*G1(0) - 8*G2(0))*pow(kP,2) + 
                     2*(12*G1(0) - G2(0))*pow(kq,2)))*pow(qe2,2))) + 
         Md*G1(qe2)*(Mp*(4*kq*nu*
                (-(kK*(4*kP - 5*kq)*kq*nu*(2*G1(0) - G2(0))) + 
                  Mp2*((-4*kP*G1(0) + 4*kq*G1(0) + kP*G2(0) + 
                        3*kq*G2(0))*pow(kK,2) + 
                     (4*kP*G1(0) + 4*kq*G1(0) - kP*G2(0) - 7*kq*G2(0))*
                      pow(kq,2))) - 
               4*me2*(kK*(kP + kq)*qe2*
                   (kP*kq*(34*G1(0) - 21*G2(0)) + 
                     (6*G1(0) - G2(0))*pow(kP,2) + 
                     2*(-34*G1(0) + 11*G2(0))*pow(kq,2)) + 
                  2*kq*(kq*(kP + kq)*nu*
                      (-14*kP*G1(0) + 22*kq*G1(0) + 5*kP*G2(0) - 
                        7*kq*G2(0)) - 
                     4*kK*Mp2*
                      (-3*kP*kq*G2(0) + (3*G1(0) - G2(0))*pow(kP,2) + 
                        (-12*G1(0) + G2(0))*pow(kq,2)))) + 
               qe2*(-4*kK*Mp2*
                   (-2*(kP*G2(0) + kq*(-11*G1(0) + G2(0)))*pow(kK,2) + 
                     kq*(kP*kq*(-10*G1(0) + 8*G2(0)) + 
                        (-7*G1(0) + 4*G2(0))*pow(kP,2) + 
                        2*(-11*G1(0) + G2(0))*pow(kq,2))) + 
                  nu*((kP + kq)*
                      (-22*kP*G1(0) + 54*kq*G1(0) + 5*kP*G2(0) - 
                        11*kq*G2(0))*pow(kq,2) + 
                     pow(kK,2)*
                      (kP*kq*(64*G1(0) - 38*G2(0)) + 6*G2(0)*pow(kP,2) + 
                        (-74*G1(0) + 25*G2(0))*pow(kq,2)))) + 
               kK*(kq*(30*G1(0) - 13*G2(0))*pow(kP,2) + 
                  (-6*G1(0) + 5*G2(0))*pow(kP,3) + 
                  4*kq*(15*G1(0) - 4*G2(0))*(pow(kK,2) - pow(kq,2)) - 
                  2*kP*(-7*G2(0)*pow(kK,2) + 
                     (12*G1(0) + G2(0))*pow(kq,2)))*pow(qe2,2)) + 
            Md*(2*nu*(-3*kK*kP*nu*G2(0) + 
                  4*(kP + kq)*me2*
                   (4*kP*G1(0) - 8*kq*G1(0) + kP*G2(0) + kq*G2(0)))*
                pow(kq,2) + qe2*
                (-4*kK*kP*(kP + kq)*me2*
                   (4*kP*G1(0) - 8*kq*G1(0) + kP*G2(0) - 5*kq*G2(0)) + 
                  nu*((kP + kq)*
                      (3*kP*(-4*G1(0) + G2(0)) + 
                        kq*(8*G1(0) + 17*G2(0)))*pow(kq,2) + 
                     3*pow(kK,2)*
                      (-4*kP*kq*G1(0) + (4*G1(0) - 2*G2(0))*pow(kP,2) + 
                        (-8*G1(0) + 3*G2(0))*pow(kq,2)))) - 
               4*Mp2*(kK*qe2*(2*(kP*G1(0) + 3*kq*(G1(0) - G2(0)))*
                      pow(kK,2) + 
                     kq*(6*kP*kq*(-G1(0) + G2(0)) + 
                        (G1(0) - 4*G2(0))*pow(kP,2) + 
                        6*(-G1(0) + G2(0))*pow(kq,2))) + 
                  kq*(-8*kK*me2*
                      (2*kP*kq*(-3*G1(0) + G2(0)) + 
                        (G1(0) + G2(0))*pow(kP,2) + 
                        (-4*G1(0) + G2(0))*pow(kq,2)) + 
                     nu*((2*kP*G1(0) - 4*kq*G1(0) + kP*G2(0) + 
                         3*kq*G2(0))*pow(kK,2) + 
                        (-2*kP*G1(0) + 4*kq*G1(0) - kP*G2(0) - kq*G2(0))*
                         pow(kq,2)))) + 
               kK*(2*(2*kP*G1(0) + 2*kq*G1(0) - kP*G2(0) - 7*kq*G2(0))*
                   pow(kK,2) + 
                  (kP + kq)*(-24*kP*kq*(G1(0) - G2(0)) + 
                     (4*G1(0) - 5*G2(0))*pow(kP,2) + 
                     2*(-2*G1(0) + 7*G2(0))*pow(kq,2)))*pow(qe2,2)))) + 
      qe2*F2(qp2)*(2*(Md - Mp)*Mp*(G1(0) - G2(0))*G3(qe2)*
          (2*kK*kP*qe2*(8*me2*Mp2 + qe2*(2*Mp2 + qe2)) - 
            kK*(kP + kq)*qe2*pow(nu,2) - pow(kq,2)*pow(nu,3) + 
            nu*(4*kq*((kP + kq)*me2 + kq*Mp2)*qe2 + 
               16*me2*Mp2*pow(kq,2) + 
               (kq*(-kP + kq) + pow(kK,2))*pow(qe2,2))) - 
         Md*G2(qe2)*(-4*pow(Mp,3)*
             (kK*kP*qe2*(4*me2 + 3*qe2)*(G1(0) - G2(0)) + 
               nu*(-4*kq*me2*(kq*G1(0) + kP*(G1(0) - G2(0))) + 
                  qe2*(kq*(kP*G1(0) + 5*kq*G1(0) - kP*G2(0) - 
                        4*kq*G2(0)) + (-2*G1(0) + G2(0))*pow(kK,2))) + 
               kK*kq*(G1(0) - G2(0))*pow(nu,2)) + 
            4*Md*Mp2*(4*me2*(kK*kP*qe2*(G1(0) - G2(0)) + 
                  kq*nu*(kq*G1(0) + kP*(-G1(0) + G2(0)))) + 
               nu*qe2*(kP*kq*(G1(0) - G2(0)) + G2(0)*pow(kK,2) + 
                  (3*G1(0) - 4*G2(0))*pow(kq,2)) + 
               kK*kq*(G1(0) - G2(0))*pow(nu,2) + 
               3*kK*kP*(G1(0) - G2(0))*pow(qe2,2)) + 
            Md*nu*(G1(0) - 2*G2(0))*
             ((4*kq*(kP + kq)*me2 + kK*(kP - kq)*nu)*qe2 + 
               pow(kq,2)*pow(nu,2) + 
               (-(kq*(kP + kq)) + pow(kK,2))*pow(qe2,2)) - 
            Mp*nu*(3*G1(0) - 2*G2(0))*
             ((4*kq*(kP + kq)*me2 + kK*(kP - kq)*nu)*qe2 + 
               pow(kq,2)*pow(nu,2) + 
               (-(kq*(kP + kq)) + pow(kK,2))*pow(qe2,2))) + 
         Md*G1(qe2)*(Mp*(-4*Mp2*
                (nu*(4*kq*me2*
                      (kq*(-4*G1(0) + G2(0)) + kP*(-G1(0) + G2(0))) + 
                     qe2*(kq*
                        (kP*G1(0) + 6*kq*G1(0) - kP*G2(0) - 3*kq*G2(0)) \
+ (-4*G1(0) + G2(0))*pow(kK,2))) + kK*kq*(G1(0) - G2(0))*pow(nu,2) + 
                  2*kK*kP*(G1(0) - G2(0))*pow(qe2,2)) - 
               3*nu*(2*G1(0) - G2(0))*
                ((4*kq*(kP + kq)*me2 + kK*(kP - kq)*nu)*qe2 + 
                  pow(kq,2)*pow(nu,2) + 
                  (-(kq*(kP + kq)) + pow(kK,2))*pow(qe2,2))) + 
            Md*(4*Mp2*(nu*(4*kq*me2*
                      (kP*(-G1(0) + G2(0)) + kq*(2*G1(0) + G2(0))) + 
                     qe2*(kq*(kP*G1(0) - kP*G2(0) - 3*kq*G2(0)) + 
                        (2*G1(0) + G2(0))*pow(kK,2))) + 
                  kK*kq*(G1(0) - G2(0))*pow(nu,2) + 
                  2*kK*kP*(G1(0) - G2(0))*pow(qe2,2)) - 
               3*nu*G2(0)*((4*kq*(kP + kq)*me2 + kK*(kP - kq)*nu)*qe2 + 
                  pow(kq,2)*pow(nu,2) + 
                  (-(kq*(kP + kq)) + pow(kK,2))*pow(qe2,2))))) + 
      2*Mp*F1(qp2)*(-2*Md*(Md - Mp)*Mp*(G1(0) - G2(0))*G2(qe2)*
          (nu*qe2*(4*kq*((kP + 2*kq)*me2 + 2*kq*Mp2) + 
               qe2*(-(kP*kq) + pow(kK,2))) + 
            kK*(2*kP - kq)*qe2*pow(nu,2) + 2*pow(kq,2)*pow(nu,3) + 
            kK*kP*(4*me2 + 8*Mp2 + qe2)*pow(qe2,2)) + 
         4*Md*Mp*G1(qe2)*(Mp*(4*Mp2*
                (nu*(4*me2*G1(0)*pow(kq,2) + 
                     qe2*(G1(0)*pow(kK,2) + 
                        (-2*G1(0) + G2(0))*pow(kq,2))) + 
                  kK*kP*(-G1(0) + G2(0))*pow(qe2,2)) - 
               nu*(2*G1(0) - G2(0))*
                ((4*kq*(kP + kq)*me2 + kK*(kP - kq)*nu)*qe2 + 
                  pow(kq,2)*pow(nu,2) + 
                  (-(kq*(kP + kq)) + pow(kK,2))*pow(qe2,2))) + 
            Md*(4*Mp2*(nu*(qe2*G1(0)*pow(kK,2) + 
                     4*me2*G1(0)*pow(kq,2) - qe2*G2(0)*pow(kq,2)) + 
                  kK*kP*(G1(0) - G2(0))*pow(qe2,2)) - 
               nu*G2(0)*((4*kq*(kP + kq)*me2 + kK*(kP - kq)*nu)*qe2 + 
                  pow(kq,2)*pow(nu,2) + 
                  (-(kq*(kP + kq)) + pow(kK,2))*pow(qe2,2)))) - 
         qe2*G3(qe2)*(Md*(4*Mp2*
                (-(kK*kP*qe2*(4*me2 + 3*qe2)*(G1(0) - G2(0))) + 
                  nu*(4*kq*me2*
                      (-2*kP*G1(0) - 4*kq*G1(0) + 2*kP*G2(0) + 
                        3*kq*G2(0)) + 
                     qe2*(kq*(2*kP*G1(0) - 2*kP*G2(0) + kq*G2(0)) + 
                        (-3*G1(0) + 2*G2(0))*pow(kK,2))) + 
                  2*kK*kq*(G1(0) - G2(0))*pow(nu,2)) + 
               nu*(2*G1(0) - G2(0))*
                ((4*kq*(kP + kq)*me2 + kK*(kP - kq)*nu)*qe2 + 
                  pow(kq,2)*pow(nu,2) + 
                  (-(kq*(kP + kq)) + pow(kK,2))*pow(qe2,2))) + 
            Mp*(-4*Mp2*(-(kK*kP*qe2*(4*me2 + 3*qe2)*(G1(0) - G2(0))) + 
                  nu*(4*kq*me2*
                      (-2*kP*G1(0) - 2*kq*G1(0) + 2*kP*G2(0) + 
                        3*kq*G2(0)) + 
                     qe2*(kq*
                         (2*kP*G1(0) - 2*kq*G1(0) - 2*kP*G2(0) + 
                         kq*G2(0)) - (G1(0) - 2*G2(0))*pow(kK,2))) + 
                  2*kK*kq*(G1(0) - G2(0))*pow(nu,2)) + 
               nu*G2(0)*((4*kq*(kP + kq)*me2 + kK*(kP - kq)*nu)*qe2 + 
                  pow(kq,2)*pow(nu,2) + 
                  (-(kq*(kP + kq)) + pow(kK,2))*pow(qe2,2))))) - 
      F2(qp2)*(qe2*G3(qe2)*(8*
             (kq*nu*(G1(0) - G2(0))*(2*pow(kK,2) - 3*pow(kq,2)) + 
               kK*(4*kq*me2*(-2*kP*G1(0) - 3*kq*G1(0) + 2*kP*G2(0) + 
                     kq*G2(0)) + 
                  qe2*(2*kP*kq*(-G1(0) + G2(0)) + 
                     G2(0)*(-pow(kK,2) + pow(kq,2)))))*pow(Mp,4) - 
            2*Mp2*(-(nu*(8*(kP + kq)*me2*G1(0) + 
                    kK*nu*(5*G1(0) - 3*G2(0)))*pow(kq,2)) - 
               qe2*(4*kK*(kP + kq)*me2*
                   (5*kP*G1(0) - 6*kP*G2(0) + kq*G2(0)) + 
                  nu*(2*(2*kq*(-G1(0) + G2(0)) + kP*(G1(0) + G2(0)))*
                      pow(kK,2) + 
                     kq*(2*kP*kq*(6*G1(0) - 7*G2(0)) + 
                        (-G1(0) + G2(0))*pow(kP,2) + 
                        (3*G1(0) - 5*G2(0))*pow(kq,2)))) + 
               kK*(kP*kq*(3*G1(0) - G2(0)) + 
                  (G1(0) - 3*G2(0))*pow(kK,2) - 
                  8*(G1(0) - G2(0))*pow(kP,2) - 
                  (G1(0) - 3*G2(0))*pow(kq,2))*pow(qe2,2)) + 
            2*Md*Mp*(4*Mp2*(kK*qe2*
                   (-(kq*(-2*kP*G1(0) + 4*kq*G1(0) + 2*kP*G2(0) + 
                        kq*G2(0))) + (4*G1(0) + G2(0))*pow(kK,2)) + 
                  kq*(4*kK*me2*
                      (2*kP*G1(0) + 5*kq*G1(0) - 2*kP*G2(0) - kq*G2(0)) \
- nu*(G1(0) - G2(0))*(2*pow(kK,2) - 3*pow(kq,2)))) + 
               nu*(8*(kP + kq)*me2*(G1(0) - G2(0)) + 
                  kK*nu*(-7*G1(0) + 3*G2(0)))*pow(kq,2) - 
               qe2*(4*kK*(kP + kq)*me2*
                   (kP*G1(0) + 10*kq*G1(0) - 4*kP*G2(0) - 5*kq*G2(0)) + 
                  nu*(2*(5*kq*(-G1(0) + G2(0)) + kP*(3*G1(0) + G2(0)))*
                      pow(kK,2) + 
                     kq*(2*kP*kq*(5*G1(0) - 6*G2(0)) + 
                        (-G1(0) + G2(0))*pow(kP,2) + 
                        (7*G1(0) - 9*G2(0))*pow(kq,2)))) + 
               kK*(kP*kq*(5*G1(0) - 3*G2(0)) + 
                  (-7*G1(0) + 3*G2(0))*pow(kK,2) + 
                  (-8*G1(0) + 6*G2(0))*pow(kP,2) + 
                  (7*G1(0) - 3*G2(0))*pow(kq,2))*pow(qe2,2)) - 
            (kP + kq)*(nu*qe2*
                (-8*kq*(kP + kq)*me2*(G1(0) - 2*G2(0)) + 
                  qe2*(kq*(2*kP*G1(0) - 2*kq*G1(0) - 4*kP*G2(0) + 
                        kq*G2(0)) - 2*(G1(0) - 2*G2(0))*pow(kK,2))) + 
               kK*qe2*(2*kP*G1(0) + 2*kq*G1(0) - kP*G2(0) - 4*kq*G2(0))*
                pow(nu,2) + (2*G1(0) - G2(0))*pow(kq,2)*pow(nu,3) + 
               kK*kP*(-4*G1(0) + 5*G2(0))*pow(qe2,3))) + 
         Md*G1(qe2)*(Md*(-(nu*qe2*
                  (4*kq*(kP + kq)*me2*
                     (-(kP*(4*G1(0) + G2(0))) + kq*(8*G1(0) + 5*G2(0))) \
+ qe2*(kq*(kP + kq)*(4*kP*G1(0) - 2*kP*G2(0) - 15*kq*G2(0)) + 
                       (8*kP*G1(0) + 8*kq*G1(0) - kP*G2(0) + 
                        5*kq*G2(0))*pow(kK,2)))) + 
               16*kK*G1(0)*(qe2*(pow(kK,2) - pow(kq,2)) + 
                  8*me2*pow(kq,2))*pow(Mp,4) - 
               kK*kq*qe2*(4*kP*G1(0) + 4*kq*G1(0) + 7*kP*G2(0) - 
                  8*kq*G2(0))*pow(nu,2) - 6*G2(0)*pow(kq,3)*pow(nu,3) + 
               kK*(12*kP*(kP + kq)*me2*G2(0) + 
                  qe2*(kP*kq*(-8*G1(0) + 10*G2(0)) + 
                     (-8*G1(0) + 7*G2(0))*pow(kP,2) + 
                     3*G2(0)*(-pow(kK,2) + pow(kq,2))))*pow(qe2,2) - 
               4*Mp2*(2*nu*(kK*nu*G2(0) + 
                     4*me2*(kP*(G1(0) - G2(0)) - kq*(2*G1(0) + G2(0))))*
                   pow(kq,2) + 
                  qe2*(-4*kK*me2*
                      (-2*kP*kq*G1(0) + (G1(0) + G2(0))*pow(kP,2) - 
                        G2(0)*pow(kq,2)) + 
                     nu*((3*kP*G1(0) - kq*(6*G1(0) + G2(0)))*
                        pow(kK,2) + 
                        (-3*kP*G1(0) + 2*kq*G1(0) + kP*G2(0) + 
                        5*kq*G2(0))*pow(kq,2))) + 
                  kK*(6*kP*kq*(-G1(0) + G2(0)) + 
                     (G1(0) - 2*G2(0))*pow(kK,2) + 
                     (G1(0) - 2*G2(0))*pow(kP,2) - 
                     (G1(0) - 2*G2(0))*pow(kq,2))*pow(qe2,2))) + 
            Mp*(kq*nu*qe2*(kK*nu*
                   (-38*kP*G1(0) + 40*kq*G1(0) + 17*kP*G2(0) - 
                     22*kq*G2(0)) + 
                  4*(kP + kq)*me2*
                   (14*kP*G1(0) - 34*kq*G1(0) - 5*kP*G2(0) + 13*kq*G2(0))\
) + 16*kK*G1(0)*(qe2*(pow(kK,2) - pow(kq,2)) + 8*me2*pow(kq,2))*
                pow(Mp,4) + 12*(-2*G1(0) + G2(0))*pow(kq,3)*pow(nu,3) + 
               (nu*(kq*(kP + kq)*
                      (-8*kP*G1(0) + 42*kq*G1(0) + 2*kP*G2(0) - 
                        21*kq*G2(0)) + 
                     (26*kP*G1(0) - 46*kq*G1(0) - 17*kP*G2(0) + 
                        19*kq*G2(0))*pow(kK,2)) - 
                  12*kK*me2*(2*G1(0) - G2(0))*
                   (-3*kP*kq + pow(kP,2) - 4*pow(kq,2)))*pow(qe2,2) - 
               4*Mp2*(2*nu*(kK*nu*(2*G1(0) - G2(0)) + 
                     4*me2*(-10*kq*G1(0) + 4*kq*G2(0) + 
                        kP*(-G1(0) + G2(0))))*pow(kq,2) + 
                  qe2*(4*kK*me2*
                      (2*kP*kq*G1(0) + (-3*G1(0) + G2(0))*pow(kP,2) + 
                        (20*G1(0) - 7*G2(0))*pow(kq,2)) + 
                     nu*((3*kP*G1(0) - 20*kq*G1(0) + 7*kq*G2(0))*
                        pow(kK,2) + 
                        (kq*(18*G1(0) - 5*G2(0)) - kP*(G1(0) + G2(0)))*
                         pow(kq,2))) + 
                  kK*((15*G1(0) - 4*G2(0))*pow(kK,2) + 
                     (-3*G1(0) + 2*G2(0))*pow(kP,2) + 
                     (-15*G1(0) + 4*G2(0))*pow(kq,2))*pow(qe2,2)) + 
               kK*((kP + kq)*
                   (kP*(14*G1(0) - 11*G2(0)) + 9*kq*(-2*G1(0) + G2(0))) \
+ 9*(2*G1(0) - G2(0))*pow(kK,2))*pow(qe2,3))) + 
         Md*G2(qe2)*(4*Md*Mp2*
             (kq*nu*qe2*(-(kq*(kP + 8*kq)*(G1(0) - G2(0))) + 
                  (2*G1(0) - 3*G2(0))*pow(kK,2)) + 
               4*me2*(2*(kP + kq)*nu*(G1(0) - G2(0))*pow(kq,2) + 
                  kK*qe2*(2*kP*kq*(-G1(0) + G2(0)) + 
                     (G1(0) - G2(0))*pow(kP,2) + 
                     (-3*G1(0) + 2*G2(0))*pow(kq,2))) + 
               2*kK*(-G1(0) + G2(0))*pow(kq,2)*pow(nu,2) + 
               kK*kP*(-9*kq*G1(0) + 2*kP*(G1(0) - G2(0)) + 8*kq*G2(0))*
                pow(qe2,2)) + 4*pow(Mp,3)*
             (2*nu*(kK*nu*(G1(0) - G2(0)) + 
                  4*me2*(kq*(-4*G1(0) + G2(0)) + kP*(-G1(0) + G2(0))))*
                pow(kq,2) + qe2*
                (kq*nu*(kq*(kP*G1(0) + 14*kq*G1(0) - kP*G2(0) - 
                        8*kq*G2(0)) + (-10*G1(0) + 3*G2(0))*pow(kK,2)) - 
                  4*kK*me2*(2*kP*kq*(-G1(0) + G2(0)) + 
                     (G1(0) - G2(0))*pow(kP,2) + 
                     (-7*G1(0) + 2*G2(0))*pow(kq,2))) + 
               kK*(kP*kq*(7*G1(0) - 8*G2(0)) + 6*G1(0)*pow(kK,2) - 
                  2*(G1(0) - G2(0))*pow(kP,2) - 6*G1(0)*pow(kq,2))*
                pow(qe2,2)) + Md*
             (kq*nu*qe2*(kK*nu*
                   (-7*kP*G1(0) + 6*kq*G1(0) + 8*kP*G2(0) - 6*kq*G2(0)) \
+ 4*me2*(G1(0) - 2*G2(0))*(pow(kP,2) - pow(kq,2))) + 
               6*(-G1(0) + G2(0))*pow(kq,3)*pow(nu,3) + 
               (4*kK*kP*(kP + kq)*me2*(G1(0) - 2*G2(0)) + 
                  nu*((5*kP*G1(0) - kq*G1(0) - 4*kP*G2(0) + 
                        2*kq*G2(0))*pow(kK,2) + 
                     3*(kP + kq)*(3*G1(0) - 4*G2(0))*pow(kq,2)))*
                pow(qe2,2) + kK*
                ((kP + kq)*(7*kP*G1(0) + kq*G1(0) - 8*kP*G2(0) - 
                     2*kq*G2(0)) - (G1(0) - 2*G2(0))*pow(kK,2))*pow(qe2,3)\
) + Mp*(kq*nu*qe2*(-4*(kP + kq)*me2*
                   (5*kP*G1(0) - 13*kq*G1(0) + 4*kq*G2(0)) + 
                  kK*nu*(17*kP*G1(0) - 22*kq*G1(0) - 6*kP*G2(0) + 
                     14*kq*G2(0))) + 
               6*(2*G1(0) - G2(0))*pow(kq,3)*pow(nu,3) + 
               (-4*kK*(kP + kq)*me2*
                   (kP*G1(0) + 12*kq*G1(0) - 4*kP*G2(0) - 6*kq*G2(0)) + 
                  nu*(kq*(kP + kq)*
                      (-29*kq*G1(0) + 22*kq*G2(0) + 
                        2*kP*(G1(0) + G2(0))) + 
                     (-13*kP*G1(0) + 23*kq*G1(0) + 6*kP*G2(0) - 
                        12*kq*G2(0))*pow(kK,2)))*pow(qe2,2) + 
               kK*(-((kP + kq)*
                     (15*kP*G1(0) - 9*kq*G1(0) - 14*kP*G2(0) + 
                       4*kq*G2(0))) + (-9*G1(0) + 4*G2(0))*pow(kK,2))*
                pow(qe2,3)))) - 
      2*Mp*F1(qp2)*(qe2*G3(qe2)*
          (Mp*(nu*(16*kq*(kP + kq)*me2*
                   (2*kP*G1(0) + 2*kq*G1(0) - 2*kP*G2(0) - 3*kq*G2(0)) + 
                  qe2*(kq*(kP + kq)*
                      (-8*kP*G1(0) + 8*kq*G1(0) + 9*kP*G2(0) - 
                        5*kq*G2(0)) + 
                     2*(2*kP*G1(0) + 2*kq*G1(0) - 3*kP*G2(0) - 
                        5*kq*G2(0))*pow(kK,2))) + 
               kK*qe2*(4*(kP + kq)*me2*
                   (4*kP*(G1(0) - G2(0)) + kq*G2(0)) + 
                  qe2*(3*kP*kq*(4*G1(0) - 5*G2(0)) + 
                     2*(6*G1(0) - 7*G2(0))*pow(kP,2) + 
                     G2(0)*(pow(kK,2) - pow(kq,2)))) - 
               4*Mp2*(4*kK*me2*
                   (kP*kq*G1(0) - 2*(G1(0) - G2(0))*pow(kP,2) + 
                     (3*G1(0) - G2(0))*pow(kq,2)) + 
                  kK*qe2*(kP*kq*(5*G1(0) - 4*G2(0)) + 
                     3*(G1(0) - G2(0))*pow(kK,2) - 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     3*(-G1(0) + G2(0))*pow(kq,2)) + 
                  nu*(2*kP*(G1(0) - G2(0))*pow(kK,2) + 
                     (-6*kP*G1(0) + 2*kq*G1(0) + 6*kP*G2(0) - kq*G2(0))*
                      pow(kq,2))) + 
               32*kK*(G1(0) - G2(0))*(pow(kK,2) - pow(kq,2))*pow(Mp,4) + 
               kK*kq*(-8*kP*G1(0) - 8*kq*G1(0) + 8*kP*G2(0) + 
                  11*kq*G2(0))*pow(nu,2)) + 
            Md*(-(kq*nu*(8*(kP + kq)*me2*
                     (kP*G1(0) + 2*kq*G1(0) - 2*kP*G2(0) - kq*G2(0)) + 
                    kK*nu*(-2*kP*G1(0) + 4*kq*G1(0) + 4*kP*G2(0) + 
                       kq*G2(0)))) + 
               qe2*(4*kK*(kP + kq)*me2*
                   (5*kq*(-2*G1(0) + G2(0)) + 2*kP*(G1(0) + G2(0))) + 
                  nu*(kq*(kP + kq)*
                      (4*kP*G1(0) - 2*kq*G1(0) - 5*kP*G2(0) + 7*kq*G2(0)) \
+ (-14*kP*G1(0) + 2*kq*G1(0) + 4*kP*G2(0) - 4*kq*G2(0))*pow(kK,2))) + 
               4*Mp2*(-4*kK*me2*
                   (kP*kq*G1(0) + 2*(G1(0) - G2(0))*pow(kP,2) + 
                     (-5*G1(0) + G2(0))*pow(kq,2)) + 
                  kK*qe2*(kP*kq*(3*G1(0) - 4*G2(0)) + 
                     (9*G1(0) - 3*G2(0))*pow(kK,2) - 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     3*(-3*G1(0) + G2(0))*pow(kq,2)) + 
                  nu*(2*kP*(G1(0) - G2(0))*pow(kK,2) - 
                     (6*kP*(G1(0) - G2(0)) + kq*G2(0))*pow(kq,2))) - 
               32*kK*(G1(0) - G2(0))*(pow(kK,2) - pow(kq,2))*pow(Mp,4) + 
               kK*(-((kP + kq)*
                     (2*kP*(5*G1(0) - 4*G2(0)) + 5*kq*(-2*G1(0) + G2(0)))\
) + 5*(-2*G1(0) + G2(0))*pow(kK,2))*pow(qe2,2))) + 
         Md*G2(qe2)*(-8*kq*me2*nu*qe2*G1(0)*pow(kP,2) + 
            4*kq*me2*nu*qe2*G2(0)*pow(kP,2) + 
            12*kP*me2*nu*qe2*G1(0)*pow(kq,2) + 
            20*me2*nu*qe2*G1(0)*pow(kq,3) - 4*me2*nu*qe2*G2(0)*pow(kq,3) - 
            16*kK*(G1(0) - G2(0))*
             (-(kK*kq*nu) + qe2*(-(kq*(2*kP + kq)) + pow(kK,2)))*pow(Mp,4) \
+ 7*kK*kP*kq*qe2*G1(0)*pow(nu,2) - 2*kK*kP*kq*qe2*G2(0)*pow(nu,2) - 
            kK*qe2*G1(0)*pow(kP,2)*pow(nu,2) + 
            2*kK*qe2*G2(0)*pow(kP,2)*pow(nu,2) - 
            10*kK*qe2*G1(0)*pow(kq,2)*pow(nu,2) + 
            5*kK*qe2*G2(0)*pow(kq,2)*pow(nu,2) - 
            kP*G1(0)*pow(kq,2)*pow(nu,3) + 
            2*kP*G2(0)*pow(kq,2)*pow(nu,3) + 5*G1(0)*pow(kq,3)*pow(nu,3) - 
            G2(0)*pow(kq,3)*pow(nu,3) - 20*kK*kP*kq*me2*G1(0)*pow(qe2,2) + 
            16*kK*kP*kq*me2*G2(0)*pow(qe2,2) - 
            8*kP*nu*G1(0)*pow(kK,2)*pow(qe2,2) + 
            10*kq*nu*G1(0)*pow(kK,2)*pow(qe2,2) + 
            4*kP*nu*G2(0)*pow(kK,2)*pow(qe2,2) - 
            5*kq*nu*G2(0)*pow(kK,2)*pow(qe2,2) + 
            4*kK*me2*G1(0)*pow(kP,2)*pow(qe2,2) + 
            2*kq*nu*G1(0)*pow(kP,2)*pow(qe2,2) + 
            4*kK*me2*G2(0)*pow(kP,2)*pow(qe2,2) - 
            kq*nu*G2(0)*pow(kP,2)*pow(qe2,2) - 
            24*kK*me2*G1(0)*pow(kq,2)*pow(qe2,2) - 
            7*kP*nu*G1(0)*pow(kq,2)*pow(qe2,2) + 
            12*kK*me2*G2(0)*pow(kq,2)*pow(qe2,2) + 
            5*kP*nu*G2(0)*pow(kq,2)*pow(qe2,2) - 
            9*nu*G1(0)*pow(kq,3)*pow(qe2,2) + 
            6*nu*G2(0)*pow(kq,3)*pow(qe2,2) + 
            2*Md*Mp*(2*kq*(kP + kq)*nu*
                (kK*nu*(-G1(0) + G2(0)) + 
                  2*me2*(2*kP*G1(0) + 3*kq*G1(0) - 2*kP*G2(0) - 
                     4*kq*G2(0))) + 
               8*kK*Mp2*(G1(0) - G2(0))*
                (-(kK*kq*nu) + qe2*(-(kq*(2*kP + kq)) + pow(kK,2))) + 
               qe2*(4*kK*me2*
                   (-(kP*kq*G1(0)) + 3*(G1(0) - G2(0))*pow(kP,2) + 
                     (-3*G1(0) + 2*G2(0))*pow(kq,2)) + 
                  nu*((kP*G1(0) + 5*kq*G1(0) - 2*kP*G2(0) - 7*kq*G2(0))*
                      pow(kK,2) + 
                     (5*kP*G1(0) + 3*kq*G1(0) - 5*kP*G2(0) - 4*kq*G2(0))*
                      pow(kq,2))) + 
               kK*(kP*kq*(7*G1(0) - 10*G2(0)) + 
                  (-3*G1(0) + 4*G2(0))*pow(kK,2) + 
                  (7*G1(0) - 8*G2(0))*pow(kP,2) + 
                  (3*G1(0) - 4*G2(0))*pow(kq,2))*pow(qe2,2)) + 
            2*Mp2*(2*kq*nu*(kK*(kP + kq)*nu*(G1(0) - G2(0)) + 
                  2*me2*(kP*kq*(-7*G1(0) + 6*G2(0)) - 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     (-11*G1(0) + 4*G2(0))*pow(kq,2))) + 
               qe2*(-4*kK*me2*
                   (kP*kq*G1(0) + 3*(G1(0) - G2(0))*pow(kP,2) + 
                     (-7*G1(0) + 2*G2(0))*pow(kq,2)) + 
                  nu*((-3*kP*G1(0) - 15*kq*G1(0) + 2*kP*G2(0) + 
                        7*kq*G2(0))*pow(kK,2) + 
                     (-7*kP*G1(0) + 5*kq*G1(0) + 9*kP*G2(0) + 
                        2*kq*G2(0))*pow(kq,2))) + 
               kK*(kP*kq*(-9*G1(0) + 8*G2(0)) + 
                  (11*G1(0) - 4*G2(0))*pow(kK,2) + 
                  (-11*G1(0) + 12*G2(0))*pow(kP,2) + 
                  (-11*G1(0) + 4*G2(0))*pow(kq,2))*pow(qe2,2)) + 
            kK*kP*kq*G1(0)*pow(qe2,3) + kK*kP*kq*G2(0)*pow(qe2,3) - 
            6*G1(0)*pow(kK,3)*pow(qe2,3) + 3*G2(0)*pow(kK,3)*pow(qe2,3) - 
            5*kK*G1(0)*pow(kP,2)*pow(qe2,3) + 
            4*kK*G2(0)*pow(kP,2)*pow(qe2,3) + 
            6*kK*G1(0)*pow(kq,2)*pow(qe2,3) - 
            3*kK*G2(0)*pow(kq,2)*pow(qe2,3)) + 
         Md*G1(qe2)*(8*Md*(kK*
                (8*kq*(-kP + kq)*me2*G1(0) + 
                  qe2*(2*kP*kq*(G1(0) - 2*G2(0)) + 
                     3*G2(0)*(pow(kK,2) - pow(kq,2)))) - 
               2*kq*nu*((-G1(0) + G2(0))*pow(kK,2) + G1(0)*pow(kq,2)))*
             pow(Mp,3) + 8*(kK*qe2*
                (kq*(-6*kP*G1(0) - 6*kq*G1(0) + 4*kP*G2(0) + 
                     3*kq*G2(0)) + (6*G1(0) - 3*G2(0))*pow(kK,2)) - 
               2*kq*(4*kK*(kP - kq)*me2*G1(0) + 
                  nu*((G1(0) - G2(0))*pow(kK,2) + G1(0)*pow(kq,2))))*
             pow(Mp,4) - 4*Md*Mp*
             (2*kK*qe2*(qe2*(kP*kq*(2*G1(0) - 3*G2(0)) + 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     G2(0)*(pow(kK,2) - pow(kq,2))) + 
                  me2*G2(0)*(-3*kP*kq - 4*pow(kP,2) + pow(kq,2))) + 
               2*nu*(kq*(kP + kq)*me2*
                   (8*kq*G1(0) - 2*kP*G2(0) - 3*kq*G2(0)) + 
                  qe2*((2*kP*G1(0) + 2*kq*G1(0) - kq*G2(0))*pow(kK,2) - 
                     2*(kP + kq)*G2(0)*pow(kq,2))) + 
               kK*kq*(kP + kq)*G2(0)*pow(nu,2)) - 
            4*Mp2*(kq*nu*(kK*(kP + kq)*nu*(2*G1(0) - G2(0)) + 
                  2*me2*(7*kP*kq*G2(0) + (-4*G1(0) + 2*G2(0))*pow(kP,2) + 
                     (-14*G1(0) + 11*G2(0))*pow(kq,2))) + 
               qe2*(2*kK*me2*(3*kP*kq*(-2*G1(0) + G2(0)) + 
                     (-8*G1(0) + 4*G2(0))*pow(kP,2) + 
                     (20*G1(0) - 7*G2(0))*pow(kq,2)) + 
                  nu*((5*kP*G1(0) - 8*kq*G1(0) + kP*G2(0) + 6*kq*G2(0))*
                      pow(kK,2) + 
                     (-10*kP*G1(0) - kq*G1(0) + 5*kP*G2(0) + 2*kq*G2(0))*
                      pow(kq,2))) + 
               kK*(kP*kq*(-9*G1(0) + 8*G2(0)) + 
                  (13*G1(0) - 5*G2(0))*pow(kK,2) + 
                  (-5*G1(0) + 6*G2(0))*pow(kP,2) + 
                  (-13*G1(0) + 5*G2(0))*pow(kq,2))*pow(qe2,2)) + 
            (2*G1(0) - G2(0))*(nu*qe2*
                (-4*kq*me2*(kP*kq - 4*pow(kP,2) + 5*pow(kq,2)) + 
                  kK*nu*(-9*kP*kq + pow(kP,2) + 8*pow(kq,2))) + 
               (kP - 5*kq)*pow(kq,2)*pow(nu,3) + 
               (12*kK*me2*(kP*kq - pow(kP,2) + 2*pow(kq,2)) + 
                  nu*(2*(5*kP - 4*kq)*pow(kK,2) + 
                     kq*(kP*kq - 4*pow(kP,2) + 5*pow(kq,2))))*pow(qe2,2) + 
               3*kK*(-(kP*kq) + 2*pow(kK,2) + pow(kP,2) - 2*pow(kq,2))*
                pow(qe2,3)))))*pow(qp2,-1)*
    pow(Gd2*Md2 + pow(kP + kq - Md2 + Mp2,2),-1))/36.
);
}
double iMdiMe2(double nu, double qe2, double qp2, double kK, double kP, double kq) {
return ((e6*Z3*pow(kK - kq,-1)*pow(kK + kq,-1)*pow(Md,-5)*pow(Mp,-1)*
    pow(-kP + kq - Md2 + Mp2,-1)*pow(qe2,-1)*
    (kK*(-2*kP*Md*Mp*F1(qp2)*(G2(qe2)*
             (kq*(-(kP*(G1(0) - 2*G2(0))) + kq*(G1(0) + G2(0)))*
                pow(kP - kq,2) + 
               pow(kK,2)*(kP*kq*(11*G1(0) - 4*G2(0)) - 
                  6*(G1(0) - G2(0))*pow(kP,2) + 
                  (-11*G1(0) + G2(0))*pow(kq,2))) + 
            G1(qe2)*(pow(kK,2)*
                (kP*kq*(2*G1(0) - G2(0)) - 6*G2(0)*pow(kP,2) + 
                  (10*G1(0) + G2(0))*pow(kq,2)) + 
               (kP - kq)*kq*(2*kP*kq*(4*G1(0) + G2(0)) + 
                  (2*G1(0) - G2(0))*pow(kP,2) + 
                  (-10*G1(0) + 11*G2(0))*pow(kq,2)))) + 
         F2(qp2)*(-((kP - kq)*qe2*G3(qe2)*
               (2*pow(kq,2)*(3*kP*kq*G2(0) + 
                    2*(G1(0) - 2*G2(0))*pow(kP,2) + 
                    (10*G1(0) - 11*G2(0))*pow(kq,2)) + 
                 pow(kK,2)*(4*kP*kq*(G1(0) - 2*G2(0)) + 
                    (6*G1(0) - 3*G2(0))*pow(kP,2) + 
                    2*(-10*G1(0) + 11*G2(0))*pow(kq,2)))) + 
            2*kP*Md*G1(qe2)*(2*Mp*
                (-((kP - kq)*
                     (4*kP*G1(0) - 4*kq*G1(0) - kP*G2(0) + 7*kq*G2(0))*
                     pow(kq,2)) + 
                  pow(kK,2)*(kP*kq*(-2*G1(0) + 5*G2(0)) + 
                     G2(0)*pow(kP,2) - (4*G1(0) + 3*G2(0))*pow(kq,2))) + 
               Md*(-((kP - kq)*
                     (4*kP*G1(0) + 8*kq*G1(0) + kP*G2(0) - kq*G2(0))*
                     pow(kq,2)) + 
                  pow(kK,2)*(kP*kq*(4*G1(0) - 2*G2(0)) + 
                     (4*G1(0) - 2*G2(0))*pow(kP,2) + 
                     (-8*G1(0) + 7*G2(0))*pow(kq,2)))) - 
            2*kP*Md*G2(qe2)*(Mp*
                (-((2*G1(0) - 3*G2(0))*pow(kP - kq,2)*pow(kq,2)) + 
                  pow(kK,2)*(2*kP*kq*(G1(0) + 3*G2(0)) - 
                     2*(G1(0) - G2(0))*pow(kP,2) - 
                     (6*G1(0) + 5*G2(0))*pow(kq,2))) + 
               Md*((G1(0) - 2*G2(0))*pow(kP - kq,2)*pow(kq,2) + 
                  pow(kK,2)*(2*kP*kq*(G1(0) - 2*G2(0)) + 
                     (2*G1(0) - G2(0))*pow(kP,2) + 
                     (-7*G1(0) + 8*G2(0))*pow(kq,2)))))) + 
      F2(qp2)*(qe2*G3(qe2)*(2*kK*Md*Mp*
             (-2*(kP - kq)*(6*kP*G1(0) + 5*kq*G1(0) - 2*kP*G2(0) + 
                  3*kq*G2(0))*pow(kq,2) + 
               pow(kK,2)*(4*kP*kq*(G1(0) + 2*G2(0)) - 
                  3*(G1(0) - G2(0))*pow(kP,2) - 
                  2*(5*G1(0) + 3*G2(0))*pow(kq,2))) + 
            2*kK*Mp2*(2*pow(kq,2)*
                (kP*kq*(-6*G1(0) + 8*G2(0)) + 
                  (G1(0) - 2*G2(0))*pow(kP,2) + 
                  (11*G1(0) - 12*G2(0))*pow(kq,2)) + 
               pow(kK,2)*(2*kP*kq*(5*G1(0) - 6*G2(0)) + 
                  (7*G1(0) - 5*G2(0))*pow(kP,2) + 
                  2*(-11*G1(0) + 12*G2(0))*pow(kq,2))) + 
            (kP - kq)*(16*kK*kP*(kP - kq)*kq*me2*(G1(0) - 2*G2(0)) + 
               kK*qe2*((-14*kq*G1(0) + 2*kP*(G1(0) - 2*G2(0)) + 
                     19*kq*G2(0))*pow(kK,2) + 
                  kq*(3*kP*kq*(4*G1(0) - 3*G2(0)) + 
                     18*(G1(0) - G2(0))*pow(kP,2) + 
                     (14*G1(0) - 19*G2(0))*pow(kq,2))) + 
               nu*(pow(kK,2)*(-4*kP*kq*(G1(0) - 2*G2(0)) + 
                     (-2*G1(0) + G2(0))*pow(kP,2) + 
                     2*(4*G1(0) - 5*G2(0))*pow(kq,2)) - 
                  2*pow(kq,2)*
                   (-3*kP*kq*(4*G1(0) - 5*G2(0)) + 
                     (-2*G1(0) + G2(0))*pow(kP,2) + 
                     2*(G1(0) - 2*G2(0))*pow(kq,2))))) + 
         Md*G2(qe2)*(Mp*(-4*kK*kP*Mp2*
                ((kP*(G1(0) - G2(0)) - 3*kq*(G1(0) + G2(0)))*
                   pow(kK,2) + (kP - kq)*(G1(0) - G2(0))*pow(kq,2)) + 
               kK*qe2*(pow(kK,2)*
                   (kP*kq*(-21*G1(0) + 16*G2(0)) + 
                     12*G2(0)*pow(kP,2) - 12*(G1(0) + G2(0))*pow(kq,2)) \
- (kP - kq)*kq*(3*kP*kq*(9*G1(0) - 8*G2(0)) + 
                     (13*G1(0) - 18*G2(0))*pow(kP,2) + 
                     12*(G1(0) + G2(0))*pow(kq,2))) - 
               2*kq*(4*kK*(kP - kq)*me2*
                   (kP*kq*(-15*G1(0) + 8*G2(0)) + 
                     (G1(0) - 4*G2(0))*pow(kP,2) + 
                     2*(G1(0) - 2*G2(0))*pow(kq,2)) + 
                  nu*((2*G1(0) - 3*G2(0))*pow(kP - kq,2)*pow(kq,2) + 
                     pow(kK,2)*
                      (2*kP*kq*(4*G1(0) - 7*G2(0)) + 4*G1(0)*pow(kP,2) + 
                        (6*G1(0) + 5*G2(0))*pow(kq,2))))) + 
            Md*(4*kK*kP*Mp2*(G1(0) - G2(0))*
                ((kP + 3*kq)*pow(kK,2) + (kP - kq)*pow(kq,2)) + 
               kK*qe2*((kP - kq)*kq*
                   (3*kP*kq*(5*G1(0) - 6*G2(0)) + 
                     (11*G1(0) - 10*G2(0))*pow(kP,2) + 
                     12*(-G1(0) + G2(0))*pow(kq,2)) + 
                  pow(kK,2)*(kP*kq*(-7*G1(0) + 8*G2(0)) + 
                     6*(G1(0) - G2(0))*pow(kP,2) + 
                     12*(-G1(0) + G2(0))*pow(kq,2))) + 
               2*kq*(4*kK*kP*me2*(G1(0) - 2*G2(0))*pow(kP - kq,2) + 
                  nu*((G1(0) - 2*G2(0))*pow(kP - kq,2)*pow(kq,2) + 
                     pow(kK,2)*
                      (2*kP*kq*(G1(0) - 2*G2(0)) + 
                        (2*G1(0) - G2(0))*pow(kP,2) + 
                        (-7*G1(0) + 8*G2(0))*pow(kq,2)))))) + 
         Md*G1(qe2)*(Md*(-4*kK*kP*Mp2*
                ((2*kP*G1(0) + 4*kq*G1(0) - kP*G2(0) - 3*kq*G2(0))*
                   pow(kK,2) - 
                  (2*kP*G1(0) + 4*kq*G1(0) + kP*G2(0) - kq*G2(0))*
                   pow(kq,2)) + 
               kK*qe2*(-((kP - kq)*kq*
                     (kP*kq*(8*G1(0) - 19*G2(0)) + 
                       (4*G1(0) - 11*G2(0))*pow(kP,2) - 
                       8*(G1(0) - 2*G2(0))*pow(kq,2))) + 
                  pow(kK,2)*(-3*kP*kq*G2(0) + 
                     (-8*G1(0) + 4*G2(0))*pow(kP,2) + 
                     8*(G1(0) - 2*G2(0))*pow(kq,2))) + 
               2*kq*(4*kK*kP*(kP - kq)*me2*
                   (4*kP*G1(0) + 8*kq*G1(0) + kP*G2(0) - kq*G2(0)) + 
                  nu*((kP - kq)*
                      (4*kP*G1(0) + 8*kq*G1(0) + kP*G2(0) - kq*G2(0))*
                      pow(kq,2) + 
                     pow(kK,2)*
                      (2*kP*kq*(-2*G1(0) + G2(0)) + 
                        (-4*G1(0) + 2*G2(0))*pow(kP,2) + 
                        (8*G1(0) - 7*G2(0))*pow(kq,2))))) + 
            Mp*(kK*qe2*((kP - kq)*kq*
                   (kP*kq*(6*G1(0) + 17*G2(0)) + 
                     (22*G1(0) - 13*G2(0))*pow(kP,2) + 
                     4*(12*G1(0) + G2(0))*pow(kq,2)) + 
                  pow(kK,2)*(kP*kq*(10*G1(0) - 17*G2(0)) + 
                     4*(2*G1(0) - 5*G2(0))*pow(kP,2) + 
                     4*(12*G1(0) + G2(0))*pow(kq,2))) + 
               4*(-(kK*kP*Mp2*
                     ((4*kq*G1(0) + kP*G2(0) + 3*kq*G2(0))*pow(kK,2) + 
                       (-4*kP*G1(0) + 4*kq*G1(0) + kP*G2(0) - 7*kq*G2(0))*
                        pow(kq,2))) + 
                  kq*(2*kK*(kP - kq)*me2*
                      (kP*kq*(-22*G1(0) + 15*G2(0)) + 
                        (6*G1(0) - G2(0))*pow(kP,2) - 
                        2*(10*G1(0) + G2(0))*pow(kq,2)) + 
                     nu*((kP - kq)*
                         (4*kP*G1(0) - 4*kq*G1(0) - kP*G2(0) + 
                         7*kq*G2(0))*pow(kq,2) + 
                        pow(kK,2)*
                         (2*kP*kq*(6*G1(0) - 5*G2(0)) + 
                          2*(G1(0) - G2(0))*pow(kP,2) + 
                          (4*G1(0) + 3*G2(0))*pow(kq,2)))))))) + 
      2*Mp*F1(qp2)*(-(kK*qe2*G3(qe2)*
            (Md*(2*(kP - kq)*kq*
                  (3*kP*kq*G1(0) + (G1(0) - 2*G2(0))*pow(kP,2) + 
                    (8*G1(0) - G2(0))*pow(kq,2)) + 
                 pow(kK,2)*(kP*kq*(-14*G1(0) + 4*G2(0)) + 
                    (8*G1(0) - 7*G2(0))*pow(kP,2) + 
                    2*(8*G1(0) - G2(0))*pow(kq,2))) + 
              Mp*(-2*(kP - kq)*kq*
                  (3*kP*kq*G2(0) + 4*(G1(0) - G2(0))*pow(kP,2) + 
                    (12*G1(0) - 13*G2(0))*pow(kq,2)) + 
                 pow(kK,2)*(kP*kq*(32*G1(0) - 38*G2(0)) + 
                    (-8*G1(0) + 11*G2(0))*pow(kP,2) + 
                    2*(-12*G1(0) + 13*G2(0))*pow(kq,2))))) + 
         Md*G1(qe2)*(-4*kK*kP*Md*Mp*
             (-((kP - kq)*kq*(4*kq*G1(0) + kP*G2(0))) + 
               (4*kP*G1(0) - 4*kq*G1(0) - 2*kP*G2(0) + 4*kq*G2(0))*
                pow(kK,2)) - 2*kP*kq*qe2*G1(0)*pow(kK,3) - 
            5*kP*kq*qe2*G2(0)*pow(kK,3) + 
            4*kq*nu*G1(0)*pow(kK,2)*pow(kP,2) - 
            8*kq*nu*G2(0)*pow(kK,2)*pow(kP,2) + 
            6*qe2*G1(0)*pow(kK,3)*pow(kP,2) - 
            9*qe2*G2(0)*pow(kK,3)*pow(kP,2) + 
            16*kK*kq*me2*G1(0)*pow(kP,3) + 10*kK*kq*qe2*G1(0)*pow(kP,3) + 
            16*kK*kq*me2*G2(0)*pow(kP,3) - 11*kK*kq*qe2*G2(0)*pow(kP,3) + 
            8*kK*me2*G1(0)*pow(kP,4) + 2*kK*qe2*G1(0)*pow(kP,4) - 
            4*kK*me2*G2(0)*pow(kP,4) - kK*qe2*G2(0)*pow(kP,4) + 
            22*kP*nu*G1(0)*pow(kK,2)*pow(kq,2) - 
            11*kP*nu*G2(0)*pow(kK,2)*pow(kq,2) + 
            32*qe2*G1(0)*pow(kK,3)*pow(kq,2) - 
            4*qe2*G2(0)*pow(kK,3)*pow(kq,2) - 
            120*kK*me2*G1(0)*pow(kP,2)*pow(kq,2) - 
            22*kK*qe2*G1(0)*pow(kP,2)*pow(kq,2) + 
            36*kK*me2*G2(0)*pow(kP,2)*pow(kq,2) + 
            29*kK*qe2*G2(0)*pow(kP,2)*pow(kq,2) + 
            2*nu*G1(0)*pow(kP,3)*pow(kq,2) - 
            nu*G2(0)*pow(kP,3)*pow(kq,2) + 
            4*kK*kP*Mp2*((-5*kq*G1(0) + kP*(G1(0) - 4*G2(0)) + 
                  3*kq*G2(0))*pow(kK,2) + 
               kq*(kP*kq*(3*G1(0) + 2*G2(0)) + 
                  (2*G1(0) - G2(0))*pow(kP,2) + 
                  (-5*G1(0) + 2*G2(0))*pow(kq,2))) + 
            16*kK*kP*me2*G1(0)*pow(kq,3) + 42*kK*kP*qe2*G1(0)*pow(kq,3) - 
            56*kK*kP*me2*G2(0)*pow(kq,3) - 21*kK*kP*qe2*G2(0)*pow(kq,3) + 
            10*nu*G1(0)*pow(kK,2)*pow(kq,3) + 
            nu*G2(0)*pow(kK,2)*pow(kq,3) + 
            6*nu*G1(0)*pow(kP,2)*pow(kq,3) + 
            3*nu*G2(0)*pow(kP,2)*pow(kq,3) + 80*kK*me2*G1(0)*pow(kq,4) - 
            18*kP*nu*G1(0)*pow(kq,4) - 32*kK*qe2*G1(0)*pow(kq,4) + 
            8*kK*me2*G2(0)*pow(kq,4) + 9*kP*nu*G2(0)*pow(kq,4) + 
            4*kK*qe2*G2(0)*pow(kq,4) + 10*nu*G1(0)*pow(kq,5) - 
            11*nu*G2(0)*pow(kq,5)) + 
         Md*G2(qe2)*(-2*kK*kP*Mp2*
             ((kP - kq)*kq*(2*kP*G1(0) - 2*kq*G1(0) - 2*kP*G2(0) + 
                  kq*G2(0)) + 
               (8*kP*G1(0) - 12*kq*G1(0) - 9*kP*G2(0) + 7*kq*G2(0))*
                pow(kK,2)) + 2*kK*kP*Md*Mp*
             ((kP - kq)*kq*(2*kP*G1(0) - 2*kP*G2(0) + kq*G2(0)) + 
               (4*kP*G1(0) - 8*kq*G1(0) - 5*kP*G2(0) + 9*kq*G2(0))*
                pow(kK,2)) - 11*kP*kq*qe2*G1(0)*pow(kK,3) + 
            7*kP*kq*qe2*G2(0)*pow(kK,3) - 
            8*kq*nu*G1(0)*pow(kK,2)*pow(kP,2) + 
            4*kq*nu*G2(0)*pow(kK,2)*pow(kP,2) + 
            qe2*G1(0)*pow(kK,3)*pow(kP,2) + 
            4*qe2*G2(0)*pow(kK,3)*pow(kP,2) + 
            16*kK*kq*me2*G1(0)*pow(kP,3) - 9*kK*kq*qe2*G1(0)*pow(kP,3) - 
            8*kK*kq*me2*G2(0)*pow(kP,3) + 12*kK*kq*qe2*G2(0)*pow(kP,3) - 
            4*kK*me2*G1(0)*pow(kP,4) - kK*qe2*G1(0)*pow(kP,4) + 
            8*kK*me2*G2(0)*pow(kP,4) + 2*kK*qe2*G2(0)*pow(kP,4) + 
            kP*nu*G1(0)*pow(kK,2)*pow(kq,2) + 
            4*kP*nu*G2(0)*pow(kK,2)*pow(kq,2) - 
            8*qe2*G1(0)*pow(kK,3)*pow(kq,2) - 
            2*qe2*G2(0)*pow(kK,3)*pow(kq,2) + 
            36*kK*me2*G1(0)*pow(kP,2)*pow(kq,2) + 
            5*kK*qe2*G1(0)*pow(kP,2)*pow(kq,2) - 
            24*kK*me2*G2(0)*pow(kP,2)*pow(kq,2) - 
            13*kK*qe2*G2(0)*pow(kP,2)*pow(kq,2) - 
            nu*G1(0)*pow(kP,3)*pow(kq,2) + 
            2*nu*G2(0)*pow(kP,3)*pow(kq,2) - 
            56*kK*kP*me2*G1(0)*pow(kq,3) - 3*kK*kP*qe2*G1(0)*pow(kq,3) + 
            40*kK*kP*me2*G2(0)*pow(kq,3) - 3*kK*kP*qe2*G2(0)*pow(kq,3) - 
            11*nu*G1(0)*pow(kK,2)*pow(kq,3) + 
            nu*G2(0)*pow(kK,2)*pow(kq,3) + 
            3*nu*G1(0)*pow(kP,2)*pow(kq,3) - 
            3*nu*G2(0)*pow(kP,2)*pow(kq,3) + 8*kK*me2*G1(0)*pow(kq,4) - 
            3*kP*nu*G1(0)*pow(kq,4) + 8*kK*qe2*G1(0)*pow(kq,4) - 
            16*kK*me2*G2(0)*pow(kq,4) + 2*kK*qe2*G2(0)*pow(kq,4) + 
            nu*G1(0)*pow(kq,5) + nu*G2(0)*pow(kq,5))) + 
      2*Mp*F1(qp2)*(qe2*G3(qe2)*
          (Mp*(-8*kP*kq*nu*G1(0)*pow(kK,2) + 
               10*kP*kq*nu*G2(0)*pow(kK,2) + 
               16*kK*kq*me2*G1(0)*pow(kP,2) - 
               32*kK*kq*me2*G2(0)*pow(kP,2) + 
               8*nu*G1(0)*pow(kK,2)*pow(kP,2) - 
               9*nu*G2(0)*pow(kK,2)*pow(kP,2) - 
               32*kK*me2*G1(0)*pow(kP,3) + 32*kK*me2*G2(0)*pow(kP,3) + 
               kK*qe2*((12*kP*G1(0) - 12*kq*G1(0) - 13*kP*G2(0) + 
                     11*kq*G2(0))*pow(kK,2) - 
                  (8*kP*G1(0) + 12*kq*G1(0) - 8*kP*G2(0) - 11*kq*G2(0))*
                   (pow(kP,2) - pow(kq,2))) + 
               64*kK*kP*me2*G1(0)*pow(kq,2) - 
               16*kK*kP*me2*G2(0)*pow(kq,2) - 
               4*nu*G2(0)*pow(kK,2)*pow(kq,2) - 
               24*nu*G1(0)*pow(kP,2)*pow(kq,2) + 
               26*nu*G2(0)*pow(kP,2)*pow(kq,2) + 
               4*kK*Mp2*((-10*kP*G1(0) + 14*kq*G1(0) + 11*kP*G2(0) - 
                     15*kq*G2(0))*pow(kK,2) + 
                  kq*(2*kP*kq*(4*G1(0) - 5*G2(0)) - 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     (-14*G1(0) + 15*G2(0))*pow(kq,2))) - 
               48*kK*me2*G1(0)*pow(kq,3) + 16*kP*nu*G1(0)*pow(kq,3) + 
               16*kK*me2*G2(0)*pow(kq,3) - 24*kP*nu*G2(0)*pow(kq,3) + 
               8*nu*G1(0)*pow(kq,4) - 2*nu*G2(0)*pow(kq,4)) + 
            Md*(12*kP*kq*nu*G1(0)*pow(kK,2) - 
               24*kK*kq*me2*G1(0)*pow(kP,2) - 
               4*nu*G1(0)*pow(kK,2)*pow(kP,2) + 
               5*nu*G2(0)*pow(kK,2)*pow(kP,2) + 
               8*kK*me2*G1(0)*pow(kP,3) - 16*kK*me2*G2(0)*pow(kP,3) - 
               64*kK*kP*me2*G1(0)*pow(kq,2) + 
               56*kK*kP*me2*G2(0)*pow(kq,2) + 
               10*nu*G1(0)*pow(kK,2)*pow(kq,2) - 
               14*nu*G2(0)*pow(kK,2)*pow(kq,2) + 
               10*nu*G1(0)*pow(kP,2)*pow(kq,2) - 
               14*nu*G2(0)*pow(kP,2)*pow(kq,2) + 
               kK*qe2*((-8*kP*G1(0) + 28*kq*G1(0) + kP*G2(0) - 
                     11*kq*G2(0))*pow(kK,2) + 
                  (kP - kq)*(3*kP*kq*(8*G1(0) - 5*G2(0)) + 
                     2*(G1(0) - 2*G2(0))*pow(kP,2) + 
                     (28*G1(0) - 11*G2(0))*pow(kq,2))) + 
               4*kK*Mp2*((kP*(6*G1(0) - 7*G2(0)) + 
                     5*kq*(-2*G1(0) + G2(0)))*pow(kK,2) + 
                  kq*(6*kP*kq*(-G1(0) + G2(0)) + 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     5*(2*G1(0) - G2(0))*pow(kq,2))) + 
               80*kK*me2*G1(0)*pow(kq,3) + 4*kP*nu*G1(0)*pow(kq,3) - 
               40*kK*me2*G2(0)*pow(kq,3) - 2*kP*nu*G2(0)*pow(kq,3) - 
               14*nu*G1(0)*pow(kq,4) + 16*nu*G2(0)*pow(kq,4))) + 
         Md*G1(qe2)*(-40*kP*kq*nu*qe2*G1(0)*pow(kK,2) + 
            14*kP*kq*nu*qe2*G2(0)*pow(kK,2) + 
            112*kK*kq*me2*qe2*G1(0)*pow(kP,2) - 
            56*kK*kq*me2*qe2*G2(0)*pow(kP,2) + 
            4*nu*qe2*G1(0)*pow(kK,2)*pow(kP,2) + 
            4*nu*qe2*G2(0)*pow(kK,2)*pow(kP,2) + 
            8*kq*me2*nu*G1(0)*pow(kP,3) - 16*kK*me2*qe2*G1(0)*pow(kP,3) - 
            4*kq*me2*nu*G2(0)*pow(kP,3) + 8*kK*me2*qe2*G2(0)*pow(kP,3) + 
            40*kK*kP*me2*qe2*G1(0)*pow(kq,2) + 
            4*kK*kP*me2*qe2*G2(0)*pow(kq,2) - 
            36*nu*qe2*G1(0)*pow(kK,2)*pow(kq,2) + 
            18*nu*qe2*G2(0)*pow(kK,2)*pow(kq,2) + 
            48*me2*nu*G1(0)*pow(kP,2)*pow(kq,2) - 
            20*nu*qe2*G1(0)*pow(kP,2)*pow(kq,2) + 
            10*nu*qe2*G2(0)*pow(kP,2)*pow(kq,2) + 
            4*Md*Mp*(kq*nu*(-((kP - kq)*kq*(4*kq*G1(0) + kP*G2(0))) + 
                  (4*kP*G1(0) - 4*kq*G1(0) - kP*G2(0) + 4*kq*G2(0))*
                   pow(kK,2)) + 
               4*kK*kP*Mp2*((G1(0) - G2(0))*pow(kK,2) - 
                  G1(0)*pow(kq,2)) + 
               kK*(-2*(kP - kq)*me2*
                   (kP*kq*(8*G1(0) + G2(0)) + 2*G2(0)*pow(kP,2) + 
                     2*(4*G1(0) - G2(0))*pow(kq,2)) - 
                  qe2*(kq*(-4*G1(0) + 5*G2(0))*pow(kP,2) + 
                     G2(0)*pow(kP,3) + 
                     2*kP*(2*G1(0) - 5*G2(0))*pow(kq,2) + 
                     4*kq*G2(0)*(-pow(kK,2) + pow(kq,2))))) + 
            24*kP*me2*nu*G1(0)*pow(kq,3) - 
            136*kK*me2*qe2*G1(0)*pow(kq,3) + 
            4*kP*nu*qe2*G1(0)*pow(kq,3) - 36*kP*me2*nu*G2(0)*pow(kq,3) + 
            44*kK*me2*qe2*G2(0)*pow(kq,3) - 8*kP*nu*qe2*G2(0)*pow(kq,3) - 
            4*Mp2*(-(kq*nu*((kP*G1(0) + 5*kq*G1(0) + 3*kP*G2(0) - 
                       3*kq*G2(0))*pow(kK,2) + 
                    kq*(-(kP*kq*(3*G1(0) + 2*G2(0))) + 
                       (-2*G1(0) + G2(0))*pow(kP,2) + 
                       (5*G1(0) - 2*G2(0))*pow(kq,2)))) + 
               kK*qe2*((-20*kq*G1(0) + 3*kP*(G1(0) - G2(0)) + 
                     7*kq*G2(0))*pow(kK,2) + 
                  kq*(9*G1(0) - 8*G2(0))*pow(kP,2) + 
                  (2*G1(0) - G2(0))*pow(kP,3) + 
                  2*kP*(-11*G1(0) + 8*G2(0))*pow(kq,2) + 
                  (20*G1(0) - 7*G2(0))*pow(kq,3)) + 
               2*kK*me2*(kq*(8*G1(0) + 3*G2(0))*pow(kP,2) + 
                  (4*G1(0) - 2*G2(0))*pow(kP,3) + 
                  3*kP*(-2*G1(0) + G2(0))*pow(kq,2) + 
                  2*(-12*G1(0) + G2(0))*pow(kq,3))) - 
            80*me2*nu*G1(0)*pow(kq,4) + 16*nu*qe2*G1(0)*pow(kq,4) + 
            40*me2*nu*G2(0)*pow(kq,4) - 2*nu*qe2*G2(0)*pow(kq,4) - 
            16*kK*kP*((G1(0) - G2(0))*pow(kK,2) + G1(0)*pow(kq,2))*
             pow(Mp,4) - 2*kK*kq*G1(0)*pow(kP,2)*pow(nu,2) + 
            kK*kq*G2(0)*pow(kP,2)*pow(nu,2) - 
            14*kK*kP*G1(0)*pow(kq,2)*pow(nu,2) + 
            7*kK*kP*G2(0)*pow(kq,2)*pow(nu,2) - 
            20*kK*G1(0)*pow(kq,3)*pow(nu,2) + 
            10*kK*G2(0)*pow(kq,3)*pow(nu,2) + 
            4*kP*G1(0)*pow(kK,3)*pow(qe2,2) - 
            40*kq*G1(0)*pow(kK,3)*pow(qe2,2) + 
            4*kP*G2(0)*pow(kK,3)*pow(qe2,2) + 
            14*kq*G2(0)*pow(kK,3)*pow(qe2,2) - 
            6*kK*kq*G1(0)*pow(kP,2)*pow(qe2,2) - 
            9*kK*kq*G2(0)*pow(kP,2)*pow(qe2,2) - 
            4*kK*G1(0)*pow(kP,3)*pow(qe2,2) + 
            8*kK*G2(0)*pow(kP,3)*pow(qe2,2) - 
            30*kK*kP*G1(0)*pow(kq,2)*pow(qe2,2) + 
            15*kK*kP*G2(0)*pow(kq,2)*pow(qe2,2) + 
            40*kK*G1(0)*pow(kq,3)*pow(qe2,2) - 
            14*kK*G2(0)*pow(kq,3)*pow(qe2,2)) + 
         Md*G2(qe2)*(6*kP*kq*nu*qe2*G1(0)*pow(kK,2) - 
            3*kP*kq*nu*qe2*G2(0)*pow(kK,2) - 
            32*kK*kq*me2*qe2*G1(0)*pow(kP,2) + 
            16*kK*kq*me2*qe2*G2(0)*pow(kP,2) + 
            6*nu*qe2*G1(0)*pow(kK,2)*pow(kP,2) - 
            3*nu*qe2*G2(0)*pow(kK,2)*pow(kP,2) - 
            4*kq*me2*nu*G1(0)*pow(kP,3) + 8*kq*me2*nu*G2(0)*pow(kP,3) - 
            12*kK*me2*qe2*G2(0)*pow(kP,3) - 
            12*kK*kP*me2*qe2*G1(0)*pow(kq,2) + 
            12*kK*kP*me2*qe2*G2(0)*pow(kq,2) + 
            24*nu*qe2*G1(0)*pow(kK,2)*pow(kq,2) - 
            12*nu*qe2*G2(0)*pow(kK,2)*pow(kq,2) + 
            4*nu*qe2*G1(0)*pow(kP,2)*pow(kq,2) - 
            24*me2*nu*G2(0)*pow(kP,2)*pow(kq,2) - 
            5*nu*qe2*G2(0)*pow(kP,2)*pow(kq,2) - 
            36*kP*me2*nu*G1(0)*pow(kq,3) + 44*kK*me2*qe2*G1(0)*pow(kq,3) + 
            16*kP*nu*qe2*G1(0)*pow(kq,3) + 24*kP*me2*nu*G2(0)*pow(kq,3) - 
            16*kK*me2*qe2*G2(0)*pow(kq,3) - 5*kP*nu*qe2*G2(0)*pow(kq,3) - 
            2*Mp2*(kq*nu*(-((kP - kq)*kq*
                     (2*kP*(G1(0) - G2(0)) + kq*(-2*G1(0) + G2(0)))) + 
                  (-6*kP*G1(0) + 12*kq*G1(0) + 7*kP*G2(0) - 7*kq*G2(0))*
                   pow(kK,2)) + 
               4*kK*me2*(kq*(3*G1(0) - 2*G2(0))*pow(kP,2) - 
                  2*(G1(0) - G2(0))*pow(kP,3) + 3*kP*G1(0)*pow(kq,2) + 
                  2*G1(0)*pow(kq,3)) + 
               kK*qe2*((3*kP*G1(0) + 12*kq*G1(0) + 4*kP*G2(0) - 
                     6*kq*G2(0))*pow(kK,2) - 
                  3*kq*(5*G1(0) - 6*G2(0))*pow(kP,2) - 
                  2*(G1(0) - G2(0))*pow(kP,3) + 
                  kP*(17*G1(0) - 19*G2(0))*pow(kq,2) + 
                  6*(-2*G1(0) + G2(0))*pow(kq,3))) + 
            40*me2*nu*G1(0)*pow(kq,4) - 20*nu*qe2*G1(0)*pow(kq,4) - 
            8*me2*nu*G2(0)*pow(kq,4) + 10*nu*qe2*G2(0)*pow(kq,4) - 
            2*Md*Mp*(2*kP*kq*nu*G1(0)*pow(kK,2) - 
               3*kP*kq*nu*G2(0)*pow(kK,2) + 
               8*kP*Mp2*(G1(0) - G2(0))*pow(kK,3) - 
               4*kK*kq*me2*G1(0)*pow(kP,2) + 
               8*kK*kq*me2*G2(0)*pow(kP,2) + 8*kK*me2*G1(0)*pow(kP,3) - 
               8*kK*me2*G2(0)*pow(kP,3) - 12*kK*kP*me2*G1(0)*pow(kq,2) - 
               8*nu*G1(0)*pow(kK,2)*pow(kq,2) + 
               9*nu*G2(0)*pow(kK,2)*pow(kq,2) + 
               2*nu*G1(0)*pow(kP,2)*pow(kq,2) - 
               2*nu*G2(0)*pow(kP,2)*pow(kq,2) + 
               8*kK*me2*G1(0)*pow(kq,3) - 2*kP*nu*G1(0)*pow(kq,3) + 
               3*kP*nu*G2(0)*pow(kq,3) - 
               kK*qe2*((kP*G1(0) + 6*kq*G1(0) - 8*kq*G2(0))*pow(kK,2) + 
                  kq*(-9*G1(0) + 10*G2(0))*pow(kP,2) - 
                  2*(G1(0) - G2(0))*pow(kP,3) + 
                  kP*(15*G1(0) - 19*G2(0))*pow(kq,2) + 
                  2*(-3*G1(0) + 4*G2(0))*pow(kq,3)) - nu*G2(0)*pow(kq,4)) \
+ 16*kP*(G1(0) - G2(0))*pow(kK,3)*pow(Mp,4) + 
            kK*kq*G1(0)*pow(kP,2)*pow(nu,2) - 
            2*kK*kq*G2(0)*pow(kP,2)*pow(nu,2) + 
            7*kK*kP*G1(0)*pow(kq,2)*pow(nu,2) + 
            kK*kP*G2(0)*pow(kq,2)*pow(nu,2) + 
            10*kK*G1(0)*pow(kq,3)*pow(nu,2) - 
            8*kK*G2(0)*pow(kq,3)*pow(nu,2) + 
            2*kP*G1(0)*pow(kK,3)*pow(qe2,2) + 
            16*kq*G1(0)*pow(kK,3)*pow(qe2,2) - 
            4*kP*G2(0)*pow(kK,3)*pow(qe2,2) - 
            5*kq*G2(0)*pow(kK,3)*pow(qe2,2) + 
            7*kK*kq*G1(0)*pow(kP,2)*pow(qe2,2) + 
            kK*kq*G2(0)*pow(kP,2)*pow(qe2,2) + 
            6*kK*G1(0)*pow(kP,3)*pow(qe2,2) - 
            9*kK*G2(0)*pow(kP,3)*pow(qe2,2) + 
            3*kK*kP*G1(0)*pow(kq,2)*pow(qe2,2) + 
            3*kK*kP*G2(0)*pow(kq,2)*pow(qe2,2) - 
            16*kK*G1(0)*pow(kq,3)*pow(qe2,2) + 
            5*kK*G2(0)*pow(kq,3)*pow(qe2,2))) - 
      F2(qp2)*(qe2*G3(qe2)*(2*Mp2*
             (8*kK*(kP - kq)*kq*me2*
                (5*kP*G1(0) - 6*kq*G1(0) - 6*kP*G2(0) + 2*kq*G2(0)) + 
               kK*qe2*((-7*kq*G1(0) + kP*(G1(0) - 5*G2(0)) + 
                     15*kq*G2(0))*pow(kK,2) + 
                  kq*(-2*kP*kq*(G1(0) - 5*G2(0)) + 
                     17*(G1(0) - G2(0))*pow(kP,2) + 
                     (7*G1(0) - 15*G2(0))*pow(kq,2))) + 
               nu*(2*(G1(0) - G2(0))*
                   (12*kP*kq + pow(kP,2) - 7*pow(kq,2))*pow(kq,2) + 
                  pow(kK,2)*(-2*kP*kq*(5*G1(0) - 6*G2(0)) + 
                     (-G1(0) + G2(0))*pow(kP,2) + 
                     2*(6*G1(0) - 5*G2(0))*pow(kq,2)))) - 
            2*Md*Mp*(4*kP*kq*nu*G1(0)*pow(kK,2) + 
               8*kP*kq*nu*G2(0)*pow(kK,2) + 
               8*kK*kq*me2*G1(0)*pow(kP,2) - 
               32*kK*kq*me2*G2(0)*pow(kP,2) - 
               nu*G1(0)*pow(kK,2)*pow(kP,2) + 
               nu*G2(0)*pow(kK,2)*pow(kP,2) - 
               88*kK*kP*me2*G1(0)*pow(kq,2) + 
               72*kK*kP*me2*G2(0)*pow(kq,2) + 
               14*nu*G1(0)*pow(kK,2)*pow(kq,2) - 
               18*nu*G2(0)*pow(kK,2)*pow(kq,2) + 
               2*nu*G1(0)*pow(kP,2)*pow(kq,2) - 
               2*nu*G2(0)*pow(kP,2)*pow(kq,2) + 
               kK*qe2*(-((3*kP*G1(0) - 19*kq*G1(0) + 5*kP*G2(0) + 
                       3*kq*G2(0))*pow(kK,2)) + 
                  kq*(6*kP*kq*G2(0) + (23*G1(0) - 13*G2(0))*pow(kP,2) + 
                     (-19*G1(0) + 3*G2(0))*pow(kq,2))) + 
               80*kK*me2*G1(0)*pow(kq,3) + 18*kP*nu*G1(0)*pow(kq,3) - 
               40*kK*me2*G2(0)*pow(kq,3) - 18*kP*nu*G2(0)*pow(kq,3) + 
               4*kK*Mp2*((kP*(G1(0) - G2(0)) - 3*kq*(G1(0) + G2(0)))*
                   pow(kK,2) + 3*(G1(0) + G2(0))*pow(kq,3)) - 
               20*nu*G1(0)*pow(kq,4) + 20*nu*G2(0)*pow(kq,4)) + 
            8*kK*(G1(0) - G2(0))*((kP + 3*kq)*pow(kK,2) - 3*pow(kq,3))*
             pow(Mp,4) + (kP - kq)*
             (nu*(kK*nu*(10*G1(0) - 11*G2(0)) - 
                  16*(kP - kq)*me2*(G1(0) - 2*G2(0)))*pow(kq,2) + 
               qe2*(8*kK*kP*(kP - kq)*me2*(G1(0) - 2*G2(0)) + 
                  nu*((4*kP*G1(0) + 8*kq*G1(0) + kP*G2(0) - 
                        13*kq*G2(0))*pow(kK,2) + 
                     kq*(kP*kq*(16*G1(0) - 23*G2(0)) + 
                        (2*G1(0) - G2(0))*pow(kP,2) + 
                        2*(G1(0) + G2(0))*pow(kq,2)))) + 
               kK*(2*kP*kq*(7*G1(0) - 8*G2(0)) - 
                  2*(G1(0) - 2*G2(0))*pow(kK,2) + 
                  (8*G1(0) - 7*G2(0))*pow(kP,2) + 
                  2*(G1(0) - 2*G2(0))*pow(kq,2))*pow(qe2,2))) + 
         Md*G2(qe2)*(Md*(kK*qe2*
                (4*kP*me2*(G1(0) - 2*G2(0))*(pow(kP,2) - pow(kq,2)) + 
                  qe2*(kq*(17*G1(0) - 22*G2(0))*pow(kP,2) + 
                     (5*G1(0) - 4*G2(0))*pow(kP,3) - 
                     6*kP*(5*G1(0) - 6*G2(0))*pow(kq,2) + 
                     2*kq*(4*G1(0) - 5*G2(0))*(-pow(kK,2) + pow(kq,2)))) \
+ 4*Mp2*(kq*(8*kK*(kP - kq)*me2*(-(kq*G1(0)) + kP*(G1(0) - G2(0))) - 
                     nu*(G1(0) - G2(0))*
                      ((kP - 3*kq)*pow(kK,2) + (-kP + kq)*pow(kq,2))) + 
                  kK*qe2*((4*kq*G1(0) + kP*(G1(0) - 2*G2(0)) - 
                        2*kq*G2(0))*pow(kK,2) + 
                     kq*(kP*kq*(7*G1(0) - 5*G2(0)) + 
                        4*(G1(0) - G2(0))*pow(kP,2) + 
                        2*(-2*G1(0) + G2(0))*pow(kq,2)))) + 
               nu*(-8*me2*(G1(0) - 2*G2(0))*pow(kP - kq,2)*pow(kq,2) + 
                  qe2*((kP - kq)*
                      (kP*G1(0) + 13*kq*G1(0) - 2*kP*G2(0) - 
                        14*kq*G2(0))*pow(kq,2) + 
                     pow(kK,2)*
                      (4*kP*kq*(G1(0) - 2*G2(0)) + 
                        (4*G1(0) - 2*G2(0))*pow(kP,2) + 
                        3*(-3*G1(0) + 4*G2(0))*pow(kq,2)))) + 
               6*kK*kP*(G1(0) - G2(0))*pow(kq,2)*pow(nu,2)) + 
            Mp*(2*nu*(kK*nu*(-10*kq*G1(0) + 8*kq*G2(0) + 
                     kP*(-8*G1(0) + G2(0))) + 
                  4*me2*G1(0)*(2*kP*kq + 5*pow(kP,2) - 7*pow(kq,2)))*
                pow(kq,2) + 4*Mp2*
                (kq*(8*kK*me2*
                      (kP*kq*(3*G1(0) - G2(0)) + 
                        (-G1(0) + G2(0))*pow(kP,2) + G1(0)*pow(kq,2)) + 
                     nu*((kP*(G1(0) - G2(0)) + 3*kq*(G1(0) + G2(0)))*
                        pow(kK,2) - (kP - kq)*(G1(0) - G2(0))*pow(kq,2))\
) + kK*qe2*((3*kP*G1(0) + 6*kq*G1(0) + 2*kP*G2(0) + 2*kq*G2(0))*
                      pow(kK,2) - 
                     kq*(kP*kq*(3*G1(0) - 5*G2(0)) + 
                        4*(G1(0) - G2(0))*pow(kP,2) + 
                        2*(3*G1(0) + G2(0))*pow(kq,2)))) + 
               qe2*(-4*kK*(kP - kq)*me2*
                   (-13*kP*kq*G1(0) + (G1(0) - 4*G2(0))*pow(kP,2) + 
                     2*(-11*G1(0) + 4*G2(0))*pow(kq,2)) + 
                  nu*(-((kP - kq)*
                        (5*kP*G1(0) + 47*kq*G1(0) - 30*kq*G2(0))*
                        pow(kq,2)) + 
                     pow(kK,2)*
                      (kP*kq*(-26*G1(0) + 18*G2(0)) + 
                        (-6*G1(0) + 2*G2(0))*pow(kP,2) + 
                        (-37*G1(0) + 14*G2(0))*pow(kq,2)))) + 
               kK*((-6*kP*(G1(0) - 2*G2(0)) + 2*kq*(-12*G1(0) + G2(0)))*
                   pow(kK,2) - 
                  (kP - kq)*(kP*kq*(46*G1(0) - 40*G2(0)) + 
                     (5*G1(0) - 8*G2(0))*pow(kP,2) + 
                     2*(12*G1(0) - G2(0))*pow(kq,2)))*pow(qe2,2))) + 
         Md*G1(qe2)*(Mp*(-4*kq*nu*
                (-(kK*kq*(4*kP + 5*kq)*nu*(2*G1(0) - G2(0))) + 
                  Mp2*((4*kP*G1(0) + 4*kq*G1(0) - kP*G2(0) + 
                        3*kq*G2(0))*pow(kK,2) + 
                     (-4*kP*G1(0) + 4*kq*G1(0) + kP*G2(0) - 7*kq*G2(0))*
                      pow(kq,2))) + 
               4*me2*(kK*(kP - kq)*qe2*
                   (kP*kq*(-34*G1(0) + 21*G2(0)) + 
                     (6*G1(0) - G2(0))*pow(kP,2) + 
                     2*(-34*G1(0) + 11*G2(0))*pow(kq,2)) + 
                  2*kq*(-((kP - kq)*kq*nu*
                        (14*kP*G1(0) + 22*kq*G1(0) - 5*kP*G2(0) - 
                         7*kq*G2(0))) + 
                     4*kK*Mp2*
                      (3*kP*kq*G2(0) + (3*G1(0) - G2(0))*pow(kP,2) + 
                        (-12*G1(0) + G2(0))*pow(kq,2)))) - 
               qe2*(4*kK*Mp2*
                   (2*(11*kq*G1(0) + kP*G2(0) - kq*G2(0))*pow(kK,2) + 
                     kq*(2*kP*kq*(5*G1(0) - 4*G2(0)) + 
                        (-7*G1(0) + 4*G2(0))*pow(kP,2) + 
                        2*(-11*G1(0) + G2(0))*pow(kq,2))) + 
                  nu*(-((kP - kq)*
                        (22*kP*G1(0) + 54*kq*G1(0) - 5*kP*G2(0) - 
                        11*kq*G2(0))*pow(kq,2)) + 
                     pow(kK,2)*
                      (kP*kq*(-64*G1(0) + 38*G2(0)) + 
                        6*G2(0)*pow(kP,2) + 
                        (-74*G1(0) + 25*G2(0))*pow(kq,2)))) + 
               kK*(kq*(30*G1(0) - 13*G2(0))*pow(kP,2) + 
                  (6*G1(0) - 5*G2(0))*pow(kP,3) + 
                  4*kq*(15*G1(0) - 4*G2(0))*(pow(kK,2) - pow(kq,2)) + 
                  2*kP*(-7*G2(0)*pow(kK,2) + 
                     (12*G1(0) + G2(0))*pow(kq,2)))*pow(qe2,2)) + 
            Md*(2*nu*(3*kK*kP*nu*G2(0) - 
                  4*(kP - kq)*me2*
                   (4*kP*G1(0) + 8*kq*G1(0) + kP*G2(0) - kq*G2(0)))*
                pow(kq,2) + qe2*
                (4*kK*kP*(kP - kq)*me2*
                   (4*kP*G1(0) + 8*kq*G1(0) + kP*G2(0) + 5*kq*G2(0)) - 
                  nu*(-((kP - kq)*
                        (12*kP*G1(0) + 8*kq*G1(0) - 3*kP*G2(0) + 
                         17*kq*G2(0))*pow(kq,2)) + 
                     3*pow(kK,2)*
                      (4*kP*kq*G1(0) + (4*G1(0) - 2*G2(0))*pow(kP,2) + 
                        (-8*G1(0) + 3*G2(0))*pow(kq,2)))) + 
               4*Mp2*(kK*qe2*(2*(kP*G1(0) + 3*kq*(-G1(0) + G2(0)))*
                      pow(kK,2) + 
                     kq*(6*kP*kq*(-G1(0) + G2(0)) - 
                        (G1(0) - 4*G2(0))*pow(kP,2) + 
                        6*(G1(0) - G2(0))*pow(kq,2))) + 
                  kq*(8*kK*me2*
                      (kP*kq*(6*G1(0) - 2*G2(0)) + 
                        (G1(0) + G2(0))*pow(kP,2) + 
                        (-4*G1(0) + G2(0))*pow(kq,2)) + 
                     nu*(-((2*kP*G1(0) + 4*kq*G1(0) + kP*G2(0) - 
                         3*kq*G2(0))*pow(kK,2)) + 
                        (2*kP*G1(0) + 4*kq*G1(0) + kP*G2(0) - kq*G2(0))*
                         pow(kq,2)))) + 
               kK*((-4*kP*G1(0) + 4*kq*G1(0) + 2*kP*G2(0) - 14*kq*G2(0))*
                   pow(kK,2) - 
                  (kP - kq)*(24*kP*kq*(G1(0) - G2(0)) + 
                     (4*G1(0) - 5*G2(0))*pow(kP,2) + 
                     2*(-2*G1(0) + 7*G2(0))*pow(kq,2)))*pow(qe2,2)))) + 
      qe2*F2(qp2)*(Md*G2(qe2)*
          (-4*pow(Mp,3)*(-(kK*kP*qe2*(4*me2 + 3*qe2)*(G1(0) - G2(0))) + 
               nu*(4*kq*me2*(kq*G1(0) + kP*(-G1(0) + G2(0))) + 
                  qe2*(kq*(kP*G1(0) - 5*kq*G1(0) - kP*G2(0) + 
                        4*kq*G2(0)) + (2*G1(0) - G2(0))*pow(kK,2))) + 
               kK*kq*(G1(0) - G2(0))*pow(nu,2)) - 
            4*Md*Mp2*(4*me2*(kq*nu*(kq*G1(0) + kP*(G1(0) - G2(0))) + 
                  kK*kP*qe2*(G1(0) - G2(0))) + 
               nu*qe2*(kP*kq*(-G1(0) + G2(0)) + G2(0)*pow(kK,2) + 
                  (3*G1(0) - 4*G2(0))*pow(kq,2)) + 
               kK*kq*(-G1(0) + G2(0))*pow(nu,2) + 
               3*kK*kP*(G1(0) - G2(0))*pow(qe2,2)) - 
            Md*nu*(G1(0) - 2*G2(0))*
             ((4*kq*(-kP + kq)*me2 + kK*(kP + kq)*nu)*qe2 + 
               pow(kq,2)*pow(nu,2) + 
               ((kP - kq)*kq + pow(kK,2))*pow(qe2,2)) + 
            Mp*nu*(3*G1(0) - 2*G2(0))*
             ((4*kq*(-kP + kq)*me2 + kK*(kP + kq)*nu)*qe2 + 
               pow(kq,2)*pow(nu,2) + ((kP - kq)*kq + pow(kK,2))*pow(qe2,2)\
)) + 2*(Md - Mp)*Mp*(G1(0) - G2(0))*G3(qe2)*
          (2*kK*kP*qe2*(8*me2*Mp2 + qe2*(2*Mp2 + qe2)) + 
            kK*(-kP + kq)*qe2*pow(nu,2) - pow(kq,2)*pow(nu,3) + 
            nu*(4*kq*((-kP + kq)*me2 + kq*Mp2)*qe2 + 
               16*me2*Mp2*pow(kq,2) + 
               (kq*(kP + kq) + pow(kK,2))*pow(qe2,2))) + 
         Md*G1(qe2)*(Mp*(4*Mp2*
                (nu*(4*kq*me2*
                      (-(kP*G1(0)) + 4*kq*G1(0) + kP*G2(0) - kq*G2(0)) \
+ qe2*(kq*(kP*G1(0) - 6*kq*G1(0) - kP*G2(0) + 3*kq*G2(0)) + 
                        (4*G1(0) - G2(0))*pow(kK,2))) + 
                  kK*kq*(G1(0) - G2(0))*pow(nu,2) + 
                  2*kK*kP*(-G1(0) + G2(0))*pow(qe2,2)) - 
               3*nu*(2*G1(0) - G2(0))*
                ((4*kq*(-kP + kq)*me2 + kK*(kP + kq)*nu)*qe2 + 
                  pow(kq,2)*pow(nu,2) + 
                  ((kP - kq)*kq + pow(kK,2))*pow(qe2,2))) + 
            Md*(-4*Mp2*(-(nu*(4*kq*me2*
                        (kP*(G1(0) - G2(0)) + kq*(2*G1(0) + G2(0))) + 
                       qe2*(kq*(-(kP*G1(0)) + kP*G2(0) - 3*kq*G2(0)) + 
                         (2*G1(0) + G2(0))*pow(kK,2)))) + 
                  kK*kq*(G1(0) - G2(0))*pow(nu,2) + 
                  2*kK*kP*(-G1(0) + G2(0))*pow(qe2,2)) - 
               3*nu*G2(0)*((4*kq*(-kP + kq)*me2 + kK*(kP + kq)*nu)*qe2 + 
                  pow(kq,2)*pow(nu,2) + 
                  ((kP - kq)*kq + pow(kK,2))*pow(qe2,2))))) + 
      2*Mp*F1(qp2)*(-2*Md*(Md - Mp)*Mp*(G1(0) - G2(0))*G2(qe2)*
          (nu*qe2*(4*kq*(-((kP - 2*kq)*me2) + 2*kq*Mp2) + 
               qe2*(kP*kq + pow(kK,2))) + kK*(2*kP + kq)*qe2*pow(nu,2) + 
            2*pow(kq,2)*pow(nu,3) + kK*kP*(4*me2 + 8*Mp2 + qe2)*pow(qe2,2)\
) + 4*Md*Mp*G1(qe2)*(Mp*(4*Mp2*
                (nu*(4*me2*G1(0)*pow(kq,2) + 
                     qe2*(G1(0)*pow(kK,2) + 
                        (-2*G1(0) + G2(0))*pow(kq,2))) + 
                  kK*kP*(-G1(0) + G2(0))*pow(qe2,2)) - 
               nu*(2*G1(0) - G2(0))*
                ((4*kq*(-kP + kq)*me2 + kK*(kP + kq)*nu)*qe2 + 
                  pow(kq,2)*pow(nu,2) + 
                  ((kP - kq)*kq + pow(kK,2))*pow(qe2,2))) + 
            Md*(4*Mp2*(nu*(qe2*G1(0)*pow(kK,2) + 
                     4*me2*G1(0)*pow(kq,2) - qe2*G2(0)*pow(kq,2)) + 
                  kK*kP*(G1(0) - G2(0))*pow(qe2,2)) - 
               nu*G2(0)*((4*kq*(-kP + kq)*me2 + kK*(kP + kq)*nu)*qe2 + 
                  pow(kq,2)*pow(nu,2) + 
                  ((kP - kq)*kq + pow(kK,2))*pow(qe2,2)))) - 
         qe2*G3(qe2)*(-(Md*(4*Mp2*
                  (kK*kP*qe2*(4*me2 + 3*qe2)*(G1(0) - G2(0)) + 
                    nu*(4*kq*me2*
                        (-2*kP*G1(0) + 4*kq*G1(0) + 2*kP*G2(0) - 
                        3*kq*G2(0)) + 
                       qe2*(-(kq*
                        (-2*kP*G1(0) + 2*kP*G2(0) + kq*G2(0))) + 
                         (3*G1(0) - 2*G2(0))*pow(kK,2))) + 
                    2*kK*kq*(G1(0) - G2(0))*pow(nu,2)) - 
                 nu*(2*G1(0) - G2(0))*
                  ((4*kq*(-kP + kq)*me2 + kK*(kP + kq)*nu)*qe2 + 
                    pow(kq,2)*pow(nu,2) + 
                    ((kP - kq)*kq + pow(kK,2))*pow(qe2,2)))) + 
            Mp*(4*Mp2*(kK*kP*qe2*(4*me2 + 3*qe2)*(G1(0) - G2(0)) + 
                  nu*(4*kq*me2*
                      (-2*kP*G1(0) + 2*kq*G1(0) + 2*kP*G2(0) - 
                        3*kq*G2(0)) + 
                     qe2*(kq*
                         (2*kP*G1(0) + 2*kq*G1(0) - 2*kP*G2(0) - 
                         kq*G2(0)) + (G1(0) - 2*G2(0))*pow(kK,2))) + 
                  2*kK*kq*(G1(0) - G2(0))*pow(nu,2)) + 
               nu*G2(0)*((4*kq*(-kP + kq)*me2 + kK*(kP + kq)*nu)*qe2 + 
                  pow(kq,2)*pow(nu,2) + 
                  ((kP - kq)*kq + pow(kK,2))*pow(qe2,2))))) + 
      F2(qp2)*(qe2*G3(qe2)*(-8*
             (kK*(4*kq*me2*(-2*kP*G1(0) + 3*kq*G1(0) + 2*kP*G2(0) - 
                     kq*G2(0)) + 
                  qe2*(2*kP*kq*(-G1(0) + G2(0)) + 
                     G2(0)*(pow(kK,2) - pow(kq,2)))) + 
               kq*nu*(G1(0) - G2(0))*(2*pow(kK,2) - 3*pow(kq,2)))*
             pow(Mp,4) - 2*Mp2*
             (nu*(8*(-kP + kq)*me2*G1(0) + kK*nu*(-5*G1(0) + 3*G2(0)))*
                pow(kq,2) - qe2*
                (4*kK*(kP - kq)*me2*
                   (5*kP*G1(0) - 6*kP*G2(0) - kq*G2(0)) + 
                  nu*(2*(2*kq*(G1(0) - G2(0)) + kP*(G1(0) + G2(0)))*
                      pow(kK,2) + 
                     kq*(2*kP*kq*(6*G1(0) - 7*G2(0)) + 
                        (G1(0) - G2(0))*pow(kP,2) + 
                        (-3*G1(0) + 5*G2(0))*pow(kq,2)))) + 
               kK*(kP*kq*(-3*G1(0) + G2(0)) + 
                  (G1(0) - 3*G2(0))*pow(kK,2) - 
                  8*(G1(0) - G2(0))*pow(kP,2) - 
                  (G1(0) - 3*G2(0))*pow(kq,2))*pow(qe2,2)) + 
            2*Md*Mp*(4*Mp2*(kK*qe2*
                   (-(kq*(2*kP*G1(0) + 4*kq*G1(0) - 2*kP*G2(0) + 
                        kq*G2(0))) + (4*G1(0) + G2(0))*pow(kK,2)) + 
                  kq*(-4*kK*me2*
                      (2*kP*G1(0) - 5*kq*G1(0) - 2*kP*G2(0) + kq*G2(0)) \
+ nu*(G1(0) - G2(0))*(2*pow(kK,2) - 3*pow(kq,2)))) + 
               nu*(8*(kP - kq)*me2*(G1(0) - G2(0)) + 
                  kK*nu*(-7*G1(0) + 3*G2(0)))*pow(kq,2) - 
               qe2*(4*kK*(kP - kq)*me2*
                   (kP*(G1(0) - 4*G2(0)) + 5*kq*(-2*G1(0) + G2(0))) + 
                  nu*(2*(5*kq*(G1(0) - G2(0)) + kP*(3*G1(0) + G2(0)))*
                      pow(kK,2) + 
                     kq*(2*kP*kq*(5*G1(0) - 6*G2(0)) + 
                        (G1(0) - G2(0))*pow(kP,2) + 
                        (-7*G1(0) + 9*G2(0))*pow(kq,2)))) + 
               kK*(kP*kq*(-5*G1(0) + 3*G2(0)) + 
                  (-7*G1(0) + 3*G2(0))*pow(kK,2) + 
                  (-8*G1(0) + 6*G2(0))*pow(kP,2) + 
                  (7*G1(0) - 3*G2(0))*pow(kq,2))*pow(qe2,2)) + 
            (kP - kq)*(nu*qe2*
                (-8*(kP - kq)*kq*me2*(G1(0) - 2*G2(0)) + 
                  qe2*(kq*(2*kP*G1(0) + 2*kq*G1(0) - 4*kP*G2(0) - 
                        kq*G2(0)) + 2*(G1(0) - 2*G2(0))*pow(kK,2))) + 
               kK*qe2*(2*kq*(G1(0) - 2*G2(0)) + kP*(-2*G1(0) + G2(0)))*
                pow(nu,2) + (-2*G1(0) + G2(0))*pow(kq,2)*pow(nu,3) + 
               kK*kP*(4*G1(0) - 5*G2(0))*pow(qe2,3))) + 
         Md*G1(qe2)*(Md*(nu*qe2*
                (-4*(kP - kq)*kq*me2*
                   (4*kP*G1(0) + 8*kq*G1(0) + kP*G2(0) + 5*kq*G2(0)) + 
                  qe2*((kP - kq)*kq*
                      (4*kP*G1(0) - 2*kP*G2(0) + 15*kq*G2(0)) + 
                     (-8*kP*G1(0) + 8*kq*G1(0) + kP*G2(0) + 5*kq*G2(0))*
                      pow(kK,2))) + 
               16*kK*G1(0)*(qe2*(pow(kK,2) - pow(kq,2)) + 
                  8*me2*pow(kq,2))*pow(Mp,4) + 
               kK*kq*qe2*(4*kP*G1(0) - 4*kq*G1(0) + 7*kP*G2(0) + 
                  8*kq*G2(0))*pow(nu,2) + 6*G2(0)*pow(kq,3)*pow(nu,3) + 
               kK*(12*kP*(kP - kq)*me2*G2(0) + 
                  qe2*(2*kP*kq*(4*G1(0) - 5*G2(0)) + 
                     (-8*G1(0) + 7*G2(0))*pow(kP,2) + 
                     3*G2(0)*(-pow(kK,2) + pow(kq,2))))*pow(qe2,2) - 
               4*Mp2*(2*nu*(kK*nu*G2(0) + 
                     4*me2*(kP*(G1(0) - G2(0)) + kq*(2*G1(0) + G2(0))))*
                   pow(kq,2) + 
                  qe2*(-4*kK*me2*
                      (2*kP*kq*G1(0) + (G1(0) + G2(0))*pow(kP,2) - 
                        G2(0)*pow(kq,2)) + 
                     nu*((3*kP*G1(0) + 6*kq*G1(0) + kq*G2(0))*
                        pow(kK,2) + 
                        (-3*kP*G1(0) - 2*kq*G1(0) + kP*G2(0) - 
                        5*kq*G2(0))*pow(kq,2))) + 
                  kK*(6*kP*kq*(G1(0) - G2(0)) + 
                     (G1(0) - 2*G2(0))*pow(kK,2) + 
                     (G1(0) - 2*G2(0))*pow(kP,2) - 
                     (G1(0) - 2*G2(0))*pow(kq,2))*pow(qe2,2))) + 
            Mp*(kq*nu*qe2*(kK*nu*
                   (38*kP*G1(0) + 40*kq*G1(0) - 17*kP*G2(0) - 
                     22*kq*G2(0)) - 
                  4*(kP - kq)*me2*
                   (14*kP*G1(0) + 34*kq*G1(0) - 5*kP*G2(0) - 13*kq*G2(0))\
) + 16*kK*G1(0)*(qe2*(pow(kK,2) - pow(kq,2)) + 8*me2*pow(kq,2))*
                pow(Mp,4) + 12*(2*G1(0) - G2(0))*pow(kq,3)*pow(nu,3) + 
               (nu*((kP - kq)*kq*
                      (8*kP*G1(0) + 42*kq*G1(0) - 2*kP*G2(0) - 
                        21*kq*G2(0)) + 
                     (26*kP*G1(0) + 46*kq*G1(0) - 17*kP*G2(0) - 
                        19*kq*G2(0))*pow(kK,2)) - 
                  12*kK*me2*(2*G1(0) - G2(0))*
                   (3*kP*kq + pow(kP,2) - 4*pow(kq,2)))*pow(qe2,2) - 
               4*Mp2*(2*nu*(kK*nu*(2*G1(0) - G2(0)) + 
                     4*me2*(-(kP*G1(0)) + 10*kq*G1(0) + kP*G2(0) - 
                        4*kq*G2(0)))*pow(kq,2) + 
                  qe2*(4*kK*me2*
                      (-2*kP*kq*G1(0) + (-3*G1(0) + G2(0))*pow(kP,2) + 
                        (20*G1(0) - 7*G2(0))*pow(kq,2)) + 
                     nu*((3*kP*G1(0) + 20*kq*G1(0) - 7*kq*G2(0))*
                        pow(kK,2) - 
                        (kP*G1(0) + 18*kq*G1(0) + kP*G2(0) - 
                        5*kq*G2(0))*pow(kq,2))) + 
                  kK*((15*G1(0) - 4*G2(0))*pow(kK,2) + 
                     (-3*G1(0) + 2*G2(0))*pow(kP,2) + 
                     (-15*G1(0) + 4*G2(0))*pow(kq,2))*pow(qe2,2)) + 
               kK*((kP - kq)*
                   (14*kP*G1(0) + 18*kq*G1(0) - 11*kP*G2(0) - 
                     9*kq*G2(0)) + 9*(2*G1(0) - G2(0))*pow(kK,2))*
                pow(qe2,3))) + 
         Md*G2(qe2)*(4*Md*Mp2*
             (kq*nu*qe2*(-((kP - 8*kq)*kq*(G1(0) - G2(0))) + 
                  (-2*G1(0) + 3*G2(0))*pow(kK,2)) + 
               4*me2*(2*(kP - kq)*nu*(G1(0) - G2(0))*pow(kq,2) + 
                  kK*qe2*(2*kP*kq*(G1(0) - G2(0)) + 
                     (G1(0) - G2(0))*pow(kP,2) + 
                     (-3*G1(0) + 2*G2(0))*pow(kq,2))) + 
               2*kK*(-G1(0) + G2(0))*pow(kq,2)*pow(nu,2) + 
               kK*kP*(2*kP*G1(0) + 9*kq*G1(0) - 2*kP*G2(0) - 8*kq*G2(0))*
                pow(qe2,2)) + 4*pow(Mp,3)*
             (2*nu*(kK*nu*(G1(0) - G2(0)) - 
                  4*me2*(kP*G1(0) - 4*kq*G1(0) - kP*G2(0) + kq*G2(0)))*
                pow(kq,2) + qe2*
                (kq*nu*(kq*(kP*G1(0) - 14*kq*G1(0) - kP*G2(0) + 
                        8*kq*G2(0)) + (10*G1(0) - 3*G2(0))*pow(kK,2)) - 
                  4*kK*me2*(2*kP*kq*(G1(0) - G2(0)) + 
                     (G1(0) - G2(0))*pow(kP,2) + 
                     (-7*G1(0) + 2*G2(0))*pow(kq,2))) + 
               kK*(-7*kP*kq*G1(0) + 8*kP*kq*G2(0) + 6*G1(0)*pow(kK,2) - 
                  2*G1(0)*pow(kP,2) + 2*G2(0)*pow(kP,2) - 
                  6*G1(0)*pow(kq,2))*pow(qe2,2)) + 
            Md*(kq*nu*qe2*(kK*nu*
                   (7*kP*G1(0) + 6*kq*G1(0) - 8*kP*G2(0) - 6*kq*G2(0)) - 
                  4*me2*(G1(0) - 2*G2(0))*(pow(kP,2) - pow(kq,2))) + 
               6*(G1(0) - G2(0))*pow(kq,3)*pow(nu,3) + 
               (4*kK*kP*(kP - kq)*me2*(G1(0) - 2*G2(0)) + 
                  nu*((5*kP*G1(0) + kq*G1(0) - 4*kP*G2(0) - 
                        2*kq*G2(0))*pow(kK,2) + 
                     3*(kP - kq)*(3*G1(0) - 4*G2(0))*pow(kq,2)))*
                pow(qe2,2) + kK*
                ((kP - kq)*(7*kP*G1(0) - kq*G1(0) - 8*kP*G2(0) + 
                     2*kq*G2(0)) - (G1(0) - 2*G2(0))*pow(kK,2))*pow(qe2,3)\
) + Mp*(kq*nu*qe2*(4*(kP - kq)*me2*
                   (5*kP*G1(0) + 13*kq*G1(0) - 4*kq*G2(0)) + 
                  kK*nu*(-17*kP*G1(0) - 22*kq*G1(0) + 6*kP*G2(0) + 
                     14*kq*G2(0))) + 
               6*(-2*G1(0) + G2(0))*pow(kq,3)*pow(nu,3) - 
               (4*kK*(kP - kq)*me2*
                   (kP*(G1(0) - 4*G2(0)) + 6*kq*(-2*G1(0) + G2(0))) + 
                  nu*((kP - kq)*kq*
                      (kq*(29*G1(0) - 22*G2(0)) + 2*kP*(G1(0) + G2(0))) \
+ (13*kP*G1(0) + 23*kq*G1(0) - 6*kP*G2(0) - 12*kq*G2(0))*pow(kK,2)))*
                pow(qe2,2) + kK*
                (-((kP - kq)*(15*kP*G1(0) + 9*kq*G1(0) - 14*kP*G2(0) - 
                       4*kq*G2(0))) + (-9*G1(0) + 4*G2(0))*pow(kK,2))*
                pow(qe2,3)))) + 
      2*Mp*F1(qp2)*(qe2*G3(qe2)*
          (Mp*(nu*(-16*(kP - kq)*kq*me2*
                   (-2*kq*G1(0) + 2*kP*(G1(0) - G2(0)) + 3*kq*G2(0)) + 
                  qe2*((kP - kq)*kq*
                      (8*kP*G1(0) + 8*kq*G1(0) - 9*kP*G2(0) - 
                        5*kq*G2(0)) + 
                     2*(2*kP*G1(0) - 2*kq*G1(0) - 3*kP*G2(0) + 
                        5*kq*G2(0))*pow(kK,2))) + 
               kK*qe2*(4*(kP - kq)*me2*
                   (4*kP*(G1(0) - G2(0)) - kq*G2(0)) + 
                  qe2*(-3*kP*kq*(4*G1(0) - 5*G2(0)) + 
                     2*(6*G1(0) - 7*G2(0))*pow(kP,2) + 
                     G2(0)*(pow(kK,2) - pow(kq,2)))) - 
               4*Mp2*(-4*kK*me2*
                   (kP*kq*G1(0) + 2*(G1(0) - G2(0))*pow(kP,2) + 
                     (-3*G1(0) + G2(0))*pow(kq,2)) + 
                  kK*qe2*(kP*kq*(-5*G1(0) + 4*G2(0)) + 
                     3*(G1(0) - G2(0))*pow(kK,2) - 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     3*(-G1(0) + G2(0))*pow(kq,2)) + 
                  nu*(2*kP*(G1(0) - G2(0))*pow(kK,2) + 
                     (-6*kP*G1(0) - 2*kq*G1(0) + 6*kP*G2(0) + kq*G2(0))*
                      pow(kq,2))) + 
               32*kK*(G1(0) - G2(0))*(pow(kK,2) - pow(kq,2))*pow(Mp,4) + 
               kK*kq*(8*kP*G1(0) - 8*kq*G1(0) - 8*kP*G2(0) + 11*kq*G2(0))*
                pow(nu,2)) + Md*
             (-(kq*nu*(-8*(kP - kq)*me2*
                     (kP*(G1(0) - 2*G2(0)) + kq*(-2*G1(0) + G2(0))) + 
                    kK*nu*(2*kP*(G1(0) - 2*G2(0)) + kq*(4*G1(0) + G2(0))))\
) + qe2*(4*kK*(kP - kq)*me2*(5*kq*(2*G1(0) - G2(0)) + 
                     2*kP*(G1(0) + G2(0))) + 
                  nu*(-((kP - kq)*kq*
                        (4*kP*G1(0) + 2*kq*G1(0) - 5*kP*G2(0) - 
                         7*kq*G2(0))) - 
                     2*(7*kP*G1(0) + kq*G1(0) - 2*kP*G2(0) - 2*kq*G2(0))*
                      pow(kK,2))) + 
               4*Mp2*(4*kK*me2*
                   (kP*kq*G1(0) - 2*(G1(0) - G2(0))*pow(kP,2) + 
                     (5*G1(0) - G2(0))*pow(kq,2)) + 
                  kK*qe2*(kP*kq*(-3*G1(0) + 4*G2(0)) + 
                     (9*G1(0) - 3*G2(0))*pow(kK,2) - 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     3*(-3*G1(0) + G2(0))*pow(kq,2)) + 
                  nu*(2*kP*(G1(0) - G2(0))*pow(kK,2) + 
                     (-6*kP*G1(0) + 6*kP*G2(0) + kq*G2(0))*pow(kq,2))) - 
               32*kK*(G1(0) - G2(0))*(pow(kK,2) - pow(kq,2))*pow(Mp,4) + 
               kK*(-((kP - kq)*
                     (10*kP*G1(0) + 10*kq*G1(0) - 8*kP*G2(0) - 
                       5*kq*G2(0))) + 5*(-2*G1(0) + G2(0))*pow(kK,2))*
                pow(qe2,2))) + 
         Md*G2(qe2)*(8*kq*me2*nu*qe2*G1(0)*pow(kP,2) - 
            4*kq*me2*nu*qe2*G2(0)*pow(kP,2) + 
            12*kP*me2*nu*qe2*G1(0)*pow(kq,2) - 
            20*me2*nu*qe2*G1(0)*pow(kq,3) + 4*me2*nu*qe2*G2(0)*pow(kq,3) - 
            16*kK*(G1(0) - G2(0))*
             (kK*kq*nu + qe2*((2*kP - kq)*kq + pow(kK,2)))*pow(Mp,4) - 
            7*kK*kP*kq*qe2*G1(0)*pow(nu,2) + 
            2*kK*kP*kq*qe2*G2(0)*pow(nu,2) - 
            kK*qe2*G1(0)*pow(kP,2)*pow(nu,2) + 
            2*kK*qe2*G2(0)*pow(kP,2)*pow(nu,2) - 
            10*kK*qe2*G1(0)*pow(kq,2)*pow(nu,2) + 
            5*kK*qe2*G2(0)*pow(kq,2)*pow(nu,2) - 
            kP*G1(0)*pow(kq,2)*pow(nu,3) + 
            2*kP*G2(0)*pow(kq,2)*pow(nu,3) - 5*G1(0)*pow(kq,3)*pow(nu,3) + 
            G2(0)*pow(kq,3)*pow(nu,3) + 20*kK*kP*kq*me2*G1(0)*pow(qe2,2) - 
            16*kK*kP*kq*me2*G2(0)*pow(qe2,2) - 
            8*kP*nu*G1(0)*pow(kK,2)*pow(qe2,2) - 
            10*kq*nu*G1(0)*pow(kK,2)*pow(qe2,2) + 
            4*kP*nu*G2(0)*pow(kK,2)*pow(qe2,2) + 
            5*kq*nu*G2(0)*pow(kK,2)*pow(qe2,2) + 
            4*kK*me2*G1(0)*pow(kP,2)*pow(qe2,2) - 
            2*kq*nu*G1(0)*pow(kP,2)*pow(qe2,2) + 
            4*kK*me2*G2(0)*pow(kP,2)*pow(qe2,2) + 
            kq*nu*G2(0)*pow(kP,2)*pow(qe2,2) - 
            24*kK*me2*G1(0)*pow(kq,2)*pow(qe2,2) - 
            7*kP*nu*G1(0)*pow(kq,2)*pow(qe2,2) + 
            12*kK*me2*G2(0)*pow(kq,2)*pow(qe2,2) + 
            5*kP*nu*G2(0)*pow(kq,2)*pow(qe2,2) + 
            9*nu*G1(0)*pow(kq,3)*pow(qe2,2) - 
            6*nu*G2(0)*pow(kq,3)*pow(qe2,2) + 
            2*Md*Mp*(2*(kP - kq)*kq*nu*
                (kK*nu*(G1(0) - G2(0)) + 
                  2*me2*(-2*kP*G1(0) + 3*kq*G1(0) + 2*kP*G2(0) - 
                     4*kq*G2(0))) + 
               8*kK*Mp2*(G1(0) - G2(0))*
                (kK*kq*nu + qe2*((2*kP - kq)*kq + pow(kK,2))) + 
               qe2*(4*kK*me2*
                   (kP*kq*G1(0) + 3*(G1(0) - G2(0))*pow(kP,2) + 
                     (-3*G1(0) + 2*G2(0))*pow(kq,2)) + 
                  nu*((kP*G1(0) - 5*kq*G1(0) - 2*kP*G2(0) + 7*kq*G2(0))*
                      pow(kK,2) + 
                     (5*kP*G1(0) - 3*kq*G1(0) - 5*kP*G2(0) + 4*kq*G2(0))*
                      pow(kq,2))) + 
               kK*(kP*kq*(-7*G1(0) + 10*G2(0)) + 
                  (-3*G1(0) + 4*G2(0))*pow(kK,2) + 
                  (7*G1(0) - 8*G2(0))*pow(kP,2) + 
                  (3*G1(0) - 4*G2(0))*pow(kq,2))*pow(qe2,2)) + 
            2*Mp2*(2*kq*nu*(-(kK*(kP - kq)*nu*(G1(0) - G2(0))) + 
                  2*me2*(kP*kq*(-7*G1(0) + 6*G2(0)) + 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     (11*G1(0) - 4*G2(0))*pow(kq,2))) + 
               qe2*(4*kK*me2*
                   (kP*kq*G1(0) - 3*(G1(0) - G2(0))*pow(kP,2) + 
                     (7*G1(0) - 2*G2(0))*pow(kq,2)) - 
                  nu*((3*kP*G1(0) - 15*kq*G1(0) - 2*kP*G2(0) + 
                        7*kq*G2(0))*pow(kK,2) + 
                     (7*kP*G1(0) + 5*kq*G1(0) - 9*kP*G2(0) + 2*kq*G2(0))*
                      pow(kq,2))) + 
               kK*(kP*kq*(9*G1(0) - 8*G2(0)) + 
                  (11*G1(0) - 4*G2(0))*pow(kK,2) + 
                  (-11*G1(0) + 12*G2(0))*pow(kP,2) + 
                  (-11*G1(0) + 4*G2(0))*pow(kq,2))*pow(qe2,2)) - 
            kK*kP*kq*G1(0)*pow(qe2,3) - kK*kP*kq*G2(0)*pow(qe2,3) - 
            6*G1(0)*pow(kK,3)*pow(qe2,3) + 3*G2(0)*pow(kK,3)*pow(qe2,3) - 
            5*kK*G1(0)*pow(kP,2)*pow(qe2,3) + 
            4*kK*G2(0)*pow(kP,2)*pow(qe2,3) + 
            6*kK*G1(0)*pow(kq,2)*pow(qe2,3) - 
            3*kK*G2(0)*pow(kq,2)*pow(qe2,3)) + 
         Md*G1(qe2)*(-8*Md*(-(kK*
                  (8*kq*(kP + kq)*me2*G1(0) + 
                    qe2*(-2*kP*kq*(G1(0) - 2*G2(0)) + 
                       3*G2(0)*(pow(kK,2) - pow(kq,2))))) - 
               2*kq*nu*((-G1(0) + G2(0))*pow(kK,2) + G1(0)*pow(kq,2)))*
             pow(Mp,3) + 8*(kK*qe2*
                (kq*(6*kP*G1(0) - 6*kq*G1(0) - 4*kP*G2(0) + 
                     3*kq*G2(0)) + (6*G1(0) - 3*G2(0))*pow(kK,2)) + 
               2*kq*(4*kK*(kP + kq)*me2*G1(0) + 
                  nu*((G1(0) - G2(0))*pow(kK,2) + G1(0)*pow(kq,2))))*
             pow(Mp,4) + 4*Md*Mp*
             (-2*kK*qe2*(qe2*(kP*kq*(-2*G1(0) + 3*G2(0)) + 
                     2*(G1(0) - G2(0))*pow(kP,2) + 
                     G2(0)*(pow(kK,2) - pow(kq,2))) + 
                  me2*G2(0)*(3*kP*kq - 4*pow(kP,2) + pow(kq,2))) - 
               2*nu*((kP - kq)*kq*me2*
                   (8*kq*G1(0) + 2*kP*G2(0) - 3*kq*G2(0)) + 
                  qe2*((2*kP*G1(0) - 2*kq*G1(0) + kq*G2(0))*pow(kK,2) + 
                     2*(-kP + kq)*G2(0)*pow(kq,2))) + 
               kK*(kP - kq)*kq*G2(0)*pow(nu,2)) - 
            4*Mp2*(kq*nu*(-(kK*(kP - kq)*nu*(2*G1(0) - G2(0))) + 
                  2*me2*(7*kP*kq*G2(0) + (4*G1(0) - 2*G2(0))*pow(kP,2) + 
                     (14*G1(0) - 11*G2(0))*pow(kq,2))) + 
               qe2*(2*kK*me2*(kP*kq*(6*G1(0) - 3*G2(0)) + 
                     (-8*G1(0) + 4*G2(0))*pow(kP,2) + 
                     (20*G1(0) - 7*G2(0))*pow(kq,2)) + 
                  nu*((8*kq*G1(0) - 6*kq*G2(0) + kP*(5*G1(0) + G2(0)))*
                      pow(kK,2) + 
                     (kq*(G1(0) - 2*G2(0)) + 5*kP*(-2*G1(0) + G2(0)))*
                      pow(kq,2))) + 
               kK*(kP*kq*(9*G1(0) - 8*G2(0)) + 
                  (13*G1(0) - 5*G2(0))*pow(kK,2) + 
                  (-5*G1(0) + 6*G2(0))*pow(kP,2) + 
                  (-13*G1(0) + 5*G2(0))*pow(kq,2))*pow(qe2,2)) + 
            (2*G1(0) - G2(0))*(nu*qe2*
                (-4*kq*me2*(kP*kq + 4*pow(kP,2) - 5*pow(kq,2)) + 
                  kK*nu*(9*kP*kq + pow(kP,2) + 8*pow(kq,2))) + 
               (kP + 5*kq)*pow(kq,2)*pow(nu,3) + 
               (nu*(2*(5*kP + 4*kq)*pow(kK,2) + 
                     kq*(kP*kq + 4*pow(kP,2) - 5*pow(kq,2))) - 
                  12*kK*me2*(kP*kq + pow(kP,2) - 2*pow(kq,2)))*pow(qe2,2) + 
               3*kK*(kP*kq + 2*pow(kK,2) + pow(kP,2) - 2*pow(kq,2))*
                pow(qe2,3)))))*pow(qp2,-1))/36.
);
}

