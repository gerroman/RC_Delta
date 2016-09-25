#include "integrand.h"

#include "constants.h"
#include "kinematic.h"
#include "iMd2.h"
#include <TLorentzVector.h>
#include <cmath>

double evaluate(double W2, double ThetaG, double PhiG, const double Ebeam,
    const double Q2, const double /*dE*/, const double dTheta, const double dPhi,
    double (*func)(double nu, double qe2, double qp2, double kK, double kP, double kq)) {
  static const TLorentzVector p1 = TLorentzVector(0, 0, sqrt(Ebeam*Ebeam - me2), Ebeam);
  static const TLorentzVector p2 = TLorentzVector(0, 0, 0, Mp);

  static const double eta = (2*Ebeam*Mp) / (2 * Ebeam * Mp - Q2); // me << E1, E3el
  static const double Theta3 = 2 * asin(sqrt( Q2 * Mp / (2 * Ebeam * (2 * Ebeam * Mp - Q2)))); // me << E1, E3el
  static const double E3el = getE3(Ebeam, Mp2, Theta3, eta);
  static const double p3vel = sqrt(E3el*E3el - me2);
  static const TLorentzVector p3el(p3vel * sin(Theta3), 0, p3vel * cos(Theta3), E3el);

  static const TLorentzVector p4el = p1 + p2 - p3el;
  static const double Theta4 = p4el.Theta();
  static const double Phi4 = p4el.Phi();

  /* E3, assuming me << E1, E3 */
  /* double E3 = p1.E() / eta - (W2 - Mp2) / (2 * Mp * eta); */
  /* E3, exact */
  double E3 = getE3(Ebeam, W2, Theta3, eta);
  double p3v = sqrt(E3 * E3 - me2); 
  TLorentzVector p3(p3v * sin(Theta3), 0, p3v * cos(Theta3), E3); 

  double omega =  (W2 - Mp2)/(2 * sqrt(W2));
  TLorentzVector k(
      omega * sin(ThetaG) * cos(PhiG), 
      omega * sin(ThetaG) * sin(PhiG), 
      omega * cos(ThetaG), 
      omega
      );
  TLorentzVector p4(
      -omega * sin(ThetaG) * cos(PhiG), 
      -omega * sin(ThetaG) * sin(PhiG), 
      -omega * cos(ThetaG), 
      sqrt(Mp2 + omega * omega)
      );

  TLorentzVector t1 = p1 + p2 - p3;
  double tv = sqrt(t1.E() * t1.E() - W2);
  double ThetaW = atan2(t1.Px(), t1.Pz());
  k.Boost(0, 0, tv/t1.E());
  k.RotateY(ThetaW);
  p4.Boost(0, 0, tv/t1.E());
  p4.RotateY(ThetaW);

  if ((fabs(p4.Theta() - Theta4) < dTheta) && (Phi4 - fabs(p4.Phi()) < dPhi)) {
    double nu = p1 * p2 + p3 * p2 + p1 * p4 + p3 * p4; // (p1 + p3) * (p2 + p4)
    double qe2 = 2*me2 - 2*(p1*p3);
    double qp2 = 2*Mp2 - 2*(p2*p4);
    double kK = k*p1 + k*p3;
    double kP = k*p2 + k*p4;
    double kq = (qe2 - qp2) / 2;
    return (
        (1 - (W2 - Mp2) / (2*Mp*Ebeam)) * (W2 - Mp2) / (4 * W2) / (4*M_PI) * sin(ThetaG) *
        func(nu, qe2, qp2, kK, kP, kq)
      );
  } else {
    return 0;
  }
}
