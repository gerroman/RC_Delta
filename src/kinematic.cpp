#include "kinematic.h"
#include "constants.h"

#include <TLorentzVector.h>

#include <iostream>
#include <iomanip>

std::ostream& operator<<(std::ostream& os, const TLorentzVector& p) {
  std::setprecision(6);
  os  << std::setw(10) << p.M()  << "\t("
      << std::setw(10) << p.E()  << ",\t"
      << std::setw(10) << p.Px() << ",\t"
      << std::setw(10) << p.Py() << ",\t"
      << std::setw(10) << p.Pz() << ")";
  return os;
}

void printKinematic(const TLorentzVector& p1, const TLorentzVector& p2,
    const TLorentzVector& p3, const TLorentzVector& p4, const TLorentzVector& k) {
  std::cout << p1 << std::endl
            << p2 << std::endl
            << p3 << std::endl
            << p4 << std::endl
            << k  << std::endl
            << (p1 + p2 - p3 - p4 - k) << std::endl;
}

double getE3(double Ebeam, double W2, double theta, double eta) {
  double x = W2 - Mp2;
  double e3approx = (Ebeam - x/(2*Mp))/eta;
  double cosTheta2 = pow(cos(theta), 2);
  double disc = sqrt(cosTheta2*(Ebeam - me)*(Ebeam + me)*(4*(Ebeam - me)*(Ebeam + me)*(
          (-1 + cosTheta2)*me2 + Mp2) - 4*(me2 + Ebeam*Mp)*x + x*x));
  double b = (Ebeam + Mp)*(2*me2 + 2*Ebeam*Mp - x);
  double c = (2*(1 - cosTheta2)*Ebeam*Ebeam + 4*Ebeam*Mp + 2*(cosTheta2*me2 + Mp2));
  double e3solution1 = (b - disc)/c;
  double e3solution2 = (b + disc)/c;
  if (fabs(e3approx - e3solution1) < fabs(e3approx - e3solution2)) {
    /* std::cerr << "Another solution" << std::endl; */
    return e3solution1;
  } else {
    return e3solution2;
  }
}
