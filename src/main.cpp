#include "constants.h"
#include "iMdummy2.h"
#include "iMd2.h"
#include "iMp2.h"
#include "iMe2.h"
#include "iMpiMe.h"
#include "iMdiMe.h"
#include "iMdiMp.h"
#include "iMdiMeSoft.h"
#include "iMdiMeHard.h"
#include "born.h"
#include "integrand.h"
#include "integrate.h"

#include <cmath>
#include <iostream>
#include <fstream>

int calculate_1() {
  const double Ebeam = 1.594;
  const double Q2 = 1.;
  const double MeV = 1e-3;
  const double wmin = 0 * MeV;
  const double wmax = 100 * MeV;
  const double Wmin2 = pow(wmin + sqrt(wmin*wmin + Mp2), 2);
  const double Wmax2 = pow(wmax + sqrt(wmax*wmax + Mp2), 2);
  const double RESULT = (
      2*Mp2*(2*Ebeam + Mp)*log(Mp2/Wmax2) - (Mp2-Wmax2)*(4*Ebeam + 3 * Mp - Wmax2/Mp)
      ) / (16*Ebeam) / (born(Ebeam, Q2) * 4 * M_PI * M_PI);
  const double dE = 0.25;
  const double dTheta = M_PI;
  const double dPhi = 2*M_PI;
  Integrand integrand(&iMdummy2);
  integrate(integrand, Ebeam, Q2, dE, dTheta, dPhi, Wmin2, Wmax2, RESULT);
  return 0;
}

int calculate_2() {
  const double Ebeam = 1.594;
  const double Q2 = 1.;
  const double MeV = 1e-3;
  const double wmin = 1 * MeV;
  const double wmax = 100 * MeV;
  const double Wmin2 = pow(wmin + sqrt(wmin*wmin + Mp2), 2);
  const double Wmax2 = pow(wmax + sqrt(wmax*wmax + Mp2), 2);
  const double dE = 0.25;
  const double dTheta = M_PI;
  const double dPhi = 2*M_PI;
  const double eta = (2*Ebeam*Mp) / (2 * Ebeam * Mp - Q2);
  const double RESULT = Z * e2 / (M_PI*M_PI) * log(wmax/wmin) * log(eta);
  Integrand integrand(&iMpiMe);
  integrate(integrand, Ebeam, Q2, dE, dTheta, dPhi, Wmin2, Wmax2, RESULT);
  return 0;
}

int calculate_3() {
  const double Ebeam = 1.594;
  const double Q2 = 1.;
  const double MeV = 1e-3;
  const double wmin = 0.01 * MeV;
  const double wmax = 1 * MeV;
  const double Wmin2 = pow(wmin + sqrt(wmin*wmin + Mp2), 2);
  const double Wmax2 = pow(wmax + sqrt(wmax*wmax + Mp2), 2);
  const double dE = 0.25;
  const double dTheta = M_PI;
  const double dPhi = 2*M_PI;
  const double eta = (2*Ebeam*Mp) / (2 * Ebeam * Mp - Q2);
  const double RESULT = Z * e2 / (M_PI*M_PI) * log(wmax/wmin) * log(eta);
  Integrand integrand(&iMpiMe);
  integrate(integrand, Ebeam, Q2, dE, dTheta, dPhi, Wmin2, Wmax2, RESULT);
  return 0;
}

int calculate_4(const double dTheta, const double dPhi, const char* filename) {
  std::ofstream output;
  output.open(filename);
  const double Ebeam = 1.594;
  const double Q2 = 1.51;
  const double Wmin2 = Mp2;
  const double dE = 0.;
  const unsigned int ncall = 10000; // sufficient for 1% accuracy

  const double Wmax2l = 1.0;
  const double Wmax2r = 2.4;
  const int nPoints = 50;
  const double dWmax2 = (Wmax2r - Wmax2l) / nPoints;

  Integrand integrand11(&iMd211);
  Integrand integrand12(&iMd212);
  Integrand integrand22(&iMd222);

  for (double Wmax2 = Wmax2l; Wmax2 <= Wmax2r; Wmax2+=dWmax2) {
    double v1 = integrate(integrand11, Ebeam, Q2, dE, dTheta, dPhi, Wmin2, Wmax2, 0, ncall);
    double v2 = integrate(integrand12, Ebeam, Q2, dE, dTheta, dPhi, Wmin2, Wmax2, 0, ncall);
    double v3 = integrate(integrand22, Ebeam, Q2, dE, dTheta, dPhi, Wmin2, Wmax2, 0, ncall);
    output << Wmax2 << "\t" << v1 << "\t" << v2 << "\t" << v3 << std::endl;
  }
  output.close();
  return 0;
}

int calculate_5(const double Ebeam, const double Q2, const double dE, const double dTheta, const double dPhi, const unsigned int ncall) {
  const double Wmin2 = Mp2;
  const double Wmax2 = Mp2 + 2 * Mp * Ebeam * dE;
  const double RESULT = 0;
  std::cout << "Ebeam = " << Ebeam << "    Q2 = " << Q2 << std::endl;
  std::cout << "ncall = " << ncall << std::endl;

  std::cout << " == Full == " << std::endl;
  Integrand integrand1(&iMdiMe);
  double delta_full = integrate(integrand1, Ebeam, Q2, dE, dTheta, dPhi, Wmin2, Wmax2, RESULT, ncall);

  std::cout << " == (s1) == " << std::endl;
  Integrand integrand2(&iMdiMeSoft1);
  double delta_s1 = integrate(integrand2, Ebeam, Q2, dE, dTheta, dPhi, Wmin2, Wmax2, RESULT, ncall);

  std::cout << " == (h1) == " << std::endl;
  Integrand integrand3(&iMdiMeHard1);
  double delta_h1 = integrate(integrand3, Ebeam, Q2, dE, dTheta, dPhi, Wmin2, Wmax2, RESULT, ncall);

  std::cout << " == (s2) == " << std::endl;
  Integrand integrand4(&iMdiMeSoft2);
  double delta_s2 = integrate(integrand4, Ebeam, Q2, dE, dTheta, dPhi, Wmin2, Wmax2, RESULT, ncall);

  std::cout << " == (h2) == " << std::endl;
  Integrand integrand5(&iMdiMeHard2);
  double delta_h2 = integrate(integrand5, Ebeam, Q2, dE, dTheta, dPhi, Wmin2, Wmax2, RESULT, ncall);

  std::cout << " == CHECK == " << std::endl;
  std::cout << " delta_full = " <<  delta_full  << std::endl
            << " delta_sum  = " << delta_s1 + delta_s2 + delta_h1 + delta_h2 << std::endl;
  std::cout << " =========== " << std::endl;
  return 0;
}

int main() {
  std::cout << "1 - Self test 1: sunstitute matrix element to 1 and integrate without angle cuts" << std::endl;
  std::cout << "2 - Selt test 2: integrate interference 2 Re[Me Mp^+] and compare with Maximon&Tjon" << std::endl;
  std::cout << "3 - integrate |Md|^2, tabulate and save results without and with angle cuts to output_{1,2}.txt" << std::endl;
  std::cout << "4 - integrate 2 Re[Md (Me)^+] for VEPP-3 experiment" << std::endl;
  std::cout << "Choose mode: ";
  char c;
  std::cin >> c;
  switch (c) {
    case '1':
      calculate_1();
      break;
    case '2':
      calculate_2();
      calculate_3();
      break;
    case '3':
      calculate_4(180 * Deg, 360 * Deg, "output_1.txt");
      calculate_4(3 * Deg, 3 * Deg, "output_2.txt");
      break;
    case '4':
      char c2;
      std::cout << "1 - Run I,  No.1: Ebeam = 1.594, Q2 = 1.51 , dE/E3 = 0.25, dTheta = dPhi = 3 Degree" << std::endl;
      std::cout << "2 - Run I,  No.2: Ebeam = 1.594, Q2 = 0.298, dE/E3 = 0.45, dTheta = dPhi = 5 Degree" << std::endl;
      std::cout << "3 - Run II, No.1: Ebeam = 0.998, Q2 = 0.976, dE/E3 = 0.29, dTheta = dPhi = 3 Degree" << std::endl;
      std::cout << "4 - Run II, No.2: Ebeam = 0.998, Q2 = 0.840, dE/E3 = 0.29, dTheta = dPhi = 3 Degree" << std::endl;
      std::cout << "Choose point: ";
      std::cin >> c2;
      switch (c2) {
        case '1':
          calculate_5(1.594, 1.51,  0.25, 3.0 * Deg, 3.0 * Deg, 1E+7);
          break;
        case '2':
          calculate_5(1.594, 0.298, 0.45, 5.0 * Deg, 5.0 * Deg, 1E+7);
          break;
        case '3':
          calculate_5(0.998, 0.976, 0.29, 3.0 * Deg, 3.0 * Deg, 1E+7);
          break;
        case '4':
          calculate_5(0.998, 0.830, 0.29, 3.0 * Deg, 3.0 * Deg, 1E+7);
          break;
        default:
          break;
      }
    default:
      std::cout << "Wrong choice ... " << std::endl;
  }
  return 0;
}
