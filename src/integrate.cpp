#include "constants.h"
#include "born.h"
#include "integrand.h"

#include "Math/IntegratorMultiDim.h"

#include <iostream>
#include <cmath>

double integrate(Integrand& integrand,
    const double Ebeam, const double Q2,
    const double dE, const double dTheta, const double dPhi,
    const double Wmin2, const double Wmax2,
    const double RESULT = 0, const unsigned int ncall = 0) {

  std::cerr << "[Integrate] " 
            << " Ebeam = " << Ebeam
            << " Q2 = "    << Q2
            << " Wmin2 = " << Wmin2
            << " Wmax2 = " << Wmax2
            << " dTheta = "<< dTheta
            << " dPhi = "  << dPhi
            << " ncall = " << ncall
            << std::endl;
  double p[5] = {Ebeam, Q2, dE, dTheta, dPhi};
  integrand.SetParameters(p);

  double a[3] = {Wmin2, 0, -M_PI};
  double b[3] = {Wmax2, M_PI, M_PI};

  ROOT::Math::IntegratorMultiDim ig(integrand, ROOT::Math::IntegrationMultiDim::kVEGAS, 0, 0, ncall); 
  double val = ig.Integral(a, b) / (born(Ebeam, Q2) * 4 * M_PI * M_PI);
  double error = ig.Error() / (born(Ebeam, Q2) * 4 * M_PI * M_PI);

  if (fabs(RESULT) >  0) {
    std::cerr << "expect = " << RESULT << std::endl;
  }
  std::cerr << "result = " << val << " +/- " << error << std::endl;

  return val;
}
