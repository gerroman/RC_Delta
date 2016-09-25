#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "integrand.h"

double integrate(Integrand& integrand,
    const double Ebeam, const double Q2,
    const double dE, const double dTheta, const double dPhi,
    const double Wmin2, const double Wmax2,
    const double RESULT = 0, const unsigned int ncall = 0);
#endif
