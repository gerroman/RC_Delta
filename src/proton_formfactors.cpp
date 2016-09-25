#include "constants.h"
#include "proton_formfactors.h"
#include <cmath>

double GE(double q2) 
{ 
    double tau = (-q2)/(4.*Mp2);
    return (1. + a11*tau)/(1. + b11*tau + b12*pow(tau, 2) + b13*pow(tau, 3));
}
double GM(double q2)
{
    double tau = (-q2)/(4.*Mp2);
    return mu*(1. + a21*tau)/(1. + b21*tau + b22*pow(tau, 2) + b23*pow(tau, 3));
}

double F1(double q2)
  {
  double tau = (-q2)/(4.*Mp2);
  return (GE(q2) + tau*GM(q2))/(1. + tau);
  }
double F2(double q2)
{
    double tau = (-q2)/(4.*Mp2);
    return (GM(q2) - GE(q2))/(1. + tau);
}
