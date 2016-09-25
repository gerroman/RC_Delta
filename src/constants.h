#ifndef CONSTANTSH
#define CONSTANTSH

#include <cmath>

const double me  = 0.51099893e-3;
const double Mp  = 0.938272046;
const double e2  = 4 * M_PI / 137.036;
const double mu = 2.79284736;
const double me2 = pow(me, 2);
const double me4 = pow(me, 4);
const double Mp2 = pow(Mp, 2);
const double Mp4 = pow(Mp, 4);
const double e4 = pow(e2, 2);
const double e6 = pow(e2, 3);
const double Z   = 1.;
const double Z2  = 1.;
const double Z3  = 1.;
const double Z4  = 1.;

const double Deg = M_PI / 180.0;

//parametrization of the proton form factors (dipole form)
const double a11 = 0.; // Electric form factor
const double b11 = 9.92;
const double b12 = 24.6;
const double b13 = 0.;

const double a21 = 0.; // Magnetic form factor
const double b21 = 9.92;
const double b22 = 24.6;
const double b23 = 0.;

//parametrization of the delta form factors (Zhou & Yang):
const double g1constant = 6.59;
const double g2constant = 9.08;
const double g3constant = 7.12;
const double Lambda1 = sqrt(0.71);
const double Lambda2 = 2.;
const double Lambda3 = sqrt(2.);
const double Lambda4 = 0.2;
const double aValue = -0.3;
const double Md = 1.232;
const double Md2 = pow(Md, 2);
const double Gd = 0.117;
const double Gd2 = pow(Gd, 2);

#endif
