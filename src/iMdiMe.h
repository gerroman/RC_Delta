#ifndef IMDIME_H
#define IMDIME_H

double iMdiMe1(double nu, double qe2, double qp2, double kK, double kP, double kq);
double iMdiMe2(double nu, double qe2, double qp2, double kK, double kP, double kq);
double iMdiMe(double nu, double qe2, double qp2, double kK, double kP, double kq) {
  return (iMdiMe1(nu, qe2, qp2, kK, kP, kq) + iMdiMe2(nu, qe2, qp2, kK, kP, kq));
}

#endif
