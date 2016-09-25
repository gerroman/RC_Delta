#ifndef INTEGRAND_H
#define INTEGRAND_H

#include <Math/IParamFunction.h>

double evaluate(double W2, double ThetaG, double PhiG, const double Ebeam,
    const double Q2, const double /*dE*/, const double dTheta, const double dPhi,
    double (*func)(double nu, double qe2, double qp2, double kK, double kP, double kq));

class Integrand: public ROOT::Math::IParametricFunctionMultiDim
{
  private:
    const double* pars;
    double (*func)(double nu, double qe2, double qp2, double kK, double kP, double kq);
  public:
    Integrand(double (*func)(double nu, double qe2, double qp2, double kK, double kP, double kq)) {
        this->func = func;
    };
    double DoEvalPar(const double* x, const double* p) const {
      return evaluate(x[0], x[1], x[2], p[0], p[1], p[2], p[3], p[4], func);
    }
    unsigned int NDim() const{
      return 3;
    }
    ROOT::Math::IParametricFunctionMultiDim* Clone() const{
      return new Integrand(func);
    }
    const double* Parameters() const{
      return pars;
    }
    void SetParameters(const double* p){
      pars = p;
    }
    unsigned int NPar() const{
      return 5;
    }
};


#endif
