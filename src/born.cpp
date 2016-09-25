#include "born.h"
#include "constants.h"
#include "kinematic.h"
#include "iMborn2.h"
#include <TLorentzVector.h>

double born(const double Ebeam, const double Q2) {
  const TLorentzVector p1 = TLorentzVector(0, 0, sqrt(Ebeam * Ebeam - me2), Ebeam);
  const TLorentzVector p2 = TLorentzVector(0, 0, 0, Mp);
  /* theta, assuming me << E1, E3el */
  const double theta = 2 * asin(sqrt(
        Q2 * Mp / (2 * Ebeam * (2 * Ebeam * Mp - Q2))
        ));
  /* eta, assuming me << E1, E3el */
  const double eta = (2 * Ebeam * Mp) / (2 * Ebeam * Mp - Q2);
  const double E3el = getE3(Ebeam, Mp2, theta, eta);
  const double p3vel = sqrt(E3el * E3el - me2);
  const TLorentzVector p3el(p3vel * sin(theta), 0, p3vel * cos(theta), E3el);
  const TLorentzVector p4el = p1 + p2 - p3el;
  return iMborn2((p1+p3el)*(p2+p4el), -Q2);
}
