#ifndef KINEMATIC_H
#define KINEMATIC_H

#include <TLorentzVector.h>
#include <iostream>

std::ostream& operator<<(std::ostream& os, const TLorentzVector& p);

void printKinematic(const TLorentzVector& p1, const TLorentzVector& p2,
    const TLorentzVector& p3, const TLorentzVector& p4, const TLorentzVector& k);

double getE3(double Ebeam, double W2, double theta, double eta);

#endif
