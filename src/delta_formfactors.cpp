#include "constants.h"
#include "delta_formfactors.h"
#include <cmath>

double G1(double q2) {
    return g1constant * pow(1. - q2 / pow(Lambda1, 2), -2) * pow(1. - q2 / pow(Lambda3, 2), -1);
}
double G2(double q2) {
    return g2constant * pow(1. - q2 / pow(Lambda1, 2), -2) * pow(1. - q2 / pow(Lambda3, 2), -1);
}
double G3(double q2) {
    return g3constant * pow(1. - q2 / pow(Lambda1, 2), -2) * pow(1. - q2 / pow(Lambda3, 2), -1) * (
        aValue * pow(1. - q2 / pow(Lambda2, 2), -1) + (1.-aValue) * pow(1. - q2 / pow(Lambda4, 2), -1)
        );
}

double GMstar(double q2) {
    return Mp/(3.0 * (Md + Mp)) * ( (pow(Md + Mp, 2) - q2)/Md2 * G1(q2) - (Md2 - Mp2 + q2)/(2.0*Md2) * (G1(q2) - G2(q2)) + q2/Md2 * G3(q2) );
}
double GEstar(double q2) {
    return Mp/(3.0 * (Md + Mp)) * (- (Md2 - Mp2 + q2)/(2.0*Md2) * (G1(q2) - G2(q2)) + q2/Md2 * G3(q2) );
}
double GCstar(double q2) {
    return 2.0*Mp/(3.0 * (Md + Mp)) * ( -G1(q2) + G2(q2) + (Md2 - Mp2 + q2)/(2.0*Md2) * G3(q2) );
}
