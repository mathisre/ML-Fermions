#pragma once
#include "hamiltonian.h"
#include <vector>
using namespace std;

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega, double omega_z);
    double computeLocalEnergy(bool interaction, vector<double> X, vector<double> Hidden, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double> > w);

    std::vector<double> omega() const;
    void setOmega(const std::vector<double> &omega);

    double computeInteractionPotential();
    double computePotentialEnergy(vector<double> X);
private:
    std::vector<double> m_omega = std::vector<double>();
};

