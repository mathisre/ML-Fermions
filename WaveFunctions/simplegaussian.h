#pragma once
#include "wavefunction.h"
#include "vector"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha, double beta);
    double evaluate(std::vector<double>& m_X);

    double computeDoubleDerivative(std::vector<class Particle>& particles);
    std::vector<std::vector<double>> QuantumForce              (std::vector<class Particle>& particles);

protected:
    class WaveFunction* m_wavefunction = nullptr;
};
