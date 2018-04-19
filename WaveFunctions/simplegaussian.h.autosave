#pragma once
#include "wavefunction.h"
#include "vector"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double numberOfHiddenNodes, double numberOfVisibleNodes, std::vector<double> X, std::vector<double> Hidden, std::vector<double> a_bias, std::vector<double> b_bias, std::vector<std::vector<double>> w);
    double evaluate(std::vector<double> X, std::vector<double> Hidden, std::vector<double> a_bias, std::vector<double> b_bias, std::vector<std::vector<double> > w)=0;

    double computeDoubleDerivative(std::vector<double> X, std::vector<double> Hidden, std::vector<double> a_bias, std::vector<double> b_bias, std::vector<std::vector<double> > w);
    std::vector<std::vector<double>> QuantumForce              (std::vector<class Particle>& particles);

protected:
    class WaveFunction* m_wavefunction = nullptr;
};
