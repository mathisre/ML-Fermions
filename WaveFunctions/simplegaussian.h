#pragma once
#include "wavefunction.h"
#include "vector"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system);//, std::vector<double> a_bias, std::vector<double> b_bias, std::vector<std::vector<double>> w);
    double evaluate(double GibbsValue, std::vector<double> X, std::vector<double> Hidden, std::vector<double> a_bias, std::vector<double> b_bias, std::vector<std::vector<double> > w);

    double computeDoubleDerivative(double GibbsValue, std::vector<double> X, std::vector<double> Hidden, std::vector<double> a_bias, std::vector<double> b_bias, std::vector<std::vector<double> > w);
    //std::vector<std::vector<double>> QuantumForce              (std::vector<class Particle>& particles);

    std::vector<double> QuantumForce(double GibbsValue, std::vector<double> X, std::vector<double> a_bias, std::vector<double> b_bias, std::vector<std::vector<double> > w);
protected:
    class WaveFunction* m_wavefunction = nullptr;
};
