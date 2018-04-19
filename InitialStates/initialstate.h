#pragma once
#include <vector>
#include "../particle.h"

class InitialState {
public:
    InitialState(class System* system);
    virtual void setupInitialState() = 0;
    std::vector<class Particle> getParticles() { return m_particles; }

    double getSigma() const;
    void setSigma(double sigma);

    double getSigma_squared() const;
    void setSigma_squared(double sigma_squared);

    int getNumberOfVisibleNodes() const;
    void setNumberOfVisibleNodes(int numberOfVisibleNodes);

    int getNumberOfHiddenNodes() const;

    void setNumberOfHiddenNodes(int numberOfHiddenNodes);

protected:
    class System* m_system = nullptr;
    class Random* m_random = nullptr;
    std::vector<Particle> m_particles; //= std::vector<Particle>();
    int m_numberOfDimensions = 0;
    int m_numberOfParticles = 0;
    int m_numberOfVisibleNodes =0;
    int m_numberOfHiddenNodes=0;
    double m_sigma =0;
    double m_sigma_squared =0;
    std::vector<double> m_X=std::vector<double>();
    std::vector<double> m_Hidden=std::vector<double>();
    std::vector<double> m_a_bias=std::vector<double>();
    std::vector<double> m_b_bias=std::vector<double>();
    std::vector<std::vector<double>> m_w=std::vector<std::vector<double>>;

};

