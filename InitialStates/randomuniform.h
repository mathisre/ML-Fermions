#pragma once
#include "initialstate.h"

class RandomUniform : public InitialState {
public:
    RandomUniform(System* system, int numberOfDimensions, int numberOfParticles, double sigma, vector<double> m_X, vector<double> m_Hidden, vector<double> m_a_bias, vector<double> m_b_bias, vector<std::vector<double>> m_w, double interactionSize, double timeStep, int bins, double bucketSize);
    void setupInitialState(vector<double> &m_X, vector<double> &m_Hidden, vector<double> &m_a_bias, vector<double> &m_b_bias, vector<std::vector<double> > &m_w);
  //  void setupInitialStateWithInteraction();

};

