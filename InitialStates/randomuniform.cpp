#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include "../Math/random.h"
#include "../particle.h"
#include "../system.h"
#include <cmath>
#include <random>
#include <vector>

using std::cout;
using std::endl;

using namespace std;

RandomUniform::RandomUniform(System* system, int numberOfHiddenNodes, int numberOfVisibleNodes, double sigma,vector<double>& m_X, vector<double>& m_Hidden, vector<double>& m_a_bias, vector<double>& m_b_bias, vector<std::vector<double>>& m_w,  double interactionSize, double timeStep, int bins, double bucketSize)   :
        InitialState(system) {
    assert(numberOfHiddenNodes > 0 && numberOfVisibleNodes > 0);


    setSigma(sigma);
    setSigma_squared(sigma*sigma);

   setNumberOfHiddenNodes(numberOfHiddenNodes);
   setNumberOfVisibleNodes(numberOfVisibleNodes);

    //m_system->setNumberOfParticles(numberOfParticles);
    m_system->setinteractionSize(interactionSize);
    m_system->setTimeStep(timeStep);
    m_system->setSqrtTimeStep(sqrt(timeStep));
    m_system->setBins(bins);
    m_system->setBucketSize(bucketSize);

    setupInitialState(m_X, m_Hidden, m_a_bias, m_b_bias, m_w); //delete m_hidden
}

void RandomUniform::setupInitialState(vector<double>& m_X, vector<double>& m_Hidden, vector<double>& m_a_bias, vector<double>& m_b_bias, vector<std::vector<double>>& m_w) {

   double sigma_0=0.001;

    for (int i=0; i < m_numberOfVisibleNodes; i++) {
     m_X.push_back(Random::nextDouble()-0.5);
    }

    for(int i=0; i<m_numberOfVisibleNodes; i++){
        m_a_bias.push_back(Random::nextGaussian(0,sigma_0));
    }

    for(i=0; i<m_numberOfHiddenNodes; i++){
        m_b_bias.push_back(Random::nextGaussian(0,sigma_0));
    }

    for(int i=0; i<m_numberOfVisibleNodes; i++){
        for(j=0;j<m_numberOfHiddenNodes; j++){
            m_w.push_back(Random::nextGaussian(0,sigma_0));
        }
    }


}






/*
void RandomUniform::setupInitialStateWithInteraction() {
    int placedParticles = 0;
    double R_ki = 0;
    vector<vector<double>> distancematrix(m_numberOfParticles, vector<double>(m_numberOfParticles));

    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();
        for (int d=0; d < m_numberOfDimensions; d++) {
            position.push_back(Random::nextDouble()-0.5);
        }
        for (int k = 0; k < m_numberOfParticles; k++){
            for (int d=0; d < m_numberOfDimensions; d++) {
                R_ki += (m_particles.at(k).getPosition()[d] - m_particles.at(i).getPosition()[d]) *
                        (m_particles.at(k).getPosition()[d] - m_particles.at(i).getPosition()[d]);
            }
            R_ki = sqrt(R_ki);
            //if (R_ki < m_system->getinteractionSize()){
                Particle p;
                m_particles.push_back(p);
                m_particles.at(i).setNumberOfDimensions(m_numberOfDimensions);
                m_particles.at(i).setPosition(position);

                distancematrix[k][i] = R_ki;
                distancematrix[i][k] = R_ki;
                placedParticles++;
          //  }
        }
    }
    m_system->setDistanceMatrix(distancematrix);
}
*/
