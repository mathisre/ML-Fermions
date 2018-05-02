#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include <vector>

using namespace std;
using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega, double omega_z) :
        Hamiltonian(system) {
    assert(omega > 0);
    assert(omega_z > 0);
    m_omega.reserve(3);
    m_omega.push_back(omega);
    m_omega.push_back(omega);
    m_omega.push_back(omega_z);

}

double HarmonicOscillator::computeLocalEnergy(bool interaction, vector<double> X, vector<double> Hidden, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w) {

    double potentialEnergy = 0;
    double kineticEnergy   = 0;
    double interaction_potential=0;

    kineticEnergy = -0.5*m_system->getWaveFunction()->computeDoubleDerivative(X,Hidden,a_bias,b_bias,w);   //analytical double derivative
//    kineticEnergy = -0.5*m_system->getHamiltonian()->computeNumericalDoubleDerivative(particles);   //numerical double derivative

//compute the Potential, from the choices made by setting beta and omega_z at the beginning the potential is spherical or elliptical
   potentialEnergy=computePotentialEnergy();

if (interaction==true){
    interaction_potential=computeInteractionPotential();
}

    //cout<<"Kinetic energy = "<<kineticEnergy<<endl;
    return kineticEnergy + potentialEnergy+interaction_potential;
}

std::vector<double> HarmonicOscillator::omega() const
{
    return m_omega;
}

void HarmonicOscillator::setOmega(const std::vector<double> &omega)
{
    m_omega = omega;
}

double HarmonicOscillator::computeInteractionPotential(){
    double interaction_potential=0;
    for(int i=0; i<m_system->getNumberOfParticles(); i++){
        for(int j=0; j<i; j++){
            interaction_potential+=m_system->getDistanceMatrixij(i,j);
        }
    }
    return interaction_potential;
}


double HarmonicOscillator::computePotentialEnergy(){
    double potentialEnergy =0;
    for (int k = 0; k < m_system->getNumberOfVisibleNodes(); k++ ){
        potentialEnergy += m_omega[d]*m_omega[d]*X[d]*X[d];
    }
    return potentialEnergy * 0.5;
}
