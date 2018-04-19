#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>

using namespace std;
using std::vector;

SimpleGaussian::SimpleGaussian(System* system, double numberOfHiddenNodes, double numberOfVisibleNodes, vector<double> X, vector<double> Hidden, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w) :
    WaveFunction(system) {
    assert(alpha >= 0);
    assert(beta >= 0);
    m_numberOfParameters = 3;
    m_parameters.reserve(3);
    m_parameters.push_back(alpha);
    m_parameters.push_back(alpha);
    m_parameters.push_back(alpha*beta);
    //m_parameters.push_back(beta);
}

double SimpleGaussian::evaluate(vector<double> X, vector<double> Hidden, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w) {

    // WF is just the product of the individual wavefunctions
    // Divide by the wf of the old particle and multiply in the wf of the new one

















    double r_squared = 0;
    double f=1;

    for(int j=0; j<m_system->getNumberOfParticles();j++){

        if(m_system->getNumberOfParticles()==1) break;
        for(int i=0; i<j; i++){
            if(m_system->getDistanceMatrixij(i,j)<m_system->getinteractionSize()){f=0.0;cout<<"Hey"<<endl; break;}
            f *= 1-m_system->getinteractionSize()/(m_system->getDistanceMatrixij(i,j));
    }
}

for(int i=0;i < m_system->getNumberOfParticles();i++){

        for(int d=0; d < m_system->getNumberOfDimensions();d++){
            r_squared += particles.at(i).getPosition()[d]*particles.at(i).getPosition()[d]*m_parameters[d];
                  }
}

return exp(-r_squared)*f;

}


double SimpleGaussian::computeDoubleDerivative(vector<double> X, vector<double> Hidden, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w) {
//function to compute tha analytical double drivative
    double one=0;
    double interaction=0;

    double a=m_system->getinteractionSize();

    //Interaction terms
    for(int i=0; i<m_system->getNumberOfParticles(); i++){
        double r_i_square=0;
        for(int d = 0; d < m_system->getNumberOfDimensions() - 1; d++){
         r_i_square += particles.at(i).getPosition()[d]*
                       particles.at(i).getPosition()[d];
         }
         int d = m_system->getNumberOfDimensions()-1;
         r_i_square += particles.at(i).getPosition()[d]*
                       particles.at(i).getPosition()[d]*m_parameters[2]/(m_parameters[0]);

       double second=0;
       double third=0;
       double fourth=0;
       double fifth=0;
       double temp;

       for(int j=0; j < i; j++) {

           double r_ij = m_system->getDistanceMatrixij(i,j);

           temp= a / ( (r_ij-a) * r_ij );

           second += temp;

           double r_ir_j = 0;
           for(int d = 0; d < m_system->getNumberOfDimensions() - 1; d++){
               r_ir_j += particles.at(i).getPosition()[d]*
                         particles.at(j).getPosition()[d];
           }
           int d = m_system->getNumberOfDimensions() - 1;
               r_ir_j += particles.at(i).getPosition()[d]*
                         particles.at(j).getPosition()[d]*
                         m_parameters[2]/(m_parameters[0]);

           fourth-= temp * temp;

           fifth -= 4 * m_parameters[0]  * (r_i_square - r_ir_j) * temp/
                   ( r_ij );

       }
       for(int j = j+1; j < m_system->getNumberOfParticles(); j++){

           double r_ij = m_system->getDistanceMatrixij(i,j);

           temp = a / ( (r_ij-a) * r_ij );

           second += temp;

           double r_ir_j = 0;
           for(int d = 0; d < m_system->getNumberOfDimensions() - 1; d++){
               r_ir_j += particles.at(i).getPosition()[d]*
                         particles.at(j).getPosition()[d];
           }
           int d = m_system->getNumberOfDimensions() - 1;
               r_ir_j += particles.at(i).getPosition()[d]*
                         particles.at(j).getPosition()[d]*
                         m_parameters[2]/(m_parameters[0]);

           fourth-= temp * temp;

           fifth -= 4 * m_parameters[0]  * (r_i_square - r_ir_j) * temp/
                   ( r_ij );
       }

       third=second*second;

       interaction+=second+third+fourth+fifth;

    }


//One body term
    for(int i = 0; i < m_system->getNumberOfParticles(); i++){
        for(int d = 0; d < m_system->getNumberOfDimensions(); d++){
            one += m_parameters[d]*m_parameters[d]*
                   particles.at(i).getPosition()[d]*
                   particles.at(i).getPosition()[d];
        }
    }

    one*=4.0;
    one-= 2 * ( (m_system->getNumberOfDimensions() - 1) * m_parameters[0] + m_parameters[2])
            * m_system->getNumberOfParticles(); //constant term

    return one+interaction;
}

std::vector<vector<double>> SimpleGaussian::QuantumForce(std::vector<class Particle>& particles) {
//Function to comput the Quantum Force for the Importance Sampling method

    double a = m_system->getinteractionSize() ;
    double constant;
    double R_kj;
    double dimension=m_system->getNumberOfDimensions();
    double number =m_system->getNumberOfParticles();
    std::vector<std::vector<double>> QuantumForce(dimension,vector<double>(number));
    for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
        for (int k = 0; k < m_system->getNumberOfParticles(); k++){
    QuantumForce[d][k] = -2 * (m_parameters[d]*particles.at(k).getPosition()[d]);
        for (int j = 0; j < k; j++){
                R_kj = m_system->getDistanceMatrixij(k,j);
                constant = 2*a / (R_kj*R_kj*(R_kj-a));

                    QuantumForce[d][k] += (particles.at(k).getPosition()[d] - particles.at(j).getPosition()[d]) * constant;

            }
        }
    }
    return QuantumForce;
}


