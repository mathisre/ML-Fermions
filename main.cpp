#include <iostream>
#include <random>
#include <cmath>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include <chrono>
#include <string>
#include <vector>

using namespace std;

int main(){

    int numberOfParticles   = 10;        //Number of particles of the system considered
    int numberOfDimensions  = 3;         // NUmber of dimensions

    int numberOfHiddenNodes = 2;
    int numberOfVisibleNodes =numberOfDimensions*numberOfParticles;

   double sigma =0.1;

int numberOfParameters=numberOfVisibleNodes+numberOfHiddenNodes+numberOfVisibleNodes*numberOfHiddenNodes;
    double learning_rate = 0.2;

        std::vector<double> X=std::vector<double>();
        std::vector<double> Hidden=std::vector<double>();
        std::vector<double> a_bias=std::vector<double>();
        std::vector<double> b_bias=std::vector<double>();
        std::vector<std::vector<double>> w=std::vector<std::vector<double>>;


    X.resize(numberOfVisibleNodes);
    a_bias.resize(numberOfVisibleNodes);
    Hidden.resize(numberOfHiddenNodes);
    b_bias.resize(numberOfHiddenNodes);
    for(int i=0; i<numberOfVisibleNodes; i++){
            w[i].resize(numberOfHiddenNodes);
    }
    double alpha            = 0.50;      // Variational parameter.
    double beta             = 1.0;            // for interacting case: beta=2.82843
    int numberOfSteps       = (int) 1e+6;   // NUmber of Monte Carlo steps
    double interactionSize  = 0.0; // for interacting case: interactionSize=0.0043;

    double timeStep         = 0.01;        // Importance sampling time step
    double stepLength       = 1.0;          // Metropolis step length.
    double omega            = 1.0;          // Oscillator frequency.
    double omega_z = beta;                  // Oscillator frequency in z-direction
    double equilibration    = 0.2;          // Amount of the total steps used for equilibration.

    // Parameters for onebody density histogram
    double bucketSize = 0.01;
    int bins = ceil(4 / bucketSize);

//    string filename = "0";
    string filename         = "Alpha_" + to_string(alpha) + "_dim_" + to_string(numberOfDimensions) + "_particles_" + to_string(numberOfParticles)  + ".dat";
    // Set filename to "0" to stop from writing to file

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z));
    system->setInitialState             (new RandomUniform(system, numberOfParticles, numberOfDimensions, sigma, X, Hidden, a_bias, b_bias, w, interactionSize, timeStep, bins, bucketSize,learning_rate,numberOfParameters));
    system->setWaveFunction             (new SimpleGaussian(system, numberOfHiddenNodes, numberOfVisibleNodes, X, Hidden, a_bias, b_bias, w,));
    system->openDataFile                (filename);
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);



    // Gradient descent method to find energy minimum

//    int maxIterations = 200;
//    double initialAlpha = 0.60;
//    string minFilename = "find_minimum_" + to_string(initialAlpha) +"_N_" + to_string(numberOfParticles) + "_doublecheck.dat";
//    alpha = system->gradientDescent(initialAlpha, minFilename, maxIterations);
//    cout << "Optimal alpha found by steepest descent: " << alpha << endl;
//    vector<double> parameters(3);
//    parameters[0] = alpha;
//    parameters[1] = alpha;
//    parameters[2] = alpha*beta;
//    system->getWaveFunction()->setParameters(parameters);


    auto start = std::chrono::system_clock::now();
    system->runMetropolisSteps          (numberOfSteps, X, Hidden, a_bias, b_bias, w,);
    system->printOut();


    string densityFilename = "density_alpha_" + to_string(alpha) + "_beta_" + to_string(beta) + "_inter.dat";
//    system->printOneBodyDensity(densityFilename);

    //time
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end-start;
    std::cout << " Computation time = " << diff.count() << " s\n" << endl; //display run time
    return 0;
}
