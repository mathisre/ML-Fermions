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

std::vector<double> explist(double first, double last, double size)
{
    //if(first>last) std::swap(first,last);
    double expfirst = exp(first);
    double explast = exp(last);
    double step = -(explast-expfirst)/(size-1);
    std::vector<double> out;
    for(double x=explast; x<=expfirst; x+=step)
    {
        double a = log(x);
        out.push_back(a);
    }
    return out;
}
int main(){

    int numberOfParticles   = 2;        //Number of particles of the system considered
    int numberOfDimensions  = 2;         // NUmber of dimensions
    int numberOfHiddenNodes = 2;
    int numberOfVisibleNodes =numberOfDimensions*numberOfParticles;

    double sigma =1.0;

    int numberOfParameters = numberOfVisibleNodes + numberOfHiddenNodes
            + numberOfVisibleNodes * numberOfHiddenNodes;


    bool interaction = false;

    std::vector<double> X(numberOfVisibleNodes);

    std::vector<double> Hidden=std::vector<double>(numberOfHiddenNodes);

    std::vector<double> a_bias=std::vector<double>(numberOfVisibleNodes);

    std::vector<double> b_bias=std::vector<double>(numberOfHiddenNodes);

    std::vector<std::vector<double>> w(numberOfVisibleNodes, vector<double>(numberOfHiddenNodes));//=std::vector<std::vector<double>>();


    double alpha            = 0.50;      // Variational parameter.
    double beta             = 1.0;            // for interacting case: beta=2.82843
    int numberOfSteps       = (int) 10000;   // NUmber of Monte Carlo steps
    double interactionSize  = 0.0; // for interacting case: interactionSize=0.0043;

    double timeStep         = 0.3;        // Importance sampling time step
    double stepLength       = 0.45;          // Metropolis step length.
    double omega            = 1.0;          // Oscillator frequency.
    double omega_z          = beta;         // Oscillator frequency in z-direction
    double equilibration    = 0.2;          // Amount of the total steps used for equilibration.

    // Parameters for onebody density histogram
    double bucketSize = 0.01;
    int bins = ceil(4 / bucketSize);

    //    string filename = "0";
    string filename         = "0" + to_string(alpha) + "_dim_" + to_string(numberOfDimensions) + "_particles_" + to_string(numberOfParticles)  + ".dat";
    // Set filename to "0" to stop from writing to file


    double learning_rate = 0.02;//list[i];

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z));
    system->setInitialState             (new RandomUniform(system, numberOfParticles, numberOfDimensions, numberOfHiddenNodes, numberOfVisibleNodes, sigma, X, Hidden, a_bias, b_bias, w, interactionSize, timeStep, bins, bucketSize,numberOfParameters));
    system->setWaveFunction             (new SimpleGaussian(system));//, a_bias, b_bias, w));
    system->openDataFile                (filename);
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);


    int TotalNumberOfCycles = 1000;
  // string method = "MetropolisBruteForce";
    string method = "MetropolisImportance";
   // string method = "Gibbs";

    vector <double> Gradient(numberOfParameters);
    //vector<double> list(TotalNumberOfCycles);
    //list=explist(0.001,0.0005,1000);

    //for(int j=0; j<list.size(); j++){
    //    cout<<list[j]<<endl;
    //}
    auto start = std::chrono::system_clock::now();
    //int i=TotalNumberOfCycles;

    for(int cycles = 0; cycles < TotalNumberOfCycles; cycles ++){

        system->setLearningRate           (learning_rate);

        system->setNumberOfParameters     (numberOfParameters);

        system->runMetropolisSteps        (method, Gradient,numberOfSteps,interaction, X, Hidden, a_bias, b_bias, w);

        system->StochasticGradientDescent (Gradient,X,a_bias,b_bias,w);

        system->printOut                  (cycles);



        string densityFilename = "density_alpha_" + to_string(alpha) + "_beta_" + to_string(beta) + "_inter.dat";
        //    system->printOneBodyDensity(densityFilename);

        //learning_rate-=1.995*0.0001;
        //i--;
    }
    //time

    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> diff = end - start;

    std::cout << " Computation time = " << diff.count() / 60.0 << " min\n" << endl; //display run time

    return 0;
}


