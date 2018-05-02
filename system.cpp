#include "system.h"
#include <cassert>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <iostream>
#include <fstream>
#include <time.h>
#//include "conjugategradient.h"

using namespace std;

void System::metropolisStepBruteForce(bool interaction,vector<double> &X, vector<double> Hidden, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w) { //Brute Force Metropolis method
    int randparticle=Random::nextInt(getNumberOfParticles());

    vector <double> X_old=X;
    vector <double> X_new(getNumberOfVisibleNodes());



//choose a new move
    int init =randparticle*m_numberOfDimensions;
   for(int d = init; d < init+3; d++){
        X_new[d] = X_old[d] + m_stepLength*(Random::nextDouble()-0.5);
    }

//    m_particles.at(randparticle).setPosition(r_new); UPDATE FUNCTION DISTANCE MATRIX
//    updateDistanceMatrix(m_particles, randparticle);

    double psi_new = m_waveFunction->evaluate(X_new, Hidden, a_bias, b_bias, w);


    if (Random::nextDouble() <= psi_new * psi_new / (m_psiOld * m_psiOld)){ // Accept the new move
        m_psiOld = psi_new;

        X=X_new;
        getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(interaction,X_new,Hidden,a_bias,b_bias,w));

    }
    else{ // Don't accept accept the new move
        //updateDistanceMatrix(m_particles, randparticle); //update
    }
}


void System::metropolisStepImportance(bool interaction,vector<double> &X, vector<double> Hidden, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w) { //Importance Sampling method
    int randparticle=Random::nextInt(m_initialState->getNumberOfVisibleNodes());
    vector <double> X_old=X;
    vector <double> X_new(m_initialState->getNumberOfVisibleNodes());

    vector <vector<double>> QF_old(m_numberOfDimensions,vector<double>(m_numberOfParticles)); //UPDATE
    vector <vector<double>> QF_new(m_numberOfDimensions,vector<double>(m_numberOfParticles));//UPDATE
    QF_old=m_QuantumForce;

//Choose a new move
   for(int d = randparticle; d < m_numberOfDimensions; d++){
        X_new[d] = X_old[d] + 0.5*m_QuantumForce[d][randparticle]*m_timeStep + m_sqrtTimeStep*(Random::nextGaussian(0, 1));
//        cout << r_new[d] - r_old[d] << endl;
//        cout << "R_old = " << r_old[d] << endl;
//        cout << "R_new = " << r_new[d] << endl;
    }

   // m_particles.at(randparticle).setPosition(r_new);
    setQuantumForce(m_waveFunction->QuantumForce(m_particles)); //UPDATE
    QF_new = m_QuantumForce;
    updateDistanceMatrix(m_particles, randparticle);

// Compute Green function
    double GreensFunction =0.0;
    for(int d = 0; d < m_numberOfDimensions; d++){
        int j = randparticle;
        GreensFunction += 0.5*(QF_old[d][j ] + QF_new[d][j])
                *(0.5 * 0.5 * m_timeStep * (QF_old[d][j] - QF_new[d][j]) - r_new[d] + r_old[d]);
        }

    GreensFunction = exp(GreensFunction);
    double psi_new = m_waveFunction->evaluate(X,Hidden,a_bias,b_bias, w);

    // Accept  new move
    if (Random::nextDouble() <= GreensFunction*psi_new * psi_new / (m_psiOld * m_psiOld)){
        m_psiOld = psi_new;
        X=X_new;
        getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(interaction,X,Hidden,a_bias,b_bias, w));
    }

    // Don't accept new move
    else{
       // m_particles.at(randparticle).setPosition(r_old);
        X=X_old; //useless?
        setQuantumForce(QF_old);
        updateDistanceMatrix(m_particles, randparticle); //UPDATE
    }
}

void System::runMetropolisSteps(int numberOfMetropolisSteps,bool interaction, vector<double> X, vector<double> Hidden, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w) {
    //Principal function of the whole code. Here the Monte Carlo method is
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    //cout<<getSampler()->getEnergy()<<endl;
    getSampler()->setStepNumber(0);
    getSampler()->setAcceptedNumber(0);
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    // Initial values
  // do it again  setDistanceMatrix(computematrixdistance(m_particles));
    m_psiOld = m_waveFunction->evaluate(X,Hidden, a_bias, b_bias,w);
    getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(interaction, X,Hidden, a_bias, b_bias,w));
  // update  setQuantumForce(m_waveFunction->QuantumForce(m_particles));

    setHistogram();

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStepBruteForce(interaction, X,Hidden, a_bias, b_bias,w);   //run the system with Brute Force Metropolis
    //    bool acceptedStep = metropolisStepImportance(X,Hidden, a_bias, b_bias,w);     //run the system with Importance Sampling

        m_sampler->sample(acceptedStep,interaction,X,Hidden,a_bias,b_bias,w);                    //sample results and write them to file
// INSERT STOCHASTIC GRADIENT DESCENT HERE
        m_sampler->writeToFile();
    }
}


void System::StochasticGradientDescent(vector<double>& a_bias, vector<double>& b_bias, vector<std::vector<double>>& w){

    for(int i=0; i<m_numberOfVisibleNodes; i++){
       a_bias[i]-=m_learningRate*gradiente;
    }

    for(int i=0; i<m_numberOfHiddenNodes; i++){
       b_bias[i]-=m_learningRate*gradiente;
    }

    for(int i=0; i<m_numberOfVisibleNodes; i++){
        for(j=0;j<m_numberOfHiddenNodes; j++){
            w[i][j]-=m_learningRate*gradiente;
        }
    }
}

int System::getNumberOfParameters() const
{
    return m_numberOfParameters;
}

void System::setNumberOfParameters(int numberOfParameters)
{
    m_numberOfParameters = numberOfParameters;
}


vector<double> System::GradientParameters(){
    vector<double> Gradient(m_numberOfParameters);
    vector<double> First(m_numberOfVisibleNodes);


    for(int i=0; i<m_numberOfVisibleNodes; i++){
    First+=(X[i]-a_bias[i])/getSigma_squared();
    Gradient[i]=2*((m_sampler->getCumulativeEnergy()*First[i])/m_sampler->getAcceptedNumber()-m_sampler->getEnergy())
    }



}




void System::oneBodyDensity(){
    //Function to made the histrograms needed to compute the one body density
    vector<int> histogram(m_bins);
    double r2 = 0;
    for (int j = 0; j < getNumberOfParticles(); j++){
        r2 = 0;
        for (int d = 0; d < getNumberOfDimensions(); d++){
            r2 += m_particles.at(j).getPosition()[d]*m_particles.at(j).getPosition()[d];
        }
        r2 = sqrt(r2);
        int bucket = (int)floor(r2 / m_bucketSize);
        histogram[bucket] += 1;
    }

    // Update histogram
    for (int k = 0; k < m_bins; k++){
       m_histogram[k] += histogram[k];
    }
}

int System::getNumberOfVisibleNodes() const
{
    return m_numberOfVisibleNodes;
}

void System::setNumberOfVisibleNodes(int numberOfVisibleNodes)
{
    m_numberOfVisibleNodes = numberOfVisibleNodes;
}

int System::getNumberOfHiddenNodes() const
{
    return m_numberOfHiddenNodes;
}

void System::setNumberOfHiddenNodes(int numberOfHiddenNodes)
{
    m_numberOfHiddenNodes = numberOfHiddenNodes;
}

int System::getNumberOfVisibleNodes() const
{
    return m_numberOfVisibleNodes;
}

void System::setNumberOfVisibleNodes(int numberOfVisibleNodes)
{
    m_numberOfVisibleNodes = numberOfVisibleNodes;
}

int System::getNumberOfHiddenNodes() const
{
    return m_numberOfHiddenNodes;
}

void System::setNumberOfHiddenNodes(int numberOfHiddenNodes)
{
    m_numberOfHiddenNodes = numberOfHiddenNodes;
}

double System::getSigma() const
{
    return m_sigma;
}

void System::setSigma(double sigma)
{
    m_sigma = sigma;
}

double System::getSigma_squared() const
{
    return m_sigma_squared;
}

void System::setSigma_squared(double sigma_squared)
{
    m_sigma_squared = sigma_squared;
}

double System::getLearningRate() const
{
    return m_learningRate;
}

void System::setLearningRate(double learningRate)
{
    m_learningRate = learningRate;
}

void System::printOneBodyDensity(string filename){
    ofstream myFile;
    myFile.open(filename);
    for (int j = 0; j < m_bins; j++){
        myFile <<(double) m_histogram[j] /(getNumberOfParticles()*getNumberOfMetropolisSteps()*getEquilibrationFraction()*m_bucketSize) << "    " << j*m_bucketSize<< endl;
    }
    cout << "Printed ob density! " << endl;
}
void System::setHistogram()
{
    vector<int> histogram(getBins());
    m_histogram = histogram;
}

//double System::gradientDescent(double initialAlpha, string filename, int maxIterations){
//Gradient descent method to find the optimal variational parameter alpha given an initial parameter initialAlpha
//    int steepestDescentSteps = (int) 1e+5;
//    double alpha = initialAlpha;
//    double beta = getWaveFunction()->getParameters()[2] / getWaveFunction()->getParameters()[0];
//    double lambda = -0.001;
//    int iterations = 0;
//    double energyDerivative = 100;
//    double cumulativeAlpha = 0;
//    double tol = 1e-10;
//    double percentAlphasToSave = 0.3;
//    ofstream myFile;
//    myFile.open(filename);

//    while (iterations < maxIterations && fabs(energyDerivative) > tol){
//        vector<double> parameters(3);
//        parameters[0] = alpha;
//        parameters[1] = alpha;
//        parameters[2] = alpha*beta;
//        getWaveFunction()->setParameters(parameters);
//        runMetropolisSteps(steepestDescentSteps);
//        printOut();
//        energyDerivative = findEnergyDerivative();

//        // Make sure we accept enough moves (with interaction can get stuck)
//        if ((double)m_sampler->getAcceptedNumber() / steepestDescentSteps > 0.90){
//            alpha += lambda*energyDerivative;
//            iterations++;
//        }

//        cout << " New alpha = "  << alpha <<  endl;
//        cout << " Energy derivative = " << energyDerivative << endl;
//        cout << " Iterations = " << iterations << endl;

//        // Write alpha, mean local energy and st dev to file
//        myFile << alpha << "   "  << getSampler()->getEnergy() << "  " <<
//                  sqrt(getSampler()->getCumulativeEnergySquared() - getSampler()->getEnergy()*getSampler()->getEnergy())/getNumberOfMetropolisSteps() << endl;

//        if ((double) iterations / maxIterations > 1-percentAlphasToSave){
//            cumulativeAlpha += alpha;
//        }
//    }
//    myFile.close();

//    alpha = cumulativeAlpha / (maxIterations*percentAlphasToSave);
//    return alpha;
//}


int System::getBins() const
{
    return m_bins;
}

void System::setBins(int bins)
{
    m_bins = bins;
}

double System::getBucketSize() const
{
    return m_bucketSize;
}

void System::setBucketSize(double bucketSize)
{
    m_bucketSize = bucketSize;
}

void System::printOut()
{
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
}



double System::computedistance(int i){
    double temp=0;
    for(int j=0;j<m_numberOfDimensions;j++){
        temp+=m_particles.at(i).getPosition()[j]*m_particles.at(i).getPosition()[j];
    }
    return sqrt(temp);
}


bool System::updateDistanceMatrix( std::vector<class Particle> &particles, int randparticle){
    double temp = 0;
    for (int j = 0; j < randparticle; j++){ //need to understand what is randparticle
        temp = 0;
        for (int d = 0; d < m_numberOfDimensions; d++){
            temp += (particles.at(randparticle).getPosition()[d] - particles.at(j).getPosition()[d]) *
                    (particles.at(randparticle).getPosition()[d] - particles.at(j).getPosition()[d]);
        }
        m_distanceMatrix[randparticle][j] = sqrt(temp);
        if (m_distanceMatrix[randparticle][j] < getinteractionSize()){
            return true;
        }
        m_distanceMatrix[j][randparticle] = m_distanceMatrix[randparticle][j];
    }
//    for (int j = randparticle+1; j < m_numberOfParticles; j++){
//        temp = 0;
//        for (int d = 0; d < m_numberOfDimensions; d++){
//            temp += (particles.at(randparticle).getPosition()[d] - particles.at(j).getPosition()[d]) *
//                    (particles.at(randparticle).getPosition()[d] - particles.at(j).getPosition()[d]);
//        }
//        m_distanceMatrix[randparticle][j] = sqrt(temp);
//        if (m_distanceMatrix[randparticle][j] < getinteractionSize()){
//            return true;
//        }
//        m_distanceMatrix[j][randparticle] = m_distanceMatrix[randparticle][j];

//    }
    return false;
}

std::vector<vector<double>> System::computematrixdistance(vector<double>& m_X){

    vector<vector<double>> distancematrix(m_numberOfParticles, vector<double>(m_numberOfParticles));
    double temp=0;
    int j=0;
    int z;
    int k=0;
    while(j < m_numberOfVisibleNodes){
        temp = 0;
        z=0;
            for(int i=0;i<j-2;i+=3){
                temp=(m_X[i]-m_X[j])*(m_X[i]-m_X[j])+(m_X[i+1]-m_X[j+1])*(m_X[i+1]-m_X[j+1])+(m_X[i+2]-m_X[j+2])*(m_X[i+2]-m_X[j+2]);
                distancematrix[z][k]=sqrt(temp);
                distancematrix[k][z]=distancematrix[z][k];
                z++;
            }
        j+=3;
        k++;
    }
    return distancematrix;
}

void System::updateQuantumForce(std::vector<std::vector<double> > deltaQuantumForce, bool subtract){
    if (subtract == false){
        for (int d = 0; d < m_numberOfDimensions; d++){
        for(int i=0; i<=m_numberOfParticles; i++){
            m_QuantumForce[d][i] += deltaQuantumForce[d][i];
        }
        }
    }
    else {
        for (int d = 0; d < m_numberOfDimensions; d++){
        for(int i=0; i<=m_numberOfParticles; i++){
            m_QuantumForce[d][i] += deltaQuantumForce[d][i];
        }
        }
    }
}



double System::computedistanceABS(int i, int j){
    double temp=0;
    for(int k=0;k<m_numberOfDimensions;k++){
        temp+=(m_particles.at(i).getPosition()[k] - m_particles.at(j).getPosition()[k]) *
                (m_particles.at(i).getPosition()[k] - m_particles.at(j).getPosition()[k]);
    }
    return sqrt(temp);{
    }
}

void System::openDataFile(string filename){
    m_sampler->openDataFile(filename);
}

double System::getinteractionSize() const
{
    return m_interactionSize;
}

void System::setinteractionSize(double interactionSize)
{
    m_interactionSize = interactionSize;
}

double System::getTimeStep() const
{
    return m_timeStep;
}

void System::setTimeStep(double timeStep)
{
    m_timeStep = timeStep;
}

double System::getSqrtTimeStep() const
{
    return m_sqrtTimeStep;
}

void System::setSqrtTimeStep(double sqrtTimeStep)
{
    m_sqrtTimeStep = sqrtTimeStep;
}

std::vector<vector<double> > System::getDistanceMatrix() const
{
    return m_distanceMatrix;
}

double System::getDistanceMatrixij(int i, int j) const
{
    return m_distanceMatrix[i][j];
}

double System::getPsiOld() const
{
    return m_psiOld;
}

void System::setPsiOld(double psiOld)
{
    m_psiOld = psiOld;
}

std::vector<vector<double>> System::getQuantumForce() const
{
    return m_QuantumForce;
}

void System::setQuantumForce(const std::vector<vector<double>> &QuantumForce)
{
    m_QuantumForce = QuantumForce;
}

void System::setDistanceMatrix(const std::vector<vector<double> > &distanceMatrix)
{
    m_distanceMatrix = distanceMatrix;
}


void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}


double System::findEnergyDerivative()
{

    double meanEnergy      = getSampler()->getCumulativeEnergy() / (m_numberOfMetropolisSteps*getEquilibrationFraction());
    double meanWFderiv     = getSampler()->getCumulativeWFderiv() / (m_numberOfMetropolisSteps*getEquilibrationFraction());
    double meanWFderivEloc =  getSampler()->getCumulativeWFderivMultEloc() / (m_numberOfMetropolisSteps*getEquilibrationFraction());




    return 2 * (meanWFderivEloc - meanEnergy*meanWFderiv);
}

