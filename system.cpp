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

bool System::metropolisStepBruteForce(double GibbsValue, bool interaction,vector<double> &X, vector<double> Hidden, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w) { //Brute Force Metropolis method
    //    int randparticle=Random::nextInt(getNumberOfParticles());
    int randparticle=Random::nextInt(getNumberOfVisibleNodes());

    vector <double> X_old (getNumberOfVisibleNodes());
    vector <double> X_new (getNumberOfVisibleNodes());

    X_old = X;
    X_new = X;

    //    for(int i=0; i<getNumberOfVisibleNodes(); i++){
    //        cout<<"ora"<<X[i]<<endl;
    //    }
    vector <vector<double>> OldDistanceMatrix;
    OldDistanceMatrix=getDistanceMatrix();
    //choose a new move
    int init =randparticle*m_numberOfDimensions;

    //  cout<<"random index: "<<randparticle<<endl;
    // for(int d = init; d < init+m_numberOfDimensions; d++){

    X_new[randparticle] = X_old[randparticle] + m_stepLength *( Random::nextDouble() - 0.5 );

    //   cout<<"x_new: ["<<randparticle<<"]"<<X_new[randparticle]<<endl;
    //  cout<<"x_old: ["<<randparticle<<"]"<<X_old[randparticle]<<endl;
    //   }

    //    for(int i=0; i<getNumberOfVisibleNodes(); i++){
    //        cout<<"after"<<X_new[i]<<endl;
    //    }
    //    m_particles.at(randparticle).setPosition(r_new); UPDATE FUNCTION DISTANCE MATRIX
    //    updateDistanceMatrix(m_particles, randparticle);
    setDistanceMatrix (computematrixdistance(X_new));

    double psi_new = m_waveFunction->evaluate(GibbsValue, X_new, Hidden, a_bias, b_bias, w);

    //   cout<<"new"<<psi_new<<endl;

    double prob = psi_new * psi_new / ( m_psiOld * m_psiOld );

    if (Random::nextDouble() < prob||1.0<prob){ // Accept the new move

        m_psiOld = psi_new;
        X        = X_new;

        //  getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(interaction,X_new,Hidden,a_bias,b_bias,w));
        // cout<<"ac"<<getSampler()->getEnergy()<<endl;

        return true;

    } else { // Don't accept accept the new move

        X = X_old;
        //  cout<<m_psiOld<<endl;
        // cout<<"ehi"<<endl;
        //updateDistanceMatrix(m_particles, randparticle); //update
        //    getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(interaction,X_old,Hidden,a_bias,b_bias,w));
        //  cout<<"notac"<<getSampler()->getEnergy()<<endl;
setDistanceMatrix(OldDistanceMatrix);
        return false;

    }
}


bool System::metropolisStepImportance(double GibbsValue, bool interaction,vector<double> &X, vector<double> Hidden, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w) { //Importance Sampling method

    int randparticle = Random::nextInt (getNumberOfVisibleNodes());

    vector <double> X_old (getNumberOfVisibleNodes());
    vector <double> X_new (getNumberOfVisibleNodes());
    X_old = X;
    X_new = X;

    vector <double> QF_old (m_numberOfVisibleNodes);
    vector <double> QF_new (m_numberOfVisibleNodes);
    QF_old = getQuantumForce();

    vector <vector<double>> OldDistanceMatrix;
    OldDistanceMatrix=getDistanceMatrix();
//    for(int i=0; i<m_numberOfParticles; i++){
//        for(int j=0; j<m_numberOfParticles; j++){
//            cout<<m_distanceMatrix[i][j]<<endl;
//        }
//    }
    //int init = randparticle * m_numberOfDimensions;


    //randparticle=2;

    //Choose a new move
    //   for(int d = init; d < init+ m_numberOfDimensions; d++){

    X_new[randparticle] = X_old[randparticle] + 0.5 * m_QuantumForce[randparticle] * m_timeStep + m_sqrtTimeStep * ( Random::nextGaussian(0.0, 1.0) );


    // cout << X_new[randparticle]-X_old[randparticle]<< endl;
    //        cout << "R_old = " << r_old[d] << endl;
    //        cout << "R_new = " << r_new[d] << endl;
    //   }
    //X_new[randparticle]=1.0;
    //for(int i=0; i<X_new.size(); i++){
    //cout<<X_new[i]<<endl;
    //}
    //  cout<<"old"<<QF_old[randparticle]<<endl;
    // m_particles.at(randparticle).setPosition(r_new);

    setQuantumForce (m_waveFunction->QuantumForce(GibbsValue,X_new,a_bias,b_bias,w)); //UPDATE
    QF_new = m_QuantumForce;

    //   cout<<X_new[randparticle]<<endl;
    updateDistanceMatrix(X_new, randparticle);

//    for(int i=0; i<m_numberOfParticles; i++){
//        for(int j=0; j<m_numberOfParticles; j++){
//            cout<<m_distanceMatrix[i][j]<<endl;
//        }
//    }

  //  setDistanceMatrix (computematrixdistance(X_new));

//    for(int i=0; i<m_numberOfParticles; i++){
//        for(int j=0; j<m_numberOfParticles; j++){
//            cout<<m_distanceMatrix[i][j]<<endl;
//        }
//    }
    // Compute Green function
    double GreensFunction = 0.0;
    //  for(int d = 0; d < m_numberOfDimensions; d++){
    //        GreensFunction += 0.5*(QF_old[init+d] + QF_new[init+d])
    //                *(0.5 * 0.5 * m_timeStep * (QF_old[init+d] - QF_new[init+d]) - X_new[init+d] + X_old[init+d]);
    //    }
    //   cout<<"QF_old: "<<QF_old[randparticle]<<endl;
    //  cout<<"QF_new: "<<QF_new[randparticle]<<endl;
    GreensFunction = 0.5* (QF_old[randparticle] + QF_new[randparticle])
                   * (0.5 * 0.5 * m_timeStep * (+QF_old[randparticle] - QF_new[randparticle])
                       - X_new[randparticle] + X_old[randparticle] );

    GreensFunction = exp( GreensFunction );

    double psi_new = m_waveFunction->evaluate(GibbsValue,X_new,Hidden,a_bias,b_bias, w);
    //    cout<<"step: "<<m_timeStep<<endl
    ;
    //    cout<<"old: "<<m_psiOld<<endl;
    //    cout<<"new: "<<psi_new<<endl;
    //    cout<<"green: "<<GreensFunction<<endl;

    double prob = GreensFunction * psi_new * psi_new / (m_psiOld * m_psiOld);
    //cout<<"prob: "<<prob<<endl;
    // Accept  new move

    if ( (Random::nextDouble() < prob) || ( 1.0<prob )){

        m_psiOld = psi_new;
        X        = X_new;

        //cout<<"ehi"<<endl;
        // getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(interaction,X_new,Hidden,a_bias,b_bias, w)); //?? need?
        return true;
    }

    // Don't accept new move
    else{
        // cout<<"ehi"<<endl;
        // m_particles.at(randparticle).setPosition(r_old);
        X = X_old;

        setQuantumForce   (QF_old);
  //setDistanceMatrix(OldDistanceMatrix);
updateDistanceMatrix(X_old,randparticle);
        //setDistanceMatrix (OldDistanceMatrix);

        //     getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(interaction,X_old,Hidden,a_bias,b_bias, w));
        return false;
    }
}

bool System::GibbsSampling(vector<double> &X, vector<double> Hidden, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w){

    int N = getNumberOfHiddenNodes();
    int M = getNumberOfVisibleNodes();

    //vector <double> argument(N);

    double sum;
    double sum2;
    double mean;
    double argument;

    for(int j = 0; j < N; j++){
        sum = 0;

        for(int i = 0; i < M; i++){
            sum += X[i] * w[i][j] / getSigma_squared();
        // cout<<"old X["<<i<<"]"<<X[i]<<endl;
        }

        argument = exp( - b_bias[j] - sum);
        Hidden[j]=1/(1+argument);
    }

    for( int i = 0; i < M; i++){

        sum2 = 0;

            for( int j = 0; j < N; j++){

            sum2 += Hidden[j] * w[i][j];

        }

        mean = a_bias[i] + sum2;
        X[i] = Random::nextGaussian(mean,getSigma());
       // cout<<"new X["<<i<<"]"<<X[i]<<endl;

    }

    return true;
}

void System::runMetropolisSteps(string method, vector<double>&Gradient,int numberOfMetropolisSteps,bool interaction, vector<double> X, vector<double> Hidden, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w) {
    //Principal function of the whole code. Here the Monte Carlo method is
    // m_particles                 = m_initialState->getParticles();
    double GibbsValue = 1.0;
    if(method=="Gibbs") GibbsValue=0.5;
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    //cout<<getSampler()->getEnergy()<<endl;
    getSampler()->setStepNumber(0);
    getSampler()->setAcceptedNumber(0);
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    m_sampler->setDimensionOfGradient(m_numberOfParameters);
    // Initial values
    setDistanceMatrix(computematrixdistance(X));
//    for(int i=0; i<m_numberOfParticles; i++){
//        for(int j=0; j<m_numberOfParticles; j++){
//            cout<<m_distanceMatrix[i][j]<<endl;
//        }
//    }
    m_psiOld = m_waveFunction->evaluate(GibbsValue, X,Hidden, a_bias, b_bias,w);
    //  getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(interaction, X,Hidden, a_bias, b_bias,w));
    //cout<<"ehi"<<endl;
    setQuantumForce(m_waveFunction->QuantumForce(GibbsValue, X,a_bias,b_bias,w));

    setHistogram();
    vector<double> temp (m_numberOfParameters);
    setCumulativeEnGradient(temp);
    setCumulativeGradient(temp);
    setGradient(temp);
    setEnGradient_average(temp);
bool acceptedStep;
    for (int i=0; i < numberOfMetropolisSteps; i++) {
        if( method == "MetropolisBruteForce" ) acceptedStep =  metropolisStepBruteForce(GibbsValue,interaction, X,Hidden, a_bias, b_bias,w);   //run the system with Brute Force Metropolis
        if( method == "MetropolisImportance" ) acceptedStep =  metropolisStepImportance(GibbsValue, interaction,X,Hidden, a_bias, b_bias,w);    //run the system with Importance Sampling
        if( method == "Gibbs"                ) acceptedStep =  GibbsSampling(X,Hidden,a_bias,b_bias,w);

        m_sampler->sample(GibbsValue, acceptedStep,interaction,X,Hidden,a_bias,b_bias,w);                    //sample results and write them to file
        m_sampler->writeToFile();
    }
    //StochasticGradientDescent(Gradient,X,a_bias,b_bias,w);
    m_sampler->computeAverages(Gradient);
}


void System::StochasticGradientDescent(vector<double> Gradient, vector<double> X, vector<double>& a_bias, vector<double>& b_bias, vector<std::vector<double>>& w){

    for( int i = 0; i < m_numberOfVisibleNodes; i++){

        a_bias[i] -= m_learningRate * Gradient[i];
        // cout<<"a:["<<i<<"]"<<a_bias[i]<<endl;
    }

    for(int i = 0; i < m_numberOfHiddenNodes; i++){

        b_bias[i] -= m_learningRate * Gradient[i + m_numberOfVisibleNodes];
        //cout<<"b:["<<i<<"]"<<b_bias[i]<<endl;
    }

    int z = m_numberOfVisibleNodes + m_numberOfHiddenNodes;
    for(int i = 0; i < m_numberOfVisibleNodes; i++){

        for(int j = 0; j < m_numberOfHiddenNodes; j++){

            w[i][j] -= m_learningRate * Gradient[z];
            z++;
            //      cout<<"w:["<<i<<"]["<<j<<"]"<< w[i][j]<<endl;
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


vector<double> System::GradientParameters(double GibbsValue, vector<double> X, vector<double>& a_bias, vector<double>& b_bias, vector<std::vector<double>>& w){
    vector<double> GradientPsi (m_numberOfParameters);
    vector<double> argument    (m_numberOfHiddenNodes);

    double sum;

    for(int j = 0; j < m_numberOfHiddenNodes; j++){

        sum=0;

        for(int i = 0; i < m_numberOfVisibleNodes; i++){

            sum += X[i] * w[i][j] / getSigma_squared();

        }

        argument[j] = b_bias[j] + sum;
    }

    for(int i = 0; i < m_numberOfVisibleNodes; i++){
        GradientPsi[i] = ( X[i] - a_bias[i] ) * GibbsValue / getSigma_squared();
    }


    for(int i = m_numberOfVisibleNodes; i < m_numberOfHiddenNodes + m_numberOfVisibleNodes; i++){
        GradientPsi[i] = GibbsValue / ( 1 + exp ( - argument[i - m_numberOfVisibleNodes] ) );
    }

    int z = m_numberOfVisibleNodes + m_numberOfHiddenNodes;
    for(int i = 0; i < m_numberOfVisibleNodes; i++){

        for(int j=0; j<m_numberOfHiddenNodes; j++){

            GradientPsi[z] = GibbsValue * X[i] / ( getSigma_squared() * ( 1 + exp ( - argument[j] ) ) );
            z++;

        }

    }

    return GradientPsi;
}

std::vector<double> System::getCumulativeGradient() const
{
    return m_cumulativeGradient;
}

void System::setCumulativeGradient(const std::vector<double> &cumulativeGradient)
{
    m_cumulativeGradient = cumulativeGradient;
}

std::vector<double> System::getCumulativeEnGradient() const
{
    return m_cumulativeEnGradient;
}

void System::setCumulativeEnGradient(const std::vector<double> &cumulativeEnGradient)
{
    m_cumulativeEnGradient = cumulativeEnGradient;
}

std::vector<double> System::getGradient() const
{
    return m_Gradient;
}

void System::setGradient(const std::vector<double> &Gradient)
{
    m_Gradient = Gradient;
}

std::vector<double> System::getEnGradient_average() const
{
    return m_EnGradient_average;
}

void System::setEnGradient_average(const std::vector<double> &EnGradient_average)
{
    m_EnGradient_average = EnGradient_average;
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

void System::printOut(int cycle)
{
    // m_sampler->computeAverages();

    m_sampler->printOutputToTerminal(cycle);
}



double System::computedistance(int i){
    double temp=0;
    for(int j=0;j<m_numberOfDimensions;j++){
        temp+=m_particles.at(i).getPosition()[j]*m_particles.at(i).getPosition()[j];
    }
    return sqrt(temp);
}

int System::computeIndex(int index){
    int init=index;
   //cout<<"begin"<<index<<endl;
    if((index%m_numberOfDimensions)==0) {
   //cout<<"end"<<init<<endl;
        return index;
    }

    while(index>m_numberOfDimensions){
    index -=m_numberOfDimensions;
    }
    init -=index;
  // cout<<"end"<<init<<endl;
   return init;
}

void System::updateDistanceMatrix(std::vector<double> m_X, int randparticle){
    double temp = 0;
  double  init=computeIndex(randparticle);
//cout<<"ehi"<<endl;
    double part=init/m_numberOfDimensions;

            //int init=randparticle*m_numberOfDimensions;
    for (int j = 0; j < init; j+=m_numberOfDimensions){
        temp = 0;
        for (int d = 0; d < m_numberOfDimensions; d++){
            temp += (m_X[j+d] - m_X[init+d]) *
                    (m_X[j+d] - m_X[init+d]);
           // cout<<temp<<endl;
        }
        m_distanceMatrix[part][j/m_numberOfDimensions] = sqrt(temp);
       // cout<<part<<endl;
        //cout<<j/m_numberOfDimensions<<endl;
//        if (m_distanceMatrix[randparticle][j] < getinteractionSize()){
//            return true;
//        }
        m_distanceMatrix[j/m_numberOfDimensions][part] = m_distanceMatrix[part][j/m_numberOfDimensions];
    }
    for (int j = init+m_numberOfDimensions; j < m_numberOfVisibleNodes; j+=m_numberOfDimensions){
        temp = 0;
     //   cout<<"ehi"<<endl;
        for (int d = 0; d < m_numberOfDimensions; d++){
            temp += (m_X[j+d] - m_X[init+d]) *
                    (m_X[j+d] - m_X[init+d]);
        }
        m_distanceMatrix[part][j/m_numberOfDimensions] = sqrt(temp);
//        if (m_distanceMatrix[randparticle][j] < getinteractionSize()){
//            return true;
//        }
        m_distanceMatrix[j/m_numberOfDimensions][part] = m_distanceMatrix[part][j/m_numberOfDimensions];

    }

}

std::vector<vector<double>> System::computematrixdistance(vector<double>& m_X){

    vector<vector<double>> distancematrix(m_numberOfParticles, vector<double>(m_numberOfParticles));
    double temp=0;
    int j=0;
    int z;
    int k=0;

    while(j < m_numberOfVisibleNodes){

        temp = 0;
        z    = 0;

        for(int i = 0;i < j; i += m_numberOfDimensions){

            for(int q = 0; q < m_numberOfDimensions; q++){

                temp += ( m_X[i + q] - m_X[j + q] ) * ( m_X[i + q] - m_X[j + q] );

            }
            //temp=(m_X[i]-m_X[j])*(m_X[i]-m_X[j])+(m_X[i+1]-m_X[j+1])*(m_X[i+1]-m_X[j+1])+(m_X[i+2]-m_X[j+2])*(m_X[i+2]-m_X[j+2]);

            distancematrix[z][k] = sqrt(temp);

            distancematrix[k][z] = distancematrix[z][k];

            z++;
        }

        j += m_numberOfDimensions;
        k++;
    }
    return distancematrix;
}

void System::updateQuantumForce(std::vector<std::vector<double> > deltaQuantumForce, bool subtract){
    //    if (subtract == false){
    //        for (int d = 0; d < m_numberOfDimensions; d++){
    //        for(int i=0; i<=m_numberOfParticles; i++){
    //            m_QuantumForce[d][i] += deltaQuantumForce[d][i];
    //        }
    //        }
    //    }
    //    else {
    //        for (int d = 0; d < m_numberOfDimensions; d++){
    //        for(int i=0; i<=m_numberOfParticles; i++){
    //            m_QuantumForce[d][i] += deltaQuantumForce[d][i];
    //        }
    //        }
    //    }
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

std::vector<double> System::getQuantumForce() const
{
    return m_QuantumForce;
}

void System::setQuantumForce(const std::vector<double> &QuantumForce)
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

