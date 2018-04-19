#pragma once
#include <string>
class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();

    void openDataFile(std::string filename);
    void writeToFile();

    void computeAverages();
    double getEnergy()          { return m_energy; }
    int getStepNumber() const;

    int getNumberOfMetropolisSteps() const;

    double getWFderivMultELoc() const;
    void setWFderivMultELoc(double WFderivMultELoc);

    double getCumulativeWF() const;
    void setCumulativeWF(double cumulativeWF);

    double getWFderiv() const;
    void setWFderiv(double WFderiv);

    double getCumulativeWFderiv() const;
    void setCumulativeWFderiv(double cumulativeWFderiv);

    double getCumulativeWFderivMultEloc() const;
    void setCumulativeWFderivMultEloc(double cumulativeWFderivMultEloc);

    int getAcceptedNumber() const;
    void setAcceptedNumber(int acceptedNumber);

    void setStepNumber(int stepNumber);

    void setEnergy(double energy);

    void updateEnergy(double dE);
    double getCumulativeEnergy() const;

    double getCumulativeEnergySquared() const;
    void setCumulativeEnergySquared(double cumulativeEnergySquared);

private:
    int     m_numberOfMetropolisSteps  = 0;
    double     m_stepNumber               = 0;
    double     m_acceptedNumber           = 0;
    double  m_energy                   = 0;
    double  m_energySquared            = 0;
    double  m_cumulativeEnergySquared  = 0;
    double  m_cumulativeEnergy         = 0;
    double  m_cumulativeWF             = 0;
    double  m_WFderiv                  = 0;
    double  m_WFderivMultELoc          = 0;
    double  m_cumulativeWFderiv        = 0;
    double  m_cumulativeWFderivMultEloc= 0;

    std::string  m_filename;
    class System* m_system = nullptr;
};
