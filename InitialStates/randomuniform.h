#pragma once
#include "initialstate.h"

class RandomUniform : public InitialState {
public:
    RandomUniform(System* system, int numberOfDimensions, int numberOfParticles, double sigma, double interactionSize, double timeStep, int bins, double bucketSize);
    void setupInitialState();
  //  void setupInitialStateWithInteraction();

};

