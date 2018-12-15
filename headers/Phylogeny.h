#include <random>
#include <vector>

#include "treeVertex.h"

class Phylogeny
{
  public:
    Phylogeny(double alpha, double beta);
    std::vector<treeVertex> phylogenesy(std::vector<protein> &initialSequence, int epochs, double timeGeneratorMean, double speciationProb);

  private:
    double transversionProbability;
    double transpositionProbability;
    std::default_random_engine generator;

    protein singleMutate(const protein &initialProtein, double randValue);
    std::vector<protein> mutate(const std::vector<protein> &initialSequence, std::default_random_engine &generator);
};