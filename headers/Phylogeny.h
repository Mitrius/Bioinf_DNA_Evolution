#include <random>
#include <vector>
#include <map>

#include "treeVertex.h"

class Phylogeny
{
public:
  Phylogeny(double alpha, double beta);
  std::vector<treeVertex> phylogenesy(std::vector<treeVertex> &tree, int epochs, double timeGeneratorMean, double speciationProb);
  void writeBranchesIntoFile(std::string filename, std::vector<treeVertex> &tree);
  void printTree(std::vector<treeVertex> &tree);
  std::vector<treeVertex> reversePhylogeny(const std::string &filename);

private:
  double transversionProbability;
  double transpositionProbability;
  std::default_random_engine generator;

  base singleMutate(const base &initialProtein, double randValue);
  ProteinSequence mutate(ProteinSequence &initialSequence, std::default_random_engine &generator);
  std::string getReversedOrderOfAncestors(std::vector<treeVertex> &tree, treeVertex &vertex);
};