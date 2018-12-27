#include <random>
#include <vector>
#include <map>

#include "treeVertex.h"

class Phylogeny
{
public:
  Phylogeny(double alpha, double beta);
  std::vector<treeVertex> phylogenesy(std::vector<treeVertex> &tree, int epochs, double timeGeneratorMean);
  void writeBranchesIntoFile(std::string filename, std::vector<treeVertex> &tree);
  void printTree(std::vector<treeVertex> &tree);
  std::vector<treeVertex> reversePhylogeny(const std::string &filename);

private:
  double alpha;
  double beta;

  std::default_random_engine generator;
  std::vector<std::vector<double>> generateProbMatrix(double time);
  ProteinSequence mutate(ProteinSequence &initialSequence, double &time,std::default_random_engine& generator);
  std::string getReversedOrderOfAncestors(std::vector<treeVertex> &tree, treeVertex &vertex);
};