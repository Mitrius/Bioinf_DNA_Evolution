#include <string>
#include <vector>

#include "../headers/Phylogeny.h"

int main(int argc, char const *argv[])
{
    const double alpha = std::stod(argv[1]);
    const double beta = std::stod(argv[2]);

    Phylogeny phylogeny(alpha, beta);

    const double generatorMean = std::stod(argv[3]);
    const double speciationProbability = std::stod(argv[4]);
    const int epochs = std::stoi(argv[5]);

    const std::string sequence = std::string(argv[6]);

    ProteinSequence initialSequence(sequence);
    treeVertex initialVertex(initialSequence);
    initialVertex.root = -1;

    std::vector<treeVertex> initialTree;
    initialTree.push_back(initialVertex);

    std::vector<treeVertex> outputTree = phylogeny.phylogenesy(initialTree, epochs, generatorMean, speciationProbability);

    phylogeny.printTree(outputTree);

    phylogeny.writeBranchesIntoFile("tree.csv", outputTree);

    std::vector<treeVertex> loadedTree = phylogeny.reversePhylogeny("tree.csv");

    phylogeny.phylogenesy(loadedTree, 1, generatorMean, speciationProbability);
    
    phylogeny.printTree(loadedTree);
    return 0;
}
