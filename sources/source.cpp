#include <string>
#include <vector>

#include "../headers/Phylogeny.h"

int main(int argc, char const *argv[])
{
    const double alpha = std::stod(argv[1]);
    const double beta = std::stod(argv[2]);
    const double generatorMean = std::stod(argv[3]);
    const int epochs = std::stoi(argv[4]);
    const std::string sequence = std::string(argv[5]);

    Phylogeny phylogeny(alpha, beta);   

    ProteinSequence initialSequence(sequence);
    treeVertex initialVertex(initialSequence);
    initialVertex.root = -1;

    std::vector<treeVertex> initialTree;
    initialTree.push_back(initialVertex);

    std::vector<treeVertex> outputTree = phylogeny.phylogenesy(initialTree, epochs,generatorMean);

    phylogeny.printTree(outputTree);

    phylogeny.writeBranchesIntoFile("tree.csv", outputTree);

    std::vector<treeVertex> loadedTree = phylogeny.reversePhylogeny("tree.csv");

    phylogeny.phylogenesy(loadedTree, 1, generatorMean);
    
    phylogeny.printTree(loadedTree);
    return 0;
}
