#include <string>
#include <vector>

#include "../headers/Phylogeny.h"

int main(int argc, char const *argv[])
{
    const double alpha = std::stod(argv[1]);
    const double beta = std::stod(argv[2]);
    const double generatorMean = std::stod(argv[3]);
    const double speciationProbability = std::stod(argv[4]);
    const int epochs = std::stod(argv[5]);


    std::vector<protein> initialSequence = {A,T,C,G};
    Phylogeny phylogeny(alpha,beta);

    std::vector<treeVertex> outputTree = phylogeny.phylogenesy(initialSequence,epochs,generatorMean,speciationProbability); 
    phylogeny.writeBranchesIntoFile("tree.csv",outputTree);
    return 0;
}
