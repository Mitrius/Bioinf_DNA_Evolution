#include <string>
#include <vector>

#include "../headers/Phylogeny.h"

int main(int argc, char const *argv[])
{
    const double alpha = std::stod(argv[1]);
    const double beta = std::stod(argv[2]);
    const double generatorMean = std::stod(argv[3]);
    const int epochs = std::stod(argv[4]);

    std::vector<protein> initialSequence = {A,T,C,G};
    Phylogeny phylogeny(alpha,beta);

    std::vector<treeVertex> outputTree = phylogeny.phylogenesy(initialSequence,epochs,generatorMean,0.5); 
    return 0;
}
