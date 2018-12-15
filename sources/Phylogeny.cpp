#include <chrono>
#include <queue>
#include "../headers/Phylogeny.h"

Phylogeny::Phylogeny(double alpha, double beta)
{
    this->transversionProbability = beta;
    this->transpositionProbability = alpha;
}
std::vector<treeVertex> Phylogeny::phylogenesy(std::vector<protein> &initialSequence, int epochs, double timeGeneratorMean, double speciationProb)
{
    int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::exponential_distribution<double> distribution(timeGeneratorMean);

    std::vector<treeVertex> tree;
    std::queue<treeVertex> vertexQueue;

    treeVertex initialVertex(initialSequence);

    initialVertex.id = 0;
    initialVertex.root = -1;

    tree.push_back(initialVertex);
    vertexQueue.push(initialVertex);

    for (int epoch = 0; epoch < epochs; epoch++)
    {
        int vertexAmount = vertexQueue.size();
        for (int i = 0; i < vertexAmount; i++)
        {
            treeVertex vertex = vertexQueue.front();
            vertexQueue.pop();
            double prob = (distribution(generator) - distribution.min()) / distribution.max() - 1; // Squash value in 0-1 interval
            if (prob < speciationProb)
            {
                std::vector<protein> nextSpecie = mutate(vertex.sequence, generator);
                treeVertex nextVertex(vertex);

                nextVertex.root = vertex.id;
                nextVertex.sequence = nextSpecie;
                vertex.right = nextVertex.id = tree.size();

                tree.push_back(nextVertex);
                vertexQueue.push(nextVertex);

                treeVertex dummyVertex(vertex.sequence);
                vertex.left = dummyVertex.id = tree.size();
                dummyVertex.root = vertex.id;

                tree.push_back(dummyVertex);
                vertexQueue.push(dummyVertex);
            }
            else
            {
                vertexQueue.push(vertex);
            }
        }
    }
    return tree;
}
protein Phylogeny::singleMutate(const protein &initialProtein, double randValue)
{
    const int protId = static_cast<int>(initialProtein);

    if (randValue <= transversionProbability) // Transversion occurs
        return static_cast<protein>((protId + 1) % 4);
    else if (randValue <= transpositionProbability) // Transposition occurs
        return static_cast<protein>((protId + 2) % 4);
    else //nothing changes
        return initialProtein;
}

std::vector<protein> Phylogeny::mutate(const std::vector<protein> &initialSequence, std::default_random_engine &generator)
{
    std::vector<protein> mutatedSequence;
    double val;

    for (const protein &prot : initialSequence)
    {
        val = static_cast<double>(generator()) / generator.max();
        mutatedSequence.push_back(singleMutate(prot, val));
    }

    return mutatedSequence;
}
