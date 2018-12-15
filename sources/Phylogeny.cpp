#include <chrono>
#include <queue>
#include <fstream>
#include <iostream>
#include <map>

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
    std::queue<int> vertexQueue;

    treeVertex initialVertex(initialSequence);

    initialVertex.id = 0;
    initialVertex.root = -1;

    tree.push_back(initialVertex);
    vertexQueue.push(initialVertex.id);

    for (int epoch = 0; epoch < epochs; epoch++)
    {
        int vertexAmount = vertexQueue.size();
        for (int i = 0; i < vertexAmount; i++)
        {
            int vertex = vertexQueue.front();
            vertexQueue.pop();
            double prob = (distribution(generator) - distribution.min()) / distribution.max() - 1; // Squash value into 0-1 interval
            if (prob < speciationProb)
            {
                std::vector<protein> nextSpecie = mutate(tree[vertex].sequence, generator);
                treeVertex nextVertex(tree[vertex]);

                nextVertex.root = tree[vertex].id;
                nextVertex.sequence = nextSpecie;
                tree[vertex].right = nextVertex.id = tree.size();

                tree.push_back(nextVertex);
                vertexQueue.push(nextVertex.id);

                treeVertex dummyVertex(tree[vertex].sequence);
                tree[vertex].left = dummyVertex.id = tree.size();
                dummyVertex.root = tree[vertex].root;

                tree.push_back(dummyVertex);
                vertexQueue.push(dummyVertex.id);
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

void Phylogeny::writeBranchesIntoFile(std::string filename, std::vector<treeVertex> &tree)
{
    std::ofstream targetFile;
    std::vector<treeVertex> leaves;

    for (treeVertex &vertex : tree)
    {
        if (vertex.right == 0 && vertex.left == 0)
            leaves.push_back(vertex);
    }

    targetFile.open(filename, std::ios::out);
    
    if (targetFile.is_open())
    {

        for (treeVertex &vertex : leaves)
        {
            targetFile << this->getReversedOrderOfAncestors(tree, vertex) << "\n";
        }
    }
}

std::string Phylogeny::getReversedOrderOfAncestors(std::vector<treeVertex> &tree, treeVertex &vertex)
{
    std::vector<std::string> sequences;
    std::map<protein, std::string> codings;

    codings[A] = "A";
    codings[T] = "T";
    codings[C] = "C";
    codings[G] = "G";

    int nextVertex = vertex.id;
    std::vector<protein> sequence;
    std::vector<std::string> stringVector;

    while (nextVertex >= 0)
    {
        sequence = tree[nextVertex].sequence;

        std::string sequenceStringified;
        for (protein &prot : sequence)
        {
            sequenceStringified += codings[prot];
        }
        stringVector.push_back(sequenceStringified);

        nextVertex = tree[nextVertex].root;
    }

    std::string output;

    for (int i = stringVector.size() - 1; i > 0; i--)
    {
        output += stringVector[i] + ",";
    }
    output += stringVector[0];
    return output;
}
