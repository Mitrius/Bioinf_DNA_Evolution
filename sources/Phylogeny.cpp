#include <chrono>
#include <queue>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include "../headers/Phylogeny.h"

Phylogeny::Phylogeny(double alpha, double beta)
{
    transversionProbability = beta;
    transpositionProbability = alpha;
}
std::vector<treeVertex> Phylogeny::phylogenesy(ProteinSequence &initialSequence, int epochs, double timeGeneratorMean, double speciationProb)
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
            double prob = std::fmod(distribution(generator), 1.0); // Squash value into 0-1 interval
            if (prob <= speciationProb)
            {
                ProteinSequence nextSpecie = mutate(tree[vertex].sequence, generator);
                if (nextSpecie.baseVec2String() == tree[vertex].sequence.baseVec2String())
                    continue;
                treeVertex dummyVertex(tree[vertex].sequence);
                tree[vertex].left = dummyVertex.id = tree.size();
                dummyVertex.root = tree[vertex].id;
                dummyVertex.depth = tree[vertex].depth + 1;

                tree.push_back(dummyVertex);
                vertexQueue.push(dummyVertex.id);

                treeVertex nextVertex(nextSpecie);
                nextVertex.root = tree[vertex].id;
                nextVertex.sequence = nextSpecie;
                nextVertex.depth = tree[vertex].depth + 1;

                tree[vertex].right = nextVertex.id = tree.size();

                tree.push_back(nextVertex);
                vertexQueue.push(nextVertex.id);
            }
            else
            {
                vertexQueue.push(vertex);
            }
        }
    }
    return tree;
}
base Phylogeny::singleMutate(const base &initialProtein, double randValue)
{
    const int protId = static_cast<int>(initialProtein);

    if (randValue <= transversionProbability) // Transversion occurs
        return static_cast<base>((protId + 1) % 4);
    else if (randValue <= transpositionProbability) // Transposition occurs
        return static_cast<base>((protId + 2) % 4);
    else //nothing changes
        return initialProtein;
}

ProteinSequence Phylogeny::mutate(ProteinSequence &initialSequence, std::default_random_engine &generator)
{
    ProteinSequence mutatedSequence;
    double val;

    for (base &prot : initialSequence)
    {
        val = static_cast<double>(generator()) / generator.max();
        mutatedSequence.pushBack(singleMutate(prot, val));
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

    int nextVertex = vertex.id;
    ProteinSequence sequence;
    std::vector<std::string> stringVector;
    while (nextVertex >= 0)
    {
        sequence = tree[nextVertex].sequence;

        stringVector.push_back(sequence.baseVec2String());

        nextVertex = tree[nextVertex].root;
    }

    std::string output;

    for (int i = stringVector.size() - 1; i >= 0; i--)
    {
        output += stringVector[i] + ",";
    }
    output += std::to_string(vertex.depth);
    return output;
}
void Phylogeny::printTree(std::vector<treeVertex> &tree)
{
    treeVertex maxDepthVertex = *std::max_element(tree.begin(), tree.end(), [](treeVertex &a, treeVertex &b) { return a.depth < b.depth; });
    for (int depth = 0; depth <= maxDepthVertex.depth; depth++)
    {

        for (treeVertex &vertex : tree)
        {
            if (vertex.depth == depth)
            {
                std::cout << vertex.sequence.baseVec2String() << " (d=" << vertex.depth << "|id=" << vertex.id << "|r=" << vertex.root << ")   ";
            }
        }
        std::cout << '\n'
                  << '\n';
    }
}

std::vector<treeVertex> Phylogeny::reversePhylogeny(const std::string &filename)
{
    std::ifstream inputStream;
    inputStream.open(filename);
    std::vector<treeVertex> tree;

    if (inputStream.is_open())
    {
        std::vector<std::vector<std::string>> readTree;
        // read to map
        std::string line;
        std::vector<std::string> lineSplitted;
        int maxLength = 0;
        do
        {
            std::getline(inputStream, line);
            boost::split(lineSplitted, line, boost::is_any_of(","));
            if (lineSplitted.back() != "")
            {
                int length = std::stoi(lineSplitted.back()) +1;
                if (length > maxLength)
                    maxLength = length;
            }
            lineSplitted.pop_back();
            readTree.push_back(lineSplitted);
            lineSplitted.clear();

        } while (!inputStream.eof());
        readTree.pop_back();
        // set root

        treeVertex rootVertex;
        rootVertex.sequence = ProteinSequence(readTree[0][0]);
        rootVertex.id = 0;
        rootVertex.root = -1;
        rootVertex.depth = 0;

        tree.push_back(rootVertex);

        for (int i = 1; i < maxLength; i++)
        {
            for (std::vector<std::string> &stringVec : readTree)
            {

                for (int epoch = 1; epoch <= i && epoch < stringVec.size(); epoch++)
                {
                    std::string sequence = stringVec[epoch];
                    std::string ancestor = stringVec[epoch - 1];
                    std::vector<int> possibleAncestors;
                    for (int i = 0; i < tree.size(); i++)
                    {
                        if (tree[i].sequence.baseVec2String() == ancestor && tree[i].depth == epoch - 1)
                            possibleAncestors.push_back(i);
                    }
                    for (int vertex : possibleAncestors) // possible immediate ancestors, now we need to find way to the root
                    {
                        bool goodWay = true;
                        int iter = 1;
                        int nextVertex = vertex;
                        while (goodWay && tree[nextVertex].root != -1 && epoch > iter)
                        {
                            if (tree[nextVertex].sequence.baseVec2String() != stringVec[epoch - iter])
                            {
                                goodWay = false;
                                break;
                            }
                            iter++;
                            nextVertex = tree[nextVertex].root;
                        }
                        if (!goodWay)
                            continue;
                        else
                        {
                            treeVertex newVertex;
                            newVertex.sequence = ProteinSequence(sequence);
                            newVertex.id = tree.size();
                            newVertex.root = tree[vertex].id;
                            newVertex.depth = epoch;

                            ancestor = tree[vertex].sequence.baseVec2String();
                            if (sequence == ancestor) //left child (no left child present)
                            {
                                if (tree[vertex].left == 0)
                                {
                                    tree.push_back(newVertex);
                                    tree[vertex].left = newVertex.id;
                                }
                                else
                                {
                                    break;
                                }
                            }
                            else if (tree[vertex].right == 0) //right child
                            {
                                tree.push_back(newVertex);
                                tree[vertex].right = newVertex.id;
                            }

                            break;
                        }
                    }
                }
            }
        }
    }
    return tree;
}