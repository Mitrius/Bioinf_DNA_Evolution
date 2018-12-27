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
    this->alpha = alpha;
    this->beta = beta;
}
std::vector<std::vector<double>> Phylogeny::generateProbMatrix(double time)
{
    double s = (1 - std::exp(-4 * beta * time)) / 4;
    double u = (1 + std::exp(-4 * beta * time) - 2 * std::exp(-2 * (alpha + beta) * time)) / 4;
    double r = 1 - 2 * s - u;

    std::vector<std::vector<double>> retArray = {
        {r, s, u, s},
        {s, r, s, u},
        {u, s, r, s},
        {s, u, s, r}};

    return retArray;
}
std::vector<treeVertex> Phylogeny::phylogenesy(std::vector<treeVertex> &tree, int epochs, double timeGeneratorMean)
{
    int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::exponential_distribution<double> distribution(timeGeneratorMean);
    int actualEpoch = 0;
    std::queue<int> vertexQueue;

    for (treeVertex &vertex : tree)
    {
        if (vertex.left == 0 && vertex.right == 0) //only push leaves
            vertexQueue.push(vertex.id);
    }
    while (!vertexQueue.empty())
    {
        int vertex = vertexQueue.front();
        vertexQueue.pop();

        double time = tree[vertex].timeDepth + distribution(generator);

        ProteinSequence nextSpecie = mutate(tree[vertex].sequence, time, generator);

        if (nextSpecie.baseVec2String() == tree[vertex].sequence.baseVec2String())
            continue;
        treeVertex dummyVertex(tree[vertex].sequence);
        tree[vertex].left = dummyVertex.id = tree.size();

        dummyVertex.root = tree[vertex].id;
        dummyVertex.depth = tree[vertex].depth + 1;
        dummyVertex.timeDepth = time;

        tree.push_back(dummyVertex);

        treeVertex nextVertex(nextSpecie);

        nextVertex.root = tree[vertex].id;
        nextVertex.sequence = nextSpecie;
        nextVertex.depth = tree[vertex].depth + 1;
        nextVertex.timeDepth = time;

        tree[vertex].right = nextVertex.id = tree.size();

        tree.push_back(nextVertex);
        actualEpoch++;
        if (actualEpoch < epochs)
        {
            vertexQueue.push(dummyVertex.id);
            vertexQueue.push(nextVertex.id);
        }
    }
    return tree;
}
//TODO: mutatis mutandis
ProteinSequence Phylogeny::mutate(ProteinSequence &initialSequence, double &time, std::default_random_engine &generator)
{
    ProteinSequence mutatedSequence;

    std::vector<std::vector<double>> probMatrix = generateProbMatrix(time);
    double val;

    for (base &prot : initialSequence)
    {
        mutatedSequence.pushBack(prot);

        int protCode = static_cast<int>(prot);

        for (int i = 0; i < 4; i++)
        {
            if (i == protCode)
                continue;
            double prob = generator() / generator.max();
            if (prob < probMatrix[protCode][i])
            {
                mutatedSequence.popBack();
                mutatedSequence.pushBack(static_cast<base>(i));
                break;
            }
        }
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
                std::cout << vertex.sequence.baseVec2String() << " (d=" << vertex.time << "|id=" << vertex.id << "|r=" << vertex.root << ")   ";
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
                int length = std::stoi(lineSplitted.back()) + 1;
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