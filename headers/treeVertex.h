#include <vector>
#include "protein.h"

struct treeVertex
{
    std::vector<protein> sequence;
    int left = 0;
    int right = 0;
    int id = 0;
    int root = 0;

    treeVertex(std::vector<protein> &sequence)
    {
        this->sequence = sequence;
    }
};