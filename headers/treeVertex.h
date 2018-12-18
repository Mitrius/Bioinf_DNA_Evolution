#include "ProteinSequence.h"

struct treeVertex
{
    ProteinSequence sequence;
    int left = 0;
    int right = 0;
    int id = 0;
    int root = 0;
    int depth = 0;

    treeVertex(ProteinSequence &sequence)
    {
        this->sequence = sequence;
    }
    treeVertex() {}
};