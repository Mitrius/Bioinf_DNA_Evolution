#include "../headers/ProteinSequence.h"
ProteinSequence::ProteinSequence()
{
    codings[A] = 'A';
    codings[T] = 'T';
    codings[C] = 'C';
    codings[G] = 'G';

    deCodings['A'] = A;
    deCodings['T'] = T;
    deCodings['C'] = C;
    deCodings['G'] = G;
}
ProteinSequence::ProteinSequence(const std::string &initialSequence)
{
    codings[A] = 'A';
    codings[T] = 'T';
    codings[C] = 'C';
    codings[G] = 'G';

    deCodings['A'] = A;
    deCodings['T'] = T;
    deCodings['C'] = C;
    deCodings['G'] = G;

    for (const char &base : initialSequence)
    {
        inside.push_back(deCodings[base]);
    }
}
std::string ProteinSequence::baseVec2String()
{
    std::string sequenceStringified;
    for (base &prot : inside)
    {
        sequenceStringified += codings[prot];
    }
    return sequenceStringified;
}
base &ProteinSequence::operator[](int index)
{
    return inside[index];
}

void ProteinSequence::pushBack(const base &newBase)
{
    inside.push_back(newBase);
}
void ProteinSequence::popBack()
{
    inside.pop_back();
}