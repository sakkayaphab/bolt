#include "alignment.h"
#include <iostream>
#include "readparser.h"
#include <algorithm>

Alignment::Alignment(std::string *reference)
{
    setSVType(svtype);
    setReference(reference);
}

void Alignment::setReference(std::string *reference_m)
{
    reference = reference_m;
}

std::string *Alignment::getReference()
{
    return reference;
}

void Alignment::setSVType(std::string svtype_m)
{
    svtype = svtype_m;
}

void Alignment::genarateMatrix(std::string svtype)
{
    if (svtype == "DELEND")
    {
        SmithWaterman smw(getReference(), 0, true);
        swmMatrix = smw;
    }
    else if (svtype == "DELSTART")
    {
        SmithWaterman smw(getReference(), 0, false);
        swmMatrix = smw;
    }
    else if (svtype == "DUPSTART")
    {
        SmithWaterman smw(getReference(), 0, true);
        swmMatrix = smw;
    }
    else if (svtype == "DUPEND")
    {
        SmithWaterman smw(getReference(), 0, false);
        swmMatrix = smw;
    }
    else if (svtype == "TRASTART")
    {
        SmithWaterman smw(getReference(), 0, true);
        swmMatrix = smw;
        SmithWaterman smwSecond(getReference(), 0, false);
        swmMatrixSecond = smwSecond;
    }
    else if (svtype == "TRAEND")
    {
        SmithWaterman smw(getReference(), 0, true);
        swmMatrix = smw;
    }
    else if (svtype == "INVSTART")
    {
        SmithWaterman smw(getReference(), 0, true);
        swmMatrix = smw;
        SmithWaterman smwSecond(getReference(), 0, false);
        swmMatrixSecond = smwSecond;
    }
     else if (svtype == "INVEND")
    {
        SmithWaterman smw(getReference(), 0, true);
        swmMatrix = smw;
        SmithWaterman smwSecond(getReference(), 0, false);
        swmMatrixSecond = smwSecond;
    }
}

void Alignment::setPosReference(int32_t pos)
{
    posReference = pos;
}

int32_t Alignment::getPosReference()
{
    return posReference;
}

std::vector<SmithWaterman::ScoreAlignment> Alignment::alignDeletionTargetAtEnd(std::string *seq)
{
    // std::cout << *seq << std::endl;
    std::vector<SmithWaterman::ScoreAlignment> data = swmMatrix.findStartToEnd(*seq);
    std::vector<SmithWaterman::ScoreAlignment> result;
    for (auto n : data)
    {
        n.pos = n.pos+getPosReference();
        n.end = n.end+getPosReference();

        // std::cout << " n.pos : " << n.pos
        //           << " n.end : " << n.end
        //           << " n.posseq : " << n.posseq
        //           << " n.endseq : " << n.endseq
        //           << " n.pattern : " << n.pattern << std::endl;
        
    

        result.push_back(n);
    }

    return result;
}

std::vector<SmithWaterman::ScoreAlignment> Alignment::alignDeletionTargetAtStart(std::string *seq)
{
    // std::cout << *seq << std::endl;
    std::vector<SmithWaterman::ScoreAlignment> data = swmMatrix.findEndToStart(*seq);
     std::vector<SmithWaterman::ScoreAlignment> result;
    for (auto n : data)
    {
        
        n.pos = n.pos+getPosReference();
        n.end = n.end+getPosReference();

        // std::cout << " n.pos : " << n.pos
        //           << " n.end : " << n.end
        //           << " n.posseq : " << n.posseq
        //           << " n.endseq : " << n.endseq
        //           << " n.pattern : " << n.pattern << std::endl;

        result.push_back(n);
    }

    return result;
}

std::vector<SmithWaterman::ScoreAlignment> Alignment::alignDuplicationTargetAtStart(std::string *seq)
{
    return alignDeletionTargetAtEnd(seq);
}

std::vector<SmithWaterman::ScoreAlignment> Alignment::alignDuplicationTargetAtEnd(std::string *seq)
{
    return alignDeletionTargetAtStart(seq);
}

std::vector<SmithWaterman::ScoreAlignment> Alignment::alignTranslocationTargetAtStartSCS( std::string *seq)
{
    return alignDeletionTargetAtEnd(seq);
}

std::vector<SmithWaterman::ScoreAlignment> Alignment::alignTranslocationTargetAtStartSCE(std::string *seq)
{
    std::vector<SmithWaterman::ScoreAlignment> data = swmMatrixSecond.findEndToStart(*seq);
     std::vector<SmithWaterman::ScoreAlignment> result;
    for (auto n : data)
    {
        // n.pos = n.pos-1;
        // n.end = n.end-1;

        // int32_t end = getReference()->size()-n.pos;
        // int32_t pos = getReference()->size()-n.end;
        n.pos = n.pos+getPosReference();
        n.end = n.end+getPosReference();
        // std::cout << " n.pos : " << n.pos
        //           << " n.end : " << n.end
        //           << " n.posseq : " << n.posseq
        //           << " n.endseq : " << n.endseq
        //           << " n.pattern : " << n.pattern << std::endl;

        result.push_back(n);
    }

    return result;
//     return alignDeletionTargetAtStart(seq);
}

std::vector<SmithWaterman::ScoreAlignment> Alignment::alignTranslocationTargetAtEndSCS(std::string *seq)
{
    return alignDeletionTargetAtEnd(seq);
}

std::vector<SmithWaterman::ScoreAlignment> Alignment::alignTranslocationTargetAtEndSCE(std::string *seq)
{
    return alignTranslocationTargetAtStartSCE(seq);
}

std::vector<SmithWaterman::ScoreAlignment> Alignment::alignInversionTargetAtStartSCE( std::string *seq)
{
    ReadParser readparser;
    readparser.replaceToReverseComplement(seq);
    // reverse(seq->begin(), seq->end());
    std::vector<SmithWaterman::ScoreAlignment> result =  alignDeletionTargetAtEnd(seq);
    std::vector<SmithWaterman::ScoreAlignment> data;
    for (auto n:result) {
        int32_t endseq = seq->size()-n.posseq;
        int32_t posseq = seq->size()-n.endseq;
        n.posseq = posseq;
        n.endseq = endseq;

        // std::cout << " n.pos : " << n.pos
        //           << " n.end : " << n.end
        //           << " n.posseq : " << n.posseq
        //           << " n.endseq : " << n.endseq
        //           << " n.pattern : " << n.pattern << std::endl;
                  
        data.push_back(n);
    }

    return data;
}

std::vector<SmithWaterman::ScoreAlignment> Alignment::alignInversionTargetAtStartSCS(std::string *seq)
{
    ReadParser readparser;
    readparser.replaceToReverseComplement(seq);

    
    std::vector<SmithWaterman::ScoreAlignment> data = swmMatrixSecond.findEndToStart(*seq);
     std::vector<SmithWaterman::ScoreAlignment> result;
    for (auto n : data)
    {
        n.pos = n.pos+getPosReference()-1;
        n.end = n.end+getPosReference()-1;
        int32_t seqend = seq->size()- n.posseq;
        int32_t seqpos = seq->size()- n.endseq;

        n.posseq = seqpos;
        n.endseq = seqend;

        // std::cout << " n.pos : " << n.pos
        //           << " n.end : " << n.end
        //           << " n.posseq : " << n.posseq
        //           << " n.endseq : " << n.endseq
                  
        //           << " n.pattern : " << n.pattern << std::endl;

        result.push_back(n);
    }

    return result;
}

std::vector<SmithWaterman::ScoreAlignment> Alignment::alignInversionTargetAtEndSCE(std::string *seq)
{
    return alignInversionTargetAtStartSCE(seq);
}

std::vector<SmithWaterman::ScoreAlignment> Alignment::alignInversionTargetAtEndSCS(std::string *seq)
{
    return alignInversionTargetAtStartSCS(seq);
}