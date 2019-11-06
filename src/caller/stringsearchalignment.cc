#include "stringsearchalignment.h"
#include <iostream>

StringSearchAlignment::StringSearchAlignment()
{
}

void StringSearchAlignment::setReference(std::string reference)
{
    StringSearchAlignment::reference = reference;
}

void StringSearchAlignment::setSVType(std::string svtype_m)
{
    svtype = svtype_m;
}

void StringSearchAlignment::setPosReference(int32_t pos)
{
    posReference = pos;
}

int32_t StringSearchAlignment::getPosReference()
{
    return posReference;
}

void StringSearchAlignment::buildReference()
{
    ss.setReference(reference.c_str());
    ss.buildHashTable();
}

std::vector<StringSearch::Score> StringSearchAlignment::alignDeletionTargetAtStart(std::string *seq,StringSearchConfig *ssc)
{
    std::vector<StringSearch::Score> scoreresult = ss.searchEndToStart(seq,ssc);
    std::vector<StringSearch::Score> result;
    for (auto m : scoreresult)
    {
        m.pos += getPosReference();
        m.end += getPosReference();
        result.push_back(m);
        // std::cout
        //     << "pos : " << m.pos << "\t"
        //     << "end : " << m.end << "\t"
        //     << "posseq : " << m.posseq << "\t"
        //     << "endseq : " << m.endseq << "\t"
        //     << "pattern : ";

        // for (auto n : m.matchSeqPattern)
        // {
        //     std::cout << n;
        // }

        // std::cout << std::endl;
    }

    return result;
}

std::vector<StringSearch::Score> StringSearchAlignment::alignDeletionTargetAtEnd(std::string *seq,StringSearchConfig *ssc)
{
    std::vector<StringSearch::Score> scoreresult = ss.searchStartToEnd(seq,ssc);
    std::vector<StringSearch::Score> result;
    for (auto m : scoreresult)
    {
        m.pos += getPosReference();
        m.end += getPosReference();
        result.push_back(m);
        // std::cout
        //     << "pos : " << m.pos << "\t"
        //     << "end : " << m.end << "\t"
        //     << "posseq : " << m.posseq << "\t"
        //     << "endseq : " << m.endseq << "\t"
        //     << "pattern : ";

        // for (auto n : m.matchSeqPattern)
        // {
        //     std::cout << n;
        // }

        // std::cout << std::endl;
    }

    return result;
}

std::vector<StringSearch::Score> StringSearchAlignment::alignDuplicationTargetAtStart(std::string *seq,StringSearchConfig *ssc)
{
    return alignDeletionTargetAtEnd(seq,ssc);
}

std::vector<StringSearch::Score> StringSearchAlignment::alignDuplicationTargetAtEnd(std::string *seq,StringSearchConfig *ssc)
{
    return alignDeletionTargetAtStart(seq,ssc);
}

std::vector<StringSearch::Score> StringSearchAlignment::alignTranslocationTargetAtStartSCS(std::string *seq,StringSearchConfig *ssc)
{
    return alignDeletionTargetAtEnd(seq,ssc);
}

std::vector<StringSearch::Score> StringSearchAlignment::alignTranslocationTargetAtStartSCE(std::string *seq,StringSearchConfig *ssc)
{
    return alignDeletionTargetAtStart(seq,ssc);
}
std::vector<StringSearch::Score> StringSearchAlignment::alignTranslocationTargetAtEndSCS(std::string *seq,StringSearchConfig *ssc)
{
    return alignDeletionTargetAtEnd(seq,ssc);
}
std::vector<StringSearch::Score> StringSearchAlignment::alignTranslocationTargetAtEndSCE(std::string *seq,StringSearchConfig *ssc)
{
    return alignDeletionTargetAtStart(seq,ssc);
}

std::vector<StringSearch::Score> StringSearchAlignment::alignInversionTargetAtStartSCS(std::string *seq,StringSearchConfig *ssc)
{
    ReadParser readparser;
    readparser.replaceToReverseComplement(seq);
    std::vector<StringSearch::Score> result;
    std::vector<StringSearch::Score> scorelist = alignDeletionTargetAtStart(seq,ssc);
    for (auto m : scorelist)
    {
        int32_t tempPos = seq->size() - m.endseq;
        int32_t tempEnd = seq->size() - m.posseq;
        m.posseq = tempPos;
        m.endseq = tempEnd;
        result.push_back(m);

        // std::cout
        //     << "pos : " << m.pos << "\t"
        //     << "end : " << m.end << "\t"
        //     << "posseq : " << m.posseq << "\t"
        //     << "endseq : " << m.endseq << "\t"
        //     << "pattern : ";

        // for (auto n : m.matchSeqPattern)
        // {
        //     std::cout << n;
        // }
    }

    return result;
}

std::vector<StringSearch::Score> StringSearchAlignment::alignInversionTargetAtStartSCE(std::string *seq,StringSearchConfig *ssc)
{
    ReadParser readparser;
    readparser.replaceToReverseComplement(seq);
    std::vector<StringSearch::Score> result;
    std::vector<StringSearch::Score> scorelist = alignDeletionTargetAtEnd(seq,ssc);
    for (auto m : scorelist)
    {
        int32_t tempPos = seq->size() - m.endseq;
        int32_t tempEnd = seq->size() - m.posseq;
        m.posseq = tempPos;
        m.endseq = tempEnd;
        result.push_back(m);

        // std::cout
        //     << "pos : " << m.pos << "\t"
        //     << "end : " << m.end << "\t"
        //     << "posseq : " << m.posseq << "\t"
        //     << "endseq : " << m.endseq << "\t"
        //     << "pattern : ";

        // for (auto n : m.matchSeqPattern)
        // {
        //     std::cout << n;
        // }
    }

    return result;
}

std::vector<StringSearch::Score> StringSearchAlignment::alignInversionTargetAtEndSCE(std::string *seq,StringSearchConfig *ssc)
{
    return alignInversionTargetAtStartSCE(seq,ssc);
}

std::vector<StringSearch::Score> StringSearchAlignment::alignInversionTargetAtEndSCS(std::string *seq,StringSearchConfig *ssc)
{
    return alignInversionTargetAtStartSCS(seq,ssc);
}