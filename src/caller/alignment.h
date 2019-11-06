#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "smithwaterman.h"

class Alignment
{
private:
    std::string svtype;
    std::string *reference;
    SmithWaterman swmMatrix;
    SmithWaterman swmMatrixSecond;
    int32_t posReference = 0;

public:
    Alignment(std::string *reference);
    void setReference(std::string *reference);
    void setSVType(std::string svtype);
    std::vector<SmithWaterman::ScoreAlignment> alignDeletionTargetAtEnd(std::string *seq);
    std::vector<SmithWaterman::ScoreAlignment> alignDeletionTargetAtStart(std::string *seq);

    std::vector<SmithWaterman::ScoreAlignment> alignDuplicationTargetAtStart(std::string *seq);
    std::vector<SmithWaterman::ScoreAlignment> alignDuplicationTargetAtEnd(std::string *seq);

    std::vector<SmithWaterman::ScoreAlignment> alignTranslocationTargetAtStartSCS(std::string *seq);
    std::vector<SmithWaterman::ScoreAlignment> alignTranslocationTargetAtStartSCE(std::string *seq);
    std::vector<SmithWaterman::ScoreAlignment> alignTranslocationTargetAtEndSCS(std::string *seq);
    std::vector<SmithWaterman::ScoreAlignment> alignTranslocationTargetAtEndSCE(std::string *seq);

    std::vector<SmithWaterman::ScoreAlignment> alignInversionTargetAtStartSCS(std::string *seq);
    std::vector<SmithWaterman::ScoreAlignment> alignInversionTargetAtStartSCE(std::string *seq);

    std::vector<SmithWaterman::ScoreAlignment> alignInversionTargetAtEndSCS(std::string *seq);
    std::vector<SmithWaterman::ScoreAlignment> alignInversionTargetAtEndSCE(std::string *seq);


    void replaceToReverseComplement(std::string *nucs);

    void genarateMatrix(std::string svtype);
    void setPosReference(int32_t pos);
    int32_t getPosReference();
    std::string *getReference();
};

#endif