#ifndef STRINGSEARCHALIGNMENT_H
#define STRINGSEARCHALIGNMENT_H
#include <string>
#include "stringsearch.h"
#include "readparser.h"
#include "stringsearchconfig.h"

class StringSearchAlignment
{
private:
    ;
    std::string svtype;
    std::string reference;
    StringSearch ss;
    int32_t posReference = 0;

public:
    StringSearchAlignment();
    void setReference(std::string reference);
    void setSVType(std::string svtype);
    void buildReference();
    void setPosReference(int32_t pos);
    int32_t getPosReference();

    std::vector<StringSearch::Score> alignDeletionTargetAtStart(std::string *seq,StringSearchConfig *ssc);
    std::vector<StringSearch::Score> alignDeletionTargetAtEnd(std::string *seq,StringSearchConfig *ssc);
    std::vector<StringSearch::Score> alignDuplicationTargetAtStart(std::string *seq,StringSearchConfig *ssc);
    std::vector<StringSearch::Score> alignDuplicationTargetAtEnd(std::string *seq,StringSearchConfig *ssc);
    std::vector<StringSearch::Score> alignTranslocationTargetAtStartSCS(std::string *seq,StringSearchConfig *ssc);
    std::vector<StringSearch::Score> alignTranslocationTargetAtStartSCE(std::string *seq,StringSearchConfig *ssc);
    std::vector<StringSearch::Score> alignTranslocationTargetAtEndSCS(std::string *seq,StringSearchConfig *ssc);
    std::vector<StringSearch::Score> alignTranslocationTargetAtEndSCE(std::string *seq,StringSearchConfig *ssc);
    std::vector<StringSearch::Score> alignInversionTargetAtStartSCS(std::string *seq,StringSearchConfig *ssc);
    std::vector<StringSearch::Score> alignInversionTargetAtStartSCE(std::string *seq,StringSearchConfig *ssc);

     std::vector<StringSearch::Score> alignInversionTargetAtEndSCE(std::string *seq,StringSearchConfig *ssc);
     std::vector<StringSearch::Score> alignInversionTargetAtEndSCS(std::string *seq,StringSearchConfig *ssc);
};

#endif