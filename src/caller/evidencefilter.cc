#include "evidencefilter.h"

EvidenceFilter::EvidenceFilter()
{
}

bool EvidenceFilter::passFilterEvidence(Evidence *e)
{

    // if (e->getVariantType()=="DEL")
    // {
    //     return passDeletion(e);
    // }

    return true;
}

bool EvidenceFilter::passDeletion(Evidence *e)
{
    uint8_t maxMapQ;
    uint8_t avgMapQ;
    uint8_t minMapQ=255;
    int count;
    for (auto mapq:*e->getMapQVector())
    {
        count++; 
        avgMapQ+=mapq;
        
        if (maxMapQ<mapq)
        {
            maxMapQ = mapq;
        }

        if (minMapQ>mapq)
        {
            minMapQ = mapq;
        }
    }

    avgMapQ = avgMapQ/count;

    if (minMapQ==255)
    {
        minMapQ = 0;
    }

    if (avgMapQ==0)
    {
        return false;
    }else
    {
        return true;
    }

}