#ifndef EVIDENCEFILTER_H
#define EVIDENCEFILTER_H
#include "evidence.h"

class EvidenceFilter
{
  private:
  public:
    EvidenceFilter();
    bool passFilterEvidence(Evidence *e);
    bool passDeletion(Evidence *e);
};

#endif