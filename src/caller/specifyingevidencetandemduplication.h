#ifndef DEFINITIONEVIDENCETANDEMDUPLICATION_H
#define DEFINITIONEVIDENCETANDEMDUPLICATION_H
#include "specifyingevidence.h"
#include "evidence.h"

class SpecifyingEvidenceTandemDuplication : public SpecifyingEvidence
{
private:
  void checkRange();

  uint32_t currentPos = 0;
  uint32_t currentMPos = 0;
  void proveEvidence(int index);
  void checkProveEvidence();
  bool filterEvidence(Evidence *evidence);
  void calculateVCF(Evidence *evidence);
protected:
public:
  SpecifyingEvidenceTandemDuplication();
  void updateRead();
  void done();
};

#endif