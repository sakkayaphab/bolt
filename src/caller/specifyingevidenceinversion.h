#ifndef DEFINITIONEVIDENCEINVERSION_H
#define DEFINITIONEVIDENCEINVERSION_H
#include "specifyingevidence.h"
#include "evidence.h"

class SpecifyingEvidenceInversion : public SpecifyingEvidence
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
  SpecifyingEvidenceInversion();
  void updateRead();
  void done();
};

#endif