#ifndef SPECIFYINGEVIDENCETRANSLOCATION_H
#define SPECIFYINGEVIDENCETRANSLOCATION_H
#include "specifyingevidence.h"
#include "evidence.h"

class SpecifyingEvidenceTranslocation : public SpecifyingEvidence
{
private:
  void checkRange();

  uint32_t currentPos = 0;
  uint32_t currentMPos = 0;
  void proveEvidence(int index);
  void checkProveEvidence();
  bool filterEvidence(Evidence *evidence);
  void filterEvidenceFinal();
  void calculateVCF(Evidence *evidence);
  bool incrementSVFreq(int32_t overlappedpos,int32_t overlappedsvlength, int32_t pos, int32_t mpos);
  
protected:
public:
  SpecifyingEvidenceTranslocation();
  void updateRead();
  void done();
};

#endif