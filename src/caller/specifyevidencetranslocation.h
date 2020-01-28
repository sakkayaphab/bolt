#ifndef SPECIFYINGEVIDENCETRANSLOCATION_H
#define SPECIFYINGEVIDENCETRANSLOCATION_H
#include "specifyevidence.h"
#include "evidence.h"

class SpecifyEvidenceTranslocation : public SpecifyEvidence
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
  SpecifyEvidenceTranslocation();
  void updateRead();
  void done();
};

#endif