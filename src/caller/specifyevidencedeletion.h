#ifndef SPECIFYINGEVIDENCEDELETION_H
#define SPECIFYINGEVIDENCEDELETION_H
#include "specifyevidence.h"
#include "evidence.h"

class SpecifyEvidenceDeletion : public SpecifyEvidence
{
private:
  int32_t currentPos = 0;
  int32_t currentMPos = 0;

  void checkRange();
  void checkErrorEvidence();
  void proveEvidence(int index);
  bool filterEvidence(Evidence *evidence);
  //  void addErrorToEvidence();
  void calculateVCF(Evidence *evidence);
  bool incrementSVFreq(int32_t overlappedpos,int32_t overlappedsvlength, int32_t pos, int32_t mpos);
protected:
public:
  SpecifyEvidenceDeletion();
  void updateRead();
  void done();
  void checkProveEvidence();
  int32_t getSVLength();
  void removeDuplicateFinalEvidence();
  
};

#endif