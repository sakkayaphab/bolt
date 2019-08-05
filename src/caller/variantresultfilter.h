#ifndef VARIANTRESULTFILTER_H
#define VARIANTRESULTFILTER_H
#include "evidence.h"
class VariantResultFilter
{
private:
    
    
public:
    VariantResultFilter();
    bool passFilterSV(Evidence *variantresult);
    bool passFilterDeletion(Evidence *variantresult);
    bool passFilterTandemDuplication(Evidence *variantresult);
    bool passFilterInsertion(Evidence *variantresult);
    bool passFilterInversion(Evidence *variantresult);
    bool passFilterTranslocation(Evidence *variantresult);
};

#endif