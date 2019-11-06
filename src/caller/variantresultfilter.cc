#include "variantresultfilter.h"

VariantResultFilter::VariantResultFilter()
{
}

bool VariantResultFilter::passFilterSV(Evidence *variantresult)
{
    if (!variantresult->isQuailtyPass())
    {
        return false;
    }

    //  std::cout << variantresult->getResultVcfFormatString() << std::endl;

    

    if (variantresult->getChr() == "")
    {
        return false;
    }

    if (variantresult->getPos() == 0)
    {
        return false;
    }

    // if (variantresult->getVariantType() == "INS")
    // {
    //     return passFilterInsertion(variantresult);
    // }

    // if (variantresult->getVariantType() == "DEL")
    // {
    //     return passFilterDeletion(variantresult);
    // }

    // if (variantresult->getVariantType() == "INV")
    // {
    //     return passFilterInversion(variantresult);
    // }

    // if (variantresult->getVariantType() == "DUP")
    // {
    //     return passFilterTandemDuplication(variantresult);
    // }

    // if (variantresult->getVariantType() == "BND")
    // {
    //     return passFilterTranslocation(variantresult);
    // }

    return true;
}

bool VariantResultFilter::passFilterDeletion(Evidence *variantresult)
{
    return true;
}

bool VariantResultFilter::passFilterTandemDuplication(Evidence *variantresult)
{
    // if (variantresult->getFrequency() < 8)
    // {
    //     return false;
    // }

    // if (variantresult->getMaxMapQ() < 40)
    // {
    //     return false;
    // }

    return true;
}

bool VariantResultFilter::passFilterInsertion(Evidence *variantresult)
{
    if (variantresult->getFrequency() <= 2)
    {
        return false;
    }

    return true;
}

bool VariantResultFilter::passFilterInversion(Evidence *variantresult)
{
    return true;
}

bool VariantResultFilter::passFilterTranslocation(Evidence *variantresult)
{
    return true;
}
