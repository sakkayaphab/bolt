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

    if (variantresult->getVariantType() == "INS")
    {
        return passFilterInsertion(variantresult);
    }

    if (variantresult->getChr() == "")
    {
        return false;
    }

    if (variantresult->getPos() == 0)
    {
        return false;
    }

    if (variantresult->getVariantType() == "DEL")
    {
        return passFilterDeletion(variantresult);
    }

    if (variantresult->getVariantType() == "INV")
    {
        return passFilterInversion(variantresult);
    }

    if (variantresult->getVariantType() == "DUP")
    {
        return passFilterTandemDuplication(variantresult);
    }

    if (variantresult->getVariantType() == "BND")
    {
        return passFilterTranslocation(variantresult);
    }

    return false;
}

bool VariantResultFilter::passFilterDeletion(Evidence *variantresult)
{
    // if (variantresult->getEnd() - variantresult->getPos() < 20)
    // {
    //     return false;
    // }

    // if (variantresult->getEnd() - variantresult->getPos() > 50000)
    // {
    //     return false;
    // }

    // if (variantresult->getSvLength() < 1000)
    // {
    //     // if (variantresult->getMaxMapQ() == 0)
    //     // {
    //     //     return false;
    //     // }

    //     // if (variantresult->getM() < 3)
    //     // {
    //     //     return false;
    //     // }

    //     // if (variantresult->LNGMATCH>100) {
    //     //     return false;
    //     // }
    // }
    // else
    // {
    // }

    // if (variantresult->getFrequency() > 40)
    // {
    //     return false;
    // }

    // if (variantresult->getMaxMapQ() > 45 && variantresult->getFrequency() >= 3)
    // {
    //     return true;
    // }

    // if (variantresult->getMaxMapQ() >= 0 && variantresult->getFrequency() >= 6)
    // {
    //     return true;
    // }

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
    if (variantresult->getFrequency() < 2)
    {
        return false;
    }

    return true;
}

bool VariantResultFilter::passFilterInversion(Evidence *variantresult)
{
    // if (variantresult->getMaxMapQ() < 10)
    // {
    //     return false;
    // }

    // if (variantresult->getFrequency() < 3)
    // {
    //     return false;
    // }

    return true;
}

bool VariantResultFilter::passFilterTranslocation(Evidence *variantresult)
{
    return true;
}
