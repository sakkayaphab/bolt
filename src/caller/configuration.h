#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "configuration.h"
#include <stdint.h>
class Configuration
{
private:
    

public:
    Configuration();

public:
    struct PE_DELETION
    {
        int32_t MININUM_READ_SMALL = 0;
        int32_t MININUM_READ_MEDIUM = 0;
        int32_t MININUM_READ_LARGE = 0;
        int8_t MININUM_MAPQ_SMALL = 0;
        int8_t MININUM_MAPQ_MEDIUM = 0;
        int8_t MININUM_MAPQ_LARGE = 0;
    };
    
};

#endif