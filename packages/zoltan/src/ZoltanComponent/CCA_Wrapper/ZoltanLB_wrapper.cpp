#include "PartitionerFactory.h"
#include "CONV_NS.h"

extern "C" {

CONV_NS(Component) *create_PartitionerFactory_JR() 
{
        CONV_NS(Component) *wanker;
        ::ZoltanSpace::PartitionerFactory_JR *component;
        component = new ::ZoltanSpace::PartitionerFactory_JR();
        wanker = dynamic_cast< CONV_NS(Component) *>(component);
        return wanker;
}    

char **getComponentList() 
{
        static char *list[2];
        list[0] = "create_PartitionerFactory_JR PartitionerFactory_JR";
        list[1] = 0;
        return( list );
}
}
