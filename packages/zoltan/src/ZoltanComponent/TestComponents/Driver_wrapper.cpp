#include "Driver.h"

#include "CONV_NS.h"

extern "C" {

CONV_NS(Component) *create_Driver() 
{
        CONV_NS(Component) *wanker;
        ::ZoltanTestSpace::Driver *component;
        component = new ::ZoltanTestSpace::Driver();
        wanker = dynamic_cast< CONV_NS(Component) *>(component);
        return wanker;
}    

char **getComponentList() 
{
        static char *list[2];
        list[0] = "create_Driver Driver";
        list[1] = 0;
        return( list );
}
}
