#include "Mesh.h"
#include "CONV_NS.h"

extern "C" {

CONV_NS(Component) *create_Mesh() 
{
        CONV_NS(Component) *wanker;
        ::ZoltanTestSpace::Mesh *component;
        component = new ::ZoltanTestSpace::Mesh();
        wanker = dynamic_cast< CONV_NS(Component) *>(component);
        return wanker;
}    

char **getComponentList() 
{
        static char *list[2];
        list[0] = "create_Mesh Mesh";
        list[1] = 0;
        return( list );
}
}
