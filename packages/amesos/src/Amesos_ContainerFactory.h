#ifndef AMESOS_CONTAINTERFACTORY_H
#define AMESOS_CONTAINTERFACTORY_H

#include "Amesos_ConfigDefs.h"
class Amesos_Container;

class Amesos_ContainerFactory {

public:

  Amesos_Container* Create(const string Type);


};  

#endif // AMESOS_CONTAINTERFACTORY_H
