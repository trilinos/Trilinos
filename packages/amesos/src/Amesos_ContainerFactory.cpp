#include "Amesos_ConfigDefs.h"
#include "Amesos_ContainerFactory.h"
#include "Amesos_Container.h"
#include "Amesos_ContainerLAPACK.h"
#include "Amesos_ContainerEpetraCrs.h"

Amesos_Container* Amesos_ContainerFactory::Create(const string Type)
{

  if (Type == "LAPACK") {
    return(new Amesos_ContainerLAPACK);
  }
  else if (Type == "EpetraCrs") {
    return(new Amesos_ContainerEpetraCrs);
  }

  return(0);

}
			 

