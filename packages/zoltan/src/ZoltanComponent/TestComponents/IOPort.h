#ifndef IOPortHSeen
#define  IOPortHSeen

// CCA standard definition of Port and Component
#include "cca.h"
#include "CONV_NS.h"

namespace ZoltanTestSpace
{
  class IOPort : public virtual CONV_NS(Port)
  {
    public :

   IOPort() {}

   virtual ~IOPort() {}

   virtual int ProduceGnuplotOutput(char *tag) = 0 ;

   virtual int ProduceFileOutput(char *tag) = 0 ;

   virtual int DeleteAllAndRestart() = 0 ;

   virtual int SetCmdFilename(char *) = 0 ;
  } ;
} ;
#endif
