#ifndef SIERRA_Ioss_EntityType_H
#define SIERRA_Ioss_EntityType_H

namespace Ioss {
  enum EntityType {NODEBLOCK    =   1,
		   ELEMENTBLOCK =   2,
		   NODESET      =   4,
		   EDGESET      =   8,
		   FACESET      =  16,
		   SURFACE      =  16, //: Same as faceset
		   COMMSET      =  32,
		   EDGEBLOCK    =  64,
		   FACEBLOCK    = 128,
                   REGION       = 256,
                   SUPERELEMENT = 512};
}
#endif
