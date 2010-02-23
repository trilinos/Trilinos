INCLUDE(TPLDeclareLibraries)

TPL_DECLARE_LIBRARIES( Peano
  REQUIRED_HEADERS SundanceCellInterface.h SundancePeanoInterface2D.h SundancePeanoInterface3D.h SundanceVertexInterface.h
  REQUIRED_LIBS_NAMES peano2D peano3D
  )
