SET(SacadoExample_SOURCES
  example.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  SacadoExample 
  SOURCES ${SacadoExample_SOURCES}
  )
