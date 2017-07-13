

ADD_SUBDIRECTORIES(
  ParameterList # Put this first to define a test-only library
  Conditions
  Dependencies
  FancyOutputting
  FunctionObjects
  ObjectBuilder
  Validators
  XML
  XMLParser
  )
IF(TeuchosCore_ENABLE_yaml-cpp) 
  ADD_SUBDIRECTORIES(yaml)
ENDIF()
