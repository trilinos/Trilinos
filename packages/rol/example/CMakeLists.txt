# Always built.
ADD_SUBDIRECTORY(rosenbrock)
ADD_SUBDIRECTORY(poisson-inversion)
ADD_SUBDIRECTORY(simple-eq-constr)
ADD_SUBDIRECTORY(diode-circuit)
ADD_SUBDIRECTORY(diode-circuit-simopt)
ADD_SUBDIRECTORY(topology-optimization)
ADD_SUBDIRECTORY(dual-spaces)
ADD_SUBDIRECTORY(zakharov)
ADD_SUBDIRECTORY(gross-pitaevskii)
ADD_SUBDIRECTORY(json)
ADD_SUBDIRECTORY(intrepid)
ADD_SUBDIRECTORY(minimax)
ADD_SUBDIRECTORY(wave-calibration)
ADD_SUBDIRECTORY(quadratic)
ADD_SUBDIRECTORY(arrayfire)
ADD_SUBDIRECTORY(complex)
ADD_SUBDIRECTORY(binary-design)
ADD_SUBDIRECTORY(paraboloid-circle)
ADD_SUBDIRECTORY(putting)
ADD_SUBDIRECTORY(PinT)
ADD_SUBDIRECTORY(dense-hessian)


IF(${PACKAGE_NAME}_ENABLE_Sacado)
ADD_SUBDIRECTORY(sacado)
ENDIF()

# Examples which require use of Teuchos
GET_PROPERTY( PTR_STRING GLOBAL PROPERTY PTR_IMPL )
GET_PROPERTY( PARAMETERLIST_STRING GLOBAL PROPERTY PARAMETERLIST_IMPL )

#IF( PTR_STRING STREQUAL "Teuchos::RCP" )
#IF( PARAMETERLIST_STRING STREQUAL "Teuchos::ParameterList" )
IF( ${USING_TEUCHOS_ALL} )

MESSAGE("-- Enabling ROL examples with required dependency on Teuchos::RCP or Teuchos::ParameterList")

ADD_SUBDIRECTORY(burgers-control)
ADD_SUBDIRECTORY(parabolic-control)
ADD_SUBDIRECTORY(poisson-control)
ADD_SUBDIRECTORY(tensor-opt)
ADD_SUBDIRECTORY(PDE-OPT)

IF(${PROJECT_NAME}_ENABLE_Tempus)
  ADD_SUBDIRECTORY(tempus)
ENDIF()

ENDIF( ${USING_TEUCHOS_ALL} )

#ENDIF()
#ENDIF()



