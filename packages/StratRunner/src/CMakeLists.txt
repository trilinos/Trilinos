INCLUDE(PackageAddExecutable)
INCLUDE(PackageLibraryMacros)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})


INCLUDE( FindQt4 )

IF( NOT QT4_FOUND)
	MESSAGE(SEND_ERROR "Could not find QT 4")
ENDIF( NOT QT4_FOUND)

SET( QT_USE_QTXML true)

INCLUDE( ${QT_USE_FILE} )

#INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/packages/epetra )

#APPEND_GLOBAL_SET( ${PACKAGE_NAME}_INCLUDE_DIRS
#${CMAKE_SOURCE_DIR}/packages/epetra/src/
#${CMAKE_SOURCE_DIR}/packages/thyra
#${CMAKE_SOURCE_DIR}/pacakges/teuchos
#${CMAKE_SOURCE_DIR}/packages/epetraext
#)

SET(MOC_HDRS metawindow.hpp runwindow.hpp solvertree.hpp startupwidget.hpp stratrunner.hpp)

SET(SRCS amesositem.cpp amesossettingsitem.cpp aztecooitem.cpp aztecoosettingsitem.cpp belositem.cpp 
belossolveritem.cpp blockcgitem.cpp blockgmresitem.cpp booltreeitem.cpp customstreambuff.cpp
doubletreeitem.cpp filenametreeitem.cpp freestringtreeitem.cpp ifpackitem.cpp ifpacksettingsitem.cpp main.cpp 
metawindow.cpp mlitem.cpp mlsettingsitem.cpp nopreconitem.cpp parameterlisttreeitem.cpp parametertreeitem.cpp
preconitem.cpp psuedoblockgmresitem.cpp runwindow.cpp solveritem.cpp solvertree.cpp startupwidget.cpp
stratroot.cpp stratrunner.cpp stringtreeitem.cpp verboseobjectitem.cpp )

QT4_WRAP_CPP(MOC_SRCS ${MOC_HDRS} )

PACKAGE_ADD_EXECUTABLE(
	StratRunner 
#	HEADERS ${MOC_SRCS} 
	SOURCES ${SRCS} ${MOC_SRCS} 
	DEPLIBS
	stratimikos
	#stratimikosamesos
	#stratimikosbelos
	#stratimikosifpack
	#stratimikosml
	#stratimikosaztecoo
	thyra
	teuchos
	epetra
	epetraext
	)

#ADD_EXECUTABLE(StratRunner ${SRCS} ${MOC_SRCS} )
#TARGET_LINK_LIBRARIES(StratRunner ${QT_LIBRARIES} )
#INSTALL(TARGETS StratRunner DESTINATION ${CMAKE_BINARY_DIR} )
