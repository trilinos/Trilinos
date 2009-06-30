INCLUDE(PackageLibraryMacros)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

FIND_PACKAGE(Qt4)

IF( NOT QT4_FOUND)
	MESSAGE(SEND_ERROR "Could not find Qt 4")
ENDIF( NOT QT4_FOUND)
IF(QT_VERSION_MINOR LESS 5)
	MESSAGE(SEND_ERROR "You must have Qt 4.5 or higher installed")
ENDIF(QT_VERSION_MINOR LESS 5)


INCLUDE( ${QT_USE_FILE} )

SET(HDRS TivaBuena_Dependency.hpp 
TivaBuena_DependencySheet.hpp TivaBuena_ArrayHelperFunctions.hpp 
TivaBuena_GUI.hpp TivaBuena_InvalidDependencyException.hpp
TivaBuena_SpecificParameterEntryValidators.hpp TivaBuena_StandardDependencies.hpp
TivaBuena_treeitem.hpp TivaBuena_Types.hpp)


SET(MOC_HDRS TivaBuena_metawindow.hpp TivaBuena_treemodel.hpp
TivaBuena_delegate.hpp TivaBuena_ArrayWidget.hpp TivaBuena_FilenameWidget.hpp
TivaBuena_treeview.hpp)

SET(SRCS TivaBuena_delegate.cpp TivaBuena_Dependency.cpp 
TivaBuena_DependencySheet.cpp TivaBuena_ArrayHelperFunctions.cpp 
TivaBuena_FilenameWidget.cpp TivaBuena_GUI.cpp TivaBuena_metawindow.cpp 
TivaBuena_StandardDependencies.cpp TivaBuena_treeitem.cpp 
TivaBuena_treemodel.cpp TivaBuena_treeview.cpp TivaBuena_Types.cpp)
 
QT4_WRAP_CPP(MOC_SRCS ${MOC_HDRS} )

PACKAGE_ADD_LIBRARY(
	TivaBuena 
	HEADERS ${HDRS} 
	SOURCES ${SRCS} ${MOC_SRCS} 
	DEPLIBS
	teuchos ${QT_LIBRARIES}
	)

