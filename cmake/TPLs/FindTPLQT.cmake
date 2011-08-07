INCLUDE(TPLDeclareLibraries)

SET(TPL_QT_QMAKE_EXECUTABLE "" CACHE STRING "A Trilinos specific variable that defines where the Qt Qmake Execuatable is")

IF(TPL_QT_QMAKE_EXECUTABLE AND QT_QMAKE_EXECUTABLE)
	IF(NOT(${TPL_QT_QMAKE_EXECUTABLE} EQUAL ${QT_QMAKE_EXECUTABLE}))
		MESSAGE(FATAL_ERROR "Uh oh. Looks like you set both the TPL_QT_QMAKE_EXECUTABLE and QT_QMAKE_EXECUTABLE variables and set them differently. You only need to set one.")
	ENDIF()
ENDIF()

IF(TPL_QT_QMAKE_EXECUTABLE)
	SET(QT_QMAKE_EXECUTABLE ${TPL_QT_QMAKE_EXECUTABLE})
ENDIF()

IF(QT_LIBRARY_DIRS)
	SET(QT_LIBRARY_DIR ${QT_LIBRARY_DIRS})
ENDIF()

IF(TPL_QT_LIBRARY_DIRS)
	SET(QT_LIBRARY_DIR ${TPL_QT_LIBRARY_DIRS})
ENDIF()

IF(QT_INCLUDE_DIRS)
	SET(QT_INCLUDE_DIR ${QT_INCLUDE_DIRS})
ENDIF()

IF(TPL_QT_INCLUDE_DIRS)
	SET(QT_INCLUDE_DIR ${TPL_QT_INCLUDE_DIRS})
ENDIF()

IF(TPL_QT_LIBRARIES)
	SET(QT_LIBRARIES ${TPL_QT_LIBRARIES})
ENDIF()


FIND_PACKAGE(Qt4 4.5.0 COMPONENTS QtCore QtGui QtTest QtXml REQUIRED)
if(NOT QT4_FOUND)
	message("                             ____")
	message("                     __,-~~/~    `---.")
	message("                   _/_,---(      ,    )")
	message("               __ /        <    /   )  \\___")
	message("- ------===;;;'====------------------===;;;===----- -  -")
	message("                 \\/  ~\"~\"~\"~\"~\"~\\~\"~)~\"/")
	message("                 (_ (   \\  (     >    \\)")
	message("                  \\_( _ <         >_>'")
	message("                      ~ `-i' ::>|--\"")
	message("                          I;|.|.|")
	message("                         <|i::|i|`.")
	message("                         (` ^'"`-' ")")
	MESSAGE(FATAL_ERROR "Couldn't find Qt4.5 or Greater. This causes explosions.")
endif()

IF(NOT(QT_INCLUDE_DIRS))
	SET(QT_INCLUDE_DIRS ${QT_INCLUDE_DIR} ${QT_QTCORE_INCLUDE_DIR} ${QT_QTGUI_INCLUDE_DIR})
ENDIF()

IF(NOT(QT_LIBRARY_DIRS))
	SET(QT_LIBRARY_DIRS ${QT_LIBRARY_DIR})
ENDIF()

IF(NOT(TPL_QT_INCLUDE_DIRS))
	SET(TPL_QT_INCLUDE_DIRS ${QT_INCLUDE_DIRS})
ENDIF()

IF(NOT(TPL_QT_LIBRARY_DIRS))
	SET(TPL_QT_LIBRARY_DIRS ${QT_LIBRARY_DIRS})
ENDIF()

IF(NOT(TPL_QT_LIBRARIES))
	SET(TPL_QT_LIBRARIES ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTTEST_LIBRARY} ${QT_QTXML_LIBRARY})
ENDIF()

TPL_DECLARE_LIBRARIES( QT
  REQUIRED_HEADERS QMainWindow QDialog QAbstratItemModel QTreeView QItemDelegate QPushButton QGridLayout QSpinBox QComboBox QLineEdit QLabel QScrollArea QDir QXmlStreamWriter QXmlStreamReader QStringList QDomElement
  REQUIRED_LIBS_NAMES QtCore QtGui QtXml
  )

