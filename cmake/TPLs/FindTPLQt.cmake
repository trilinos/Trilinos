INCLUDE(TPLDeclareLibraries)

FIND_PACKAGE(Qt4 4.5.0 COMPONENTS QtCore QtGui)

TPL_DECLARE_LIBRARIES( Qt
  REQUIRED_HEADERS QMainWindow QDialog QAbstratItemModel QTreeView QItemDelegate QPushButton QGridLayout QSpinBox QComboBox QLineEdit QLabel QScrollArea QDir QXmlStreamWriter QXmlStreamReader QStringList
  REQUIRED_LIBS_NAMES QtCore QtGui
  LIBRARY_DIR_HINT ${QT_LIBRARY_DIR} 
  )

