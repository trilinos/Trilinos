# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER

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

