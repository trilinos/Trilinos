// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2010) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%define %isorropia_epetra_docstring
"
PyTrilinos.Isorropia.Epetra is the python interface to namespace Epetra for
the Trilinos package Isorropia:

    http://trilinos.sandia.gov/packages/isorropia

The purpose of Isorropia.Epetra is to ....
"
%enddef

%module(package      = "PyTrilinos.Isorropia",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %isorropia_epetra_docstring) NestedEpetra

%{
// Configuration
#include "PyTrilinos_config.h"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_LocalMap.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_InvOperator.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_MapColoring.h"
#include "Epetra_SerialDistributor.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_OffsetIndex.h"
#include "Epetra_Time.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif

// Isorropia includes
#include "Isorropia_EpetraOperator.hpp"
#include "Isorropia_EpetraColorer.hpp"
#include "Isorropia_EpetraPartitioner.hpp"
#include "Isorropia_EpetraPartitioner2D.hpp"
#include "Isorropia_EpetraRedistributor.hpp"
#include "Isorropia_EpetraCostDescriber.hpp"
#include "Isorropia_EpetraOrderer.hpp"
#include "Isorropia_EpetraLevelScheduler.hpp"

// Local includes
#define NO_IMPORT_ARRAY
#define SWIG_FILE_WITH_INIT
#include "numpy_include.h"
#include "Epetra_NumPyFEVector.h"
#include "Epetra_NumPyIntSerialDenseMatrix.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPySerialSymDenseMatrix.h"
#include "Epetra_NumPySerialDenseVector.h"
#include "Epetra_NumPyVector.h"
#include "Teuchos_PythonParameter.h"
%}

// Configuration
%include "Epetra_DLLExportMacro.h"

// General exception handling
%include "exception.i"

%feature("director:except")
{
  if ($error != NULL)
  {
    throw Swig::DirectorMethodException();
  }
}

%exception
{
  try
  {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(Teuchos::Exceptions::InvalidParameterType & e)
  {
    SWIG_exception(SWIG_TypeError, e.what());
  }
  catch(Teuchos::Exceptions::InvalidParameter & e)
  {
    PyErr_SetString(PyExc_KeyError, e.what());
    SWIG_fail;
  }
  catch(Swig::DirectorException &e)
  {
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// Include Isorropia documentation (same as for Isorropia.__init__.i)
%include "Isorropia_dox.i"

// General ignore directives
%ignore *::operator=;
%ignore *::operator<<;

// Teuchos interface import
%import "Teuchos.i"

// numpy interface import
%include "numpy.i"


// Epetra interface import
%import "Epetra.i"

// Isorropia import (let SWIG know about the base classes that will be
// needed for the derived classes below)
%import "Isorropia.__init__.i"

/////////////////////////////////////////
// Isorropia::Epetra::Operator support //
/////////////////////////////////////////
%teuchos_rcp(Isorropia::Epetra::Operator)
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* elementList, int len)};
%apply (int DIM1, int* ARGOUT_ARRAY1) {(int len, int* array)};
%include "Isorropia_EpetraOperator.hpp"

////////////////////////////////////////
// Isorropia::Epetra::Colorer support //
////////////////////////////////////////
%teuchos_rcp(Isorropia::Epetra::Colorer)
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* elementList, int len)};
%apply (int DIM1, int* ARGOUT_ARRAY1) {(int len, int* array)};
%include "Isorropia_EpetraColorer.hpp"

////////////////////////////////////////////
// Isorropia::Epetra::Partitioner support //
////////////////////////////////////////////
%teuchos_rcp(Isorropia::Epetra::Partitioner)
%ignore Isorropia::Epetra::Partitioner::createNewMap(Epetra_Map *&);
%include "Isorropia_EpetraPartitioner.hpp"

//////////////////////////////////////////////
// Isorropia::Epetra::Partitioner2D support //
//////////////////////////////////////////////
//%teuchos_rcp(Isorropia::Epetra::Partioner2D)
//%include "Isorropia_EpetraPartitioner2D.hpp"

//////////////////////////////////////////////
// Isorropia::Epetra::Redistributor support //
//////////////////////////////////////////////
//%teuchos_rcp(Isorropia::Epetra::Redistributor)
%ignore Isorropia::Epetra::Redistributor::redistribute(Epetra_CrsMatrix const &, Epetra_CrsMatrix *&, bool);
%ignore Isorropia::Epetra::Redistributor::redistribute(Epetra_RowMatrix const &, Epetra_CrsMatrix *&, bool);
%ignore Isorropia::Epetra::Redistributor::redistribute(Epetra_Vector const &, Epetra_Vector *&);
%ignore Isorropia::Epetra::Redistributor::redistribute(Epetra_MultiVector const &, Epetra_MultiVector *&);
%include "Isorropia_EpetraRedistributor.hpp"

//////////////////////////////////////////////
// Isorropia::Epetra::CostDescriber support //
//////////////////////////////////////////////
%teuchos_rcp(Isorropia::Epetra::CostDescriber)
%apply (int DIM1, float* ARGOUT_ARRAY1) {(int len, float* weights)};
%apply (int DIM1, int* ARGOUT_ARRAY1) {(int len, int* global_ids)};
%include "Isorropia_EpetraCostDescriber.hpp"

////////////////////////////////////////
// Isorropia::Epetra::Orderer support //
////////////////////////////////////////
%teuchos_rcp(Isorropia::Epetra::Orderer)
%apply (int DIM1, int* ARGOUT_ARRAY1) {(int len, int* array)};
%include "Isorropia_EpetraOrderer.hpp"

///////////////////////////////////////////////
// Isorropia::Epetra::LevelScheduler support //
///////////////////////////////////////////////
%teuchos_rcp(Isorropia::Epetra::LevelScheduler)
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* elementList, int len)};
%include "Isorropia_EpetraLevelScheduler.hpp"


// Isorropia Visualizer Code
%pythoncode
%{
try:
    import sys
    import random
    import colorsys
    import math

    try:
        from PIL import Image, ImageDraw
    except ImportError:
        import Image, ImageDraw

    Image_loaded = True
except ImportError:
    Image_loaded = False 

class Visualizer:
    """
    The Visualizer class provides a means of simply visualizing partitions
    formed by Isorropia.Epetra.Partitioner or Partitioner2D. It can also draw
    Epetra.CrsMatrix objects directly.
    
    The Visualizer relies on the Python Imaging Library and ImageMagick.
    
    To generate images of multiple sizes given the same partitioner, you'll
    need to call generateImage() several times, calling getPILImage() after
    each call and saving the result. You can then save or show those images as
    necessary.

    """

    data     = None  # The matrix/partition we are visualizing. After
                     # redistribution this will always be a CrsMatrix.
    dataType = None  # The type of the original data passed in.
                     # Currently supported: CrsMatrix, Partitioner,
                     # Partitioner2D
    
    numProcs        = 0
    verbosity       = 0
    iAmRoot         = False
    loadThreshold   = 0

    img             = None

    background      = 0x000000
    matrixHeight    = -1
    matrixWidth     = -1
    imageHeight     = -1
    imageWidth      = -1
        
    pointTL         = (0,0)
    pointBR         = (-1,-1)
    regionHeight    = -1
    regionWidth     = -1

    def __init__(self, data, numParts=-1, loadThreshold=0, verbosity=0):
        """
        Visualizer constructor.
        Inputs: data      - CrsMatrix or Partitioner
                numParts  - The number of parts to distribute into.
                            Defaults to the number of processes in MPI.
                loadThreshold - The threshold ratio that must be exceeded
                            for the visualizer to display something in that
                            pixel. Zero is a safe choice.
                verbosity - Print debugging output. High number increases
                            quantity of output.

       The initialization first reassigns the values in the matrix to represent
       the processor they reside on (or have been partitioned to), then
       redistributes the matrix onto the root processor for use in show().
        """

        if Image_loaded == False:
            print "Python Imaging Library, ImageMagick and colorsys libraries required for Visualizer function"
            sys.exit(1)

        if isinstance(data, PyTrilinos.Epetra.CrsMatrix):
            self.comm = data.Comm()
        elif isinstance(data, PyTrilinos.Isorropia.Epetra.Partitioner) or isinstance(data, PyTrilinos.Isorropia.Epetra.Partitioner2D):
            self.comm = data.InputMatrix().Comm()
        else:
            sys.exit(1)

        iAmRoot = self.iAmRoot = self.comm.MyPID() == 0

        if numParts == -1:
            self.numProcs = comm.NumProc()
        else:
            self.numProcs = numParts

        self.verbosity = verbosity
        self.data = data
        self.loadThreshold = loadThreshold

        if isinstance(data, PyTrilinos.Epetra.CrsMatrix):
            if iAmRoot and verbosity > 1: print "Proceeding with CrsMatrix"
            
            self.dataType = "CrsMatrix"
            
            if iAmRoot and verbosity > 2: print "Redistributing for visualization"
        
            # Reset the values of the matrix to their respective processors
            self.centralizeMatrix()
        
        elif isinstance(data, PyTrilinos.Isorropia.Epetra.Partitioner):
            if iAmRoot and verbosity > 1: print "Proceeding with Partitioner"
            
            self.dataType = "Partitioner"

            if iAmRoot and verbosity > 2: print "Redistributing for visualization"
        
            # Reset the values of the matrix to their partitioned processors
            self.centralizePartitioner()

        elif isinstance(data, PyTrilinos.Isorropia.Epetra.Partitioner2D):
            if iAmRoot and verbosity > 1: print "Proceeding with Partitioner2D"
            
            self.dataType = "Partitioner2D"

            if iAmRoot and verbosity > 2: print "Redistributing for visualization"
        
            # Reset the values of the matrix to their partitioned processors
            self.centralizePartitioner2D()

        # Only support CrsMatrix and Partitioner objects currently
        else:
            sys.exit(1)


    def processRow(self, row, draw):
        """
        processRow() is a private helper to draw an individual row of pixels of
        the image when simpleGenerate() can't be used.

        """

        if not self.iAmRoot: return

        if self.verbosity > 3: print "Processing row " + str(row)

        widthRatio = float(self.regionWidth) / float(self.imageWidth)
        heightRatio = float(self.regionHeight) / float(self.imageHeight)

        # Calculate starting and ending positions in the matrix which will be 
        # drawn into this row of pixels
        startY = int(math.ceil(row * heightRatio)) + self.pointTL[1]
        endY = min(int(math.floor((row + 1) * heightRatio)) + self.pointTL[1], self.matrixHeight - 1)
        #startY = int(math.ceil(row * heightRatio))
        #endY = min(int(math.floor((row + 1) * heightRatio)), self.matrixHeight - 1)

        # Extract the rows to be represented
        nonZeroes = []
        for j in range(startY, endY + 1):
            nonZeroes += [self.data.ExtractGlobalRowCopy(j)]

        # procCount: List of lists. One list per pixel, whose values represent
        # the number of nonzeroes for the i-th processor found in this pixel so far
        procCount = [[0]*(self.numProcs + 1) for i in range(self.imageWidth)]
        for k in range(len(nonZeroes)):
                                                # x = arbitrary id. Ideally would be dictionary.
            procNums = nonZeroes[k][0]          # Maps x -> processor ID
            colIndices = nonZeroes[k][1]        # Maps x -> column index.

            # Go through the nonzeroes in this row and increment the appropriate
            # counter
            for l in range(len(procNums)): # (len(colIndices) would be the same)
                # Only show the columns we care about
                if colIndices[l] < self.pointTL[0]:
                    continue
                if colIndices[l] >= self.pointBR[0]:
                    break

                # Take the matrix column and show us which pixel column it should go into
                pixelColumn = int((colIndices[l] - self.pointTL[0]) / widthRatio)
                if (pixelColumn >= len(procCount)):
                    print "Invalid value found for pixelColumn. Quitting."
                    print pixelColumn, len(procCount)
                    sys.exit(1)

                # Increment the counter for whomever owns that nonzero
                if (int(procNums[l]) >= len(procCount[pixelColumn])) or (int(procNums[l]) < 0):
                    print "Invalid value found for part value. Quitting."
                    sys.exit(1)
                else:
                    procCount[pixelColumn][int(procNums[l])] += 1

        # Draw each pixel in the row, once we have counted up the nonzeroes
        for col in range(self.imageWidth):
            # Ensure there were enough nonzeroes in this pixel to warrant drawing.
            # Currently the threshold is 0, so any nonzero is displayed, but it
            # supports change.
            numElementsInPixel = int(widthRatio * (endY - startY))
            if numElementsInPixel == 0: loadRatio = 0
            else: loadRatio = float(sum(procCount[col])) / numElementsInPixel
            if(max(procCount[col]) != 0 and loadRatio > self.loadThreshold):
                draw.rectangle([col, row,
                                col + 1,
                                row + 1],
                               fill=self.fillColors[procCount[col].index(max(procCount[col])) - 1])

    
    def generateImage(self):
        """
        Opens an image window representing the matrix stored in this object's
        data field. This opened window allows saving the image.

        The colors chosen to represent each processor are evenly distributed
        along the hue and presently can not be customized.

        """

        if not self.iAmRoot: return

        if self.dataType != "CrsMatrix" and self.dataType != "Partitioner" and self.dataType != "Partitioner2D":
            raise "Error: Only CrsMatrix and Partitioner[2D] data formats supported at present"
            return

        self.matrixHeight = self.data.NumGlobalRows()
        self.matrixWidth  = self.data.NumGlobalCols()
        
        if self.pointBR[0] == -1 or self.pointBR[1] == -1:
            self.pointBR = (self.matrixWidth, self.matrixHeight)

        self.regionHeight = self.pointBR[1] - self.pointTL[1]
        self.regionWidth  = self.pointBR[0] - self.pointTL[0]

        # Generate an image that maintains aspect ratio if we have not been given a height.
        if self.imageHeight == -1:
            self.imageHeight = int(float(self.regionHeight) / self.regionWidth * self.imageWidth)
            if self.verbosity > 3: print "Defining height as", self.imageHeight

        # Create Image and ImageDraw objects
        im = Image.new("RGB", (self.imageWidth, self.imageHeight), self.background)
        draw = ImageDraw.Draw(im)


        if (self.pointBR[0] > self.matrixWidth or
            self.pointBR[1] > self.matrixWidth or
            self.pointBR[0] < self.pointTL[0] or
            self.pointBR[1] < self.pointTL[1] or
            self.pointTL[0] < 0 or
            self.pointTL[1] < 0 or
            self.regionHeight <= 0 or
            self.regionWidth <= 0):
            raise Exception("Invalid region of interest")
        
        cellWidth = self.imageWidth / self.regionWidth
        cellHeight = self.imageHeight / self.regionHeight

        self.fillColors = self.fillColorArray(self.numProcs)

        """ If we have more matrix rows or columns than we have pixels, can't use simpleGenerate() """
        if cellWidth == 0 and cellHeight == 0:    
            if self.verbosity > 2: print "NOT using simpleGenerate()"
            for i in range(self.imageHeight): # Go through the rows
                self.processRow(i, draw)
        elif cellWidth == 0 or cellHeight == 0:
            sys.stderr.write("""
                ERROR: There exists a known bug where nothing prints if you try to
                print part of your matrix such that you have more rows in your
                image than in your matrix but fewer columns than in your
                matrix. Or vice versa. We're working on fixing this. In the
                meantime, please re-size your print region, your image output,
                or both!
                """)
        else:
            if self.verbosity > 2: print "Using simpleGenerate()"
            self.simpleGenerate(draw)

        
        self.img = im

    def show(self):
        """
        Shows the pre-generated image. Must be preceded by a call to generateImage().
        """
        if not self.iAmRoot: return

        if self.img == None:
            print "Image must be computed before call to show(). Call generateImage() first."
            return
        
        self.img.show()

    def save(self, filename):
        """
        Saves the pre-generated image to the given file. The format is derived
        from the filename. Must be preceded by a call to generateImage().
        
        """
        if not self.iAmRoot: return

        if self.img == None:
            print "Image must be computed before call to save(). Call generateImage() first."
            return
        
        self.img.save(filename)

    def getPILImage(self):
        """
        Returns the PIL Image object generated by generateImage(). Returns None
        if no image has yet been generated.
        """
        return self.img

    def simpleGenerate(self, draw):
        """
        Draws the matrix to the provided image object, provided that the width
        and height of the image exceed the number of columns and rows of the matrix.
        
        Used by generateImage().
        """
        if not self.iAmRoot: return
        
        regionWidth = self.pointBR[0] - self.pointTL[0]
        regionHeight = self.pointBR[1] - self.pointTL[1]

        fillColors = self.fillColorArray(self.numProcs)

        cellWidth = float(self.imageWidth) / regionWidth
        cellHeight = float(self.imageHeight) / regionHeight

        # Loop through rows of matrix
        for i in range(self.pointTL[1], self.pointBR[1]):
            nonZeroes = self.data.ExtractGlobalRowCopy(i)
                    
            procNums = nonZeroes[0]
            colIndices = nonZeroes[1]
   
            startY = (i - self.pointTL[1]) * cellHeight
            
            # Loop through nonzeroes, drawing appropriate rectangles
            for j in range(len(colIndices)):
                if colIndices[j] < self.pointTL[0]: continue
                elif colIndices[j] >= self.pointBR[0]: break
                startX = (colIndices[j] - self.pointTL[0]) * cellWidth
                draw.rectangle([int(startX), int(startY),
                                int(startX + cellWidth),
                                int(startY + cellHeight)],
                               fill=fillColors[int(procNums[j]) - 1])
        
    def setBgColor(self, background):
        self.background = background
    
    def setImageHeight(self, height):
        self.imageHeight = height

    def setImageWidth(self, width):
        self.imageWidth = width

    def setViewingWindow(self, pointTL, pointBR=(-1,-1)):
        self.pointTL = pointTL
        self.pointBR = pointBR

    def setVerbosity(self, verbosity):
        self.verbosity = verbosity
    
    def fillColorArray(self, numParts):
        """
        Given the number of partitions, returns an array of at least that many
        distinct color values
        """

        HSV_tuples = [(x*1.0/numParts, 0.9, 0.8) for x in range(numParts)]
        RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
        hex_tuples = map(lambda x: (int(255*x[2]) << 16)+(int(255*x[1]) << 8)+(int(255*x[0])), RGB_tuples)

        random.shuffle(hex_tuples)

        return hex_tuples


    def centralizeMatrix(self):
        """
        Performs two operations:
            1. Replaces all non-zero values in the matrix with a value
               representing that non-zero's global part ID
            2. Repartitions such that all elements in the matrix are owned by
               the main process.
        """

        """ Reset the values in each row of the matrix to this processor's ID """
        for i in range(self.data.NumMyRows()):
            rowView = self.data.ExtractMyRowCopy(i) # Values first, indices second
            self.data.ReplaceMyValues(i, [self.comm.MyPID() + 1]*len(rowView[1]), rowView[1])

        globalNumRows = self.data.NumGlobalRows()

        # Create a map to move all rows to the root processor
        if self.iAmRoot:
            localNumRows = globalNumRows
        else:
            localNumRows = 0
       
        targetMap = PyTrilinos.Epetra.Map(globalNumRows, localNumRows, 0, self.comm)
        redist = PyTrilinos.Isorropia.Epetra.Redistributor(targetMap)
        # Redistribute to root processor
        self.data = redist.redistribute(self.data)
    
    def centralizePartitioner(self):
        """
        Performs three operations:
            1. Replaces the partitioner with a new CrsMatrix, maintaining
               non-zero distribution (with the same RowMap).
            1. Replaces all non-zero values in the matrix with a value
               representing that non-zero's global part ID
            2. Repartitions such that all elements in the matrix are owned by
               the main process.

        """

        rowMatrix = self.data.InputMatrix()
        map = rowMatrix.RowMatrixRowMap()

        # Get information needed to create a duplicate CrsMatrix
        rowLengths = []
        rows = []
        for i in range(rowMatrix.NumMyRows()):
            rowView = rowMatrix.ExtractMyRowCopy(i) # Values first, indices second
            rows += [rowView[1]]
            rowLengths += [len(rowView[1])]
        
        # Create a new CrsMatrix to duplicate the Partitioner's RowMatrix
        crsm = PyTrilinos.Epetra.CrsMatrix(PyTrilinos.Epetra.Copy, map, rowMatrix.RowMatrixColMap(), rowLengths)
        for gid in range(map.MaxLID()+1):
            values = [1]*len(rows[gid])
            crsm.InsertMyValues(gid, values, rows[gid])
        crsm.FillComplete()

        # Reassign the CrsMatrix's values to their mapped processor IDs
        newMap = self.data.createNewMap()

        matrixHeight = rowMatrix.NumGlobalRows()

        for i in range(crsm.NumMyRows()):
            rowView = rowMatrix.ExtractMyRowCopy(i) # Values first, indices second
            gid = map.GID(i)
            newLoc = self.data.index(i) ###########################################################
            if newLoc < 0 or newLoc > matrixHeight:
                print "newLoc is wrong, quitting", newLoc
                sys.exit(1)
            crsm.ReplaceMyValues(i, [newLoc + 1]*len(rowView[1]), rowView[1])
        crsm.FillComplete()

        # Set self.data to the CrsMatrix representation
        self.data = crsm

        globalNumRows = self.data.NumGlobalRows()

        # Create a map to move all rows to the root processor
        if self.iAmRoot:
            localNumRows = globalNumRows
        else:
            localNumRows = 0
       
        targetMap = PyTrilinos.Epetra.Map(globalNumRows, localNumRows, 0, self.comm)
        redist = PyTrilinos.Isorropia.Epetra.Redistributor(targetMap)
        # Redistribute to root processor
        self.data = redist.redistribute(self.data)

    def centralizePartitioner2D(self):
        """
        Exactly the same as centralizePartitioner(), but accepting a Partitioner2D.

        Creates a new CrsMatrix to represent the Partitioner2D's matrix
        structure. This new CrsMatrix has the same distribution as the
        Partitioner2D's matrix, but the values of the nonzeroes represent their
        'destination' processors in the Partitioner's new Map
        """
        
        rowMatrix = self.data.InputMatrix()
        map = rowMatrix.RowMatrixRowMap()

        # Get information needed to create a duplicate CrsMatrix
        rowLengths = []
        rows = []
        for i in range(rowMatrix.NumMyRows()):
            rowView = rowMatrix.ExtractMyRowCopy(i) # Values first, indices second
            rows += [rowView[1]]
            rowLengths += [len(rowView[1])]
        
        # Create a new CrsMatrix to duplicate the Partitioner's RowMatrix
        crsm = PyTrilinos.Epetra.CrsMatrix(PyTrilinos.Epetra.Copy, map, rowMatrix.RowMatrixColMap(), rowLengths)
        for gid in range(map.MaxLID()+1):
            values = [0]*len(rows[gid])
            crsm.InsertMyValues(gid, values, rows[gid])
        crsm.FillComplete()
        
        # Reassign the CrsMatrix's values to their mapped processor IDs
        matrixHeight = rowMatrix.NumGlobalRows()

        numMyNZs = rowMatrix.NumMyNonzeros()
        counter = 0

        for i in range(crsm.NumMyRows()):
            rowView = rowMatrix.ExtractMyRowCopy(i) # Values first, indices second
            gid = map.GID(i)
            newProcs = [0] * len(rowView[1])
            for j in range(len(newProcs)):
                newProcs[j] = self.data.index(counter)
                counter += 1
            crsm.SumIntoMyValues(i, newProcs, rowView[1])

        crsm.FillComplete()

        # Set self.data to the CrsMatrix representation
        self.data = crsm

        globalNumRows = self.data.NumGlobalRows()

        # Create a map to move all rows to the root processor
        if self.iAmRoot:
            localNumRows = globalNumRows
        else:
            localNumRows = 0
       
        targetMap = PyTrilinos.Epetra.Map(globalNumRows, localNumRows, 0, self.comm)
        redist = PyTrilinos.Isorropia.Epetra.Redistributor(targetMap)
        # Redistribute to root processor
        self.data = redist.redistribute(self.data)
        
%}
