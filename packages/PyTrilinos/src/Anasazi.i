// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

// This documentation string will be the python help facility help
// string
%define %anasazi_docstring
"
PyTrilinos.Anasazi is the python interface to the Trilinos eigensolver
package Anasazi:

    http://trilinos.sandia.gov/packages/anasazi

Anasazi is a collection of eigensolver technologies.  The C++ version
makes liberal use of templates to abstract out the scalar type,
multivector type and operator type.  Currently, the python version
supports only an Epetra interface, so that the scalar type is float
(C++ double), the multivector type is Epetra.MultiVector and the
operator type is Epetra.Operator.  The following python factories
produce python wrappers to underlying C++ instances of the given C++
names with concrete Epetra implementations:

    * SortManager
    * BasicSort
    * Eigenproblem
    * BasicEigenproblem
    * StatusTest
    * StatusTestCombo
    * StatusTestMaxIters
    * StatusTestOutput
    * StatusTestResNorm
    * OrthoManager
    * MatOrthoManager
    * SVQBOrthoManager
    * Eigensolver
    * SolverManager
    * BlockDavidson
    * BlockDavidsonSolMgr
    * BlockKrylovSchur
    * BlockKrylovSchurSolMgr
    * LOBPCG
    * LOBPCGSolMgr
    * Eigensolution

For an example of usage of Anasazi to solve an eigenproblem, see the
following script in the example subdirectory of the PyTrilinos
package:

    * exAnasazi_BlockDavidson.py
"
%enddef

// Define the module name, its package and documentation string
%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %anasazi_docstring) Anasazi

%{
// System include files
#include <ostream>
#include <sstream>

// Configuration include files
#include "PyTrilinos_config.h"
#include "AnasaziConfigDefs.hpp"

#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// Local include files
#include "PyTrilinos_FILEstream.hpp"

// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"

// Epetra include files
#ifdef HAVE_PYTRILINOS_EPETRA
#include "PyTrilinos_Epetra_Headers.hpp"
#endif

// Anasazi include files
#include "PyTrilinos_Anasazi_Headers.hpp"
%}

// Include configuration
%include "PyTrilinos_config.h"
%include "Anasazi_config.h"
#ifdef HAVE_PYTRILINOS_EPETRA
%include "Anasaziepetra_DLLExportMacro.h"
#endif

// Standard exception handling
%include "exception.i"

// Auto-documentation feature
%feature("autodoc", "1");

// Include Anasazi documentation
%include "Anasazi_dox.i"

// C++ STL support
%include "stl.i"

// General ignore directives
%ignore *::operator=;
%ignore *::print;

// General rename directive
%rename(_global) global;

// Support for other Trilinos packages
%import "Teuchos.i"

#ifdef HAVE_PYTRILINOS_EPETRA
%include "numpy.i"
%import "Epetra.i"
#endif

// General exception handling
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

//////////////////////////////////////////////
// Support these classes, encapsulated in a //
// Teuchos::RCP<...>, as function arguments //
//////////////////////////////////////////////
%teuchos_rcp(std::basic_ostream)
%teuchos_rcp(std::ostream)
%teuchos_rcp(std::vector<int, std::allocator<int> >)
%teuchos_rcp(Anasazi::MultiVec< double >)
%teuchos_rcp(Anasazi::OutputManager< double >)
%teuchos_rcp(Anasazi::BasicOutputManager< double >)
#ifdef HAVE_PYTRILINOS_EPETRA
%teuchos_rcp(Anasazi::EpetraMultiVecAccessor)
%teuchos_rcp(Anasazi::EpetraMultiVec)
%teuchos_rcp(Anasazi::EpetraGenOp)
%teuchos_rcp(Anasazi::EpetraSymOp)
%teuchos_rcp(Anasazi::EpetraOp)
%teuchos_rcp(Anasazi::EpetraSymMVOp)
%teuchos_rcp(Anasazi::EpetraWSymMVOp)
%teuchos_rcp(Anasazi::EpetraW2SymMVOp)
%teuchos_rcp(Anasazi::StatusTest< double, Epetra_MultiVector, Epetra_Operator >)
%teuchos_rcp(Anasazi::StatusTestCombo< double, Epetra_MultiVector, Epetra_Operator >)
%teuchos_rcp(Anasazi::StatusTestMaxIters< double, Epetra_MultiVector, Epetra_Operator >)
%teuchos_rcp(Anasazi::StatusTestOutput< double, Epetra_MultiVector, Epetra_Operator >)
%teuchos_rcp(Anasazi::StatusTestResNorm< double, Epetra_MultiVector, Epetra_Operator >)
%teuchos_rcp(Anasazi::SortManager< double >)
%teuchos_rcp(Anasazi::BasicSort< double >)
%teuchos_rcp(Anasazi::Eigenproblem< double, Epetra_MultiVector, Epetra_Operator >)
%teuchos_rcp(Anasazi::BasicEigenproblem< double, Epetra_MultiVector, Epetra_Operator >)
%teuchos_rcp(Anasazi::OrthoManager< double, Epetra_MultiVector >)
%teuchos_rcp(Anasazi::MatOrthoManager< double, Epetra_MultiVector, Epetra_Operator >)
%teuchos_rcp(Anasazi::BasicOrthoManager< double, Epetra_MultiVector, Epetra_Operator >)
%teuchos_rcp(Anasazi::SVQBOrthoManager< double, Epetra_MultiVector, Epetra_Operator >)
#endif

/////////////////////////////////////////////////////////////////////////
// Anasazi returns eigenvalues in a std::vector< Anasazi::Value<       //
// ScalarType > > container.  This is support for converting them to a //
// more convenient numpy array.                                        //
/////////////////////////////////////////////////////////////////////////
%define %anasazi_eigenvalues_typemap(ScalarType, NumPyType)
%typemap(out,fragment="NumPy_Backward_Compatibility,NumPy_Macros")
  (std::vector< Anasazi::Value< ScalarType > >)
{
  npy_intp dims[1] = { static_cast< npy_intp >($1.size()) };
  PyObject * array = PyArray_SimpleNew(1, dims, NumPyType);
  ScalarType * data = (ScalarType*) array_data(array);
  for (npy_intp i=0; i<dims[0]; ++i)
  {
    data[2*i  ] = $1[i].realpart;
    data[2*i+1] = $1[i].imagpart;
  }
  return array;
}
%enddef
%anasazi_eigenvalues_typemap(float , NPY_CFLOAT )
%anasazi_eigenvalues_typemap(double, NPY_CDOUBLE)

/////////////////////////////
// Anasazi Version support //
/////////////////////////////
%include "AnasaziVersion.cpp"
%pythoncode
%{
__version__ = Anasazi_Version().split()[2]
%}

///////////////////////////
// Anasazi Types support //
///////////////////////////
%extend Anasazi::Eigensolution
{
  std::vector< Anasazi::Value< ScalarType > > evals()
  {
    return self->Evals;
  }
  Teuchos::RCP<MV> evecs()
  {
    return self->Evecs;
  }
  Teuchos::RCP<MV> espace()
  {
    return self->Espace;
  }
}
%ignore Anasazi::Eigensolution::Evals;
%ignore Anasazi::Eigensolution::Evecs;
%ignore Anasazi::Eigensolution::Espace;
%include "AnasaziTypes.hpp"

%extend Anasazi::Value
{
  std::string __str__()
  {
    std::stringstream output;
    output << self->realpart << "+" << self->imagpart << "j";
    return output.str();
  }
}
%template (ValueDouble)
  Anasazi::Value<double>;

///////////////////////////////////
// Anasazi OutputManager support //
///////////////////////////////////
%include "AnasaziOutputManager.hpp"
%template (OutputManagerDouble)
  Anasazi::OutputManager<double>;

////////////////////////////////////////
// Anasazi BasicOutputManager support //
////////////////////////////////////////
namespace std { struct ostream; }
%include "AnasaziBasicOutputManager.hpp"
%template (BasicOutputManagerDouble)
  Anasazi::BasicOutputManager< double >;

/////////////////////////////////
// Anasazi SortManager support //
/////////////////////////////////
%include "AnasaziSortManager.hpp"

///////////////////////////////
// Anasazi BasicSort support //
///////////////////////////////
%include "AnasaziBasicSort.hpp"

////////////////////////////////////
// Anasazi MultiVecTraits support //
////////////////////////////////////
%include "AnasaziMultiVecTraits.hpp"

//////////////////////////////
// Anasazi MultiVec support //
//////////////////////////////
%include "AnasaziMultiVec.hpp"
%template (MultiVecDouble)
  Anasazi::MultiVec<double>;

////////////////////////////////////
// Anasazi OperatorTraits support //
////////////////////////////////////
%include "AnasaziOperatorTraits.hpp"

//////////////////////////////
// Anasazi Operator support //
//////////////////////////////
%include "AnasaziOperator.hpp"
%teuchos_rcp(Anasazi::Operator< double >)
%template (OperatorDouble)
  Anasazi::Operator<double>;

//////////////////////////////////
// Anasazi Eigenproblem support //
//////////////////////////////////
%include "AnasaziEigenproblem.hpp"

///////////////////////////////////////
// Anasazi BasicEigenproblem support //
///////////////////////////////////////
%include "AnasaziBasicEigenproblem.hpp"

////////////////////////////////
// Anasazi StatusTest support //
////////////////////////////////
%include "AnasaziStatusTest.hpp"

/////////////////////////////////////
// Anasazi StatusTestCombo support //
/////////////////////////////////////
%include "AnasaziStatusTestCombo.hpp"

////////////////////////////////////////
// Anasazi StatusTestMaxIters support //
////////////////////////////////////////
%include "AnasaziStatusTestMaxIters.hpp"

//////////////////////////////////////
// Anasazi StatusTestOutput support //
//////////////////////////////////////
%include "AnasaziStatusTestOutput.hpp"

///////////////////////////////////////
// Anasazi StatusTestResNorm support //
///////////////////////////////////////
%include "AnasaziStatusTestResNorm.hpp"

//////////////////////////////////
// Anasazi OrthoManager support //
//////////////////////////////////
%ignore Anasazi::OrthoManager::project;
%ignore Anasazi::OrthoManager::projectAndNormalize;
%ignore Anasazi::OrthoManager::projectMat;
%ignore Anasazi::OrthoManager::projectAndNormalizeMat;
%include "AnasaziOrthoManager.hpp"

/////////////////////////////////////
// Anasazi MatOrthoManager support //
/////////////////////////////////////
%include "AnasaziMatOrthoManager.hpp"

///////////////////////////////////////
// Anasazi BasicOrthoManager support //
///////////////////////////////////////
%ignore Anasazi::BasicOrthoManager::normalizeMat;
%include "AnasaziBasicOrthoManager.hpp"

//////////////////////////////////////
// Anasazi SVQBOrthoManager support //
//////////////////////////////////////
%ignore Anasazi::SVQBOrthoManager::normalizeMat;
%include "AnasaziSVQBOrthoManager.hpp"

/////////////////////////////////
// Anasazi Eigensolver support //
/////////////////////////////////
%include "AnasaziEigensolver.hpp"

///////////////////////////////////
// Anasazi SolverManager support //
///////////////////////////////////
%include "AnasaziSolverManager.hpp"

///////////////////////////////////
// Anasazi BlockDavidson support //
///////////////////////////////////
%include "AnasaziBlockDavidson.hpp"

/////////////////////////////////////////
// Anasazi BlockDavidsonSolMgr support //
/////////////////////////////////////////
%include "AnasaziBlockDavidsonSolMgr.hpp"

//////////////////////////////////////
// Anasazi BlockKrylovSchur support //
//////////////////////////////////////
%include "AnasaziBlockKrylovSchur.hpp"

////////////////////////////////////////////
// Anasazi BlockKrylovSchurSolMgr support //
////////////////////////////////////////////
%include "AnasaziBlockKrylovSchurSolMgr.hpp"

////////////////////////////
// Anasazi LOBPCG support //
////////////////////////////
%include "AnasaziLOBPCG.hpp"

//////////////////////////////////
// Anasazi LOBPCGSolMgr support //
//////////////////////////////////
%include "AnasaziLOBPCGSolMgr.hpp"

///////////////////////////////////
// Anasazi EpetraAdapter support //
///////////////////////////////////
#ifdef HAVE_PYTRILINOS_EPETRA
%include "AnasaziEpetraAdapter.hpp"
%ignore
Anasazi::MultiVecTraits< double,
			 Epetra_MultiVector >::CloneView(const Epetra_MultiVector &,
							 const std::vector< int > &);
%template (SortManagerEpetra)
  Anasazi::SortManager< double >;
%template (BasicSortEpetra)
  Anasazi::BasicSort< double >;
%template (MultiVecTraitsEpetra)
  Anasazi::MultiVecTraits< double, Epetra_MultiVector >;
%template (OperatorTraitsEpetra)
  Anasazi::OperatorTraits< double, Epetra_MultiVector, Epetra_Operator >;
%template (EigenproblemEpetra)
  Anasazi::Eigenproblem< double, Epetra_MultiVector, Epetra_Operator >;
%template (BasicEigenproblemEpetra)
  Anasazi::BasicEigenproblem< double, Epetra_MultiVector, Epetra_Operator >;
%template (StatusTestEpetra)
  Anasazi::StatusTest< double, Epetra_MultiVector, Epetra_Operator >;
%template (StatusTestComboEpetra)
  Anasazi::StatusTestCombo< double, Epetra_MultiVector, Epetra_Operator >;
%template (StatusTestMaxItersEpetra)
  Anasazi::StatusTestMaxIters< double, Epetra_MultiVector, Epetra_Operator >;
%template (StatusTestOutputEpetra)
  Anasazi::StatusTestOutput< double, Epetra_MultiVector, Epetra_Operator >;
%template (StatusTestResNormEpetra)
  Anasazi::StatusTestResNorm< double, Epetra_MultiVector, Epetra_Operator >;
%template (OrthoManagerEpetra)
  Anasazi::OrthoManager< double, Epetra_MultiVector >;
%template (MatOrthoManagerEpetra)
  Anasazi::MatOrthoManager< double, Epetra_MultiVector, Epetra_Operator >;
%template (BasicOrthoManagerEpetra)
  Anasazi::BasicOrthoManager< double, Epetra_MultiVector, Epetra_Operator >;
%template (SVQBOrthoManagerEpetra)
  Anasazi::SVQBOrthoManager< double, Epetra_MultiVector, Epetra_Operator >;
%template (EigensolverEpetra)
  Anasazi::Eigensolver<  double, Epetra_MultiVector, Epetra_Operator >;
%template (SolverManagerEpetra)
  Anasazi::SolverManager<  double, Epetra_MultiVector, Epetra_Operator >;
%template (BlockDavidsonEpetra)
  Anasazi::BlockDavidson< double, Epetra_MultiVector, Epetra_Operator >;
%template (BlockDavidsonSolMgrEpetra)
  Anasazi::BlockDavidsonSolMgr< double, Epetra_MultiVector, Epetra_Operator >;
%template (BlockKrylovSchurEpetra)
  Anasazi::BlockKrylovSchur< double, Epetra_MultiVector, Epetra_Operator >;
%template (BlockKrylovSchurSolMgrEpetra)
  Anasazi::BlockKrylovSchurSolMgr< double, Epetra_MultiVector, Epetra_Operator >;
%template (LOBPCGEpetra)
  Anasazi::LOBPCG< double, Epetra_MultiVector, Epetra_Operator >;
%template (LOBPCGSolMgrEpetra)
  Anasazi::LOBPCGSolMgr< double, Epetra_MultiVector, Epetra_Operator >;
%template(EigensolutionEpetra)
  Anasazi::Eigensolution< double, Epetra_MultiVector >;
%pythoncode
%{
EigensolutionEpetra.Evals  = EigensolutionEpetra.evals
EigensolutionEpetra.Evecs  = EigensolutionEpetra.evecs
EigensolutionEpetra.Espace = EigensolutionEpetra.espace
%}
#endif

/////////////////////////
// std::vector support //
/////////////////////////
%template (VectorValueDouble) std::vector< Anasazi::Value< double > >;
%template (VectorInt        ) std::vector< int >;

//////////////////////////////
// Generic python interface //
//////////////////////////////
%define %anasazi_scalartype_factory(ClassName)
%pythoncode
%{
def ClassName(*args):
    """
    Factory function for class ClassName.  Currently, this returns a python
    wrapper around class ClassName<double>.
    """
    return ClassName##Double(*args)
%}
%enddef

%define %anasazi_factory(ClassName)
%pythoncode
%{
def ClassName(*args):
    """
    Factory function for class ClassName.  Currently, this returns a python
    wrapper around class ClassName<double,Epetra_MultiVector,Epetra_Operator>.
    """
    return ClassName##Epetra(*args)
%}
%enddef

// Implement the python factories for classes templated purely on ScalarType
%anasazi_scalartype_factory(Value             )
%anasazi_scalartype_factory(OutputManager     )
%anasazi_scalartype_factory(BasicOutputManager)
%anasazi_scalartype_factory(MultiVec          )
%anasazi_scalartype_factory(Operator          )

// Implement the python factories for classes templated on MV and/or OP 
%anasazi_factory(SortManager           )
%anasazi_factory(BasicSort             )
%anasazi_factory(MultiVecTraits        )
%anasazi_factory(OperatorTraits        )
%anasazi_factory(Eigenproblem          )
%anasazi_factory(BasicEigenproblem     )
%anasazi_factory(StatusTest            )
%anasazi_factory(StatusTestCombo       )
%anasazi_factory(StatusTestMaxIters    )
%anasazi_factory(StatusTestOutput      )
%anasazi_factory(StatusTestResNorm     )
%anasazi_factory(OrthoManager          )
%anasazi_factory(MatOrthoManager       )
%anasazi_factory(SVQBOrthoManager      )
%anasazi_factory(Eigensolver           )
%anasazi_factory(SolverManager         )
%anasazi_factory(BlockDavidson         )
%anasazi_factory(BlockDavidsonSolMgr   )
%anasazi_factory(BlockKrylovSchur      )
%anasazi_factory(BlockKrylovSchurSolMgr)
%anasazi_factory(LOBPCG                )
%anasazi_factory(LOBPCGSolMgr          )
%anasazi_factory(Eigensolution         )

// Turn off the exception handling
%exception;
