/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_Broker_hpp_
#define _snl_fei_Broker_hpp_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>

namespace fei {
  class VectorSpace;
  class MatrixGraph;
  class Vector;
  class Matrix;
  class LinearSystem;
}//namespace fei

namespace snl_fei {

  /** Internal interface. Similar to a factory, for creating
      Matrix/Vector/etc instances which are views of LinearSystemCore
      implementations. One of these will be instantiated and
      used internally by snl_fei::Factory.
  */
  class Broker {
  public:
    /** Usual virtual destructor. */
    virtual ~Broker(){}

    /** Produce an instance of an fei::Vector. This overloading of the
	create() method is for use by Broker implementations that are
	dispensing 'views' of vectors that reside in LinearSystemCore or
	FiniteElementData container implementations. In those cases, there is
	a distinction that must be made between solution-vectors and
	rhs-vectors.

	@param isSolutionVector
     */
    virtual fei::SharedPtr<fei::Vector> createVector(bool isSolutionVector=false) = 0;

    /** Produce an instance of an fei::Matrix
     */
    virtual fei::SharedPtr<fei::Matrix> createMatrix() = 0;

    /** Produce an instance of an fei::LinearSystem
     */
    virtual fei::SharedPtr<fei::LinearSystem> createLinearSystem() = 0;

    /** Set the MatrixGraph object used by this broker. */
    virtual void setMatrixGraph(fei::SharedPtr<fei::MatrixGraph> matrixGraph) = 0;
  };//class Broker
}//namespace snl_fei

#endif // _snl_fei_Broker_hpp_
