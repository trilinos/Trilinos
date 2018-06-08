// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_ARRAYFIREVECTOR_H
#define ROL_ARRAYFIREVECTOR_H

#include "ROL_Vector.hpp"
#include "arrayfire.h"

/** \class ROL::ArrayFireVector
    \brief Provides the ArrayFire vector implementation of the ROL::Vector interface.
*/


namespace ROL {


template <class Real, class Element = Real>
class ArrayFireVector : public ROL::Vector<Real> {

private:
    
    const ROL::Ptr<af::array> afVector;

public:
    
    ArrayFireVector(const ROL::Ptr<af::array> &initVector) : afVector(initVector) {}

    void set (const ROL::Vector<Real> &newValue)  {
        const ArrayFireVector<Real, Element> &newVector = dynamic_cast<const ArrayFireVector<Real, Element>&>(newValue);
        const af::array &newArray = *(newVector.getVector());
        *afVector = newArray;
    }

    void plus (const ROL::Vector<Real> &operand)  {
        const ArrayFireVector<Real, Element> &operandVector = dynamic_cast<const ArrayFireVector<Real, Element>&>(operand);
        const af::array &operandArray = *(operandVector.getVector());
        
        dim_t currentDim = afVector->dims(0);
        dim_t operandDim = operandArray.dims(0);
        
        if (currentDim == operandDim) {
            *afVector += operandArray;
        }
    }

    void axpy (const Real alpha, const ROL::Vector<Real> &operand)   {
        const ArrayFireVector<Real, Element> &operandVector = dynamic_cast<const ArrayFireVector<Real, Element>&>(operand);
        const af::array &operandArray = *(operandVector.getVector());

        dim_t currentDim = afVector->dims(0);
        dim_t operandDim = operandArray.dims(0);

        if (currentDim == operandDim) {
            *afVector += operandArray * static_cast<Element>(alpha);
        }
    }

    void scale (const Real alpha)   {
        *afVector *= static_cast<Element>(alpha);
    }

    Real dot (const ROL::Vector<Real> &operand) const   {
        const ArrayFireVector<Real, Element> &operandVector = dynamic_cast<const ArrayFireVector<Real, Element>&>(operand);
        const af::array &operandArray = *(operandVector.getVector());
        
        dim_t currentDim = afVector->dims(0);
        dim_t operandDim = operandArray.dims(0);

        Real value = 0;

        if (currentDim == operandDim) {
            af::array dotProduct = af::dot(*afVector, operandArray);
            value = static_cast<Real>(dotProduct.scalar<Element>());
        }

        return value;
    }

    Real norm () const   {
        return static_cast<Real>(af::norm(*afVector));
    }

    ROL::Ptr<ROL::Vector<Real> > clone () const   {
        return(ROL::makePtr<ArrayFireVector<Real, Element>(ROL::makePtr<af::array(afVector->dims(0), afVector->type>>())));
    }

    ROL::Ptr<const af::array> getVector () const {
        return afVector;
    }

    ROL::Ptr< af::array> getVector()  {
        return afVector;
    }

/*  to be implemented later
    ROL::Ptr<ROL::Vector<Real> > basis (const int index) const   {
		ROL::Ptr<af::array> basisArrary = af::constant(0, afVector.dims(0), afVector.type());
        basisArrary(index) = 1;
        
        ArrayFireVector<Real, Element> *basisVector = new ArrayFireVector<Real, Element>(basisArrary);
        
        return basisVector;
    }
   
    
    int dimension() const   {
        dim_t currentDim = afVector.dims(0);
        return static_cast<int>(currentDim);
    }
	
    void applyUnary (const ROL::Elementwise::UnaryFunction<Real> &function)   {
        Element * vectorData = afVector.host<Element>();
        
        dim_t currentDim = afVector.dims(0);
        
        for (uint32_t i = 0; i < currentDim; i++) {
            vectorData[i] = function.apply(vectorData[i]);
        }
        
        af::array arrayFromHost = af::array(currentDim, vectorData);
        afVector = arrayFromHost;
    }
    
    
    void applyBinary (const ROL::Elementwise::BinaryFunction<Real> &function,
                      const ROL::Vector<Real> &operand)   {
        const ArrayFireVector<Real, Element> &operandVector = dynamic_cast<const ArrayFireVector&&>(operand);
        const af::array &operandArray = operandVector.getVector();
        
        Element * vectorData = afVector.host<Element>();
        Element * operandData = operandArray.host<Element>();
        dim_t currentDim = afVector.dims(0);
        
        for (uint32_t i = 0; i < currentDim; i++) {
            vectorData[i] = function.apply(vectorData[i], operandData[i]);
        }
        
        af::array arrayFromHost = af::array(currentDim, vectorData);
        afVector = arrayFromHost;
    }
    
    
    Real reduce (const ROL::Elementwise::ReductionOp<Real> &reducer) const   {
        Real result = reducer.initialValue();
        
        Element * hostData = afVector.host<Element>();
        dim_t currentDim = afVector.dims(0);
        
        for (uint32_t i = 0; i < currentDim; i++) {
            reducer.reduce(hostData[i], result);
        }
        
        return result;
    }
*/

};

} // namespace ROL

#endif
