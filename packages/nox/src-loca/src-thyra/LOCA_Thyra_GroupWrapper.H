// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Glen Hansen (gahanse@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_THYRA_GROUPWRAPPER_H
#define LOCA_THYRA_GROUPWRAPPER_H

#include "LOCA_Thyra_Group.H"                   // base class

// Forward declares
namespace Teuchos {
  class ParameterList;
}

namespace LOCA {

  //! %LOCA %Thyra support
  namespace Thyra {

    //! Extension of the LOCA::Thyra::Group to provide access to solution data
    class GroupWrapper : public LOCA::Thyra::Group {

    public:

      //! Constructor
      GroupWrapper(const Teuchos::RCP<LOCA::GlobalData>& global_data,
        const NOX::Thyra::Vector& initial_guess,
        const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model,
        const LOCA::ParameterVector& p,
        int p_index,
        bool implement_dfdp = false);

      /*!
       * \brief Copy constructor. If type is DeepCopy, takes ownership of
       * valid shared Jacobian.
       */
      GroupWrapper(const GroupWrapper& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~GroupWrapper();

      //! Assignment operator
      virtual GroupWrapper& operator=(const GroupWrapper& source);

      /*!
       * @name Overloaded NOX::Thyra::Group  methods.
       */
      //@{

      //! Assignment operator.
      virtual NOX::Abstract::Group&
      operator=(const NOX::Abstract::Group& source);

      //! Assignment operator.
      virtual NOX::Abstract::Group&
      operator=(const NOX::Thyra::Group& source);

      //! Cloning function
      virtual Teuchos::RCP<NOX::Abstract::Group>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Access the solution
      const Teuchos::RCP<const NOX::Thyra::Vector> getNOXThyraVecRCPX(){ return x_vec_; }

    };

  } // namespace Thyra
} // namespace LOCA


#endif
