// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_MULTICONTINUATION_NATURALGROUP_H
#define LOCA_MULTICONTINUATION_NATURALGROUP_H

#include "LOCA_MultiContinuation_ExtendedGroup.H"       // base class

namespace LOCA { 

  namespace MultiContinuation { 

    /*! 
     * \brief Specialization of LOCA::MultiContinuation::ExtendedGroup
     * to natural continuation.
     */
    /*!
     * Natural continuation corresponds to a continuation equation
     * \f$g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=0\f$ with \f$g\f$ given by
     * \f[
     *     g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=
     *         p-p_0-v_p \Delta s
     * \f]
     * where \f$v_p\f$ is the parameter component of the predictor direction 
     * \f$v\f$. This corresponds geometrically to constraining the nonlinear 
     * solver steps used in calculating \f$F(x,p)=0\f$ to be orthogonal to the 
     * parameter axis.  The natural constraint \f$g\f$ is
     * represented by a LOCA::MultiContinuation::NaturalConstraint object.
     */
    class NaturalGroup : 
      public virtual LOCA::MultiContinuation::ExtendedGroup {

    public:
  
      //! Constructor
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list.
       * \param continuationParams [in] Continuation parameters.
       * \param grp [in] Group representing \f$F\f$.
       * \param pred [in] Predictor strategy.
       * \param paramIDs [in] Parameter IDs of continuation parameters.
       */
      NaturalGroup(
      const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
      const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& continuationParams,
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const vector<int>& paramIDs);

      //! Copy constructor
      NaturalGroup(const NaturalGroup& source, 
		     NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~NaturalGroup();

      /*! 
       * @name Implementation of NOX::Abstract::Group virtual methods 
       */
      //@{

      //! Assignment operator
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);

      //! Clone function
      virtual Teuchos::RefCountPtr<NOX::Abstract::Group>
      clone(NOX::CopyType type = NOX::DeepCopy) const;
      
      //@}

      /*! 
       * @name Implementation of LOCA::MultiContinuation::AbstractStrategy virtual methods 
       */
      //@{

      //! Copy
      virtual void copy(const NOX::Abstract::Group& source);

      //@}

    private:

      //! Prohibit generation and use of operator=()
      NaturalGroup& operator=(const NaturalGroup& source);

    };

  } // namespace MultiContinuation
} // namespace LOCA

#endif
