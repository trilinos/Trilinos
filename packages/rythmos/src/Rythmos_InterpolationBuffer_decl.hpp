//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_INTERPOLATION_BUFFER_DECL_H
#define Rythmos_INTERPOLATION_BUFFER_DECL_H

#include "Rythmos_InterpolationBufferBase.hpp"
#include "Rythmos_Types.hpp"
#include "Rythmos_DataStore.hpp"
#include "Rythmos_InterpolatorAcceptingObjectBase.hpp"



namespace Rythmos {

enum IBPolicy {
  BUFFER_POLICY_INVALID = 0,
  BUFFER_POLICY_STATIC = 1,
  BUFFER_POLICY_KEEP_NEWEST = 2
};


/** \brief concrete class for interpolation buffer functionality. */
template<class Scalar> 
class InterpolationBuffer : 
  virtual public InterpolationBufferBase<Scalar>,
  virtual public InterpolatorAcceptingObjectBase<Scalar>
{
public:

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /// Redefined from Rythmos::InterpolationBufferBase
  /** \brief. */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
    
  /** \brief. */
  InterpolationBuffer();

  /// Initialize the buffer:
  void initialize( const RCP<InterpolatorBase<Scalar> >& interpolator, int storage );

  /** \brief Redefined from Rythmos::InterpolatorAcceptingObjectBase */
  //@{
  
  /// Set the interpolator for this buffer
  void setInterpolator(const RCP<InterpolatorBase<Scalar> >& interpolator);
  
  /** \brief . */
  RCP<InterpolatorBase<Scalar> >
    getNonconstInterpolator();

  /** \brief . */
  RCP<const InterpolatorBase<Scalar> >
    getInterpolator() const;
  
  /// Unset the interpolator for this buffer
  RCP<InterpolatorBase<Scalar> > unSetInterpolator();

  //@}

  /// Set the maximum storage of this buffer
  void setStorage( int storage );
    
  /// Get the maximum storage of this buffer
  int getStorage() const;

  /** \brief . */
  IBPolicy getIBPolicy();
        
  /// Destructor
  ~InterpolationBuffer() {};

  /// Add point to buffer
  void addPoints(
    const Array<Scalar>& time_vec
    ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec
    ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec);

  /// Get value from buffer
  void getPoints(
    const Array<Scalar>& time_vec
    ,Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec
    ,Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec
    ,Array<ScalarMag>* accuracy_vec
    ) const;

  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;

  /// Get interpolation nodes
  void getNodes(Array<Scalar>* time_vec) const;

  /// Get order of interpolation
  int getOrder() const;
    
  /// Remove interpolation nodes
  void removeNodes(Array<Scalar>& time_vec);

  /// Redefined from Teuchos::Describable
  /** \brief . */
  std::string description() const;

  /** \brief . */
  void describe(
    Teuchos::FancyOStream       &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ) const;

  /// Redefined from Teuchos::ParameterListAcceptor
  /** \brief . */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);

  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();

  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();

  RCP<const Teuchos::ParameterList> getValidParameters() const;
    
private:

  RCP<InterpolatorBase<Scalar> > interpolator_;
  int storage_limit_;
  RCP<typename DataStore<Scalar>::DataStoreVector_t> data_vec_;

  RCP<Teuchos::ParameterList> paramList_;

  IBPolicy policy_;


  // Private member functions:
  void defaultInitializeAll_();

};

// Nonmember constructor
template<class Scalar>
RCP<InterpolationBuffer<Scalar> > interpolationBuffer( 
  const RCP<InterpolatorBase<Scalar> >& interpolator = Teuchos::null,
  int storage = 0 
  );


} // namespace Rythmos


#endif // Rythmos_INTERPOLATION_BUFFER_DECL_H
