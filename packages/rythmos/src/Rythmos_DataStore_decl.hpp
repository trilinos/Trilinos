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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_DATA_STORE_DECL_H
#define Rythmos_DATA_STORE_DECL_H

#include "Rythmos_Types.hpp"
#include "Thyra_VectorBase.hpp"
#include "Teuchos_Describable.hpp"

namespace Rythmos {

template<class Scalar>
class DataStore : virtual public Teuchos::Describable
{

  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /// Destructor
    ~DataStore() {};

    /** \brief. */
    DataStore();

    /** \brief. */
    // This is a shallow copy constructor, use clone for a deep copy
    DataStore(Scalar& time_
      ,const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& x_
      ,const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdot_
      ,ScalarMag& accuracy_);

    /** \brief. */
    // This is a shallow copy constructor, use clone for a deep copy
    DataStore(const DataStore<Scalar>& ds_in);
    
    /** \brief. */
    // This is a deep clone and copies the underlying vectors
    RCP<DataStore<Scalar> > clone() const;

    /// Time value of data:
    Scalar time;

    /// Solution value of data at above time:
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x;

    /// Solution dot value of data at above time:
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot;

    /// Accuracy of x data.  This is the accuracy of interpolations
    ScalarMag accuracy;

    /// Less than comparison for sorting:
    bool operator< (const DataStore<Scalar>& ds) const;

    /// Less than comparison for sorting:
    bool operator<= (const DataStore<Scalar>& ds) const;

    /// Less than comparison for sorting:
    bool operator< (const Scalar& t) const;

    /// Less than comparison for sorting:
    bool operator<= (const Scalar& t) const;

    /// Less than comparison for sorting:
    bool operator> (const DataStore<Scalar>& ds) const;

    /// Less than comparison for sorting:
    bool operator>= (const DataStore<Scalar>& ds) const;

    /// Less than comparison for sorting:
    bool operator> (const Scalar& t) const;

    /// Less than comparison for sorting:
    bool operator>= (const Scalar& t) const;

    /// Equality comparison for matching:
    bool operator== (const DataStore<Scalar>& ds) const;

    /// Equality comparison for matching:
    bool operator== (const Scalar& t) const;

    /// typedef for creating vectors of DataStore objects.
    typedef Array<DataStore<Scalar> > DataStoreVector_t;

    /// typedef for creating vectors of DataStore objects.
    typedef Array<const DataStore<Scalar> > constDataStoreVector_t;

    /// typedef for creating lists of DataStore objects.
    typedef std::list<DataStore<Scalar> > DataStoreList_t;

    /// typedef for creating lists of DataStore objects.
    typedef std::list<const DataStore<Scalar> > constDataStoreList_t;

    /// Inherited from Describable:
    /** \brief . */
    std::string description() const;

    /** \brief . */
    /** \brief . */
    void describe(
      Teuchos::FancyOStream       &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const;
};


// This is a helper function to convert a vector of DataStore objects to vectors of t,x,xdot,accuracy
template<class Scalar>
void dataStoreVectorToVector(
      const typename DataStore<Scalar>::DataStoreVector_t &ds
      ,Array<Scalar> *time_vec
      ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > *x_vec
      ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > *xdot_vec
      ,Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> *accuracy_vec);

// This is a helper function to convert vectors of t,x,xdot,accuracy to a vector of DataStore objects
template<class Scalar>
void vectorToDataStoreVector(
      const Array<Scalar> &time_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &x_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &xdot_vec
      ,const Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &accuracy_vec
      ,typename DataStore<Scalar>::DataStoreVector_t *ds);

// This is a helper function to convert vectors of t,x,xdot,[accuracy] to a list of DataStore objects
template<class Scalar>
void vectorToDataStoreList(
      const Array<Scalar> &time_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &x_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &xdot_vec
      ,const Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &accuracy_vec
      ,typename DataStore<Scalar>::DataStoreList_t *ds);

template<class Scalar>
void vectorToDataStoreList(
      const Array<Scalar> &time_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &x_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &xdot_vec
      ,typename DataStore<Scalar>::DataStoreList_t *ds);

} // namespace Rythmos

#endif // Rythmos_DATA_STORE_DECL_H

