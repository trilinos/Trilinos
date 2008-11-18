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

#ifndef Rythmos_DATA_STORE_H
#define Rythmos_DATA_STORE_H

#include "Rythmos_Types.hpp"
#include "Thyra_VectorBase.hpp"
#include "Teuchos_RCP.hpp"
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

// DataStore definitions:
template<class Scalar>
DataStore<Scalar>::DataStore()
  :time(-1),
   accuracy(-1)
{}

template<class Scalar>
DataStore<Scalar>::DataStore(
  Scalar &time_
  ,const Teuchos::RCP<const Thyra::VectorBase<Scalar> > &x_
  ,const Teuchos::RCP<const Thyra::VectorBase<Scalar> > &xdot_
  ,ScalarMag &accuracy_)
{
  time = time_;
  x = x_;
  xdot = xdot_;
  accuracy = accuracy_;
}

template<class Scalar>
DataStore<Scalar>::DataStore(
    const DataStore<Scalar>& ds_in
    )
{
  time = ds_in.time;
  x = ds_in.x;
  xdot = ds_in.xdot;
  accuracy = ds_in.accuracy;
}

template<class Scalar>
RCP<DataStore<Scalar> > DataStore<Scalar>::clone() const
{
  Scalar t_out = time;
  RCP<VectorBase<Scalar> > x_out;
  if (!Teuchos::is_null(x)) {
    x_out = x->clone_v();
  }
  RCP<VectorBase<Scalar> > xdot_out;
  if (!Teuchos::is_null(xdot)) {
    xdot_out = xdot->clone_v();
  }
  ScalarMag accuracy_out = accuracy;
  RCP<DataStore<Scalar> > ds_out = rcp(new DataStore<Scalar>(t_out,x_out,xdot_out,accuracy_out));
  return ds_out;
}

template<class Scalar>
bool DataStore<Scalar>::operator< (const DataStore<Scalar>& ds) const
{ 
  return( this->time < ds.time ); 
}

template<class Scalar>
bool DataStore<Scalar>::operator<= (const DataStore<Scalar>& ds) const
{ 
  return( this->time <= ds.time ); 
}

template<class Scalar>
bool DataStore<Scalar>::operator< (const Scalar& t) const
{ 
  return( this->time < t ); 
}

template<class Scalar>
bool DataStore<Scalar>::operator<= (const Scalar& t) const
{ 
  return( this->time <= t ); 
}

template<class Scalar>
bool DataStore<Scalar>::operator> (const DataStore<Scalar>& ds) const
{ 
  return( this->time > ds.time ); 
}

template<class Scalar>
bool DataStore<Scalar>::operator>= (const DataStore<Scalar>& ds) const
{ 
  return( this->time >= ds.time ); 
}

template<class Scalar>
bool DataStore<Scalar>::operator> (const Scalar& t) const
{ 
  return( this->time > t ); 
}

template<class Scalar>
bool DataStore<Scalar>::operator>= (const Scalar& t) const
{ 
  return( this->time >= t ); 
}

template<class Scalar>
bool DataStore<Scalar>::operator== (const DataStore<Scalar>& ds) const
{ 
  return( this->time == ds.time ); 
}

template<class Scalar>
bool DataStore<Scalar>::operator== (const Scalar& t) const
{ 
  return( this->time == t ); 
}

template<class Scalar>
std::string DataStore<Scalar>::description() const
{
  std::string name = "Rythmos::DataStore";
  return(name);
}

template<class Scalar>
void DataStore<Scalar>::describe(
      Teuchos::FancyOStream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const
{
  if (verbLevel == Teuchos::VERB_EXTREME) {
    out << description() << "::describe:" << std::endl;
    out << "time = " << time << std::endl;
    out << "x = " << std::endl;
    x->describe(out,verbLevel);
    if (xdot != Teuchos::null) {
      out << "xdot = " << std::endl;
      xdot->describe(out,verbLevel);
    }
    out << "accuracy = " << accuracy << std::endl;
  }
}

// DataStore Helper Function definitions:
template<class Scalar>
void dataStoreVectorToVector(
      const typename DataStore<Scalar>::DataStoreVector_t &ds
      ,Array<Scalar> *time_vec
      ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > *x_vec
      ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > *xdot_vec
      ,Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> *accuracy_vec)
{
  if(time_vec)
    time_vec->clear();
  if(x_vec)
    x_vec->clear();
  if(xdot_vec)
    xdot_vec->clear();
  if(accuracy_vec)
    accuracy_vec->clear();
  int N = ds.size();
  for (int i=0; i<N ; ++i) {
    if(time_vec)
      time_vec->push_back(ds[i].time);
    if(x_vec)
      x_vec->push_back(ds[i].x);
    if(xdot_vec)
      xdot_vec->push_back(ds[i].xdot);
    if(accuracy_vec)
      accuracy_vec->push_back(ds[i].accuracy);
  }
}

template<class Scalar>
void vectorToDataStoreVector(
      const Array<Scalar> &time_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &x_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &xdot_vec
      ,const Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &accuracy_vec
      ,typename DataStore<Scalar>::DataStoreVector_t *ds
      ) 
{
  int N = time_vec.size();
  int Nx = x_vec.size();
  int Nxdot = xdot_vec.size();
  int Nacc = accuracy_vec.size();
  if ( (N != Nx) || (N != Nxdot) || (N != Nacc) ) {
    ds = NULL;
    return;
  }
  ds->clear();
  for (int i=0; i<N ; ++i) {
    Scalar time_temp = time_vec[i];
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x_temp = x_vec[i];
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot_temp = xdot_vec[i];
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType accuracy_temp = accuracy_vec[i];
    DataStore<Scalar> ds_tmp(time_temp,x_temp,xdot_temp,accuracy_temp);
    ds->push_back(ds_tmp);
  }
}

template<class Scalar>
void vectorToDataStoreList(
      const Array<Scalar> &time_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &x_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &xdot_vec
      ,const Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &accuracy_vec
      ,typename DataStore<Scalar>::DataStoreList_t *ds) 
{
  int N = time_vec.size();
  int Nx = x_vec.size();
  int Nxdot = xdot_vec.size();
  int Nacc = accuracy_vec.size();
  if ( (N != Nx) || (N != Nxdot) || (N != Nacc) ) {
    ds = NULL;
    return;
  }
  ds->clear();
  for (int i=0; i<N ; ++i) {
    Scalar time_temp = time_vec[i];
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x_temp = x_vec[i];
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot_temp = xdot_vec[i];
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType accuracy_temp = accuracy_vec[i];
    DataStore<Scalar> ds_tmp(time_temp,x_temp,xdot_temp,accuracy_temp);
    ds->push_back(ds_tmp);
  }
}

template<class Scalar>
void vectorToDataStoreList(
      const Array<Scalar> &time_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &x_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &xdot_vec
      ,typename DataStore<Scalar>::DataStoreList_t *ds) 
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> accuracy_vec;
  int N = time_vec.size();
  accuracy_vec.reserve(N);
  for (int i=0 ; i<N ; ++i) {
    accuracy_vec.push_back(ST::zero());
  }
  vectorToDataStoreList(time_vec,x_vec,xdot_vec,accuracy_vec,ds);
}

} // namespace Rythmos

#endif // Rythmos_DATA_STORE_H

