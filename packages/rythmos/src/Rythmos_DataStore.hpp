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

#include "Thyra_VectorBase.hpp"
#include "Teuchos_RefCountPtr.hpp"
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
    DataStore() {};
    DataStore(Scalar &time_
      ,const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &x_
      ,const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &xdot_
      ,ScalarMag &accuracy_);

    /// Time value of data:
    Scalar time;

    /// Solution value of data at above time:
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x;

    /// Solution dot value of data at above time:
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdot;

    /// Accuracy of x data.  This is the accuracy of interpolations
    ScalarMag accuracy;

    /// Less than comparison for sorting:
    bool operator< (const DataStore<Scalar>& ds) const;

    /// Greather than comparison for sorting:
    bool operator> (const DataStore<Scalar>& ds) const;

    /// Equality comparison for matchin:
    bool operator== (const DataStore<Scalar>& ds) const;

    /// typedef for creating vectors of DataStore objects.
    typedef std::vector<DataStore<Scalar> > DataStoreVector_t;

    /// typedef for creating lists of DataStore objects.
    typedef std::list<DataStore<Scalar> > DataStoreList_t;

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
void DataStoreVectorToVector(
      const typename DataStore<Scalar>::DataStoreVector_t &ds
      ,std::vector<Scalar> *time_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *x_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *xdot_vec
      ,std::vector<typename DataStore<Scalar>::ScalarMag> *accuracy_vec);

// This is a helper function to convert vectors of t,x,xdot,accuracy to a vector of DataStore objects
template<class Scalar>
void VectorToDataStoreVector(
      const std::vector<Scalar> &time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &xdot_vec
      ,const std::vector<typename DataStore<Scalar>::ScalarMag> &accuracy_vec
      ,typename DataStore<Scalar>::DataStoreVector_t *ds);

// This is a helper function to convert vectors of t,x,xdot,[accuracy] to a list of DataStore objects
template<class Scalar>
void VectorToDataStoreList(
      const std::vector<Scalar> &time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &xdot_vec
      ,const std::vector<typename DataStore<Scalar>::ScalarMag> &accuracy_vec
      ,typename DataStore<Scalar>::DataStoreList_t *ds);

template<class Scalar>
void VectorToDataStoreList(
      const std::vector<Scalar> &time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &xdot_vec
      ,typename DataStore<Scalar>::DataStoreList_t *ds);

// DataStore definitions:
template<class Scalar>
DataStore<Scalar>::DataStore(
  Scalar &time_
  ,const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &x_
  ,const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &xdot_
  ,ScalarMag &accuracy_)
{
  time = time_;
  x = x_;
  xdot = xdot_;
  accuracy = accuracy_;
}

template<class Scalar>
bool DataStore<Scalar>::operator< (const DataStore<Scalar>& ds) const
{ 
  return( this->time < ds.time ); 
}

template<class Scalar>
bool DataStore<Scalar>::operator> (const DataStore<Scalar>& ds) const
{ 
  return( this->time > ds.time ); 
}

template<class Scalar>
bool DataStore<Scalar>::operator== (const DataStore<Scalar>& ds) const
{ 
  return( this->time == ds.time ); 
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
  if (verbLevel == Teuchos::VERB_EXTREME)
  {
    out << description() << "::describe:" << std::endl;
    out << "time = " << time << std::endl;
    out << "x = " << std::endl;
    x->describe(out,verbLevel);
    out << "xdot = " << std::endl;
    xdot->describe(out,verbLevel);
    out << "accuracy = " << accuracy << std::endl;
  }
}

// DataStore Helper Function definitions:
template<class Scalar>
void DataStoreVectorToVector(
      const typename DataStore<Scalar>::DataStoreVector_t &ds
      ,std::vector<Scalar> *time_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *x_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *xdot_vec
      ,std::vector<typename DataStore<Scalar>::ScalarMag> *accuracy_vec)
{
  time_vec->clear();
  x_vec->clear();
  xdot_vec->clear();
  accuracy_vec->clear();
  int N = ds.size();
  for (int i=0; i<N ; ++i)
  {
    time_vec->push_back(ds[i].time);
    x_vec->push_back(ds[i].x);
    xdot_vec->push_back(ds[i].xdot);
    accuracy_vec->push_back(ds[i].accuracy);
  }
}

template<class Scalar>
void VectorToDataStoreVector(
      const std::vector<Scalar> &time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &xdot_vec
      ,const std::vector<typename DataStore<Scalar>::ScalarMag> &accuracy_vec
      ,typename DataStore<Scalar>::DataStoreVector_t *ds
      ) 
{
  int N = time_vec.size();
  int Nx = x_vec.size();
  int Nxdot = xdot_vec.size();
  int Nacc = accuracy_vec.size();
  if ( (N != Nx) || (N != Nxdot) || (N != Nacc) )
  {
    ds = NULL;
    return;
  }
  ds->clear();
  for (int i=0; i<N ; ++i)
  {
    Scalar time_temp = time_vec[i];
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x_temp = x_vec[i];
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdot_temp = xdot_vec[i];
    typename DataStore<Scalar>::ScalarMag accuracy_temp = accuracy_vec[i];
    DataStore<Scalar> ds_tmp(time_temp,x_temp,xdot_temp,accuracy_temp);
    ds->push_back(ds_tmp);
  }
}

template<class Scalar>
void VectorToDataStoreList(
      const std::vector<Scalar> &time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &xdot_vec
      ,const std::vector<typename DataStore<Scalar>::ScalarMag> &accuracy_vec
      ,typename DataStore<Scalar>::DataStoreList_t *ds) 
{
  int N = time_vec.size();
  int Nx = x_vec.size();
  int Nxdot = xdot_vec.size();
  int Nacc = accuracy_vec.size();
  if ( (N != Nx) || (N != Nxdot) || (N != Nacc) )
  {
    ds = NULL;
    return;
  }
  ds->clear();
  for (int i=0; i<N ; ++i)
  {
    Scalar time_temp = time_vec[i];
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x_temp = x_vec[i];
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdot_temp = xdot_vec[i];
    typename DataStore<Scalar>::ScalarMag accuracy_temp = accuracy_vec[i];
    DataStore<Scalar> ds_tmp(time_temp,x_temp,xdot_temp,accuracy_temp);
    ds->push_back(ds_tmp);
  }
}

template<class Scalar>
void VectorToDataStoreList(
      const std::vector<Scalar> &time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &xdot_vec
      ,typename DataStore<Scalar>::DataStoreList_t *ds) 
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  std::vector<typename DataStore<Scalar>::ScalarMag> accuracy_vec;
  int N = time_vec.size();
  accuracy_vec.reserve(N);
  for (int i=0 ; i<N ; ++i)
    accuracy_vec.push_back(ST::zero());
  VectorToDataStoreList(time_vec,x_vec,xdot_vec,accuracy_vec,ds);
}


} // namespace Rythmos

#endif // Rythmos_DATA_STORE_H

