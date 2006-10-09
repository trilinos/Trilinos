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
#include "Teuchos_Describable.hpp"

namespace Rythmos {

template<class Scalar>
class DataStore : virtual public Teuchos::Describable
{
  public:

    /// Destructor
    ~DataStore();

    /** \brief. */
    DataStore();
    DataStore(ScalarMag &time
      ,Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &x
      ,Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &xdot);

    /// Time value of data:
    ScalarMag time;

    /// Solution value of data at above time:
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x;

    /// Solution dot value of data at above time:
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdot;

    /// Accuracy of x data.  This is the accuracy of interpolations
    ScalarMag accuracy;

    /// Less than comparison for sorting:
    bool operator< (
        const DataStore<Scalar>& d1
        ,const DataStore<Scalar>& d2) const
    { return( d1.time < d2.time ); }

    /// Redefined from describable
    /** \brief . */
    std::string description() const;

    /** \brief . */
    void describe(
      Teuchos::FancyOStream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const;

}

// This is a helper function to convert a sequence of DataStore objects to vectors of t,x,xdot,accuracy
void DataStoreSeqToVector(
      const std::sequence<Teuchos::RefCountPtr<DataStore<Scalar> > > &ds
      ,std::vector<ScalarMag> *time_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *x_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *xdot_vec
      ,std::vector<ScalarMag> *accuracy_vec) const

// This is a helper function to convert vectors of t,x,xdot,accuracy to a sequence of DataStore objects
void VectorToDataStoreSeq(
      const std::vector<ScalarMag> &time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &xdot_vec
      ,const std::vector<ScalarMag> &accuracy_vec
      ,std::sequence<Teuchos::RefCountPtr<DataStore<Scalar> > > *ds) const

// DataStore definitions:
template<class Scalar>
DataStore<Scalar>::DataStore(ScalarMag &time_
  ,Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &x_
  ,Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &xdot_
  ,ScalarMag &accuracy_)
{
  time = time_;
  x = x_;
  xdot = xdot_;
  accuracy = accuracy_;
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
    out << description() << "::describe:";
}

// DataStore Helper Function definitions:
template<class Scalar>
void DataStoreSeqToVector(
      const std::sequence<Teuchos::RefCountPtr<DataStore<Scalar> > > &ds
      ,std::vector<ScalarMag> *time_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *x_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *xdot_vec
      ,std::vector<ScalarMag> *accuracy_vec) const
{
  std::sequence<DataStore<Scalar> >::iterator ds_it = ds.front();
  for (; ds_it != ds.end() ; ds_it++)
  {
    time_vec->push_back(ds_it->time);
    x_vec->push_back(ds_it->x);
    xdot_vec->push_back(ds_it->xdot);
    accuracy_vec->push_back(ds_it->accuracy);
  }
}

template<class Scalar>
void VectorToDataStoreSeq(
      const std::vector<ScalarMag> &time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &xdot_vec
      ,const std::vector<ScalarMag> &accuracy_vec
      ,std::sequence<Teuchos::RefCountPtr<DataStore<Scalar> > > *ds) const
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
    ds->insert(ds.end(),Teuchos::rcp(new DataStore(time_vec[i],x_vec[i],xdot_vec[i],accuracy_vec[i])));
  }
}

} // namespace Rythmos

#endif // Rythmos_DATA_STORE_H

