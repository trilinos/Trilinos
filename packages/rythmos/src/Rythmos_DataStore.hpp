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

    /// Less than comparison for sorting:
    bool operator< (const DataStore<Scalar>& d1, const DataStore<Scalar>& d2) const
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

// This is a helper function to convert vectors of t,x,xdot to a list of DataStore objects
void VectorToDataStoreList(
  std::list<Teuchos::RefCountPtr<DataStore<Scalar> > > *list_ds
  ,const std::vector<ScalarMag> &time_vec
  ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &x_vec
  ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &xdot_vec) const;

// This is a helper function to convert a list of DataStore objects to vectors of t,x,xdot
void DataStoreListToVector(
  std::vector<ScalarMag> *time_vec
  ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *x_vec
  ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *xdot_vec
  ,const std::list<Teuchos::RefCountPtr<DataStore<Scalar> > > &list_ds) const;


// DataStore definitions:
template<class Scalar>
DataStore<Scalar>::DataStore(ScalarMag &time_
  ,Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &x_
  ,Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &xdot_)
{
  time = time_;
  x = x_;
  xdot = xdot_;
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
void VectorToDataStoreList(
      std::list<Teuchos::RefCountPtr<DataStore<Scalar> > > *list_ds
      ,const std::vector<ScalarMag> &time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > &xdot_vec) const
{
  int N = time_vec.size();
  int Nx = x_vec.size();
  int Nxdot = xdot_vec.size();
  if ((N != Nx) || (N != Nxdot))
  {
    list_ds = NULL;
    return;
  }
  list_ds->clear();
  for (int i=0; i<N ; ++i)
  {
    list_ds->push_back(Teuchos::rcp(new DataStore(time_list[i],x_list[i],xdot_list[i])));
  }
}

template<class Scalar>
void DataStoreListToVector(
      std::vector<ScalarMag> *time_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *x_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *xdot_vec
      ,const std::list<Teuchos::RefCountPtr<DataStore<Scalar> > > &list_ds) const
{
  int N = list_ds.size();
  time_list->reserve(N); time_list->clear();
  x_list->reserve(N);    x_list->clear();
  xdot_list->reserve(N); xdot_list->clear();
  std::list<DataStore<Scalar> >::iterator list_it = list_ds.begin();
  for (; list_it != list_ds.end() ; list_it++)
  {
    time_list->push_back(list_it->time);
    x_list->push_back(list_it->x);
    x_dot_list->push_back(list_it->xdot);
  }
}


} // namespace Rythmos

#endif // Rythmos_DATA_STORE_H

