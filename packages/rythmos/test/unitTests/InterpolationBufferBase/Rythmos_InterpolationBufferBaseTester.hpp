//@HEADER

// ***********************************************************************
//
//                     Rythmos Package
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

#ifndef Rythmos_INTERPOLATION_BUFFER_BASE_TESTER_H
#define Rythmos_INTERPOLATION_BUFFER_BASE_TESTER_H

#include "Rythmos_InterpolationBufferBase.hpp"

namespace Rythmos {

template<class Scalar>
class InterpolationBufferBaseTester
  : virtual public Teuchos::Describable
  , virtual public Teuchos::ParameterListAcceptor
  , virtual public Teuchos::VerboseObject<InterpolationBufferBaseTester<Scalar> >
{
  public:

    ~InterpolationBufferBaseTester() {};

    InterpolationBufferBaseTester(); 

    bool check(const InterpolationBufferBase<Scalar> &IB
        ,const Thyra::VectorSpaceBase<Scalar>& vs
        ) const;

    bool checkSetPoints(const InterpolationBufferBase<Scalar> &IB
        ,const Thyra::VectorSpaceBase<Scalar>& vs
        ) const;
    bool checkGetPoints(const InterpolationBufferBase<Scalar> &IB
        ,const Thyra::VectorSpaceBase<Scalar>& vs
        ) const;
    bool checkSetRange(const InterpolationBufferBase<Scalar> &IB
        ,const Thyra::VectorSpaceBase<Scalar>& vs
        ) const;
    bool checkGetNodes(const InterpolationBufferBase<Scalar> &IB
        ,const Thyra::VectorSpaceBase<Scalar>& vs
        ) const;
    bool checkRemoveNodes(const InterpolationBufferBase<Scalar> &IB
        ,const Thyra::VectorSpaceBase<Scalar>& vs
        ) const;
    bool checkGetOrder(const InterpolationBufferBase<Scalar> &IB
        ,const Thyra::VectorSpaceBase<Scalar>& vs
        ) const;

    /// Redefined from describable
    /** \brief . */
    std::string description() const;

    /** \brief . */
    void describe(
      Teuchos::FancyOStream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const;

    /// Redefined from Teuchos::ParameterListAcceptor
    /** \brief . */
    void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();


  private:
    int num_rand_vectors;
    Teuchos::RefCountPtr<Teuchos::ParameterList> parameterList;
};

template<class Scalar>
InterpolationBufferBaseTester<Scalar>::InterpolationBufferBaseTester()
{
  num_rand_vectors=10;
}

template<class Scalar>
bool InterpolationBufferBaseTester<Scalar>::check(
    const InterpolationBufferBase<Scalar>& IB
    ,const Thyra::VectorSpaceBase<Scalar>& vs
    ) const
{
  bool status = false;
  status = checkSetPoints(IB,vs);
  if (!status) return(status);
  status = checkGetPoints(IB,vs);
  if (!status) return(status);
  status = checkSetRange(IB,vs);
  if (!status) return(status);
  status = checkGetNodes(IB,vs);
  if (!status) return(status);
  status = checkRemoveNodes(IB,vs);
  if (!status) return(status);
  status = checkGetOrder(IB,vs);
  return(status);
}

template<class Scalar>
bool InterpolationBufferBaseTester<Scalar>::checkSetPoints(
    const InterpolationBufferBase<Scalar> &IB
    ,const Thyra::VectorSpaceBase<Scalar>& vs
    ) const
{
  bool status=false;
  /*
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  typename Teuchos::ScalarTraits<Scalar> ST;
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  // insert valid points then GetNodes, check that t values are the same
  // insert valid points then GetNodes then GetPoints at node values, verify data comes back out
  Thyra::UniversalMultiVectorRandomizer<Scalar> randomizer;
  Thyra::seed_randomize<Scalar>(12345);
  Teuchos::ScalarTraits<Scalar>::seedrandom(12345);

  std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > x_vec;
  std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > xdot_vec;
  std::vector<Scalar> t_vec;
  std::vector<ScalarMag> accuracy_vec;
  int N = num_rand_vectors;
  for (int i=0 ; i<N ; ++i)
  {
    x_vec.push_back(Thyra::createMember(vs));
    xdot_vec.push_back(Thyra::createMember(vs));
  }
  for (int i=0 ; i<N ; ++i)
  {
    randomizer(&*x_vec[i]);
    randomizer(&*xdot_vec[i]);
    t_vec.push_back(Scalar(i/N)); // what about negative times?
    accuracy_vec.push_back(Teuchos::ScalarTraits<Scalar>::random()+ST::one());
  }
  IB.SetPoints(t_vec,x_vec,xdot_vec,accuracy_vec);
  std::vector<Scalar> t_vec_out;
  IB.GetNodes(&t_out_vec);
  for (int i=0 ; i<N ; ++i)
  {
    if (t_out_vec[i] != t_vec[i])
    {
      *out << "InterpolationBufferBaseTester IB = " << IB.description() << std::endl;
      *out << "Error:  IB.SetPoints() did not set a time value correctly" << std::endl;
    }
  }
  std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > x_vec_out;
  std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > xdot_vec_out;
  */




      // call with empty time_vec input
      // call with Teuchos::null input
      // call wtih x valid but xdot = Teuchos::null
      // call with x valid but random xdot = Teuchos::null
      // call with non-empty time_vec input and empty vectors 
      // call with non-empty time_vec input and vectors of different sizes
      // call with replacement data at a node and check that it was replaced correctly with GetNodes and GetPoints
      // call with both replacement data and new data and verify it works correctly
  return(status);
}

template<class Scalar>
bool InterpolationBufferBaseTester<Scalar>::checkGetPoints(
    const InterpolationBufferBase<Scalar> &IB
    ,const Thyra::VectorSpaceBase<Scalar>& vs
    ) const
{
      // call with empty time_vec input
      // call with values between nodes and verify that time value is between nodes
      //   verify that if x != null then x_out != null (same for xdot)
      // call with NULLs for output
      // call with partially filled vectors
  return(false);
}
template<class Scalar>
bool InterpolationBufferBaseTester<Scalar>::checkSetRange(
    const InterpolationBufferBase<Scalar> &IB
    ,const Thyra::VectorSpaceBase<Scalar>& vs
    ) const
{
      // call with another IB and then check with GetNodes and GetPoints
      // test time_lower and time_upper values for correctness of bound checking
      //   use valid and invalid ranges
      // test empty IB passed in
  return(false);
}
template<class Scalar>
bool InterpolationBufferBaseTester<Scalar>::checkGetNodes(
    const InterpolationBufferBase<Scalar> &IB
    ,const Thyra::VectorSpaceBase<Scalar>& vs
    ) const
{
      // call with null ptr
      // call with partially filled vector
  return(false);
}
template<class Scalar>
bool InterpolationBufferBaseTester<Scalar>::checkRemoveNodes(
    const InterpolationBufferBase<Scalar> &IB
    ,const Thyra::VectorSpaceBase<Scalar>& vs
    ) const
{
      // call and then then verify nodes are gone with GetNodes
      // call with empty vector
  return(false);
}
template<class Scalar>
bool InterpolationBufferBaseTester<Scalar>::checkGetOrder(
    const InterpolationBufferBase<Scalar> &IB
    ,const Thyra::VectorSpaceBase<Scalar>& vs
    ) const
{
  bool status = false;
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  // call and verify that a reasonable number comes out, i.e. larger than zero
  int order = IB.GetOrder();
  if (order >= 0) 
    status = true;
  else
  {
    *out << "InterpolationBufferBaseTester IB = " << IB.description() << std::endl;
    *out << "Error:  IB.GetOrder() returns a negative number" << std::endl;
  }
  return(status);
}

template<class Scalar>
std::string InterpolationBufferBaseTester<Scalar>::description() const
{
  std::string name = "Rythmos::InterpolationBufferBaseTester";
  return(name);
}

template<class Scalar>
void InterpolationBufferBaseTester<Scalar>::describe(
      Teuchos::FancyOStream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const
{
  if ( (static_cast<int>(verbLevel) == static_cast<int>(Teuchos::VERB_DEFAULT) ) ||
       (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)     )
     )
  {
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
  {
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM))
  {
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_HIGH))
  {
  }
}

template<class Scalar>
void InterpolationBufferBaseTester<Scalar>::setParameterList(
    Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList
    )
{
  parameterList = paramList;
}

template<class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> InterpolationBufferBaseTester<Scalar>::getParameterList()
{
  return(parameterList);
}

template<class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> InterpolationBufferBaseTester<Scalar>::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> temp_param_list = parameterList;
  parameterList = Teuchos::null;
  return(temp_param_list);
}


} // namespace Rythmos

#endif // Rythmos_INTERPOLATION_BUFFER_BASE_TESTER_H


