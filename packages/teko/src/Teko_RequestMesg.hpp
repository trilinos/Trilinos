// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_RequestMesg_hpp__
#define __Teko_RequestMesg_hpp__

#include <string>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Teko {

class RequestMesg {
 public:
  /** Construct a request messages specifing the
   * details of the request.
   *
   * \param[in] name Name of request to be satisfied.
   * \param[in] tag Optional tag describing other information
   *                about this tag.
   */
  explicit RequestMesg(const std::string& name, unsigned int tag = 0) : name_(name), tag_(tag) {}

  /** Construct a parameter list message. This sets the
   * name to "Parameter List" and the tag to 0. But now the
   * parameter list can completely describe the request.
   *
   * \param[in] pl Parameter list describing the request
   */
  explicit RequestMesg(const Teuchos::RCP<const Teuchos::ParameterList>& pl) : paramList_(pl) {
    fromParameterList(*paramList_);
  }

  //! Simple base class destructor
  virtual ~RequestMesg() {}

  //! Get the name for this request
  std::string getName() const { return name_; }

  //! Get the tag for this request
  unsigned int getTag() const { return tag_; }

  //! Get parameter list for this request
  const Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const {
    return paramList_.getConst();
  }

 protected:
  void fromParameterList(const Teuchos::ParameterList& pl) {
    name_ = "Parameter List";
    tag_  = 0;
    if (pl.isParameter("Name")) name_ = pl.get<std::string>("Name");
    if (pl.isParameter("Tag")) tag_ = pl.get<unsigned int>("Tag");
  }

  std::string name_;
  unsigned int tag_;
  Teuchos::RCP<const Teuchos::ParameterList> paramList_;
};

// simple stream interface for RequestMesg
inline std::ostream& operator<<(std::ostream& os, const Teko::RequestMesg& rm) {
  os << "RequestMesg <"
     << "name = \"" << rm.getName() << "\", "
     << "tag = " << rm.getTag() << ">";

  return os;
}

}  // end namespace Teko

#endif
