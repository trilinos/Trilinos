// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_LineSearch_Utils_Printing.H"
#include "NOX_Utils.H"
#include "NOX_Common.H"

NOX::LineSearch::Utils::Printing::
Printing(const Teuchos::RCP<NOX::Utils>& u) :
  NOX::Utils(*u)
{

}

NOX::LineSearch::Utils::Printing::~Printing()
{

}

void NOX::LineSearch::Utils::Printing::
reset(const Teuchos::RCP<NOX::Utils>& u)
{
  static_cast<NOX::Utils&>(*this) = *u;
}

void NOX::LineSearch::Utils::Printing::
printOpeningRemarks(const std::string& lineSearchName) const
{
  if (this->isPrintType(NOX::Utils::InnerIteration))
    {
      this->out() << "\n" << NOX::Utils::fill(72) << "\n"
          << "-- " << lineSearchName << " -- \n";
    }
}


void NOX::LineSearch::Utils::Printing::
printStep(int n, double step, double oldf, double newf, const std::string s,
      bool unscaleF) const
{
  if (isPrintType(NOX::Utils::InnerIteration))
  {
    this->out() << std::setw(3) << n << ":";
    this->out() << NOX::Utils::fill(1,' ') << "step = " << sciformat(step);
    if (unscaleF == true) {
      this->out() << NOX::Utils::fill(1,' ') << "old f = "
          << sciformat(std::sqrt(2. * oldf));
      this->out() << NOX::Utils::fill(1,' ') << "new f = "
          << sciformat(std::sqrt(2. * newf));
    }
    else {
      this->out() << NOX::Utils::fill(1,' ') << "old f = " << sciformat(oldf);
      this->out() << NOX::Utils::fill(1,' ') << "new f = " << sciformat(newf);
    }
    if (!s.empty())
    {
      this->out() << " " << s << "\n";
      this->out() << NOX::Utils::fill(72);
    }
    this->out() << std::endl;
  }
}
