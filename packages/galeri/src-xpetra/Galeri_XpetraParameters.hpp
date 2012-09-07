// @HEADER
//
// ***********************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_GALLERYPARAMETERS_HPP
#define MUELU_GALLERYPARAMETERS_HPP

#include <Teuchos_Describable.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>


namespace Galeri {
  
  namespace Xpetra {
    
    // TODO nx/ny/nz == GO or global_size_t ? But what is the best to do?

    template<typename GO>
    class Parameters 
      : public Teuchos::VerboseObject<Parameters<GO> >, public Teuchos::Describable
    {

    public:
      
      Parameters(Teuchos::CommandLineProcessor& clp, GO nx=16, GO ny=-1, GO nz=-1, const std::string &
      matrixType="Laplace1D", int keepBCs=0): nx_(nx), ny_(ny), nz_(nz), matrixType_(matrixType), keepBCs_(keepBCs){
        clp.setOption("nx", &nx_, "mesh points in x-direction.");
        clp.setOption("ny", &ny_, "mesh points in y-direction.");
        clp.setOption("nz", &nz_, "mesh points in z-direction.");
        clp.setOption("matrixType", &matrixType_, "matrix type: Laplace1D, Laplace2D, Laplace3D"); //TODO: Star2D, numGlobalElements=...
        clp.setOption("keepBCs", &keepBCs_, "keep Dirichlet boundary rows in matrix (0=false,1=true)");
      }
      
      void check() const {
        //if (nx < 0) ...
      }
      
      GO GetNumGlobalElements() const {
        check();

        GO numGlobalElements=-1;
        if (matrixType_ == "Laplace1D")
          numGlobalElements = static_cast<GO>(nx_);
        else if (matrixType_ == "Laplace2D")
          numGlobalElements = static_cast<GO>(nx_*ny_);
        else if (matrixType_ == "Laplace3D")
          numGlobalElements = static_cast<GO>(nx_*ny_*nz_);
        //TODO else throw

        if (numGlobalElements < 0) throw Exceptions::Overflow("Gallery: numGlobalElements < 0 (did you forget --ny for 2D problems?)");

        return numGlobalElements;
      }

      const std::string & GetMatrixType() const {
        check();
        return matrixType_;
      }

      Teuchos::ParameterList & GetParameterList() const {
        paramList_ = Teuchos::ParameterList();

        check();

        paramList_.set("nx", static_cast<GO>(nx_));
        paramList_.set("ny", static_cast<GO>(ny_));
        paramList_.set("nz", static_cast<GO>(nz_));
        paramList_.set("keepBCs", static_cast<bool>(keepBCs_));

        return paramList_;
      }

      //! @name Overridden from Teuchos::Describable 
      //@{
    
      //! Return a simple one-line description of this object.
      std::string description() const {
        std::ostringstream out;
        out << Teuchos::Describable::description();
        out << "{type = "  << matrixType_ << ", size = " << GetNumGlobalElements() << "} ";
        return out.str();
      }
    
      //! Print the object with some verbosity level to an FancyOStream object.
      void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = verbLevel_default) const {
        using std::endl;
        int vl = (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;
        if (vl == VERB_NONE) return;
      
        if (vl == VERB_LOW) { out << description() << endl; } else { out << Teuchos::Describable::description() << endl; }
      
        if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
          Teuchos::OSTab tab1(out);

          out << "Matrix type: " << matrixType_ << endl
              << "Problem size: " << GetNumGlobalElements();

          if (matrixType_ == "Laplace2D")       out << " (" << nx_ << "x" << ny_ << ")";
          else if (matrixType_ == "Laplace3D")  out << " (" << nx_ << "x" << ny_ << "x" << nz_ << ")";

          out << endl;
        }
      
      }

      //@}

    private:
      // See Teuchos BUG 5249: https://software.sandia.gov/bugzilla/show_bug.cgi?id=5249
      double nx_;
      double ny_;
      double nz_;
      // GO nx_;
      // GO ny_;
      // GO nz_;

      std::string matrixType_;

      int keepBCs_;

      mutable Teuchos::ParameterList paramList_; // only used by GetParameterList(). It's temporary data. TODO: bad design...
    };
  
  }
}

#endif


//TODO: add capability to read from file ??

//GOAL: link between InputReader and ParameterList + Hide everything from examples
