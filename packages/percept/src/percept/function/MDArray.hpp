// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef stk_encr_MDArray_hpp
#define stk_encr_MDArray_hpp

#include <string>
#include <stdexcept>
#include <vector>
#include <Kokkos_Core.hpp>
#include <Kokkos_DynRankView.hpp>
#include <percept/Util.hpp>

  namespace percept
  {

    typedef Kokkos::DynRankView<double,Kokkos::HostSpace> MDArray;
    typedef Kokkos::DynRankView<int,Kokkos::HostSpace> MDArrayInt;
    typedef Kokkos::DynRankView<unsigned,Kokkos::HostSpace> MDArrayUInt;

    template<class MDA>
    inline std::string
    printForMathematica(MDA& mda, bool convert=true, int precision=6)
    {
      std::ostringstream str;
      if (!convert)
        str << std::setprecision(precision);

      if (mda.rank() == 1)
        {
          str << "{";
          for (int i=0; i < mda.extent_int(0); ++i)
            {
              if (convert)
                str << Util::convert_to_mm(mda(i),precision) << (i == mda.extent_int(0) - 1 ?  "}" : ",");
              else
                str << mda(i) << (i == mda.extent_int(0) - 1 ?  "}" : ",");
            }
          return str.str();
        }
      else if (mda.rank() == 2)
        {
          str << "{";
          int nij[2] = {mda.extent_int(0), mda.extent_int(1)};
          int ind[2] = {0,0};
          for (ind[0]=0; ind[0] < nij[0]; ++ind[0])
            {
              str << "{";
              for (ind[1]=0; ind[1] < nij[1]; ++ind[1])
                {
                  if (convert)
                    str << Util::convert_to_mm(mda(ind[0],ind[1]),precision) << (ind[1] == nij[1] - 1 ? "}" : ",");
                  else
                    str << mda(ind[0],ind[1]) << (ind[1] == nij[1] - 1 ? "}" : ",");
                }
              str << (ind[0] == nij[0] - 1 ? "}" : ",");
            }
          return str.str();
        }
      else
        return "not ready";
    }

    template<class MDA>
    inline std::string printContainer(const MDA& mda) {

      // Save the format state of the original ostream os.
      std::ostringstream os;

      os.setf(std::ios_base::scientific, std::ios_base::floatfield);
      os.setf(std::ios_base::right);
      int myprec = os.precision();

      int size = mda.size();
      int rank = mda.rank();
      Teuchos::Array<int> multiIndex(rank);

      os<< "===============================================================================\n"\
        << "\t Container size = " << size << "\n"
        << "\t Container rank = " << rank << "\n" ;

      if( (rank == 0 ) && (size == 0) ) {
        os<< "====================================================================================\n"\
          << "|                        *** This is an empty container ****                       |\n";
      }
      else {
        os<< "\t extents     = ";

        for(int r = 0; r < rank; r++){
          os << " (" << mda.extent_int(r) <<") ";
        }
        os << "\n";

        os<< "============================================================\n"\
          << "|              Multi-index          Value                  |\n"\
          << "============================================================\n";
      }

      for(int i=0; i< mda.extent_int(0); ++i ) 
      {
        for(int j=0; j< mda.extent_int(1); ++j )
          for(int k=0; k< mda.extent_int(2); ++k )
            for(int l=0; l< mda.extent_int(3); ++l ) 
            {
              std::ostringstream mistring;
              mistring << i << std::dec << " ";
              if(rank > 1) mistring << j << std::dec << " ";
              if(rank > 2) mistring << k << std::dec << " ";
              if(rank > 3) mistring << l << std::dec << " ";
              os.setf(std::ios::right, std::ios::adjustfield);
              os << std::setw(27) << mistring.str();
              os.setf(std::ios::left, std::ios::adjustfield);
              os << std::setw(myprec+8) << mda.access(i,j,k,l) << "\n";
            }
      }

      os<< "====================================================================================\n\n";

      return os.str();
    }


    class MDArrayString
    {
      typedef std::vector<std::string > VecOfString;
      int m_rank;
      VecOfString m_array_1;  // rank one
      std::vector<VecOfString > m_array_2; // rank two
      //...

      void clone(const MDArrayString& mda)
      {
        this->m_rank = mda.m_rank;
        if (this->m_rank == 1)
          {
            this->m_array_1 = mda.m_array_1;
          }
        else
          {
            int dim0 = mda.extent_int(0);
            this->m_array_2.resize(dim0);
            for (int j = 0; j < dim0; j++)
              {
                this->m_array_2[j] = mda.m_array_2[j];
              }
          }
      }

    public:
      MDArrayString() : m_rank(1) { m_array_1.resize(0); }
      MDArrayString(int dim) : m_rank(1) { m_array_1.resize(dim); }
      MDArrayString(int dim0, int dim1) : m_rank(2)
      {
        m_array_2.resize(dim0);
        for (int j = 0; j < dim0; j++)
          {
            m_array_2[j].resize(dim1);
          }
      }

      MDArrayString(const MDArrayString& mda)
      {
        clone(mda);
      }
      MDArrayString& operator=(const MDArrayString& mda)
      {
        clone(mda);
        return *this;
      }

      void resize(int dim)
      {
        m_rank = 1;
        m_array_1.resize(dim);
      }
      void resize(int dim0, int dim1)
      {
        m_rank = 2;
        m_array_2.resize(dim0);
        for (int j = 0; j < dim0; j++)
          {
            m_array_2[j].resize(dim1);
          }

      }
      int rank() const { return m_rank; }

      void setValues(std::string *data)
      {
        if (m_rank == 1)
          {
            for (unsigned j = 0; j < m_array_1.size(); j++)
              m_array_1[j] = data[j];
          }
        else if (m_rank == 2)
          {
            unsigned k = 0;
            for (unsigned i = 0; i < m_array_2.size(); i++)
              for (unsigned j = 0; j < m_array_2[0].size(); j++)
                {
                  m_array_2[i][j] = data[k++];
                }
          }
      }

      //MDArrayString(std::string array[][]) : m_rank(2)
      std::string& operator()(int i1) {
        if (m_rank==2) throw std::runtime_error("MDArrayString:: rank 2 but asking for 1 dim");
        return m_array_1[i1];
      }
      std::string& operator()(int i1, int i2) {
        if (m_rank==1) throw std::runtime_error("MDArrayString:: rank 1 but asking for 2 dim");
        return m_array_2[i1][i2];
      }

      inline std::string& access(int i1) {
        return (*this)(i1);
      }

      inline std::string& access(int i1, int i2) {
        return (*this)(i1,i2);
      }

      inline std::string& access(int i1, int i2, int i3) {
        if (i3==0) 
          return access(i1,i2);
        else
          throw std::runtime_error("MDArrayString:: Rank at most 2 but requiring a rank 3 container");
      }

      inline std::string& access(int i1, int i2, int i3, int i4) {
        if (i4==0) 
          return access(i1,i2,i3);
        else
          throw std::runtime_error("MDArrayString:: Rank at most 2 but requiring a rank 4 container");
      }

      int extent_int(int i1) const {
        if (m_rank==1 && i1 > 0) throw std::runtime_error("MDArrayString:: rank 1 but asking for 2nd dim");
        if (m_rank == 1)
          {
            return m_array_1.size();
          }
        else if (m_rank ==2)
          {
            if (i1 == 0)
              {
                return m_array_2.size();
              }
            else if (i1 == 1)
              {
                return m_array_2[0].size();
              }
            else
              {
                throw std::runtime_error("MDArrayString:: asking for dimension greater than max rank of 2");
              }
          }
        return 0;
      }

    };

  }

#endif

