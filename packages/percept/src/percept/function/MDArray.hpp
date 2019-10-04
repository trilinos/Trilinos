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
#include "Intrepid_FieldContainer.hpp"
#include <percept/Util.hpp>

  namespace percept
  {

    //     class MDArray : public FieldContainer<double>
    //     {
    //       public:
    //       typedef FieldContainer<double> base;
    //       MDArray(std::vector<int> dimensions) : FieldContainer<double>( Teuchos::Array<int>(dimensions.begin(), dimensions.end()) ) {}
    //     };

    typedef Intrepid::FieldContainer<double> MDArray;
    typedef Intrepid::FieldContainer<int> MDArrayInt;
    typedef Intrepid::FieldContainer<unsigned> MDArrayUInt;

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
          for (int i=0; i < mda.dimension(0); ++i)
            {
              if (convert)
                str << Util::convert_to_mm(mda(i),precision) << (i == mda.dimension(0) - 1 ?  "}" : ",");
              else
                str << mda(i) << (i == mda.dimension(0) - 1 ?  "}" : ",");
            }
          return str.str();
        }
      else if (mda.rank() == 2)
        {
          str << "{";
          int nij[2] = {mda.dimension(0), mda.dimension(1)};
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


    //typedef Intrepid::FieldContainer<std::string> MDArrayString;
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
            int dim0 = mda.dimension(0);
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
      const int rank() const { return m_rank; }

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

      int dimension(int i1) const {
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

