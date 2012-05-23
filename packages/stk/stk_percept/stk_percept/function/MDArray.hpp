#ifndef stk_encr_MDArray_hpp
#define stk_encr_MDArray_hpp

#include <string>
#include <stdexcept>
#include <vector>
#include "Intrepid_FieldContainer.hpp"

namespace stk
{
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

    //typedef Intrepid::FieldContainer<std::string> MDArrayString;
    class MDArrayString
    {
      typedef std::vector<std::string > VecOfString;
      int m_rank;
      VecOfString m_array_1;  // rank one
      std::vector<VecOfString > m_array_2; // rank two
      //...
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

      int dimension(int i1) { 
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
}

#endif

