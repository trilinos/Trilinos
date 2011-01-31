#ifndef stk_percept_hashtable_hpp
#define  stk_percept_hashtable_hpp

#include <stk_percept/Util.hpp>

#include <stk_percept/Teuchos_My_Hashtable.hpp>


namespace stk
{
  namespace percept 
  {
    template<class Key, class Value> 
    class Hashtable;

    template<class Key, class Value> 
    class HT_Iterator
    {
    public:
      typedef typename Teuchos::HashPair<Key, Value> * iterator;
      //typedef typename Hashtable<Key, Value> HT;

    private:
      Hashtable<Key, Value> * m_ht;
      int m_bucket;
      int m_item_in_bucket;
      int num_buckets() 
      { 
        VERIFY_OP_ON(m_ht, !=, 0, "000");
        return m_ht->data().size(); 
      }
      int num_in_bucket(int bucket) 
      { 
        VERIFY_OP_ON(m_ht, !=, 0, "001");
        return m_ht->data()[bucket].size(); 
      }
    public:


      HT_Iterator(Hashtable<Key, Value> * ht=0, bool begin=true) : m_ht(ht), m_bucket(0), m_item_in_bucket(0)
      {
        if (ht && begin)
          {
            m_bucket=0;
            m_item_in_bucket=0;

            int ibs = m_bucket;
            m_bucket = num_buckets();
            for (int ib = ibs; ib < num_buckets(); ib++)
              {
                if (num_in_bucket(ib) > 0)
                  {
                    m_bucket = ib;
                    break;
                  }
              }
          }
        if (ht && !begin)
          {
            m_bucket = num_buckets();
            m_item_in_bucket = 0;
          }
      }

      HT_Iterator& operator++()
      {
        //std::cout << "m_bucket = " << m_bucket << " num_buckets = " << num_buckets() << " m_item_in_bucket= " << m_item_in_bucket 
        //          << " num_in_bucket = " << num_in_bucket(m_bucket)
        //          << std::endl;

        VERIFY_OP_ON(num_buckets(), >, 0, "operator++ 0");

        VERIFY_OP_ON(m_bucket, < , num_buckets(), "operator++ 1");

        if (m_item_in_bucket == num_in_bucket(m_bucket) - 1)
          {
            m_item_in_bucket = 0;
            ++m_bucket;
            if (m_bucket == num_buckets()) // we're at the end
              {
              }
            else
              {
                int ibs = m_bucket;
                m_bucket = num_buckets();
                m_item_in_bucket = 0;
                for (int ib = ibs; ib < num_buckets(); ib++)
                  {
                    if (num_in_bucket(ib) > 0)
                      {
                        m_bucket = ib;
                        break;
                      }
                  }
              }
          }
        else
          {
            ++m_item_in_bucket;        
          }
        //std::cout << "m_bucket = " << m_bucket << " num_buckets = " << num_buckets() << " m_item_in_bucket= " << m_item_in_bucket  << std::endl;
        VERIFY_OP_ON( ((m_item_in_bucket == 0 && m_bucket == num_buckets()) || (m_bucket < num_buckets() && m_item_in_bucket < num_in_bucket(m_bucket)) ) , == , true, "operator++ 2");
        return *this;
      }

#if 0      
      // postfix
      HT_Iterator operator++(int)
      {
        HT_Iterator temp = *this;

        ++*this;
        return temp;
      }
#endif

      Teuchos::HashPair<Key, Value> operator*() { return m_ht->data()[m_bucket][m_item_in_bucket]; }
    
      bool operator==(const HT_Iterator& iter)
      {
        if (m_bucket == iter.m_bucket && m_item_in_bucket == iter.m_item_in_bucket)
          return true;
        return false;
      }

      bool operator!=(const HT_Iterator& iter)
      {
        return !(this->operator==(iter));
      }

    };

    template<class Key, class Value> 
    class Hashtable : public Teuchos::Teuchos_Hashtable<Key, Value>
    {
    public:
      typedef Teuchos::Teuchos_Hashtable<Key, Value> base_type;
      typedef HT_Iterator<Key, Value> iterator;
      typedef const HT_Iterator<Key, Value> const_iterator;

      typedef Teuchos::Array<Key> KeyArray;
      typedef Teuchos::Array<Value> ValueArray;

      KeyArray m_keys;
      ValueArray m_values;

      Hashtable(int capacity=101, double rehashDensity = 0.8) : base_type(capacity, rehashDensity) {}


      void initialize_iterators()
      {
      }

      void clear_iterators()
      {
      }
      
      iterator begin() 
      {
        return iterator(this, true);
      }

      iterator end() 
      {
        return iterator(this, false);
      }

    };
  }
}
#endif
