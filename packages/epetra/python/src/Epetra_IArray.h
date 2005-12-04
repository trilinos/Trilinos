class IArray
{
  public:
    IArray(int size)
    {
      ptr_ = new int[size];
      size_ = size;
      ownership_ = true;
    }

    IArray(int* ptr, int size)
    {
      ptr_ = ptr;
      size_ = size;
      ownership_ = false;
    }

    ~IArray()
    {
      if (ownership_)
        delete[] ptr_;
    }

    int& operator[](int i)
    {
      return(ptr_[i]);
    }

    int* Values()
    {
      return(ptr_);
    }

  private:
    int* ptr_;
    int size_;
    bool ownership_;
};
