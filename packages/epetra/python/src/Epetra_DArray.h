class DArray
{
  public:
    DArray(int size)
    {
      ptr_ = new double[size];
      size_ = size;
      ownership_ = true;
    }

    DArray(double* ptr, int size)
    {
      ptr_ = ptr;
      size_ = size;
      ownership_ = false;
    }

    ~DArray()
    {
      if (ownership_)
        delete[] ptr_;
    }

    double& operator[](int i)
    {
      return(ptr_[i]);
    }

    double* Values()
    {
      return(ptr_);
    }

  private:
    double* ptr_;
    int size_;
    bool ownership_;
};
