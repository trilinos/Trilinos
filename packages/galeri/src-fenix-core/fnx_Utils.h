#ifndef FNX_UTILS_H
#define FNX_UTILS_H

namespace fnx 
{
  class Error
  {
    public:
      Error(const string& what, const string& file, int line) :
        what_(what),
        file_(file),
        line_(line)
      {}

      void Print(ostream& os)
      {
        os << "fnx::error::" << what_ << endl;
        os << "fnx::error:: File " << file_ << ", line " << line_ << endl;
      }

    private:
      const string& what_, file_;
      int line_;
  };

  template<class T>
  class DistArray : public Epetra_DistObject
  {
    public:
      DistArray(Epetra_Map& Map, const int RowSize) :
        Epetra_DistObject(Map)
      {
        // only Epetra_Map's with constant element size of 1 are allowed
        if (Map_.MaxElementSize() != 1) 
          FNX_THROW("Map.MaxElementSize() != 1");
        if (!Map_.ConstantElementSize()) 
          FNX_THROW("Map.ConstantElementSize() != true");

        MyLength_     = Map.NumMyElements();
        GlobalLength_ = Map.NumGlobalElements();
        RowSize_      = RowSize;
        count_        = 0;
        values_.resize(MyLength_ * RowSize_);
      }

      inline int MyLength() const
      {
        return(MyLength_);
      }

      inline int GlobalLength() const
      {
        return(GlobalLength_);
      }

      inline int RowSize() const
      {
        return(RowSize_);
      }

      inline T& operator()(const int LEID, const int ID, const bool isLocal)
      {
        assert (isLocal ==true);
        assert (ID <= RowSize_);
        return(values_[LEID * RowSize_ + ID]);
      }

      inline T& operator()(const int GEID, const int ID)
      {
        int LEID = Map().LID(GEID);
        assert (LEID != -1);
        assert (ID <= RowSize_);
        return(values_[LEID * RowSize_ + ID]);
      }

      void Print(ostream& os) const
      {
        os << "DistArray object, label   = " << Label() << endl;
        os << "Number of local elements  = " << NumMyElements() << endl;
        os << "Number of global elements = " << NumGlobalElements() << endl;
        os << endl;

        for (int iproc=0; iproc < Comm().NumProc(); iproc++) 
        {
          if (iproc == 0)
          {
            os << "GEID\t";
            for (int i = 0; i < RowLength(); ++i) os << "V\t";
            os << endl;
          }

          if (Comm().MyPID() == iproc) 
          {
            for (int i = 0; i < NumMyElements(); ++i)
            {
              os << Map().GID(i) << '\t';
              for (int j = 0; j < RowSize_; ++j)
                os << values_[i * RowSize_ + j] << '\t';
              os << endl;
            }
          }
        }
        os << endl;
      }

      int NextGID()
      {
        ++count_;
        if (count_ < NumMyElements_)
          return(Map().GID(count_));
        else
          return(-1);
      }

      int FirstGID()
      {
        count_ = 0;
        return(Map().GID(0));
      }

    private:
      virtual int CheckSizes(const Epetra_SrcDistObject& Source)
      {
        return(0);
      }

      virtual int CopyAndPermute(const Epetra_SrcDistObject& Source,
                                 int NumSameIDs,
                                 int NumPermuteIDs,
                                 int * PermuteToLIDs,
                                 int * PermuteFromLIDs,
                                 const Epetra_OffsetIndex * Indexor)
      {
        const DistArray& S = dynamic_cast<const DistArray&>(Source);
        const vector<T>& From = S.ExtractView();

        vector<T>& To = values_;

        int * ToFirstPointInElementList = 0;
        int * FromFirstPointInElementList = 0;
        int * FromElementSizeList = 0;

        int j, jj, jjj, k;

        int NumSameEntries;

        NumSameEntries = NumSameIDs;

        // Short circuit for the case where the source and target vector is the same.
        if (To==From) NumSameEntries = 0;

        // Do copy first
        if (NumSameIDs>0)
          if (To!=From) {
            for (j=0; j<NumSameEntries * RowSize_; j++)
            {
              To[j] = From[j];
            }
          }

        // Do local permutation next
        if (NumPermuteIDs>0) {

          for (j=0; j<NumPermuteIDs * RowSize_; j++) 
            To[PermuteToLIDs[j]] = From[PermuteFromLIDs[j]];
          // constant element size case
        }

        return(0);
      }

      virtual int PackAndPrepare(const Epetra_SrcDistObject& Source,
                                 int NumExportIDs,
                                 int* ExportLIDs,
                                 int& LenExports,
                                 char*& Exports,
                                 int& SizeOfPacket,
                                 int* Sizes,
                                 bool & VarSizes,
                                 Epetra_Distributor& Distor)
      {
        const DistArray& S = dynamic_cast<const DistArray&>(Source);
        const vector<T>& From = S.ExtractView();

        vector<T> To = values_;

        int * FromFirstPointInElementList = 0;
        int * FromElementSizeList = 0;

        SizeOfPacket = RowSize_ * sizeof(T); 

        if(NumExportIDs*SizeOfPacket>LenExports) {
          if (LenExports>0) delete [] Exports;
          LenExports = NumExportIDs*SizeOfPacket;
          Exports = new char[LenExports];
        }

        T* ptr;

        if (NumExportIDs>0) {
          ptr = (T*) Exports;

          // Point entry case
          for (int j=0; j<NumExportIDs; j++) 
            for (int k = 0; k < RowSize_ ; ++k)
              *ptr++ = From[ExportLIDs[j] * RowSize_ + k];
        }

        return(0);
      }

      virtual int UnpackAndCombine(const Epetra_SrcDistObject& Source,
                                   int NumImportIDs,
                                   int* ImportLIDs,
                                   int LenImports,
                                   char* Imports,
                                   int& SizeOfPacket,
                                   Epetra_Distributor& Distor,
                                   Epetra_CombineMode CombineMode,
                                   const Epetra_OffsetIndex * Indexor)
      {
        int j, jj, k;

        if (CombineMode != Insert)
          EPETRA_CHK_ERR(-1); //Unsupported CombinedMode, will default to Zero

        cout << NumImportIDs << endl;
        if (NumImportIDs<=0) return(0);

        T* To = &values_[0];
        int * ToFirstPointInElementList = 0;
        int * ToElementSizeList = 0;

        T* ptr;
        // Unpack it...

        ptr = (T*) Imports;

        for (j=0; j<NumImportIDs; j++) 
          for (int k = 0; k < RowSize_ ; ++k)
            To[ImportLIDs[j] * RowSize_ + k] = *ptr++;

        return(0);
      }

      const vector<T>& ExtractView() const
      {
        return(values_);
      }

      vector<T> values_;
      int MyLength_;
      int GlobalLength_;
      int RowSize_;
      int count_;
      int last_;

  };

} // namespace fnx

#endif
