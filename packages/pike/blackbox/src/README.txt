1. For hierarchichal solves, add a call to the ME interface that passes down a vector for global numbering of the system.  During the setup phase, we assign all gids to mes.  If registered in a large system, a call from the top level solve will reassign gids appropriately.

void assignGIDs(int numApps);

Discuss with Ross - unique strings fine for now.  no more than 100 apps so string lookups should not be that important.  Only really used during setup phase anyway.

2. Build a DAG from global numbers of MEs based on data transfers that return source and target MEs. Add get source and get target mes.  Is this even needed???

3. Possible change parameter and responses so that solver returns a ParameterLisbrary.  Something to handle active parameters vs inactive parameters?

Talk with Ross and Brian.  Going to experiment since this is not really needed yet.  Have a ME equivalent interface that Ross and I worked on:

    /** \brief Np */
    int getNumParamsVecs() const = 0;

    Teuchos::ArrayView<const std::string> getParamVecName(int l) const = 0;

    /** \brief 0 <= l < Np */
    void setParamVec(const Teuchos::ArrayView<const double> &p_l, int l) = 0;

    Teuchos::Ordinal getParamVecSize(int l) const = 0;

    /** \brief get Ng */
    int getNumResponseVecs() const  = 0;

    Teuchos::ArrayView<const std::string> getResponseName(int j) const = 0;

    /** \brief 0 <= j < Ng */
    void getResponse(int j, Teuchos::ArrayView<double> &g_j) const = 0;

    Teuchos::Ordinal getResponseSize(int j) const = 0;


Going to try experiment using a modified boost::any for now.

4. Rename DataTransfer to TransferOperator????  Does this matter?

5. Need to add time monitor support.

6. Multicore/threading support - need fence() operation on solve and transfer.
