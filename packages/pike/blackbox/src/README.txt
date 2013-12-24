1. For hierarchichal solves, add a call to the ME interface that passes down a vector for global numbering of the system.  During the setup phase, we assign all gids to mes.  If registered in a large system, a call from the top level solve will reassign gids appropriately.

void assignGIDs(int numApps);

2. Build a DAG from global numbers of MEs based on data transfers that return source and target MEs. Add get source and get target mes.

3. Possible change parameter and responses so that solver returns a PArameterLisbrary.  Something to handle active parameters vs inactive parameters?

4. Rename DataTransfer to TransferOperator????  Does this matter?

