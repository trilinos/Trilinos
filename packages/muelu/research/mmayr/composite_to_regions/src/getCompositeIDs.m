function [myList] = getCompositeIDs(myNodes,myRank) 

%
% Returns the global IDs of nodes that myRank owns in the composite map.
%
% Note: the first entry in any sublist procIds(1) corresponds to the owner
%       of any shared node.
   myList = [];
   count = 0;
   for i=1:length(myNodes),
      if  myNodes(i).procIds(1) == myRank % I own this guy
         count = count + 1;
         myList(count) = myNodes(i).ID;
      end;
   end;
