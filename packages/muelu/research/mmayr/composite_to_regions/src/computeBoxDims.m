function [procsCorner,globalDims,localDims,regionCorner] = computeBoxDims(rp1,pp1,i,j,...
          procsCorner,globalDims,localDims,regionCorner)
         % rp1 is region plus one
         % pp1 is proc plus one

         % this function computes some box information associated with a grid point,
         % region, and processor. This information corresponds to the lower left corner,
         % global dimensions of the region, and the local dimension of the region that
         % this processor owns. For the lower left corner, we compute the true regionCorner 
         % (within the absolute global 2D mesh that spans all regions). That is, the
         % location of a region's lower-left corner, regardless of what processor owns it.
         % We also compute the lower left corner that is held by a processor. That is,
         % when multiple processors share a region, we need to know the location where
         % this processor's piece of the region starts (again within the absolute 2D
         % mesh that spans all regoins). Later, in mk2DRegionFile, we'll compute relative
         % locations. That is, where does a processor's portion of a sub-region start
         % (using relative locations within just the region).

         % NOTE: this function DEFINITELY assumes that it is invoked in a lexicographical
         % fashion. So, it is first invoked for i=0,j=0, then for i=1,j=0, etc. 

% A corner has not yet been assigned, so this i,j must be the lower left corner.
% The first 'if statement' below does not depend on processors (only regions)
% The second if is processor specific (so this must be the lower corner of this 
% processors portion of the region rp1-1. 
if regionCorner(rp1,1) == -1,
   regionCorner(rp1,1)=i; regionCorner(rp1,2)=j;
end;
if procsCorner(rp1,pp1,1) == -1,
   procsCorner(rp1,pp1,1)=i; procsCorner(rp1,pp1,2)=j; 
end;

inds = find( procsCorner(rp1,:,1) ~= -1);
cornerx = min(procsCorner(rp1,inds,1));
cornery = min(procsCorner(rp1,inds,2));
if cornerx ~= regionCorner(rp1,1), fprintf('ugh\n'); keyboard; end;
if cornery ~= regionCorner(rp1,2), fprintf('UGH\n'); keyboard; end;
if i-cornerx+1 > globalDims(rp1,1),
   globalDims(rp1,1)= i-cornerx+1;
end
if j-cornery+1 > globalDims(rp1,2),
   globalDims(rp1,2)= j-cornery+1;
end
if i-procsCorner(rp1,pp1,1)+1 > localDims(rp1,pp1,1),
   localDims(rp1,pp1,1)= i-procsCorner(rp1,pp1,1)+1;
end
if j-procsCorner(rp1,pp1,2)+1 > localDims(rp1,pp1,2),
   localDims(rp1,pp1,2)= j-procsCorner(rp1,pp1,2)+1;
end       
