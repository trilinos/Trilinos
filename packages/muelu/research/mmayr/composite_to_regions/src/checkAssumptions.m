function [procsRegionAssignment,regionsProcAssignment] = checkAssumptions(whichCase,procsRegionAssignment,regionsProcAssignment,pAssign,region)


if whichCase == 'RegionsSpanProcs',
  if procsRegionAssignment(pAssign+1) == -1,
    procsRegionAssignment(pAssign+1) = region;
  else if procsRegionAssignment(pAssign+1) ~= region,
          fprintf('violates one region per proc assumption\n'); keyboard;
       end
  end
else
  if regionsProcAssignment(region+1) == -1,
    regionsProcAssignment(region+1) = pAssign;
  else if regionsProcAssignment(region+1) ~= pAssign,
          fprintf('violates one proc assigned to each region assumption\n'); keyboard;
       end
  end
end;

