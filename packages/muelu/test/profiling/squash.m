function [newvec,newlabels] = squash(vec,labels,thresh,maxEntries)
  if (exist('thresh') ~= 1)
    thresh = 0.05;
  end
  if (exist('maxEntries') ~= 1)
    maxEntries = 10;
  end
  newleng = length(find(vec>thresh))+1;
  % limit number of entries in vector (as this is the number of slices in pie chart)
  if newleng > maxEntries, newleng=maxEntries; end

  newvec = zeros(newleng,1);
  newvec(1:newleng-1) = vec(1:newleng-1);
  newvec(newleng) = sum(vec(newleng:length(vec)));
  newlabels = cell(newleng,1);
  newlabels(1:newleng-1) = labels(1:newleng-1);
  newlabels(newleng) = {'misc'};
  total = sum(newvec);
  checksum = 0;
  for ii=1:length(newlabels)
    perc = sprintf('%2.1f%% (%3.1f s)',100*newvec(ii)/total,newvec(ii));
    newlabels(ii) = {[char(newlabels(ii)) ' ' perc]};
    checksum = checksum + newvec(ii);
  end
  if (abs(checksum - total) > 1e-12)
    msg = sprintf('sum of slices (%4.3f) do not add up to total (%4.3f)',checksum,total);
    warning(msg);
  end
  % labelsPlusTimings = [char(newlabels)  num2str(newT)];
