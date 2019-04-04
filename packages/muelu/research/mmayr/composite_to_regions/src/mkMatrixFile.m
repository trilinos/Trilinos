% writeMatrixFile
%
% Write matrix to file 'Amat.mm'
%
% Input:
%   A       matrix to be written
%   outDir  path to output directory
%

% now send matrix to C++

function [] = writeMatrixFile(A, outDir)

% open file
fp = fopen(sprintf('%s/Amat.mm', outDir), 'w');
if (fp == -1)
  error('Cannot open file.');
end

% write to file
fprintf(fp,'%%%%MatrixMarket matrix coordinate real general\n');
fprintf(fp,'%d %d %d\n',size(A,1),size(A,2),nnz(A));
[aaa,bbb,ccc] = find(A);
[~,ii] = sort(aaa);
aaa = aaa(ii); bbb = bbb(ii); ccc = ccc(ii);
row = aaa(1);
startRow = 1;  endRow = startRow;
for i=2:length(aaa)
  if aaa(i) ~= row % sort and print out the previous row
    aa = aaa(startRow:i-1);
    bb = bbb(startRow:i-1);
    cc = ccc(startRow:i-1);
    [~,ii] = sort(bb);
    aa = aa(ii);
    bb = bb(ii);
    cc = cc(ii);
    for k=1:length(aa)
      fprintf(fp,'%d %d %20.13e\n',aa(k),bb(k),cc(k));
    end
    startRow=i;   % indicates that we are now storing a new row
  end
end
% print the last row
aa = aaa(startRow:end);
bb = bbb(startRow:end);
cc = ccc(startRow:end);
[~,ii] = sort(bb);
aa = aa(ii);
bb = bb(ii);
cc = cc(ii);
for k=1:length(aa)
  fprintf(fp,'%d %d %20.13e\n',aa(k),bb(k),cc(k));
end

% close file
fclose(fp);

end