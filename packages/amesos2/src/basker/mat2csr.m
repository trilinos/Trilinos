function out = mat2csr(A)

[n,m] = size(A);
col_ptr = [1];
row_idx = [];
vals =[];
count = 1;

fid = fopen('out.csc', 'w');
fprintf(fid, '%d %d \n',n, nnz(A)); 
fprintf(1, '%d %d \n',n, nnz(A)); 


for i=1:n

    [i,j,val] = find(A(:,i));
    count = count + length(val);
    col_ptr = [col_ptr; count];
    row_idx = [row_idx; i];
    vals = [vals; val];

end
%check
%col_ptr
size(col_ptr)
size(row_idx)
size(vals)

for i=1:length(col_ptr)
    fprintf(fid,'%d ', col_ptr(i));
    fprintf(1,'%d ', col_ptr(i));
end
fprintf('\n');

for i=1:length(row_idx)
    fprintf(fid, '%d ' , row_idx(i));
     fprintf(1, '%d ' , row_idx(i));
end
fprintf('\n');
for i = 1:length(vals)
    fprintf(fid, '%f ', vals(i));
end


    
