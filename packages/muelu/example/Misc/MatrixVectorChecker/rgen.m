for kk = 1:10
n = 1001; m = n; den = 1/n;
mat = sprand(n,m,den);
[aaa,bbb,ccc] = find(mat);
[aaa,perm] = sort(aaa);
bbb = bbb(perm);
ccc = ccc(perm);
first = 1;
row   = aaa(first);
for i=2:length(aaa)
   if  aaa(i) ~= row,
      last = i-1;
      [ttt,perm] = sort(bbb(first:last));
      temp = aaa(first:last); aaa(first:last) = temp(perm);
      temp = ccc(first:last); ccc(first:last) = temp(perm);
      first = i; row = aaa(first);
   end;
end;
last = length(aaa);
[ttt,perm] = sort(bbb(first:last));
temp = aaa(first:last); aaa(first:last) = temp(perm);
temp = ccc(first:last); ccc(first:last) = temp(perm);

fid = fopen('mat_example.mm','w');
fprintf(fid,'%%%%MatrixMarket matrix coordinate real general\n');
fprintf(fid,'%d %d %d\n',n,m,length(aaa));
for i=1:length(aaa)
   fprintf(fid,'%6d %6d %20.13e\n',aaa(i),bbb(i),ccc(i));
end;
fclose(fid);
!mpirun -np 11 ml_read_MatrixMarket.exe -f mat_example.mm
!rm -f mat_example.mm
end;
