function dm2hb(file,A,rhs,title,key,type,ifmt,job) 
%
% Purpose: 
% =======
%
%  Writes a sparse matrix in MATLAB into a standard Harwell-Boeing format 
%  file. The matrix is stored in compressed column format. 
%
%  This is the real double precision version.
%
% Usage: 
% ======
%
%    dm2hb('file',A,rhs,'title','key','type',ifmt,job) 
%
%    file : ASCII filename
%    A    : sparse matrix in Matlab
%    rhs  : array containing the right-hand side(s). Accessed only
%	    if job > 2 (see job)
%    title: title of matrix (72 characters)
%    key  : key of matrix   (8 characters)
%    type : type of matrix  (3 characters)
%    ifmt : specifies output format of the numerical values
%	    if ifmt < 100, then the format chosen is Dxx.yy, in which yy 
%              is precisely the integer ifmt (and xx is ifmt+7)
%	    if ifmt > 100, then the format chosen is Fxx.yy, in which 
%              the length of the mantissa yy is the integer mod(ifmt,100)
%              and the length of the integer part is ifmt/100.
%	    For examples,   
%             	   ifmt= 4   means  E11.4   [-]x.xxxxE+ee    
%	     	   ifmt=104  means  F7.4    [-]x.xxxx
%    job  : indicates whether A or A and rhs are to be out put
%           job = 1,   write srtucture only, i.e., the arrays ja and ia.
%           job = 2,   write matrix including values, i.e., a, ja, ia
%           job = 3,   write matrix and one right hand side: a,ja,ia,rhs.
% 	    job = k+2, write matrix and k successive right-hand sides.
%
% Short form: m2hb('file',A) -- use default printing format:
%	      title = 'No title';
%	      key   = 'No key';
%	      type  = 'RUA';
%	      ifmt  = 8;
%	      job   = 2;
%
% This M-file was modified and finished by Xiaoye Li at UC Berkeley. 
%
% =======================================================================

[nrow,ncol] = size( A );
nnzeros = nnz(A);
n = ncol;

% Coordinate form --> compressed column format: ja,ia,a
k = 0;
ja( 1 )  = k + 1;
for j = 1:n,  
  [rows,temp,vals] = find( A(:,j) );
  sz = size( rows ); 
  for k2 = 1:sz,
     k = k + 1;
     ia( k ) = rows(k2);
     a ( k ) = vals(k2);
  end
  ja( j+1 )  = k + 1;
end

fid = fopen(file,'w+');
if fid < 0; 
    error(['Can''t open file "' file '" for writing.']); 
end;

% Default format
if nargin < 5
    key = 'key';
    if nargin < 4
     	title = 'title';
    end;
end;
if nargin < 6
    type = 'RUA';
end;
if nargin < 7
    ifmt = 8;
end;
if nargin < 8
    job = 2;
end;

% Compute column pointer format
len = ceil( log10(0.1 + nnzeros + 1) ) + 1;
nperline = min(floor(80/len), ncol+1);
ptr_len = len;
ptr_nperline = nperline;
ptrcrd = floor(ncol / nperline) + 1;
ptrfmt = ['(' int2str(nperline) 'I' int2str(len) ')'];

% Compute row index format
nperline = min(floor(80/len), nnzeros);
ind_len = len;
ind_nperline = nperline;
indcrd = floor((nnzeros-1) / nperline) + 1;
indfmt = ['(' int2str(nperline) 'I' int2str(len) ')'];

% Compute values and rhs format (same format for both)
valcrd = 0;
rhscrd = 0;
c_valfmt = [];
if job > 1 
    if ifmt >= 100
	ihead = floor(ifmt / 100);
	ifmt = ifmt - 100*ihead;
	len = ihead + ifmt + 2;
	nperline = floor(80 / len);
	c_len = len;   			% 80 / nperline;
	for i = 1:nperline
	    c_valfmt = [c_valfmt '%' int2str(c_len) '.' int2str(ifmt) 'f'];
	end;
	valfmt = [int2str(nperline) 'F' int2str(len) '.' int2str(ifmt)];
    else
	len = ifmt + 7;
	nperline = floor(80 / len);
	c_len = len;   			% 80 / nperline;
	for i = 1:nperline
	    c_valfmt = [c_valfmt '%' int2str(c_len) '.' int2str(ifmt) 'E'];
	end;
	valfmt = [int2str(nperline) 'E' int2str(len) '.' int2str(ifmt)];
    end;
    valcrd = floor((nnzeros-1) / nperline) + 1;
    valfmt = ['(' valfmt ')'];
    c_valfmt = [c_valfmt '\n'];
end;

nrhs = job - 2;
if nrhs >= 1
    rhscrd = floor((nrhs*nrow - 1) / nperline) + 1;
end;

totcrd = ptrcrd + indcrd + valcrd + rhscrd;


% Now write 4-line header 

% Line 1
t = title; m = size(t,2);
for i = m+1:72, t = [t ' ']; end;
fprintf(fid, '%72s', t);
t = key; m = size(t,2);
for i = m+1:8, t = [t ' ']; end;
fprintf(fid,'%8s\n', t);

% Line 2
fprintf(fid,'%14d%14d%14d%14d%14d\n', totcrd, ptrcrd, indcrd, valcrd, rhscrd);

% Line 3
t = type; m = size(t,2);
for i=m+1:14, t = [t ' ']; end;
fprintf(fid, '%14s', t);
fprintf(fid, '%14i%14i%14i%14i\n', nrow, ncol, nnzeros, nrhs);

% Line 4
t = ptrfmt; m = size(t,2);
for i=m+1:16, t = [t ' ']; end;
fprintf(fid, '%16s', t);
t = indfmt; m = size(t,2);
for i=m+1:16, t = [t ' ']; end;
fprintf(fid, '%16s', t);
t = valfmt; m = size(t,2);
for i=m+1:20, t = [t ' ']; end;
fprintf(fid, '%20s', t);
fprintf(fid, '%20s\n', t);

% column pointers
t = [];
for j = 1:ptr_nperline, t = [t '%' int2str(ptr_len) 'd']; end;
t = [t '\n'];
fprintf(fid, t, ja(1:ncol+1)); 
if ncol+1 > ptr_nperline & ...
	(ncol+1)/ptr_nperline ~= floor((ncol+1)/ptr_nperline)
    fprintf(fid,'\n');
end;

% row indices
t = [];
for j = 1:ind_nperline, t = [t '%' int2str(ind_len) 'd']; end;
t = [t '\n']; 
fprintf(fid, t, ia(1:nnzeros)); 
if nnzeros > ind_nperline & ...
	nnzeros/ind_nperline ~= floor(nnzeros/ind_nperline)
    fprintf(fid,'\n');
end;

% numerical values of nonzero elements of the matrix
if job >= 2
    if job == 2
    	fprintf(fid, c_valfmt, a(1:nnzeros));
    else
    	fprintf(fid, c_valfmt, a(1:nnzeros));
	if nnzeros > nperline & ...
		nnzeros/nperline ~= floor(nnzeros/nperline)
	    fprintf(fid,'\n');
	end;
	fprintf(fid, c_valfmt, rhs(1:nrhs*nrow));
    end;
end;

fclose(fid);
