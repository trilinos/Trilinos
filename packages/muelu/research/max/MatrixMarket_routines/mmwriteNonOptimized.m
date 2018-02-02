function [ err ] = mmwriteNonOptimized(filename,A,comment,field,precision)
%
% Function: mmwrite(filename,A,comment,field,precision)
%
%    Writes the sparse or dense matrix A to a Matrix Market (MM)
%    formatted file.
%
% Required arguments:
%
%                 filename  -  destination file
%
%                 A         -  sparse or full matrix
%
% Optional arguments:
%
%                 comment   -  matrix of comments to prepend to
%                              the MM file.  To build a comment matrix,
%                              use str2mat. For example:
%
%                              comment = str2mat(' Comment 1' ,...
%                                                ' Comment 2',...
%                                                ' and so on.',...
%                                                ' to attach a date:',...
%                                                [' ',date]);
%                              If ommitted, a single line date stamp comment
%                              will be included.
%
%                 field     -  'real'
%                              'complex'
%                              'integer'
%                              'pattern'
%                              If ommitted, data will determine type.
%
%                 precision -  number of digits to display for real
%                              or complex values
%                              If ommitted, full working precision is used.
%

if ( nargin == 5)
  precision = 16;
elseif ( nargin == 4)
  precision = 16;
elseif ( nargin == 3)
  mattype = 'real'; % placeholder, will check after FIND-ing A
  precision = 16;
elseif ( nargin == 2)
  comment = '';
  % Check whether there is an imaginary part:
  mattype = 'real'; % placeholder, will check after FIND-ing A
  precision = 16;
end

mmfile = fopen([filename],'w');
if ( mmfile == -1 )
 error('Cannot open file for output');
end;


[M,N] = size(A);

%%%%%%%%%%%%%       This part for sparse matrices     %%%%%%%%%%%%%%%%
if ( issparse(A) )

  [I,J,V] = find(A);
  if ( sum(abs(imag(nonzeros(V)))) > 0 )
    Vreal = 0;
  else
    Vreal = 1;
  end

  if ( ~ strcmp(mattype,'pattern') & Vreal )
    mattype = 'real';
  elseif ( ~ strcmp(mattype,'pattern') )
    mattype = 'complex';
  end
%
% Determine symmetry:
%

    symm = 'general';
    issymm = 0;
    NZ = length(V);

% Sparse coordinate format:

  rep = 'coordinate';


  fprintf(mmfile,'%%%%MatrixMarket matrix %s %s %s\n',rep,mattype,symm);
%  [MC,NC] = size(comment);
%  if ( MC == 0 )
%    fprintf(mmfile,'%% Generated %s\n',[date]);
%  else
%    for i=1:MC,
%      fprintf(mmfile,'%%%s\n',comment(i,:));
%    end
%  end
  fprintf(mmfile,'%d %d %d\n',M,N,NZ);
  cplxformat = sprintf('%%d %%d %% .%dg %% .%dg\n',precision,precision);
  realformat = sprintf('%%d %%d %% .%dg\n',precision);
  if ( strcmp(mattype,'real') )
     for i=1:NZ
        fprintf(mmfile,realformat,I(i),J(i),V(i));
     end;
  elseif ( strcmp(mattype,'complex') )
  for i=1:NZ
     fprintf(mmfile,cplxformat,I(i),J(i),real(V(i)),imag(V(i)));
  end;
  elseif ( strcmp(mattype,'pattern') )
     for i=1:NZ
        fprintf(mmfile,'%d %d\n',I(i),J(i));
     end;
  else
     err = -1;
     disp('Unsupported mattype:')
     mattype
  end;

%%%%%%%%%%%%%       This part for dense matrices      %%%%%%%%%%%%%%%%
else
  if ( sum(abs(imag(nonzeros(A)))) > 0 )
    Areal = 0;
  else
    Areal = 1;
  end
  if ( ~strcmp(mattype,'pattern') & Areal )
    mattype = 'real';
  elseif ( ~strcmp(mattype,'pattern')  )
    mattype = 'complex';
  end
%
% Determine symmetry:
%

    issymm = 0;
    symm = 'general';


% Dense array format:

  rep = 'array';
  [MC,NC] = size(comment);
  fprintf(mmfile,'%%%%MatrixMarket matrix %s %s %s\n',rep,mattype,symm);
  for i=1:MC,
    fprintf(mmfile,'%%%s\n',comment(i,:));
  end;
  fprintf(mmfile,'%d %d\n',M,N);
  cplxformat = sprintf('%% .%dg %% .%dg\n', precision,precision);
  realformat = sprintf('%% .%dg\n', precision);
  if ( ~ strcmp(symm,'general') )
     rowloop = 'j';
  else
     rowloop = '1';
  end
  if ( strcmp(mattype,'real') )
     for j=1:N
       for i=eval(rowloop):M
          fprintf(mmfile,realformat,A(i,j));
       end
     end
  elseif ( strcmp(mattype,'complex') )
     for j=1:N
       for i=eval(rowloop):M
          fprintf(mmfile,cplxformat,real(A(i,j)),imag(A(i,j)));
       end
     end
  elseif ( strcmp(mattype,'pattern') )
     err = -2
     disp('Pattern type inconsistent with dense matrix')
  else
     err = -2
     disp('Unknown matrix type:')
     mattype
  end
end

fclose(mmfile);
