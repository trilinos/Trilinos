fname = 'vdp.out'

fid = fopen(fname)

data = fscanf(fid,'%g',[3,Inf])';

plot(data(:,2),data(:,3),'b-*')
