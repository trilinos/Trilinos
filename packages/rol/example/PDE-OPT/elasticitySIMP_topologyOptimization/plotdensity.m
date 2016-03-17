data_obj = importdata('density.txt', ' ', 2);  %% we need to skip the first two lines
density = data_obj.data;

xsize = round(sqrt(size(density, 1)));
ysize = xsize;

density = reshape(density, xsize, ysize);

imagesc(density.^3);
colorbar
