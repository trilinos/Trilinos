function plotresults(M,N,T,width,height)

data = load('design.txt');
W    = [2.70, 8.96, 7.75];
X = [];
Y = [];
Z = [];
for i=1:M
  A  = data(((i-1)*(M+N)+1):(i*(M+N)),2:4);
  x0 = (i-1)*width/M;
  x1 = i*width/M;
  for j=1:N
    B  = A(((j-1)*T+1):(j*T),3);
    y0 = (j-1)*height/N;
    y1 = j*height/N;
    X  = [X, [x0;x0;x1;x1]];
    Y  = [Y, [y0;y1;y1;y0]];
    Z  = [Z, W*B];
  end
end

figure,
patch(X,Y,Z)
axis('equal','tight')
colorbar
