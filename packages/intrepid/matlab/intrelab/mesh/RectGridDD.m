
function [mesh] = RectGridDD(nsdx, nsdy, t, p, e)
%
%  AUTHOR:  Matthias Heinkenschloss and Denis Ridzal
%           Department of Computational and Applied Mathematics
%           Rice University
%           November 23, 2005

xmin = min(p(:,1));
xmax = max(p(:,1));
ymin = min(p(:,2));
ymax = max(p(:,2));

elempart = zeros(size(t,1),1);

for k = 1:size(t,1)
    x = p(t(k,:),1);  % x coordinates of all vertives in triangle k
    y = p(t(k,:),2);  % y coordinates of all vertives in triangle k
    isd = 0;  % subdomain number
    for i = 1:nsdx
        for j = 1:nsdy
            isd = isd+1;
            if( all( x >= xmin+(i-1)*(xmax-xmin)/nsdx-1.e-10 ) ...
                & all( x <= xmin+i*(xmax-xmin)/nsdx+1.e-10 ) ...
                & all( y >= ymin+(j-1)*(ymax-ymin)/nsdy-1.e-10 ) ...
                & all( y <= ymin+j*(ymax-ymin)/nsdy+1.e-10 ) )
                % all vertices of triangle k are inside subdomain k,
                % i.e. triangke k is inside subdomain k
                elempart(k) = isd-1;
            end
        end
    end
end
       
[mesh] = meshpart(t, p, e, elempart);


% plot the partition
color = ['w','r','g','b','c','m','y'];
nc    = length(color);

if 0
hold on;
axis('equal','tight');
for i = 1:size(t,1)                  % plot each element
    x = p(t(i,1:3),1);
    y = p(t(i,1:3),2);
    c = mod(elempart(i),nc) + 1;
    fill(x,y,color(c));
    plot(x, y, '-k');
end
hold off;
end
