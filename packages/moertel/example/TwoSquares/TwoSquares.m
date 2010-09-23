% creates the matlab data structure for the following problem:
%
%  +------+
%  |  S2  |
%  +------+
%  +------+
%  |  S1  |
%  +------+
%
% where S1 = (-1,1) x (-1,1) and S2 = (-1, 1) x (1, 3).
% BC label is set to 0 everywhere except for the upper side of S1
% and the lower side of S2

% Creates S1
[p1,e1,t1] = poimesh('squareg', 10, 10);
% put `10' on the upper side, 0 otherwise
for i=1:size(e1,2)
  n1 = e1(1,i);
  n2 = e1(2,i);
  if (p1(2,n1) == 1.) && (p1(2,n2) == 1.0)
    e1(5,i) = 10;
  else
    e1(5,i) = 0;
  end
end

% Creates S2
[p2,e2,t2] = poimesh('squareg', 9, 9);
% translate by a factor 2 on the Y-axix
p2(2,:) = p2(2,:) + 2;

% put `10' on the lower side, 0 otherwise
for i=1:size(e2,2)
  n1 = e2(1,i);
  n2 = e2(2,i);
  if (p2(2,n1) == 1.) && (p2(2,n2) == 1.0)
    e2(5,i) = 20;
  else
    e2(5,i) = 0;
  end
end

% finally glue them together
p = [p1(1,:), p2(1,:); p1(2,:), p2(2,:)];
t2(1:3,:) = t2(1:3,:) + size(p1,2);
t2(4,:) = 2; 
t = [t1(1,:), t2(1,:);
     t1(2,:), t2(2,:);
     t1(3,:), t2(3,:);
     t1(4,:), t2(4,:)];
e2(1:2,:) = e2(1:2,:) + size(p1,2);
e = [e1(1,:), e2(1,:);
     e1(2,:), e2(2,:);
     e1(3,:), e2(3,:);
     e1(4,:), e2(4,:);
     e1(5,:), e2(5,:);
     e1(6,:), e2(6,:);
     e1(7,:), e2(7,:)];

pdeplot(p,e,t)
convert(p,e,t,'TwoSquares.grid');
