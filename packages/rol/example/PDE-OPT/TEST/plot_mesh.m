S.Vertices = load('nodes.txt');
S.Faces = load('cells.txt') + 1;
S.FaceColor = 'white';
S.EdgeColor = 'red';
S.LineWidth = 3;

patch(S);

set(gca, 'Ylim', [-0.1 1.1]);
axis fill;

%idx = [0:1:size(S.Faces, 1)-1];
%idx_str = num2str(idx);
%idx_cay = strsplit(idx_str);
%text(S.Vertices(S.Faces(:,1), 1), S.Vertices(S.Faces(:,1), 2), idx_cay, 'FontSize', 14)

idx = load('sideset.txt');
idx_str = num2str(idx');
idx_cay = strsplit(idx_str);

text(S.Vertices(S.Faces(idx+1,1), 1), S.Vertices(S.Faces(idx+1,1), 2), idx_cay, 'FontSize', 14)

