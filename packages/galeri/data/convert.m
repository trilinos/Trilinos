function convert(p, e, t, output)

fid = fopen(output, 'w');

fprintf(fid, 'NumDimensions: 2\n');
fprintf(fid, 'NumVerticesPerElement: 3\n');
fprintf(fid, 'NumVerticesPerFace: 2\n');
fprintf(fid, 'NumFacesPerElement: 3\n');
fprintf(fid, 'NumMyElements: %d\n', size(t, 2));
fprintf(fid, 'NumMyVertices: %d\n', size(p, 2));
fprintf(fid, 'NumMyBoundaryFaces: %d\n', size(e, 2));
fprintf(fid, 'ElementType: GALERI_TRIANGLE\n');

% print coord, last is the tag
fprintf(fid, 'VertexCoord:\n');
for i=1:size(p, 2)
  fprintf(fid, '%e %e 0\n', p(1,i), p(2,i));
end

% print elements, last is the tag
fprintf(fid, 'ElementVertices:\n');
for i=1:size(t, 2)
  fprintf(fid, '%d %d %d %d\n', t(1,i)-1, t(2,i)-1, t(3,i)-1, t(4,i));
end

% print faces, last is the tag/patch
fprintf(fid, 'FaceVertices:\n');
for i=1:size(e, 2)
  fprintf(fid, '%d %d %d\n', e(1,i)-1, e(2,i)-1, e(5,i));
end

fprintf(fid, 'End\n');

fclose(fid);
