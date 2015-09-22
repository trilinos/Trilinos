function [mesh] = exoread(filename)
%
% Cubit / Exodus reader
%
% INPUT:   filename              = name of ExodusII file generated through Cubit
%
% OUTPUT: 'mesh'                 = mesh data structure
%               .t               = element connectivity
%               .p               = node coordinates
%               .elem_ss{i}      = global element IDs on sideset (boundary) i
%               .side_ss{i}      = local side (subcell) IDs on sideset (boundary) i
%               .sidenodes_ss{i} = global node ID tuples on each side of sideset i
%               .e               = global node ID tuples for all sides in all sidesets
%               .cellType        = string for cell type ('Triangle', etc.)
%               .sideType        = string for side subcell type ('Line', etc.)
%               .sidesPerCell    = number of sides per parent cell
%
% Author: Denis Ridzal
%         Sandia National Labs

  ncid=netcdf.open(filename,'NC_NOWRITE');
  [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
  exomesh = containers.Map();
  for i=1:numvars-1
      [varname vartype vardimIDs varatts] = netcdf.inqVar(ncid,i);
      varid = netcdf.inqVarID(ncid,varname);
      exomesh(varname) = netcdf.getVar(ncid, varid);
  end

  % Set root strings for cell choice, coordinates, connectivity, etc.
  str_coord    = 'coord';
  str_connect  = 'connect1';
  str_elem_ss  = 'elem_ss';
  str_side_ss  = 'side_ss';

  % Define local connectivities; this is ExodusII format.
  tricon  = int32(  [1 2; 2 3; 3 1]  );
  quadcon = int32(  [1 2; 2 3; 3 4; 4 1]  );
  tetcon  = int32(  [1 2 4; 2 3 4; 1 4 3; 1 3 2]  );
  hexcon  = int32(  [1 2 6 5; 2 3 7 6; 3 4 8 7; 1 5 8 4; 1 4 3 2; 5 6 7 8]  );

  % Grab nodes.
  mesh.p = exomesh(str_coord);

  % Grab node connectivity.
  mesh.t = double(exomesh(str_connect))';

  % Determine cell type.
  if (size(mesh.p,2) == 2) & (size(mesh.t,2) == 3)
    str_cell = 'Triangle';  str_sidecell = 'Line';
  elseif (size(mesh.p,2) == 2) & (size(mesh.t,2) == 4)
    str_cell = 'Quadrilateral';  str_sidecell = 'Line';
  elseif (size(mesh.p,2) == 3) & (size(mesh.t,2) == 4)
    str_cell = 'Tetrahedron';  str_sidecell = 'Triangle';
  elseif (size(mesh.p,2) == 3) & (size(mesh.t,2) == 8)
    str_cell = 'Hexahedron';  str_sidecell = 'Quadrilateral';
  else
    error('\nexoread: Unknown cell type!\n')
  end
  mesh.cellType = str_cell;
  mesh.sideType = str_sidecell;

  switch lower(str_cell)
    case {'triangle'}
      cellcon = tricon;
    case {'quadrilateral'}
      cellcon = quadcon;
    case {'tetrahedron'}
      cellcon = tetcon;
    case {'hexahedron'}
      cellcon = hexcon;
    otherwise
      error('\nexoread: Unknown cell type string!\n')
  end
  mesh.sidesPerCell = size(cellcon, 1);

  % Grab sidesets (elem_ss) with local boundary side IDs (side_ss).
  % Also compute global side nodes associated with sidesets.
  nss = 1;
  mesh.e = [];
  elem_ss_key = strcat(str_elem_ss, '1');
  side_ss_key = strcat(str_side_ss, '1');
  while exomesh.isKey(elem_ss_key)
    mesh.elem_ss{nss} = double(exomesh(elem_ss_key));
    mesh.side_ss{nss} = double(exomesh(side_ss_key))-1; % subtract 1 for zero-based local indexing

    local_sides = cellcon(mesh.side_ss{nss}+1, :);
    local_sides_flat = double(reshape(local_sides, numel(local_sides), 1));
    elem_repeat = repmat(mesh.elem_ss{nss}, size(cellcon,2), 1);
    linearInd = sub2ind(size(mesh.t), elem_repeat, local_sides_flat);
    global_side_flat = mesh.t(linearInd);
    mesh.sidenodes_ss{nss} = reshape(global_side_flat, size(mesh.side_ss{nss},1), size(cellcon,2));
    mesh.e = [mesh.e; mesh.sidenodes_ss{nss}];

    nss = nss+1;
    elem_ss_key = strcat(str_elem_ss, num2str(nss));
    side_ss_key = strcat(str_side_ss, num2str(nss));
  end

  % Close netCDF file.
  netcdf.close(ncid);

end % function exoread
