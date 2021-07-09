clear all;

mesh = 0; wave = 5; damp = 0; impd = 1; fullObs = false; fullCtrl = false; solver = 'PRDW';
%% mesh     = 0, 1, 2, 3, 4 is mesh resolution from coarsest to finest
%%            NOTE: Meshes are available only for the cases 0 and 1 in the repository.
%%            For finer meshes, consult the /mesh/README.md file.
%% wave     = 5, 10, 20, 40, 80 is wavenumber
%%            NOTE: The minimum mesh resolutions that support these wavenumbers
%%            are 0, 1, 2, 3, and 4, respectively.
%% damp     = 0, 1 is to turn off/on damping
%% impd     = 0, 1 is to select reflective/radiating boundary conditions
%% fullObs  = false, true is to force full state observations, i.e., everywhere in the domain
%% fullCtrl = false, true is to force distributed controls, i.e., everywhere in the domain
%% solver   = 'DIRECT'   for direct KKT solve
%%            'PRDW'     for (perturbed) Rees-Dollar-Wathen preconditioner with SYMMLQ
%%            'IMSHIFT'  for imaginary shift preconditioner with SYMMLQ

file = ['input_mesh' int2str(mesh) '_wave' int2str(wave) '_damp' int2str(damp) '_imp' int2str(impd) '.xml'];
command = ['cp inputs/' file ' input_ex01.xml'];
fprintf('Executing system command: $ %s\n', command);
system(command);

fprintf("\nGenerating matrix files ...\n");
tic
!./ROL_example_PDE-OPT_published_Helmholtz_KouriRidzalTuminaro2020_example_01.exe
toc

outfile = ['results/Results_' file]
fileID  = fopen(outfile, 'w');

% Run study.
for i=0:6
  beta  = 10^(-i)
  niter = run(beta, wave, fullObs, fullCtrl, solver)
  fprintf(fileID, '%e  %d\n', beta, niter);
end

fclose(fileID);
