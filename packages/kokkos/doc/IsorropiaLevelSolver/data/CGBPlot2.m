function CGBPlot2(ALLDATA);
  ts = size(ALLDATA,1) / 4;
  TARGET = ALLDATA(0*ts+1:1*ts,:);
  NOLNCH = ALLDATA(1*ts+1:2*ts,:);
  NOPFOR = ALLDATA(2*ts+1:3*ts,:);
  NOCOMP = ALLDATA(3*ts+1:4*ts,:);
  
  % want three bits of data:
  % epetra time
  % kokkos::crsmatrix time vs. numthread
  % levelsolver::applyinverse time vs. numthread vs. capability
  % the first two come from target, the last one from all four
  LEVSOL_LBL = 0;
  KOKKOS_LBL = 1;
  EPETRA_LBL = 2;

  % individual data set (TARGET,NOLNCH,NOPFOR,NOCOMP) looks like
  % testlabel numthreads time
  EpetraTime = TARGET( min(find(TARGET(:,1) == EPETRA_LBL)), 3);
  fprintf('Epetra time is %d\n',EpetraTime);
  KokkosTime = TARGET( find(TARGET(:,1) == KOKKOS_LBL), 2:end);
  Threads = KokkosTime(:,1);
  numThreads = size(KokkosTime,1);
  fprintf('Number of thread tests is %d: ',numThreads);
  for i=1:numThreads, 
      fprintf('%d ',KokkosTime(i,1)); 
  end
  fprintf('\n');

  % get LevelSolver times
  LS_TARGET = TARGET( find(TARGET(:,1) == LEVSOL_LBL), 2:end);
  LS_NOLNCH = NOLNCH( find(NOLNCH(:,1) == LEVSOL_LBL), 2:end);
  LS_NOPFOR = NOPFOR( find(NOPFOR(:,1) == LEVSOL_LBL), 2:end);
  LS_NOCOMP = NOCOMP( find(NOCOMP(:,1) == LEVSOL_LBL), 2:end);

  % plot
  lbls = {};
  fig = figure;
  loglog(KokkosTime(:,1),KokkosTime(:,2),'-+r'); 
  lbls = {lbls{:} 'Kokkos'};
  hold on;
  loglog(LS_TARGET(:,1),LS_TARGET(:,2),'-ob');
  lbls = {lbls{:} 'LevSol (target)'};
  loglog(LS_NOCOMP(:,1),LS_NOCOMP(:,2),'-.db');
  lbls = {lbls{:} 'LevSol (no comp)'};
  loglog(LS_NOPFOR(:,1),LS_NOPFOR(:,2),'--sb');
  lbls = {lbls{:} 'LevSol (no pfor)'};
  loglog(LS_NOLNCH(:,1),LS_NOLNCH(:,2),':*b');
  lbls = {lbls{:} 'LevSol (no launch)'};
  eline = line(get(gca,'XLim'),[EpetraTime EpetraTime]);
  lbls = {lbls{:} 'Epetra'};
  set(eline,'Color',[0 0 0]);
  set(gca,'XTickLabel',Threads);
  set(gca,'XTick',Threads);
  set(gca,'XLim',[min(Threads)-eps,max(Threads)+eps]);
  set(gca,'YLim',[1e-8 1]);
  legend(lbls{:},'Location','SouthEast');
  xlabel('Number of Threads');
  ylabel('Average runtime (s)');
  grid off;
  fn = inputname(1);
  title(strcat(fn, ' (TPI)'));
  saveas(fig,strcat(fn,'_tpi.png'),'png');