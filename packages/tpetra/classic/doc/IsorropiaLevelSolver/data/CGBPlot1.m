function CGBPlot2(fn,ALLDATA);
  ts = size(ALLDATA,1) / 4;
  TARGET = ALLDATA(1:3,:);
  NOLNCH = ALLDATA(4:6,:);
  NOPFOR = ALLDATA(7:9,:);
  NOCOMP = ALLDATA(10:12,:);
  
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
  EpetraTime = TARGET( find(TARGET(:,1) == EPETRA_LBL), end);
  fprintf('Epetra time is %d\n',EpetraTime);
  KokkosTime = TARGET( find(TARGET(:,1) == KOKKOS_LBL), end);
  fprintf('Kokkos time is %d\n',KokkosTime);

  % get LevelSolver times
  LS_TARGET = TARGET( find(TARGET(:,1) == LEVSOL_LBL), end);
  LS_NOLNCH = NOLNCH( find(NOLNCH(:,1) == LEVSOL_LBL), end);
  LS_NOPFOR = NOPFOR( find(NOPFOR(:,1) == LEVSOL_LBL), end);
  LS_NOCOMP = NOCOMP( find(NOCOMP(:,1) == LEVSOL_LBL), end);

  % plot
  lbls = {};
  fig = figure;
  res = semilogy([1 2 3 4],[LS_TARGET LS_NOCOMP LS_NOPFOR LS_NOLNCH],'-ob');
  set(res,'LineWidth',1.0);
  hold on;
  xlbls = {'' 'Target' 'No Comp' 'No PFor' 'No Launch' ''};
  lbls = {lbls{:} 'LevSol'};
  kline = line([0 5],[KokkosTime KokkosTime]);
  lbls = {lbls{:} 'Kokkos'};
  set(kline,'Color',[1 0 0]);
  set(kline,'LineWidth',1.0);
  eline = line([0 5],[EpetraTime EpetraTime]);
  lbls = {lbls{:} 'Epetra'};
  set(eline,'Color',[0 0 0]);
  set(eline,'LineStyle','--');
  set(eline,'LineWidth',1.0);
  set(gca,'XTickLabel',xlbls);
  set(gca,'XTick',[0:5])
  set(gca,'XLim',[0 5]);
  set(gca,'YLim',[1e-8 1]);
  set(gca,'FontSize',13);
  set(gca,'LineWidth',1.0);
  l = legend(lbls{:},'Location','SouthWest');
  set(l,'FontSize',13);
  %ylabel('Average runtime (s)','FontSize',14);
  grid off;
  title(deTeX(strcat(fn, ' (Serial)')),'FontSize',16);
  saveas(fig,strcat(fn,'_serial.png'),'png');