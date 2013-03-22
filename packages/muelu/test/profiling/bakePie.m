% disables interpreting underscore as LaTeX
set(0,'defaulttextinterpreter','none')
% load data, N-by-2 cell array, with data in 1st col, labels in 2nd col
timings
[T,P] = sort(cell2mat(data(:,2)),'descend');
labels = data(:,1);
% group all timings less than "thresh" into a single entry
thresh = 0.075;
maxSlices=13;
[newT,newlabels] = squash(T,labels(P),thresh,maxSlices);
% plot results
pie(newT,newlabels);
