F = {'F1','bmw3_2','bmw7','bmwcra','boneS10','inline_1'};
for i=1:length(F),
    DAT = load(strcat(F{i},'.serial.alllog'));
    CGBPlot1(F{i},DAT); close;
end

for i=1:length(F),
    DAT = load(strcat(F{i},'.tbb.alllog'));
    CGBPlot2(F{i},DAT,'tbb'); close;
end

for i=1:length(F),
    DAT = load(strcat(F{i},'.alllog'));
    CGBPlot2(F{i},DAT,'tpi'); close;
end