function [Rating_per_SE,Rating_per_RE]=Rating_EvaGame

version=1;
%This function is used to calculate The Rating Subjects provided for each
%condition [1= 'Extremely Unlikely',2='Slightly Unlikely',3='Slightly Likely',4='Extremely Likely']

%check the input arguments
% subject list
if nargin<2
    sublist_name=sprintf('sublist_ver%d.txt',version);
    sublist = textread(sublist_name,'%s','delimiter','\n');
    sub_exclude = [];
    sublist(sub_exclude)=[];
    
    
else
    if ~iscell(sublist)
        sublist = {sublist};
    end
end

nSubj = numel(sublist);

% directory
dir_fig = '../figures';
dir_rebuilt_data='../organized';

%index of self-efficacy and response efficacy
SE=[0,1,2];
RE=[0,1,2];
Rating_per_SE=zeros(size(SE,2),nSubj);
Rating_per_RE=zeros(size(RE,2),nSubj);

for s=1:nSubj
    dataFile=fullfile(dir_rebuilt_data,sprintf('%s_StructData.mat',sublist{s}));
    load(dataFile);   
    % Take the data from evaluation game
    for b=1:size(SE,2)
        idx_include_pair_SE= (StructData.Only_block_pair_SE==SE(b) &(StructData.Only_block_Rating~=-1)&(StructData.Only_block_RatingRT~=-1));
        Rating_per_SE(b,s) = mean(StructData.Only_block_Rating(idx_include_pair_SE));
    end
    for b=1:size(RE,2)
        idx_include_pair_RE= (StructData.Only_block_pair_RE==RE(b)&(StructData.Only_block_Rating~=-1)&(StructData.Only_block_RatingRT~=-1));
        Rating_per_RE(b,s) = mean(StructData.Only_block_Rating(idx_include_pair_RE));
    end
    
end




figure
fig_setting_default
Subj_M_bar=nanmean(Rating_per_SE,2);
Subj_M_bar_se=nanstd(Rating_per_SE,0,2)/sqrt(nSubj);
hold on
b=bar(Subj_M_bar,'grouped');
b.FaceColor=[0.8,0.8,0.8];
drawnow;

bar_xpos = NaN(size(Rating_per_SE));
nbars = size(Subj_M_bar, 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',Subj_M_bar,Subj_M_bar_se,'k','linestyle','none','LineWidth',2);
title('Self-efficacy','Interpreter', 'none');
title_name=sprintf('Rating_OnlySE_ver%d',version);
%title(title_name,'Interpreter', 'none');
set(gca, 'XTick', [1,2, 3], 'XTickLabel', {'Low','Mid' 'High'});
set(gca, 'YTick', [1,2,3,4], 'YTickLabel', {'Extremely Unlikely','Slightly Unlikely','Slightly Likely','Extremely Likely'});
%legend('Low RE','Mid RE','High RE','Location','best');
ylim([1,5]);
hold off
% output figure
outputfile = fullfile(dir_fig,title_name);
print(outputfile,'-dpng');

figure
fig_setting_default
Subj_M_bar=nanmean(Rating_per_RE,2);
Subj_M_bar_se=nanstd(Rating_per_RE,0,2)/sqrt(nSubj);
hold on
b=bar(Subj_M_bar,'grouped');
b.FaceColor=[0.8,0.8,0.8];
drawnow;

bar_xpos = NaN(size(Rating_per_RE));
nbars = size(Subj_M_bar, 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',Subj_M_bar,Subj_M_bar_se,'k','linestyle','none','LineWidth',2);
title_name=sprintf('Rating_OnlyRE_ver%d',version);
%title(title_name,'Interpreter', 'none');
title('Response-efficacy','Interpreter', 'none');
set(gca, 'XTick', [1,2, 3], 'XTickLabel', {'Low','Mid' 'High'});
set(gca, 'YTick', [1,2,3,4], 'YTickLabel', {'Extremely Unlikely','Slightly Unlikely','Slightly Likely','Extremely Likely'});
%legend('Low RE','Mid RE','High RE','Location','best');
ylim([1,5]);
hold off
% output figure
outputfile = fullfile(dir_fig,title_name);
print(outputfile,'-dpng');




%%%%%%%%%%%%%%%%%%
%%%%%%stats%%%%%%%
%%%%%%%%%%%%%%%%%%
%Rating_per_SE
dataTable = array2table(Rating_per_SE', 'VariableNames', {'Group1', 'Group2', 'Group3'});

% Create a grouping variable
group = table({'Group1'; 'Group2'; 'Group3'},'VariableNames',{'Group'});


% Perform one-way ANOVA using ranova
rm = fitrm(dataTable, 'Group1-Group3~1', 'WithinDesign', group);
ranovaResults = ranova(rm);

% Display the one-way ANOVA table using ranova
disp(ranovaResults);


%Rating_per_RE
dataTable = array2table(Rating_per_RE', 'VariableNames', {'Group1', 'Group2', 'Group3'});

% Create a grouping variable
group = table({'Group1'; 'Group2'; 'Group3'},'VariableNames',{'Group'});


% Perform one-way ANOVA using ranova
rm = fitrm(dataTable, 'Group1-Group3~1', 'WithinDesign', group);
ranovaResults = ranova(rm);

% Display the one-way ANOVA table using ranova
disp(ranovaResults);


