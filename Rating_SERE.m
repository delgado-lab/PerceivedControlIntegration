function [Rating_per_con_SE,Rating_per_con_RE]=Rating_SERE(version)

%This function is used to calculate The Rating Subjects provided for each
%condition [1='Extremely Unconfident',2='Slightly Unconfident',3='Slightly Confident',4='Extremely Confident']
version=1;
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
dir_fig = 'figures';
dir_rebuilt_data='organized';

%index of self-efficacy and response efficacy
SE=[0,1,2];
RE=[0,1,2];
Rating_per_con_SE=zeros(size(SE,2),size(SE,2),nSubj);
Rating_per_con_RE=zeros(size(RE,2),size(RE,2),nSubj);

for s=1:nSubj
    %Load data for each participant
    dataFile=fullfile(dir_rebuilt_data,sprintf('%s_StructData.mat',sublist{s}));
    load(dataFile);
    
    %Put two blocks in together
    StructData.All_block_pair_SE= [StructData.S1_block_pair_SE;StructData.S2_block_pair_SE];
    StructData.All_block_pair_RE= [StructData.S1_block_pair_RE;StructData.S2_block_pair_RE];
    StructData.All_block_typeQ= [StructData.S1_typeQ;StructData.S2_typeQ];
    StructData.All_block_Question= [StructData.S1_Question;StructData.S2_Question];
    StructData.All_block_Rating= [StructData.S1_block_Rating;StructData.S2_block_Rating];
    StructData.All_block_RatingRT= [StructData.S1_block_RatingRT;StructData.S2_block_RatingRT];
    
    %Two type of questions, but only take the one focused on SE and have
    %respond
    for b=1:size(SE,2)
        idx_include_pair_SE= (StructData.All_block_pair_SE==SE(b) &(StructData.All_block_typeQ==0) &(StructData.All_block_Rating~=0));
        Rating_per_SE(b,s) = mean(StructData.All_block_Rating(idx_include_pair_SE));
        for p=1:size(RE,2)
            idx_include_pair=(StructData.All_block_pair_SE==SE(b) ...
                & StructData.All_block_pair_RE==RE(p) ...
                &(StructData.All_block_typeQ==0)...
                & StructData.All_block_RatingRT~=0);
            Rating_per_con_SE(b,p,s)=mean(StructData.All_block_Rating(idx_include_pair));
            
            
        end
    end
    %Two type of questions, but only take the one focused on RE and have
    %respond
    for p=1:size(RE,2)
        idx_include_pair_RE= (StructData.All_block_pair_RE==RE(p) &(StructData.All_block_typeQ==1) &(StructData.All_block_Rating~=0));
        Rating_per_RE(p,s) = mean(StructData.All_block_Rating(idx_include_pair_RE));
        for b=1:size(SE,2)
            idx_include_pair=(StructData.All_block_pair_SE==SE(b) ...
                & StructData.All_block_pair_RE==RE(p) ...
                &(StructData.All_block_typeQ==1)...
                & StructData.All_block_Rating~=0);
            Rating_per_con_RE(b,p,s)=mean(StructData.All_block_Rating(idx_include_pair));
            
            
        end
    end
end



figure
fig_setting_default
Subj_M_bar=nanmean(Rating_per_SE,2);
Subj_M_bar_se=nanstd(Rating_per_SE,0,2)/sqrt(nSubj);
hold on
b=bar(Subj_M_bar,'grouped','FaceColor','none');
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
title_name=sprintf('Rating_SE_ver%d',version);
title('Self-efficacy','Interpreter', 'none')
set(gca, 'XTick', [1,2, 3], 'XTickLabel', {'Low','Mid' 'High'});

set(gca, 'YTick', [1,2,3,4], 'YTickLabel', {'Extremely Unconfident','Slightly Unconfident','Slightly Confident','Extremely Confident'});
ylim([1,5]);
hold off
% output figure
outputfile = fullfile(dir_fig,title_name);
print(outputfile,'-dpng');



%%statistic test
%%%%% comparison %%%%%
fprintf('# comparison\n');
nVar = size(SE,2)*size(RE,2);
factorNames = {'SE', 'RE'};
model_within = 'SE*RE';
list_SE={'Low_SE','Mid_SE','High_SE'};
list_RE={'Low_RE','Mid_RE','High_RE'};

groupData = NaN(nSubj, nVar);
list_var = cell(nVar,1);
var_SE = cell(nVar,1);
var_RE = cell(nVar,1);
idx_var = 0;
for v = 1:size(SE,2)
    for k = 1:size(RE,2)
        idx_var = idx_var + 1;
        groupData(:,idx_var) = Rating_per_con_SE(v,k,:);
        list_var{idx_var} = sprintf('V%d', idx_var);
        var_SE{idx_var} = list_SE{v};
        var_RE{idx_var} = list_RE{k};
    end
end
table_data = array2table(groupData, 'VariableNames', list_var);
table_within = table(var_SE, var_RE, 'VariableNames', factorNames);

% ANOVA
rm = fitrm(table_data, sprintf('V1-V%d~1', nVar), 'WithinDesign', table_within);
anova_table = ranova(rm, 'WithinModel', model_within)

%mc_table = multcompare(rm, 'SE', 'ComparisonType', 'lsd')
%mc_table = multcompare(rm, 'SE', 'by', 'RE', 'ComparisonType', 'lsd')
%mc_table = multcompare(rm, 'RE', 'ComparisonType', 'lsd')
%mc_table = multcompare(rm, 'RE', 'by', 'SE', 'ComparisonType', 'lsd')

fprintf('\n');



%%%%%%RE FIGURE%%%%%%
figure
fig_setting_default
Subj_M_bar=nanmean(Rating_per_RE,2);
Subj_M_bar_se=nanstd(Rating_per_RE,0,2)/sqrt(nSubj);
hold on
b=bar(Subj_M_bar,'grouped','w');
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
title_name=sprintf('Rating_RE_ver%d',version);
title('Response-efficacy','Interpreter', 'none')
set(gca, 'XTick', [1,2, 3], 'XTickLabel', {'Low','Mid' 'High'});
set(gca, 'YTick', [1,2,3,4], 'YTickLabel', {'Extremely Unconfident','Slightly Unconfident','Slightly Confident','Extremely Confident'});
ylim([1,5]);
hold off
% output figure
outputfile = fullfile(dir_fig,title_name);
print(outputfile,'-dpng');







%%statistic test
%%%%% comparison %%%%%
fprintf('# comparison\n');
nVar = size(SE,2)*size(RE,2);
factorNames = {'SE', 'RE'};
model_within = 'SE*RE';
list_SE={'Low_SE','Mid_SE','High_SE'};
list_RE={'Low_RE','Mid_RE','High_RE'};

groupData = NaN(nSubj, nVar);
list_var = cell(nVar,1);
var_SE = cell(nVar,1);
var_RE = cell(nVar,1);
idx_var = 0;
for v = 1:size(SE,2)
    for k = 1:size(RE,2)
        idx_var = idx_var + 1;
        groupData(:,idx_var) = Rating_per_con_RE(v,k,:);
        list_var{idx_var} = sprintf('V%d', idx_var);
        var_SE{idx_var} = list_SE{v};
        var_RE{idx_var} = list_RE{k};
    end
end
table_data = array2table(groupData, 'VariableNames', list_var);
table_within = table(var_SE, var_RE, 'VariableNames', factorNames);

% ANOVA
rm = fitrm(table_data, sprintf('V1-V%d~1', nVar), 'WithinDesign', table_within);
anova_table = ranova(rm, 'WithinModel', model_within)

%mc_table = multcompare(rm, 'SE', 'ComparisonType', 'lsd')
%mc_table = multcompare(rm, 'SE', 'by', 'RE', 'ComparisonType', 'lsd')
%mc_table = multcompare(rm, 'RE', 'ComparisonType', 'lsd')
%mc_table = multcompare(rm, 'RE', 'by', 'SE', 'ComparisonType', 'lsd')

fprintf('\n');





