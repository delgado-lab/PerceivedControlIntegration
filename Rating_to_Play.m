function Rating_to_Play

%This function is used to calculate The Rating Subjects provided for each
%condition [1= Definitely Not, 2= Not, 3=Yes, 4=Definitely Yes]
version=1;
%check the input arguments
% subject list
if nargin<2
    sublist_name=sprintf('sublist_ver%d.txt',version);
    sublist = textread(sublist_name,'%s','delimiter','\n');
    sub_exclude = [39]; %This participant did not finish the Play Game
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
Rating_per_con=zeros(size(SE,2),size(SE,2),nSubj);

for s=1:nSubj
    dataFile=fullfile(dir_rebuilt_data,sprintf('%s_StructData.mat',sublist{s}));
    load(dataFile);
    
    %organized data based on each condition
    for b=1:size(SE,2)
        for p=1:size(RE,2)
            idx_include_pair=(StructData.Play_block_pair_SE==SE(b)&StructData.Play_block_pair_RE==RE(p) & (StructData.Play_block_Play~=0));
            Rating_per_con(b,p,s)=mean(StructData.Play_block_Play(idx_include_pair));
           
        end
    end

end

figure
fig_setting_default
Subj_M_bar=mean(Rating_per_con,3);
Subj_M_bar_se=std(Rating_per_con,0,3)/sqrt(nSubj);
hold on
b=bar(Subj_M_bar,'grouped');
b(1).FaceColor=[0.8,0.8,0.8];
b(2).FaceColor=[0.7,0.7,0.7];
b(3).FaceColor=[0.5,0.5,0.5];
drawnow;

bar_xpos = NaN(size(Rating_per_con));
nbars = size(Subj_M_bar, 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',Subj_M_bar,Subj_M_bar_se,'k','linestyle','none');
title_name=sprintf('Rating_play_ver%d',version);
set(gca, 'XTick', [1,2, 3], 'XTickLabel', {'Low SE','Mid SE' 'High SE'});
set(gca, 'YTick', [1,2,3,4], 'YTickLabel', {'Extremely No','Slightly No','Slightly Yes','Extremely Yes'});

legend('Low RE','Mid RE','High RE','Location','best');
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
        groupData(:,idx_var) = Rating_per_con(v,k,:);
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





