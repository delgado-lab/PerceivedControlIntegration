function [allData] = group_mle_condition(sublist, fig_on)

version=1;
%check the input arguments
% subject list
if nargin<1
    sublist_name=sprintf('sublist_ver%d.txt',version);
    sublist = textread(sublist_name,'%s','delimiter','\n');
    sub_exclude = [39]; %This participant did not finish the study
    sublist(sub_exclude)=[];
    
else
    if ~iscell(sublist)
        sublist = {sublist};
    end
end

nSubj = numel(sublist);

if nargin<2
    fig_on = 1;
end

% directory
dir_fig = 'figures';
dir_reg = 'mle';
dir_rebuilt_data='organized';

% which mle model
logitModel ='SERE'; %Since the SERE model wins the model comparision, we are checking data from this model 
%condition
nameSE={'LowSE','MidSE','HighSE'};
nameRE={'LowRE','MidRE','HighRE'};

switch logitModel

    case {'SERE'}
        title_parameter = {'b', 's', 'r'};

        
end
nParameter = numel(title_parameter);

% initialize variables
allData = zeros(nSubj,nParameter);

% start to group data
for s = 1:nSubj
    
    subID = sublist{s};
    regFile = fullfile(dir_reg, sprintf('mle_ver%d_%s_%s.mat',version, logitModel, sublist{s}));
    
    load(regFile);
    for i = 1:nParameter
        allData(s,:) = x_bestfit;
    end
    
end
switch logitModel

    case {'SERE'}
        SR_Data=allData(:,[2,3]);

end



figure
fig_setting_default
Subj_M_bar=mean(SR_Data,1);
Subj_M_bar_se=std(SR_Data,0,1)/sqrt(nSubj);
hold on
b=bar(Subj_M_bar,'grouped','w');
b.FaceColor=[0.8,0.8,0.8];
e = errorbar([1:size(SR_Data,2)], Subj_M_bar, Subj_M_bar_se,'LineWidth',2);
set(e, 'linestyle', 'none', 'linewidth', 2, 'color', 'k');
drawnow;

hold off

title_name=sprintf('Estimated_parameter_%s_ver%d',logitModel,version);
%title(title_name,'Interpreter', 'none');
switch logitModel
    case {'SERE'}
        set(gca, 'XTick', [1,2], 'XTickLabel', {'SE','RE'});
end

ylabel('Beta')
ylim([0,14]);
hold off
% output figure
outputfile = fullfile(dir_fig,title_name);
print(outputfile,'-dpng');

[h1,p1]=signrank(allData(:,2),allData(:,3)); 
[h1,p1]=signrank(allData(:,2)); %Check whether Wse is significant different than 0
[h1,p1]=signrank(allData(:,3)); %Check whether Wre is significant different than 0

