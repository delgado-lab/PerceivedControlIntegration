function allData = group_mle_condition_BIC(sublist, fig_on)

% group results of BIC from multiple subjects
%

version=1;
%check the input arguments
% subject list
if nargin<1
    sublist_name=sprintf('sublist_ver%d.txt',version);
    sublist = textread(sublist_name,'%s','delimiter','\n');
    sub_exclude = [39];
    sublist(sub_exclude)=[];
    
else
    if ~iscell(sublist)
        sublist = {sublist};
    end
end

nSubj = numel(sublist);
% whether to draw figures
if nargin<2
    fig_on = 1;
end

% which mle model
logitModel = {'SE','RE','SERE','EP'};
nModel = numel(logitModel);

% directory
dir_reg = 'mle';
dir_fig = 'figures';


% initialize variables
allData = cell(1, nModel);
for i = 1:nModel
    allData{i} = zeros(nSubj,1);
end

nTrials=102;

% start to group data
for i = 1:nModel
    nameModel = logitModel{i};
    for s = 1:nSubj
        subID = sublist{s};
        regFile = fullfile(dir_reg, sprintf('mle_ver%d_%s_%s.mat',version, logitModel{i}, sublist{s}));
        load(regFile);
        % BIC = -2*LLE+k*ln(nTrials)
        nParameter = numel(x_bestfit);
        BIC = -2*logmle_bestfit+nParameter*log(nTrials);
        allData{i}(s,1) = BIC; % In all trials condition
        
    end
end
subjData = zeros(nSubj,1);
for i = 1:nModel
    subjData(:,i) = mean(allData{i},2);
    %     subjData(:,i,:) = allData{i};
end

% bar graph
if fig_on
    
    barData = mean(subjData,1);
    stdData = std(subjData,1);
    seData = stdData./sqrt(nSubj);
    figure;
    fg = fig_setting_default();
    hold on
    nrow = size(barData,1);
    ncol = size(barData,2);
    h = bar(barData);
    h.FaceColor=[0.8,0.8,0.8];
e = errorbar([1:size(barData,2)], barData, seData,'LineWidth',2);
set(e, 'linestyle', 'none', 'linewidth', 2, 'color', 'k');
drawnow;
    
        set(gca, 'XTick', [1:size(barData,2)], 'XTickLabel', {'SE-only','RE-only','Integrated','ExpectedValue'});
    %legend(logitModel, 'interpreter', 'none', 'Location', 'NorthWest');
    ylabel('BIC');
    
    %ylim([0,80]);
    hold off
    
  
    
    % output figure
    print(fullfile(dir_fig,'BIC'), '-depsc');
    
end


% % statistics
model_comparison = [1 2;1 3;1 4;2 3;2 4;3 4];


for m = 1:size(model_comparison, 1)
    
    model_1 = model_comparison(m,1);
    model_2 = model_comparison(m,2);
    fprintf('%s vs %s\n', logitModel{model_1}, logitModel{model_2});
    data_1 = subjData(:,model_1);
    data_2 = subjData(:,model_2);
    [P,H]=signtest(data_1,data_2)

    
end



