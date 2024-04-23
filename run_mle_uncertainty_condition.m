function run_mle_uncertainty_condition(sublist, logitModel)

% run model fitting through maximum likelihood estimation
% implement on data under each condition
%

%Four models could be fit to using the following code.
% 'SE' model -only using SE information, and completely ignore RE
% 'RE' model -only using RE information, and completely ignore SE
% 'SERE' model- using both SE and RE information to make decision
% 'EP' model - In addition to SERE inforamtion, adding expected value in the model 

% YYY, 20200309
version=1;
% check the input arguments
% subject list
if nargin<1
    sublist_name=sprintf('sublist_ver%d.txt',version);
    sublist = textread(sublist_name,'%s','delimiter','\n');
    sub_exclude = [39]; % This participant did not finish the study
    sublist(sub_exclude)=[];
else
    if ~iscell(sublist)
        sublist = {sublist};
    end
end

nSubj = numel(sublist);



% which mle model 
logitModel ='SERE'; %Change the model name while switching models

% directory
dir_rebuilt_data='organized';
dir_reg = 'mle';
mkdir(dir_reg);

% blocks

% start to model fitting
for s = 1:nSubj
    
    dataFile=fullfile(dir_rebuilt_data,sprintf('%s_StructData.mat',sublist{s}));
    load(dataFile);
    data = StructData;

    % model fitting
    switch logitModel
        case {'SE','RE'}
            nameParameter = {'b','s'}; % name of paramenters
            lb = [0.0001,  0.0001]; % low bound
            ub = [    50,   20]; % high bound
            cutinterval = [5,  5]; % grid
            nParameter = numel(cutinterval);
            % multiple start seeds
            x0 = zeros(prod(cutinterval), nParameter);
            idx = 0;
            for i = 1:cutinterval(1)
                for j = 1:cutinterval(2)
                    idx = idx + 1;
                    x0(idx, :) = [...
                        lb(1)+(i-1)*(ub(1)-lb(1))/cutinterval(1);...
                        lb(2)+(j-1)*(ub(2)-lb(2))/cutinterval(2)
                        ];
                end
            end
            
            mleModel = @(x)mle_modelFitting_simple(x, logitModel, data);
            
        case {'SERE'}
            
            nameParameter = {'b','s','r'}; % name of paramenters
            lb = [0.0001,  0.0001, 0.0001]; % low bound
            ub = [    50,   20, 20]; % high bound
            cutinterval = [5,  5, 5]; % grid
            nParameter = numel(cutinterval);
            % multiple start seeds
            x0 = zeros(prod(cutinterval), nParameter);
            idx = 0;
            for i = 1:cutinterval(1)
                for j = 1:cutinterval(2)
                    for k = 1:cutinterval(3)
                        idx = idx + 1;
                        x0(idx, :) = [...
                            lb(1)+(i-1)*(ub(1)-lb(1))/cutinterval(1);...
                            lb(2)+(j-1)*(ub(2)-lb(2))/cutinterval(2);...
                            lb(3)+(k-1)*(ub(3)-lb(3))/cutinterval(3)
                            ];
                    end
                end
            end
            
            mleModel = @(x)mle_modelFitting_simple(x, logitModel, data);
            
        case {'EP'}
            
            nameParameter = {'b','s','r','k'}; % name of paramenters
            lb = [0.0001,  0.0001, 0.0001, 0.0001]; % low bound
            ub = [    50,   20, 20, 20]; % high bound
            cutinterval = [5,  5, 5, 5]; % grid
            nParameter = numel(cutinterval);
            % multiple start seeds
            x0 = zeros(prod(cutinterval), nParameter);
            idx = 0;
            for i = 1:cutinterval(1)
                for j = 1:cutinterval(2)
                    for k = 1:cutinterval(3)
                        for l= 1:cutinterval(4)
                            idx = idx + 1;
                            x0(idx, :) = [...
                                lb(1)+(i-1)*(ub(1)-lb(1))/cutinterval(1);...
                                lb(2)+(j-1)*(ub(2)-lb(2))/cutinterval(2);...
                                lb(3)+(k-1)*(ub(3)-lb(3))/cutinterval(3);...
                                lb(4)+(l-1)*(ub(4)-lb(4))/cutinterval(4)
                                ];
                        end
                    end
                end
            end
            
            mleModel = @(x)mle_modelFitting_simple(x, logitModel, data);

            
    end % end of model fitting
    
    
    
    
    % run multiple fitting at different start seeds
    sum_logmle_all = zeros(idx, 1);
    x_all = zeros(idx, nParameter);
    options = optimoptions('fmincon',...
        'Display','off',...
        'Algorithm','interior-point');
    tic
    for i = 1:idx
        [x, neg_sum_logmle] = fmincon(@(x)mleModel(x), x0(i,:), [], [], [], [], lb, ub,[],options);
        sum_logmle = -1*neg_sum_logmle;
        sum_logmle_all(i,1) = sum_logmle;
        x_all(i,:) = x;
        fprintf('%d\t', i);
        fprintf('%.3f\t', sum_logmle);
        fprintf('%.3f\t', x);
        fprintf('\n');
    end
    toc
    
    % choose the best fitting value and the corresponding parament set
    [logmle_bestfit, idx_bestfit] = max(sum_logmle_all);
    x_bestfit = x_all(idx_bestfit, :);
    
    fittingFile = fullfile(dir_reg, sprintf('mle_ver%d_%s_%s.mat',version, logitModel, sublist{s}));
    save(fittingFile,'nameParameter','x_bestfit', 'logmle_bestfit');
    
    
    
end % end fo subject


