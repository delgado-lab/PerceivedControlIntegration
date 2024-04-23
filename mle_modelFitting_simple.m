function neg_sum_logmle = mle_modelFitting_simple(parameter_list, type, decisionProblems)

% run model fitting through maximum likelihood estimation
%
% YYY, 20200309

% choose the model
switch type
    case {'SE','RE'}
        b = parameter_list(1);
        s = parameter_list(2);
    case {'SERE'}
        b = parameter_list(1);
        s = parameter_list(2);
        r = parameter_list(3);
    case {'EP'}
        b = parameter_list(1);
        s = parameter_list(2);
        r = parameter_list(3);
        k = parameter_list(4);
        
end

SE=[0.2,0.5,0.8];
RE=[0.2,0.5,0.8];
%PersonalSE=decisionProblems.PersonalSE;


% initialize variables
idx_inclde= (decisionProblems.Choice_block_pair1_chosen~=-1);
realChoice = decisionProblems.Choice_block_pair1_chosen(idx_inclde);
decisionProblems.pair1SE=decisionProblems.Choice_block_pair1SE(idx_inclde);
decisionProblems.pair2SE=decisionProblems.Choice_block_pair2SE(idx_inclde);
decisionProblems.pair1RE=decisionProblems.Choice_block_pair1RE(idx_inclde);
decisionProblems.pair2RE=decisionProblems.Choice_block_pair2RE(idx_inclde);

nTrial = sum(decisionProblems.Choice_block_pair1_chosen~=-1); %Need to change later
pChoice_all = zeros(nTrial, 1);
logit_choice = @(x) 1/(1+exp(-1*x));

% run model
for t = 1:nTrial
    switch type
        case {'EP'}
            EP1=SE(decisionProblems.pair1SE(t)+1)*RE(decisionProblems.pair1RE(t)+1);
            EP2=SE(decisionProblems.pair2SE(t)+1)*RE(decisionProblems.pair2RE(t)+1);
            Option1=s*SE(decisionProblems.pair1SE(t)+1)+r*RE(decisionProblems.pair1RE(t)+1)+k*EP1;
            Option2=s*SE(decisionProblems.pair2SE(t)+1)+r*RE(decisionProblems.pair2RE(t)+1)+k*EP2;
            decisionValue = b*(Option1-Option2);
        case {'SE'}
            Option1=s*SE(decisionProblems.pair1SE(t)+1);
            Option2=s*SE(decisionProblems.pair2SE(t)+1);
            decisionValue = b*(Option1-Option2);
        case {'RE'}
            Option1=s*RE(decisionProblems.pair1RE(t)+1);
            Option2=s*RE(decisionProblems.pair2RE(t)+1);
            decisionValue = b*(Option1-Option2);
        case {'SERE'}
            Option1=s*SE(decisionProblems.pair1SE(t)+1)+r*RE(decisionProblems.pair1RE(t)+1);
            Option2=s*SE(decisionProblems.pair2SE(t)+1)+r*RE(decisionProblems.pair2RE(t)+1);
            decisionValue = b*(Option1-Option2);
        
    end
    
    pChoice_t = logit_choice(decisionValue);
    pChoice_t = min(max(pChoice_t,0.0000001),0.9999999); % constrain the probability between [0.0001, 0.9999]
    switch realChoice(t)
        case 1
            pChoice_all(t,1) = pChoice_t;
        case 0
            pChoice_all(t,1) = 1-pChoice_t;
    end
    
end

% negative sum of log maximum likelihood estimation
% because fmincon could only determine the minimum fitting value,
% we have to transform the sum of log likelihood estimation to negative.
neg_sum_logmle = -1*sum(log(pChoice_all));




