function fMRI_organized_onset_new(subID)

% organize input and data for overall fMRI onset files
% the output onset files could be used to produce files for GLM analysis
% under different GLM models.
%
% 2022,03,02
version=1;
% subject list
if nargin<2
    sublist_name=sprintf('sublist_ver%d.txt',version);
    sublist = textread(sublist_name,'%s','delimiter','\n');
    
    %sub_exclude = [];
    %sublist(sub_exclude)=[];
else
    if ~iscell(sublist)
        sublist = {sublist};
    end
end
nSubj = numel(sublist);


% directory
dirInput = 'input';
dirData = 'organized';
dir_data= 'data';
dirParam = 'mle';
dirOnset = 'fMRI_onset';
mkdir(dirOnset);


% setting
list_SE=[1,2,3];
list_RE=[1,2,3];
Block_type={'EvaGame','SERE1','SERE2','Play'}; %Sub049 doesn't do Play block in MRI ,'Play'
nSE = numel(list_SE);
nRE = numel(list_RE);
nBlock=numel(Block_type);

for sub=1:nSubj
    
    % load files that is organized
    subID=sublist{sub};
    fprintf(subID)
    
    dataFile = fullfile(dirData, sprintf('%s_StructData.mat', subID));
    load(dataFile);
    data=StructData;
    
    % which types of data will be organized
    onsetFile = fullfile(dirOnset, sprintf('%s_fmri_onset.mat', subID));
    
    %load the original data (Time data)
    %Evaluation Game info
    dataFile=fullfile(dir_data,sprintf('%s_2_1_MRI_Rating_SE_RE.csv',subID));
    Timedata_Eva = readtable(dataFile, 'delimiter', ',');
    Ready_text_line=find(~isnan(Timedata_Eva.Ready_text_started));
    Eva_started=max(Timedata_Eva.center_OnsetTime_started(Ready_text_line),Timedata_Eva.PairONSET_card_prob_2_started(Ready_text_line));
    
    %SERE info
    if strcmp(sublist{sub},'Sub044') %Sub006;Sub049
       %read RESE1 first
       dataFile=fullfile(dir_data,sprintf('%s_2_2_Rating_MRI_SEREcombine_Eyetracking_session%d.csv',sublist{sub},1));
       Timedata = readtable(dataFile, 'delimiter', ',');
       Ready_text_line=find(~isnan(Timedata.Ready_text_started));
       idx_CombineB1_started=Ready_text_line;
    else
    dataFile=fullfile(dir_data,sprintf('%s_2_2_Rating_MRI_SEREcombine_Eyetracking.csv',subID));
    Timedata = readtable(dataFile, 'delimiter', ',');
    Ready_text_line=find(~isnan(Timedata.Ready_text_started));
    idx_CombineB1_started=Ready_text_line(1);
    idx_CombineB2_started=Ready_text_line(2);
    end
    
    %Play info
     dataFile=fullfile(dir_data,sprintf('%s_2_3_MRI_Rating_Play_Eyetracking.csv',subID));
     Timedata_Play = readtable(dataFile, 'delimiter', ',');
     idx_Play_started=find(~isnan(Timedata_Play.Ready_text_started));
     
    
    
    %Time started for the block
    CombineB1_started=max(Timedata.center_OnsetTime_2_started(idx_CombineB1_started),Timedata.center_OnsetTime_started(idx_CombineB1_started));
    if strcmp(sublist{sub},'Sub044') %Sub006;Sub049
    else 
    CombineB2_started=max(Timedata.center_OnsetTime_2_started(idx_CombineB2_started),Timedata.center_OnsetTime_started(idx_CombineB2_started));
    end
    Play_Block_started=Timedata_Play.center_OnsetTime_started(idx_Play_started);
    
    
    % time
    t_stim_present = 2;
    
    
    
    for b=1:nBlock
        
        switch Block_type{b}
            case'EvaGame'
                nTrial=90;
                onsetTitle = {...
                    'block_info', 'trial','typeQ','SE_info','RE_info','Stimulus_onset','Stimulus_duration',...
                    'ISI_onset','ISI_duration','Question','Rating','Rating_onset','Rating_duration',...
                    'ITI_onset','ITI_duration'
                    };
                nTitle = numel(onsetTitle);
                fmri_onset.onsetTitle{b} = onsetTitle;
                
                % start to organize data
                %allOnset = zeros(nBlock*nTrial, nTitle);
                allOnset{b} = zeros(1*nTrial, nTitle);
                current_trial = 0;
                
                %Calculated Time data
                %Stimulus_onset
                idx_included=find(~isnan(max(Timedata_Eva.center_OnsetTime_started,Timedata_Eva.PairONSET_card_prob_2_started)));
                n = find(strcmp(onsetTitle, 'Stimulus_onset'));
                StimOnsetTime=max(Timedata_Eva.center_OnsetTime_started(idx_included),Timedata_Eva.PairONSET_card_prob_2_started(idx_included));
                allOnset{b}(:, n) = StimOnsetTime-Eva_started;
                
                
                %ISI_onset
                n = find(strcmp(onsetTitle, 'ISI_onset'));
                switch StructData.Run_hand
                    case 0
                        ISIOnsetTime = max(Timedata_Eva.fixOnset_text_ISI_SE_2_started(idx_included),Timedata_Eva.fixOnset_text_ISI_RE_2_started(idx_included));
                    case 1
                         ISIOnsetTime = max(Timedata_Eva.fixOnset_text_ISI_SE_started(idx_included),Timedata_Eva.fixOnset_text_ISI_RE_started(idx_included));
                end
                allOnset{b}(:, n) =ISIOnsetTime-Eva_started;
                
                %Rating_onset  ->Not all trials have question
                n = find(strcmp(onsetTitle, 'Rating_onset'));
                switch StructData.Run_hand
                    case 0
                        QuestionOnsetTime=max(Timedata_Eva.key_resp_rating_SE_L_started(idx_included),Timedata_Eva.key_resp_rating_RE_L_started(idx_included));
                        
                    case 1
                        QuestionOnsetTime=max(Timedata_Eva.key_resp_rating_SE_R_started(idx_included),Timedata_Eva.key_resp_rating_RE_R_started(idx_included));
                end
                RatingOnset=QuestionOnsetTime-Eva_started;
                RatingOnset(isnan(RatingOnset))=0;
                allOnset{b}(:, n) = RatingOnset;
                
                %ITI_onset
                n = find(strcmp(onsetTitle, 'ITI_onset'));
                ITIOnset=Timedata_Eva.fixOnset_text_ITI_started(idx_included)-Eva_started;
                ITIOnset(isnan(ITIOnset))=0;
                allOnset{b}(:, n) = ITIOnset;
                
                
                %Stimulus duration
                n_Stim = find(strcmp(onsetTitle, 'Stimulus_onset'));
                n_ISI = find(strcmp(onsetTitle, 'ISI_onset'));
                n = find(strcmp(onsetTitle, 'Stimulus_duration'));
                allOnset{b}(:, n) = allOnset{b}(:, n_ISI)-allOnset{b}(:, n_Stim);
                
                %Rating_duration
                n = find(strcmp(onsetTitle, 'Rating_duration'));
                n_ITI=find(strcmp(onsetTitle, 'ITI_onset'));
                n_Rating = find(strcmp(onsetTitle, 'Rating_onset'));
                allOnset{b}(:, n) = allOnset{b}(:, n_ITI)-allOnset{b}(:, n_Rating);
                
                
                
                for t = 1:nTrial
                    
                    current_trial = current_trial + 1;
                    
                    % Block
                    n = find(strcmp(onsetTitle, 'block_info'));
                    allOnset{b}(current_trial, n) = b;
                    
                    % Trial
                    n = find(strcmp(onsetTitle, 'trial'));
                    allOnset{b}(current_trial, n) = t;
                    
                    % Type_Q
                    n = find(strcmp(onsetTitle, 'typeQ'));
                    allOnset{b}(current_trial, n) = data.Only_block_Q(t);%SE=0;RE=1
                    
                    % SE Info
                    n = find(strcmp(onsetTitle, 'SE_info'));
                    allOnset{b}(current_trial, n) = data.Only_block_pair_SE(t)+1;
                    
                    % RE Info
                    n = find(strcmp(onsetTitle, 'RE_info'));
                    allOnset{b}(current_trial, n) = data.Only_block_pair_RE(t)+1;
                    
                    %ISI duration ->Some trials don't have questions
                    n = find(strcmp(onsetTitle, 'ISI_duration'));
                    if data.Only_block_Q(t)==0 %No question trial
                        if current_trial==nTrial
                            current_ISI=data.Only_block_ISI(current_trial);
                        else
                            current_ISI=allOnset{b}(current_trial+1, n_Stim)-allOnset{b}(current_trial, n_ISI);
                        end
                    else
                        current_ISI=allOnset{b}(current_trial, n_Rating)-allOnset{b}(current_trial, n_ISI);
                    end
                    allOnset{b}(current_trial, n) =current_ISI;
                    
                    %Question
                    n = find(strcmp(onsetTitle, 'Question'));
                    allOnset{b}(current_trial, n)=data.Only_block_Q(t);
                    
                    %Rating
                    n = find(strcmp(onsetTitle, 'Rating'));
                    allOnset{b}(current_trial, n) = data.Only_block_Rating(t);
                    
                    %ITI duration
                    n = find(strcmp(onsetTitle, 'ITI_duration'));
                    if data.Only_block_Q(t)==0 %No question trial
                        current_ITI=0;
                    else
                        if current_trial==nTrial
                            current_ITI=data.Only_block_ITI(current_trial);
                        else
                            current_ITI=allOnset{b}(current_trial+1, n_Stim)-allOnset{b}(current_trial, n_ITI);
                        end
                    end
                    allOnset{b}(current_trial, n) =current_ITI;
                    
                end
                
            case 'SERE1'
                nTrial=72;
                
                onsetTitle = {...
                    'block_info', 'trial','typeQ','SE_info','RE_info','Stimulus_onset','Stimulus_duration',...
                    'ISI_onset','ISI_duration','Question','Rating','Rating_onset','Rating_duration',...
                    'ITI_onset','ITI_duration'
                    };
                nTitle = numel(onsetTitle);
                fmri_onset.onsetTitle{b} = onsetTitle;
                
                % start to organize data
                %allOnset = zeros(nBlock*nTrial, nTitle);
                allOnset{b} = zeros(1*nTrial, nTitle);
                current_trial = 0;
                
                %Calculated Time data
                %Stimulus_onset   
                if strcmp(sublist{sub},'Sub044') %Sub006;Sub049
                    idx_included=find(~isnan(max(Timedata.center_OnsetTime_started,Timedata.center_OnsetTime_2_started)));
                else 
                idx_included=find(~isnan(max(Timedata.center_OnsetTime_started,Timedata.center_OnsetTime_2_started)));
                idx_included=idx_included(idx_included<idx_CombineB2_started);
                end
                n = find(strcmp(onsetTitle, 'Stimulus_onset'));
                StimOnsetTime=max(Timedata.center_OnsetTime_started(idx_included),Timedata.center_OnsetTime_2_started(idx_included));
                allOnset{b}(:, n) = StimOnsetTime-CombineB1_started;
                
                
                %ISI_onset
                n = find(strcmp(onsetTitle, 'ISI_onset'));
                switch StructData.Run_hand
                    case 0
                        ISIOnsetTime = max(Timedata.fixOnset_text_ISI_SE_2_started(idx_included),Timedata.fixOnset_text_ISI_RE_2_started(idx_included));
                    case 1
                        ISIOnsetTime = max(Timedata.fixOnset_text_ISI_SE_started(idx_included),Timedata.fixOnset_text_ISI_RE_started(idx_included));
                end
                allOnset{b}(:, n) =ISIOnsetTime-CombineB1_started;
                
                %Rating_onset  ->Not all trials have question
                n = find(strcmp(onsetTitle, 'Rating_onset'));
                switch StructData.Run_hand
                    case 0
                        QuestionOnsetTime=max(Timedata.key_resp_rating_SE_L_started(idx_included),Timedata.key_resp_rating_RE_L_started(idx_included));
                        
                    case 1
                        QuestionOnsetTime=max(Timedata.key_resp_rating_SE_R_started(idx_included),Timedata.key_resp_rating_RE_R_started(idx_included));
                end
                RatingOnset=QuestionOnsetTime-CombineB1_started;
                RatingOnset(isnan(RatingOnset))=0;
                allOnset{b}(:, n) = RatingOnset;
                
                %ITI_onset
                n = find(strcmp(onsetTitle, 'ITI_onset'));
                ITIOnset=Timedata.fixOnset_text_ITI_started(idx_included)-CombineB1_started;
                ITIOnset(isnan(ITIOnset))=0;
                allOnset{b}(:, n) = ITIOnset;
                
                
                %Stimulus duration
                n_Stim = find(strcmp(onsetTitle, 'Stimulus_onset'));
                n_ISI = find(strcmp(onsetTitle, 'ISI_onset'));
                n = find(strcmp(onsetTitle, 'Stimulus_duration'));
                allOnset{b}(:, n) = allOnset{b}(:, n_ISI)-allOnset{b}(:, n_Stim);
                
                %Rating_duration
                n = find(strcmp(onsetTitle, 'Rating_duration'));
                n_ITI=find(strcmp(onsetTitle, 'ITI_onset'));
                n_Rating = find(strcmp(onsetTitle, 'Rating_onset'));
                allOnset{b}(:, n) = allOnset{b}(:, n_ITI)-allOnset{b}(:, n_Rating);
                
                
                
                for t = 1:nTrial
                    
                    current_trial = current_trial + 1;
                    
                    % Block
                    n = find(strcmp(onsetTitle, 'block_info'));
                    allOnset{b}(current_trial, n) = b;
                    
                    % Trial
                    n = find(strcmp(onsetTitle, 'trial'));
                    allOnset{b}(current_trial, n) = t;
                    
                    % Type_Q
                    n = find(strcmp(onsetTitle, 'typeQ'));
                    allOnset{b}(current_trial, n) = data.S1_typeQ(t);%SE=0;RE=1
                    
                    % SE Info
                    n = find(strcmp(onsetTitle, 'SE_info'));
                    allOnset{b}(current_trial, n) = data.S1_block_pair_SE(t)+1;
                    
                    % RE Info
                    n = find(strcmp(onsetTitle, 'RE_info'));
                    allOnset{b}(current_trial, n) = data.S1_block_pair_RE(t)+1;
                    
                    %ISI duration ->Some trials don't have questions
                    n = find(strcmp(onsetTitle, 'ISI_duration'));
                    if data.S1_Question(t)==0 %No question trial
                        if current_trial==nTrial
                            current_ISI=data.S1_ISI(current_trial);
                        else
                            current_ISI=allOnset{b}(current_trial+1, n_Stim)-allOnset{b}(current_trial, n_ISI);
                        end
                    else
                        current_ISI=allOnset{b}(current_trial, n_Rating)-allOnset{b}(current_trial, n_ISI);
                    end
                    allOnset{b}(current_trial, n) =current_ISI;
                    
                    

                    
                    %Question
                    n = find(strcmp(onsetTitle, 'Question'));
                    allOnset{b}(current_trial, n)=data.S1_Question(t);
                    
                    %Rating
                    n = find(strcmp(onsetTitle, 'Rating'));
                    allOnset{b}(current_trial, n) = data.S1_block_Rating(t);
                    
                    %ITI duration
                    n = find(strcmp(onsetTitle, 'ITI_duration'));
                    if data.S1_Question(t)==0 %No question trial
                        current_ITI=0;
                    else
                        if current_trial==nTrial
                            current_ITI=data.S1_block_ITI(current_trial);
                        else
                            current_ITI=allOnset{b}(current_trial+1, n_Stim)-allOnset{b}(current_trial, n_ITI);
                        end
                    end
                    allOnset{b}(current_trial, n) =current_ITI;
                    
                end
            case 'SERE2'
                nTrial=72;
                
                onsetTitle = {...
                    'block_info', 'trial','typeQ','SE_info','RE_info','Stimulus_onset','Stimulus_duration',...
                    'ISI_onset','ISI_duration','Question','Rating','Rating_onset','Rating_duration',...
                    'ITI_onset','ITI_duration'
                    };
                nTitle = numel(onsetTitle);
                fmri_onset.onsetTitle{b} = onsetTitle;
                
                % start to organize data
                %allOnset = zeros(nBlock*nTrial, nTitle);
                allOnset{b} = zeros(1*nTrial, nTitle);
                current_trial = 0;
                
                %Calculated Time data
                %Stimulus_onset
                if strcmp(sublist{sub},'Sub044') %Sub006;Sub044
                    dataFile=fullfile(dir_data,sprintf('%s_2_2_Rating_MRI_SEREcombine_Eyetracking_session%d.csv',sublist{sub},2));
                    Timedata = readtable(dataFile, 'delimiter', ',');
                    Ready_text_line=find(~isnan(Timedata.Ready_text_started));
                    idx_CombineB2_started=Ready_text_line;
                    CombineB2_started=max(Timedata.center_OnsetTime_2_started(idx_CombineB2_started),Timedata.center_OnsetTime_started(idx_CombineB2_started));
                    idx_included=find(~isnan(max(Timedata.center_OnsetTime_started,Timedata.center_OnsetTime_2_started)));
                else
                idx_included=find(~isnan(max(Timedata.center_OnsetTime_started,Timedata.center_OnsetTime_2_started)));
                idx_included=idx_included(idx_included>=idx_CombineB2_started);
                end
                n = find(strcmp(onsetTitle, 'Stimulus_onset'));
                StimOnsetTime=max(Timedata.center_OnsetTime_started(idx_included),Timedata.center_OnsetTime_2_started(idx_included));
                allOnset{b}(:, n) = StimOnsetTime-CombineB2_started;
                
                
                %ISI_onset
                n = find(strcmp(onsetTitle, 'ISI_onset'));
                switch StructData.Run_hand
                    case 0
                        ISIOnsetTime = max(Timedata.fixOnset_text_ISI_SE_2_started(idx_included),Timedata.fixOnset_text_ISI_RE_2_started(idx_included));
                    case 1
                        ISIOnsetTime = max(Timedata.fixOnset_text_ISI_SE_started(idx_included),Timedata.fixOnset_text_ISI_RE_started(idx_included));
                end
                allOnset{b}(:, n) =ISIOnsetTime-CombineB2_started;
                
                %Rating_onset  ->Not all trials have question
                n = find(strcmp(onsetTitle, 'Rating_onset'));
                switch StructData.Run_hand
                    case 0
                        QuestionOnsetTime=max(Timedata.key_resp_rating_SE_L_started(idx_included),Timedata.key_resp_rating_RE_L_started(idx_included));
                        
                    case 1
                        QuestionOnsetTime=max(Timedata.key_resp_rating_SE_R_started(idx_included),Timedata.key_resp_rating_RE_R_started(idx_included));
                end
                RatingOnset=QuestionOnsetTime-CombineB2_started;
                RatingOnset(isnan(RatingOnset))=0;
                allOnset{b}(:, n) = RatingOnset;
                
                %ITI_onset
                n = find(strcmp(onsetTitle, 'ITI_onset'));
                ITIOnset=Timedata.fixOnset_text_ITI_started(idx_included)-CombineB2_started;
                ITIOnset(isnan(ITIOnset))=0;
                allOnset{b}(:, n) = ITIOnset;
                
                
                %Stimulus duration
                n_Stim = find(strcmp(onsetTitle, 'Stimulus_onset'));
                n_ISI = find(strcmp(onsetTitle, 'ISI_onset'));
                n = find(strcmp(onsetTitle, 'Stimulus_duration'));
                allOnset{b}(:, n) = allOnset{b}(:, n_ISI)-allOnset{b}(:, n_Stim);
                
                %Rating_duration
                n = find(strcmp(onsetTitle, 'Rating_duration'));
                n_ITI=find(strcmp(onsetTitle, 'ITI_onset'));
                n_Rating = find(strcmp(onsetTitle, 'Rating_onset'));
                allOnset{b}(:, n) = allOnset{b}(:, n_ITI)-allOnset{b}(:, n_Rating);
                
                
                
                for t = 1:nTrial
                    
                    current_trial = current_trial + 1;
                    
                    % Block
                    n = find(strcmp(onsetTitle, 'block_info'));
                    allOnset{b}(current_trial, n) = b;
                    
                    % Trial
                    n = find(strcmp(onsetTitle, 'trial'));
                    allOnset{b}(current_trial, n) = t;
                    
                    % Type_Q
                    n = find(strcmp(onsetTitle, 'typeQ'));
                    allOnset{b}(current_trial, n) = data.S2_typeQ(t);%SE=0;RE=1
                    
                    % SE Info
                    n = find(strcmp(onsetTitle, 'SE_info'));
                    allOnset{b}(current_trial, n) = data.S2_block_pair_SE(t)+1;
                    
                    % RE Info
                    n = find(strcmp(onsetTitle, 'RE_info'));
                    allOnset{b}(current_trial, n) = data.S2_block_pair_RE(t)+1;
                    
                    %ISI duration ->Some trials don't have questions
                    n = find(strcmp(onsetTitle, 'ISI_duration'));
                   if data.S2_Question(t)==0 %No question trial
                        if current_trial==nTrial
                            current_ISI=data.S2_ISI(current_trial);
                        else
                            current_ISI=allOnset{b}(current_trial+1, n_Stim)-allOnset{b}(current_trial, n_ISI);
                        end
                    else
                        current_ISI=allOnset{b}(current_trial, n_Rating)-allOnset{b}(current_trial, n_ISI);
                    end
                    allOnset{b}(current_trial, n) =current_ISI;
                    
                    %Question
                    n = find(strcmp(onsetTitle, 'Question'));
                    allOnset{b}(current_trial, n)=data.S2_Question(t);
                    
                    %Rating
                    n = find(strcmp(onsetTitle, 'Rating'));
                    allOnset{b}(current_trial, n) = data.S2_block_Rating(t);
                    
                    %ITI duration
                    n = find(strcmp(onsetTitle, 'ITI_duration'));
                    if data.S2_Question(t)==0 %No question trial
                        current_ITI=0;
                    else
                        if current_trial==nTrial
                            current_ITI=data.S2_block_ITI(current_trial);
                        else
                            current_ITI=allOnset{b}(current_trial+1, n_Stim)-allOnset{b}(current_trial, n_ITI);
                        end
                    end
                    allOnset{b}(current_trial, n) =current_ITI;
                    
                end
            case 'Play'
                nTrial=72;
                onsetTitle = {...
                    'block_info', 'trial','SE_info','RE_info','Stimulus_onset','Stimulus_duration',...
                    'ISI1_onset','ISI1_duration','Question','Rating_Play','Play_onset','Play_duration',...
                    'ISI2_onset','ISI2_duration','Rating_Attri','Attri_onset','Attri_duration',...
                    'ITI_onset','ITI_duration'
                    };
                nTitle = numel(onsetTitle);
                fmri_onset.onsetTitle{b} = onsetTitle;
                
                % start to organize data
                %allOnset = zeros(nBlock*nTrial, nTitle);
                allOnset{b} = zeros(1*nTrial, nTitle);
                current_trial = 0;
                %Calculated Time data
                %Stimulus_onset
                idx_included=find(~isnan(Timedata_Play.center_OnsetTime_started));
                n = find(strcmp(onsetTitle, 'Stimulus_onset'));
                allOnset{b}(:, n) = Timedata_Play.center_OnsetTime_started(idx_included)-Play_Block_started;
                
                
                %ISI1_onset
                n = find(strcmp(onsetTitle, 'ISI1_onset'));
                switch StructData.Run_hand
                    case 0
                        allOnset{b}(:, n) = Timedata_Play.fixOnset_text_5_started(idx_included)-Play_Block_started;
                    case 1
                        allOnset{b}(:, n) = Timedata_Play.fixOnset_text_4_started(idx_included)-Play_Block_started;
                end
                
                %Stimulus duration
                n = find(strcmp(onsetTitle, 'Stimulus_duration'));
                n_ISI=find(strcmp(onsetTitle, 'ISI1_onset'));
                n_Stim = find(strcmp(onsetTitle, 'Stimulus_onset'));
                allOnset{b}(:, n) = allOnset{b}(:, n_ISI)-allOnset{b}(:, n_Stim);
                
                %Play_onset
                n = find(strcmp(onsetTitle, 'Play_onset'));
                switch StructData.Run_hand
                    case 0
                        PlayOnset=Timedata_Play.PairQuestionONSET_text_2_started(idx_included)-Play_Block_started;
                    case 1
                        PlayOnset=Timedata_Play.PairQuestionONSET_text_4_started(idx_included)-Play_Block_started;
                end
                PlayOnset(isnan(PlayOnset))=0;
                allOnset{b}(:, n) =PlayOnset;
                
                %ISI2_onset
                n = find(strcmp(onsetTitle, 'ISI2_onset'));
                ISI2Onset=Timedata_Play.fixOnset_text_started(idx_included)-Play_Block_started;
                ISI2Onset(isnan(ISI2Onset))=0;
                allOnset{b}(:, n) =ISI2Onset;
                
                %Play_duration
                n = find(strcmp(onsetTitle, 'Play_duration'));
                switch StructData.Run_hand
                    case 0
                        allOnset{b}(:, n) = Timedata_Play.fixOnset_text_started(idx_included)-Timedata_Play.PairQuestionONSET_text_2_started(idx_included);
                    case 1
                        allOnset{b}(:, n) = Timedata_Play.fixOnset_text_started(idx_included)-Timedata_Play.PairQuestionONSET_text_4_started(idx_included);
                end
                
                
                %Attri_onset
                n = find(strcmp(onsetTitle, 'Attri_onset'));
                switch StructData.Run_hand
                    case 0
                        AttriOnset=Timedata_Play.PairQuestionONSET_attri_2_started(idx_included)-Play_Block_started;
                        
                    case 1
                        
                        AttriOnset=Timedata_Play.PairQuestionONSET_attri_started(idx_included)-Play_Block_started;
                end
                AttriOnset(isnan(AttriOnset))=0;
                allOnset{b}(:, n) = AttriOnset;
                

                %ISI2_duration
                n_Attri = find(strcmp(onsetTitle, 'Attri_onset'));
                n_ISI2=find(strcmp(onsetTitle, 'ISI2_onset'));
                n = find(strcmp(onsetTitle, 'ISI2_duration'));
                allOnset{b}(:, n) = allOnset{b}(:, n_Attri)-allOnset{b}(:, n_ISI2);
                
                %ITI_onset
                n = find(strcmp(onsetTitle, 'ITI_onset'));
                switch StructData.Run_hand
                    case 0
                        ITIOnset=Timedata_Play.fixOnset_text_ITI_started(idx_included)-Play_Block_started;
                    case 1
                        ITIOnset=Timedata_Play.fixOnset_text_ITI_started(idx_included)-Play_Block_started;
                end
                ITIOnset(isnan(ITIOnset))=0;
                allOnset{b}(:, n) = ITIOnset;
                
                
                %Attri_duration
                n_Attri = find(strcmp(onsetTitle, 'Attri_onset'));
                n_ITI=find(strcmp(onsetTitle, 'ITI_onset'));
                n = find(strcmp(onsetTitle, 'Attri_duration'));
                allOnset{b}(:, n) = allOnset{b}(:, n_ITI)-allOnset{b}(:, n_Attri);
                
                for t = 1:nTrial
                    
                    current_trial = current_trial + 1;
                    
                    % Block
                    n = find(strcmp(onsetTitle, 'block_info'));
                    allOnset{b}(current_trial, n) = b;
                    
                    % Trial
                    n = find(strcmp(onsetTitle, 'trial'));
                    allOnset{b}(current_trial, n) = t;
                    
                    % SE Info
                    n = find(strcmp(onsetTitle, 'SE_info'));
                    allOnset{b}(current_trial, n) = data.Play_block_pair_SE(t)+1;
                    
                    % RE Info
                    n = find(strcmp(onsetTitle, 'RE_info'));
                    allOnset{b}(current_trial, n) = data.Play_block_pair_RE(t)+1;
                    
                    
                    %ISI1 duration->Some trials don't have questions
                    n = find(strcmp(onsetTitle, 'ISI1_duration'));
                    n_Stim = find(strcmp(onsetTitle, 'Stimulus_onset'));
                    n_ISI=find(strcmp(onsetTitle, 'ISI1_onset'));
                    n_PlayOnset = find(strcmp(onsetTitle, 'Play_onset'));
                    if data.Play_block_Q(t)==0 %No question trial
                        if current_trial==nTrial
                            current_ISI=data.Play_block_ISI1(current_trial);
                        else
                            current_ISI=allOnset{b}(current_trial+1, n_Stim)-allOnset{b}(current_trial, n_ISI);
                        end
                    else
                        current_ISI=allOnset{b}(current_trial, n_PlayOnset)-allOnset{b}(current_trial, n_ISI);
                    end
                    allOnset{b}(current_trial, n) =current_ISI;
                    
                    
                    
                    %Question
                    n = find(strcmp(onsetTitle, 'Question'));
                    allOnset{b}(current_trial, n)=data.Play_block_Q(t);
                    
                    %Rating_Play
                    n = find(strcmp(onsetTitle, 'Rating_Play'));
                    allOnset{b}(current_trial, n) = data.Play_block_Play(t);
                    
                    
                    %Rating_Attri
                    n = find(strcmp(onsetTitle, 'Rating_Attri'));
                    allOnset{b}(current_trial, n) = data.Play_block_Attri(t);
                    
                    %ITI_duration
                    n = find(strcmp(onsetTitle, 'ITI_duration'));
                    n_ITI = find(strcmp(onsetTitle, 'ITI_onset'));
                    if data.Play_block_Q(t)==0 %No question trial
                        current_ITI=0;
                    else
                        if current_trial==nTrial
                            current_ITI=data.Play_block_ITI(current_trial);
                        else
                            current_ITI=allOnset{b}(current_trial+1, n_Stim)-allOnset{b}(current_trial, n_ITI);
                        end
                    end
                    allOnset{b}(current_trial, n) =current_ITI;
                    
                    
                    
                end
                
        end
    end
    
    
    % save onsets
    fmri_onset.allOnset = allOnset;
    save(onsetFile,'fmri_onset');
end



