function reorganized_data(version,Sublist)

% check the input arguments
% subject list
if nargin<3
    sublist_name=sprintf('sublist_ver%d.txt',version);
    sublist = textread(sublist_name,'%s','delimiter','\n');
    
else
    if ~iscell(sublist)
        sublist = {sublist};
    end
end

nSubj = numel(sublist);

% directory
dir_fig = 'figures';
dir_data='data';
dir_rebuilt_data='organized';
mkdir(dir_rebuilt_data);

speed=[7];
target_size_ori=[10,20,30];
prob_card=[20,50,80];
block_sequence=[0,1,2];
Rating_pair={'Extremely Unlikely','Slightly Unlikely','Slightly Likely','Extremely Likely'};
Rating_Choice={'Definitely No','No','Yes','Definitely Yes'};
Attri_Choice={'Definitely Card','Card','Bar','Definitely Bar'};
for s=1:nSubj
    %% for practice data
    
    dataFile=fullfile(dir_data,sprintf('%s_1trainingbehavior_beforeMRI.csv',sublist{s}));
    fprintf('%s',sublist{s})
    data = readtable(dataFile, 'delimiter', ',');
    
    StructData.Run_hand=data.run_R(find(~isnan(data.run_R))); %if Run_hand=0: Left hand; if Run_hand =1: Right hand;
    StructData.target_size_practice=data.Target_Size(find(~isnan(data.Target_Size)));
    
    %Calculate the information from Practice Data
    StructData.target_size_ori=target_size_ori;
    StructData.sigma=data.sigma_indi(find(~isnan(data.sigma_indi)));
    
    for t=1:size(target_size_ori,2)
        StructData.practice_hit(1,t)=mean(data.key_resp_bar_hit(find(data.Current_Target_Size==target_size_ori(t))));
        StructData.practice_position(:,t)=data.key_resp_bar_position(find(data.Current_Target_Size==target_size_ori(t)));
        StructData.practice_RT(:,t)=data.key_resp_bar_rt(find(data.Current_Target_Size==target_size_ori(t)));
        if strcmp(sublist{s},'Sub006')
            switch StructData.Run_hand
                case 0
                    Q_row=find(~isnan(data.Slider_LikelyBar_L_response));
                    idx_Q = (data.Target_Size(find(~isnan(data.Slider_LikelyBar_L_response)))==target_size_ori(t));
                    StructData.practiceRating(1,t)=data.Slider_LikelyBar_L_response(Q_row(idx_Q));
                    StructData.practiceRatingRT(1,t)=data.Slider_LikelyBar_L_rt(Q_row(idx_Q));
                case 1
                    Q_row=find(~isnan(data.Slider_LikelyBar_R_response));
                    idx_Q = (data.Target_Size(find(~isnan(data.Slider_LikelyBar_R_response)))==target_size_ori(t));
                    StructData.practiceRating(1,t)=(data.Slider_LikelyBar_R_response(Q_row(idx_Q)));
                    StructData.practiceRatingRT(1,t)=data.Slider_LikelyBar_R_rt(Q_row(idx_Q));
            end
        else
            
            
            switch StructData.Run_hand
                case 0
                    Q_row=find(~isnan(data.key_resp_LikeQ_bar_L_keys));
                    idx_Q = (data.Target_Size(find(~isnan(data.key_resp_LikeQ_bar_L_keys)))==target_size_ori(t));
                    StructData.practiceRating(1,t)=data.key_resp_LikeQ_bar_L_keys(Q_row(idx_Q))-4;
                    StructData.practiceRatingRT(1,t)=data.key_resp_LikeQ_bar_L_rt(Q_row(idx_Q));
                case 1
                    Q_row=find(~isnan(data.key_resp_LikeQ_bar_R_keys));
                    idx_Q = (data.Target_Size(find(~isnan(data.key_resp_LikeQ_bar_R_keys)))==target_size_ori(t));
                    StructData.practiceRating(1,t)=(data.key_resp_LikeQ_bar_R_keys(Q_row(idx_Q))-4)*1;
                    StructData.practiceRatingRT(1,t)=(data.key_resp_LikeQ_bar_R_rt(Q_row(idx_Q)));
            end
            
        end
        
        
    end
    
    %Organized data for Card Practice
    for q=1:size(prob_card,2)
        if strcmp(sublist{s},'Sub006')
            switch StructData.Run_hand
                case 0
                    StructData.likely_card_Rating(1,q)= data.Slider_LikelyCard_L_response(find(data.Q_card_prob==prob_card(q)));
                    StructData.likely_card_RatingRT(1,q)= data.Slider_LikelyCard_L_rt(find(data.Q_card_prob==prob_card(q)));
                case 1
                    StructData.likely_card_Rating(1,q)= data.Slider_LikelyCard_R_response(find(data.Q_card_prob==prob_card(q)));
                    StructData.likely_card_RatingRT(1,q)= data.Slider_LikelyCard_R_rt(find(data.Q_card_prob==prob_card(q)));
                    
            end
        else
            switch StructData.Run_hand
                case 0
                    Q_row=find(~isnan(data.key_resp_LikeQ_card_L_keys));
                    idx_Q = (data.Q_card_prob(find(~isnan(data.key_resp_LikeQ_card_L_keys)))==prob_card(q));
                    StructData.likely_card_Rating(1,q)=data.key_resp_LikeQ_card_L_keys(Q_row(idx_Q))-4;
                    StructData.likely_card_RatingRT(1,q)=data.key_resp_LikeQ_card_L_rt(Q_row(idx_Q));
                case 1
                    Q_row=find(~isnan(data.key_resp_LikeQ_card_R_keys));
                    idx_Q = (data.Q_card_prob(find(~isnan(data.key_resp_LikeQ_card_R_keys)))==prob_card(q));
                    StructData.likely_card_Rating(1,q)=(data.key_resp_LikeQ_card_R_keys(Q_row(idx_Q))-4)*1;
                    StructData.likely_card_RatingRT(1,q)=(data.key_resp_LikeQ_card_R_rt(Q_row(idx_Q)));
            end
        end
        
        
    end
    
    
    %Practice in MRI
    
    
    dataFile=fullfile(dir_data,sprintf('%s_1trainingbehavior_inMRI.csv',sublist{s}));
    data = readtable(dataFile, 'delimiter', ',');
    
    %last stage
    n_practice=20;
    StructData.Practice2_Indi_target_size=unique(data.Current_Target_Size(~isnan(data.Current_Target_Size)));
    StructData.last_position=zeros(n_practice,size(StructData.Practice2_Indi_target_size,1));
    StructData.last_RT=zeros(n_practice,size(StructData.Practice2_Indi_target_size,1));
    StructData.PracMRI_Bar_Rating=zeros(size(StructData.Practice2_Indi_target_size,1),1);
    
    for t=1:size(StructData.Practice2_Indi_target_size,1)
        StructData.last_hit(1,t)=mean(data.key_resp_bar_hit(find(data.Current_Target_Size==StructData.Practice2_Indi_target_size(t))));
        StructData.last_position(:,t)=data.key_resp_bar_position(find(data.Current_Target_Size==StructData.Practice2_Indi_target_size(t)));
        StructData.last_RT(:,t)=data.key_resp_bar_rt(find(data.Current_Target_Size==StructData.Practice2_Indi_target_size(t)));
        
        
        switch StructData.Run_hand
            case 0
                Q_row=find(~isnan(data.key_resp_likely_bar_L_rt));
                idx_Q = (data.Target_Size(find(~isnan(data.key_resp_likely_bar_L_rt)))==target_size_ori(t));
                for n=1:size(Rating_pair,2)
                    if ~isempty(find(strcmp(cell2mat(data.Rating_Likely_bar(Q_row(idx_Q))),Rating_pair{n})))
                        StructData.PracMRI_Bar_Rating(t,1)=n;
                    end
                end
                StructData.Prac_Bar_RatingRT(1,t)=(data.key_resp_likely_bar_L_rt(Q_row(idx_Q)));
            case 1
                Q_row=find(~isnan(data.key_resp_likely_bar_R_rt));
                idx_Q = (data.Target_Size(find(~isnan(data.key_resp_likely_bar_R_rt)))==target_size_ori(t));
                for n=1:size(Rating_pair,2)
                    if ~isempty(find(strcmp(cell2mat(data.Rating_Likely_bar(Q_row(idx_Q))),Rating_pair{n})))
                        StructData.PracMRI_Bar_Rating(t,1)=n;
                    end
                end
                StructData.PracMRI_Bar_RatingRT(1,t)=(data.key_resp_likely_bar_R_rt(Q_row(idx_Q)));
        end
    end
    
    
    %Organized data for Card Practice
    StructData.PracMRI_Card_Rating=zeros(size(prob_card,2),1);
    for q=1:size(prob_card,2)
        switch StructData.Run_hand
            case 0
                Q_row=find(~isnan(data.key_resp_likely_card_L_rt));
                idx_Q=(data.Q_card_prob(find(~isnan(data.Q_card_prob)))==prob_card(q));
                for n=1:size(Rating_pair,2)
                    if ~isempty(find(strcmp(cell2mat(data.Rating_Likely_card(Q_row(idx_Q))),Rating_pair{n})))
                        StructData.PracMRI_Card_Rating(q,1)=n;
                    end
                end
                StructData.PracMRI_Card_RatingRT(q,1)=(data.key_resp_likely_card_L_rt(Q_row(idx_Q)));
                
            case 1
                Q_row=find(~isnan(data.key_resp_likely_card_R_rt));
                idx_Q=(data.Q_card_prob(find(~isnan(data.Q_card_prob)))==prob_card(q));
                for n=1:size(Rating_pair,2)
                    if ~isempty(find(strcmp(cell2mat(data.Rating_Likely_card(Q_row(idx_Q))),Rating_pair{n})))
                        StructData.PracMRI_Card_Rating(q,1)=n;
                    end
                end
                StructData.PracMRI_Card_RatingRT(q,1)=(data.key_resp_likely_card_R_rt(Q_row(idx_Q)));
        end
    end
    
    
    %% for MRI data
    
    %For the SE RE only session
    dataFile=fullfile(dir_data,sprintf('%s_2_1_MRI_Rating_SE_RE.csv',sublist{s}));
    fprintf('%s',sublist{s})
    data = readtable(dataFile, 'delimiter', ',');
    
    %organized data
    StructData.Only_block_trialN=[1:90]';
    %for bar trial
    idx_bar=find(~isnan(data.target_pair));
    %for card trial
    idx_card=find(~isnan(data.card_pair));
    
    %idx_all
    idx_all=sort([idx_bar;idx_card]);
    
    switch StructData.Run_hand
        case 0
            TrialN=data.Bar_Card_L_thisTrialN;
        case 1
            TrialN=data.Bar_Card_R_thisTrialN;
    end
    %Q_type: 0:bar 1:card
    StructData.Only_block_type=zeros(size(StructData.Only_block_trialN,1),1);
    StructData.Only_block_type(TrialN(~isnan(data.card_pair))+1,1)=1;
    %Question or not
    StructData.Only_block_Q=zeros(size(StructData.Only_block_trialN,1),1);
    trial_barQ=TrialN(data.Q_bar==1 & ~isnan(data.ISI_SE))+1;
    trial_cardQ=TrialN(data.Q_card==1 & ~isnan(data.ISI_RE))+1;
    StructData.Only_block_Q(sort([trial_barQ,trial_cardQ]),1)=1;
    %IF no SE trial, current pair SE=-1;
    StructData.Only_block_pair_SE=ones(size(StructData.Only_block_trialN,1),1)*-1;
    Bar_trial=TrialN(~isnan(data.target_pair))+1;
    StructData.Only_block_pair_SE(Bar_trial,1)=data.current_pair_SE(~isnan(data.ISI_SE));
    %IF no RE trial, current pair RE=-1;
    StructData.Only_block_pair_RE=ones(size(StructData.Only_block_trialN,1),1)*-1;
    Card_trial=TrialN(~isnan(data.card_pair))+1;
    StructData.Only_block_pair_RE(Card_trial,1)=data.current_pair_RE(~isnan(data.ISI_RE));
    %ISI
    StructData.Only_block_ISI=ones(size(StructData.Only_block_trialN,1),1)*-1;
    StructData.Only_block_ISI(Bar_trial,1)=data.ISI_SE(~isnan(data.ISI_SE));
    StructData.Only_block_ISI(Card_trial,1)=data.ISI_RE(~isnan(data.ISI_RE));
    %ITI
    StructData.Only_block_ITI(:,1)=data.current_ITI(idx_all);
    StructData.Only_block_ITI(isnan(StructData.Only_block_ITI),1)=-1;
    %
    StructData.Only_block_Rating=ones(size(StructData.Only_block_trialN,1),1)*-1;
    for n=1:size(Rating_pair,2)
        StructData.Only_block_Rating(find(strcmp(data.Rating_pair(idx_all),Rating_pair{n})),1)=n;
    end
    StructData.Only_block_RatingRT=ones(size(StructData.Only_block_trialN,1),1)*-1;
    switch StructData.Run_hand
        case 0
            StructData.Only_block_RatingRT=max(data.key_resp_rating_SE_L_rt(idx_all),data.key_resp_rating_RE_L_rt(idx_all));
            StructData.Only_block_RatingRT(find(isnan(StructData.Only_block_RatingRT)),1)=-1;
        case 1
            StructData.Only_block_RatingRT=max(data.key_resp_rating_SE_R_rt(idx_all),data.key_resp_rating_RE_R_rt(idx_all));
            StructData.Only_block_RatingRT(find(isnan(StructData.Only_block_RatingRT)),1)=-1;
    end
    
    
    
    
    
    
    %For the SE_RE combine session
    if strcmp(sublist{s},'Sub044') %Sub006; Sub049
        dataFile=fullfile(dir_data,sprintf('%s_2_2_Rating_MRI_SEREcombine_Eyetracking_session%d.csv',sublist{s},1));
    else
        dataFile=fullfile(dir_data,sprintf('%s_2_2_Rating_MRI_SEREcombine_Eyetracking.csv',sublist{s}));
        
    end
    fprintf('%s',sublist{s})
    data = readtable(dataFile, 'delimiter', ',');
    %organozed Data
    %session1
    StructData.S1_block_trialN=[1:72]';
    info_block=find(~isnan(data.Ready_text_started));
    idx=sort([find(~isnan(data.center_OnsetTime_2_started));find(~isnan(data.center_OnsetTime_started))]);
    if strcmp(sublist{s},'Sub044') %Sub006;Sub049
        idx=idx;
    else
        idx=idx(idx<info_block(2));
    end
    StructData.S1_block_pair_combi=data.current_pair(idx);
    StructData.S1_block_pair_SE= data.current_pair_SE(idx);
    StructData.S1_block_pair_RE= data.current_pair_RE(idx);
    StructData.S1_typeQ=data.type_Q(idx); %0 for SE; 1 for RE
    StructData.S1_typeQ(find(isnan(StructData.S1_typeQ)),1)=1;
    StructData.S1_ISI=zeros(size(StructData.S1_block_trialN,1),1);
    StructData.S1_ISI(find(StructData.S1_typeQ==0),1)=data.ISI_SE(idx(find(StructData.S1_typeQ==0)),1);
    StructData.S1_ISI(find(StructData.S1_typeQ==1),1)=data.ITI_RE(idx(find(StructData.S1_typeQ==1)),1);
    StructData.S1_Question=zeros(size(StructData.S1_block_trialN,1),1);
    StructData.S1_Question(find(StructData.S1_typeQ==0),1)=data.SE_block_Q(idx(find(StructData.S1_typeQ==0)),1);
    StructData.S1_Question(find(StructData.S1_typeQ==1),1)=data.RE_block_Q(idx(find(StructData.S1_typeQ==1)),1);
    StructData.S1_block_ITI=data.current_ITI(idx);
    StructData.S1_block_ITI(isnan(StructData.S1_block_ITI),1)=0;
    
    StructData.S1_block_Rating=zeros(size(StructData.S1_block_trialN,1),1);
    for n=1:size(Rating_pair,2)
        StructData.S1_block_Rating(find(strcmp(data.Rating_pair(idx),Rating_pair{n})),1)=n;
    end
    StructData.S1_block_RatingRT=zeros(size(StructData.S1_block_trialN,1),1);
    switch StructData.Run_hand
        case 0
            StructData.S1_block_RatingRT(find(StructData.S1_typeQ==0),1)=data.key_resp_rating_SE_L_rt(idx(find(StructData.S1_typeQ==0)),1);
            StructData.S1_block_RatingRT(find(StructData.S1_typeQ==1),1)=data.key_resp_rating_RE_L_rt(idx(find(StructData.S1_typeQ==0)),1);
        case 1
            StructData.S1_block_RatingRT(find(StructData.S1_typeQ==0),1)=data.key_resp_rating_SE_R_rt(idx(find(StructData.S1_typeQ==0)),1);
            StructData.S1_block_RatingRT(find(StructData.S1_typeQ==1),1)=data.key_resp_rating_RE_R_rt(idx(find(StructData.S1_typeQ==0)),1);
    end
    StructData.S1_block_RatingRT(isnan(StructData.S1_block_RatingRT),1)=0;
    
    
    %session2
    if strcmp(sublist{s},'Sub044') %Sub006; Sub049
        dataFile=fullfile(dir_data,sprintf('%s_2_2_Rating_MRI_SEREcombine_Eyetracking_session%d.csv',sublist{s},2));
    end
    fprintf('%s',sublist{s})
    data = readtable(dataFile, 'delimiter', ',');
    %organozed Data
    %session1
    StructData.S2_block_trialN=[1:72]';
    idx=sort([find(~isnan(data.center_OnsetTime_2_started));find(~isnan(data.center_OnsetTime_started))]);
    if strcmp(sublist{s},'Sub044') %Sub006 ;Sub049
        idx=idx;
    else
        idx=idx(idx>=info_block(2));
    end
    StructData.S2_block_pair_combi=data.current_pair(idx);
    StructData.S2_block_pair_SE= data.current_pair_SE(idx);
    StructData.S2_block_pair_RE= data.current_pair_RE(idx);
    StructData.S2_typeQ=data.type_Q(idx); %0 for SE; 1 for RE
    StructData.S2_typeQ(find(isnan(StructData.S2_typeQ)),1)=1;
    StructData.S2_ISI=zeros(size(StructData.S2_block_trialN,1),1);
    StructData.S2_ISI(find(StructData.S2_typeQ==0),1)=data.ISI_SE(idx(find(StructData.S2_typeQ==0)),1);
    StructData.S2_ISI(find(StructData.S2_typeQ==1),1)=data.ITI_RE(idx(find(StructData.S2_typeQ==1)),1);
    StructData.S2_Question=zeros(size(StructData.S2_block_trialN,1),1);
    StructData.S2_Question(find(StructData.S2_typeQ==0),1)=data.SE_block_Q(idx(find(StructData.S2_typeQ==0)),1);
    StructData.S2_Question(find(StructData.S2_typeQ==1),1)=data.RE_block_Q(idx(find(StructData.S2_typeQ==1)),1);
    StructData.S2_block_ITI=data.current_ITI(idx);
    StructData.S2_block_ITI(isnan(StructData.S2_block_ITI),1)=0;
    
    StructData.S2_block_Rating=zeros(size(StructData.S2_block_trialN,1),1);
    for n=1:size(Rating_pair,2)
        StructData.S2_block_Rating(find(strcmp(data.Rating_pair(idx),Rating_pair{n})),1)=n;
    end
    StructData.S2_block_RatingRT=zeros(size(StructData.S2_block_trialN,1),1);
    switch StructData.Run_hand
        case 0
            StructData.S2_block_RatingRT(find(StructData.S2_typeQ==0),1)=data.key_resp_rating_SE_L_rt(idx(find(StructData.S2_typeQ==0)),1);
            StructData.S2_block_RatingRT(find(StructData.S2_typeQ==1),1)=data.key_resp_rating_RE_L_rt(idx(find(StructData.S2_typeQ==0)),1);
        case 1
            StructData.S2_block_RatingRT(find(StructData.S2_typeQ==0),1)=data.key_resp_rating_SE_R_rt(idx(find(StructData.S2_typeQ==0)),1);
            StructData.S2_block_RatingRT(find(StructData.S2_typeQ==1),1)=data.key_resp_rating_RE_R_rt(idx(find(StructData.S2_typeQ==0)),1);
    end
    StructData.S2_block_RatingRT(isnan(StructData.S2_block_RatingRT),1)=0;
    
    
    
    
    %For the play session
    dataFile=fullfile(dir_data,sprintf('%s_2_3_MRI_Rating_Play_Eyetracking.csv',sublist{s}));
    fprintf('%s',sublist{s})
    data = readtable(dataFile, 'delimiter', ',');
    %organozed Data
    StructData.Play_block_trialN=[1:72]';
    idx=find(~isnan(data.center_OnsetTime_started));
    StructData.Play_block_pair_combi=data.current_pair(idx);
    StructData.Play_block_pair_SE= data.current_pair_SE(idx);
    StructData.Play_block_pair_RE= data.current_pair_RE(idx);
    StructData.Play_block_ISI1=data.ISI1_Play(idx);
    StructData.Play_block_ISI2=data.ISI2_Play(idx);
    StructData.Play_block_Play=zeros(size(StructData.Play_block_trialN,1),1);
    for n=1:size(Rating_Choice,2)
        StructData.Play_block_Play(find(strcmp(data.Rating_pair_Play(idx),Rating_Choice{n})),1)=n;
    end
    StructData.Play_block_Attri=zeros(size(StructData.Play_block_trialN,1),1);
    for n=1:size(Attri_Choice,2)
        StructData.Play_block_Attri(find(strcmp(data.Attri_Ans(idx),Attri_Choice{n})),1)=n;
    end
    switch StructData.Run_hand
        case 0
            %%%%%ITI need to modify after chang the exp code
            StructData.Play_block_ITI=data.ITI_Play(idx)+2-data.key_resp_rating_Attri_L2_rt(idx);
            %%%%%
            StructData.Play_block_PlayRT=data.key_resp_rating_play_L_rt(idx);
            StructData.Play_block_AttriRT=data.key_resp_rating_Attri_L2_rt(idx);
        case 1
            %%%%%ITI need to modify after chang the exp code
            StructData.Play_block_ITI=data.ITI_actual_Play(idx);
            %%%%%
            StructData.Play_block_PlayRT=data.key_resp_rating_play_R_rt(idx);
            StructData.Play_block_AttriRT=data.key_resp_rating_Attri_R_rt(idx);
    end
    
    
    %Question or not
    StructData.Play_block_Q=data.Play_block_Q(idx);
    
    
%         %% Final behavior
        dataFile=fullfile(dir_data,sprintf('%s_3final_behavior_ButtonPress.csv',sublist{s}));
        fprintf('%s',sublist{s})
        data = readtable(dataFile, 'delimiter', ',');


    StructData.Indi_target_size=unique(data.target_pair_1(find(~isnan(data.target_pair_1))));
    %Choice
    ChoiceIndex=find(~isnan(data.ChoiceTrial_thisN));
    StructData.Choice_block_pair1SE=data.current_pair_1_SE(ChoiceIndex);
    StructData.Choice_block_pair1RE=data.current_pair_1_RE(ChoiceIndex);
    StructData.Choice_block_pair2SE=data.current_pair_2_SE(ChoiceIndex);
    StructData.Choice_block_pair2RE=data.current_pair_2_RE(ChoiceIndex);
    StructData.Choice_block_pair1_chosen=data.pair1_chosen(ChoiceIndex);
    StructData.Choice_block_pair2_chosen=data.pair2_chosen(ChoiceIndex);
    StructData.Choice_block_pair_chosen_RT=data.key_resp_choice_rt(ChoiceIndex);
    
     StructData.Indi_target_size=unique(data.Target_size_Confident(find(~isnan(data.Target_size_Confident))));

        %Questions
        idx=0;
        for i=1:size(StructData.Indi_target_size,1)
            for p=1:size(prob_card,2)


            % %Confident Question
            idx=idx+1;
            idx_confi=((data.Target_size_Confident==StructData.Indi_target_size(i))&(data.Prob_Confident==prob_card(p)));
            StructData.Confi_SE(i,p)=i-1; % Use [0,1,2] same as python
            StructData.Confi_RE(i,p)=p-1;
            
            %Likely Question
            idx_likely=((data.Target_size_likely==StructData.Indi_target_size(i))&(data.Prob_Likely==prob_card(p)));
            StructData.Likely_SE(i,p)=i-1; % Use [0,1,2] same as python
            StructData.Likely_RE(i,p)=p-1;
            % Control Question
            idx_control=((data.Target_size_Control==StructData.Indi_target_size(i))&(data.Prob_Control==prob_card(p)));
            StructData.Control_SE(i,p)=i-1; % Use [0,1,2] same as python
            StructData.Control_RE(i,p)=p-1;
   

                
                switch StructData.Run_hand
                    case 0
                        StructData.Confi_Rating(i,p)=(data.key_resp_Confi_L_keys(idx_confi)-4);
                        StructData.Confi_RatingRT(i,p)=data.key_resp_Confi_L_rt(idx_confi);
                        StructData.Likely_Rating(i,p)=(data.key_resp_Likely_L_keys(idx_likely)-4);
                        StructData.Likely_RatingRT(i,p)=data.key_resp_Likely_L_rt(idx_likely);
                        StructData.Control_Rating(i,p)=(data.key_resp_Control_L_keys(idx_control)-4);
                        StructData.Control_RatingRT(i,p)=data.key_resp_Control_L_rt(idx_control);
                    case 1
                        StructData.Confi_Rating(i,p)=(data.key_resp_Confi_R_keys(idx_confi)-4)*-1;
                        StructData.Confi_RatingRT(i,p)=data.key_resp_Confi_R_rt(idx_confi);
                        StructData.Likely_Rating(i,p)=(data.key_resp_Likely_R_keys(idx_likely)-4)*-1;
                        StructData.Likely_RatingRT(i,p)=data.key_resp_Likely_R_rt(idx_likely);
                        StructData.Control_Rating(i,p)=(data.key_resp_Control_R_keys(idx_control)-4)*-1;
                        StructData.Control_RatingRT(i,p)=data.key_resp_Control_R_rt(idx_control);
    
                end
    
            end
        end
    
        %betting
        BettingIndex=find(~isnan(data.current_pair));
        StructData.Betting_SE=data.current_pair_SE(BettingIndex);
        StructData.Betting_RE=data.current_pair_RE(BettingIndex);
        StructData.Betting=data.Slider_rating_betting_response(BettingIndex);
        StructData.Betting_rt=data.Slider_rating_betting_rt(BettingIndex);
    
        
    

    
        idx_isempty = cellfun(@isempty, data.key_resp_finalQ_keys, 'UniformOutput', false);
        idx_isempty = cell2mat(idx_isempty);
        StructData.key_resp_finalQ1=(data.key_resp_finalQ_keys{~idx_isempty});
    
        idx_isempty = cellfun(@isempty, data.key_resp_finalQ_2_keys, 'UniformOutput', false);
        idx_isempty = cell2mat(idx_isempty);
        StructData.key_resp_finalQ2=(data.key_resp_finalQ_2_keys{~idx_isempty});
    
        idx_isempty = cellfun(@isempty, data.key_resp_finalQ_3_keys, 'UniformOutput', false);
        idx_isempty = cell2mat(idx_isempty);
        StructData.key_resp_finalQ3=(data.key_resp_finalQ_3_keys{~idx_isempty});
%     
%     
    
    Filename=sprintf('%s_StructData',sublist{s});
    save(fullfile(dir_rebuilt_data,Filename),'StructData');
    
end



