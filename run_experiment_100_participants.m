% for each parameter set we run the 100 participants
% they share the internal model, they receive different ordering of the
% stimuli.

data_panel = [];
n_panel = [];
n2_panel = [];
params_pooled = 0; params_sig_pooled = 0;
nsubs = 100;
recog_all = [];

for irep = 1:nsubs % number of "subjects"/experiments
    % -------------------------------------------------
    % definition of model parameters
    % they define the listener/agent
    th_lip_k = Alip(1:3);
    th_f2_k = Af2(1:3);
    
    % initial values of model parameters. 
    % They will be updated after each perceptual trial. 
    th_lip_k0 = th_lip_k;    
    th_f2_k0 = th_f2_k;
    sig_Alip_k = sig_Alip*ones(3,1);
    sig_Alip_k0 = sig_Alip_k;
    sig_Af2_k = sig_Af2*ones(3,1);
    sig_Af2_k0 = sig_Af2_k;
       
    
    % exp_params (EXPERIMENTAL PARAMETERS)
    % exp_params(itrial,:) =  [s_vis s_ac Alip_val Af2_val]
    % s_vis = 0: trial with no visual stimulus, 1: with visual stimulus
    % s_ac  = 0: trial with no auditory stimulus, 1: with auditory stimulus    
   
    exp_params = experimental_paradigm(irep).exp_params;
    ind_exp = experimental_paradigm(irep).ind_exp;
    lip_randn = experimental_paradigm(irep).lip_randn;
    f2_randn = experimental_paradigm(irep).f2_randn;    
    
    % --------------------------------------------
    %nstim = size(exp_params,1);
    Ntrials = size(exp_params,1);  
    
    
    %vtem2 = zeros(Ntrials,2);
    params_adapt = zeros(Ntrials,6);
    params_adapt_sig = zeros(Ntrials,6);
    n_recog = [0 0 0];
    ind_c_over = zeros(1,Ntrials);
    
    clear Hsame pk 
                    
    % START THE EXPERIMENT
    for itrial = 1:Ntrials
        close all
        % read the stimulus features for the current trial:
        ind_stim = ind_exp(itrial);
        s_vis = exp_params(ind_stim,1); 
        s_ac = exp_params(ind_stim,2);
        Alip_exp = exp_params(ind_stim,3);
        Af2_exp = exp_params(ind_stim,4);
        
        % internal estimates of feature amplitudes at the end of the
        % stimulus presentation:
        s_bar_lip = Alip_exp + sig_s_lip*lip_randn(itrial);
        s_bar_f2 = Af2_exp + sig_s_f2*f2_randn(itrial);
        
        Aest(itrial,:) = [s_vis*s_bar_lip, s_bar_f2];
        
        % calculation of posterior probability
        % assuming common source:
        Hsame(:,itrial) =  exp(-s_vis*(th_lip_k-s_bar_lip).^2/2./(sig_s_lip^2+sig_Alip_k.^2) ...
            -s_ac*(th_f2_k-s_bar_f2).^2/2./(sig_s_f2^2+sig_Af2_k.^2))...
            ./(sqrt(sig_s_lip^2+sig_Alip_k.^2)*s_vis+(1-s_vis)) ...
            ./(sqrt(sig_s_f2^2+sig_Af2_k.^2)*s_ac+(1-s_ac))/3;
        
        % posterior probability for each of the three tokens
        pk(itrial,:) = softmax(log(Hsame(:,itrial)))';
        
        
        % ------------------------------------------------------
        % PERCEPTUAL DECISION BASED ON THE POSTERIOR pk(itrial,:)
        % a: the current trial's posterior probability
        % b: the current trial's percept (1: aba, 2: ada, 3: aga).
                
        [a b] = max(pk(end,:));
        a = pk(itrial,:);
            
        % vector keeping the percepts for the current listener/agent
        b_est(itrial) = b;
        % counter for number of times category 'b' has been perceived for a
        % given listener/agent.
        n_recog(b) = n_recog(b)+1;
        
        
        
        % ------------------------------------------------------
        % PARAMETER UPDATES        
        % ------------------------------------------------------
        % Standard Bayesian update
        if strcmp(update_method, 'standard Bayes')
            Alip_est = s_bar_lip;
            Af2_est = s_bar_f2;
            % a_th: "strength" of the prior for location parameters
            % a_sig:  "strength" of the prior for width parameters
            a_th = up_meth_params(isim,1);
            a_sig = up_meth_params(isim,2);
            % update for the spread parameters corrected on November 21st, 2017
            sig_Alip_k(b) = sqrt(sig_Alip_k(b)^2 + ...
                (a_th+n_recog(b))/(a_th+n_recog(b)+1)/(a_sig+n_recog(b))*s_vis*(Alip_est-th_lip_k(b))^2 ...
                -1/(a_sig+n_recog(b))*s_vis*(sig_Alip_k(b)^2));
            sig_Af2_k(b) = sqrt(sig_Af2_k(b)^2 + ...
                (a_th+n_recog(b))/(a_th+n_recog(b)+1)/(a_sig+n_recog(b))*s_ac*(Af2_est-th_f2_k(b))^2 ...
                -1/(a_sig+n_recog(b))*s_ac*(sig_Af2_k(b)^2));
            %update of the location parameters
            th_lip_k(b) = th_lip_k(b)+ 1/(a_th+n_recog(b))*s_vis*(Alip_est-th_lip_k(b));
            th_f2_k(b)  = th_f2_k(b)+ 1/(a_th+n_recog(b))*s_ac*(Af2_est-th_f2_k(b));     
      
        elseif strcmp(update_method, 'decaying accumulation')
            % -----------------------------------
            %  update_method = 'decaying accumulation';
            Alip_est = s_bar_lip;
            Af2_est = s_bar_f2;
            
            r_lip_1 = up_meth_params(isim,1);
            r_lip_2 = up_meth_params(isim,2);
            r_f2_1 = up_meth_params(isim,3);
            r_f2_2 = up_meth_params(isim,4);
            
            th_lip_k1 = th_lip_k + 1*r_lip_1*(th_lip_k0-th_lip_k); %decay
            th_lip_k(b) = th_lip_k1(b)+ 1*r_lip_2*s_vis*(Alip_est-th_lip_k(b))*a(b);
            th_lip_k(~([1:3]==b)) = th_lip_k1(~([1:3]==b));
                       
            th_f2_k1 = th_f2_k + 1*r_f2_1*(th_f2_k0-th_f2_k); % decay
            th_f2_k(b)  = th_f2_k1(b)+ 1*r_f2_2*s_ac*(Af2_est-th_f2_k(b))*a(b);
            th_f2_k(~([1:3]==b)) = th_f2_k1(~([1:3]==b));
            
        elseif strcmp(update_method, 'pure accumulation')
            % -----------------------------------
            %  update_method = 'pure accumulation'  (delta rule)
            Alip_est = s_bar_lip;
            Af2_est = s_bar_f2;
            r_del_lip = up_meth_params(isim, 1);
            r_del_f2 = up_meth_params(isim, 2);
            
            th_lip_k(b) = th_lip_k(b)+ 1*r_del_lip*s_vis*(Alip_est-th_lip_k(b))*a(b);            
            th_f2_k(b)  = th_f2_k(b)+ 1*r_del_f2*s_ac*(Af2_est-th_f2_k(b))*a(b);
                  
        elseif strcmp(update_method, 'no updates')
            % do nothing
        else
            disp('please define an update method')
            return
        end
        
        % ------------------------------------------------------
        % prevent cross-overs:
        
        if (th_f2_k(2) < th_f2_k(1)) & cr_over == 0
            th_f2_k([1 2]) = th_f2_k([2 1]);
            ind_c_over(itrial) = 1;
        end
        
        if (th_lip_k(2) > th_lip_k(3)) & cr_over == 0
            th_lip_k([2 3]) = th_lip_k([3 2]);
            ind_c_over(itrial) = 1;
        end
        
        % ------------------------------------------------------
        % save parameter updates
        params_adapt(itrial,:) = [th_lip_k'  th_f2_k'];
        params_adapt_sig(itrial,:) = [sig_Alip_k'  sig_Af2_k'];
    end
    % END OF CURRENT EXPERIMENT/LISTENER
    
    % CONFUSION MATRIX AND CONTRASTS
    exp_stim = exp_params(ind_exp,:);
    % separate trials according to experimental stimuli
    % ind_Ab: trials with auditory only, 'aba'
    % ind_Ad: trials with auditory only, 'ada'
    % ind_Ag: trials with auditory only, 'aga'
    % ind_VbAb: trials with congruent audiovisual 'aba'
    % ind_VgAb: trials with incongruent McGurk stimuli visual 'aga',
    %      acoustic 'aba'
    % ind_VgAg: trials with congruent audiovisual 'aga'
    ind_Ab = intersect(find(exp_stim(:,1)==0),find(exp_stim(:,4)==Af2(1)));
    ind_Ad = intersect(find(exp_stim(:,1)==0),find(exp_stim(:,4)==Af2(2)));
    ind_Ag = intersect(find(exp_stim(:,1)==0),find(exp_stim(:,4)==Af2(3)));
    ind_VbAb = intersect(find(exp_stim(:,1)==1),find(exp_stim(:,3)==Alip(1)));
    ind_VgAg = intersect(find(exp_stim(:,1)==1),find(exp_stim(:,4)==Af2(3)));
    ind_VgAb = intersect(find(exp_stim(:,3)==Alip(3)),find(exp_stim(:,4)==Af2(1)));
    ind_mcgurk = ind_VgAb;
    
    % after running all the trials for the current participant
    % calculate the relevant quantities

        b = b_est;
        % identify the fused McGurk stimili:
        ind_mcgurk_fused = intersect(ind_VgAb,find(b==2));
        % control stimuli for fused McGurk
        ind_mc_other = setdiff([1:length(b)],ind_mcgurk);
        % ------------------------------------------------------
        
        % 
        ind_tem = intersect(ind_Ab-1,ind_mcgurk_fused)+1;
        n1_tem = hist(b(ind_tem),[1:3]);
        [n1_tem, n1_tem(2)/sum(n1_tem)*100];
        N1 = sum(n1_tem);
        n1 = n1_tem(2)/sum(n1_tem)*100;
        
        ind_tem = intersect(ind_Ab-1,ind_mc_other)+1;
        n2_tem = hist(b(ind_tem),[1:3]);
        [n2_tem, n2_tem(2)/sum(n2_tem)*100];
        N2 = sum(n2_tem);
        n2 = n2_tem(2)/sum(n2_tem)*100;
        
        ind_tem = intersect(ind_Ad-1,ind_mcgurk_fused)+1;
        n3_tem = hist(b(ind_tem),[1:3]);
        [n3_tem, n3_tem(2)/sum(n3_tem)*100];
        N3 = sum(n3_tem);
        n3 = n3_tem(2)/sum(n3_tem)*100;
        
        ind_tem = intersect(ind_Ad-1,ind_mc_other)+1;
        n4_tem = hist(b(ind_tem),[1:3]);
        [n4_tem, n4_tem(2)/sum(n4_tem)*100];
        N4 = sum(n4_tem);
        n4 = n4_tem(2)/sum(n4_tem)*100;
        
        % control
        ind_Ab_aba = ind_Ab(b(ind_Ab)==1);
        ind_Ad_ada = ind_Ad(b(ind_Ad)==2);
        ind_Ag_aga = ind_Ag(b(ind_Ag)==3);
        ind_aud_other2 = [ind_Ab_aba; ind_Ag_aga];
        
        ind_tem = intersect(ind_Ab-1,ind_Ad_ada)+1;
        n5_tem = hist(b(ind_tem),[1:3]);
        [n5_tem, n5_tem(2)/sum(n5_tem)*100];
        N5 = sum(n5_tem);
        n5 = n5_tem(2)/sum(n5_tem)*100;
        
        ind_tem = intersect(ind_Ab-1,ind_aud_other2)+1;
        n6_tem = hist(b(ind_tem),[1:3]);
        [n6_tem, n6_tem(2)/sum(n6_tem)*100];
        N6 = sum(n6_tem);
        n6 = n6_tem(2)/sum(n6_tem)*100;
        
        ind_tem = intersect(ind_Ad-1,ind_Ad_ada)+1;
        n7_tem = hist(b(ind_tem),[1:3]);
        [n7_tem, n7_tem(2)/sum(n7_tem)*100];
        N7 = sum(n7_tem);
        n7 = n7_tem(2)/sum(n7_tem)*100;
        
        ind_tem = intersect(ind_Ad-1,ind_aud_other2)+1;
        n8_tem = hist(b(ind_tem),[1:3]);
        [n8_tem, n8_tem(2)/sum(n8_tem)*100];
        N8 = sum(n8_tem);
        n8 = n8_tem(2)/sum(n8_tem)*100;
                
        
    
    n_panel(irep,:) = [N2,N1,N4,N3 N6 N5 N8 N7];
    n2_panel(irep,:) = [n2_tem(2) n1_tem(2) n4_tem(2) n3_tem(2) ...
        n6_tem(2) n5_tem(2) n8_tem(2) n7_tem(2)];

    % confusion matrix
    recog_all(irep,:) = [hist(b(ind_Ab),[1 2 3])    hist(b(ind_Ad),[1 2 3])       hist(b(ind_Ag),[1 2 3]) ...
        hist(b(ind_VbAb),[1 2 3])  hist(b(ind_mcgurk),[1 2 3])   hist(b(ind_VgAg),[1 2 3])];
    
    n_prev(irep) = sum(ind_c_over);
end
% each row contains data for multiple contrasts for a single participant. 
data_panel = n2_panel./n_panel*100;


return
