%% MODEL PARAMETERS 
clear all; close all;
% the columns stand for
% sim_s_lip, sim_s_f2, sim_C_lip, sim_C_f2
sig_vals_1 = [0.1,   0.1,   0.1,  0.1 ; ...
            0.12,  0.12,  0.1,  0.1 ;...
            0.15,  0.15,  0.1,  0.1; ...
            0.12,  0.15,  0.1,  0.1; ...
            0.15,  0.12,  0.1,  0.1];
        
sig_vals_2 = [0.1,   0.1,   0.2,  0.1 ; ...
            0.12,  0.12,  0.2,  0.1 ;...
            0.15,  0.15,  0.2,  0.1; ...
            0.12,  0.15,  0.2,  0.1; ...
            0.15,  0.12,  0.2,  0.1];
        
sig_vals_3 = [0.1,   0.1,   0.2,  0.2 ; ...
            0.12,  0.12,  0.2,  0.2 ;...
            0.15,  0.15,  0.2,  0.2; ...
            0.12,  0.15,  0.2,  0.2; ...
            0.15,  0.12,  0.2,  0.2];
        
sig_vals_4= [0.1,   0.1,   0.1,  0.2 ; ...
            0.12,  0.12,  0.1,  0.2 ;...
            0.15,  0.15,  0.1,  0.2; ...
            0.12,  0.15,  0.1,  0.2; ...
            0.15,  0.12,  0.1,  0.2];

        
sig_vals = [sig_vals_1; sig_vals_2; sig_vals_3; sig_vals_4];



% -------------------------------------------------
% update method parameters: choose one
%update_method = 'pure accumulation';  % delta rule
update_method = 'decaying accumulation';
%update_method = 'standard Bayes';
%update_method = 'no updates';

clear up_meth_params
if strcmp(update_method, 'pure accumulation')
    [a_vals, b_vals] = meshgrid([0.1 0.2], [0.1:0.1:0.6]); 
    up_meth_params = [a_vals(1:end)'  b_vals(1:end)'];
    [a_vals, b_vals] = meshgrid([0.4:0.1:0.6]); 
    up_meth_params = [up_meth_params; a_vals(1:end)'  b_vals(1:end)'];
    [a_vals, b_vals] = meshgrid([0.7:0.1:0.8]); 
    up_meth_params = [up_meth_params; a_vals(1:end)'  b_vals(1:end)'];
    [a_vals, b_vals] = meshgrid([0.05 0.1], [0.02:0.02:0.08]); 
    up_meth_params = [up_meth_params; a_vals(1:end)'  b_vals(1:end)'];
elseif strcmp(update_method, 'decaying accumulation')
    % [r_lip_1 r_lip_2 r_f2_1 r_f2_2]
    % up_meth_params =[r_lip_1 r_lip_2 r_f2_1 r_f2_2]; 
    [a_vals, b_vals] = meshgrid([0.05, 0.1:0.02:0.2],[0.05 0.1 0.2 0.3 0.4]);
    up_meth_params = [a_vals(1:end)'  b_vals(1:end)'  a_vals(1:end)'  b_vals(1:end)'];
elseif strcmp(update_method, 'standard Bayes')
    [a_vals, b_vals] = meshgrid([1 5 10]);  
    up_meth_params = [a_vals(1:end)'  b_vals(1:end)'];
elseif strcmp(update_method, 'no updates')
    % no parameters to set
    up_meth_params =[ ];
else
    disp('you need to define an update method')
    return
end

n_up_meth_params = size(up_meth_params,1);
n_sig_vals = size(sig_vals,1);


% combine model and update parameters
[iup isig] = meshgrid([1:n_up_meth_params],[1:n_sig_vals]);
iup = iup(1:end);
isig = isig(1:end);
sig_vals = sig_vals(isig,:);
up_meth_params = up_meth_params(iup,:);

[size(up_meth_params)  size(sig_vals)]


%% RUN THE SIMULATIONS
vals_1 = sig_vals(:,1);
vals_2 = sig_vals(:,2);
vals_3 = sig_vals(:,3);
vals_4 = sig_vals(:,4);


disp(['number of simus: ' num2str(length(vals_1))])

% run the simus 
disp('************       ***************')


% determine the name of the file
% define the output directory
sim_dir = [pwd '/simulation_output']
try
    date_str = char(datetime);
catch
    date_str = datestr(clock);
end

date_str(findstr(date_str,':'))='-';
date_str(findstr(date_str,' '))='-';

date_str


%>>>>>>>>>>>>>>>>>>>>>>>
% Normalized location parameters for categories 
% visual feature:  (/aba/, /ada/, /aga/)
Alip = [1  0.6  0.37]';
% acoustic feature;  (/aba/, /ada/, /aga/)
Af2 = [0.1  0.4  1]';
cr_over = 0; %  0: don't allow cross-overs, 1: allow cross-overs

% file containing the 100 experimental paradigms
% each has 69 repetitions of the 6 stimuli in random order with sensory 
% noise
load experimental_stimuli

tic
for isim = 1:3 %length(vals_1)
    sig_s_lip = vals_1(isim);
    sig_s_f2 = vals_2(isim);
    sig_Alip = vals_3(isim);
    sig_Af2 = vals_4(isim);
    
    % this runs the 100 subjects with the current parameter values
    % it also derives the most relevant quantities
    clear n_prev 
    run_experiment_100_participants;
      
    update_parameters = up_meth_params(isim,:);
    sig_parameters = sig_vals(isim,:);
    % save 
    save([sim_dir '/luttke_params_' date_str '_iparam_' num2str(isim) '.mat'],...
    'data_panel', 'recog_all','params_adapt','params_adapt_sig', ...
    'update_method','sig_parameters','update_parameters');
    disp(['done with parameter set ' num2str(isim) ' of ' num2str(length(vals_1))])
end

disp('ALL DONE')
return

