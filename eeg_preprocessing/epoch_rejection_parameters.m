
% give epoch rejection values/parameters here subjects.
% main_preprocessing will use the same value on all of the sessions of a
% subject. and save it in the end of processing for book keeping.

function [epoch_rej_parameters] = epoch_rejection_parameters(number_of_subjects)


%% parameters for all subjects
% to tune parameters go to EEGLAB gui > load existing dataset >
% Marked Epochs - Fir > -Epoched.set > ERPLAB > Artifact detection in
% epoched data > select parameter > modify parameters values & select Mark Flag 1 +
% # corresponding to selected parameter (see below) > Accept > see Matlab command window output & scroll
% data > once tuned then update parameter values below
for i = 1:number_of_subjects

    epoch_rej_parameters(i) = struct('max_voltage_threshold', [], ...
        'pp_threshold', [], ...
        'step_like',[], ...
        'sample_difference', [], ...
        'flatline', []);
    
    % extreme/max voltage threshold, positive and negative. Mark flag 1 & 2 
    % help pop_artextval
    epoch_rej_parameters(i).max_voltage_threshold.value = 100; % original 100 
    
    % peak-to-peak voltage threshold in a given window. Mark flag 1 & 3.
    % help pop_artmwppth
    epoch_rej_parameters(i).pp_threshold.value = 150 % original 150
    epoch_rej_parameters(i).pp_threshold.window_size = 200;
    epoch_rej_parameters(i).pp_threshold.window_step = 100;
    
    % step-like artifacts. Mark flag 1 & 4.
    % help pop_artstep
    epoch_rej_parameters(i).step_like.value = 100; % original 100
    epoch_rej_parameters(i).step_like.window_size = 200;
    epoch_rej_parameters(i).step_like.window_step = 50;
    
    % sample-to-sample difference artifacts. Mark flag 1 & 5.
    % help pop_artdiff
    epoch_rej_parameters(i).sample_difference.value = 50; % original 50 
    
    % flatline artifacts. Mark flag 1 & 6.
    % help pop_artflatline
    epoch_rej_parameters(i).flatline.value = 1;
    epoch_rej_parameters(i).flatline.duration = 150;
    
end


%% parameters per subject



% %% subject 2
% epoch_rej_parameters(2).max_voltage_threshold.value = 88;
% 
% % peak-to-peak volage threshold in a given window. Marker bits 1
% % and 3.
% % help pop_artmwppth
% epoch_rej_parameters(2).pp_threshold.value = 99;
% epoch_rej_parameters(2).pp_threshold.window_size = 201;
% epoch_rej_parameters(2).pp_threshold.window_step = 101;
% 
% % step-like artifacts. Marker bits 1 and 4.
% % help pop_artstep
% epoch_rej_parameters(2).step_like.value = 100;
% epoch_rej_parameters(2).step_like.window_size = 200;
% epoch_rej_parameters(2).step_like.window_step = 50;
% 
% % sample-to-sample difference artifacts. Marker bits 1 and 5.
% % help pop_artdiff
% epoch_rej_parameters(2).sample_difference.value = 50;
% 
% % flatline artifacts. Marker bits 1 and 6.
% % help pop_artflatline
% epoch_rej_parameters(2).flatline.value = 1;
% epoch_rej_parameters(2).flatline.duration = 150;




% %% subject 6
% epoch_rej_parameters(6).max_voltage_threshold.value = 100;
% 
% % peak-to-peak volage threshold in a given window. Marker bits 1
% % and 3.
% % help pop_artmwppth
% epoch_rej_parameters(6).pp_threshold.value = 150;
% epoch_rej_parameters(6).pp_threshold.window_size = 200;
% epoch_rej_parameters(6).pp_threshold.window_step = 100;
% 
% % step-like artifacts. Marker bits 1 and 4.
% % help pop_artstep
% epoch_rej_parameters(6).step_like.value = 100;
% epoch_rej_parameters(6).step_like.window_size = 200;
% epoch_rej_parameters(6).step_like.window_step = 50;
% 
% % sample-to-sample difference artifacts. Marker bits 1 and 5.
% % help pop_artdiff
% epoch_rej_parameters(6).sample_difference.value = 50;
% 
% % flatline artifacts. Marker bits 1 and 6.
% % help pop_artflatline
% epoch_rej_parameters(6).flatline.value = 1;
% epoch_rej_parameters(6).flatline.duration = 150;

end