

function [EEG, epoch_length] = add_events (EEG, EEG_length)

epoch_length = 1;      %epoch length in seconds

%% remove redudant data in the beginning and end (resting-state data)
% % not required/gives wrong output for already removed redundant data, as
% % EEG.xmax can be 600.002 to begin with, leading to removal of more data
% % than EEG_length = 600.
% if (EEG.xmax > EEG_length)
%     data_start = (EEG.xmax - EEG_length)/2;
%     data_end = data_start + EEG_length;
% else
%     data_start = 0;
%     data_end = floor((EEG.xmax-epoch_length)/epoch_length)*epoch_length + epoch_length;
% end
% 
% EEG = pop_select( EEG,'time',[data_start, data_end] );

%%
EEG.event = struct('type', {}, 'latency',{}, 'urevent',{});
nevents = length(EEG.event);


%% put markers in resting state data

event_type = {'1'};
event_latency = (1/EEG.srate:epoch_length:(1/EEG.srate)+(EEG.xmax)-epoch_length)*EEG.srate;

[EEG.event, nevents] = func_add_events (EEG.event, event_type, event_latency, nevents);

EEG = eeg_checkset( EEG );

end


%%
function [EEG_event, nevents] = func_add_events (EEG_event, event_type, event_latency, nevents)

event_latency = num2cell(event_latency);

event_type = repmat (event_type, size(event_latency));

Urevent = num2cell(nevents+1:nevents+length(event_latency));

[EEG_event(1, nevents+1:nevents+length(event_type)).type] = event_type{:};
[EEG_event(1, nevents+1:nevents+length(event_type)).latency] = event_latency{:};
[EEG_event(1, nevents+1:nevents+length(event_type)).urevent] = Urevent{:};

nevents = length(EEG_event);

end