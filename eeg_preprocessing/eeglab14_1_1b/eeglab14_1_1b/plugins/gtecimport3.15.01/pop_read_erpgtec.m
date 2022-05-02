% gtec import function
function [EEG, ALLEEG, command, ERP, ALLERP] = pop_read_erpgtec(ALLEEG)

ALLERP = buildERPstruct([]);
CURRENTERP = 0;

% import data to eeglab
[EEG, ALLEEG, command] = pop_read_gtec(ALLEEG);

% select continous dataset
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'retrieve',1,'study',0);

% create eventlist
EEG = pop_editeventlist(EEG);

% epoch data according to bin
EEG = pop_epochbin( EEG , [-200.0  800.0],  'none');

% calc averages
ERP = pop_averager( EEG , 'Criterion', 'all', 'DSindex',  1, 'Stdev', 'on', 'Warning', 'on' );

% save created ERPset
ERP.erpname = 'gtec_erpimport';  % name for erpset menu
pop_savemyerp(ERP, 'erpname', ERP.erpname, 'filename', [ERP.erpname '.erp']);

CURRENTERP = CURRENTERP + 1;
ALLERP(CURRENTERP) = ERP;
updatemenuerp(ALLERP);
% plot averages
pop_ploterps(ERP, 1:2, 1:ERP.nchan);
return;
