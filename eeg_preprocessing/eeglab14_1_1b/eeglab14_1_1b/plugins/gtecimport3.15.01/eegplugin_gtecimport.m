% eegplugin_gtecimport() - EEGLAB plugin for importing gtec data files.

function vers = eegplugin_gtecimport(fig, trystrs, catchstrs)

    vers = 'gtecimport 3.15.01';
    
    if nargin < 3
        error('eegplugin_gtecimport requires 3 arguments');
    end;
  
    % add folder to path
    % ------------------
    if ~exist('pop_loadeep')
        p = which('eegplugin_gtecimport.m');
        p = p(1:findstr(p,'eegplugin_gtecimport.m')-1);
        addpath( p );
    end;
    
    % find import data menu
    % ---------------------
    menu = findobj(fig, 'tag', 'import data');

    % menu callbacks
    % --------------
 	cb_read_gtec  = [ trystrs.no_check '[EEG ALLEEG LASTCOM] = pop_read_gtec(ALLEEG);' catchstrs.new_non_empty ]; 
    uimenu( menu, 'label','From g.tec file (Sample Data)', 'CallBack', cb_read_gtec, 'Separator', 'on');
    
 	cb_read_gtecposition  = [ trystrs.no_check '[EEG LASTCOM] = pop_read_gtecposition(EEG);' catchstrs.new_non_empty ]; 
    uimenu( menu, 'label','From g.tec file (Location Data)', 'CallBack', cb_read_gtecposition, 'Separator', 'off');
    
    % not available in eeglab 13.3.2
    %menu = findobj(fig, 'tag', 'Exerp');
    %cb_read_gtecerp  = [ trystrs.no_check '[EEG ALLEEG LASTCOM ERP ALLERP] = pop_read_erpgtec(ALLEEG);' catchstrs.new_non_empty ]; 
    %uimenu( menu, 'label','<html>Import <b>ERP</b> from g.tec file<html>', 'CallBack', cb_read_gtecerp, 'Separator', 'on');