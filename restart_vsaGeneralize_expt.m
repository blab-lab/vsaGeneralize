function expPath = restart_vsaGeneralize_expt(expt)
%RESTART_EXPT  Restart experiment after a crash.

if nargin < 1, expt = []; end

if ~isfield(expt,'name'), expt.name = 'vsaGeneralize'; end
if ~isfield(expt,'snum'), expt.snum = get_snum; end

expFun = get_experiment_function(expt.name);

% find all temp trial dirs
subjPath = get_acoustSavePath(expt.name,expt.snum);
tempdirs = regexp(genpath(subjPath),'[^;]*temp_trials','match')';
if isempty(tempdirs)
    fprintf('No unfinished experiments to restart.\n')
    expPath = [];
    return;
end

% prompt for restart
for d = 1:length(tempdirs)
    %find last trial saved
    trialnums = get_sortedTrials(tempdirs{d});
    if trialnums
        lastTrial = trialnums(end);
    else 
        lastTrial = 1;
    end
    
    %check to see if experiment completed. only prompt to rerun if
    %incomplete.
    dataPath = fileparts(strip(tempdirs{d},'right',filesep));
    load(fullfile(dataPath,'expt.mat'),'expt') % get expt file 
    if lastTrial ~= expt.ntrials
        startName = regexp(dataPath,expt.snum);
        expName = dataPath(startName:end);
        q = sprintf('Restart experiment "%s" at trial %d? [y/n] ', expName, lastTrial+1);
        q = strrep(q,'\','\\'); %add extra \ to string to display correctly in "input" command
        response = input(q, 's');
        if strcmp(strip(response),'y')
            % restart expt
            expt.startTrial = lastTrial+1;      % set starting trial
            expt.isRestart = 1;
            if ~isfield(expt,'crashTrials')
                expt.crashTrials = [];
            end
            expt.crashTrials = [expt.crashTrials expt.startTrial];
            save(fullfile(dataPath,'expt.mat'),'expt')            
            
            % set current condition
            conds2run = expt.conds(expt.allConds(expt.startTrial));

            % run
            expFun(expt,conds2run)

            expPath = fileparts(strip(tempdirs{d},'right',filesep));
            fprintf('%s phase completed.\n',conds2run{end})
            break;
        else
            fprintf('Restart canceled.\n')
        end
    end
    fprintf('All %d trials completed.\n',expt.ntrials)
    expPath = [];
end



end

