function expt = run_vsaGeneralize_expt(expt,bTestMode)
%RUN_VSAADAPT_SETUP  Run VSA pilot adaptation experiment.
%   RUN_VSAADAPT_SETUP(EXPT,BTESTMODE)

if nargin < 1, expt = []; end
if nargin < 2 || isempty(bTestMode), bTestMode = 0; end

%% set up stimuli
expt.name = 'vsaGeneralize';
if ~isfield(expt,'snum'), expt.snum = get_snum; end
if ~isfield(expt,'gender'), expt.gender = get_gender; end
expt.dataPath = get_acoustSavePath(expt.name,expt.snum);


% stimuli
expt.conds = {'baselineGeneralize' 'baselineTrain' 'train' 'generalization'};
words2train = {'bead' 'bad' 'booed' 'bod'}; 
% words2generalize = {'bid', 'bayed', 'bed', 'bode'}; 
words2generalize = {'bid', 'bayed', 'bed', 'bode', 'bud'}; 


% timing
expt.timing.stimdur = 1.5;         % time stim is on screen, in seconds
expt.timing.interstimdur = .75;    % minimum time between stims, in seconds
expt.timing.interstimjitter = .75; % maximum extra time between stims (jitter)


nwords2train = length(words2train);
nwords2generalize = length(words2generalize);
nwords = nwords2train + nwords2generalize;

if bTestMode
    testModeReps = 1;
    nPre = 1*nwords2train;
    nBaseline_generalize = testModeReps*nwords2generalize;
    nBaseline_train = testModeReps*nwords2train;
    nTrain = testModeReps*nwords2train;
    nGeneralization = testModeReps*nwords2generalize;
    expt.breakFrequency = testModeReps*nwords;
else
    nPre = 5*nwords2train; 
    nBaseline_generalize = 15*nwords2generalize;
    nBaseline_train =  10*nwords2train;
    nTrain = 70*nwords2train;
    nGeneralization = 10*nwords2generalize;
    expt.breakFrequency = 30;
end

%% set up calibration phase
exptpre = expt;
exptpre.dataPath = fullfile(expt.dataPath,'pre');
exptpre.conds = {'baselineTrain'};
exptpre.ntrials = nPre;
exptpre.words = words2train;

%% set up main experiment
% ntrials
expt.ntrials = nBaseline_generalize + nBaseline_train + nTrain + nGeneralization;
expt.breakTrials = expt.breakFrequency:expt.breakFrequency:expt.ntrials;

% conds
expt.allConds = [1*ones(1,nBaseline_generalize) 2*ones(1,nBaseline_train) 3*ones(1,nTrain) 4*ones(1,nGeneralization)];

groupStart = 1;
groupEnd = ceil(nBaseline_generalize/nwords2generalize);
indStart = -nwords2generalize; indEnd = 0;
for i = groupStart:groupEnd
    if (i == groupStart)
        indStart = indEnd + 1;
    else
        indStart = indStart + nwords2generalize;
    end
    indEnd = indStart + nwords2generalize - 1;
%     indStart = (i-1)*nwords2generalize+1;
%     indEnd = min(i*nwords2generalize,expt.ntrials);
    rp = nwords2train + randperm(nwords2generalize);
    allWords(indStart:indEnd) = rp(1:length(indStart:indEnd));
end
groupStart = groupStart + ceil(nBaseline_generalize/nwords2generalize);
groupEnd = groupEnd + ceil(nBaseline_train/nwords2train);
for i = groupStart:groupEnd
    if (i == groupStart)
        indStart = indEnd + 1;
    else
        indStart = indStart + nwords2train;
    end
    indEnd = indStart + nwords2train - 1;
%     indStart = (i-1)*nwords2train+1;
%     indEnd = min(i*nwords2train,expt.ntrials);
    rp = randperm(nwords2train);
    allWords(indStart:indEnd) = rp(1:length(indStart:indEnd));
end
groupStart = groupStart + ceil(nBaseline_train/nwords2train);
groupEnd = groupEnd + ceil(nTrain/nwords2train);
for i = groupStart:groupEnd
    if (i == groupStart)
        indStart = indEnd + 1;
    else
        indStart = indStart + nwords2train;
    end
    indEnd = indStart + nwords2train - 1;
%     indStart = (i-1)*nwords2train+1;
%     indEnd = min(i*nwords2train,expt.ntrials);
    rp = randperm(nwords2train);
    allWords(indStart:indEnd) = rp(1:length(indStart:indEnd));
end
groupStart = groupStart + ceil(nTrain/nwords2train);
groupEnd = groupEnd + ceil(nGeneralization/nwords2generalize);
for i = groupStart:groupEnd
    if (i == groupStart)
        indStart = indEnd + 1;
    else
        indStart = indStart + nwords2generalize;
    end
    indEnd = indStart + nwords2generalize - 1;
%     indStart = (i-1)*nwords2generalize+1;
%     indEnd = min(i*nwords2generalize,expt.ntrials);
    rp = nwords2train + randperm(nwords2generalize);
    allWords(indStart:indEnd) = rp(1:length(indStart:indEnd));
end
expt.allWords = allWords;

expt.words = [words2train words2generalize];

% shifts
fieldDim = 257;
p.F1Min = 200;
p.F1Max = 1500;
p.F2Min = 500;
p.F2Max = 3500;
p.pertf1 = floor(p.F1Min:(p.F1Max-p.F1Min)/(fieldDim-1):p.F1Max);
p.pertf2 = floor(p.F2Min:(p.F2Max-p.F2Min)/(fieldDim-1):p.F2Max);
% p.pertAmp = zeros(fieldDim,fieldDim); % define dummy pert field
p.pertAmp = zeros(1, fieldDim);
% p.pertPhi = zeros(fieldDim,fieldDim); % define dummy pert field
p.pertPhi = zeros(1, fieldDim);
expt.audapterParams = p;

maxScaleFact = .5;
% shiftScaleFact is a scalar between 0 and maxScaleFact that is used to scale shiftMag matrix
expt.shiftScaleFact = [zeros(1,nBaseline_generalize) zeros(1, nBaseline_train) maxScaleFact*ones(1,nTrain) zeros(1,nGeneralization)];

%% save expt
if ~exist(expt.dataPath,'dir')
    mkdir(expt.dataPath)
end
exptfile = fullfile(expt.dataPath,'expt.mat');
bSave = savecheck(exptfile);
if bSave
    save(exptfile, 'expt')
    fprintf('Saved expt file: %s.\n',exptfile);
end

h_fig = setup_exptFigs;
get_figinds_audapter; % names figs: stim = 1, ctrl = 2, dup = 3;
% h_sub = get_subfigs_audapter(h_fig(ctrl),1);

pronouce_txt = sprintf(['Please read the following words out [y/n]: \n\n' ...
    'bad, bayed, bead, bed \n' ...
    'bid, bod, bode, booed bud\n']);
txtparams.Color = 'White';
txtparams.FontSize = 45;
txtparams.HorizontalAlignment = 'Center';
h_pronouce = draw_exptText(h_fig, 0.5,0.5,pronouce_txt, txtparams);
while (1)
    if strcmp(input(pronouce_txt, 's'), 'y')
        break;
    end
end
delete_exptText(h_fig, h_pronouce);
close(h_fig);


pertFieldOK = 'no';
% measure vowel space
while strcmp(pertFieldOK, 'no')
    % if ~exist(fullfile(exptpre.dataPath,'data.mat'),'file')
    exptpre = run_measureFormants_audapter(exptpre);
    % end
    
    %check LPC order
    check_audapterLPC(exptpre.dataPath)
    hGui = findobj('Tag','check_LPC');
    waitfor(hGui);
    
    load(fullfile(exptpre.dataPath,'nlpc'),'nlpc')
    
    %% calibrate 2D pert field
    fmtMeans = calc_vowelMeans(exptpre.dataPath);
    [p,h_pertField] = calc_pertField('in',fmtMeans,1);
    pertFieldOK = askNChoiceQuestion('Is the perturbation field OK?', {'yes', 'no'});
    try
        close(h_pertField)
    catch
    end
    
    %set lpc order
    p.nLPC = nlpc;
end
expt.audapterParams = add2struct(expt.audapterParams,p);

% resave expt
save(exptfile, 'expt');
fprintf('Saved pertfield to expt file: %s.\n',exptfile);



%% run adaptation experiment

% conds2run = {'baselineGeneralize', 'baselineTrain', 'train', 'generalization'};
% expt = run_vsaGeneralize_audapter(expt,conds2run);

conds2run = {'baselineGeneralize'};
% expt.words = words2generalize;
expt = run_vsaGeneralize_audapter(expt,conds2run);

conds2run = {'baselineTrain', 'train'};
% expt.words = words2train;
expt = run_vsaGeneralize_audapter(expt,conds2run);

conds2run = {'generalization'};
% expt.words = words2generalize;
expt = run_vsaGeneralize_audapter(expt,conds2run);

