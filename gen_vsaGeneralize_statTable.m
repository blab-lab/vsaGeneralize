function gen_vsaGeneralize_statTable

%% get the processed data
dataPath = get_acoustLoadPath('vsaGeneralize');
load(fullfile(dataPath,'processedData_25-75.mat'))

%% generate data table for AVS
nSubs = size(data2plot.vsaGeneralizeAbsNormalized,2);

phases = {'adapt','gen'};
phaseNames = {'tPdistPercentNormalized'  'gPdistPercentNormalized'};
nPhases = length(phases);

stab = cell(1,nSubs);
for s = 1:nSubs
    ptab = cell(1,nPhases);
    for p = 1:nPhases
        phase = phases{p};
        name = phaseNames{p};
        dat.AVS = data2plot.(name)(s);

        fact.phase = phase;
        fact.subj = s;

        ptab{p} = get_datatable(dat,fact);
        clear dat;
    end
    stab{s} = vertcat(ptab{:});
end
datatable = vertcat(stab{:});

%% saving AVS datatable
writetable(datatable,fullfile(dataPath,'datatable_AVS'))
save(fullfile(dataPath,'datatable_AVS.mat'),'datatable')
clear datatable
clear stab
clear ptab
clear dat
clear fact
%% generate data table for vowels for changes from baseline
nSubs = size(data2plot.vsaGeneralizeAbsNormalized,2);

vowels = fieldnames(data2plot.dists{1});
nVowels = length(vowels);

stab = cell(1,nSubs);
for s = 1:nSubs
    vtab = cell(1,nVowels);
    for v = 1:nVowels
        vow = vowels{v};
        %get adaptation magnitude, change in distance from center
        dat.mag = data2plot.dists{s}.(vow);
        %get adaptation angle, relative to baseline
        dat.ang = get_angle([0,0], data2plot.compensations{s}.(vow));
        
        %calculate single value for normalized distance to each vowel
        if any(strcmp(fieldnames(distances_all{s}),vow))
            mags = [data2plot.dists{s}.iy data2plot.dists{s}.ae...
                data2plot.dists{s}.uw data2plot.dists{s}.aa];
            angs = [get_angle([0,0], data2plot.compensations{s}.iy) ...
                get_angle([0,0], data2plot.compensations{s}.ae) ...
                get_angle([0,0], data2plot.compensations{s}.uw) ...
                get_angle([0,0], data2plot.compensations{s}.aa)];
            dat.trainDistMag = sum(1.*mags./distances_all{1}.(vow));
            dat.trainDistAng = circ_mean(angs,1./distances_all{1}.(vow));
        else
            dat.trainDistMag = NaN;
            dat.trainDistAng = NaN;
        end
        fact.vowel = vow;
        fact.vowNUM = v;
        fact.subj = s;

        vtab{v} = get_datatable(dat,fact);
        clear dat;
    end
    stab{s} = vertcat(vtab{:});
end
datatable = vertcat(stab{:});

%% saving
writetable(datatable,fullfile(dataPath,'datatable_vowels'))
save(fullfile(dataPath,'datatable_vowels.mat'),'datatable')
clear datatable
clear stab
clear ptab
clear dat
clear fact

%% generate data table for vowels, absolute magnitudes of distance from center
nSubs = size(data2plot.vsaGeneralizeAbsNormalized,2);

vowels = fieldnames(data2plot.dists{1});
nVowels = length(vowels);
phases = {'Base';'End'};
nPhases = length(phases);

stab = cell(1,nSubs);
for s = 1:nSubs
    vtab = cell(1,nVowels);
    for v = 1:nVowels
        vow = vowels{v};
        ptab = cell(1,nPhases);
        for p = 1:nPhases
            phaseName = strcat('dists',phases{p});

            %get distance from center
            dat.mag = data2plot.(phaseName){s}.(vow);
            fact.phase = phases{p};
            fact.vowel = vow;
            fact.vowNUM = v;
            fact.subj = s;
    
            ptab{p} = get_datatable(dat,fact);
            clear dat;
        end
        vtab{v} = vertcat(ptab{:});
    end
    stab{s} = vertcat(vtab{:});
end
datatable = vertcat(stab{:});

%% saving
writetable(datatable,fullfile(dataPath,'datatable_vowels_byphase'))
save(fullfile(dataPath,'datatable_vowels_byphase.mat'),'datatable')
clear datatable
clear stab
clear ptab
clear dat
clear fact

%% generate data table for clear speech variables
nSubs = size(data2plot.vsaGeneralizeAbsNormalized,2);

vowels = fieldnames(data2plot.dists{1});
nVowels = length(vowels);
phases = {'Base';'End'};
nPhases = length(phases);

stab = cell(1,nSubs);
for s = 1:nSubs
    vtab = cell(1,nVowels);
    for v = 1:nVowels
        vow = vowels{v};
        ptab = cell(1,nPhases);
        for p = 1:nPhases
            phaseNamef0Max = strcat('f0Max',phases{p});
            phaseNamef0Range = strcat('f0Range',phases{p});
            phaseNameInt = strcat('int',phases{p});
            phaseNameDur = strcat('dur',phases{p});

            %get distance from center
            dat.f0Max = data2plot.(phaseNamef0Max){s}.(vow);
            dat.f0Range = data2plot.(phaseNamef0Range){s}.(vow);
            dat.int = data2plot.(phaseNameInt){s}.(vow);
            dat.dur= data2plot.(phaseNameDur){s}.(vow);
            fact.phase = phases{p};
            fact.vowel = vow;
            fact.vowNUM = v;
            fact.subj = s;
    
            ptab{p} = get_datatable(dat,fact);
            clear dat;
        end
        vtab{v} = vertcat(ptab{:});
    end
    stab{s} = vertcat(vtab{:});
end
datatable = vertcat(stab{:});

%% saving
writetable(datatable,fullfile(dataPath,'datatable_clearSpeech'))
save(fullfile(dataPath,'datatable_clearSpeech.mat'),'datatable')
clear datatable
clear stab
clear ptab
clear dat
clear fact


end
function phi = get_angle(loc1,loc2)
    phi = atan((loc2(2)-loc1(2))/(loc2(1)-loc1(1)));
    if (loc2(1)-loc1(1)) < 0
        phi = phi + pi;
    end
end

