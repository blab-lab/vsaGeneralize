function [h, data2plot,distances_euclidean,distances_defined] = analyze_vsaGeneralize(dPs,window,bPlot)
if nargin < 1 || isempty(dPs)
    dPs = get_dataPaths_vsaGeneralize;
end
if nargin < 2 || isempty(window)
    window = [25 75];
end
if ~isequal(size([25 50]),[1 2]) || any(window>100) || any (window<1)
    warning('Analysis window must be a 1 x 2 vector of percentages from 1-100. Using default window (25-50%).')
    window = [];
end

if nargin < 3 || isempty(bPlot)
    bPlot = 0;
end
    
exptName = 'vsaGeneralize';

[data2plot] = preprocess(dPs, window);
[distances_euclidean,distances_defined,nearestTrainVowel,distances_all] ...
        = measure_nearest_vowels(data2plot.posLists);
    
saveDir = get_acoustLoadPath('vsaGeneralize');
saveFile = fullfile(saveDir,sprintf('processedData_%.0f-%.0f.mat',window(1),window(2)));
bSave = savecheck(saveFile);
if bSave
    save(saveFile,'data2plot','distances_euclidean','distances_defined','nearestTrainVowel','distances_all')
end

if bPlot
    h = [];
    plotting(h, data2plot,length(dPs),distances_euclidean,distances_defined,nearestTrainVowel);
end

end

%%

function [distances_euclidean,distances_defined,nearestTrainVowel,distances_all] = measure_nearest_vowels(posLists)
    numParticipants = length(posLists);
    numGenVowels    = 5;

    % init distances to nearest vowels
    distances_euclidean = zeros(numParticipants, numGenVowels);
    nearestTrainVowel   = zeros(numParticipants, numGenVowels);
    % init distances to self-defined ones
    distances_defined.ih2iy    = zeros(numParticipants);
    distances_defined.ey2front = zeros(numParticipants);
    distances_defined.eh2ae    = zeros(numParticipants);
    distances_defined.ou2back  = zeros(numParticipants);
    distances_defined.ah2aa    = zeros(numParticipants);
    for p = 1:numParticipants
        for i = 1:numGenVowels
            % calculate distances to nearest vowels
            dists = euclidean(posLists{p}.baselineGen(i,:), posLists{p}.baselineTrain);
            vowInd = find(dists == min(dists));
            nearestTrainVowel(p, i)   = vowInd;
            distances_euclidean(p, i) = dists(vowInd);
            vow = num2vowstr(i+4);
            distances_all{p}.(vow) = dists;
        end
        % calculate distances to self-defined vowels
        distances_defined.ih2iy(p)    = euclidean(posLists{p}.baselineGen(1, :), posLists{p}.baselineTrain(1, :));
        distances_defined.ey2front(p) = euclidean(posLists{p}.baselineGen(2, :), mean(posLists{p}.baselineTrain([1 2], :)));
        distances_defined.eh2ae(p)    = euclidean(posLists{p}.baselineGen(3, :), posLists{p}.baselineTrain(2, :));
        distances_defined.ou2back(p)  = euclidean(posLists{p}.baselineGen(4, :), mean(posLists{p}.baselineTrain([3 4], :)));
        distances_defined.ah2aa(p)    = euclidean(posLists{p}.baselineGen(5, :), posLists{p}.baselineTrain(4, :));
    end
end

function h = plotting(h, data2plot,num,distances_euclidean,distances_defined,nearestTrainVowel)
%     h = plot_vsa(h, data2plot.posLists);
    h = plot_pdist(h, data2plot);
    h = plot_sp_normalized(h,data2plot,nearestTrainVowel,num);
    nwords = 9;
    [h,mag_c]        = plot_compensation(h,data2plot,num,nwords);
    [h,mag_sp]       = plot_scalar_proj(h,data2plot,num,nwords);
    [h,efficient_sp] = plot_efficient_comp(h,mag_c,mag_sp);
    [h,mag_vp]       = plot_vec_proj(h,data2plot,num,nwords);
    [h,efficient_vp] = plot_efficient_comp(h,mag_c,mag_vp);
    
%     h = plot_dist_vs_vp_by_participant(h, data2plot.vector_projs, data2plot.posLists);
    
    h = plot_dist_vs_sp_by_participant(h, data2plot.scalar_projs,data2plot.posLists,nearestTrainVowel,num);
    
%     h = plot_vsa_train_vs_gen(h,data2plot.vsaTraining,data2plot.vsaGeneralize);
    h = plot_pdist_train_vs_gen(h,data2plot.tPdistAbsNormalized,data2plot.gPdistAbsNormalized,...
            data2plot.tPdistPercentNormalized,data2plot.gPdistPercentNormalized);
end

function h = plot_vsa(h, posLists)
    nwords = 9;
    posBase = zeros(nwords, 2); posNonBase = zeros(nwords, 2);
    for p = 1:length(posLists)
        posBase(1:2,:) = posLists{p}.baselineTrain(1:2,:);
        posBase(3,:) = posLists{p}.baselineTrain(4,:);
        posBase(4,:) = posLists{p}.baselineTrain(3,:);
        posBase(5:7,:) = posLists{p}.baselineGen(1:3,:);  
        posBase(8,:) = posLists{p}.baselineGen(5,:);
        posBase(9,:) = posLists{p}.baselineGen(4,:); 
        
        posNonBase(1:2,:) = posLists{p}.train(1:2,:);  
        posNonBase(3,:) = posLists{p}.train(4,:);  
        posNonBase(4,:) = posLists{p}.train(3,:);  
        posNonBase(5:7,:) = posLists{p}.gen(1:3,:);
        posNonBase(8,:) = posLists{p}.gen(5,:);
        posNonBase(9,:) = posLists{p}.gen(4,:);
        
        h(end+1) = figure;
        scatter(posBase(1:4,1), posBase(1:4,2), 'ko');
        hold all
        scatter(posBase(5:9,1), posBase(5:9,2), 'kx');
        scatter(posNonBase(1:4,1), posNonBase(1:4,2), 'ro');
        scatter(posNonBase(5:9,1), posNonBase(5:9,2), 'rx');
        plot([posBase(1:4,1); posBase(1,1)], [posBase(1:4,2); posBase(1,2)], 'k-');
        plot([posBase(5:9,1); posBase(5,1)], [posBase(5:9,2); posBase(5,2)], 'k--');
        plot([posNonBase(1:4,1); posNonBase(1,1)], [posNonBase(1:4,2); posNonBase(1,2)], 'r-');
        plot([posNonBase(5:9,1); posNonBase(5,1)], [posNonBase(5:9,2); posNonBase(5,2)], 'r--');
        legend({'train base' 'generalize base' 'train' 'generalize'});
        xlabel('formant1')
        ylabel('formant2')
        title('vowel space')
        hold off
    end
end

function h = plot_pdist(h,data2plot)
    h(end+1) = figure;
    scatter(ones(1,length(data2plot.tPdistAbsNormalized)), data2plot.tPdistAbsNormalized, 'ko');
    hold on
    scatter(2.*ones(1,length(data2plot.gPdistAbsNormalized)), data2plot.gPdistAbsNormalized, 'ro');
    for i = 1:length(data2plot.tPdistAbsNormalized)
        plot([1 2], [data2plot.tPdistAbsNormalized(i) data2plot.gPdistAbsNormalized(i)], 'b-');
    end
    xx = 0:0.05:3;
    plot(xx, zeros(1,length(xx)), 'k.');
    xlim([0 3])
    xticks([1, 2])
    xticklabels({'train', 'generalize'})
    xlabel('condition');
    ylabel('mean pairwise distances');
    title('normalized by subtraction')
    hold off
    
    h(end+1) = figure;
    scatter(ones(1,length(data2plot.tPdistPercentNormalized)), data2plot.tPdistPercentNormalized, 'ko');
    hold on
    scatter(2.*ones(1,length(data2plot.gPdistPercentNormalized)), data2plot.gPdistPercentNormalized, 'ro');
    for i = 1:length(data2plot.tPdistPercentNormalized)
        plot([1 2], [data2plot.tPdistPercentNormalized(i) data2plot.gPdistPercentNormalized(i)], 'b-');
    end
    xx = 0:0.05:3;
    plot(xx, ones(1,length(xx)), 'k.');
    xlim([0 3])
    xticks([1, 2])
    xticklabels({'train', 'generalize'})
    xlabel('condition');
    ylabel('mean pairwise distance');
    title('normalized by division')
    hold off
    
    [h_t,p_t,ci_t,stats_t] = ttest(data2plot.tPdistPercentNormalized, 1)
    [h_g,p_g,ci_g,stats_g] = ttest(data2plot.gPdistPercentNormalized, 1)
end

function h = plot_sp_normalized(h,data2plot,nearestTrainVowel,num)
    numTrainVowels = 4;
    numGenVowels   = 5;
    % ih, eh, ey, ow, ah
    cals = zeros(num, numGenVowels);
    defs = zeros(num, numGenVowels);
    for p = 1:num
        for i = 1:numGenVowels
            genVow = num2vowstr(numTrainVowels+i); 
            trainVow = num2vowstr(nearestTrainVowel(p,i));
            % using nanmean() bc scalar_projs have badTrails
            cals(p,i) = nanmean(data2plot.scalar_projs{p}.(genVow)) / nanmean(data2plot.scalar_projs{p}.(trainVow));
        end
        defs(p,1) = nanmean(data2plot.scalar_projs{p}.ih) / nanmean(data2plot.scalar_projs{p}.iy);
        defs(p,2) = nanmean(data2plot.scalar_projs{p}.eh) / nanmean([data2plot.scalar_projs{p}.iy data2plot.scalar_projs{p}.ae]);
        defs(p,3) = nanmean(data2plot.scalar_projs{p}.ey) / nanmean(data2plot.scalar_projs{p}.ae);
        defs(p,4) = nanmean(data2plot.scalar_projs{p}.ow) / nanmean([data2plot.scalar_projs{p}.aa data2plot.scalar_projs{p}.uw]);
        defs(p,5) = nanmean(data2plot.scalar_projs{p}.ah) / nanmean(data2plot.scalar_projs{p}.aa);
    end
    h(end+1) = figure;
    xx = 0:0.05:6;
    plot(xx, ones(1,length(xx)), 'k.');
    hold on
    for i = 1:numGenVowels
        scatter(ones(1,num)*i, cals(:,i), 'ko');
        hold on
        scatter(i, nanmean(cals(:,i)), 700, 'r.');
    end
    xlim([0 6])
    xticks([1, 2, 3, 4, 5])
    xticklabels({'ih', 'eh', 'ey', 'ow', 'ah'})
    xlabel('vowel');
    ylabel('normalized scalar projections');
    title('calculated nearest training vowels');
    h(end+1) = figure;
    xx = 0:0.05:6;
    plot(xx, ones(1,length(xx)), 'k.');
    hold on
    for i = 1:numGenVowels
        scatter(ones(1,num)*i, defs(:,i), 'ko');
        hold on
        scatter(i, nanmean(defs(:,i)), 700, 'r.');
    end
    xlim([0 6])
    xticks([1, 2, 3, 4, 5])
    xticklabels({'ih', 'eh', 'ey', 'ow', 'ah'})
    xlabel('vowel');
    ylabel('normalized scalar projections');
    title('self-defined nearest trainig vowel');
    
    [h_cal,p_cal,ci_cal,stats_cal] = ttest(cals(:), 1)
    [h_def,p_def,ci_def,stats_def] = ttest(defs(:), 1)
    [h_cal_mean,p_cal_mean,ci_cal_mean,stats_cal_mean] = ttest(nanmean(cals, 1), 1)
    [h_def_mean,p_def_mean,ci_def_mean,stats_def_mean] = ttest(nanmean(defs, 1), 1)
end

function [h,magnitude] = plot_compensation(h,data2plot,num,nwords)
    magnitude = zeros(1, nwords); 
    compensations.iy = [];  compensations.ae = [];  compensations.aa = [];  compensations.uw = []; 
    compensations.ih = [];  compensations.eh = [];  compensations.ey = [];  compensations.ow = []; compensations.ah = [];
    for i = 1:num
        magnitude = magnitude + ...
            [nanmean(vector_norm(data2plot.compensations{i}.iy)), nanmean(vector_norm(data2plot.compensations{i}.ae)), ...
             nanmean(vector_norm(data2plot.compensations{i}.uw)), nanmean(vector_norm(data2plot.compensations{i}.aa)), ...
             nanmean(vector_norm(data2plot.compensations{i}.ih)), nanmean(vector_norm(data2plot.compensations{i}.eh)), ...
             nanmean(vector_norm(data2plot.compensations{i}.ey)), nanmean(vector_norm(data2plot.compensations{i}.ow)), ...
             nanmean(vector_norm(data2plot.compensations{i}.ah))]; 
        vn = vector_norm(data2plot.compensations{i}.iy); compensations.iy(:, length(compensations.iy)+1:length(compensations.iy)+length(vn)) = vn;
        vn = vector_norm(data2plot.compensations{i}.ae); compensations.ae(:, length(compensations.ae)+1:length(compensations.ae)+length(vn)) = vn;
        vn = vector_norm(data2plot.compensations{i}.aa); compensations.aa(:, length(compensations.aa)+1:length(compensations.aa)+length(vn)) = vn;
        vn = vector_norm(data2plot.compensations{i}.uw); compensations.uw(:, length(compensations.uw)+1:length(compensations.uw)+length(vn)) = vn;
        vn = vector_norm(data2plot.compensations{i}.ih); compensations.ih(:, length(compensations.ih)+1:length(compensations.ih)+length(vn)) = vn;
        vn = vector_norm(data2plot.compensations{i}.eh); compensations.eh(:, length(compensations.eh)+1:length(compensations.eh)+length(vn)) = vn;
        vn = vector_norm(data2plot.compensations{i}.ey); compensations.ey(:, length(compensations.ey)+1:length(compensations.ey)+length(vn)) = vn;
        vn = vector_norm(data2plot.compensations{i}.ow); compensations.ow(:, length(compensations.ow)+1:length(compensations.ow)+length(vn)) = vn; 
        vn = vector_norm(data2plot.compensations{i}.ah); compensations.ah(:, length(compensations.ah)+1:length(compensations.ah)+length(vn)) = vn; 
    end
    h(end+1) = figure;
    magnitude = magnitude ./ num;
    sz = 100;
    scatter(1*ones(1, length(compensations.iy)), compensations.iy, 'ko'); 
    hold all
    scatter(1, magnitude(1), sz, 'bd', 'filled');
    scatter(2*ones(1, length(compensations.ae)), compensations.ae, 'ko'); scatter(2, magnitude(2), sz, 'bd', 'filled')
    scatter(3*ones(1, length(compensations.uw)), compensations.uw, 'ko'); scatter(3, magnitude(3), sz, 'bd', 'filled')
    scatter(4*ones(1, length(compensations.aa)), compensations.aa, 'ko'); scatter(4, magnitude(4), sz, 'bd', 'filled')
    scatter(5*ones(1, length(compensations.ih)), compensations.ih, 'ro'); scatter(5, magnitude(5), sz, 'md', 'filled')
    scatter(6*ones(1, length(compensations.eh)), compensations.eh, 'ro'); scatter(6, magnitude(6), sz, 'md', 'filled')
    scatter(7*ones(1, length(compensations.ey)), compensations.ey, 'ro'); scatter(7, magnitude(7), sz, 'md', 'filled')
    scatter(8*ones(1, length(compensations.ow)), compensations.ow, 'ro'); scatter(8, magnitude(8), sz, 'md', 'filled')
    scatter(9*ones(1, length(compensations.ah)), compensations.ah, 'ro'); scatter(9, magnitude(9), sz, 'md', 'filled')
    xlim([0 10])
    xticks([1, 2, 3, 4, 5, 6, 7, 8, 9])
    xticklabels({'iy', 'ae', 'uw', 'aa', 'ih', 'eh', 'ey', 'ow', 'ah'})
    xlabel('vowel');
    ylabel('compensation (Mels)');
    title('magnitudes of compensations')
    hold off
end

function [h,mag_scalar_proj] = plot_scalar_proj(h,data2plot,num,nwords)
    mag_scalar_proj = zeros(1, nwords);
    scalar_projs.iy = [];  scalar_projs.ae = [];  scalar_projs.aa = [];  scalar_projs.uw = []; 
    scalar_projs.ih = [];  scalar_projs.eh = [];  scalar_projs.ey = [];  scalar_projs.ow = []; scalar_projs.ah = [];
    for i = 1:num
        mag_scalar_proj = mag_scalar_proj + ...
            [nanmean(data2plot.scalar_projs{i}.iy), nanmean(data2plot.scalar_projs{i}.ae), ...
            nanmean(data2plot.scalar_projs{i}.uw), nanmean(data2plot.scalar_projs{i}.aa), ...
            nanmean(data2plot.scalar_projs{i}.ih), nanmean(data2plot.scalar_projs{i}.eh), ...
            nanmean(data2plot.scalar_projs{i}.ey), nanmean(data2plot.scalar_projs{i}.ow), nanmean(data2plot.scalar_projs{i}.ah)];
        scalar_projs.iy(:, length(scalar_projs.iy)+1:length(scalar_projs.iy)+length(data2plot.scalar_projs{i}.iy)) = data2plot.scalar_projs{i}.iy;
        scalar_projs.ae(:, length(scalar_projs.ae)+1:length(scalar_projs.ae)+length(data2plot.scalar_projs{i}.ae)) = data2plot.scalar_projs{i}.ae;
        scalar_projs.aa(:, length(scalar_projs.aa)+1:length(scalar_projs.aa)+length(data2plot.scalar_projs{i}.aa)) = data2plot.scalar_projs{i}.aa;
        scalar_projs.uw(:, length(scalar_projs.uw)+1:length(scalar_projs.uw)+length(data2plot.scalar_projs{i}.uw)) = data2plot.scalar_projs{i}.uw;
        scalar_projs.ih(:, length(scalar_projs.ih)+1:length(scalar_projs.ih)+length(data2plot.scalar_projs{i}.ih)) = data2plot.scalar_projs{i}.ih;
        scalar_projs.eh(:, length(scalar_projs.eh)+1:length(scalar_projs.eh)+length(data2plot.scalar_projs{i}.eh)) = data2plot.scalar_projs{i}.eh;
        scalar_projs.ey(:, length(scalar_projs.ey)+1:length(scalar_projs.ey)+length(data2plot.scalar_projs{i}.ey)) = data2plot.scalar_projs{i}.ey;
        scalar_projs.ow(:, length(scalar_projs.ow)+1:length(scalar_projs.ow)+length(data2plot.scalar_projs{i}.ow)) = data2plot.scalar_projs{i}.ow;
        scalar_projs.ah(:, length(scalar_projs.ah)+1:length(scalar_projs.ah)+length(data2plot.scalar_projs{i}.ah)) = data2plot.scalar_projs{i}.ah;
    end
    
    h(end+1) = figure;
    mag_scalar_proj = mag_scalar_proj ./ num;
    sz = 100;
    scatter(1*ones(1, length(scalar_projs.iy)), scalar_projs.iy, 'ko'); 
    hold all
    scatter(1, mag_scalar_proj(1), sz, 'bd', 'filled');
    scatter(2*ones(1, length(scalar_projs.ae)), scalar_projs.ae, 'ko'); scatter(2, mag_scalar_proj(2), sz, 'bd', 'filled')
    scatter(3*ones(1, length(scalar_projs.uw)), scalar_projs.uw, 'ko'); scatter(3, mag_scalar_proj(3), sz, 'bd', 'filled')
    scatter(4*ones(1, length(scalar_projs.aa)), scalar_projs.aa, 'ko'); scatter(4, mag_scalar_proj(4), sz, 'bd', 'filled')
    scatter(5*ones(1, length(scalar_projs.ih)), scalar_projs.ih, 'ro'); scatter(5, mag_scalar_proj(5), sz, 'md', 'filled')
    scatter(6*ones(1, length(scalar_projs.ey)), scalar_projs.ey, 'ro'); scatter(6, mag_scalar_proj(6), sz, 'md', 'filled')
    scatter(7*ones(1, length(scalar_projs.eh)), scalar_projs.eh, 'ro'); scatter(7, mag_scalar_proj(7), sz, 'md', 'filled')
    scatter(8*ones(1, length(scalar_projs.ow)), scalar_projs.ow, 'ro'); scatter(8, mag_scalar_proj(8), sz, 'md', 'filled')
    scatter(9*ones(1, length(scalar_projs.ah)), scalar_projs.ah, 'ro'); scatter(9, mag_scalar_proj(9), sz, 'md', 'filled')
    xlim([0 10])
    xticks([1, 2, 3, 4, 5, 6, 7, 8, 9])
    xticklabels({'iy', 'ae', 'uw', 'aa', 'ih', 'ey', 'eh', 'ow', 'ah'})
    xlabel('vowel');
    ylabel('scalar projections (Mels)');
    title('magnitudes of dot projections to the direction of perturbations (dot)')
    hold off
end

function [h,effcient_comp] = plot_efficient_comp(h,mag_c,mag_p)
    effcient_comp = mag_p ./ mag_c;
    h(end+1) = figure;
    c = categorical({'iy', 'ae', 'uw', 'aa', 'ih', 'eh', 'ey', 'ow', 'ah'});
    c = reordercats(c, {'iy', 'ae', 'uw', 'aa', 'ih', 'eh', 'ey', 'ow', 'ah'});
    b = bar(c, effcient_comp);
    b.FaceColor = 'flat';
    b.CData(1:4,:) = [.5 0 .5; .5 0 .5; .5 0 .5; .5 0 .5;];
    xlabel('words')
    ylabel('efficient compensation magnitude (Mels)')
    title('magnitudes of efficient compensations')
end

function [h,mag_vec_proj] = plot_vec_proj(h,data2plot,num,nwords)
    mag_vec_proj = zeros(1, nwords);
    vec_projs.iy = [];  vec_projs.ae = [];  vec_projs.aa = [];  vec_projs.uw = []; 
    vec_projs.ih = [];  vec_projs.eh = [];  vec_projs.ey = [];  vec_projs.ow = []; vec_projs.ah = [];
    for i = 1:num
        mag_vec_proj = mag_vec_proj + ...
            [nanmean(data2plot.vector_projs{i}.iy, 2), nanmean(data2plot.vector_projs{i}.ae, 2), ...
             nanmean(data2plot.vector_projs{i}.uw, 2), nanmean(data2plot.vector_projs{i}.aa, 2), ...
             nanmean(data2plot.vector_projs{i}.ih, 2), nanmean(data2plot.vector_projs{i}.eh, 2), ...
             nanmean(data2plot.vector_projs{i}.ey, 2), nanmean(data2plot.vector_projs{i}.ow, 2), nanmean(data2plot.vector_projs{i}.ah, 2)];
    end
    h(end+1) = figure;
    mag_vec_proj = mag_vec_proj ./ num;
    mag_vec_proj = vector_norm(mag_vec_proj);
    c = categorical({'iy', 'ae', 'uw', 'aa', 'ih', 'eh', 'ey', 'ow', 'ah'});
    c = reordercats(c, {'iy', 'ae', 'uw', 'aa', 'ih', 'eh', 'ey', 'ow', 'ah'});
    b = bar(c, mag_vec_proj);
    hold on
    b.FaceColor = 'flat';
    b.CData(1:4,:) = [.5 0 .5; .5 0 .5; .5 0 .5; .5 0 .5; ];
    xlabel('words')
    ylabel('compensation (Mels)')
    title('magnitudes of compensations after projecting to the perturbation')
end

function h = plot_dist_vs_mag(h,dists_eu,dists_def,mag_sp,mag_vp,efficient_sp,efficient_vp)

    h(end+1) = figure;
    scatter([dists_def.ih2iy dists_def.ey2front dists_def.eh2ae dists_def.ou2back dists_def.ah2aa], mag_sp(5:9), 'ro');
    hold on
    dx = 4; dy = 4;
    text(dists_def.ih2iy+dx,    mag_sp(5)-dy, cellstr('ih'));
    text(dists_def.ey2front+dx, mag_sp(6)-dy, cellstr('ey'));
    text(dists_def.eh2ae+dx,    mag_sp(7)-dy, cellstr('eh'));
    text(dists_def.ou2back+dx,  mag_sp(8)-dy, cellstr('ow'));
    text(dists_def.ah2aa+dx,    mag_sp(9)-dy, cellstr('ah'));
    xlabel('distances to the nearest vowel (Mels)')
    ylabel('magnitudes of compensations (Mels)')
    title('distances vs magnitudes of scalar projection (self-defined)')
    hold off

    h(end+1) = figure;
    scatter([dists_def.ih2iy dists_def.ey2front dists_def.eh2ae dists_def.ou2back dists_def.ah2aa], mag_vp(5:9), 'ro');
    hold on
    dx = 4; dy = 4;
    text(dists_def.ih2iy+dx,    mag_vp(5)-dy, cellstr('ih'));
    text(dists_def.ey2front+dx, mag_vp(6)-dy, cellstr('ey'));
    text(dists_def.eh2ae+dx,    mag_vp(7)-dy, cellstr('eh'));
    text(dists_def.ou2back+dx,  mag_vp(8)-dy, cellstr('ow'));
    text(dists_def.ah2aa+dx,    mag_vp(9)-dy, cellstr('ah'));
    xlabel('distances to the nearest vowel (Mels)')
    ylabel('magnitudes of compensations (Mels)')
    title('distances vs magnitudes of vector projection (self-defined)')
    hold off

    h(end+1) = figure;
    scatter([dists_def.ih2iy dists_def.ey2front dists_def.eh2ae dists_def.ou2back dists_def.ah2aa], efficient_sp(5:9), 'ro');
    hold on
    dx = 4; dy = .02;
    text(dists_def.ih2iy+dx,    efficient_sp(5)-dy, cellstr('ih'));
    text(dists_def.ey2front+dx, efficient_sp(6)-dy, cellstr('ey'));
    text(dists_def.eh2ae+dx,    efficient_sp(7)-dy, cellstr('eh'));
    text(dists_def.ou2back+dx,  efficient_sp(8)-dy, cellstr('ow'));
    text(dists_def.ah2aa+dx,    efficient_sp(9)-dy, cellstr('ah'));
    xlabel('distances to the nearest vowel (Mels)')
    ylabel('magnitudes of compensations (Mels)')
    title('distances vs magnitudes of efficient scalar projection (self-defined)')
    hold off

% 
    h(end+1) = figure;
    scatter([dists_def.ih2iy dists_def.ey2front dists_def.eh2ae dists_def.ou2back dists_def.ah2aa], efficient_vp(5:9), 'ro');
    hold on
    dx = 4; dy = .02;
    text(dists_def.ih2iy+dx,    efficient_vp(5)-dy, cellstr('ih'));
    text(dists_def.ey2front+dx, efficient_vp(6)-dy, cellstr('ey'));
    text(dists_def.eh2ae+dx,    efficient_vp(7)-dy, cellstr('eh'));
    text(dists_def.ou2back+dx,  efficient_vp(8)-dy, cellstr('ow'));
    text(dists_def.ah2aa+dx,    efficient_vp(9)-dy, cellstr('ah'));
    xlabel('distances to the nearest vowel (Mels)')
    ylabel('magnitudes of compensations (Mels)')
    title('distances vs magnitudes- of efficent vector projection (self-defined)')
    hold off

    h(end+1) = figure;
    scatter(dists_eu, mag_sp(5:9), 'ro');
    hold on
    dx = 4; dy = 4;
    text(dists_eu(1)+dx, mag_sp(5)-dy, cellstr('ih'));
    text(dists_eu(2)+dx, mag_sp(6)-dy, cellstr('ey'));
    text(dists_eu(3)+dx, mag_sp(7)-dy, cellstr('eh'));
    text(dists_eu(4)+dy, mag_sp(8)-dy, cellstr('ow'));
    text(dists_eu(5)+dy, mag_sp(9)-dy, cellstr('ah'));
    xlabel('distances to the nearest vowel (Mels)')
    ylabel('magnitudes of compensations (Mels)')
    title('distances vs magnitudes of scalar projection (nearest in VSA)')
    hold off
% 
    h(end+1) = figure;
    scatter(dists_eu, mag_vp(5:9), 'ro');
    hold on
    dx = 4; dy = 4;
    text(dists_eu(1)+dx, mag_vp(5)-dy, cellstr('ih'));
    text(dists_eu(2)+dx, mag_vp(6)-dy, cellstr('ey'));
    text(dists_eu(3)+dx, mag_vp(7)-dy, cellstr('eh'));
    text(dists_eu(4)+dy, mag_vp(8)-dy, cellstr('ow'));
    text(dists_eu(5)+dy, mag_vp(9)-dy, cellstr('ah'));
    xlabel('distances to the nearest vowel (Mels)')
    ylabel('magnitudes of compensations (Mels)')
    title('distances vs magnitudes of vector projection (nearest in VSA)')
    hold off

    h(end+1) = figure;
    scatter(dists_eu, efficient_sp(5:9), 'ro');
    hold on
    dx = 4; dy = 0.02;
    text(dists_eu(1)+dx, efficient_sp(5)-dy, cellstr('ih'));
    text(dists_eu(2)+dx, efficient_sp(6)-dy, cellstr('ey'));
    text(dists_eu(3)+dx, efficient_sp(7)-dy, cellstr('eh'));
    text(dists_eu(4)+dy, efficient_sp(8)-dy, cellstr('ow'));
    text(dists_eu(5)+dy, efficient_sp(9)-dy, cellstr('ah'));
    xlabel('distances to the nearest vowel (Mels)')
    ylabel('magnitudes of compensations (Mels)')
    title('distances vs magnitudes of efficient scalar projection (nearest in VSA)')
    hold off

    h(end+1) = figure;
    scatter(dists_eu, efficient_vp(5:9), 'ro');
    hold on
    dx = 4; dy = 0.02;
    text(dists_eu(1)+dx, efficient_vp(5)-dy, cellstr('ih'));
    text(dists_eu(2)+dx, efficient_vp(6)-dy, cellstr('ey'));
    text(dists_eu(3)+dx, efficient_vp(7)-dy, cellstr('eh'));
    text(dists_eu(4)+dy, efficient_vp(8)-dy, cellstr('ow'));
    text(dists_eu(5)+dy, efficient_vp(9)-dy, cellstr('ah'));
    xlabel('distances to the nearest vowel (Mels)')
    ylabel('magnitudes of compensations (Mels)')
    title('distances vs magnitudes of efficient vector projection (nearest in VSA)')
    hold off
end

function h = plot_dist_vs_vp_by_participant(h,vps,posLists)
    NUM = length(vps);
    vp.iy = NaN(2,NUM); vp.ae = NaN(2,NUM); vp.uw = NaN(2,NUM); vp.aa = NaN(2,NUM);
    vp.ih = NaN(2,NUM); vp.eh = NaN(2,NUM); vp.ey = NaN(2,NUM); vp.ow = NaN(2,NUM); vp.ah = NaN(2,NUM);
    pos.iy = NaN(NUM,2); pos.ae = NaN(NUM,2); pos.uw = NaN(NUM,2); pos.aa = NaN(NUM,2);
    pos.ih = NaN(NUM,2); pos.eh = NaN(NUM,2); pos.ey = NaN(NUM,2); pos.ow = NaN(NUM,2); pos.ah = NaN(NUM,2);
    for p = 1:NUM
        vp.iy(:,p) = nanmean(vps{p}.iy,2); vp.ae(:,p) = nanmean(vps{p}.ae,2);
        vp.uw(:,p) = nanmean(vps{p}.uw,2); vp.aa(:,p) = nanmean(vps{p}.aa,2); 
        vp.ih(:,p) = nanmean(vps{p}.ih,2); vp.eh(:,p) = nanmean(vps{p}.ey,2); 
        vp.ey(:,p) = nanmean(vps{p}.ey,2); vp.ow(:,p) = nanmean(vps{p}.ow,2); 
        vp.ah(:,p) = nanmean(vps{p}.ah,2); 
        pos.iy(p,:) = posLists{p}.baselineTrain(1,:); pos.ae(p,:) = posLists{p}.baselineTrain(2,:); 
        pos.uw(p,:) = posLists{p}.baselineTrain(3,:); pos.aa(p,:) = posLists{p}.baselineTrain(4,:);
        pos.ih(p,:) = posLists{p}.baselineGen(1,:);   pos.ey(p,:) = posLists{p}.baselineGen(2,:); 
        pos.eh(p,:) = posLists{p}.baselineGen(3,:);   pos.ow(p,:) = posLists{p}.baselineGen(4,:); 
        pos.ah(p,:) = posLists{p}.baselineGen(5,:);   
    end
    vp_norm.iy = vector_norm(vp.iy); vp_norm.ae = vector_norm(vp.ae);
    vp_norm.uw = vector_norm(vp.uw); vp_norm.aa = vector_norm(vp.aa);
    vp_norm.ih = vector_norm(vp.ih); vp_norm.eh = vector_norm(vp.eh);
    vp_norm.ey = vector_norm(vp.ey); vp_norm.ow = vector_norm(vp.ow);
    vp_norm.ah = vector_norm(vp.ah);
    ih2iy = NaN(1, NUM); eh2ae = NaN(1, NUM); ey2front = NaN(1, NUM);
    ou2back = NaN(1, NUM); ah2aa = NaN(1, NUM);
    for p = 1:NUM
        ih2iy(1,p) = euclidean(pos.ih(p,:), pos.iy(p,:));
        eh2ae(1,p) = euclidean(pos.eh(p,:), pos.ae(p,:));
%         ey2front(1,p) = nanmean([euclidean(pos.ey(p,:),pos.iy(p,:)) euclidean(pos.ey(p,:),pos.ae(p,:))]);
%         ou2back(1,p) = nanmean([euclidean(pos.ow(p,:),pos.uw(p,:)) euclidean(pos.ow(p,:),pos.aa(p,:))]);
        ey2front(1,p) = euclidean(pos.ey(p,:), mean([pos.iy(p,:); pos.ae(p,:)]));
        ou2back(1,p) = euclidean(pos.ey(p,:), mean([pos.uw(p,:); pos.aa(p,:)]));
        ah2aa(1,p) = euclidean(pos.ah(p,:), pos.aa(p,:));
    end
    
    sz = 35;
    
    h(end+1) = figure;
    p1 = scatter(nanmean(ih2iy), nanmean(vp_norm.ih), sz, 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410]);
    hold all
    p2 = scatter(nanmean(ey2front), nanmean(vp_norm.ey), sz, 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
    p3 = scatter(nanmean(eh2ae), nanmean(vp_norm.eh), sz, 'MarkerEdgeColor', [0.9290, 0.6940, 0.1250], 'MarkerFaceColor', [0.9290, 0.6940, 0.1250]);
    p4 = scatter(nanmean(ou2back), nanmean(vp_norm.ow), sz, 'MarkerEdgeColor', [0.4940, 0.1840, 0.5560], 'MarkerFaceColor', [0.4940, 0.1840, 0.5560]);
    p5 = scatter(nanmean(ah2aa), nanmean(vp_norm.ah), sz, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor', [0.4660, 0.6740, 0.1880]);
    p6 = scatter(ih2iy, vp_norm.ih, sz, 'MarkerEdgeColor', [0, 0.4470, 0.7410]);
    p7 = scatter(ey2front, vp_norm.ey, sz, 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980]);
    p8 = scatter(eh2ae, vp_norm.eh, sz, 'MarkerEdgeColor', [0.9290, 0.6940, 0.1250]);
    p9 = scatter(ou2back, vp_norm.ow, sz, 'MarkerEdgeColor', [0.4940, 0.1840, 0.5560]);
    p10 = scatter(ah2aa, vp_norm.ah, sz, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880]);
    legend([p1,p2,p3,p4,p5], {'ih','ey','eh','ow','ah'})
    xlabel('distances to the nearest vowel (Mels)')
    ylabel('magnitudes of projectoins (Mels)')
    title('distances vs vector projection by participant (self-defined)')
end

function h = plot_dist_vs_sp_by_participant(h,sps,posLists,nearestTrainVowel,num)
    numGenVowels   = 5;
    x_cals = zeros(num, numGenVowels);
    x_defs = zeros(num, numGenVowels);
    for p = 1:num
        for i = 1:numGenVowels
            trainVow = nearestTrainVowel(p,i);
            x_cals(p,i) = euclidean(posLists{p}.baselineGen(i,:), posLists{p}.baselineTrain(trainVow, :));
        end
        x_defs(p,1) = euclidean(posLists{p}.baselineGen(1,:), posLists{p}.baselineTrain(1, :) );
        x_defs(p,2) = euclidean(posLists{p}.baselineGen(2,:), mean([posLists{p}.baselineTrain(1,:); posLists{p}.baselineTrain(2,:)], 1) );
        x_defs(p,3) = euclidean(posLists{p}.baselineGen(3,:), posLists{p}.baselineTrain(2, :) );
        x_defs(p,4) = euclidean(posLists{p}.baselineGen(4,:), mean([posLists{p}.baselineTrain(3,:); posLists{p}.baselineTrain(4,:)], 1) );
        x_defs(p,5) = euclidean(posLists{p}.baselineGen(5,:), posLists{p}.baselineTrain(4, :) );
    end
    y = zeros(num, numGenVowels);
    for p = 1:num
        for i = 1:numGenVowels
            genVow = num2vowstr(i);
            y(p,i) = nanmean(sps{p}.(genVow));
        end
    end
    
    sz = 35;
    colors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880];
    h(end+1) = figure;
    for i = 1:numGenVowels
        scatter(x_cals(:,i), y(:,i), sz, 'MarkerEdgeColor', colors(i,:));
        hold on
    end
    legend('ih', 'eh', 'ey', 'ow', 'ah');
    for i = 1:numGenVowels
        scatter(mean(x_cals(:,i)), nanmean(y(:,i)), sz, 'MarkerEdgeColor', colors(i,:), 'MarkerFaceColor', colors(i,:));
        hold on
    end
    xlabel('distances to the nearest vowel (Mels)')
    ylabel('magnitudes of projectoins (Mels)')
    title('distances vs scalar projection by participant (calculated)')
    
    h(end+1) = figure;
    for i = 1:numGenVowels
        scatter(x_defs(:,i), y(:,i), sz, 'MarkerEdgeColor', colors(i,:));
        hold on
    end
    legend('ih', 'eh', 'ey', 'ow', 'ah');
    for i = 1:numGenVowels
        scatter(mean(x_defs(:,i)), nanmean(y(:,i)), sz, 'MarkerEdgeColor', colors(i,:), 'MarkerFaceColor', colors(i,:));
        hold on
    end
    xlabel('distances to the nearest vowel (Mels)')
    ylabel('magnitudes of projectoins (Mels)')
    title('distances vs scalar projection by participant (self-defined)')
end

function h = plot_vsa_train_vs_gen(h,vsaTraining, vsaGeneralize)
    h(end+1) = figure;
    scatter(vsaTraining, vsaGeneralize, 'ro')
    xlabel('vsaTraining (Hz)')
    ylabel('vsaGeneralize (Hz)')
    title('magnitudes of training vs generalization')
end

function h = plot_pdist_train_vs_gen(h,tPdistAbsNormalized,gPdistAbsNormalized,tPdistPercentNormalized,gPdistPercentNormalized)
    h(end+1) = figure;
    scatter( tPdistAbsNormalized, gPdistAbsNormalized, 'ro')
    xlabel('training data pairwise distance differences')
    ylabel('generalization data pairwise distance diffrences')
    title('pairwise distance training vs generalization (differences)')
    [rho_abs, p_abs] = corr(tPdistAbsNormalized', gPdistAbsNormalized', 'type', 'spearman')
    
    h(end+1) = figure;
    scatter( tPdistPercentNormalized, gPdistPercentNormalized, 'ro')
    xlabel('training data pairwise distance percentage')
    ylabel('generalization data pairwise distance percentage')
    title('pairwise distance training vs generalization (percentage)')
    [rho_percent, p_percent] = corr(tPdistPercentNormalized', gPdistPercentNormalized', 'type', 'spearman')
end

%% preprocessing
% vowels: expt.vowels
% data2plot.field(i): data to plot grouped by participant
%
function [data2plot] = preprocess(dPs, window)
    for p = 1:length(dPs)
        dataPath = dPs{p};
        fprintf('Processing participant %d / %d\n',p,length(dPs))
        dataPathPre = fullfile(dataPath,'pre');
        fprintf('\t loading data\n')
        load(fullfile(dataPath, 'data.mat'));
        fprintf('\t loading expt\n')
        load(fullfile(dataPath, 'expt.mat'));
        fprintf('\t loading dataVals\n')
        load(fullfile(dataPath, 'dataVals.mat'));
        ntrials = length(expt.allWords);
        vowels = expt.vowels;

        %calculate pre-experiment vowel centeroid
        preMeans = calc_vowelMeans(dataPathPre);
        vowels = fieldnames(preMeans);
        for v = 1:length(vowels) %convert to mels
            vow = vowels{v};
            preMeans.(vow) = hz2mels(preMeans.(vow));
        end
        fCenField = calc_vowelCentroid(preMeans);
        
        % indexing
        binSizeTrain = 40; binSizeGen = 50;
        inds.bGenIndices = find(expt.allConds==1);
        inds.bGenIndices = inds.bGenIndices(26:75); % only use last 50 trials
        inds.bTrainIndices = find(expt.allConds==2); 
        inds.trainIndices = find(expt.allConds==3);
        inds.genIndices = find(expt.allConds==4);
        % last 40 (#train_vows * 10) words in the training phase
        inds.lastTrain = inds.trainIndices(length(inds.trainIndices)-binSizeTrain+1:length(inds.trainIndices));

        %calculate baseline vowel centeroid
        baseMeans = calc_vowelMeans(dataPath,{'baselineTrain'});
        vowels = fieldnames(baseMeans);
        for v = 1:length(vowels) %convert to mels
            vow = vowels{v};
            baseMeans.(vow) = hz2mels(baseMeans.(vow));
        end
        [fCen] = calc_vowelCentroid(baseMeans);
        
        [f1sIn,f2sIn,f1sOut,f2sOut,ints,f0Maxs,f0Ranges,durs] = setup_formants(dataVals, ntrials,window);
        [vows,vsa,clearVars] = setup_vowels(f1sIn,f2sIn,ints,f0Maxs,f0Ranges,durs,binSizeTrain,binSizeGen,inds,expt.allVowels);
        pertPhi = data(inds.trainIndices(1)).params.pertPhi;
        fLims.f1Min = data(inds.trainIndices(1)).params.f1Min;
        fLims.f1Max = data(inds.trainIndices(1)).params.f1Max;
        fLims.f2Min = data(inds.trainIndices(1)).params.f2Min;
        fLims.f2Max = data(inds.trainIndices(1)).params.f2Max;
        
        [compensations,scalar_projs,vector_projs,dists,distsBase,distsEnd,clearVarsSorted,perts] = measure_changes(vows,clearVars,pertPhi,fLims,fCen);        
        
        vowelList_bt = expt.allVowels(inds.bTrainIndices);
        vowelList_t  = expt.allVowels(inds.lastTrain);
        vowelList_bg = expt.allVowels(inds.bGenIndices);
        vowelList_g  = expt.allVowels(inds.genIndices);
        [btPdist,bgPdist,tPdist,gPdist,baseTrainInds,posList] = measure_pdist(f1sIn,f2sIn,inds,vowelList_bt,vowelList_t,vowelList_bg,vowelList_g,expt.vowels);

        data2plot.btPdists(p) =  btPdist;
        data2plot.bgPdists(p) =  bgPdist;
        data2plot.tPdists(p) = tPdist;
        data2plot.gPdists(p) = gPdist;

        data2plot.tPdistAbsNormalized(p) = tPdist - btPdist;
        data2plot.gPdistAbsNormalized(p) = gPdist - bgPdist;
        data2plot.tPdistPercentNormalized(p) = tPdist/btPdist;
        data2plot.gPdistPercentNormalized(p) = gPdist/bgPdist;
        data2plot.vsaTrainingAbsNormalized(p) = vsa.VSA_train - vsa.VSA_btrain;
        data2plot.vsaGeneralizeAbsNormalized(p) = vsa.VSA_gen - vsa.VSA_bgen;
        data2plot.vsaTrainingPercentNormalized(p) = vsa.VSA_train / vsa.VSA_btrain;
        data2plot.vsaGeneralizePercentNormalized(p) = vsa.VSA_gen / vsa.VSA_bgen;
        data2plot.compensations{p} = compensations;
        data2plot.scalar_projs{p} = scalar_projs;
        data2plot.vector_projs{p} = vector_projs;
        data2plot.dists{p} = dists;
        data2plot.distsEnd{p} = distsEnd;
        data2plot.distsBase{p} = distsBase;
        data2plot.intEnd{p} = clearVarsSorted.int.end;
        data2plot.intBase{p} = clearVarsSorted.int.base;
        data2plot.durEnd{p} = clearVarsSorted.dur.end;
        data2plot.durBase{p} = clearVarsSorted.dur.base;
        data2plot.f0MaxEnd{p} = clearVarsSorted.f0Max.end;
        data2plot.f0MaxBase{p} = clearVarsSorted.f0Max.base;
        data2plot.f0RangeEnd{p} = clearVarsSorted.f0Range.end;
        data2plot.f0RangeBase{p} = clearVarsSorted.f0Range.base;
        data2plot.posLists{p} = posList;
        data2plot.vows{p} = vows;
        data2plot.fCen{p} = fCen;
        data2plot.fCenField{p} = fCenField;
        data2plot.perts{p} = perts;
        if p == 1
            data2plot.VSA = NaN(length(dPs),length(vsa.VSA));
            data2plot.AVS = NaN(length(dPs),length(vsa.VSA));
        end
        data2plot.VSA(p,:) = vsa.VSA;
        data2plot.AVS(p,:) = vsa.AVS;
        
    end
end

% mean formants between 0.25T to 0.75T
function [f1sIn,f2sIn,f1sOut,f2sOut,ints,f0Maxs,f0Ranges,durs] = setup_formants(dataVals, ntrials,window)
    f1sIn = NaN(1,ntrials); f2sIn = NaN(1,ntrials);
    f1sOut = NaN(1,ntrials); f2sOut = NaN(1,ntrials);
    ints = NaN(1,ntrials); 
    f0Maxs = NaN(1,ntrials);
    f0Ranges = NaN(1,ntrials);
    durs = NaN(1,ntrials);
    for j = 1:ntrials
        % skip bad trials
        if isempty(dataVals(j).f1) || isempty(dataVals(j).f2)
            f1sIn(j) = NaN; f1sOut = NaN;
            f2sIn(j) = NaN; f2sOut = NaN;
            ints(j) = NaN;
            f0Maxs(j) = NaN;
            f0Ranges(j) = NaN;
            durs(j) = NaN;
            continue;
        end
        ftrackInSamps1 = find(~isnan(dataVals(j).f1));
        ftrackIn1 = dataVals(j).f1(ftrackInSamps1);
        ftrackInSamps2 = find(~isnan(dataVals(j).f2));
        ftrackIn2 = dataVals(j).f2(ftrackInSamps2);
        ftrackOutSamps1 = find(~isnan(dataVals(j).f2)); 
        ftrackOut1 = dataVals(j).f1(ftrackOutSamps1);
        ftrackOutSamps2 = find(~isnan(dataVals(j).f2)); 
        ftrackOut2 = dataVals(j).f2(ftrackOutSamps2);
        ftrackLength = length(ftrackIn1); % TODO
        if dataVals(j).vowel == 2 %only use specified analysis window for /ae/
            p1 = round(ftrackLength.*window(1)./100);
            p2 = round(ftrackLength.*window(2)./100);
        else %use 25-75% window for all other vowels
            p1 = round(ftrackLength.*25./100);
            p2 = round(ftrackLength.*75./100);
        end
        f1sIn(j) = hz2mels(mean(ftrackIn1(p1:p2)));
        f2sIn(j) = hz2mels(mean(ftrackIn2(p1:p2)));
        f1sOut(j) = hz2mels(mean(ftrackOut1(p1:p2))); 
        f2sOut(j) = hz2mels(mean(ftrackOut2(p1:p2)));

        ints(j) = max(dataVals(j).int);
        f0Maxs(j) = max(dataVals(j).f0);
        f0Ranges(j) = range(dataVals(j).f0);
        durs(j) = dataVals(j).dur;
    end
end

% vsa: stores all VSA of a participant 
% vsa for training vowels       - polyarea(i, a, u)
% vsa for generalization vowels - poluarea(ih, ey, eh, ah, ow) (could be
% NaN)
% vsa.VSA_bgen, vsa.VSA_btrain, vsa.VSA_train, vsa.VSA_gen
% vows: \in 2 x #total_trails; Formants In; For training vowels, it only
% contains the last 40 trials
function [vows,vsa,clearVars] = setup_vowels(f1sIn,f2sIn,ints,f0Maxs,f0Ranges,durs,binSizeTrain,binSizeGen,inds,allVowels)
    vsa.VSA_bgen = NaN(1,1);
    vsa.VSA_btrain = NaN(1,1);
    vsa.VSA_train = NaN(1,1);
    vsa.VSA_gen = NaN(1,1);
    vsa.VSA = [vsa.VSA_bgen vsa.VSA_btrain vsa.VSA_train vsa.VSA_gen];
    vsa.AVS = vsa.VSA;

    iVSA = 0; 

    % baseline generalize
    iVSA = iVSA+1;
    currTrials = inds.bGenIndices;
    vowelList = allVowels(currTrials);
    ihLoc = currTrials(vowelList == 5); % bid
    ehLoc = currTrials(vowelList == 7); % bed
    eyLoc = currTrials(vowelList == 6); % bayed
    owLoc = currTrials(vowelList == 8); % bode
    ahLoc = currTrials(vowelList == 9); % bud
    
    vows.ih_bgen(1) = nanmean(f1sIn(ihLoc)); vows.ih_bgen(2) = nanmean(f2sIn(ihLoc)); 
    vows.eh_bgen(1) = nanmean(f1sIn(ehLoc)); vows.eh_bgen(2) = nanmean(f2sIn(ehLoc)); 
    vows.ey_bgen(1) = nanmean(f1sIn(eyLoc)); vows.ey_bgen(2) = nanmean(f2sIn(eyLoc)); 
    vows.ow_bgen(1) = nanmean(f1sIn(owLoc)); vows.ow_bgen(2) = nanmean(f2sIn(owLoc)); 
    vows.ah_bgen(1) = nanmean(f1sIn(ahLoc)); vows.ah_bgen(2) = nanmean(f2sIn(ahLoc));
    
    clearVars.int.ih_bgen = nanmean(ints(ihLoc));
    clearVars.int.eh_bgen = nanmean(ints(ehLoc));
    clearVars.int.ey_bgen = nanmean(ints(eyLoc));
    clearVars.int.ow_bgen = nanmean(ints(owLoc));
    clearVars.int.ah_bgen = nanmean(ints(ahLoc));

    clearVars.f0Max.ih_bgen = nanmean(f0Maxs(ihLoc));
    clearVars.f0Max.eh_bgen = nanmean(f0Maxs(ehLoc));
    clearVars.f0Max.ey_bgen = nanmean(f0Maxs(eyLoc));
    clearVars.f0Max.ow_bgen = nanmean(f0Maxs(owLoc));
    clearVars.f0Max.ah_bgen = nanmean(f0Maxs(ahLoc));

    clearVars.f0Range.ih_bgen = nanmean(f0Ranges(ihLoc));
    clearVars.f0Range.eh_bgen = nanmean(f0Ranges(ehLoc));
    clearVars.f0Range.ey_bgen = nanmean(f0Ranges(eyLoc));
    clearVars.f0Range.ow_bgen = nanmean(f0Ranges(owLoc));
    clearVars.f0Range.ah_bgen = nanmean(f0Ranges(ahLoc));

    clearVars.dur.ih_bgen = nanmean(durs(ihLoc));
    clearVars.dur.eh_bgen = nanmean(durs(ehLoc));
    clearVars.dur.ey_bgen = nanmean(durs(eyLoc));
    clearVars.dur.ow_bgen = nanmean(durs(owLoc));
    clearVars.dur.ah_bgen = nanmean(durs(ahLoc));
    
    vsa.VSA_bgen = polyarea([vows.ih_bgen(1) vows.ey_bgen(1)...
        vows.eh_bgen(1)  vows.ah_bgen(1) vows.ow_bgen(1)], ...
        [vows.ih_bgen(2) vows.ey_bgen(2)...
        vows.eh_bgen(2)  vows.ah_bgen(2) vows.ow_bgen(2)]);
    
    vsa.VSA(iVSA) = NaN;

    %calculate AVS
    posList = [vows.ih_bgen(1) vows.ih_bgen(2); ...
        vows.ey_bgen(1) vows.ey_bgen(2); ...
        vows.eh_bgen(1) vows.eh_bgen(2); ...
        vows.ah_bgen(1) vows.ah_bgen(2); ...
        vows.ow_bgen(1) vows.ow_bgen(2)];
    vsa.AVS(iVSA) = mean(pdist(posList));
    
    % baseline training
    iVSA = iVSA+1;
    currTrials = inds.bTrainIndices;
    vowelList = allVowels(currTrials);
    uLoc = currTrials(vowelList==3); 
    iLoc = currTrials(vowelList==1);
    aLoc = currTrials(vowelList==4); 
    aeLoc = currTrials(vowelList==2); 
    
    vows.u_btrain(1) = nanmean(f1sIn(uLoc));   vows.u_btrain(2) = nanmean(f2sIn(uLoc)); 
    vows.i_btrain(1) = nanmean(f1sIn(iLoc));   vows.i_btrain(2) = nanmean(f2sIn(iLoc)); 
    vows.a_btrain(1) = nanmean(f1sIn(aLoc));   vows.a_btrain(2) = nanmean(f2sIn(aLoc)); 
    vows.ae_btrain(1) = nanmean(f1sIn(aeLoc)); vows.ae_btrain(2) = nanmean(f2sIn(aeLoc)); 
    
    clearVars.int.u_btrain = nanmean(ints(uLoc));
    clearVars.int.i_btrain = nanmean(ints(iLoc));
    clearVars.int.a_btrain = nanmean(ints(aLoc));
    clearVars.int.ae_btrain = nanmean(ints(aeLoc));

    clearVars.dur.u_btrain = nanmean(durs(uLoc));
    clearVars.dur.i_btrain = nanmean(durs(iLoc));
    clearVars.dur.a_btrain = nanmean(durs(aLoc));
    clearVars.dur.ae_btrain = nanmean(durs(aeLoc));
    
    clearVars.f0Max.u_btrain = nanmean(f0Maxs(uLoc));
    clearVars.f0Max.i_btrain = nanmean(f0Maxs(iLoc));
    clearVars.f0Max.a_btrain = nanmean(f0Maxs(aLoc));
    clearVars.f0Max.ae_btrain = nanmean(f0Maxs(aeLoc));

    clearVars.f0Range.u_btrain = nanmean(f0Ranges(uLoc));
    clearVars.f0Range.i_btrain = nanmean(f0Ranges(iLoc));
    clearVars.f0Range.a_btrain = nanmean(f0Ranges(aLoc));
    clearVars.f0Range.ae_btrain = nanmean(f0Ranges(aeLoc));

    vsa.VSA_btrain = polyarea([vows.i_btrain(1) vows.ae_btrain(1)...
        vows.a_btrain(1) vows.u_btrain(1)], ...
        [vows.i_btrain(2) vows.ae_btrain(2)...
        vows.a_btrain(2) vows.u_btrain(2)]);
    
    vsa.VSA(iVSA) = vsa.VSA_btrain;

    %calculate AVS
    posList = [vows.i_btrain(1) vows.i_btrain(2);  ...
        vows.ae_btrain(1) vows.ae_btrain(2);...
        vows.u_btrain(1) vows.u_btrain(2);...
        vows.a_btrain(1) vows.a_btrain(2)];
    vsa.AVS(iVSA) = mean(pdist(posList));
    
    % training phase
    for j = inds.trainIndices(1):binSizeTrain:inds.trainIndices(end)
        iVSA = iVSA+1;
        currTrials = j:j+binSizeTrain-1;
        vowelList = allVowels(currTrials);
        uLoc = currTrials(vowelList==3); iLoc = currTrials(vowelList==1);
        aLoc = currTrials(vowelList==4); aeLoc = currTrials(vowelList==2); 
        
        vows.u_train(1) = nanmean(f1sIn(uLoc));   vows.u_train(2) = nanmean(f2sIn(uLoc)); 
        vows.i_train(1) = nanmean(f1sIn(iLoc));   vows.i_train(2) = nanmean(f2sIn(iLoc)); 
        vows.a_train(1) = nanmean(f1sIn(aLoc));   vows.a_train(2) = nanmean(f2sIn(aLoc)); 
        vows.ae_train(1) = nanmean(f1sIn(aeLoc)); vows.ae_train(2) = nanmean(f2sIn(aeLoc)); 
            
        clearVars.int.u_train = nanmean(ints(uLoc));
        clearVars.int.i_train = nanmean(ints(iLoc));
        clearVars.int.a_train = nanmean(ints(aLoc));
        clearVars.int.ae_train = nanmean(ints(aeLoc));
    
        clearVars.dur.u_train = nanmean(durs(uLoc));
        clearVars.dur.i_train = nanmean(durs(iLoc));
        clearVars.dur.a_train = nanmean(durs(aLoc));
        clearVars.dur.ae_train = nanmean(durs(aeLoc));
        
        clearVars.f0Max.u_train = nanmean(f0Maxs(uLoc));
        clearVars.f0Max.i_train = nanmean(f0Maxs(iLoc));
        clearVars.f0Max.a_train = nanmean(f0Maxs(aLoc));
        clearVars.f0Max.ae_train = nanmean(f0Maxs(aeLoc));
    
        clearVars.f0Range.u_train = nanmean(f0Ranges(uLoc));
        clearVars.f0Range.i_train = nanmean(f0Ranges(iLoc));
        clearVars.f0Range.a_train = nanmean(f0Ranges(aLoc));
        clearVars.f0Range.ae_train = nanmean(f0Ranges(aeLoc));
        
        vsa.VSA_train = polyarea([vows.i_train(1) vows.ae_train(1)...
        vows.a_train(1) vows.u_train(1)], ...
        [vows.i_train(2) vows.ae_train(2)...
        vows.a_train(2) vows.u_train(2)]);
        
        vsa.VSA(iVSA) = vsa.VSA_train;
        
        %calculate AVS
        posList = [vows.i_train(1) vows.i_train(2);  ...
            vows.ae_train(1) vows.ae_train(2);...
            vows.u_train(1) vows.u_train(2);...
            vows.a_train(1) vows.a_train(2)];
        vsa.AVS(iVSA) = mean(pdist(posList));
    end
    
    % generalization phase
    iVSA = iVSA+1;
    currTrials = inds.genIndices;
    vowelList = allVowels(currTrials);
    ihLoc = currTrials(vowelList == 5); % bid
    ehLoc = currTrials(vowelList == 7); % bed
    eyLoc = currTrials(vowelList == 6); % bayed
    owLoc = currTrials(vowelList == 8); % bode
    ahLoc = currTrials(vowelList == 9); % bud
    
    vows.ih_gen(1) = nanmean(f1sIn(ihLoc)); vows.ih_gen(2) = nanmean(f2sIn(ihLoc));
    vows.eh_gen(1) = nanmean(f1sIn(ehLoc)); vows.eh_gen(2) = nanmean(f2sIn(ehLoc));
    vows.ey_gen(1) = nanmean(f1sIn(eyLoc)); vows.ey_gen(2) = nanmean(f2sIn(eyLoc));
    vows.ow_gen(1) = nanmean(f1sIn(owLoc)); vows.ow_gen(2) = nanmean(f2sIn(owLoc));
    vows.ah_gen(1) = nanmean(f1sIn(ahLoc)); vows.ah_gen(2) = nanmean(f2sIn(ahLoc));
    
    clearVars.int.ih_gen = nanmean(ints(ihLoc));
    clearVars.int.eh_gen = nanmean(ints(ehLoc));
    clearVars.int.ey_gen = nanmean(ints(eyLoc));
    clearVars.int.ow_gen = nanmean(ints(owLoc));
    clearVars.int.ah_gen = nanmean(ints(ahLoc));

    clearVars.f0Max.ih_gen = nanmean(f0Maxs(ihLoc));
    clearVars.f0Max.eh_gen = nanmean(f0Maxs(ehLoc));
    clearVars.f0Max.ey_gen = nanmean(f0Maxs(eyLoc));
    clearVars.f0Max.ow_gen = nanmean(f0Maxs(owLoc));
    clearVars.f0Max.ah_gen = nanmean(f0Maxs(ahLoc));

    clearVars.f0Range.ih_gen = nanmean(f0Ranges(ihLoc));
    clearVars.f0Range.eh_gen = nanmean(f0Ranges(ehLoc));
    clearVars.f0Range.ey_gen = nanmean(f0Ranges(eyLoc));
    clearVars.f0Range.ow_gen = nanmean(f0Ranges(owLoc));
    clearVars.f0Range.ah_gen = nanmean(f0Ranges(ahLoc));

    clearVars.dur.ih_gen = nanmean(durs(ihLoc));
    clearVars.dur.eh_gen = nanmean(durs(ehLoc));
    clearVars.dur.ey_gen = nanmean(durs(eyLoc));
    clearVars.dur.ow_gen = nanmean(durs(owLoc));
    clearVars.dur.ah_gen = nanmean(durs(ahLoc));

    vsa.VSA_gen = polyarea([vows.ih_gen(1) vows.ey_gen(1) vows.eh_gen(1)...
        vows.ah_gen(1) vows.ow_gen(1)], ...
        [vows.ih_gen(2) vows.ey_gen(2) vows.eh_gen(2)...
        vows.ah_gen(2) vows.ow_gen(2)]);
    vsa.VSA(iVSA) = NaN;
    
    %calculate AVS
    posList = [vows.ih_gen(1) vows.ih_gen(2); ...
        vows.ey_gen(1) vows.ey_gen(2); ...
        vows.eh_gen(1) vows.eh_gen(2); ...
        vows.ah_gen(1) vows.ah_gen(2); ...
        vows.ow_gen(1) vows.ow_gen(2)];
    vsa.AVS(iVSA) = mean(pdist(posList));
    
end

% all in Mels; perturbations are calculated per vowel per individual
% compensations: train - baseline OR generalization - baseline;
% compensations.vow \in 2*10 (first two formants * binSize)
% scalar_projs: call helper function scalar_proj(compensation, getPertPhi)
% scalar_projs.vow \in 1*10
% vector_projs: call helper function vector_proj(compensation, getPertPhi)
% vector_projs.vow \in 2*10
function [compensations,scalar_projs,vector_projs,dists,distsBase,distsEnd,clearVarsSorted,perts] = measure_changes(vows,clearVars,pertPhi,fLims,fCen)
    % setup variables
    u_train = vows.u_train;   u_btrain = vows.u_btrain;
    i_train = vows.i_train;   i_btrain = vows.i_btrain;
    a_train = vows.a_train;   a_btrain = vows.a_btrain;
    ae_train = vows.ae_train; ae_btrain = vows.ae_btrain;
    ih_gen = vows.ih_gen;     ih_bgen = vows.ih_bgen;
    eh_gen = vows.eh_gen;     eh_bgen = vows.eh_bgen;
    ey_gen = vows.ey_gen;     ey_bgen = vows.ey_bgen;
    ow_gen = vows.ow_gen;     ow_bgen = vows.ow_bgen;
    ah_gen = vows.ah_gen;     ah_bgen = vows.ah_bgen;
    
    compensations.('iy') = i_train - i_btrain;
    compensations.('ae') = ae_train - ae_btrain;
    compensations.('uw') = u_train - u_btrain;
    compensations.('aa') = a_train - a_btrain;
    compensations.('ih') = ih_gen - ih_bgen;
    compensations.('eh') = eh_gen - eh_bgen;
    compensations.('ey') = ey_gen - ey_bgen;
    compensations.('ow') = ow_gen - ow_bgen;
    compensations.('ah') = ah_gen - ah_bgen;
    
    %TODO: check getPertPhi
    perts.('iy') = getPertPhi(i_btrain', pertPhi, fLims);
    perts.('ae') = getPertPhi(ae_btrain', pertPhi, fLims);
    perts.('uw') = getPertPhi(u_btrain', pertPhi, fLims);
    perts.('aa') = getPertPhi(a_btrain', pertPhi, fLims); 
    perts.('ih') = getPertPhi(ih_bgen', pertPhi, fLims);
    perts.('eh') = getPertPhi(eh_bgen', pertPhi, fLims);
    perts.('ey') = getPertPhi(ey_bgen', pertPhi, fLims);
    perts.('ow') = getPertPhi(ow_bgen', pertPhi, fLims);
    perts.('ah') = getPertPhi(ah_bgen', pertPhi, fLims);
    
    scalar_projs.('iy') = scalar_proj(compensations.('iy'), getPertPhi(i_btrain', pertPhi, fLims));
    scalar_projs.('ae') = scalar_proj(compensations.('ae'), getPertPhi(ae_btrain', pertPhi, fLims));
    scalar_projs.('uw') = scalar_proj(compensations.('uw'), getPertPhi(u_btrain', pertPhi, fLims));
    scalar_projs.('aa') = scalar_proj(compensations.('aa'), getPertPhi(a_btrain', pertPhi, fLims)); 
    scalar_projs.('ih') = scalar_proj(compensations.('ih'), getPertPhi(ih_bgen', pertPhi, fLims));
    scalar_projs.('eh') = scalar_proj(compensations.('eh'), getPertPhi(eh_bgen', pertPhi, fLims));
    scalar_projs.('ey') = scalar_proj(compensations.('ey'), getPertPhi(ey_bgen', pertPhi, fLims));
    scalar_projs.('ow') = scalar_proj(compensations.('ow'), getPertPhi(ow_bgen', pertPhi, fLims));
    scalar_projs.('ah') = scalar_proj(compensations.('ah'), getPertPhi(ah_bgen', pertPhi, fLims));
    
    vector_projs.('iy') = vector_proj(compensations.('iy'), getPertPhi(i_btrain', pertPhi, fLims));
    vector_projs.('ae') = vector_proj(compensations.('ae'), getPertPhi(ae_btrain', pertPhi, fLims));
    vector_projs.('uw') = vector_proj(compensations.('uw'), getPertPhi(u_btrain', pertPhi, fLims));
    vector_projs.('aa') = vector_proj(compensations.('aa'), getPertPhi(a_btrain', pertPhi, fLims)); 
    vector_projs.('ih') = vector_proj(compensations.('ih'), getPertPhi(ih_bgen', pertPhi, fLims));
    vector_projs.('eh') = vector_proj(compensations.('eh'), getPertPhi(eh_bgen', pertPhi, fLims));
    vector_projs.('ey') = vector_proj(compensations.('ey'), getPertPhi(ey_bgen', pertPhi, fLims));
    vector_projs.('ow') = vector_proj(compensations.('ow'), getPertPhi(ow_bgen', pertPhi, fLims));
    vector_projs.('ah') = vector_proj(compensations.('ah'), getPertPhi(ah_bgen', pertPhi, fLims));
    
    distsBase.('iy') = euclidean(fCen, i_btrain);
    distsBase.('ae') = euclidean(fCen, ae_btrain);
    distsBase.('uw') = euclidean(fCen, u_btrain);
    distsBase.('aa') = euclidean(fCen, a_btrain);
    distsBase.('ih') = euclidean(fCen, ih_bgen);
    distsBase.('eh') = euclidean(fCen, eh_bgen);
    distsBase.('ey') = euclidean(fCen, ey_bgen);
    distsBase.('ow') = euclidean(fCen, ow_bgen);
    distsBase.('ah') = euclidean(fCen, ah_bgen);

    distsEnd.('iy') = euclidean(fCen, i_train);
    distsEnd.('ae') = euclidean(fCen, ae_train);
    distsEnd.('uw') = euclidean(fCen, u_train);
    distsEnd.('aa') = euclidean(fCen, a_train);
    distsEnd.('ih') = euclidean(fCen, ih_gen);
    distsEnd.('eh') = euclidean(fCen, eh_gen);
    distsEnd.('ey') = euclidean(fCen, ey_gen);
    distsEnd.('ow') = euclidean(fCen, ow_gen);
    distsEnd.('ah') = euclidean(fCen, ah_gen);

    dists.('iy') = euclidean(fCen, i_train) - euclidean(fCen, i_btrain);
    dists.('ae') = euclidean(fCen, ae_train) - euclidean(fCen, ae_btrain);
    dists.('uw') = euclidean(fCen, u_train) - euclidean(fCen, u_btrain);
    dists.('aa') = euclidean(fCen, a_train) - euclidean(fCen, a_btrain);
    dists.('ih') = euclidean(fCen, ih_gen) - euclidean(fCen, ih_bgen);
    dists.('eh') = euclidean(fCen, eh_gen) - euclidean(fCen, eh_bgen);
    dists.('ey') = euclidean(fCen, ey_gen) - euclidean(fCen, ey_bgen);
    dists.('ow') = euclidean(fCen, ow_gen) - euclidean(fCen, ow_bgen);
    dists.('ah') = euclidean(fCen, ah_gen) - euclidean(fCen, ah_bgen);

    varNames = fieldnames(clearVars);
    nVars = length(varNames);
    for v = 1:nVars
        var = varNames{v};

        % setup variables
        u_train = clearVars.(var).u_train;   u_btrain = clearVars.(var).u_btrain;
        i_train = clearVars.(var).i_train;   i_btrain = clearVars.(var).i_btrain;
        a_train = clearVars.(var).a_train;   a_btrain = clearVars.(var).a_btrain;
        ae_train = clearVars.(var).ae_train; ae_btrain = clearVars.(var).ae_btrain;
        ih_gen = clearVars.(var).ih_gen;     ih_bgen = clearVars.(var).ih_bgen;
        eh_gen = clearVars.(var).eh_gen;     eh_bgen = clearVars.(var).eh_bgen;
        ey_gen = clearVars.(var).ey_gen;     ey_bgen = clearVars.(var).ey_bgen;
        ow_gen = clearVars.(var).ow_gen;     ow_bgen = clearVars.(var).ow_bgen;
        ah_gen = clearVars.(var).ah_gen;     ah_bgen = clearVars.(var).ah_bgen;

        clearVarsSorted.(var).base.('iy') = i_btrain;
        clearVarsSorted.(var).base.('ae') = ae_btrain;
        clearVarsSorted.(var).base.('uw') = u_btrain;
        clearVarsSorted.(var).base.('aa') = a_btrain;
        clearVarsSorted.(var).base.('ih') = ih_bgen;
        clearVarsSorted.(var).base.('eh') = eh_bgen;
        clearVarsSorted.(var).base.('ey') = ey_bgen;
        clearVarsSorted.(var).base.('ow') = ow_bgen;
        clearVarsSorted.(var).base.('ah') = ah_bgen;

        clearVarsSorted.(var).end.('iy') = i_train;
        clearVarsSorted.(var).end.('ae') = ae_train;
        clearVarsSorted.(var).end.('uw') = u_train;
        clearVarsSorted.(var).end.('aa') = a_train;
        clearVarsSorted.(var).end.('ih') = ih_gen;
        clearVarsSorted.(var).end.('eh') = eh_gen;
        clearVarsSorted.(var).end.('ey') = ey_gen;
        clearVarsSorted.(var).end.('ow') = ow_gen;
        clearVarsSorted.(var).end.('ah') = ah_gen;
    end
end


function [btPdist,bgPdist,tPdist,gPdist,baseTrainInds,posList] = measure_pdist(f1sIn,f2sIn,inds,vowelList_bt,vowelList_t,vowelList_bg,vowelList_g,vowels)
    for v = 1:4 %hardcoded
        vow = vowels{v};
        baseTrainInds.(vow) = inds.bTrainIndices(vowelList_bt == v);
        trainInds.(vow)     = inds.lastTrain(vowelList_t == v);

    end
    for v = 5:9
        vow = vowels{v};
        baseGenInds.(vow)   = inds.bGenIndices(vowelList_bg == v);
        genInds.(vow)       = inds.genIndices(vowelList_g == v);
    end
    
    posList.baselineTrain = [nanmean(f1sIn(baseTrainInds.iy)) nanmean(f2sIn(baseTrainInds.iy)); ...
        nanmean(f1sIn(baseTrainInds.ae)) nanmean(f2sIn(baseTrainInds.ae)); ...
        nanmean(f1sIn(baseTrainInds.uw)) nanmean(f2sIn(baseTrainInds.uw)); ...
        nanmean(f1sIn(baseTrainInds.aa)) nanmean(f2sIn(baseTrainInds.aa))];
    posList.baselineGen = [nanmean(f1sIn(baseGenInds.ih)) nanmean(f2sIn(baseGenInds.ih)); ...
        nanmean(f1sIn(baseGenInds.ey)) nanmean(f2sIn(baseGenInds.ey)); ...
        nanmean(f1sIn(baseGenInds.eh)) nanmean(f2sIn(baseGenInds.eh)); ...
        nanmean(f1sIn(baseGenInds.ow)) nanmean(f2sIn(baseGenInds.ow)); ...
        nanmean(f1sIn(baseGenInds.ah)) nanmean(f2sIn(baseGenInds.ah))];
    posList.train = [nanmean(f1sIn(trainInds.iy)) nanmean(f2sIn(trainInds.iy));  ...
        nanmean(f1sIn(trainInds.ae)) nanmean(f2sIn(trainInds.ae));...
        nanmean(f1sIn(trainInds.uw)) nanmean(f2sIn(trainInds.uw)); ...
        nanmean(f1sIn(trainInds.aa)) nanmean(f2sIn(trainInds.aa))];
    posList.gen = [nanmean(f1sIn(genInds.ih)) nanmean(f2sIn(genInds.ih)); ...
        nanmean(f1sIn(genInds.ey)) nanmean(f2sIn(genInds.ey)); ...
        nanmean(f1sIn(genInds.eh)) nanmean(f2sIn(genInds.eh)); ...
        nanmean(f1sIn(genInds.ow)) nanmean(f2sIn(genInds.ow)); ...
        nanmean(f1sIn(genInds.ah)) nanmean(f2sIn(genInds.ah));];
    
    btPdist = mean(pdist(posList.baselineTrain));
    bgPdist = mean(pdist(posList.baselineGen));
    tPdist  = mean(pdist(posList.train));
    gPdist  = mean(pdist(posList.gen));
    
    
end

%% helper functions
function vowstr = num2vowstr(num)
    if (num == 1)
        vowstr = 'iy';
    elseif (num == 2)
        vowstr = 'ae';
    elseif (num == 3)
        vowstr = 'uw';
    elseif (num == 4)
        vowstr = 'aa';
    elseif (num == 5)
        vowstr = 'ih';
    elseif (num == 6)
        vowstr = 'eh';
    elseif (num == 7)
        vowstr = 'ey';
    elseif (num == 8)
        vowstr = 'ow';
    elseif (num == 9)
        vowstr = 'ah';
    end
end

function d = euclidean(from, tos)
    for i = 1:size(tos, 1)
        d(i) = pdist([from; tos(i, :)]);
    end
end

function n = vector_norm(v)
    v = v(:, all(~isnan(v)));
    if (isempty(v))
        n = NaN;
        return
    end
    for i = 1:size(v, 2)
        n(i) = norm(v(:, i));
    end
end

function p = scalar_proj(a, theta) 
    a = a(:, all(~isnan(a)));
    theta = theta(all(~isnan(a)));
    p = NaN(1, size(a, 2));
    for i = 1:size(a, 2)
        p(i) = norm(a(:, i)).*cos(theta(i)); % modified
    end
end

function p = vector_proj(a, theta)
    a = a(:, all(~isnan(a)));
    theta = theta(all(~isnan(a)));
    p = a.*cos(theta);
end

function amp = getPertAmp(fmts, pertAmp)
    fieldDim = 257;
    F1Min = 200; F1Max = 1500;
    F2Min = 500; F2Max = 3500;
    pertf1 = floor(F1Min:(F1Max-F1Min)/(fieldDim-1):F1Max);
    pertf2 = floor(F2Min:(F2Max-F2Min)/(fieldDim-1):F2Max);
    amp = NaN(1, size(fmts, 2));
    for i = 1:size(fmts, 2)
        if ( isnan(fmts(2, i)) || isnan(fmts(1, i)) ); continue; end
        binF1 = 0; binF2 = 0;
        fmt1 = fmts(1, i); fmt2 = fmts(2, i);
        for j = 1:length(pertf1)-1
            if (pertf1(j) <= fmt1 && pertf1(j+1) >= fmt1)
                binF1 = j; break;
            end
        end
        for j = 1:length(pertf2)-1
            if (pertf2(j) <= fmt2 && pertf2(j+1) >= fmt2)
                binF2 = j; break;
            end
        end
        if (binF1 == 0 || binF2 == 0)
            warning('getPertAmp: Either f1 or f2 is out of bound.');
            continue;
        end
        amp(i) = pertAmp(binF2, binF1);
    end
end

function theta = getPertPhi(fmts, pertPhi,fLims)
    fieldDim = 257;
    pertf1 = floor(fLims.f1Min:(fLims.f1Max-fLims.f1Min)/(fieldDim-1):fLims.f1Max);
    pertf2 = floor(fLims.f2Min:(fLims.f2Max-fLims.f2Min)/(fieldDim-1):fLims.f2Max);
    theta = NaN(1, size(fmts, 2));
    for i = 1:size(fmts, 2)
        if ( isnan(fmts(2, i)) || isnan(fmts(1, i)) ); continue; end
        binF1 = 0; binF2 = 0;
        fmt1 = fmts(1, i); fmt2 = fmts(2, i);
        for j = 1:length(pertf1)-1
            if (pertf1(j) <= fmt1 && pertf1(j+1) >= fmt1)
                binF1 = j; break;
            end
        end
        for j = 1:length(pertf2)-1
            if (pertf2(j) <= fmt2 && pertf2(j+1) >= fmt2)
                binF2 = j; break;
            end
        end
        if (binF1 == 0 || binF2 == 0)
            warning('getPertPhi: Either f1 or f2 is out of bound.');
            continue;
        end
        theta(i) = pertPhi(binF2, binF1);
    end
end
