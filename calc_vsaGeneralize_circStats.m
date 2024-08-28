function [pvalGlobal,testStat, pvals, sigPairs, vowList, r, p] = calc_vsaGeneralize_circStats

dataPath = get_acoustLoadPath('vsaGeneralize');
load(fullfile(dataPath,'datatable_vowels.mat'))

[pvalGlobal,~,testStat] = circ_cmtest(datatable.ang, datatable.vowNUM);
vows = unique(datatable.vowel);
i = 0;
for v1 = 1:length(vows)-1
    vow1 = vows{v1};
    for v2 = v1+1:length(vows)
        i = i+1;
        vow2 = vows{v2};
        temp = datatable(ismember(datatable.vowel,vow1)| ...
                         ismember(datatable.vowel,vow2),:);
        pvals(i) = circ_cmtest(temp.ang, temp.vowNUM);
        vowList{i} = sprintf('%s_%s',vow1,vow2);
    end
    
end

[h,~,~,pvalsAdj] = fdr_bh(pvals);

sigPairs = vowList(h);

%create subset of only generalization data for correlation with distance to
%training vowels
gen = datatable(ismember(datatable.vowel,'ih')| ...
    ismember(datatable.vowel,'ey')| ...
    ismember(datatable.vowel,'eh')| ...
    ismember(datatable.vowel,'ah')| ...
    ismember(datatable.vowel,'ow'),:);

[r,p] = circ_corrcc(gen.ang, gen.trainDistAng);