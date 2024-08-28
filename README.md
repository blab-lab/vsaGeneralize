# vsaGenerelize Analysis Guideline

This is documentation for ```analyze_vsaGeneralize.m```

Author: Taijing Chen

- [Running analysis script](#running-analysis-script)

- [Naming/Indexing](#naming-indexing)

- [Data preprocessing](#data-preprocessing)
  - [preprocess()](#preprocess)
  -  [setup_formants()](#setup-formants)	
  - [setup_vowels()](#setup-vowels)	
  - [measure_changes()](#measure-changes)
  - [measure_pdist()](#measure-pdist)
- [Nearesting Training Vowel](#nearest-vowels)
  - [Self-defined nearest vowel](#self-defined-neareast-vowels)
  - [measure_nearest_vowels()](#measure-nearest-vowels)
  - [plot_dist_vs_sp_by_participant()](#plot-dist-vs-sp-by-participant)
  - [plot_pdist_train_vs_gen()](#plot-pdist-train-vs-gen)
- [Helper Function](#helper-function)





<h2 id=running-analysis-script> Running analysis script</h2>

The function `analyze_vsaGeneralize` takes in datapaths for participants as cells arrays of character vectors, returns figure handler(`h`), some preprocessed data (`data2plot`), and information about acoustic distances on F1-F2 space (`distances_euclidean`, `distances defined`)  For example:

```matlab
dPs = get_dataPaths_vsaGeneralize;
[h, data2plot,distances_euclidean,distances_defined] = analyze_vsaGeneralize(dPs);
```





<h2 id=naming-indexing> Naming/Indexing </h2>

Here're how the variables are generally named:

- b: baseline
- t: training
- g/gen: generalization (testing)
- pdist: pairwise distances
- pos:  vowel positions on F1-F2 space, i.e., [f1 f2]
- vec_proj/vp: vector project, (with direction)
- scalar_proj/sp: scalar projection (magnitude of vector projection, w/t direction)

Here're how the vowels are indexed (check fuction `num2vowstr()` in the helper function section):

- 1 - iy/i
- 2 - ae
- 3 - uw/u
- 4 - aa/a
- 5 - ih
- 6 - eh
- 7 - ey
- 8 - ow
- 9 - ah

If the vowels are grouped by conditions, then for the generalization vowels (ih, eh, ey, ow, ah), subtract their indices by 4 (number of training vowels). (e.g., 1 - ih, 2 - eh, ...)



<h2 id=data-preprocessing> Data preprocessing </h2>

**<u>NOTE</u>: Many preprocessed fields are redundant. I kept them there in case there're changes to make in data analysis, but all of them can be deleted. **

Data preprocessing are done in function `preprocess`, which includes 4 helper functions `setup_formants`, `setup_vowels` `measure_change`, `measure_pdist`.

<h3 id=preprocess> preprocess() </h3>

```matlab
[data2plot] = preprocess(snums, exptName)
```

This function preprocess all data for plotting.

- `snums`: cells arrays of character vectors of participant numbers. For example:

  ```matlab
  snums = {'sp189', 'sp235', 'sp252', 'sp253'};
  ```

- `exptName`: by default, "vsaGeneralize"

- `data2plot`: stores all data for plotting

  - `btPdists, bgPdists, tPdists, gPdists`: pairwise distances grouped by conditions; not normalized (check [Naming/Indexing](#naming-indexing) for variable naming) (**<u>NOTE</u>: not used; can be deleted**)
  - `tPdistAbsNormalized, gPdistAbsNormalized`: pairwise distainces normalized by **subtracting** the baseline values in the training or generalization condition, respectively, i.e., `tPdists - btPdists` and `gPdists -bgPdists`; 1 x #participants
  - `tPdistPercentNormalized, gPdistPercentNormalized`: pairwise distainces normalized by **dividing** the baseline values in the training or generalization condition, respectively, i.e., `tPdists ./ btPdists` and `gPdists ./ bgPdists`; 1 x #participants

- `compensations, scalar_projs, vector_projs`: cell array of compensations/scalar projections/vector projections of the final behaviral changes for all participants; 1 x #participants. You can index into values for each vowel by `iy,ae,uw,aa,ih,eh,ey,ow,ah`. You'll get a 2x50 or 2x40 matrix that stores the corresponding [F1 F2] values, depended on conditions. For example:

  ```matlab
  data2plot.compensations{1} % all compensations of participant 1
  data2plot.scalar_projs{1}  % all scalar projections of participant 1
  
  data2plot.compensations{1}.iy  % the change vectors of participant 1 in the last 10 training trials for iy
  data2plot.compensations{1}.ih  % the change vectors of participant 1 in the testing trials for ih
  ```

  <b id=badTrails><u>NOTE</u>:  These 3 fields contain badTrials. Use <code>nanmean()</code> to calculate average. </b>

- `posLists`: cell arrays of F1-F2 position of vowels for all participants. You can access a 5x2 or 4x2 matrix for in each condition by `baselineTrain, baselineGen, train, gen`. In this matrix, each row stores [F1 F2] values for some vowel, whose index follows the [Naming/Indexing](#naming-indexing) section. For example:

  ```matlab
  data2plot.posLists{1} % vowel positions of participant 1
  data2plot.posLists{1}.baselineTrain % vowel positions of participant 1 in the baseline training phase
  data2plot.posLists{1}.baselineTrain(1,:) % iy's position of participant 1 in the baseline training phase
  data2plot.posLists{1}.baselineTrain(1,:) % ih's position of participant 1 in the baseline training phase
  ```

<h3 id=setup-formants> setup_formants </h3>

```matlab
function [f1sIn,f2sIn,f1sOut,f2sOut] = setup_formants(dataVals, ntrials)
```

This function extracts formant ins and outs from the data of each participant.

- `dataVals`: data from the  `dataVals.mat` file of some participant p
- `ntrails`: number of trials of some participant p
- `f1sIn, f2sIn, f1sOut, f2sOut`: as name suggested, ins and outs of first two formants. They are the mean formant values between 0.25 timepoint to 0.5 timepoint.  (**<u>NOTE</u>: `f1sOut, f2sOut` are never used; can be deleted**)

<h3 id=setup-vowels> setup_vowels() </h3>

```matlab
function [vows,vsa] = setup_vowels(f1sIn,f2sIn,binSizeTrain,binSizeGen,inds,allVowels)
```

From formants, this function groups vsa and formants by vowel for some participant p.

- `f1sIn, f2sIn`: as name suggested, ins and outs of first two formants. They are the mean formant values between 0.25 timepoint to 0.5 timepoint.
- `binSizeTrain, binSizeGen`: # of vowels for training and testing (binSizeTrain = 4, binSizeGen = 5). (**<u>NOTE</u>: There're places where I didn't pass this variables to function and hardcoded them as 4 and 5 inside functions**)
- `inds`:
  - `inds.bGenIndices`: vowel indices for the baseline generalization group. **There're 75 (15 x binSinzeGen) of them, but we only use last 50 (10 x binSinzeGen). This is because in this condition, participants' auditory feedbacks were masked. When piloting, participants report they felt odd about the masked noise. So we decided to extend the baseline phase for generalization vowel to 75 trials and only use the last 50 of them.**
  - `inds.bTrainIndices`: vowel indices for the baseline training group
  - `inds.trainIndices`: vowel indices for the training group (**<u>NOTE</u>: not used; can be deleted**)
  - `inds.lastTrain`:  indices for the last 40 (10 x binSizeTrain) vowels in training group
  - `inds.genIndices`: vowel indices for the generalization group
- `allVowels`: allVowels field in `expt.mat`
- `vows`: formants grouped by vowels & conditions
  - 2x50 (baseline testing) `vows.ih_bgen, vows.eh_gen, vows.ey_gen, vows.ow_gen, vows.ah_gen`
  - 2x40 (baseline training):`vows.u_btrain, vows.i_btrain, vows.a_btrain,vows.ae_btrain`
  - 2x40 (last 40 training trials):`vows.u_train, vows.i_train, vows.a_train, vows.ae_train`
  - 2x50 (testing):`vows.ih_gen, vows.eh_gen, vows.ey_gen, vows.ow_gen, vows.ah_gen`
- `vsa`:  areas returned by matlab function `polyarea` (**<u>NOTE</u>: not used and can be deleted. If deleted, remember to delete `vsaTraining` and `vsaGeneralize` in data2plot as well. **)
  - `vsa.VSA_train, vsa.VSA_btrain`: polyarea(iy, aa, uw), 2x10
  - `vsa.VSA_gen, vsa.VSA_bgen`: polyarea(ih, ey, eh, ow, ah), 2x10 (**very likely to be NaN**)

<h3 id=measure-changes> measure_changes() </h3>

```matlab
function [compensations,scalar_projs,vector_projs] = measure_changes(vows,pertPhi)
```

This function measures and calculates behavioral changes for some participant p (in mels).

- `vows`: formants grouped by vowels & conditions (check section [setup_vowels](#setup-vowels) for details)
- `pertPhi`: perturbation angles for each participant; 257x257. (check helper function `getPertPhi` for indexing details)

The following returns can be indexed by `iy,ae,uw,aa,ih,eh,ey,ow,ah`; all in *mels*

- `compensations`: formants - baseline_formants; 2x50 (testing vowels) or 2x40 (training vowels)
- `scalar_projs`: scalar projections of `compensations`, 1x50 (testing vowels) or 1x40 (training vowels)
- `vector_projs`: vector projections of `compensations`, 2x50 (testing vowels) or 2x40 (training vowels)

<h3 id=measure-pdist> measure_pdist() </h3>

```matlab
function [btPdist,bgPdist,tPdist,gPdist,baseTrainInds,posList] = measure_pdist(f1sIn,f2sIn,inds,vowelList_bt,vowelList_t,vowelList_bg,vowelList_g,vowels)
```

This function measures pairwise distances (mean of returns from Matlab function `pdist`) for some participant p (in mels). 

- `f1sIn,f2sIn`: formant Ins, check [setup_formants](#setup-formants) for details

- `inds`: check [setup_vowels](#setup-vowels) for details

- `vowelList_bt = expt.allVowels(inds.bTrainIndices);`

- `vowelList_t = expt.allVowels(inds.lastTrain);`

- `vowelList_bg = expt.allVowels(inds.bGenIndices);`

- `vowelList_g = expt.allVowels(inds.genIndices);`

- `vowels = expt.vowels`

- `btPdist,bgPdist,tPdist,gPdist`: pairwise distances grouped by conditions (check [Naming/Indexing](#naming-indexing) for details )

- `baseTrainInds`: stores vowel indices for the baseline training group (**<u>NOTE</u>: replaced by posList so not be used; can be deleted**)

  - `baseTrainInds.iy, baseTrainInds.ae, baseTrainInds.uw, baseTrainInds.aa`

- `posList`: F1-F2 position lists of vowels; grouped by conditions; used by `meansure_nearest_vowel`

  - `posList.baselineTrain, posList.baselineGen, posList.train, posList.gen`: mean F1-F2 values of some vowel formants

  - You can index into each vowel following indexing convention in [Naming/Indexing](#naming-indexing). For example:

    ```matlab
    posList.baselineTrain(1,:) % position of iy in baseline training phase
    posList.baselineGen(1,:)   % position of ih in baseline generalization phase 
    ```



<h2 id=nearest-vowels> Nearest Training Vowels </h2>

<h3 id=self-defined-neareast-vowels> Self-defined nearest vowels </h3>

- ih - iy
- ey - mean(iy, ae)
- eh - ae
- ow - mean(uw, aa)
- ah - aa

<h3 id=measure-nearest-vowels> measure_nearest_vowels() </h3>

```matlab
function [distances_euclidean,distances_defined,nearestTrainVowel] = measure_nearest_vowels(posLists)
```

This function calculates the distances to the nearest training vowels of each testing vowels for all particpants

- `posLists`: vowel positions for all participants, i.e., `data2plot.posLists`. Check the [preprocess](#preprocess) section for details. For example:

  ```matlab
  posLists{1} % vowel positions of participant 1
  posLists{1}.baselineTrain % vowel positions of participant 1 in the baseline training phase
  posLists{1}.baselineTrain(1,:) % iy's position of participant 1 in the baseline training phase
  posLists{1}.baselineGen(1,:)   % ih's position of participant 1 in the baseline training phase
  ```

- `distances_euclidean`: euclidean distances from each testing vowels to the **calculated** nearest trainig vowels for all participants; #participants x #testing-vowels. For example:

  ```matlab
  distances_euclidean(1,:) % distances to the calculated nearest training vowels of participant 1
  distances_euclidean(1,1) % distance from ih to its nearest training vowels of participant 1
  ```

- `distances_defined`: euclidean distances from each testing vowels to the **self-defined** nearest trainig vowels foor all participants. You can access values for each testing vowels by `ih2iy, ey2front, eh2ae, ou2back, ah2aa`, each containing a vector of length == #participants. 

  ```matlab
  distances_defined.ih2iy    % distances from ih to iy for all participants
  distances_defined.ih2iy(1) % distances from ih to iy for participant 1
  ```

- `nearestTrainVowel`: nearest training vowels of each testing vowels on F1-F2 space for all participants; #participants x #testing-vowels. It is an int matrix whose values correspond to [Naming/Indexing](#naming-indexing). For example:

  ```matlab
  nearestTrainVowel(1, 1) == 1 % For participant 1, the nearest training vowel for testing vowel ih is iy
  nearestTrainVowel(2, 1) == 2 % For participant 2, the nearest training vowel for testing vowel ih is ae
  ```



<h2 id=plotting> Plotting </h2>

```matlab
function h = plotting(h, data2plot,num,distances_euclidean,distances_defined,nearestTrainVowel)
```

All formant values are in Mels.

- `h`: figure handlers
- `data2plot`: check [preprocess](#preprocess)
- `num`: number of participants, 
- `distances_euclidean,distances_defined,nearestTrainVowel`: check [measure_nearest_vowels](#measure-nearest-vowels)

<h3 id=plot-pdist> plot_pdist() </h3>

```matlab
function h = plot_pdist(h,data2plot)
```

This function plots out two pairwise distances grouped by training vowels and testing conditions. The first one is normalized by **subtracting** the corresponding baseline values and the second one is normalized by **dividing** the corresponding baseline values, i.e., `data2plot.tPdistAbsNormalized` and `data2plot.gPdistAbsNormalized` & `data2plot.tPdistPercentNormalized` and `data2plt.gPdistPercentNormalized`. Check [preprocess](#preprocess) for detailed explanations of these fileds.

This function also prints out t-test results of the second figure (normalized by division) on the console:

```matlab
[h_t,p_t,ci_t,stats_t] = ttest(data2plot.tPdistPercentNormalized, 1)
[h_g,p_g,ci_g,stats_g] = ttest(data2plot.gPdistPercentNormalized, 1)
```

<h3 id=plot-sp-normalized> plot_sp_normalized() </h3>

```matlab
function h = plot_sp_normalized(h,data2plot,nearestTrainVowel,num)
```

**TODO: `data2plot` can be changed into `data2plot.scalar_projs` to reduce copying.**

This function plots out two scalar projections of the compensations of each testing vowels normalized by the ones of its nearest training vowels. 

It processes `data2plot.scalar_projs` field and store corresponding data into variable `defs` and `cals` , both of which store scalar projections of the generalization vowels normalized by the ones of their nearest traning vowels (grouped by participant). Thus, each has dimensions of #participants x 5 (#testing-vowels). `defs` stores those normalized by self-defined nearest training vowels, and `cals` stores those normalized by calculated nearest training vowels. `nanmean()` is used here to avoid [badTrials](#badTrails) in `data2plot.scalar_projs`.

This funcion also prints out t-test results of both figures on the console:

```matlab
% all scalar projection normalized by their calculated nearest training vowels vs 1
[h_cal,p_cal,ci_cal,stats_cal] = ttest(cals(:), 1) 
% all scalar projection normalized by their self-defined nearest training vowels vs 1
[h_def,p_def,ci_def,stats_def] = ttest(defs(:), 1)
% means of normalized scalar projections grouped by testing vowels (calculated) vs 1
[h_cal_mean,p_cal_mean,ci_cal_mean,stats_cal_mean] = ttest(nanmean(cals, 1), 1)
% means of normalized scalar projections grouped by testing vowels (self-defined) vs 1
[h_def_mean,p_def_mean,ci_def_mean,stats_def_mean] = ttest(nanmean(defs, 1), 1)
```

**TODO: `nanmean()` can be changed into `mean()`.**

<h3 id=plot-dist-vs-sp-by-participant> plot_dist_vs_sp_by_participant() </h3>

```matlab
function h = plot_dist_vs_sp_by_participant(h,sps,posLists,nearestTrainVowel,num)
```

This function plots out euclidean distances between the testing and training vowels vs the correspoinding scalar projections of the testing vowels, grouped by participant. 

- `sps`: `data2plot.scalar_projs ` (check [preprocess](#preprocess) for details)
- `posLists`: `data2plot.posLists` (check [preprocess](#preprocess) for details)

It processes `sps` field into variable `y` and `posLists` field into variable `x_cals, x_defs`. `y` stores scalar projections of testing vowels, grouped by participants. `x_cals` stores euclidean distances between the testing vowels to its self-defined nearest traing vowels, and `x_defs` stores euclidean distances between the testing vowels to the calculated nearest training vowels. All 3 variables have dimensions of #participants x #testing-vowels.

**TODO: may need some statistic test here.**

<h3 id=plot-pdist-train-vs-gen> plot_pdist_train_vs_gen() </h3>

``` matlab
function h = plot_pdist_train_vs_gen(h,tPdistAbsNormalized,gPdistAbsNormalized,tPdistPercentNormalized,gPdistPercentNormalized)
```

This function plots out the pairwise distances of the training vowels vs the pairwise distances of the testing vowels.

- `tPdistAbsNormalized,gPdistAbsNormalized,tPdistPercentNormalized,gPdistPercentNormalize`: corresponding fields in `data2plot`; (check [preprocess](#preprocess) for details)

Two spearman test results are printed to the consoles:

```matlab
% correlation between testing and training (pdist normalized by subtraction from the baseline)
[rho_abs, p_abs] = corr(tPdistAbsNormalized', gPdistAbsNormalized', 'type', 'spearman')
% correlation between testing and training (pdist normalized by division by the baseline)
[rho_percent, p_percent] = corr(tPdistPercentNormalized', gPdistPercentNormalized', 'type', 'spearman')
```



<h2 id=helper-functions> Helper Functions </h2>

<h3 id=num2vowstr> num2vowstr() </h3>

```matlab
function vowstr = num2vowstr(num)
```

This function convert integer idex into corresponding ipa strings used in this script. (Check [Naming/Indexing](#naming-indexing)).

- `num`: integer 1-9. If the structure/variable only contains testing vowels (instead of training + testing vowels), you may want to add 4 (i.e., #testing-vowels) before passing to this function. 
- `vowstr`: (Check [Naming/Indexing](#naming-indexing)).

<h3 id=euclidean> euclidean() </h3>

```matlab
function d = euclidean(from, tos)
```

This function calculates euclidean distances from `from` (as name suggested, a vector) to `tos` (a list of vectors). 

- `from`: a 1 x 2 row vector 
- `tos`:   a n x 2 matrix
- `d`:       a 1 x n vector

For example:

```matlab
% returns a scalar; the distance between the position of iy and the position of ih in baseline
euclidean(posLists{p}.baselineGen(1, :), posLists{p}.baselineTrain(1, :));
% returns a 1x4 vetor; distancs between the position of iy and the positions of all training vowels in baseline
euclidean(posLists{p}.baselineGen(1, :), posLists{p}.baselineTrain(:, :));
```

<h3 id=scalar-proj> scalar_proj() </h3>

```matlab
function p = scalar_proj(a, theta) 
```

This function calculates scalar projections of a 2xn matrix `a`, with a 1xn vector of angles `theta`. For example:

```matlab
>> scalar_proj([1 2 3; 0 0 0], [pi/6 pi/6 pi/6])
>> ans =  0.8660    1.7321    2.5981
```

<h3 id=vector-proj> vector_proj() </h3>

```matlab
function p = vector_proj(a, theta)
```

This function calculates vector projections of a 2xn matrix `a`, with a 1xn vector of angles `theta`. For example:

```matlab
>> vector_proj([1 2 3; 0 0 0], [pi/6 pi/6 pi/6])
>> ans = 0.8660    1.7321    2.5981
         0         0         0
```

**NOTE: scalar_projection == column_major_norm(vector_projection)** For example:

```
scalar_proj([1 ; 2], pi/6) == norm(vector_proj([1 ; 2], pi/6))
```







