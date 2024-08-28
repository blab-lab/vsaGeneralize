function plot_vsaGeneralize_vowelSpaces(dataPaths)

adaptColor = [237 28 26]./255;
controlColor = [4 75 214]./255;
colors = [adaptColor;controlColor];
fullPageWidth = 17.4+3.0; % 174 mm + margins
colWidth = 8.5+3.0; % 85 mm + margins

plotParams.Marker = '.';
plotParams.MarkerSize = 12;
plotParams.MarkerAlpha = .25;
plotParams.LineWidth = .6;
plotParams.LineColor = [.7 .7 .7 .5];
plotParams.avgMarker = 'o';
plotParams.avgMarkerSize = 4;
plotParams.avgLineWidth = 1.25;
plotParams.jitterFrac = .25;
plotParams.FontSize = 13;

subj2plot = 1:20;
nSubs2plot = length(subj2plot);
for sidx = 1:length(subj2plot)
    [~,subjID] = fileparts(dataPaths{sidx});
    fdataByVowel.baseline = get_fdataByVowel(dataPaths{sidx},26:115);
    fdataByVowel.end = get_fdataByVowel(dataPaths{sidx},356:445);
    h2(sidx) = plot_VSA(fdataByVowel,[0 0 0; adaptColor],plotParams);
    axlim = [250 1000 750 2000];
    axis(axlim)
    pbaspect([1 1 1]); pbaspect manual;%axis square
    set(gca,'YTick',axlim(3):250:axlim(4))
    set(gca,'XTickLabel','')
    if sidx > 1
        set(gca,'YTickLabel','')
        ylabel('')
    else
        ylabel('F2 (mels)')
    end
    xlabel('')
    title(subjID)
    %axis normal;
    makeFig4Printing;
end

nCols = floor(sqrt(nSubs2Plot));
nRows = ceil(nSubs2plot/nCols);
figpos_cm = [15 15 fullPageWidth fullPageWidth*.565];
h = figure('Units','centimeters','Position',figpos_cm);
copy_fig2subplot(h2,h,nCols,nRows,[],1);