function [h] = plot_vsaGeneralize_PaperFigs(figs2plot,data2plot)
%PLOT_VSAPAPERFIGS  Plot figures for the vsaAdapt2 paper.

if nargin < 1 || isempty(figs2plot), figs2plot = 1; end
if nargin < 2 || isempty(data2plot)
    loadPath = get_acoustLoadPath('vsaGeneralize');
    load(fullfile(loadPath,'processedData_25-75.mat'))
end
dataPaths = get_dataPaths_vsaGeneralize;

set(0, 'DefaultFigureRenderer', 'painters');

trainColor = [237 28 26]./255;
genColor = [4 75 214]./255;
colors = [trainColor;genColor];
fullPageWidth = 17.4+3.0; % 174 mm + margins
%colWidth = 8.5+3.0; % 85 mm + margins

plotParams.Marker = '.';
plotParams.MarkerSize = 25;
plotParams.MarkerAlpha = .4;
plotParams.LineWidth = .8;
plotParams.LineColor = [.7 .7 .7];
plotParams.avgMarker = 'o';
plotParams.avgMarkerSize = 8;
plotParams.avgLineWidth = 2;
plotParams.jitterFrac = .25;
plotParams.FontSize = 13;

vowColors.iy = [.4 .7 .06]; %[78 155 11]/255; %[.55 .85 .15];
vowColors.ae = [.8 0 .4];
vowColors.aa = [.1 0 .9];
vowColors.uw = [.1 .6 .9];

%% Fig 1: Experiment design

[bPlot,ifig] = ismember(1,figs2plot);
if bPlot
    
    %% A: pert field TODO
    sid = 'sp189';
    dataPath = get_acoustLoadPath('vsaGeneralize',sid);
    fmtMeans = calc_vowelMeans(fullfile(dataPath,'pre'));
    [~,h1(1)] = calc_pertField('in',fmtMeans,1);
    axlim = [350 1050 950 1850];
    axis(axlim)
    fmtMeans = calc_vowelMeans(dataPath,{'baselineGeneralize'});
    vows2plot = {'ih','ey','eh','ow','ah'};
    for v = 1:length(vows2plot)
        vow = vows2plot{v};
        text(hz2mel(fmtMeans.(vow)(1)),hz2mel(fmtMeans.(vow)(2)),arpabet2ipa(vow),...
            'HorizontalAlignment','center')
    end
    
    set(gca,'XTick',axlim(1)+50:200:axlim(2))
    set(gca,'YTick',axlim(3)+50:200:axlim(4))
    axis normal;
    makeFig4Printing;
    
    %% B: spectrogram TODO
    trials2plot = [118 117 116 119];
    params.fmtsColor = genColor;
    params.fmtsLineWidth = 3;
    params.sfmtsColor = trainColor;
    params.sfmtsLineWidth = 1.5;
    params.ylim = 4000;
    %params.figpos = [35 700 2510 150];
    
    fCen = calc_vowelCentroid(fmtMeans);
    params.fmtCen = fCen;
    params.thresh_gray = .65;
    params.max_gray = .75;
    
    load(fullfile(dataPath,'data.mat'),'data');
    plot_audapterFormants(data(trials2plot),params);
    h1(2) = plot_audapterFormants(data(trials2plot(1)),params);
    cla;
    set(gca,'YTick',[1000 1949 3429]); % 1000, 1500, 2000 in mels
    set(gca,'YTickLabel',[500 1000 2000])
    set(gca,'XTick',[])
    makeFig4Printing;
    
    %% C: trial timeline
    h1(3) = figure;
    hold on;
    plot([1 75],[0 0],'Color',genColor,'LineWidth',2);
    plot([76 115 116 395],[0 0 50 50],'Color',trainColor,'LineWidth',2);
    plot([396 445],[0 0],'Color',genColor,'LineWidth',2);
    
    
    axlim = [0 446 -15 75];
    axis(axlim);
    set(gca,'YTick',[0 50]) %set(gca,'YTick',axlim(3):25:axlim(4))
    set(gca,'YTickLabel',{'0' '50%'}) %set(gca,'YTickLabel',{'' '0' '' '50%' ''})
    hline(0,'k','--');
    vline(75.5,'k',':');
    vline(115.5,'k',':');
    vline(395.5,'k',':');
    set(gca,'XTick',[1 75 115 355 395])
    xlabel('trial number')
    ylabel('perturbation')
    
    ypos = 60; fontsize = 7; rot = 55;
    text(1+75/2,ypos,'baseline','FontSize',fontsize,'HorizontalAlignment','center','Rotation',rot)
    text(76+40/2,ypos,'baseline','FontSize',fontsize,'HorizontalAlignment','center','Rotation',rot)
    text(116+280/2,ypos,'exposure','FontSize',fontsize,'HorizontalAlignment','center','Rotation',rot)
    text(356+40/2,ypos,'adaptation','FontSize',fontsize,'HorizontalAlignment','center','Rotation',rot)
    text(396+50/2,ypos,'generalization','FontSize',fontsize,'HorizontalAlignment','center','Rotation',rot)
    
    ax = axis; ylims = ax([3 4]);
    h_fill = fill([356 395 395 356],[ylims(1) ylims(1) ylims(2) ylims(2)],[.9 .9 .9],'EdgeColor','none');
    uistack(h_fill,'bottom');

    makeFig4Printing;
    %set(gca,'XColor','none')

    %% ALL
    figpos_cm = [1 29 fullPageWidth fullPageWidth*(325/1244)]; %figpos = [58 970 1244 325];
    h(ifig) = figure('Units','centimeters','Position',figpos_cm);
    copy_fig2subplot(h1,h(ifig),2,3,{[1 4] [2 3] [5 6]},1);%h1,
    
end

%% Fig 2: Increases in AVS

[bPlot,ifig] = ismember(2,figs2plot);
if bPlot
    
    nSubs = length(dataPaths);
    clear ylims;
    ylims = [0.98 1.1];
    ylimsPaired = [0.9 1.25];
    ylimSpacingPaired = 0.1;
    ylimSpacing = 0.05;
    h2 = [];
    
    % normalize data
    for s = 1:nSubs
        vsTrack(s,2:9) = data2plot.AVS(s,2:9) ./ data2plot.AVS(s,2);
        vsTrack(s,[1 10]) = data2plot.AVS(s,[1 10]) ./ data2plot.AVS(s,1);
    end
    
    %% A,B: vs tracks, all subjects
    h2(end+1) = figure;    
    lineWidth = 2;
    
    errorbar(2:9,mean(vsTrack(:,2:9)),ste(vsTrack(:,2:9)),'o-','Color',trainColor,'MarkerFaceColor',trainColor,'MarkerSize',plotParams.avgMarkerSize,'LineWidth',lineWidth);
    hold on
    errorbar(1,mean(vsTrack(:,1)),ste(vsTrack(:,1)),'o-','Color',genColor,'MarkerFaceColor',genColor,'MarkerSize',plotParams.avgMarkerSize,'LineWidth',lineWidth);
    errorbar(10,mean(vsTrack(:,10)),ste(vsTrack(:,10)),'o-','Color',genColor,'MarkerFaceColor',genColor,'MarkerSize',plotParams.avgMarkerSize,'LineWidth',lineWidth);
    
    axlim = [0.5 10.5 ylims(1) ylims(2)];
    axis(axlim);
    set(gca,'YTick',axlim(3):ylimSpacing:axlim(4),...
        'XTick',[1 2 6 9 10],'YTick',[1 1.05 1.1])
    set(gca,'XTickLabel',{'baseline', 'baseline', 'exposure', 'adaptation', 'generalization'})
    xtickangle(30)
    h_lines = hline(1,'k','--');
    h_lines(end+1) = vline(1.5,'k',':');
    h_lines(end+1) = vline(2.5,'k',':');
    h_lines(end+1) = vline(9.5,'k',':');
    for hl = 1:length(h_lines)
        h_lines(hl).HandleVisibility = 'off';
    end
    h_fill = fill([8.5 9.5 9.5 8.5],[ylims(1) ylims(1) ylims(2) ylims(2)],[.9 .9 .9],'EdgeColor','none');
    uistack(h_fill,'bottom');
    ylabel('Normalized AVS')
    makeFig4Printing;
    
    %% B: paired data, adaptation/washout/retention phases

    pairedData.train = vsTrack(:,9)';
    pairedData.gen = vsTrack(:,10)';
    h2(end+1) = plot_pairedData(pairedData,colors,plotParams);
   
    set(gca,'YLim',ylimsPaired)
    set(gca,'YTick',ylimsPaired(1):ylimSpacingPaired:ylimsPaired(2))
    set(gca,'XTickLabel',{'adaptation', 'generalization'})
    xtickangle(30)
    ylabel('Normalized AVS')

    h_line = hline(1,'k','--');
    h_line.HandleVisibility = 'off';
    makeFig4Printing;
    
    %% C: correlation
    h2(end+1) = figure;
    [r,p] = corr(pairedData.train',pairedData.gen','type','Spearman');
    plot(pairedData.train,pairedData.gen,'k.','MarkerSize',plotParams.MarkerSize)
    pfit = polyfit(pairedData.train,pairedData.gen,1);
    ax = axis; xrange = ax(1)-5:ax(2)+5;
    hold on
    plot(xrange,xrange*pfit(1)+pfit(2),'k'); axis(ax);
    text(ax(2)-.3*diff(ax(1:2)),.96,{sprintf('r = %.3f',r),sprintf('p = %.3f',p)})
    xlabel('Adaptation (norm. AVS)')
    ylabel({'Generalization';'(norm. AVS)'})
    
    %% ALL
    figpos_cm = [1 15 fullPageWidth fullPageWidth*.565];
    tiledParams.TileSpacing = 'normal';
    h(ifig) = figure('Units','centimeters','Position',figpos_cm);
    copy_fig2tiledlayout(h2,h(ifig),2,4,[1 5 7],{[1 4] [1 2] [1 2]},1,tiledParams);

    

end

%% Fig 3: Vowel-specific distance from center

[bPlot,ifig] = ismember(3,figs2plot);
if bPlot
    h3 = [];
    nSubs = length(dataPaths);
    
    
    measures = {'dists'};
    ylabs = {'\Delta dist. to center','Compensation projection', };
    nMeas = length(measures);
    for m = 1:nMeas
        meas = measures{m};
        vDat.iy = NaN(nSubs,1);
        vDat.ih = NaN(nSubs,1);
        vDat.ey = NaN(nSubs,1);
        vDat.eh = NaN(nSubs,1);
        vDat.ae = NaN(nSubs,1);
        vDat.aa = NaN(nSubs,1);
        vDat.ah = NaN(nSubs,1);
        vDat.ow = NaN(nSubs,1);
        vDat.uw = NaN(nSubs,1);
        vowels = fieldnames(vDat);
        nVowels = length(vowels);
        for s = 1:nSubs
            for v = 1:nVowels
                vow = vowels{v};
                vDat.(vow)(s) = data2plot.(meas){s}.(vow);
            end
        end
        plotParams.bPaired = 0;
        vColors = [trainColor; repmat(genColor,3,1); repmat(trainColor,2,1); repmat(genColor,2,1); trainColor];
        h3(end+1) = plot_pairedData(vDat,vColors,plotParams);
        set(gca,'XTickLabel',arpabet2ipa(vowels))
        ylabel(ylabs{m})
        hline(0,'k','--');
        makeFig4Printing
    end
    
    %plot response angles
    if ~exist('circ_var','file')
        warning('circular statistics package not detected. no mean angle will be shown.')
        warning('https://github.com/circstat/circstat-matlab/blob/master/circ_var.m')
    end
    
    h5 = [];
    h52 = [];
    nSubs = length(dataPaths);
    
    vowels = {'iy','ae','aa','uw','ih','ey','eh','ah','ow'};
    altVowelNames = {'i','ae','a','u','ih','ey','eh','ah','ow'};
    
    nVowels = length(vowels);
    
    bTrain = logical([ones(1,4) zeros(1,5)]);
    
    for v = 1:nVowels
        vow = vowels{v};
        h3(end+1) = figure();
        %         h52(end+1) = figure();
        hold on
        x = NaN(1,nSubs);
        y = x;
        for s = 1:nSubs
            %compensation angle
            angle(s) = get_angle([0 0],data2plot.compensations{s}.(vow));
            
            %compensation magnitude
            %                 mag = pdist([0 0;data2plot.compensations{s}.(vow)]);
            mag = 1;
            
            %baseline vowel angle
            if bTrain(v)
                field2load = strcat(altVowelNames{v},'_btrain');
                plotColor = trainColor;
                meanLineStyle = '-';
            else
                field2load = strcat(altVowelNames{v},'_bgen');
                plotColor = genColor;
                meanLineStyle = ':';
            end
            
            bAngle(s) = get_angle(data2plot.vows{s}.(field2load),data2plot.fCenField{s});
            
            hold on
            %                 dAngle(s) = t/2*pi*(rand(1,1)-0.5); %for testing variance range
            x(s) = mag.*cos(angle(s));
            y(s) = mag.*sin(angle(s));
            h_q = quiver(0,0,x(s),y(s),0);
            set(h_q,'Color',get_lightcolor(plotColor,0.5),'LineWidth',1,'MaxHeadSize',0.5)
            
            hold on
            x(s) = mag.*cos(bAngle(s));
            y(s) = mag.*sin(bAngle(s));
            h_q = quiver(0,0,x(s),y(s),0);
            set(h_q,'Color',get_lightcolor([0 0 0],0.75),'LineWidth',1,'MaxHeadSize',0.5)
        end
        for r = .25:.25:1
            rectangle('Position',[-r -r 2*r 2*r],'Curvature',[1,1],...
                'EdgeColor',plotParams.LineColor,'FaceColor','none');
        end
        try
            mag = (1-circ_var(angle));
            x = mag.*cos(circ_mean(angle));
            y = mag.*sin(circ_mean(angle));
            h_q = quiver(0,0,x,y,0);
            set(h_q,'Color',plotColor,'LineWidth',3,'MaxHeadSize',0.6./mag)
            
            mag = (1-circ_var(bAngle));
            x = mag.*cos(circ_mean(bAngle));
            y = mag.*sin(circ_mean(bAngle));
            h_q = quiver(0,0,x,y,0);
            set(h_q,'Color','k','LineWidth',3,'MaxHeadSize',0.6./mag,'LineStyle',meanLineStyle)
                      
            
            
        catch
        end
        
        set(gca,'YTickLabel',[],'XTickLabel',[])
        title(arpabet2ipa(vow))
        axis equal
        makeFig4Printing
    end
    

    
    figpos_cm = [1 15 fullPageWidth fullPageWidth*4/5];
    tiledParams.TileSpacing = 'compact';
    h(ifig) = figure('Units','centimeters','Position',figpos_cm);
    copy_fig2tiledlayout(h3,h(ifig),4,5,[1 11 15 20 16 12 13 14 19 18],[{[2 5]} repmat({[1 1]},1,9)],1,tiledParams);
end

%% Fig 4: adaptation angles relative to training vowels, bar plot version

[bPlot,ifig] = ismember(41,figs2plot);
if bPlot
    h4 = [];
    nSubs = length(dataPaths);
    
    tVowels = {'iy','ae','aa','uw'};
    nTVowels = length(tVowels);
    
    vowels = {'ih','ey','eh','ah','ow'};
    nVowels = length(vowels);
    
    for v = 1:nVowels
        vow = vowels{v};
        for s = 1:nSubs  
            %compensation angle
            angle = get_angle([0 0],data2plot.compensations{s}.(vow));
            for t = 1:nTVowels
                tVow = tVowels{t};
                %training vowel angle
                tAngle = get_angle([0 0],data2plot.compensations{s}.(tVow));
                
                dAngle = wrapTo2Pi(angle-tAngle);

                if dAngle > pi
                    dAngle = dAngle - (2*pi);
                elseif tAngle < (-1 * pi)
                    dAngle = 2*pi + dAngle;
                end
                
                if dAngle > pi || dAngle < -1*pi
                    error('angle is out of range')
                end
                
                vDat.(vow).(tVow)(s) = abs(rad2deg(dAngle));
            end
        end
        plotParams.bPaired = 0;
        h4(end+1) = plot_pairedData(vDat.(vow),[],plotParams);
%         if v ~= nVowels
%             set(gca,'XTickLabel',[])
%         end
        if v ==1
            ylabel('angular difference')
        else
            set(gca,'YTickLabel',[])
        end
        title(vow)
%         hline(0,'k','--');
        makeFig4Printing
    end
    
    figpos_cm = [1 15 fullPageWidth .5*fullPageWidth];
    tiledParams.TileSpacing = 'compact';
    h(ifig) = figure('Units','centimeters','Position',figpos_cm);
    copy_fig2tiledlayout(h4,h(ifig),1,nVowels,[],[],1,tiledParams);
end

%% Fig 4: adaptation angles relative to training vowels, radial plot version
[bPlot,ifig] = ismember(42,figs2plot);
if ~exist('circ_var','file')
    warning('circular statistics package not detected. no mean angle will be shown.')
    warning('https://github.com/circstat/circstat-matlab/blob/master/circ_var.m')
end
if bPlot
    h4 = [];
    nSubs = length(dataPaths);
    
    tVowels = {'iy','ae','aa','uw'};
    nTVowels = length(tVowels);
    
    vowels = {'ih','ey','eh','ah','ow'};
    nVowels = length(vowels);
    
    for v = 1:nVowels
        vow = vowels{v};
        for t = 1:nTVowels
            tVow = tVowels{t};
            h4(end+1) = figure();
            hold on
            x = NaN(1,nSubs);
            y = x;
            for s = 1:nSubs
                %compensation angle
                angle = get_angle([0 0],data2plot.compensations{s}.(vow));
                
                %compensation magnitude
%                 mag = pdist([0 0;data2plot.compensations{s}.(vow)]);
                mag = 1;
                
                %training vowel angle
                tAngle = get_angle([0 0],data2plot.compensations{s}.(tVow));
                
                dAngle(s) = wrapTo2Pi(angle-tAngle);

                if dAngle(s) > pi
                    dAngle(s) = dAngle(s) - (2*pi);
                elseif dAngle < (-1 * pi)
                    dAngle(s) = 2*pi + dAngle(s);
                end
                
                x(s) = mag.*cos(dAngle(s));
                y(s) = mag.*sin(dAngle(s));
                h_q = quiver(0,0,x(s),y(s),0);
                set(h_q,'Color',plotParams.LineColor,'LineWidth',1,'MaxHeadSize',0.5)
            end
            for r = .25:.25:1
                rectangle('Position',[-r -r 2*r 2*r],'Curvature',[1,1],...
                    'EdgeColor',plotParams.LineColor,'FaceColor','none');
            end
            try
                mag = (1-circ_var(dAngle));
                x = mag.*cos(circ_mean(dAngle));
                y = mag.*sin(circ_mean(dAngle));
                h_q = quiver(0,0,x,y,0);
                set(h_q,'Color','k','LineWidth',3,'MaxHeadSize',0.5./mag)
            catch
            end
            if t ==1
                ylabel({vow,'F2'})
            else
                set(gca,'YTickLabel',[])
            end
            if v == 1
                title(tVow)
            end
            if v == nVowels
                xlabel('F1')
            else
%                 set(gca,'XTickLabel',[])
            end
            plot(0,0,'+k');
            axis square
            makeFig4Printing
        end
    end
    
    figpos_cm = [1 15 fullPageWidth fullPageWidth];
    tiledParams.TileSpacing = 'compact';
    h(ifig) = figure('Units','centimeters','Position',figpos_cm);
    copy_fig2tiledlayout(h4,h(ifig),nVowels,nTVowels,[],[],1,tiledParams);
end

%% Fig 5: adaptation angles relative to training vowels, radial plot version
[bPlot,ifig] = ismember(5,figs2plot);
if ~exist('circ_var','file')
    warning('circular statistics package not detected. no mean angle will be shown.')
    warning('https://github.com/circstat/circstat-matlab/blob/master/circ_var.m')
end
if bPlot
    h5 = [];
        
%     m = 257;
%     m1 = floor(m*0.5);
%     r = (0:m1-1)'/max(m1,1);
%     g = r;
%     r = [r; ones(m1+1,1)];
%     g = [g; 1; flipud(g)];
%     b = flipud(r);
%     rgb = [r g b]; 
%     

    rgb = flipud(brewermap(256,'RdBu'));
    nSubs = length(dataPaths);
    tVowels = {'ih','ey','eh','ah','ow','iy','ae','aa','uw'};
    nTVowels = length(tVowels);
    
    gVowels = tVowels;
    nGVowels = length(gVowels);
    
%     tVowels = {'iy','ae','aa','uw'};
%     nTVowels = length(tVowels);
%     
%     gVowels = {'ih','ey','eh','ah','ow'};
%     nGVowels = length(gVowels);
    
    allVowels = [tVowels gVowels];
    nVowels = length(allVowels);
    for s = 1:nSubs
        for v = 1:nVowels
            vow = allVowels{v};
            vDat.dists.(vow)(s) = data2plot.dists{s}.(vow);
            vDat.angs.(vow)(s) = get_angle([0 0],data2plot.compensations{s}.(vow));
        end
    end
    
    corrMatrix.dists = NaN(nGVowels,nTVowels);
    corrMatrix.ang= NaN(nGVowels,nTVowels);
    for v1 = 1:nGVowels
        vow1 = gVowels{v1};
        for v2 = 1:nTVowels
            vow2 = tVowels{v2};
            [corrMatrix.dists(v1,v2),p.dists(v1,v2)] = corr(vDat.dists.(vow1)',vDat.dists.(vow2)');
            [corrMatrix.angs(v1,v2),p.angs(v1,v2)] = circ_corrcc(vDat.angs.(vow1)',vDat.angs.(vow2)');
        end
        
    end
    
    titles = {'Adaptation magnitude (r)','Adaptation angle (r)'};
    measures = {'dists','angs'};
    for f = 1:2
        h5(f) = figure;
        plotMatrix = triu(corrMatrix.(measures{f}),1);
        cMapLength = size(rgb,1);
        image((plotMatrix+1)*cMapLength/2)
        
        for v1 = 1:nTVowels
            for v2 = 1:nGVowels
                if plotMatrix(v2,v1) ~=0
                    if p.(measures{f})(v2,v1)<=0.05
                        rText = sprintf('%.2f*',plotMatrix(v2,v1));
                        fWeight = 'bold';
                    else
                        rText = sprintf('%.2f',plotMatrix(v2,v1));
                        fWeight = 'normal';
                    end

                    text(v1,v2,rText,'HorizontalAlignment','center',...
                        'FontWeight',fWeight,'Color',[0 0 0])
                end
            end
        end
        set(gca,'XTick',1:nTVowels,'XTickLabel',arpabet2ipa(tVowels),'YTick',1:nGVowels,'YTickLabel',arpabet2ipa(gVowels))
        set(gca, 'XAxisLocation', 'top')
        xlim([1.5 9.5])
        ylim([0.5 8.5])
        h_l(1) = hline(5.5,'k','-');
        h_l(2) = vline(5.5,'k','-');
        for l = 1:length(h_l)
            set(h_l(l),'LineWidth',2)
        end
        title(titles{f})
        makeFig4Printing
    end
    
    figpos_cm = [1 15 fullPageWidth .5*fullPageWidth];
    tiledParams.TileSpacing = 'compact';
    h(ifig) = figure('Units','centimeters','Position',figpos_cm);
    copy_fig2tiledlayout(h5,h(ifig),1,2,[],[],1,tiledParams);
    colormap(rgb)
end


end

function phi = get_angle(loc1,loc2)
    phi = atan((loc2(2)-loc1(2))/(loc2(1)-loc1(1)));
    if (loc2(1)-loc1(1)) < 0
        phi = phi + pi;
    end
end