function NET=fig_multi_bars_release(varargin)

trainflag=0; logYflag=0; sparseXflag=0; medianflag=0; RMSEflag=0;meanRMSEflag=0;
floquetflag=0; floquetfracflag=0; dimflag=0;
Wnormflag=0; Wnorm2flag=0; condflag=0;
ERRoutlierflag=1; normflag=0;
saveflag=0;xlogflag=0;subspaceflag=0;stabmarginflag=0;sineprodflag=0;allsineflag=0;sparsestatsflag=0;
filenamestr='';
outlthr=1;
FMcutoff=1.0;
maxsubang=60;
alldata=[];
index=1;

%% read command line options
% varargin=varargin{1};
argidx = 1;
while argidx<=length(varargin)
    if strcmpi( varargin{argidx}, 'data' )
        ii=1;
        while( (argidx+ii<=length(varargin)) && (isstruct(varargin{argidx+ii})))
            alldata{ii} = varargin{argidx+ii};
            ii=ii+1;
        end
        argidx=argidx+ii;
    elseif strcmpi( varargin{argidx}, 'index' ) || strcmpi( varargin{argidx}, 'idx' )
        index=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'show training' ) || strcmpi( varargin{argidx}, 'training' ) || strcmpi( varargin{argidx}, 'train' )
        trainflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'log axis' )
        logYflag = 1;
        logYlims=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'show sparsity stats' ) || strcmpi( varargin{argidx}, 'show sparsity' ) ||...
            strcmpi( varargin{argidx}, 'sparsity' )
        sparsestatsflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show median' ) || strcmpi( varargin{argidx}, 'median' )
        medianflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show RMSE' ) || strcmpi( varargin{argidx}, 'RMSE' )
        RMSEflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show mean RMSE' ) || strcmpi( varargin{argidx}, 'meanRMSE' )
        meanRMSEflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show dim' ) || strcmpi( varargin{argidx}, 'dim' )
        dimflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show condition number' ) || strcmpi( varargin{argidx}, 'show cond' ) ||...
            strcmpi( varargin{argidx}, 'cond' ) || strcmpi( varargin{argidx}, 'CN' )
        condflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show norms' ) || strcmpi( varargin{argidx}, 'norms' )
        normflag = 1;
        argidx=argidx+1;
    elseif (strcmpi( varargin{argidx}, 'show Wnorm' ) || strcmpi( varargin{argidx}, 'Wnorm' ) )
        Wnormflag = 1;
        argidx=argidx+1;
    elseif (strcmpi( varargin{argidx}, 'show Wnorm  regularised' ) || strcmpi( varargin{argidx}, 'Wnorm reg' ) )
        Wnorm2flag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show Floquet' ) || strcmpi( varargin{argidx}, 'Floquet' )
        floquetflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show Floquet fraction' ) || strcmpi( varargin{argidx}, 'Floquet fraction' ) ||...
            strcmpi( varargin{argidx}, 'Floquet frac' )
        floquetfracflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'FMcutoff' ) || strcmpi( varargin{argidx}, 'FM cutoff' )
        FMcutoff=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'show subspace angle' ) || strcmpi( varargin{argidx}, 'subspace angles' ) ||...
            strcmpi( varargin{argidx}, 'subspace angle' ) || strcmpi( varargin{argidx}, 'subspace' ) ||...
            strcmpi( varargin{argidx}, 'subspace cos' )
        subspaceflag = 1;
        maxsubang=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'show sine product' ) || strcmpi( varargin{argidx}, 'sine product') ||...
            strcmpi( varargin{argidx}, 'sineproduct')
        sineprodflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show all sine product' ) || strcmpi( varargin{argidx}, 'all sine product') ||...
            strcmpi( varargin{argidx}, 'allsineproduct')
        allsineflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show stability margins' ) || strcmpi( varargin{argidx}, 'stability margins' ) ||...
            strcmpi( varargin{argidx}, 'stability' )
        stabmarginflag = 1;
        stabstr=varargin{argidx+1};
        xlinecutoff=varargin{argidx+2};
        if strcmpi(stabstr,'EVs') || strcmpi( stabstr, 'EV' )
            xlinecritflag=0;
        elseif strcmpi(stabstr,'error') || strcmpi( stabstr, 'RMSE' )
            xlinecritflag=1;
        elseif strcmpi(stabstr,'Floquet multiplier') || strcmpi( stabstr, 'Floquet' ) || strcmpi( stabstr, 'FM' )
            xlinecritflag=2;
        elseif strcmpi(stabstr,'Floquet multiplier median') || strcmpi( stabstr, 'Floquet median' ) || strcmpi( stabstr, 'FM median' )
            xlinecritflag=4;
        elseif strcmpi(stabstr,'Floquet fraction') || strcmpi( stabstr, 'Floquet frac' ) || strcmpi( stabstr, 'FM frac' )
            xlinecritflag=5;
        else
            error('Wrong stability margin string.');
        end
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'clean outliers' ) || strcmpi( varargin{argidx}, 'clean ERR outliers' )
        ERRoutlierflag = 1;
        outlthr=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'savefig' ) || strcmpi( varargin{argidx}, 'save' )
        filenamestr=varargin{argidx+1};
        saveflag = 1;
        argidx=argidx+2;
    else
        argidx=argidx+1;
    end
end

%% set parameters
NS=length(alldata);
%% extract learning performance, measures, stats etc.
for i=1:NS
    alldata{i}.parameters.ERRoutlierflag=ERRoutlierflag;
    alldata{i}.parameters.FMcutoff=FMcutoff;
    if(alldata{i}.parameters.flags.cortexonlyflag==1)
        NET{i}=extract_stats_v3 (alldata{i}.parameters,alldata{i}.parameters.Npara,alldata{i},'ctx',outlthr);
    else
        NET{i}=extract_stats_v3 (alldata{i}.parameters,alldata{i}.parameters.Npara,alldata{i},'cereb',outlthr);
    end
end

parameters=alldata{1}.parameters;
learnflag=parameters.flags.learnflag;
if(learnflag==100)
    X=unique(parameters.sigmaN);
    Xstr='\sigma_N';
    savestr='sigN';
    % replace sigmaN=0 with finite value in order to plot on log-scale
    X(1)=X(2)/10;
elseif(learnflag==101)
    X=unique(parameters.GCthreshold);
    Xstr='\theta_{GC}';
    savestr='theta';
elseif(learnflag==102)
    X=unique(parameters.Ncx);
    Xstr='N_{ctx}';
    savestr='Nctx';
elseif(learnflag==103)
    X=unique(parameters.Ngc);
    Xstr='N_{gc}';
    savestr='Ngc';
elseif(learnflag==105)
    X=unique(parameters.sigmaNgc);
    Xstr='\sigma_{Ngc}';
    savestr='sigNgc';
    % replace sigmaN=0 with finite value in order to plot on log-scale
    X(1)=X(2)/10;
elseif(learnflag==106)
    X=unique(parameters.Wgoc);
    Xstr='W_{goc}';
    savestr='Wgoc';
elseif(learnflag==107)
    X=unique(parameters.Ncomp);
    Xstr='Ncomp';
    savestr='Ncomp';
elseif(learnflag==108)
    X=unique(parameters.sigAmp);
    Xstr='Signal Amp';
    savestr='sigAmp';
    xlogflag=1;
elseif(learnflag==109)
    X=unique(parameters.lambda);
    Xstr='Fundamental freq. (Hz)';
    savestr='lambda';
    xlogflag=1;
elseif(learnflag==110)
    X=unique(parameters.Tlearn)*unique(parameters.lambda);
    Xstr='T_{learn} (cycles)';
    savestr='Tlearn';
elseif(learnflag==113)
    X=unique(parameters.g);
    Xstr='g';
    savestr='g';
elseif(learnflag==117)
    X=parameters.Ngc;
    Xstr='N_{gc}';
elseif(learnflag==119)
    X=parameters.Ampmod;
    Xstr='Ampmod';
elseif(learnflag==120)
    X=parameters.sigmaLS;
    Xstr='\sigma_{LS}';
    savestr='sigLS';
    xlogflag=1;
elseif(learnflag==123)
    X=parameters.fc_avrg;
    Xstr='Fraction of active GCs';
    savestr='fcoding';
elseif(learnflag==124)
    X=parameters.sp_alpha;
    Xstr='\alpha_{sp}';
    savestr='alphasp';
    xlogflag=1;
elseif(learnflag==128)
    X=parameters.rho_Wfb12;
    Xstr='atanh(\rho_{Wfb12})';
    savestr='rho12';
elseif(learnflag==129)
    X=unique(parameters.sigAmp_intf);
    Xstr='Intf. signal Amp';
    savestr='sigAmpintf';
    xlogflag=1;
end

Ncomp=max(parameters.Ncomp);

parameters.FMcutoff=FMcutoff;
interferenceflag=parameters.flags.interferenceflag;

%% filename
if(saveflag==1)
    % simname=split(data.parameters.simtype);
    if(interferenceflag>0)
        intfname=['_Nintf' num2str(parameters.Ncomp_intf) ];
    else
        intfname='';
    end
    filename=[savestr '_Ncomp' num2str(Ncomp) intfname '_' parameters.learningrule '_' filenamestr];
end

%% estimate stability margins using error/EV/Floquet
xintersect=[];
if(stabmarginflag==1)
    for i=1:NS
        if(xlinecritflag==1)
            xintersect{i} = findIntersect(X,NET{i}.RMSE_test,xlinecutoff*ones(numel(X),1),xlogflag,1); %#ok<*AGROW>
        elseif(xlinecritflag==2)
            xintersect{i} = findIntersect(X,NET{i}.FLOQUET_avrg,xlinecutoff*ones(numel(X),1),xlogflag,1);
        elseif(xlinecritflag==4)
            xintersect{i} = findIntersect(X,NET{i}.FLOQUET_median,xlinecutoff*ones(numel(X),1),xlogflag,1);
        elseif(xlinecritflag==5)
            xintersect{i} = findIntersect(X,NET{i}.FLOQUET_frac,xlinecutoff*ones(numel(X),1),xlogflag,0);
        end
    end
end

%% figures
cmp = get(groot,'DefaultAxesColorOrder');
ctxcolour=cmp(1,:);
cerebcolour=cmp(2,:);
% suppcolour=cmp(5,:);
for i=1:NS
    para=alldata{i}.parameters;
    if(para.flags.cortexonlyflag==1)
        clr{i}=ctxcolour;
        xTickstr{i}='RNN';
        if(para.flags.taudistrbflag==1)
            xTickstr{i}='hetRNN';
        end

    else
        clr{i}=cerebcolour;
        % xTickstr{i}='P+RNN';
        Expstr='';
        Gainstr='';
        GoCstr='';
        if(para.GCgain>1)
            Gainstr='+gain';
        end
        if(para.flags.taugcmodeflag==0)
            Pstr='FFL';
        elseif(para.flags.taugcmodeflag==1)
            Pstr='tauFFL';
        elseif(para.flags.taugcmodeflag==2)
            Pstr='hettauFFL';
        end
        if(para.flags.taudistrbflag==0)
            RNNstr='RNN';
        elseif(para.flags.taudistrbflag==1)
            RNNstr='hetRNN';
        end
        if(para.flags.gocmodeflag==0)
            GoCstr='noGoC ';
        end
        if((para.Ngc/para.Ncx) > 1)
            Expstr=['1:' num2str(para.Ngc/para.Ncx) ' '];
        end
        xTickstr{i}=[GoCstr Expstr RNNstr '+' Pstr Gainstr];
    end
end

FS=12;
LW=2;
MS=6;
set(0,'defaultAxesFontName', 'Arial');

Nflags=trainflag+dimflag+RMSEflag+meanRMSEflag+condflag...
    +sineprodflag+allsineflag+2*sparsestatsflag...
    +normflag+Wnormflag+Wnorm2flag+floquetflag...
    +floquetfracflag+subspaceflag+stabmarginflag;

sx=8;
sy=10;
nr=2;
if(saveflag==1)
    nr=1;
    sx=4.5;
    sy=sx;
    Fwidth=sx*Nflags/nr;
    tilestr='tight';
    paddstr='tight';
    FS=10;
    MS=4;
else
    Fwidth=sx+sx*Nflags/nr;
    tilestr='compact';
    paddstr='compact';
end

% plot parameters
plotpara.sparseXflag=sparseXflag;
plotpara.Xstr=Xstr;
plotpara.FS=FS;
plotpara.LW=LW;
plotpara.MS=MS;
plotpara.ctxcolour=ctxcolour;
plotpara.cerebcolour=cerebcolour;
plotpara.learnflag=learnflag;
plotpara.xintersect=xintersect;

if(Nflags>0)
    f1=figure;
    set(f1,'papertype','A4','Units', 'centimeters', 'position', [10,40,Fwidth,sy*nr]);
    tiledlayout(nr,ceil(Nflags/nr),'TileSpacing',tilestr,'Padding',paddstr);

    %% panels
    if(RMSEflag==1)
        nexttile;
        hold on;

        Y=[];
        Z=[];
        for i=1:NS
            if(medianflag)
                Y=[Y; NET{i}.RMSE_median_test(index)];
            else
                Y=[Y; NET{i}.RMSE_test(index)];
                Z=[Z; NET{i}.RMSE_SD_test(index)];
            end
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        errorbar(1:NS,Y,Z,Z,'LineStyle','none','Color','k');

        xlim([0.5 NS+0.5]);
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);

        ylabel('Test error (RMSE)','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);
        if(logYflag)
            set(gca, 'YScale', 'log');
            ylim([logYlims(1) logYlims(2)]);
            yticks(10.^(log10(logYlims(1)):log10(logYlims(end))));
        else
            ylim([0 1]);
        end
        title([Xstr ' = ' num2str(round(X(index),1))]); ax = gca; ax.TitleFontSizeMultiplier = 0.8;
    end

    if(meanRMSEflag==1)
        nexttile;
        hold on;

        Y=[];
        for i=1:NS
            Y=[Y; mean(NET{i}.RMSE_test)];
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        xlim([0.5 NS+0.5]);
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);

        ylabel('Mean error (meanRMSE)','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);
        if(logYflag)
            set(gca, 'YScale', 'log');
            ylim([0.1 1]);
            yticks(10.^(log10(0.1):log10(1)));
        else
            ylim([0 1]);
        end
        title([Xstr ' = ' num2str(round(X(index),1))]); ax = gca; ax.TitleFontSizeMultiplier = 0.8;
    end

    if(trainflag)
        nexttile;
        hold on;
        Y=[];
        for i=1:NS
            Y=[Y; NET{i}.RMSE_train(index)];
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        xlim([0.5 NS+0.5]);
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);
        ylabel('Training error','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);
        if(logYflag)
            set(gca, 'YScale', 'log');
            ylim([logYlims(1) logYlims(2)]);
            yticks(10.^(log10(logYlims(1)):log10(logYlims(end))));
        else
            ylim([0 1]);
        end

    end

    if(sparsestatsflag)
        nexttile;
        hold on;

        Y=[];
        for i=1:NS
            if(any("f_max" == string(NET{i}.STATS_avrg.Properties.VariableNames)))
                Y=[Y; NET{i}.STATS_avrg.f_mean(index)];
            else
                Y=[Y; nan];
            end
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        xlim([0.5 NS+0.5]);
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);

        ylim([0 1]);
        ylabel('f_{lt}','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);
        title('GC lifetime sparsity'); ax = gca; ax.TitleFontSizeMultiplier = 0.8;

        nexttile;
        hold on;

        Y=[];
        for i=1:NS
            if(any("fpop_max" == string(NET{i}.STATS_avrg.Properties.VariableNames)))
                Y=[Y; NET{i}.STATS_avrg.fpop_mean(index)];
            else
                Y=[Y; nan];
            end
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        xlim([0.5 NS+0.5]);
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);

        ylim([0 1]);
        ylabel('f_{pop}','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);
        title('GC population sparsity'); ax = gca; ax.TitleFontSizeMultiplier = 0.8;
    end

    if(normflag)
        nexttile;
        hold on;
        idx=1:2:(2*Ncomp);
        Y=[];
        for i=1:NS
            % if(alldata{i}.parameters.flags.cortexonlyflag==1)
            %     clr{i}=ctxcolour;
            % else
            %     clr{i}=cerebcolour;
            % end
            Y=[Y; mean(NET{i}.NORMS_avrg(index,idx+1),2)];
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        xlim([0.5 NS+0.5]);
        set(gca, 'YScale', 'log');
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);
        ylim([0.05 1]);
        ylabel('Norm (a.u.)','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);
    end

    if(condflag)
        nexttile;
        hold on;

        Y=[];
        for i=1:NS
            Y=[Y; NET{i}.COND_avrg(index)];
            % Y=[Y; NET{i}.COND_TOTAL_avrg(index)];
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        xlim([0.5 NS+0.5]);
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);

        ylim([1 inf]);
        set(gca, 'YScale', 'log');
        ylabel('CN','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);

    end

    if(dimflag)
        % errorbar(X,CTX.DIM_avrg,CTX.DIM_std,'o--','Color',ctxcolour,'Markerfacecolor',ctxcolour);
        % errorbar(X,CEREB.DIM_avrg,CEREB.DIM_std,'o--','Color',cerebcolour,'Markerfacecolor',cerebcolour);
        nexttile;
        hold on;

        Y=[];
        for i=1:NS
            Y=[Y; NET{i}.DIM_FT_avrg(index)];
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        xlim([0.5 NS+0.5]);
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);

        ylim([1 inf]);
        ylabel('Dimensionality','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);

    end

    if(subspaceflag)
        nexttile;
        hold on;

        Y=[];
        for i=1:NS
            Y=[Y; mean(NET{i}.SUBSPACE_ANGLES_avrg(index,:))];
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        xlim([0.5 NS+0.5]);
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);

        ylim([0 maxsubang]);
        ylabel('Avrg. subsp. angle (deg.)','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);
    end

    if(allsineflag)
        nexttile;
        hold on;

        Y=[];
        Z=[];
        for i=1:NS
            % Y=[Y; sqrt(NET{i}.SINEPROD_SIGNAL_avrg(index).*prod(sind(NET{i}.SUBSPACE_ANGLES_avrg(index,:)),2))];
            Y=[Y; NET{i}.ALLSINE_avrg(index)];
            Z=[Z; NET{i}.ALLSINE_std(index)];
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        errorbar(1:NS,Y,Z,Z,'LineStyle','none','Color','k');

        xlim([0.5 NS+0.5]);
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);

        ylim([1e-4 2e-2]);
        ylabel('Sine product','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);
        set(gca, 'YScale', 'log');
    end

    if(Wnormflag)
        nexttile;
        hold on;

        Y=[];
        for i=1:NS
            Y=[Y; NET{i}.WNORM_avrg(index)];
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        xlim([0.5 NS+0.5]);
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);

        % ylim([1e-4 1]);
        ylabel('||W_{out}||','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);
        set(gca, 'YScale', 'log');
    end

    if(Wnorm2flag)
        % plot(X,NET.WNORM_DENOM_NONORMS_avrg,'--','Color',ctxcolour*0.7,'Markerfacecolor',ctxcolour*0.7);
        % plot(X,CEREB.WNORM_DENOM_NONORMS_avrg,'--','Color',cerebcolour*0.7,'Markerfacecolor',cerebcolour*0.7);

        nexttile;
        hold on;

        Y=[];
        for i=1:NS
            Y=[Y; NET{i}.WNORM_DENOM_avrg(index)];
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        xlim([0.5 NS+0.5]);
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);

        ylabel('den(||W_{out}||)','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);
        set(gca, 'YScale', 'log');
    end

    if(floquetflag)
        nexttile;
        hold on;

        Y=[];
        for i=1:NS
            if(medianflag==1)
                Y=[Y; NET{i}.FLOQUET_median(index)];
            else
                Y=[Y; NET{i}.FLOQUET_avrg(index)];
            end
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        xlim([0.5 NS+0.5]);
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);

        yline(1,'k--');
        ylim([0.1 1e6]);
        ylabel('Max. Abs. FM','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);
        set(gca, 'YScale', 'log');
        yticks(10.^(log10(1):2:log10(10^6)));
    end

    if(floquetfracflag)
        nexttile;
        hold on;

        Y=[];
        for i=1:NS
            Y=[Y; NET{i}.FLOQUET_frac(index)];
        end
        b=bar(1:NS,Y,'LineStyle','none','FaceColor','flat');
        for i=1:NS
            b.CData(i,:) = clr{i};
        end
        xlim([0.5 NS+0.5]);
        set(gca,'xTick',1:NS);
        set(gca,'XTickLabel',xTickstr);
        set(gca,'XTickLabelRotation',45);

        ylim([0 1]);
        ylabel(['Fraction FM > ' num2str(FMcutoff)],'FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,plotpara);
    end

    if(stabmarginflag==1)
        nexttile;
        hold on;

        % Y=[];G=[];
        % for i=1:NS
        %     if(isscalar(xintersect{i}))
        %         xintersect{i}=[xintersect{i}; max(X)];
        %     end
        %     if(length(xintersect{i})>2)
        %         xintersect{i}=[min(xintersect{i}); max(xintersect{i})];
        %     end
        %     Y=[Y; xintersect{i}];
        %     G=[G; repmat({xTickstr{i}},length(xintersect{i}),1)]; %#ok<CCAT1>
        % end
        % boxplot(Y,G,'Symbol','','Whisker',0,'BoxStyle','filled','MedianStyle','target','Orientation','horizontal');
        %
        % a = get(get(gca,'children'),'children');   % Get the handles of all the objects
        % t = get(a,'tag');               % List the names of all the objects
        % idx=strcmpi(t,'box');           % Find Box objects
        % boxes=a(idx);                   % Get the children you need
        % set(boxes,'linewidth',10);       % Set width
        % for i=1:NS
        %     boxes(NS-i+1).Color = clr{i};
        % end

        for i=1:NS
            y_vert = [0 0 0.3 0.3]+(i-1) +1;
            if(isscalar(xintersect{i}))
                x_vert = [xintersect{i}(1) X(end) X(end) xintersect{i}(1)];
                patch(x_vert,y_vert,clr{i},'EdgeColor','none');
            else
                for k=1:(numel(xintersect{i})-1)
                    x_vert = [xintersect{i}(k) plotpara.xintersect{i}(k+1) xintersect{i}(k+1) xintersect{i}(k)];
                    patch(x_vert,y_vert,clr{i},'EdgeColor','none');
                end
            end
        end


        ylim([0.5 NS+0.5]);
        set(gca,'yTick',1:NS);
        set(gca,'yTickLabel',xTickstr);

        xlabel('Frequency (Hz)','FontSize',FS);
        axis(gca, 'square'); set(gca, 'TickDir', 'out');
        box off; set(gca,'FontSize',FS); set(gca, 'Color', 'None');
        if(xlogflag)
            set(gca, 'XScale', 'log');
            xticks(10.^(log10(0.1):log10(100)));
        end
        xlim([min(X) 100]);
    end

    if(saveflag==1)
        filenamefull=['temp_figures/bar_' filename];
        saveas(f1,filenamefull,'svg');
    end

end

end

%% functions
function OUT=extract_stats_v3 (parameters,Npara,data,typestring,outlthr)
Nreal=parameters.Nreal;     % number of repititions / network realisations
Ncomp=parameters.Ncomp(end);

ERRoutlierflag=parameters.ERRoutlierflag;

% get relevant data
RMSE_test=Get_structs(data,typestring,'RMSE_test',parameters,1,Npara);
RMSE_train=Get_structs(data,typestring,'RMSE_train',parameters,1,Npara);
ANGLES=Get_structs(data,typestring,'Angles',parameters,100,Npara);
ANGLES_RAW=ANGLES;
NORMS=Get_structs(data,typestring,'Norms',parameters,20,Npara);
% DIM=Get_structs(data,typestring,'DIM',parameters,1,Npara);
DIM_FT=Get_structs(data,typestring,'DIM_FT',parameters,1,Npara);
COND=Get_structs(data,typestring,'COND',parameters,1,Npara);
COND_TOTAL=Get_structs(data,typestring,'COND_total',parameters,1,Npara);
FLOQUET=Get_structs(data,typestring,'Max_Floquet_multiplier',parameters,1,Npara);
SUBSPACE_ANGLES=Get_structs(data,typestring,'Subspaceangles_cos',parameters,1,Npara);
SINEPROD_SIGNAL=Get_structs(data,typestring,'Sineprod_signal',parameters,1,Npara);

WNORM=Get_structs(data,typestring,'Wnorm',parameters,1,Npara);

STATS=data.('stats');

OUT.RMSE_distrb=zeros(Npara,Nreal);

OUT.RMSE_train=zeros(Npara,1);
OUT.RMSE_test=zeros(Npara,1);
OUT.RMSE_SD_test=zeros(Npara,1);
OUT.RMSE_median_test=zeros(Npara,1);

% OUT.DIM_avrg=zeros(Npara,size(DIM,2));
% OUT.DIM_std=zeros(Npara,size(DIM,2));
OUT.DIM_FT_avrg=zeros(Npara,size(DIM_FT,2));
OUT.DIM_FT_std=zeros(Npara,size(DIM_FT,2));

OUT.COND_avrg=zeros(Npara,size(COND,2));
OUT.COND_std=zeros(Npara,size(COND,2));
OUT.COND_TOTAL_avrg=zeros(Npara,size(COND_TOTAL,2));
OUT.COND_TOTAL_std=zeros(Npara,size(COND_TOTAL,2));

OUT.FLOQUET_distrb=zeros(Npara,Nreal);
OUT.FLOQUET_avrg=zeros(Npara,size(FLOQUET,2));
OUT.FLOQUET_std=zeros(Npara,size(FLOQUET,2));
OUT.FLOQUET_median=zeros(Npara,size(FLOQUET,2));
OUT.FLOQUET_frac=zeros(Npara,1);

OUT.SUBSPACE_ANGLES_avrg=zeros(Npara,size(SUBSPACE_ANGLES{1},2));
OUT.SUBSPACE_ANGLES_std=zeros(Npara,size(SUBSPACE_ANGLES{1},2));

maxsize_angles=max(cell2mat(cellfun(@size,ANGLES,'UniformOutput',false)));
for i=1:size(ANGLES,1)
    % pad with NaNs
    [m, n, k]=size(ANGLES{i});
    ANGLES{i}=[ANGLES{i} NaN(m,maxsize_angles(2)-n,k); NaN(maxsize_angles(1)-m,maxsize_angles(2),k)];

    ANGLES_RAW{i}=ANGLES{i};
    % map angles into 0-90 degree range
    ANGLES{i}(ANGLES{i}>90)=180-ANGLES{i}(ANGLES{i}>90);

    OUT.ANGLES_avrg=zeros(Npara,size(ANGLES{i},1),size(ANGLES{i},2));
end

maxsize_norms=max(cell2mat(cellfun(@size,NORMS,'UniformOutput',false)));
for i=1:size(NORMS,1)
    % pad with NaNs
    m=size(NORMS{i},2);
    NORMS{i}=[NORMS{i} NaN(maxsize_norms(1),maxsize_norms(2)-m)];

    OUT.NORMS_avrg=zeros(Npara,size(NORMS{i},2));
    OUT.NORMS_std=zeros(Npara,size(NORMS{i},2));
    OUT.NORMSratio_avrg=zeros(Npara,(Ncomp^2 -Ncomp)/2);
end

OUT.WNORM_avrg=zeros(Npara,1);
OUT.WNORM_std=zeros(Npara,1);
OUT.WNORM_NONORMS_avrg=zeros(Npara,1);
OUT.WNORM_DENOM_avrg=zeros(Npara,1);
OUT.WNORM_DENOM_NONORMS_avrg=zeros(Npara,1);

OUT.SINEPROD_SIGNAL_avrg=zeros(Npara,1);
OUT.ALLSINE_avrg=zeros(Npara,1);
OUT.ALLSINE_std=zeros(Npara,1);

OUT.success=zeros(Npara,1);

OUT.STATS_avrg=[];
OUT.STATS_std=[];

for i=1:Npara

    %% correct outliers and NaN values
    if(ERRoutlierflag)
        RMSE_test{i}(RMSE_test{i}>outlthr)=outlthr;
        RMSE_test{i}(isnan(RMSE_test{i}))=outlthr;
    end

    %% full RMSE distributions
    OUT.RMSE_distrb(i,:)=RMSE_test{i}(:,1);

    %% mean, SD and median of RMSE
    OUT.RMSE_train(i)=mean(RMSE_train{i},'all','omitnan');
    OUT.RMSE_SD_train(i)=std(RMSE_train{i},0,'all','omitnan')/sqrt(Nreal);

    OUT.RMSE_test(i)=mean(RMSE_test{i},'all','omitnan');
    OUT.RMSE_SD_test(i)=std(RMSE_test{i},0,'all','omitnan')/sqrt(Nreal);

    OUT.RMSE_median_test(i)=median(RMSE_test{i},'all','omitnan');

    %% dimensionality
    % OUT.DIM_avrg(i,:)=mean(DIM{i},'omitnan');
    % OUT.DIM_std(i,:)=std(DIM{i},'omitnan')/sqrt(Nreal);

    %% FT dimensionality
    OUT.DIM_FT_avrg(i,:)=mean(DIM_FT{i},'omitnan');
    OUT.DIM_FT_std(i,:)=std(DIM_FT{i},'omitnan')/sqrt(Nreal);

    %% Condition number
    OUT.COND_avrg(i,:)=mean(COND{i},'omitnan');
    OUT.COND_std(i,:)=std(COND{i},'omitnan')/sqrt(Nreal);
    OUT.COND_TOTAL_avrg(i,:)=mean(COND_TOTAL{i},'omitnan');
    OUT.COND_TOTAL_std(i,:)=std(COND_TOTAL{i},'omitnan')/sqrt(Nreal);

    %% full Floquet multiplier distributions
    OUT.FLOQUET_distrb(i,:)=FLOQUET{i}(:,1);

    %% maximum Floquet multiplier
    OUT.FLOQUET_avrg(i,:)=mean(FLOQUET{i},'omitnan');
    OUT.FLOQUET_std(i,:)=std(FLOQUET{i},'omitnan')/sqrt(Nreal);
    OUT.FLOQUET_median(i,:)=median(FLOQUET{i},'all','omitnan');

    %% fraction of maximum Floquet multiplier > threshold
    OUT.FLOQUET_frac(i)=mean(FLOQUET{i}(:,1)>parameters.FMcutoff);

    %% Subspace angles
    OUT.SUBSPACE_ANGLES_avrg(i,:)=mean(SUBSPACE_ANGLES{i},'omitnan');
    OUT.SUBSPACE_ANGLES_std(i,:)=std(SUBSPACE_ANGLES{i},'omitnan')/sqrt(Nreal);

    %% spanning vector angles
    OUT.ANGLES_avrg(i,:,:)=mean(ANGLES{i},3,'omitnan');

    %% spanning vector norms
    OUT.NORMS_avrg(i,:)=mean(NORMS{i},'omitnan');
    OUT.NORMS_std(i,:)=std(NORMS{i},'omitnan')/sqrt(Nreal);

    %% ration of signal spanning vector norms
    tmp=OUT.NORMS_avrg(i,2:2:(2*Ncomp+1));
    tmp=tmp'*(1./tmp); tmp=triu(tmp,1); tmp=tmp(:); tmp(tmp==0)=[];
    OUT.NORMSratio_avrg(i,:)=tmp;

    %% readout weight norms
    OUT.WNORM_avrg(i)=mean(WNORM{i},'omitnan');
    OUT.WNORM_std(i)=std(WNORM{i},'omitnan')/sqrt(Nreal);
    OUT.WNORM_NONORMS_avrg(i)=OUT.WNORM_avrg(i)*sqrt(sum(OUT.NORMS_avrg(i,2:2:(2*Ncomp+1)).^2))*sqrt(parameters.Ncx);
    % OUT.WNORM_NONORMS_avrg(i)=OUT.WNORM_avrg(i)*nthroot(prod(OUT.NORMS_avrg(i,2:2:(2*Ncomp+1))),length(2:2:(2*Ncomp+1)))*sqrt(parameters.Ncx);
    % OUT.WNORM_NONORMS_avrg(i)=OUT.WNORM_avrg(i)*prod(OUT.NORMS_avrg(i,2:2:(2*Ncomp+1)))*sqrt(parameters.Ncx);

    %% readout weight norm denominator
    signal_angles=[];
    for k=1:2:(2*Ncomp)
        for j=1:2:(k-1)
            signal_angles=[signal_angles, sind(OUT.ANGLES_avrg(i,k+1,j+1))];
        end
    end
    OUT.SINEPROD_SIGNAL_avrg2(i)=prod(signal_angles);
    OUT.SINEPROD_SIGNAL_avrg(i)=mean(SINEPROD_SIGNAL{i},'omitnan');
    OUT.SINEPROD_SIGNAL_std(i)=std(SINEPROD_SIGNAL{i},'omitnan')/sqrt(Nreal);

    OUT.WNORM_DENOM_avrg(i)=1/(prod(sind(OUT.SUBSPACE_ANGLES_avrg(i,:)))*OUT.SINEPROD_SIGNAL_avrg(i)*prod(OUT.NORMS_avrg(i,2:2:(2*Ncomp+1)))*sqrt(parameters.Ncx));
    OUT.WNORM_DENOM_NONORMS_avrg(i)=OUT.WNORM_DENOM_avrg(i) *sqrt(sum(OUT.NORMS_avrg(i,2:2:(2*Ncomp+1)).^2))*sqrt(parameters.Ncx);
    % OUT.WNORM_DENOM_NONORMS_avrg(i)=1/(prod(sind(OUT.SUBSPACE_ANGLES_avrg(i,:)))*prod(signal_angles));
    % OUT.WNORM_DENOM_avrg(i)=1/prod(sind(OUT.SUBSPACE_COS_ANGLES_avrg(i,:)));

    %% allsine
    % OUT.ALLSINE_avrg(i)=sqrt(OUT.SINEPROD_SIGNAL_avrg(i).*prod(sind(OUT.SUBSPACE_ANGLES_avrg(i,:)),2));

    for rr=1:Nreal
        sub_ang=sind(SUBSPACE_ANGLES{i}(rr,:));
        sig_ang=[];
        for k=1:2:(2*Ncomp)
            for j=1:2:(k-1)
                sig_ang=[sig_ang, sind(ANGLES{i}(k+1,j+1,rr))];
            end
        end
        prod_ang(rr)=prod([sub_ang,sig_ang]);
    end
    OUT.ALLSINE_avrg(i)=mean(prod_ang,'omitnan');
    OUT.ALLSINE_std(i)=std(prod_ang,'omitnan')/sqrt(Nreal);
    %% CONTINUE HERE

    %% general stats
    tmp= varfun(@nanmean, STATS{i},"InputVariables",STATS{i}.Properties.VariableNames); %#ok<NANMEAN>
    tmp.Properties.VariableNames = STATS{i}.Properties.VariableNames;
    OUT.STATS_avrg=[OUT.STATS_avrg; tmp];

    tmp= varfun(@nanstd, STATS{i},"InputVariables",STATS{i}.Properties.VariableNames); %#ok<NANSTD>
    tmp.Properties.VariableNames = STATS{i}.Properties.VariableNames;
    OUT.STATS_std=[OUT.STATS_std; tmp];

end

end

function Set_CCxC_scan_FigureProperties(ah,plotpara)
% ah = axis handle

line = get(ah,'Children');
for i=1:length(line)
    set(line(i),'Linewidth', plotpara.LW);
    if(strcmpi(line(i).Type,'line'))
        set(line(i),'MarkerEdgeColor', 'none');
        set(line(i),'MarkerSize', plotpara.MS);
    end
end
axis(ah, 'square');
set(ah, 'TickDir', 'out');
box off;
set(ah,'FontSize',plotpara.FS);
set(ah, 'Color', 'None');

end

function outdata=Get_structs(data,typestring,namestring,parameters,Narray2,Npara)
if(isfield(data,[namestring '_' typestring]))
    outdata=data.([namestring '_' typestring])(1:Npara);
else
    for i=1:Npara
        outdata{i,1}=nan([parameters.Nreal,Narray2]);
    end
end
end