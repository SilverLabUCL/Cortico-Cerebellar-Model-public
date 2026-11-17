function [NET, X]=fig_multi_scan_release(varargin)

trainflag=0; logYflag=0;gainflag=0; medianflag=0; RMSEflag=0;
floquetflag=0; floquetfracflag=0; dimflag=0;
Wnormflag=0; condflag=0; colourflag=0; averageCXflag=0;
legendflag=0; ERRoutlierflag=0; sNflag=0; normflag=0;
maxsineprodflag=0;saveflag=0;xlogflag=0;pranglesflag=0;stabmarginflag=0;indvsineprodflag=0;fullsineprodflag=0;
compactflag=0;labelsflag=0;
filenamestr='';
outlthr=inf;
FMcutoff=1.0;
maxprang=60;
alldata=[];

%% read command line options
argidx = 1;
while argidx<=length(varargin)
    if strcmpi( varargin{argidx}, 'data' )
        ii=1;
        while( (argidx+ii<=length(varargin)) && (isstruct(varargin{argidx+ii})))
            alldata{ii} = varargin{argidx+ii};
            ii=ii+1;
        end
        argidx=argidx+ii;
    elseif strcmpi( varargin{argidx}, 'show training' ) || strcmpi( varargin{argidx}, 'training' )
        trainflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'log axis' )
        logYflag = 1;
        logYlims=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'show gain' ) || strcmpi( varargin{argidx}, 'gain' )
        gainflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show median' ) || strcmpi( varargin{argidx}, 'median' )
        medianflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show RMSE' ) || strcmpi( varargin{argidx}, 'RMSE' )
        RMSEflag = 1;
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
    elseif strcmpi( varargin{argidx}, 'show Floquet' ) || strcmpi( varargin{argidx}, 'Floquet' )
        floquetflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show Floquet fraction' ) || strcmpi( varargin{argidx}, 'Floquet fraction' ) ||...
            strcmpi( varargin{argidx}, 'Floquet frac' ) ||  strcmpi( varargin{argidx}, 'FM frac' )
        floquetfracflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'FMcutoff' ) || strcmpi( varargin{argidx}, 'FM cutoff' )
        FMcutoff=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'show principal angles' ) || strcmpi( varargin{argidx}, 'principal angles' ) ||...
            strcmpi( varargin{argidx}, 'principal angle' ) || strcmpi( varargin{argidx}, 'principal' )
        pranglesflag = 1;
        maxprang=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'show individual sine products' ) || strcmpi( varargin{argidx}, 'individual sine product') ||...
            strcmpi( varargin{argidx}, 'individual sineproduct') || strcmpi( varargin{argidx}, 'individual sine')
        indvsineprodflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show full sine product' ) || strcmpi( varargin{argidx}, 'show sine product' ) ||...
            strcmpi( varargin{argidx}, 'sine product') || strcmpi( varargin{argidx}, 'sineproduct')
        fullsineprodflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show legends' ) || strcmpi( varargin{argidx}, 'show legend' ) || strcmpi( varargin{argidx}, 'legend')
        legendflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show maximum sine product' ) || strcmpi( varargin{argidx}, 'show max sine prod' )
        maxsineprodflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show stability margins' ) || strcmpi( varargin{argidx}, 'stability margins' ) ||...
            strcmpi( varargin{argidx}, 'stability' )
        stabmarginflag = 1;
        stabstr=varargin{argidx+1};
        xlinecutoff=varargin{argidx+2};
        if strcmpi(stabstr,'error') || strcmpi( stabstr, 'RMSE' )
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
    elseif strcmpi( varargin{argidx}, 'compact' )
        compactflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show labels' ) || strcmpi( varargin{argidx}, 'labels' )
        labelsflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'clean outliers' ) || strcmpi( varargin{argidx}, 'clean ERR outliers' )
        ERRoutlierflag = 1;
        outlthr=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'average cortex' ) || strcmpi( varargin{argidx}, 'avrg cortex' ) ||...
            strcmpi( varargin{argidx}, 'average ctx' ) || strcmpi( varargin{argidx}, 'avrg ctx' )
        averageCXflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'colour' )
        colour=varargin{argidx+1};
        if(~isnumeric(colour))
            error('Missing colour multiplier.');
        end
        colourflag = 1;
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
nctx=0;
for i=1:NS
    if(length(alldata{i}.simcount)<alldata{i}.parameters.Nreal*alldata{i}.parameters.Npara)
        Npara(i)=size(alldata{i}.time_test,1)-1;
    else
        Npara(i)=alldata{i}.parameters.Npara;
    end

    alldata{i}.parameters.ERRoutlierflag=ERRoutlierflag; %#ok<*AGROW>
    alldata{i}.parameters.FMcutoff=FMcutoff; 
    if(alldata{i}.parameters.flags.cortexonlyflag==1)
        nctx=nctx+1;
        NET{i}=extract_stats_v3 (alldata{i}.parameters,Npara(i),alldata{i},'ctx',outlthr); 
    else
        NET{i}=extract_stats_v3 (alldata{i}.parameters,Npara(i),alldata{i},'cereb',outlthr);
    end
end

%% set parameters
parameters=alldata{1}.parameters;
learnflag=parameters.flags.learnflag;
if(learnflag==100)
    X=unique(parameters.sigmaN);
    Xstr='\sigma_N';
    savestr='sigN';
    sNflag=1;
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
    sNflag=1;
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
elseif(learnflag==114)
    X=unique(parameters.eta);
    Xstr='eta';
    xlogflag=1;
elseif(learnflag==117)
    X=parameters.Ngc;
    Xstr='N_{gc}';
elseif(learnflag==118)
    X=parameters.Phi1;
    Xstr='\phi';
    savestr='Phi';
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
elseif(learnflag==131)
    X=unique(parameters.GCgain);
    Xstr='GC gain';
    xlogflag=1;
elseif(learnflag==132)
    X=parameters.rho_Wfb13;
    Xstr='atanh(\rho_{Wfb13})';
    savestr='rho13';
end

[Npara, minidx]=min(Npara);
Nreal=parameters.Nreal;
Ncomp=max(parameters.Ncomp);

parameters.ERRoutlierflag=ERRoutlierflag;
parameters.FMcutoff=FMcutoff;
interferenceflag=parameters.flags.interferenceflag;

if(parameters.flags.floquetflag==0)
    stabmarginflag=0;
    floquetflag=0;
    floquetfracflag=0;
end

if(length(alldata{minidx}.simcount)<Nreal*Npara)
    Npara=size(alldata{minidx}.time_test,1)-1;
    X=X(1:Npara);
end

%% filename
if(saveflag==1)
    if(interferenceflag>0)
        intfname=['_Nintf' num2str(parameters.Ncomp_intf) ];
    else
        intfname='';
    end
    RNNstr=['_' parameters.activation_function];
    filename=[savestr RNNstr '_Ncomp' num2str(Ncomp) intfname '_' parameters.learningrule '_' filenamestr];
end

%% average ctx data
if(averageCXflag)
    for i=1:NS
        if(alldata{i}.parameters.flags.cortexonlyflag==1)
            NET{i}=average_ERROR(parameters,Npara,NET{i},alldata{i}.RMSE_test_ctx,...
                alldata{i}.Angles_ctx,alldata{i}.DIM_FT_ctx,alldata{i}.COND_ctx,alldata{i}.COND_total_ctx,...
                alldata{i}.Wnorm_ctx,alldata{i}.Max_Floquet_multiplier_ctx,alldata{i}.Subspaceangles_cos_ctx);
        end
    end
end

%% create labels for angles
labels_str=strings((2*Ncomp+10),(2*Ncomp+10));
for i=1:(2*Ncomp+10)
    for j=1:(i-1)
        if(i==1)
            str1='0';
        else
            if (mod(i,2)==0); istr='c'; else; istr='s'; end
            str1=[istr '^' num2str(ceil((i-1)/2))];
        end
        if(j==1)
            str2='0';
        else
            if (mod(j,2)==0); jstr='c'; else; jstr='s'; end
            str2=[jstr '^' num2str(ceil((j-1)/2))];
        end
        labels_str(i,j)=[str1, ' ', str2];
        labels_str(j,i)=labels_str(i,j);
    end
end
for i=1:NS
    NET{i}.labels_str=labels_str;
    NET{i}.labels_str=labels_str;
end

%% find idx of maximum sine product
if(maxsineprodflag==1)
    for i=1:NS
        Y=sqrt(NET{i}.SINEPROD_SIGNAL_avrg.*prod(sind(NET{i}.PRINCIPAL_ANGLES_avrg),2));
        [~,ymax_idx{i}]=max(Y);
        Xmax{i}=X(ymax_idx{i});
    end
end

%% estimate stability margins using error/EV/Floquet
cmp = get(groot,'DefaultAxesColorOrder');
ctxcolour=cmp(1,:);
cerebcolour=cmp(2,:);
xlinestyle='--';
xintersect=[];

stableidx=[];
if(stabmarginflag==1)
    for i=1:NS
        if(xlinecritflag==1)
            xintersect{i} = findIntersect(X,NET{i}.RMSE_test,xlinecutoff*ones(numel(X),1),xlogflag,1);
        elseif(xlinecritflag==2)
            xintersect{i} = findIntersect(X,NET{i}.FLOQUET_avrg,xlinecutoff*ones(numel(X),1),xlogflag,1);
        elseif(xlinecritflag==4)
            xintersect{i} = findIntersect(X,NET{i}.FLOQUET_median,xlinecutoff*ones(numel(X),1),xlogflag,1);
        elseif(xlinecritflag==5)
            xintersect{i} = findIntersect(X,NET{i}.FLOQUET_frac,xlinecutoff*ones(numel(X),1),xlogflag,0);
        end
        if(isscalar(xintersect{i}))
            stableidx{i}=X> min(xintersect{i});
        else
            stableidx{i}=(X> min(xintersect{i}) & X< max(xintersect{i}));
        end
        X_stable{i}=X(stableidx{i});
    end
else
    for i=1:NS
        X_stable{i}=X;
    end
end

%% figures
FS=14;
LW=1.5;
set(0,'defaultAxesFontName', 'Arial');
if(colourflag)
    ctxcolour=ctxcolour*colour;
    cerebcolour=cerebcolour*colour;
end
clrmod_ctx=flip(linspace(0.6,1,nctx));
clrmod=flip(linspace(0.6,1,NS-nctx));
ictx=1; icereb=1;
for i=1:NS
    if(alldata{i}.parameters.flags.cortexonlyflag==1)
        clr{i}=ctxcolour*clrmod_ctx(ictx);
        ictx=ictx+1;
    else
        clr{i}=cerebcolour*clrmod(icereb);
        icereb=icereb+1;
    end
end

if(compactflag==1)
    dimflag=0;
    condflag=0;
end

Nflags=trainflag+dimflag+gainflag+RMSEflag+condflag+indvsineprodflag+fullsineprodflag...
    +normflag+Wnormflag+floquetflag+floquetfracflag+pranglesflag+compactflag;

sx=8;
sy=10;
nr=2;
Fwidth=sx+sx*Nflags/nr;
tilestr='compact';
paddstr='compact';
if(saveflag==1)
    nr=1;
    sx=4.5;
    sy=sx;
    Fwidth=sx*Nflags/nr;
    tilestr='tight';
    paddstr='tight';
    FS=10;
end

% plot parameters
plotpara.Ncomp=Ncomp;
plotpara.xlogflag=xlogflag;
plotpara.sNflag=sNflag;
plotpara.Xstr=Xstr;
plotpara.FS=FS;
plotpara.LW=LW;
plotpara.ctxcolour=ctxcolour;
plotpara.cerebcolour=cerebcolour;
plotpara.learnflag=learnflag;
plotpara.stabmarginflag=stabmarginflag;
plotpara.xintersect=xintersect;
plotpara.labelsflag=labelsflag;
plotpara.xlinecol=clr;
plotpara.xlinestyle=xlinestyle;

if(Nflags>0)
    f1=figure;
    set(f1,'papertype','A4','Units', 'centimeters', 'position', [10,40,Fwidth,sy*nr]);
    tiledlayout(nr,ceil(Nflags/nr),'TileSpacing',tilestr,'Padding',paddstr);

    %% panels
    if(RMSEflag==1)
        nexttile;
        hold on;
        for i=1:NS
            plot(X,NET{i}.RMSE_test,':','color',clr{i});
            if(stabmarginflag==1)
                Y_stable=NET{i}.RMSE_test(stableidx{i});
                plot(X_stable{i},Y_stable,'-','color',clr{i});
            end
            if(medianflag)
                plot(X,NET{i}.RMSE_median_test,':','color',clr{i});
            end
            if(maxsineprodflag==1)
                plot(Xmax{i},NET{i}.RMSE_test(ymax_idx{i}),'ko','MarkerSize',3);
            end
        end
        ylabel('Test error (RMSE)','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        if(logYflag)
            set(gca, 'YScale', 'log');
            ylim([logYlims(1) logYlims(2)]);
            yticks(10.^(log10(logYlims(1)):log10(logYlims(end))));
        else
            ylim([0 1]);
        end
    end

    if(trainflag)
        nexttile;
        hold on;
        for i=1:NS
            plot(X,NET{i}.RMSE_train,'-','color',clr{i});
        end
        ylabel('Training error','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        if(logYflag)
            set(gca, 'YScale', 'log');
        end
    end

    if(gainflag)
        nexttile;
        for i=1:NS
            errorbar(X,NET{i}.STATS_avrg.gain_mean,NET{i}.STATS_std.gain_mean,'-','color',clr{i});
        end
        ylabel('Gain','FontSize',FS);
        set(gca, 'yScale', 'log');
        % ylim([]);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
    end

    if(normflag)
        nexttile;
        hold on;
        idx=1:2:(2*Ncomp);
        
        for i=1:NS
            Y=mean(NET{i}.NORMS_avrg(:,idx+1),2);
            plot(X,Y,':','Color',clr{i});
            if(stabmarginflag==1)
                Y_stable=Y(stableidx{i});
                plot(X_stable{i},Y_stable,'-','Color',clr{i});
            end
        end
        set(gca, 'YScale', 'log');
        ylim([-Inf 1]);
        ylabel('Norm (a.u.)','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
    end

    if(condflag)
        nexttile;
        hold on;
        for i=1:NS
            plot(X,NET{i}.COND_TOTAL_avrg,'-','color',clr{i});
        end

        ylim([1 inf]);
        set(gca, 'YScale', 'log');
        ylabel('CN','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        if(legendflag)
            legend('ctx','total ctx','gcs','total gcs','Location','Best');
        end
    end

    if(dimflag)
        nexttile;
        hold on;
        for i=1:NS
            plot(X,NET{i}.DIM_avrg,'-','color',clr{i});
        end

        ylim([1 inf]);
        ylabel('Dimensionality','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        if(legendflag)
            legend('ctx','gcs','signal ctx','signal gcs','Location','Best');
        end
    end

    if(compactflag)
        nexttile;
        yyaxis left
        hold on;
        for i=1:NS
            plot(X,NET{i}.COND_TOTAL_avrg,'-','color',clr{i});
        end

        % ylim([1 1e5]);
        ylim([1 1e11]);
        yticks([1 1e5 1e10]);
        set(gca, 'YScale', 'log');
        ylabel('CN','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);

        yyaxis right
        for i=1:NS
            plot(X,NET{i}.DIM_avrg,'-','color',clr{i}*0.5);
        end

        ylim([1 2]);
        yticks([1 2]);
        ylabel('Total dim.','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);

        ax = gca();
        ax.YAxis(1).Color = clr{1};
        ax.YAxis(2).Color = clr{1}*0.5;
    end

    if(pranglesflag)
        nexttile;
        hold on;
        for i=1:NS
            plot(X,NET{i}.PRINCIPAL_ANGLES_avrg,'-','color',clr{i});
        end

        ylim([0 maxprang]);
        ylabel('Principal angle (deg.)','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
    end

    if(indvsineprodflag)
        nexttile;
        hold on;

        for i=1:NS
            % task-relevant separation
            plot(X,NET{i}.SINEPROD_SIGNAL_avrg,'-','color',clr{i});
            % task-relevant vs task-irrelevant separation
            plot(X,prod(sind(NET{i}.PRINCIPAL_ANGLES_avrg),2),'--','color',clr{i});
        end 

        ylim([1e-2 1]);
        ylabel('Indv. sine products','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        set(gca, 'YScale', 'log');
    end

    if(fullsineprodflag)
        nexttile;
        hold on;

        for i=1:NS
            Y=NET{i}.SINEPROD_SIGNAL_avrg.*prod(sind(NET{i}.PRINCIPAL_ANGLES_avrg),2);
            plot(X,Y,':','color',clr{i});
            if(stabmarginflag==1)
                Y_stable=Y(stableidx{i});
                plot(X_stable{i},Y_stable,'-','color',clr{i});
            end
        end
        

        ylim([1e-4 1]);
        ylabel('Sine product','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        set(gca, 'YScale', 'log');
    end

    if(Wnormflag)
        nexttile;
        hold on;
        for i=1:NS
            plot(X,NET{i}.WNORM_avrg,'-','color',clr{i});
        end

        ylabel('||W_{out}||','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        set(gca, 'YScale', 'log');
    end

    if(floquetflag)
        nexttile;
        hold on;
        for i=1:NS
            plot(X,NET{i}.FLOQUET_avrg,'-','color',clr{i});
            % plot(X,NET{i}.FLOQUET_median,'-','color',clr{i});
        end

        yline(1,'k--');
        ylim([0.01 1e5]);
        ylabel('Max. Abs. FM','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        set(gca, 'YScale', 'log');
        yticks(10.^(log10(0.01):2:log10(10^5)));
    end

    if(floquetfracflag)
        nexttile;
        hold on;
        for i=1:NS
            plot(X,NET{i}.FLOQUET_frac,'-','color',clr{i});
        end
        ylim([0 1]);
        ylabel(['Fraction FM > ' num2str(FMcutoff)],'FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
    end

    if(saveflag==1)
        filenamefull=['temp_figures/scan_' filename];
        saveas(f1,filenamefull,'svg');
    end

end

end

%% functions
function OUT=extract_stats_v3 (parameters,Npara,data,typestring,outlthr)
Nreal=parameters.Nreal;     % number of repititions / network realisations
% Ncomp=parameters.Ncomp(end);

ERRoutlierflag=parameters.ERRoutlierflag;

% get relevant data
RMSE_test=Get_structs(data,typestring,'RMSE_test',parameters,1,Npara);
RMSE_train=Get_structs(data,typestring,'RMSE_train',parameters,1,Npara);
ANGLES=Get_structs(data,typestring,'Angles',parameters,100,Npara);
ANGLES_RAW=ANGLES;
NORMS=Get_structs(data,typestring,'Norms',parameters,20,Npara);
R=Get_structs(data,typestring,'R',parameters,8,Npara);
DIM=Get_structs(data,typestring,'DIM',parameters,1,Npara);
DIM_FT=Get_structs(data,typestring,'DIM_FT',parameters,1,Npara);
COND=Get_structs(data,typestring,'COND',parameters,1,Npara);
COND_TOTAL=Get_structs(data,typestring,'COND_total',parameters,1,Npara);
FLOQUET=Get_structs(data,typestring,'Max_Floquet_multiplier',parameters,1,Npara);
PRINCIPAL_ANGLES=Get_structs(data,typestring,'Subspaceangles_cos',parameters,1,Npara);
SINEPROD_SIGNAL=Get_structs(data,typestring,'Sineprod_signal',parameters,1,Npara);
WNORM=Get_structs(data,typestring,'Wnorm',parameters,1,Npara);
STATS=data.('stats');

OUT.RMSE_train=zeros(Npara,1);
OUT.RMSE_test=zeros(Npara,1);
OUT.RMSE_SD_test=zeros(Npara,1);
OUT.RMSE_median_test=zeros(Npara,1);

OUT.DIM_avrg=zeros(Npara,size(DIM,2));
OUT.DIM_std=zeros(Npara,size(DIM,2));
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

OUT.PRINCIPAL_ANGLES_avrg=zeros(Npara,size(PRINCIPAL_ANGLES{1},2));
OUT.PRINCIPAL_ANGLES_std=zeros(Npara,size(PRINCIPAL_ANGLES{1},2));

OUT.R_avrg=zeros(Npara,size(R{1},1),size(R{1},2));

maxsize_angles=max(cell2mat(cellfun(@size,ANGLES,'UniformOutput',false)));
for i=1:size(ANGLES,1)
    % pad with NaNs
    [m, n, k]=size(ANGLES{i});
    ANGLES{i}=[ANGLES{i} NaN(m,maxsize_angles(2)-n,k); NaN(maxsize_angles(1)-m,maxsize_angles(2),k)];
    for j=1:size(ANGLES{i},1)
        for k=1:size(ANGLES{i},2)
            ANGLES{i}(j,k,:)=ANGLES{i}(k,j,:);
        end
    end

    ANGLES_RAW{i}=ANGLES{i};
    % map angles into 0-90 degree range
    ANGLES{i}(ANGLES{i}>90)=180-ANGLES{i}(ANGLES{i}>90);

    OUT.ANGLES_avrg=zeros(Npara,size(ANGLES{i},1),size(ANGLES{i},2));
    OUT.ANGLES_std=zeros(Npara,size(ANGLES{i},1),size(ANGLES{i},2));
end

maxsize_norms=max(cell2mat(cellfun(@size,NORMS,'UniformOutput',false)));
for i=1:size(NORMS,1)
    % pad with NaNs
    m=size(NORMS{i},2);
    NORMS{i}=[NORMS{i} NaN(maxsize_norms(1),maxsize_norms(2)-m)];

    OUT.NORMS_avrg=zeros(Npara,size(NORMS{i},2));
    OUT.NORMS_std=zeros(Npara,size(NORMS{i},2));
end

OUT.WNORM_avrg=zeros(Npara,1);
OUT.WNORM_std=zeros(Npara,1);

OUT.SINEPROD_SIGNAL_avrg=zeros(Npara,1);

OUT.STATS_avrg=[];
OUT.STATS_std=[];

for i=1:Npara

    %% correct outliers and NaN values
    if(ERRoutlierflag)
        RMSE_test{i}(RMSE_test{i}>outlthr)=outlthr;
        RMSE_test{i}(isnan(RMSE_test{i}))=outlthr;
    end

    %% mean, SD and median of RMSE
    OUT.RMSE_train(i)=mean(RMSE_train{i},'all','omitnan');
    OUT.RMSE_SD_train(i)=std(RMSE_train{i},0,'all','omitnan')/sqrt(Nreal);

    OUT.RMSE_test(i)=mean(RMSE_test{i},'all','omitnan');
    OUT.RMSE_SD_test(i)=std(RMSE_test{i},0,'all','omitnan')/sqrt(Nreal);

    OUT.RMSE_median_test(i)=median(RMSE_test{i},'all','omitnan');

    %% dimensionality
    OUT.DIM_avrg(i,:)=mean(DIM{i},'omitnan');
    OUT.DIM_std(i,:)=std(DIM{i},'omitnan')/sqrt(Nreal);

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
    OUT.PRINCIPAL_ANGLES_avrg(i,:)=mean(PRINCIPAL_ANGLES{i},'omitnan');
    OUT.PRINCIPAL_ANGLES_std(i,:)=std(PRINCIPAL_ANGLES{i},'omitnan')/sqrt(Nreal);
    
    %% spanning vector angles
    OUT.ANGLES_avrg(i,:,:)=mean(ANGLES{i},3,'omitnan');
    OUT.ANGLES_std(i,:,:)=std(ANGLES{i},0,3,'omitnan')/sqrt(Nreal);

    %% spanning vector norms
    OUT.NORMS_avrg(i,:)=mean(NORMS{i},'omitnan');
    OUT.NORMS_std(i,:)=std(NORMS{i},'omitnan')/sqrt(Nreal);

    %% correlation between spanning vector
    OUT.R_avrg(i,:,:)=mean(R{i},3,'omitnan');

    %% readout weight norms
    OUT.WNORM_avrg(i)=mean(WNORM{i},'omitnan');
    OUT.WNORM_std(i)=std(WNORM{i},'omitnan')/sqrt(Nreal);

    %% readout weight norm denominator
    % sin_signal_angles=[];
    % for k=1:2:(2*Ncomp)
    %     for j=1:2:(k-1)
    %         sin_signal_angles=[sin_signal_angles, sind(OUT.ANGLES_avrg(i,k+1,j+1))];
    %     end
    % end
    OUT.SINEPROD_SIGNAL_avrg(i)=mean(SINEPROD_SIGNAL{i},'omitnan');

    %% general stats
    tmp= varfun(@nanmean, STATS{i},"InputVariables",STATS{i}.Properties.VariableNames); %#ok<NANMEAN>
    tmp.Properties.VariableNames = STATS{i}.Properties.VariableNames;
    OUT.STATS_avrg=[OUT.STATS_avrg; tmp];

    tmp= varfun(@nanstd, STATS{i},"InputVariables",STATS{i}.Properties.VariableNames); %#ok<NANSTD>
    tmp.Properties.VariableNames = STATS{i}.Properties.VariableNames;
    OUT.STATS_std=[OUT.STATS_std; tmp];

end

end

function OUT=average_ERROR(parameters,Npara,IN,RMSE,ANGLES,DIM_FT,COND,COND_TOTAL,WNORM,FLOQUET,PRINCIPAL_ANGLES)

Nreal=parameters.Nreal;     % number of repititions / network realisations
ERRoutlierflag=parameters.ERRoutlierflag;

OUT=IN;

%% correct Error outliers and NaN values
if(ERRoutlierflag)
    thr=1;
    for i=1:Npara
        RMSE{i}(RMSE{i}>thr)=1;
        RMSE{i}(isnan(RMSE{i}))=1;
    end
end

%% RMSE
tmp_RMSE=vertcat(RMSE{:});
OUT.RMSE_test=mean(tmp_RMSE,'all','omitnan')*ones(Npara,1);
OUT.RMSE_median_test=median(tmp_RMSE,'all','omitnan')*ones(Npara,1);
OUT.RMSE_SD_test=std(tmp_RMSE,0,'all','omitnan')/sqrt(Npara*Nreal) *ones(Npara,1);

%% readout weight norms
tmp=vertcat(WNORM{:});
OUT.WNORM_avrg=repmat(mean(tmp,'omitnan'),Npara,1);
OUT.WNORM_std=repmat(std(tmp,'omitnan')/sqrt(Npara*Nreal),Npara,1);

%% DIMs
% tmp=vertcat(DIM{:});
% OUT.DIM_avrg=mean(tmp,'all','omitnan')*ones(Npara,1);
% OUT.DIM_std=std(tmp,0,'all','omitnan')/sqrt(Npara*Nreal)*ones(Npara,1);

tmp=vertcat(DIM_FT{:});
OUT.DIM_FT_avrg=mean(tmp,'all','omitnan')*ones(Npara,1);
OUT.DIM_FT_std=std(tmp,0,'all','omitnan')/sqrt(Npara*Nreal)*ones(Npara,1);

%% Conditions number
tmp=vertcat(COND{:});
OUT.COND_avrg=mean(tmp,'all','omitnan')*ones(Npara,1);
OUT.COND_std=std(tmp,0,'all','omitnan')/sqrt(Npara*Nreal)*ones(Npara,1);

tmp=vertcat(COND_TOTAL{:});
OUT.COND_TOTAL_avrg=mean(tmp,'all','omitnan')*ones(Npara,1);
OUT.COND_TOTAL_std=std(tmp,0,'all','omitnan')/sqrt(Npara*Nreal)*ones(Npara,1);

%% Angles
n = size(ANGLES,1);
if(n>Npara)
    ANGLES(n)=[];
end
matSize = size(ANGLES{1},1); 
B = reshape(cell2mat(ANGLES'),matSize,[],Nreal,Npara);
B = mean(B,3,'omitnan');
B = mean(B,4,'omitnan');
for i=1:size(OUT.ANGLES_avrg,1)
    OUT.ANGLES_avrg(i,:,:)=B;
end

%% maximum Floquet multiplier
tmp=vertcat(FLOQUET{:});
OUT.FLOQUET_avrg=mean(tmp,'all','omitnan')*ones(Npara,1);
OUT.FLOQUET_std=std(tmp,0,'all','omitnan')/sqrt(Npara*Nreal)*ones(Npara,1);

%% Principal angles
tmp=vertcat(PRINCIPAL_ANGLES{:});
for i=1:size(tmp,2)
    OUT.PRINCIPAL_ANGLES_avrg(:,i)=mean(tmp(:,i),'all','omitnan')*ones(Npara,1);
    OUT.PRINCIPAL_ANGLES_std(:,i)=std(tmp(:,i),0,'all','omitnan')/sqrt(Npara*Nreal)*ones(Npara,1);
end

end

function Set_CCxC_scan_FigureProperties(ah,X,plotpara)

if(plotpara.sNflag==1)
    set(ah, 'XScale', 'log');
    xlim(ah,[X(1)/2 X(end)*2]);
    xt=logspace(log10(X(1)),0,5);
    xticks(ah,xt);
    xtickformat(ah,'%g');
    xticklabels(ah,{'0',xt(2:end)});
end
if(plotpara.xlogflag==1)
    set(ah, 'XScale', 'log');
    xticks([0.01 0.1 1 10 100]);
    xticklabels(ah,{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}'});
    xtickformat(ah,'%g');
end
if(plotpara.stabmarginflag==1)
    for i=1:numel(plotpara.xintersect)
        for k=1:numel(plotpara.xintersect{i})
            xline(plotpara.xintersect{i}(k),plotpara.xlinestyle,'Color',plotpara.xlinecol{i},'LineWidth',plotpara.LW);
        end
    end

    for i=1:numel(plotpara.xintersect)
        y_vert = [0.1 0.1 0.2 0.2]+(i-1)*0.1;
        if(isscalar(plotpara.xintersect{i}))
            x_vert = [plotpara.xintersect{i}(1) X(end) X(end) plotpara.xintersect{i}(1)];
            patch(x_vert,y_vert,plotpara.xlinecol{i},'EdgeColor','none');
        else
            for k=1:(numel(plotpara.xintersect{i})-1)
                x_vert = [plotpara.xintersect{i}(k) plotpara.xintersect{i}(k+1) plotpara.xintersect{i}(k+1) plotpara.xintersect{i}(k)];
                patch(x_vert,y_vert,plotpara.xlinecol{i},'EdgeColor','none');
            end
        end
    end
end
xlabel(plotpara.Xstr,'FontSize',plotpara.FS);
line = get(ah,'Children');
for i=1:length(line)
    set(line(i),'Linewidth', plotpara.LW);
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
