function [CTX, X]=fig_RNN_scan_release(data,varargin)

trainflag=0; logYflag=0; traceflag=0; statsflag=0; paraflag=0; medianflag=0; RMSEflag=0; dimflag=0;
Xangleflag=0; Wnormflag=0; condflag=0; floquetflag=0; floquetfracflag=0;
legendflag=0;ERRoutlierflag=1;sNflag=0;anglesflag=0;angles_signal_noiseflag=0;angles_noise_noiseflag=0;
normflag=0;normnoise_flag=0;saveflag=0;xlogflag=0;ERRdistrbflag=0;pranglesflag=0;stabmarginflag=0;
indvsineprodflag=0;fullsineprodflag=0;compactflag=0;labelsflag=0;markeroffflag=0;
filenamestr='';
Xangleidx=[2,1];
outlthr=1;
FMcutoff=1.0;
Nang=12;
maxprang=60;

%% read command line options
% varargin=varargin{1};
argidx = 1;
while argidx<=length(varargin)
    if strcmpi( varargin{argidx}, 'show traces' ) || strcmpi( varargin{argidx}, 'traces' )
        traceflag = 1;
        tidx=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'show training' ) || strcmpi( varargin{argidx}, 'training' )
        trainflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'RMSE log axis' ) || strcmpi( varargin{argidx}, 'log axis' )
        logYflag = 1;
        logYlims=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'show stats' ) || strcmpi( varargin{argidx}, 'stats' )
        statsflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show parameters' ) || strcmpi( varargin{argidx}, 'parameters' )
        paraflag = 1;
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
    elseif strcmpi( varargin{argidx}, 'show angles' ) || strcmpi( varargin{argidx}, 'angles' )
        anglesflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show angles signal noise' ) || strcmpi( varargin{argidx}, 'angles signal noise' )
        angles_signal_noiseflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show angles noise noise' ) || strcmpi( varargin{argidx}, 'angles noise noise' )
        angles_noise_noiseflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'number of angles' ) || strcmpi( varargin{argidx}, 'Nang' )
        Nang =varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'show norms' ) || strcmpi( varargin{argidx}, 'norms' )
        normflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show noise norms' ) || strcmpi( varargin{argidx}, 'noise norms' )
        normnoise_flag = 1;
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
    elseif strcmpi( varargin{argidx}, 'show ERROR distributions' ) || strcmpi( varargin{argidx}, 'ERROR distributions' ) ||...
            strcmpi( varargin{argidx}, 'ERROR distrb')
        ERRdistrbflag = 1;
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
    elseif strcmpi( varargin{argidx}, 'X axis angle') || strcmpi( varargin{argidx}, 'X angle')
        Xangleflag=1;
        Xangleidx=varargin{argidx+1};
        if( length(Xangleidx)~=2)
            error('length(Xangleidx) must be 2');
        end
        argidx=argidx+2;

    elseif strcmpi( varargin{argidx}, 'markers off' ) || strcmpi( varargin{argidx}, 'marker off' )
        markeroffflag=1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'savefig' ) || strcmpi( varargin{argidx}, 'save' )
        filenamestr=varargin{argidx+1};
        saveflag = 1;
        argidx=argidx+2;
    else
        argidx=argidx+1;
    end
end

%% set parameters
parameters=data.parameters;
learnflag=parameters.flags.learnflag;
if(learnflag==100)
    X=unique(parameters.sigmaN);
    Xstr='\sigma_N';
    savestr='sigN';
    sNflag=1;
    % replace sigmaN=0 with finite value in order to plot on log-scale
    X(1)=X(2)/10;
elseif(learnflag==102)
    X=unique(parameters.Ncx);
    Xstr='N_{ctx}';
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
elseif(learnflag==113)
    X=unique(parameters.g);
    Xstr='g';
    savestr='g';
elseif(learnflag==114)
    X=unique(parameters.eta);
    Xstr='\eta (1/s)';
elseif(learnflag==118)
    X=parameters.Phi;
    Xstr='\phi';
    savestr='Phi';
elseif(learnflag==119)
    X=parameters.Ampmod;
    Xstr='Ampmod';
elseif(learnflag==120)
    X=parameters.sigmaLS;
    Xstr='\sigma_{LS}';
    xlogflag=1;
elseif(learnflag==128)
    X=acosd(parameters.rho_Wfb12);
    Xstr='angle Wfb12';
    savestr='rho12';
elseif(learnflag==129)
    X=unique(parameters.sigAmp_intf);
    Xstr='Intf. signal Amp';
    savestr='sigAmpintf';
    xlogflag=1;
end

signalflag=parameters.flags.signalflag;

Npara=parameters.Npara;
Nreal=parameters.Nreal;
Ncomp=max(parameters.Ncomp);

if(parameters.Ncomp<2)
    anglesflag=0;
end
if(parameters.flags.floquetflag==0)
    stabmarginflag=0;
    floquetflag=0;
    floquetfracflag=0;
end

parameters.ERRoutlierflag=ERRoutlierflag;
parameters.FMcutoff=FMcutoff;
interferenceflag=parameters.flags.interferenceflag;

if(length(data.simcount)<Nreal*Npara)
    Npara=size(data.time_test,1)-1;
    X=X(1:Npara);
end

%% filename
if(saveflag==1)
    if(interferenceflag>0)
        intfname=['_Nintf' num2str(data.parameters.Ncomp_intf) ];
    else
        intfname='';
    end
    filename=[savestr '_Ncomp' num2str(Ncomp) intfname '_' data.parameters.learningrule '_' filenamestr];
end

%% extract learning performance, measures, stats etc.
CTX=extract_stats_v3 (parameters,Npara,data,'ctx',outlthr);

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
CTX.labels_str=labels_str;
if(learnflag==128)
    X=atanh(X);
    X(isinf(X))=max(X(~isinf(X)))*2;
end

%% uses angles as x-axis
if(Xangleflag==1)
    if(Xangleidx(2)>Xangleidx(1))
        Xangleidx=flip(Xangleidx);
    end

    % X=atanh(X);
    % X(isinf(X))=max(X(~isinf(X)))*2;

    % Xtest=CEREB.ANGLES_avrg(:,(2*Ncomp+1)+1,1+1); % v+ z+
    % Xstr='angle (v_+, z_+)';

    % X=CEREB.ANGLES_avrg(:,(2*Ncomp+1)+1,1+2); % v- z+
    % Xstr='angle (v_-, z_+)';

    % X=CEREB.ANGLES_avrg(:,(2*Ncomp+1)+1,1+3); % v2+ z+
    % Xstr='angle (v^2_+, z_+)';

    X=CTX.ANGLES_avrg(:,Xangleidx(1)+1,Xangleidx(2)+1);
    Xstr=['Angle ', labels_str{Xangleidx(1)+1,Xangleidx(2)+1}];

    xlogflag=0;
    sNflag=0;
end

%% estimate stability margins using error/EV/Floquet
cmp = get(groot,'DefaultAxesColorOrder');
ctxcolour=cmp(1,:);
suppcolour=cmp(5,:);
xlinecol='b';
xlinestyle='--';
xintersect=[];
if(stabmarginflag==1)
    if(xlinecritflag==1)
        xintersect = findIntersect(X,CTX.RMSE_test,xlinecutoff*ones(numel(X),1),xlogflag,1);
    elseif(xlinecritflag==2)
        xintersect = findIntersect(X,CTX.FLOQUET_avrg,xlinecutoff*ones(numel(X),1),xlogflag,1);
    elseif(xlinecritflag==4)
        xintersect = findIntersect(X,CTX.FLOQUET_median,xlinecutoff*ones(numel(X),1),xlogflag,1);
        xlinecol=ctxcolour*0.7;
    elseif(xlinecritflag==5)
        xintersect = findIntersect(X,CTX.FLOQUET_frac,xlinecutoff*ones(numel(X),1),xlogflag,0);
    end
end

%% figures
FS=12;
LW=1.5;
MS=6;
set(0,'defaultAxesFontName', 'Arial');

if(compactflag==1)
    dimflag=0;
    condflag=0;
end

Nflags=trainflag+dimflag+statsflag+paraflag+RMSEflag+condflag...
    +anglesflag+angles_signal_noiseflag+angles_noise_noiseflag+indvsineprodflag+fullsineprodflag...
    +normflag+normnoise_flag+Wnormflag...
    +floquetflag+floquetfracflag+pranglesflag+compactflag;
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
    MS=4;
end

if(markeroffflag==1)
    MS=1;
end

% plot parameters
plotpara.Ncomp=Ncomp;
plotpara.xlogflag=xlogflag;
plotpara.sNflag=sNflag;
plotpara.Xstr=Xstr;
plotpara.FS=FS;
plotpara.LW=LW;
plotpara.MS=MS;
plotpara.ctxcolour=ctxcolour;
plotpara.learnflag=learnflag;
plotpara.stabmarginflag=stabmarginflag;
plotpara.xintersect=xintersect;
plotpara.labelsflag=labelsflag;
plotpara.xlinecol=xlinecol;
plotpara.xlinestyle=xlinestyle;

if(ERRdistrbflag)
    Nflags=0;
end

ii=0;
if(Nflags>0)
    f1=figure;
    set(f1,'papertype','A4','Units', 'centimeters', 'position', [10,40,Fwidth,sy*nr]);
    tiledlayout(nr,ceil(Nflags/nr),'TileSpacing',tilestr,'Padding',paddstr);

    %% panels
    if(RMSEflag==1)
        nexttile;
        hold on;
        plot(X,CTX.RMSE_test,'o-','Markerfacecolor',ctxcolour,'color',ctxcolour);
        if(medianflag)
            plot(X,CTX.RMSE_median_test,':','Markerfacecolor',ctxcolour,'color',ctxcolour);
            if(legendflag)
                legend('S-RNN','S-RNN median','Location','Best');
            end
        else
            if(legendflag)
                legend('ctx','Location','Best');
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
        errorbar(X,CTX.RMSE_train,CTX.RMSE_SD_train,'o-','Markerfacecolor',ctxcolour,'color',ctxcolour);
        ylabel('Training error','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        if(logYflag)
            set(gca, 'YScale', 'log');
        end
    end

    if(statsflag)
        nexttile;
        hold on;
        errorbar(X,CTX.STATS_avrg.hcx_std,CTX.STATS_std.hcx_std,'o--','Markerfacecolor',ctxcolour,'color',ctxcolour);
        errorbar(X,CTX.STATS_avrg.nucx_std,CTX.STATS_std.nucx_std,'o-','Markerfacecolor',ctxcolour,'color',ctxcolour);
        if(any("nucx_pop_std" == string(CTX.STATS_avrg.Properties.VariableNames)))
            errorbar(X,CTX.STATS_avrg.nucx_pop_std,CTX.STATS_std.nucx_pop_std,'o-','Markerfacecolor',suppcolour,'color',suppcolour);
        end
        legend('\sigma_{h}','\sigma_{\nu}','Location','Best');
        ylabel('Stats','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
    end

    if(anglesflag)
        plot_angles(X,CTX.ANGLES_avrg,labels_str,plotpara,ctxcolour,'Signal',Nang,'mode','signal');
    end

    if(angles_signal_noiseflag)
        plot_angles(X,CTX.ANGLES_avrg,labels_str,plotpara,cmp(4,:),'Overtones-to-signal',Nang,'mode','signal noise');
    end

    if(angles_noise_noiseflag)
        plot_angles(X,CTX.ANGLES_avrg,labels_str,plotpara,cmp(7,:),'Overtones',Nang,'mode','noise noise');
    end

    if(normflag)
        ii=ii+1;
        axnorm(ii)=nexttile;
        hold on;
        idx=1:2:(2*Ncomp);
        % plot(X,CTX.NORMS_avrg(:,idx+1),'o-','Color',[ctxcolour,0.3]);
        plot(X,mean(CTX.NORMS_avrg(:,idx+1),2),'o-','Color',ctxcolour,'Markerfacecolor',ctxcolour);
        idx=1:2:(2*Ncomp);
        for k=idx
            [normmax(k), maxidx]=max(CTX.NORMS_avrg(:,k+1)); %#ok<AGROW>
            if(labelsflag==1)
                text(X(maxidx),normmax(k),labels_str(k+1,1),'Color',ctxcolour);
            end
        end
        normmax=max(normmax);
        if(normmax<1)
            normmax=1;
        end
        set(gca, 'YScale', 'log');
        ylim([-Inf normmax*1.05]);
        ylabel('Norm (a.u.)','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        title('Signal'); ax = gca; ax.TitleFontSizeMultiplier = 0.8;
    	if(normnoise_flag)
            linkaxes(axnorm,'y');
    	end
        ylim(axnorm(end),[0 normmax*1.05])
    end

    if(normnoise_flag)
        ii=ii+1;
        axnorm(ii)=nexttile;
        hold on;
        maxidx=(2*Ncomp)+Nang;
        if(maxidx>size(CTX.NORMS_avrg,2))
            maxidx=size(CTX.NORMS_avrg,2)-1;
        end
        idx=[2:2:(2*Ncomp) (2*Ncomp+1):maxidx];
        for i=1:length(idx)
            if(idx(i)<=(2*Ncomp))
                plot(X,CTX.NORMS_avrg(:,idx(i)+1),'--','Color',cmp(4,:));
            else
                plot(X,CTX.NORMS_avrg(:,idx(i)+1),'-','Color',cmp(4,:));
            end
        end
        for k=idx
            [normmax(k), maxidx]=max(CTX.NORMS_avrg(:,k+1));
            if(labelsflag==1)
                text(X(maxidx),normmax(k),labels_str(k+1,1),'Color',cmp(4,:));
            end
        end

        set(gca, 'YScale', 'log');
        % ylim([-Inf normmax*1.05]);
        ylabel('Norm (a.u.)','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        title('Overtones'); ax = gca; ax.TitleFontSizeMultiplier = 0.8;
    	if(normflag)
            linkaxes(axnorm,'y');
    	end
        ylim(axnorm(end),[0 max(normmax)*1.05])
    end

    if(condflag)
        nexttile;
        hold on;
        plot(X,CTX.COND_avrg,'-','Color',cmp(1,:));
        plot(X,CTX.COND_TOTAL_avrg,'--','Color',cmp(1,:)*0.7);

        ylim([1 inf]);
        set(gca, 'YScale', 'log');
        ylabel('CN','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        % yticks([1 1e5 1e10]);
        yticks([1 1e2 1e4 1e6]);
        if(legendflag)
            legend('signal','total','Location','Best');
        end
    end

    if(dimflag)
        nexttile;
        hold on;
        plot(X,CTX.DIM_avrg,'--','Color',ctxcolour);
        plot(X,CTX.DIM_FT_avrg,'-','Color',ctxcolour);
        ylabel('Dimensionality','FontSize',FS);
        ylim([1 inf]);

        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        if(legendflag)
            legend('ctx','signal ctx','Location','Best');
        end
    end

    if(compactflag)
        nexttile;
        yyaxis left
        hold on;
        plot(X,CTX.COND_avrg,'o-','Color',cmp(1,:));
        % plot(X,CTX.COND_TOTAL_avrg,'--','Color','k');

        ylim([1 inf]);
        set(gca, 'YScale', 'log');
        ylabel('CN','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);

        yyaxis right
        plot(X,CTX.DIM_FT_avrg,'o-','Color','k');
        ylim([1 inf]);
        ylabel('Dimensionality','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);

        ax = gca();
        ax.YAxis(1).Color = cmp(1,:);
        ax.YAxis(2).Color = [0 0 0];
    end

    if(pranglesflag)
        nexttile;
        hold on;
        plot(X,CTX.PRINCIPAL_ANGLES_avrg,'o-','Color',ctxcolour,'MarkerFaceColor',ctxcolour);

        ylim([0 maxprang]);
        ylabel('Principal angle (deg.)','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
    end

    if(indvsineprodflag)
        nexttile;
        hold on;
        % task-relevant separation
        plot(X,CTX.SINEPROD_SIGNAL_avrg,'o-','Color',ctxcolour,'MarkerFaceColor',ctxcolour);
        % task-relevant vs task-irrelevant separation
        plot(X,prod(sind(CTX.PRINCIPAL_ANGLES_avrg),2),'o--','Color',ctxcolour,'MarkerFaceColor',ctxcolour);

        ylim([-Inf 1]);
        ylabel('Sine product','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        set(gca, 'YScale', 'log');
    end

    if(fullsineprodflag)
        nexttile;
        hold on;

        plot(X,CTX.SINEPROD_SIGNAL_avrg.*prod(sind(CTX.PRINCIPAL_ANGLES_avrg),2),'o-','Color',ctxcolour,'MarkerFaceColor',ctxcolour);

        ylim([-Inf 1]);
        ylabel('Global sine product','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        set(gca, 'YScale', 'log');
    end

    if(Wnormflag)
        nexttile;
        hold on;
        plot(X,CTX.WNORM_avrg,'o-','Color',ctxcolour,'Markerfacecolor',ctxcolour);
        ylabel('||W_{out}||','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        set(gca, 'YScale', 'log');
    end

    if(floquetflag)
        nexttile;
        hold on;
        plot(X,CTX.FLOQUET_avrg,'o-','Color',ctxcolour,'Markerfacecolor',ctxcolour);
        if(medianflag==1)
            plot(X,CTX.FLOQUET_median,':','Color',ctxcolour*0.7);
        end
        yline(1,'k--');
        ylim([0.01 1e5]);
        ylabel('Max. Abs. FM','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        set(gca, 'YScale', 'log');
        yticks(10.^(log10(1):2:log10(10^6)));
    end

    if(floquetfracflag)
        nexttile;
        hold on;
        plot(X,CTX.FLOQUET_frac,'o-','Color',ctxcolour,'Markerfacecolor',ctxcolour);
        ylim([0 1]);
        ylabel(['Fraction FM > ' num2str(FMcutoff)],'FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
    end

    if(paraflag)
        LRstring=parameters.learningrule;
        SIGstring=parameters.signaltype;
        RNNactstring=parameters.activation_function;
        if(parameters.flags.LRflag>=2)
            SIGMALSstring=[' with \sigma_{LS} = ', num2str(parameters.sigmaLS(end))];
        else
            SIGMALSstring='';
        end
        if(isfield(parameters.flags,'loopflag_test'))
            if(parameters.flags.loopflag_test==0)
                LOOPstring=['OPEN loop testphase, with \sigma_{Nopen} = ' num2str(parameters.openloop_test_noise)];
            elseif(parameters.flags.loopflag_test==1)
                LOOPstring='CLOSED loop testphase';
            elseif(parameters.flags.loopflag_test==-1)
                LOOPstring=['OPEN loop testphase, only noise with \sigma_{Nopen} = ' num2str(parameters.openloop_test_noise)];
            end
        else
            LOOPstring='';
        end

        taumodestring='';
    	if(parameters.flags.taudistrbflag==0)
    		taumodestring= [parameters.tau_mode,' with mean = ', num2str(round(parameters.tau_avrg(end)*1000,2))];
    	elseif(parameters.flags.taudistrbflag==1)
    		taumodestring= [parameters.tau_mode,' with mean = ', num2str(round(parameters.tau_avrg(end)*1000,2)) ', SD = ', num2str(round(parameters.tau_std(end)*1000,2)) ];
    	end

        parastring1=['N_{ctx} = ', num2str(parameters.Ncx(end)), ', \sigma_N = ', num2str(parameters.sigmaN(end))];
        parastring2=['w = ', num2str(parameters.lambda(end)), ', \sigma_{Amp} = ', num2str(parameters.sigAmp(end))];

        dtstring=sprintf('%.1e',parameters.dt(1).*parameters.lambda(1));

        nexttile;
        xpos=-0.1;
        ypos=1.1;
        ystep=0.09;
        text(xpos,ypos-1*ystep,['Learning rule: ', LRstring, SIGMALSstring, ', N_{real} = ', num2str(parameters.Nreal), ', dt = ', dtstring], 'FontWeight', 'Bold');
        text(xpos,ypos-2*ystep,LOOPstring, 'FontWeight', 'Bold');
        text(xpos,ypos-3*ystep,SIGstring, 'FontWeight', 'Bold');
        text(xpos,ypos-4*ystep,['RNN act: ', RNNactstring], 'FontWeight', 'Bold');
        text(xpos,ypos-5*ystep,['g = ', num2str(parameters.g(end))], 'FontWeight', 'Bold');
        text(xpos,ypos-7*ystep,taumodestring, 'FontWeight', 'Bold');

        text(xpos,ypos-9*ystep,parastring1, 'FontWeight', 'Bold');
        text(xpos,ypos-10*ystep,parastring2, 'FontWeight', 'Bold');

        if(interferenceflag>0)
            if(interferenceflag==1)
                intflocus= 'ctx intf: ';
            elseif(interferenceflag==2)
                intflocus= 'GC intf: ';
            end
            interferencestring=[intflocus, 'N^{intf}_{comp} = ', num2str(parameters.Ncomp_intf),...
                ', \rho_{Wfb} = ', num2str(parameters.rho_Wfb(end)) ,', \sigma_{Amp,intf} = ', num2str(parameters.sigAmp_intf(end))];
            text(xpos,ypos-11*ystep,interferencestring, 'FontWeight', 'Bold');
        end
        axis off ;

    end

    if(saveflag==1)
        filenamefull=['temp_figures/scan_' filename];
        saveas(f1,filenamefull,'svg');
    end

end

%% plot traces
if(traceflag)
    tidx=unique(tidx);
    tidx=tidx(tidx<=Npara);
    if(~isempty(tidx))
        FS=12;
        nr=5;
        nc=length(tidx);
        sx=6*nc;
        sy=5*nr;
        time=data.time_test;
        Z=data.Z_test;
        target=data.target_test;

        ftraces=figure;
        set(ftraces,'papertype','A4','Units', 'centimeters', 'position', [10,40,sx,sy]);
        tiledlayout(nr,nc,'TileSpacing','compact','Padding','compact');
        k=1;
        for i=tidx
            xmin=min(data.time_test{i});
            xmax=max(data.time_test{i})*1.05;

            % average output
            if(signalflag==2)
                % when targets are different for each realisation choose
                % random example traces
                iex=randi([1, Nreal]);
                Z_plot=Z{i}(:,iex);
            else
                Z_plot=mean(Z{i},2);
            end
            nexttile(k+4*nc);
            hold on;
            plot(time{i},Z_plot,'Color',ctxcolour*0.5);
            xlim([xmin xmax]);
            %         ylim([-ylimits ylimits]);
            %         yticks([-1,0,1]);
            set(gca, 'TickDir', 'out');box off;
            set(gca,'FontSize',FS); set(gca, 'Color', 'None');
            if(k==1)
                ylabel('Avrg trace','FontSize',FS);
            end

            iex=randperm(Nreal);
            % sample output 1
            nexttile(k+nc);
            plot(time{i},Z{i}(:,iex(1)),'Color',ctxcolour);
            xlim([xmin xmax]);
            %         ylim([-ylimits ylimits]);
            %         yticks([-1,0,1]);
            set(gca, 'TickDir', 'out');box off;
            set(gca,'FontSize',FS); set(gca, 'Color', 'None');
            if(k==1)
                ylabel('Sample trace','FontSize',FS);
            end

            % sample output 2
            nexttile(k+2*nc);
            plot(time{i},Z{i}(:,iex(2)),'Color',ctxcolour);
            xlim([xmin xmax]);
            %         ylim([-ylimits ylimits]);
            %         yticks([-1,0,1]);
            set(gca, 'TickDir', 'out');box off;
            set(gca,'FontSize',FS); set(gca, 'Color', 'None');
            if(k==1)
                ylabel('Sample trace','FontSize',FS);
            end

            % sample output 3
            nexttile(k+3*nc);
            plot(time{i},Z{i}(:,iex(3)),'Color',ctxcolour);
            xlim([xmin xmax]);
            %         ylim([-ylimits ylimits]);
            %         yticks([-1,0,1]);
            set(gca, 'TickDir', 'out');box off;
            set(gca,'FontSize',FS); set(gca, 'Color', 'None');
            if(k==1)
                ylabel('Sample trace','FontSize',FS);
            end

            % target signal
            nexttile(k);
            plot(time{i},target{i}(:,1),'Color','k');
            xlim([xmin xmax]);
            %         ylim([-ylimits ylimits]);
            yticks([-1,0,1]);
            set(gca, 'TickDir', 'out');box off;
            set(gca,'FontSize',FS); set(gca, 'Color', 'None');
            if(k==1)
                ylabel('Target signal','FontSize',FS);
            end
            if(sNflag && i==1)
                title(['Npara = ' num2str(i) ', ' Xstr ' = ' num2str(0)],'FontSize',11);
            else
                title(['Npara = ' num2str(i) ', ' Xstr ' = ' num2str(X(i))],'FontSize',11);
            end

            k=k+1;
        end
    end
end

%% Error distributions figure
if(ERRdistrbflag)

    ncerr=1+floquetflag;
    FS=12;
    sx=12;
    sy=10;
    % MS=3;
    Xdistrb=X+randn(1,Nreal)*0.01;
    if(saveflag==1)
        % MS=4;
        sx=4.5;
        sy=sx;
        % Fwidth=sx*Nflags/nr;
        tilestr='tight';
        paddstr='tight';
        FS=10;
    else
        % Fwidth=sx+sx*Nflags/nr;
        tilestr='compact';
        paddstr='compact';
    end

    ferrdistrb=figure;
    set(ferrdistrb,'papertype','A4','Units', 'centimeters', 'position', [10,40,sx*ncerr,sy]);
    tiledlayout(1,ncerr,'TileSpacing',tilestr,'Padding',paddstr);

    nexttile;
    hold on;
    plot(Xdistrb,CTX.RMSE_distrb,'o','MarkerEdgeColor','None','Markerfacecolor',ctxcolour,'color',ctxcolour);
    plot(X,CTX.RMSE_test,'-','color',ctxcolour*0.7);
    if(medianflag)
        plot(X,CTX.RMSE_median_test,'-.','color',ctxcolour*0.7);
    end

    set(gca, 'YScale', 'log');
    ylim([logYlims(1) logYlims(2)]);
    ylabel(['Test error (' 'RMSE' ')'],'FontSize',FS);
    %     title(LRstring);
    Set_CCxC_scan_FigureProperties(gca,X,plotpara);

    if(floquetflag)
        nexttile;
        hold on;
        plot(Xdistrb,CTX.FLOQUET_distrb,'o','MarkerEdgeColor','None','Markerfacecolor',ctxcolour,'color',ctxcolour);
        plot(X,CTX.FLOQUET_avrg,'-','color',ctxcolour*0.7);
        if(medianflag)
            plot(X,CTX.FLOQUET_median,'-.','color',ctxcolour*0.7);
        end

        set(gca, 'YScale', 'log');
        yline(1,'k--');
        ylim([0.001 1e6]);
        ylabel('Max. Abs. FM','FontSize',FS);
        Set_CCxC_scan_FigureProperties(gca,X,plotpara);
        yticks(10.^(log10(0.001):2:log10(10^6)));
    end

    if(saveflag==1)
        filenamefull=['temp_figures/ERRdistrb_' filename];
        saveas(ferrdistrb,filenamefull,'svg');
    end

end

end

%% functions
function OUT=extract_stats_v3 (parameters,Npara,data,typestring,outlthr)
Nreal=parameters.Nreal;     % number of repititions / network realisations

ERRoutlierflag=parameters.ERRoutlierflag;

% get relevant data
RMSE_test=Get_structs(data,typestring,'RMSE_test',parameters,1,Npara);
RMSE_train=Get_structs(data,typestring,'RMSE_train',parameters,1,Npara);
ANGLES=Get_structs(data,typestring,'Angles',parameters,100,Npara);
ANGLES_RAW=ANGLES;
NORMS=Get_structs(data,typestring,'Norms',parameters,20,Npara);
DIM=Get_structs(data,typestring,'DIM',parameters,1,Npara);
DIM_FT=Get_structs(data,typestring,'DIM_FT',parameters,1,Npara);
COND=Get_structs(data,typestring,'COND',parameters,1,Npara);
COND_TOTAL=Get_structs(data,typestring,'COND_total',parameters,1,Npara);
FLOQUET=Get_structs(data,typestring,'Max_Floquet_multiplier',parameters,1,Npara);
PRINCIPAL_ANGLES=Get_structs(data,typestring,'Subspaceangles_cos',parameters,1,Npara);
SINEPROD_SIGNAL=Get_structs(data,typestring,'Sineprod_signal',parameters,1,Npara);

WNORM=Get_structs(data,typestring,'Wnorm',parameters,1,Npara);

STATS=data.('stats');

OUT.RMSE_distrb=zeros(Npara,Nreal);

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

    %% full RMSE distributions
    OUT.RMSE_distrb(i,:)=RMSE_test{i}(:,1);

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

    %% Principal angles
    OUT.PRINCIPAL_ANGLES_avrg(i,:)=mean(PRINCIPAL_ANGLES{i},'omitnan');
    OUT.PRINCIPAL_ANGLES_std(i,:)=std(PRINCIPAL_ANGLES{i},'omitnan')/sqrt(Nreal);

    %% spanning vector angles
    OUT.ANGLES_avrg(i,:,:)=mean(ANGLES{i},3,'omitnan');
    OUT.ANGLES_std(i,:,:)=std(ANGLES{i},0,3,'omitnan')/sqrt(Nreal);

    %% spanning vector norms
    OUT.NORMS_avrg(i,:)=mean(NORMS{i},'omitnan');
    OUT.NORMS_std(i,:)=std(NORMS{i},'omitnan')/sqrt(Nreal);

    %% readout weight norms
    OUT.WNORM_avrg(i)=mean(WNORM{i},'omitnan');
    OUT.WNORM_std(i)=std(WNORM{i},'omitnan')/sqrt(Nreal);

    %% signal sine product
    % sin_signal_angles=[];
    % for k=1:2:(2*Ncomp)
    %     for j=1:2:(k-1)
    %         sin_signal_angles=[sin_signal_angles, sind(OUT.ANGLES_avrg(i,k+1,j+1))]; %#ok<AGROW>
    %     end
    % end
    % OUT.SINEPROD_SIGNAL_avrg2(i)=prod(sin_signal_angles);
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

function Set_CCxC_scan_FigureProperties(ah,X,plotpara)
% ah = axis handle

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
    for k=1:numel(plotpara.xintersect)
        xline(plotpara.xintersect(k),plotpara.xlinestyle,'Color',plotpara.xlinecol,'LineWidth',plotpara.LW);
    end
end
xlabel(plotpara.Xstr,'FontSize',plotpara.FS);
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
        outdata{i,1}=nan([parameters.Nreal,Narray2]);  %#ok<AGROW>
    end
end
end

function plot_angles(X,angles,labels_str,plotpara,colour,titlestr,Nang,varargin)

anglesmodeflag=0;
argidx = 1;
while argidx<=length(varargin)
    if strcmpi( varargin{argidx}, 'mode' )
        anglemodestr=varargin{argidx+1};
        if strcmpi(anglemodestr,'signal') % default option
            anglesmodeflag=0;
        elseif strcmpi(anglemodestr,'signal noise') || strcmpi( anglemodestr, 'noise signal' )
            anglesmodeflag=1;
        elseif strcmpi(anglemodestr,'noise noise' )
            anglesmodeflag=2;
        else
            error('Wrong anglemode string.');
        end
        argidx=argidx+2;
    else
        argidx=argidx+1;
    end
end

if(anglesmodeflag==0)
    kidx=1:2:(2*plotpara.Ncomp);
    jidx=@(k)(1:2:(k-1));
elseif(anglesmodeflag==1)
    kidx=[2:2:(2*plotpara.Ncomp) (2*plotpara.Ncomp+1):(2*plotpara.Ncomp+Nang)];
    jidx=@(k) 1:(2*plotpara.Ncomp);
    % jidx=@(k)jidx_signal_noise(k,plotpara.Ncomp,Nang);
elseif(anglesmodeflag==2)
    kidx=(2*plotpara.Ncomp+1):(2*plotpara.Ncomp+1+Nang);
    jidx=@(k)(2*plotpara.Ncomp+1):(k-1);
end
kidx(kidx>=size(angles,3))=[];
nexttile;
hold on;
for k=kidx
    for j=jidx(k)
        if(all(angles(:,k+1,j+1)~=0))
            plot(X,angles(:,k+1,j+1),'o-','Color',colour,'Markerfacecolor',colour,'Linewidth',plotpara.LW);
            if(plotpara.labelsflag==1)
                if(plotpara.learnflag==128 || plotpara.learnflag==132)
                    [angmax, maxidx]=min(angles(:,k+1,j+1));
                else
                    [angmax, maxidx]=max(angles(:,k+1,j+1));
                end
                text(X(maxidx),angmax,labels_str(k+1,j+1),'Color',colour);
            end
        end
    end
end

ylim([0 100]);
ylabel('Angle (deg.)','FontSize',plotpara.FS);
Set_CCxC_scan_FigureProperties(gca,X,plotpara);
title(titlestr); ax = gca; ax.TitleFontSizeMultiplier = 0.8;
end