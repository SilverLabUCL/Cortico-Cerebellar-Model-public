function [Dcx,Dgc]=fig_activity_release(data,varargin)

argidx = 1;
outputflag=0;
activityflag=0;
modeflag=0;
NPCA=4;
% varargin=varargin{1};
while argidx<=length(varargin)
    if strcmpi( varargin{argidx}, 'show activity' )
        activityflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show output' )
        outputflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'show modes' ) || strcmpi( varargin{argidx}, 'modes' )
        modeflag = 1;
        NPC=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'NPCA' )
        NPCA=varargin{argidx+1};
        argidx=argidx+2;
    else
        argidx=argidx+1;
    end
end

CORTEXONLYflag=data.parameters.flags.cortexonlyflag;
Nrnn=data.parameters.Ncx;
Nffl=data.parameters.Ngc;

RNN_testphase=data.Cx_test;
RNN_FTphase=data.Cx_FT;
if(CORTEXONLYflag==0)
    FFL_testphase=data.GC_test;
    FFL_FTphase=data.GC_FT;
else
    FFL_testphase=nan(size(RNN_testphase));
    FFL_FTphase=nan(size(RNN_FTphase));
end
time_learn=data.time_learn{1};
time_test=data.time_test{1};
time_FTphase=data.time_FT{1};
dt=time_test(2)-time_test(1);

learnflag=data.parameters.flags.learnflag;
if(learnflag==1)
    OUT_test=data.RNNoutput_test;
    target_test=data.target_test;
    RMSE=data.RMSE_test_ctx;
    tltstr=['RNN only', data.parameters.learningrule];
elseif(learnflag==2)
    OUT_test=data.FFLoutput_test;
    target_test=data.target_test;
    RMSE=data.RMSE_test_cereb;
    tltstr=['RNN+FFL', data.parameters.learningrule];
else
    tltstr='no learning';
end

% PCAs
opt = statset('pca'); opt.Display = 'off';
[~,scoreCx,latentCx,~,~] = pca(RNN_FTphase,'Options',opt);
Dcx=sum(latentCx)^2/sum(latentCx.^2);
if(CORTEXONLYflag==0)
    [~,scoreGC,latentGC,~,~] = pca(FFL_FTphase,'Options',opt);
    Dgc=sum(latentGC)^2/sum(latentGC.^2);
else
    scoreGC=nan(size(scoreCx));
    Dgc=nan;
end

sx=30;
sy=10;
FS=10;
LW=2;
set(0,'defaultAxesFontName', 'Arial');
cmp = get(groot,'DefaultAxesColorOrder');
if(learnflag==1)
    outcolor=cmp(1,:);
else
    outcolor=cmp(2,:);
end
%% test output signal
if(outputflag && learnflag>0)
    f1=figure;
    set(f1,'papertype','A4','Units', 'centimeters', 'position', [10,40,sx,sy]);

    hold on;
    plot(time_test,target_test{1},'k','LineWidth',LW);
    plot(time_test,OUT_test{1},'--','color',outcolor,'LineWidth',LW);
    Set_CCxC_plot_FigureProperties(gcf,gca,'Time (s)','output rate',FS,'axis','normal','ylim',[-inf inf]);
    TE=round(1000*RMSE(1))/1000;
    title([tltstr ', test error =', num2str(TE)]);
end

%% plot figure with principal component and Fourier modes
if(modeflag)
    nc=3;
    nr=3;
    sx=nc*8;
    sy=5*nr;
    FSp=10;
    
    clrfac_pca=linspace(1,0.6,NPCA);
    clrfac=linspace(1,0.6,NPC);

    freq=data.frequencies;
    % freqmax=max(freq(:));
    freqmin=min(freq(:));

    Tperiod=1.0/data.parameters.lambda;
    timeplot=0:dt:(2*Tperiod);       % 2D plots: plot over 2 periods
    Tidx_basis=1:length(timeplot);

    fmodes=figure;
    set(fmodes,'papertype','A4','Units', 'centimeters', 'position', [40,40,sx,sy]);
    tiledlayout(nr,nc,'TileSpacing','tight','Padding','compact');

    % activity
    nsamp=7; % number of sample traces
    nexttile;
    hold on;
    if(CORTEXONLYflag==0)
        nsampidx=randi([1 Nffl],nsamp,1); % number of sample traces
        plot(timeplot,FFL_FTphase(Tidx_basis,nsampidx),'-','LineWidth',LW,'Color',cmp(2,:));
    else
        nsampidx=randi([1 Nrnn],nsamp,1); % number of sample traces
        plot(timeplot,RNN_FTphase(Tidx_basis,nsampidx),'-','LineWidth',LW,'Color',cmp(1,:));
    end
    yticks(0);
    xticks([0 0.1 0.2]);
    Set_CCxC_plot_FigureProperties(gcf,gca,'Time (s)','Rate (a.u.)',FSp);

    ymax=max(abs(scoreCx(:,1)));
    % RNN principal components
    nexttile(4,[2 1]);
    hold on;
    for i=1:NPCA
        plot(timeplot,scoreCx(Tidx_basis,i)/ymax,'color',cmp(1,:)*clrfac_pca(i),'LineWidth',LW);
    end
    Xstr='Time (s)';
    title('PC RNN');
    yticks([-1 0 1]);
    xticks([0 0.1 0.2]);
    Set_CCxC_plot_FigureProperties(gcf,gca,Xstr,'PC',FSp,'ylim',[-1 1]);

    ymax=max(abs(scoreGC(:,1)));
    % FFL principal components
    nexttile(5,[2 1]);
    hold on;
    for i=1:NPCA
        plot(timeplot,scoreGC(Tidx_basis,i)/ymax,'color',cmp(2,:)*clrfac_pca(i),'LineWidth',LW);
    end
    Xstr='Time (s)';
    title('PC FFL');
    yticks([-1 0 1]);
    xticks([0 0.1 0.2]);
    Set_CCxC_plot_FigureProperties(gcf,gca,Xstr,'PC',FSp,'ylim',[-1 1]);

    % FC task-relevant spanning vector components
    nexttile(6);
    hold on;
    for i=1:NPC
        % norms of spanning vectors
        Snorm=vecnorm(data.stotal)*sqrt(Nrnn);
        maxnorm=max(Snorm(2:end));
        if(mod(i,2)==1)
            w=freqmin*floor((i +1)/2);
            signal=Snorm(i+1)*cos(2*pi*w*timeplot)/maxnorm;
            if(i<=2*data.parameters.Ncomp)
                colour=cmp(4,:);
                plot(timeplot,signal,'color',colour*clrfac(i),'LineWidth',LW);
            end
        end
    end
    Xstr=[];
    xticks([]);
    yticks([-1 0 1]);
    Set_CCxC_plot_FigureProperties(gcf,gca,Xstr,'FC',FSp,'ylim',[-1 1]);
    title('task-relevant frequency components');

    % FC noise spanning vector components
    nexttile(9);
    hold on;
    for i=1:NPC
        % norms of spanning vectors
        Snorm=vecnorm(data.stotal)*sqrt(Nrnn);
        maxnorm=max(Snorm(2:end));
        if(mod(i,2)==1)% all cosine components
            w=freqmin*floor((i +1)/2);
            signal=Snorm(i+1)*cos(2*pi*w*timeplot)/maxnorm;
            if(i<=2*data.parameters.Ncomp)
            else
                % task-irrelevant cosine components
                colour=cmp(5,:);
                plot(timeplot,signal,'-','color',colour*clrfac(i),'LineWidth',LW);
            end
        else % task-irrelevant sine components
            signal=Snorm(i+1)*sin(2*pi*w*timeplot)/maxnorm;
            colour=cmp(5,:);
            % if(i<=2*data.parameters.Ncomp) % sine target frequency components
                plot(timeplot,signal,'-','color',colour*clrfac(i),'LineWidth',LW);
            % else % sine overtone components
            %     plot(timeplot,signal,'-','color',colour*clrfac(i),'LineWidth',LW);
            % end
        end
    end
    yticks([-1 0 1]);
    xticks([0 0.1 0.2]);
    Set_CCxC_plot_FigureProperties(gcf,gca,'Time (s)','FC',FSp,'ylim',[-1 1]);
    title('task-irrelevant frequency components');

end
%% plot more detailed activity traces
if(activityflag)
    sx=30;
    sy=20;
    nsamp=7; % number of sample traces
    nsampidx=randi([1 Nrnn],nsamp,1); % number of sample traces
    %     nsampidx=1:nsamp;

    f7=figure;
    set(f7,'papertype','A4','Units', 'centimeters', 'position', [40,40,sx,sy]);

    subplot(3,3,1);
    plot(time_learn,data.target_learn{1});
    Set_CCxC_plot_FigureProperties(gcf,gca,'','target learn',FS);

    subplot(3,3,2);
    plot(time_test,data.target_test{1});
    Set_CCxC_plot_FigureProperties(gcf,gca,'','target test',FS);

    subplot(3,3,3);
    plot(time_FTphase,data.target_FT{1});
    Set_CCxC_plot_FigureProperties(gcf,gca,'','target FT',FS);

    subplot(3,3,5);
    plot(time_test,RNN_testphase(:,nsampidx));
    Set_CCxC_plot_FigureProperties(gcf,gca,'','Cx test',FS);

    subplot(3,3,6);
    plot(time_FTphase,RNN_FTphase(:,nsampidx));
    Set_CCxC_plot_FigureProperties(gcf,gca,'','Cx FT',FS);

    if(CORTEXONLYflag==0)

        subplot(3,3,8);
        plot(time_test,FFL_testphase(:,nsampidx));
        Set_CCxC_plot_FigureProperties(gcf,gca,'Time (s)','GC test',FS);

        subplot(3,3,9);
        plot(time_FTphase,FFL_FTphase(:,nsampidx));
        Set_CCxC_plot_FigureProperties(gcf,gca,'Time (s)','GC FT',FS);
    end

end

end

function Set_CCxC_plot_FigureProperties(~,ah,Xstr,Ystr,FS,varargin)
% fh = figure handle
% ah = axis handle
argidx = 1;
axisflag=0;
semilogflag=0;
YLIM=[-inf inf];
XLIM=[-inf inf];
while argidx<=length(varargin)
    if strcmpi( varargin{argidx}, 'axis' )
        axisflag = 1;
        AXISstr=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'ylim' )
        YLIM=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'xlim' )
        XLIM=varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'semilog' )
        semilogflag=1;
        argidx=argidx+1;
    else
        argidx=argidx+1;
    end
end

xlim(XLIM);
ylim(YLIM);
xlabel(Xstr,'FontSize',FS);
ylabel(Ystr,'FontSize',FS);
if(axisflag)
    axis(ah, AXISstr);
end
set(ah, 'TickDir', 'out');
box off;
set(ah,'FontSize',FS);
if(semilogflag)
    set(gca,'yscale','log')
end
end
