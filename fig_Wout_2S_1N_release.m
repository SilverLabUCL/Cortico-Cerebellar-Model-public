function fig_Wout_2S_1N_release(varargin)

NS=200;
Nangle=[15,45,60,90];
plot1Dflag=0;
vn=1.1;
saveflag=0;
argidx = 1;
while argidx<=length(varargin)
    if strcmpi( varargin{argidx}, 'plot 1D' )
        plot1Dflag = 1;
        argidx=argidx+1;
    elseif strcmpi( varargin{argidx}, 'NS' )
        NS = varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'Nangle' )
        Nangle = varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'vn' )
        vn = varargin{argidx+1};
        argidx=argidx+2;
    elseif strcmpi( varargin{argidx}, 'save' )
        saveflag = 1;
        argidx=argidx+1;
    else
        argidx=argidx+1;
    end
end

%% calculate weights
Sangle=linspace(0.1,90,NS);

if(plot1Dflag)
    for k=1:numel(Nangle)
        Wnorm_sangle{k}=zeros(NS,1); %#ok<*AGROW>

        for i=1:NS
            Wnorm_sangle{k}(i) = Wout_2S_1N_simplified(cosd(Sangle(i)), cosd(Nangle(k)), vn, 1);
        end

    end
end

%% plots 1D
if(plot1Dflag)
    FS=10;
    LW=1.5;
    nc=1;
    nr=1;

    if(saveflag==1)
        sx=4.5;
        sy=sx;
        tilestr='tight';
        paddstr='tight';
        FS=10;
        LW=1;
    else
        sx=8;
        sy=sx;
        tilestr='compact';
        paddstr='compact';
    end

    set(0,'defaultAxesFontName', 'Arial');
    cmp = get(groot,'DefaultAxesColorOrder');

    f1D=figure;
    set(f1D,'papertype','A4','Units', 'centimeters', 'position', [10,40,nc*sx,nr*sy]);
    tiledlayout(nr,nc,'TileSpacing',tilestr,'Padding',paddstr);

    % Wnorm vs theta for different xi
    nexttile;
    hold on;
    lgdstr=cell(numel(Nangle),1);
    kk=1;
    for k=1:numel(Nangle)
        plot(Sangle,Wnorm_sangle{k},'-','LineWidth',LW,'Color',cmp(k,:));
        lgdstr{kk}=['\gamma=' num2str(Nangle(k)) ', ||v_{+}||=' num2str(vn)];
        kk=kk+1;
    end
    set(gca, 'yScale', 'log');
    xlim([min(Sangle) max(Sangle)]);
    xlabel('\theta (degrees)','FontSize',FS);
    ylabel('||W_{out}||','FontSize',FS);
    axis(gca, 'square'); set(gca, 'TickDir', 'out'); box off;
    set(gca,'FontSize',FS);
    set(gca, 'Color', 'None');
    legend(lgdstr);

    if(saveflag==1)
        filenamefull='Wout_2S_1N_1D';
        saveas(f1D,filenamefull,'svg');
    end
end

end

function [Wnorm, Wnorm_reg, cp] = Wout_2S_1N_simplified(cv,cz,n1,n2)
% cz is the principal angle
% cv is the signal angle

cp2=0.5 * cz^2 *(1-cv^2)/(1-cv);
cp=sqrt(cp2);

% check positive semi-definiteness
CR= [1 cv cp;
    cv 1 cp;
    cp cp 1];
eigvals = eig(CR);
if min(eigvals) < 0
    % disp('Matrix is not positive semi-definite.');
    Wnorm=NaN;
    Wnorm_reg=NaN;
    cp=NaN;
else
    alpha=1-cp2;
    beta=1-cp2;
    gamma=cv-cp2;

    cv2 = cv*cv;
    cz2 = cz*cz;

    Wnorm2=(beta*n2^2 +alpha*n1^2 -2*gamma*n1*n2)/(n1^2 *n2^2 *(1-cz2)*(1-cv2));

    Wnorm=sqrt(Wnorm2);
    Wnorm_reg=Wnorm *n1 *n2;
end

end