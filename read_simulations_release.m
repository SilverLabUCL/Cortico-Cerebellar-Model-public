function data=read_simulations_release(path_sim)

% path_sim='simulation_results/';

%% read json file for flags
fid=fopen([path_sim '/flags.json']);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
flags=jsondecode(str);

CORTEXONLY=flags.cortexonlyflag;
if(CORTEXONLY==1)
    disp('Data from cortex-only simulation.');
end

%% read json file for parameters
fid=fopen([path_sim '/parameters.json']);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
parameters_json=jsondecode(str);

Npara=length(parameters_json.Ncx);     % for readout before scan simulation is done

% remove redundant information
fn = fieldnames(parameters_json);
for k=1:numel(fn)
    if(iscell(parameters_json.(fn{k})))
    else
        data.parameters.(fn{k})=nanUnique_stable(parameters_json.(fn{k}));
    end
end
if(flags.learnflag==131 || flags.gainflag==4)

else
    data.parameters.GCgain=parameters_json.GCgain(1:length(data.parameters.GCthreshold));
end
data.parameters.phases=data.parameters.phases(1:data.parameters.Ncomp);
data.parameters.phases_intf=data.parameters.phases_intf(1:data.parameters.Ncomp_intf);
data.parameters.amplitudes=data.parameters.amplitudes(1:data.parameters.Ncomp);
data.parameters.amplitudes_intf=data.parameters.amplitudes_intf(1:data.parameters.Ncomp_intf);
data.parameters.flags=flags;

data.simcount=readmatrix([path_sim,'simcount.dat']);

%% single example data
if(flags.learnflag<3)
    data.nuRNN_test=readmatrix([path_sim,'nuRNN_test.dat']);
    data.nuRNN_openloop=readmatrix([path_sim,'nuRNN_openloop.dat']);

    if(CORTEXONLY==0)
        data.nuFFL_test=readmatrix([path_sim,'nuFFL_test.dat']);
        data.nuFFL_openloop=readmatrix([path_sim,'nuFFL_openloop.dat']);
    end

    path_para='parameter_1/';

    % read json file for stats
    fid=fopen([path_sim,path_para 'stats.json']);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    data.stats=jsondecode(str);

    % conver stats structs into tables
    data.stats=struct2table(data.stats);

    % read json file for targets
    fid=fopen([path_sim,path_para 'targets.json']);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    data.targets=jsondecode(str);
    data.time_FT{1}=data.targets.time_FT;
    data.time_learn{1}=data.targets.time_learn;
    data.time_test{1}=data.targets.time_test;

    data.target_FT{1}=data.targets.target_FT';
    data.target_learn{1}=data.targets.target_learn';
    data.target_test{1}=data.targets.target_test';
    if(flags.interferenceflag>0)
        data.target_FT_interf{1}=data.targets.target_FT_interf';
    end

    % test data
    if(flags.ROwriteflag==1)
        data.RNNoutput_test{1}=checkopen([path_sim,path_para,'RNNoutput_test.dat']);
        if(CORTEXONLY==0)
            data.FFLoutput_test{1}=checkopen([path_sim,path_para,'FFLoutput_test.dat']);
        end
    end

    % noise_test data
    if(flags.loopflag_test<=0)
        data.noise_test{1}=checkopen([path_sim,path_para,'noise_test.dat']);
    end

    % read json file for measures
    fid=fopen([path_sim,path_para 'measures.json']);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    data.measures_raw=jsondecode(str);

    % re-order data
    fn = fieldnames(data.measures_raw);
    for k=1:length(fn)
        tmp=data.measures_raw.(fn{k});
        data.(fn{k})=tmp;
    end

    % Read JSON file for single trial data
    fid = fopen([path_sim, path_para, 'weights_single_trial.json']);

    if fid > -1
        raw = fread(fid, inf);
        str = char(raw');
        fclose(fid);
        data.single_trial = jsondecode(str);
    else
        data.single_trial = [];
    end

    % Read spanning vectors
    data.stotal = checkopen([path_sim, 'stotal_RNN.dat']);
    if CORTEXONLY == 0
        data.stotal = checkopen([path_sim, 'stotal_FFL.dat']);
    end

    % Connectivity weights
    % Wcx
    Wcx = data.single_trial.Wcx;
    Wcx(:,1:2) = Wcx(:,1:2) + 1;
    data.weights.Wcx = spconvert(Wcx);

    % Wout_RNN
    data.weights.Wout_RNN = data.single_trial.Wout_RNN';

    % Wfb
    data.weights.Wfb = data.single_trial.Wfb;

    % Wfb_interf
    if data.parameters.flags.interferenceflag > 0
        data.weights.Wfb_interf = data.single_trial.Wfb_interf;
    end

    % FFL connectivity and parameters
    if CORTEXONLY == 0
        % MFidx
        data.weights.MFidx = data.single_trial.MFidx;
        data.weights.Wctxmf = spconvert(data.weights.MFidx);
        if size(data.weights.Wctxmf,2) ~= data.parameters.Ncx
            data.weights.Wctxmf(data.parameters.Nmf, data.parameters.Ncx) = 0;
        end

        % Wmf
        Wmf = data.single_trial.Wmf;
        for i = 1:size(Wmf,1)
            Wmf(i,2) = data.weights.MFidx(Wmf(i,2));
        end
        data.weights.Wmf = spconvert(Wmf);

        % Wout_FFL
        data.weights.Wout_FFL = data.single_trial.Wout_FFL';

        % Distributions and neuron parameters
        data.weights.Wgoc_distrb = data.single_trial.Wgoc_distrb';
        data.weights.theta = data.single_trial.theta';
        data.weights.gain = data.single_trial.gain';
        data.weights.tauGC_distrb = data.single_trial.tauGC_distrb';

    end

    % Initial conditions
    data.hcx_inistate = data.single_trial.hcx_inistate;


    % clean up
    data = rmfield(data,'measures_raw');
    data = rmfield(data,'targets');
    data = rmfield(data,'single_trial');
end

%% scan data
if(flags.learnflag>=3)
    for j=1:Npara
        path_Npara=['parameter_' num2str(j) '/'];

        % read json file for measures
        fid=fopen([path_sim,path_Npara 'measures.json']);
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);
        data.measures_raw{j,1}=jsondecode(str);

        % read json file for stats
        fid=fopen([path_sim,path_Npara 'stats.json']);
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);
        data.stats{j,1}=jsondecode(str);

        % read json file for targets
        fid=fopen([path_sim,path_Npara 'targets.json']);
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);
        data.targets=jsondecode(str);
        data.time_test{j,1}=data.targets.time_test;
        data.target_test{j,1}=data.targets.target_test';

    end

    % re-order data
    fn = fieldnames(data.measures_raw{1});
    for k=1:length(fn)
        tmp=[];
        for i=1:Npara
            tmp{i,1}=data.measures_raw{i}.(fn{k});  %#ok<AGROW>
        end
        data.(fn{k})=tmp;
    end

    % clean up
    data = rmfield(data,'measures_raw');
    data = rmfield(data,'targets');

    % conver stats structs into tables
    for i=1:length(data.stats)
        data.stats{i}=struct2table(data.stats{i});
    end

    % output signals
    if(flags.ROwriteflag==1)
        for j=1:Npara
            path_Npara=['parameter_' num2str(j) '/'];
            t=dir([path_sim,path_Npara,'RNNoutput_test.dat']);
            if(~isempty(t))
                fname='RNNoutput_test.dat';
                data.RNNoutput_test{j}=readmatrix([path_sim,path_Npara,fname])';
            end
            % noise_test data
            if(flags.loopflag_test<=0)
                t=dir([path_sim,path_Npara,'noise_test.dat']);
                if(t.bytes~=0)
                    fname='noise_test.dat';
                    data.noise_test{j,1}=readmatrix([path_sim,path_Npara,fname])';
                end
            end
        end

        if(CORTEXONLY==0)
            % output signals
            for j=1:Npara
                path_Npara=['parameter_' num2str(j) '/'];
                s=dir([path_sim,path_Npara,'FFLoutput_test.dat']);
                if(s.bytes~=0)
                    fname='FFLoutput_test.dat';
                    data.FFLoutput_test{j}=readmatrix([path_sim,path_Npara,fname])';
                end
            end
        end
    end

end

%% reorder angles into matrices
if(isfield(data,'Angles_ctx'))
    data.Angles_ctx=reorder_angles(data.Angles_ctx,Npara);
% else
%     data.Angles_ctx=reorder_angles(data.Angles_cereb,Npara);
end

if(isfield(data,'Angles_cereb'))
    data.Angles_cereb=reorder_angles(data.Angles_cereb,Npara);
end

%% reorder R_full
if(isfield(data,'R_cereb'))
    tmp=cell(Npara,1);
    Nmat=sqrt(length(data.R_cereb{1}(1,:)));
    for i=1:Npara
        for j=1:size(data.R_cereb{i},1)
            tmp{i}(:,:,j)=reshape(data.R_cereb{i}(j,:),Nmat,Nmat);
        end
    end
    data.R_cereb=tmp;
end

%% read in target signal parameters
if(flags.signalflag==0 || flags.signalflag==1 || flags.signalflag==2 || flags.signalflag==22 || flags.signalflag==23 ||...
        flags.signalflag==24 || flags.signalflag==25 || flags.signalflag==99)

    data.frequencies = parameters_json.frequency;
    data.amplitudes = parameters_json.amp;
    data.phases = parameters_json.phase;

    if(parameters_json.Ncomp_intf>0)
        data.frequencies_intf = parameters_json.frequency_intf;
        data.amplitudes_intf = parameters_json.amp_intf;
        data.phases_intf = parameters_json.phase_intf;
    end

end

%% eliminate unused parameters
if(CORTEXONLY==1)
    data.parameters.CORTEXONLY='reduced simulation';
    flags.actfuncflag_GC=[];
    flags.fcodingflag=[];
    flags.gainflag=[];
    flags.gocmodeflag=[];
    flags.filterflag=[];
    flags.simplemapflag=[];
    flags.taugcmodeflag=[];
    flags.STPcerebflag=[];

    data.parameters.FACsigmaLS_gc=[];
    data.parameters.fc_avrg=[];
    data.parameters.fc_sd=[];
    data.parameters.GCgain=[];
    data.parameters.GCmax=[];
    data.parameters.GCgainmax=[];    
    data.parameters.GCthreshold=[];
    data.parameters.GCthreshold_base=[];
    data.parameters.Ngc=[];
    data.parameters.Nmf=[];
    data.parameters.Nsyn=[];
    data.parameters.sigmaNgc=[];
    data.parameters.tauGC_avrg=[];
    data.parameters.tauGC_std=[];
    data.parameters.tauGCcutoff=[];
    data.parameters.tol_GCgain=[];    
    data.parameters.Wgoc=[];
    data.parameters.Wgoc_SD=[];
else
    data.parameters.CORTEXONLY='full simulation';
end
data.parameters.flags=flags;

%% identify target signal
if(flags.signalflag==0)
    data.parameters.signaltype='quadruple sine';
elseif(flags.signalflag==22)
    if(flags.learnflag==107)
        data.parameters.signaltype=['sine wave with ' num2str(data.parameters.Ncomp(1)) ' to ' num2str(data.parameters.Ncomp(end)) ' components'];
    else
        data.parameters.signaltype=['sine wave with ' num2str(data.parameters.Ncomp) ' components'];
    end
elseif(flags.signalflag==23)
    if(flags.learnflag==107)
        data.parameters.signaltype=['cosine wave with ' num2str(data.parameters.Ncomp(1)) ' to ' num2str(data.parameters.Ncomp(end)) ' components'];
    else
        data.parameters.signaltype=['cosine wave with ' num2str(data.parameters.Ncomp) ' components'];
    end
elseif(flags.signalflag==24)
    if(flags.learnflag==107)
        data.parameters.signaltype=['cosine-sine wave with ' num2str(data.parameters.Ncomp(1)) ' to ' num2str(data.parameters.Ncomp(end)) ' components'];
    else
        data.parameters.signaltype=['cosine-sine wave with ' num2str(data.parameters.Ncomp) ' components'];
    end
elseif(flags.signalflag==99)
    data.parameters.signaltype='test signal';
end
%% identify noise mode
if(flags.noiseflag==0)
    data.parameters.noisetype='constant noise on ctx and gc inputs';
elseif(flags.noiseflag==1)
    data.parameters.noisetype='noise scaled with ctx and gc inputs';
elseif(flags.noiseflag==2)
    data.parameters.noisetype='noise scaled with ctx and gc rates';
elseif(flags.noiseflag==-1)
    data.parameters.noisetype='no noise';
end
%% identify learning rule
if(flags.LRflag==0)
    data.parameters.learningrule='RLS';
    data.parameters.sigmaLS=[];
    data.parameters.delta_eta=[];
    data.parameters.delta_r2_reg=[];
elseif(flags.LRflag==1)
    data.parameters.learningrule='DELTA';
    data.parameters.sigmaLS=[];
elseif(flags.LRflag==2)
    data.parameters.learningrule='LS';
    data.parameters.delta_eta=[];
    data.parameters.delta_r2_reg=[];
elseif(flags.LRflag==3)
    data.parameters.learningrule='RR';
    data.parameters.delta_eta=[];
    data.parameters.delta_r2_reg=[];
elseif(flags.LRflag==4)
    data.parameters.learningrule='PINV';
    data.parameters.delta_eta=[];
    data.parameters.delta_r2_reg=[];
end
%% identify CC input layer connectivity
if(flags.simplemapflag==0)
    data.parameters.MFtoGC='Nsyn mixed MF inputs per GC';
elseif(flags.simplemapflag==1)
    data.parameters.MFtoGC='1 unique MF per GC';
end
%% identify RNN activation function
if(flags.actfuncflag==0)
    data.parameters.activation_function='tanh';
elseif(flags.actfuncflag==1)
    data.parameters.activation_function='linear';
end
%% identify GC activation function
if(flags.actfuncflag_GC==0)
    data.parameters.activation_function_GC='relu';
    data.parameters.sp_alpha=[];
elseif(flags.actfuncflag_GC==1)
    data.parameters.activation_function_GC='softplus';
elseif(flags.actfuncflag_GC==2)
    data.parameters.activation_function_GC='linear';
    data.parameters.sp_alpha=[];
elseif(flags.actfuncflag_GC==3)
    data.parameters.activation_function_GC='tanh';
    data.parameters.sp_alpha=[];
end
%% identify loop mode of learning phase
if(flags.loopflag==0)
    data.parameters.loop='open loop learning';
elseif(flags.loopflag==1)
    data.parameters.loop='closed loop learning';
    data.parameters.openloop_test_noise=[];
end
%% identify loop mode of testphase
if(flags.loopflag_test==0)
    data.parameters.loop_test='open loop testphase';
elseif(flags.loopflag_test==1)
    data.parameters.loop_test='closed loop testphase';
elseif(flags.loopflag_test==-1)
    data.parameters.loop_test='open loop testphase with white noise';
end
%% identify interference mode
if(flags.interferenceflag==0)
    data.parameters.interference='OFF';
    % data.parameters.rho_Wfb=[];
    % data.parameters.rho_Wfb12=[];
    data.parameters.rho_Wfb13=[];
    data.parameters.sigAmp_intf=[];
elseif(flags.interferenceflag==1)
    data.parameters.interference='ON';
elseif(flags.interferenceflag==2)
    data.parameters.interference='GCs';
end
%% identify threshold mode
if(flags.fcodingflag==0)
    data.parameters.threshold_mode='Adjusting thresholds by setting common value for theta';
elseif(flags.fcodingflag==1)
    data.parameters.threshold_mode='Adjusting GC coding level globally';
elseif(flags.fcodingflag==2)
    data.parameters.threshold_mode='Adjusting GC coding level individually';
end
%% identify gain mode
if(flags.gainflag==0)
    data.parameters.gain_mode='Gains set globally via theory';
elseif(flags.gainflag==1)
    data.parameters.gain_mode='Gains set globally based on population SD';
elseif(flags.gainflag==2)
    data.parameters.gain_mode='Gains set locally based on individual SD';
elseif(flags.gainflag==3)
    data.parameters.gain_mode='Gains fixed globally';
elseif(flags.gainflag==4)
    data.parameters.gain_mode='Gains set to inverse norm of 1st spanning vector';
elseif(flags.gainflag==5)
    data.parameters.gain_mode='Gains set globally such as to match 1st spanning vector norms';
elseif(flags.gainflag==6)
    data.parameters.gain_mode='Gains set locally (with identical SD) such as to match 1st spanning vector norms';
elseif(flags.gainflag==7)
    data.parameters.gain_mode=['Gains set locally to adjust maximum GC rate to GCmax =' num2str(data.parameters.GCmax)];
elseif(flags.gainflag==8)
    data.parameters.gain_mode='Gains set locally to match maximum RNN and GC rates';
end
%% identify GoC mode
if(flags.gocmodeflag==0)
    data.parameters.GoC_mode='GoC OFF';
    data.parameters.Wgoc=[];
    data.parameters.Wgoc_SD=[];
elseif(flags.gocmodeflag==1)
    data.parameters.GoC_mode='GoC ON';
    data.parameters.Wgoc_SD=[];
elseif(flags.gocmodeflag==2)
    data.parameters.GoC_mode='GoC ON with heterogeneous weights';
elseif(flags.gocmodeflag==3)
    data.parameters.GoC_mode='GoC ON based on MF SD';
    data.parameters.Wgoc_SD=[];
end
%% identify RNN tau distribution mode
if(flags.taudistrbflag==0)
    data.parameters.tau_mode='uniform tau';
    data.parameters.tau_std=[];
elseif(flags.taudistrbflag==1)
    data.parameters.tau_mode='distributed tau';
end
%% identify GC tau mode
if(flags.taugcmodeflag==0)
    data.parameters.tauGC_mode='GC tau OFF';
    data.parameters.tauGC_avrg=[];
    data.parameters.tauGC_std=[];
elseif(flags.taugcmodeflag==1)
    data.parameters.tauGC_mode='uniform GC tau';
    data.parameters.tauGC_std=[];
elseif(flags.taugcmodeflag==2)
    data.parameters.tauGC_mode='distributed GC tau';
end
%% identify dynamical FM cleaning removal mode
if(isfield(flags,'FMdyncleanflag'))
    if(flags.FMdyncleanflag==0)
        data.parameters.FM_cleaning='standard';
    elseif(flags.FMdyncleanflag==1)
        data.parameters.FM_cleaning='dynamic';
    end
end
%% identify Floquet multiplier mode
if(flags.floquetflag==1)
    if(flags.floquetmodeflag==0)
        data.parameters.Floquet='open loop Floquet multipliers calculated';
    else
        data.parameters.Floquet='closed loop Floquet multipliers calculated';
    end
else
    data.parameters.Floquet='no Floquet multipliers calculated';
end

%% identify learning mode
if(flags.learnflag==1)
    data.parameters.simtype='sample simulation, cortical learning';
elseif(flags.learnflag==2)
    data.parameters.simtype='sample simulation, cerebellar learning';
elseif(flags.learnflag==0)
    data.parameters.simtype='sample simulation, no learning';
elseif(flags.learnflag==99)
    data.parameters.simtype='learning repeated over Nreal realisations';
elseif(flags.learnflag==100)
    data.parameters.simtype='scan over sigmaN';
elseif(flags.learnflag==101)
    data.parameters.simtype='scan over GC threshold';
elseif(flags.learnflag==102)
    data.parameters.simtype='scan over Nctx';
elseif(flags.learnflag==103)
    data.parameters.simtype='scan over Ngc';
elseif(flags.learnflag==105)
    data.parameters.simtype='scan over sigmaNgc';
elseif(flags.learnflag==106)
    data.parameters.simtype='scan over Wgoc';
elseif(flags.learnflag==107)
    data.parameters.simtype='scan over Ncomp';
elseif(flags.learnflag==108)
    data.parameters.simtype='scan over sigAmp';
elseif(flags.learnflag==109)
    data.parameters.simtype='scan over lambda';
elseif(flags.learnflag==110)
    data.parameters.simtype='scan over Tlearn';
elseif(flags.learnflag==112)
    data.parameters.simtype='scan over p_{recip}';
elseif(flags.learnflag==113)
    data.parameters.simtype='scan over g';
elseif(flags.learnflag==114)
    data.parameters.simtype='scan over eta';
elseif(flags.learnflag==118)
    data.parameters.simtype='scan over Phi';
elseif(flags.learnflag==119)
    data.parameters.simtype='scan over Ampmod';
elseif(flags.learnflag==120)
    data.parameters.simtype='scan over sigmaLS';
elseif(flags.learnflag==123)
    data.parameters.simtype='scan over average f_{coding}';
elseif(flags.learnflag==124)
    data.parameters.simtype='scan over sopftpus sharpness';
elseif(flags.learnflag==128)
    data.parameters.simtype='scan over feedback-weight correlation';
elseif(flags.learnflag==129)
    data.parameters.simtype='scan over interference signal amplitude';
elseif(flags.learnflag==131)
    data.parameters.simtype='scan over (fixed) GC gain';
end

end


function out=checkopen(filename)

if isfile(filename)
    out=readmatrix(filename);
else
    out=[];
end

end

function y = nanUnique_stable(x)
y = unique(x,'stable','row');
y(isnan(y(1:end-1))) = [];
end

function out=reorder_angles(Angles,Npara)
out=cell(Npara,1);
if(iscell(Angles))
    anglraw=Angles{1};
elseif(Npara==1)
    tmp=Angles;
    Angles=cell(1,1);
    Angles{1}=tmp;
    anglraw=tmp;
else
    error('Wrong angles format.')
end
aidx=find(all(anglraw==0,1));
for j=2:Npara
    anglraw=Angles{j};
    aidx=intersect(aidx,find(all(anglraw==0,1)));
end
for j=1:Npara
    anglraw=Angles{j};
    % aidx=find(all(anglraw==0));
    for k=1:size(anglraw,1)
        % aidx=find(anglraw(k,:)==0);
        mang=zeros(length(aidx),length(aidx)+1);
        for i=1:(length(aidx)-1)
            tmp=anglraw(k,aidx(i):(aidx(i+1)-1));
            mang(i+1,1:length(tmp))=tmp;
        end
        mang(:,1)=[];
        out{j,1}(:,:,k)=mang;
    end
end
end
