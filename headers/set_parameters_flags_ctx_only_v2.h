int set_parameters_flags(netpara &para, int &argc, char* argv[]) {

//------------------------------------------------------------------------
// 	default parameters
//------------------------------------------------------------------------
  	para.DT = 0.001;			// integration time constant
  	para.DT_euler = 0.0001;		// integration time constant for euler test-phase
	para.Ncx=500;				// number of RNN units
	para.lambda=5;				// signal fundamental frequency
	para.Ncomp=1;				// number of components of combined sine/cosine signal
	para.Ncomp_intf=2;			// number of components of interference signal
	para.nf=1;					// spacing between frequency components (nf=1: every overtones, nf=2: every second overtone, etc.)
	para.g=0.9;
	para.seed=-1;				// default random number seeds
	para.Nreal=10;				// number of network realisations
	para.Npara=10;				// number of parameter steps
	para.N_FMperiods=10;		// number of periods for estimation of Floquet multipliers
	para.duration_learn=10;		// length of learning phase (in periods)
	para.duration_test=20;		// length of test phase (in periods)
	para.sigmaLS=pow(10,-2);	// regularisation parameter for LS or RR learning
	para.sigmaN=5*pow(10,-3);	// noise std added to cortical units
	para.rho_Wfb12=1;			//
	para.rho_Wfb13=1;			//
	para.write_every=20;		// read out every nth time step
	para.outlier_thresh=1.1;    // factor determining outlier detection threshold

	para.delta_eta=0.0025;		// initial learning rate for DELTA rule
	para.delta_eps_decay=1e-5;	// weight decay parameter for DELTA rule
	para.delta_r2_reg=1e-6;		// norm regularisation for DELTA rule

	para.INT_delta_abs=1e-3;	// integrator absolute error
	para.INT_delta_rel=1e-3;	// integrator relative error

	para.sigAmp=1.0;		// global amplitude multiplier
	para.sigAmp_intf=1.0;	// global amplitude multiplier (interference signal)

	para.amps={1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};	// all amplitudes of signal / interference signals
	//~ para.amps=para.amps*sqrt(2);
	para.amps_intf={0.0,1.0,1.0,1.0,1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

	para.phases={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};		// all phases of signal / interference signals
	para.phases_intf={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

	// target signal: 'quadruple_sine', 'accordion', 'sine', 'cosine_fixed', 'lorenz', 'square_wave', 'sawtooth_wave', 'staircase_wave', 'test_signal' or 'zero_signal'
	para.signaltype="cosine_fixed";
	// scale noise with firing rate: 'constant', 'proportional', 'rate' or 'none'
	para.noisescaling="none";
	// learn mode: 'cortex', 'none', 'repeat' or 'scan_...'
	para.learnmode="cortex";
	// learn mode: 'scan_...'
	para.learnmode2="scan_lambda";
	// learning rule: 'RLS', 'LS', 'RR', 'PINV' or 'DELTA'
	para.learningrule="LS";
	// testphase in 'closed' or 'open' loop or with open loop with noise exclusively: 'open noise'
	para.loop_test="closed";
	// activation function of reservoir units: 'tanh', 'linear' or 'relu'
	para.activation_function="tanh";
	// add interference signal: 'ON' or 'OFF'
	para.interference_signal="OFF";
	// cortical tau distribution: 'ON' or 'OFF'
	para.tau_mode="OFF";
	// calculate Floquet multipiers: 'ON' or 'OFF'
	para.floquet="OFF";
	// calculate Floquet multipiers: 'closed' or 'open'
	para.floquetmode="open";
	// write readout unit output to file
	para.ROwrite="OFF";
	// dynamic cleanup of Fourier modes
	para.FMdynclean="OFF";

	try {
		po::options_description desc("Allowed options", 240);
		desc.add_options()
			("help", "produce help message")
			("signaltype", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.signaltype), "target signal: 'quadruple_sine', 'accordion', 'cosine_fixed'... ")
			("learnmode", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.learnmode), "learn mode: 'cortex', 'none', 'repeat' or 'scan_...' ")
			("learnmode2", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.learnmode2), "learn mode: 'scan_...' ")
			("LRule", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.learningrule), "learning rule: 'RLS', 'LS', 'RR', 'PINV' or 'DELTA' ")
			("noise", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.noisescaling), "noise (scaling with firing rate): 'constant', 'proportional' or 'none' ")
			("loop_test", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.loop_test), "loop type during test-phase: 'closed' or 'open'")
			("activation_function", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.activation_function), "activation function of RNN units: 'tanh' or 'linear' ")
			("interference_signal", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.interference_signal), "add interference signal: 'ON' or 'OFF' ")
			("tau_mode", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.tau_mode), "distribute RNN taus: 'ON' or 'OFF' ")
			("Floquet", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.floquet), "calculate Floquet multipiers: 'ON' or 'OFF' ")
			("loop_Floquet_mode", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.floquetmode), "Floquet mode: 'closed' or 'open' ")
			("ROwrite", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.ROwrite), "write reaout to file: 'ON' or 'OFF' ")
			("FMdynclean", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.FMdynclean), "dynamically clean Fourier spectrum: 'ON' or 'OFF' ")
			("Ncx", po::value<int>(&para.Ncx)->default_value(para.Ncx), "number of RNN units")
			("Ncomp", po::value<int>(&para.Ncomp)->default_value(para.Ncomp), "number of signal frequency components")
			("Ncomp_intf", po::value<int>(&para.Ncomp_intf)->default_value(para.Ncomp_intf), "number of interference signal frequency components")
			("nf", po::value<int>(&para.nf)->default_value(para.nf), "spacing between frequency components")
			("g", po::value<double>(&para.g)->default_value(para.g), "reccurent weight parameter")
			("lambda", po::value<double>(&para.lambda)->default_value(para.lambda), "signal fundamental frequency")
			("T_learn", po::value<double>(&para.duration_learn)->default_value(para.duration_learn), "length of learning phase (i.e. number of periods)")
			("T_test", po::value<double>(&para.duration_test)->default_value(para.duration_test), "length of test phase  (i.e. number of periods)")
			("sigmaLS", po::value<double>(&para.sigmaLS)->default_value(para.sigmaLS), "regularisation for LS or RR")
			("sigmaN", po::value<double>(&para.sigmaN)->default_value(para.sigmaN), "RNN unit noise")
			("sigAmp", po::value<double>(&para.sigAmp)->default_value(para.sigAmp), "global target signal multiplier")
			("sigAmp_intf", po::value<double>(&para.sigAmp_intf)->default_value(para.sigAmp_intf), "global interference signal multiplier")
			("rho_Wfb12", po::value<double>(&para.rho_Wfb12)->default_value(para.rho_Wfb12), "correlation rho_Wfb12")
			("rho_Wfb13", po::value<double>(&para.rho_Wfb13)->default_value(para.rho_Wfb13), "correlation rho_Wfb13")
			("Nreal", po::value<int>(&para.Nreal)->default_value(para.Nreal), "number of network realisations")
			("Npara", po::value<int>(&para.Npara)->default_value(para.Npara), "number of parameter steps")
			("N_FMperiods", po::value<int>(&para.N_FMperiods)->default_value(para.N_FMperiods), "number of periods for calculation of Floquet multiplier")
			("seed", po::value<int>(&para.seed)->default_value(para.seed), "random number seed")
			("write_every", po::value<int>(&para.write_every)->default_value(para.write_every), "read-out sampling rate (every xth time-step)")
			("dt", po::value<double>(&para.DT)->default_value(para.DT), "integration time step")
			("dt_euler", po::value<double>(&para.DT_euler)->default_value(para.DT_euler), "integration time step of euler test-phase")
			("outlier_threshold", po::value<double>(&para.outlier_thresh)->default_value(para.outlier_thresh), "factor determining outlier detection threshold")
			("delta_eta", po::value<double>(&para.delta_eta)->default_value(para.delta_eta), "initial learning rate for DELTA rule")
			("INT_delta_abs", po::value<double>(&para.INT_delta_abs)->default_value(para.INT_delta_abs), "target absolute error for numerical integration")
			("INT_delta_rel", po::value<double>(&para.INT_delta_rel)->default_value(para.INT_delta_rel), "target relative error for numerical integration")

		;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			std::cout << desc << "\n";
			return 1;
		}
		std::vector<std::string> tmp;
		if (vm["signaltype"].defaulted()==0)  {
			tmp=vm["signaltype"].as< std::vector<std::string> >();
			para.signaltype= tmp[0];
		}
		if (vm["learnmode"].defaulted()==0)  {
			tmp=vm["learnmode"].as< std::vector<std::string> >();
			para.learnmode= tmp[0];
		}
		if (vm["learnmode2"].defaulted()==0)  {
			tmp=vm["learnmode2"].as< std::vector<std::string> >();
			para.learnmode2= tmp[0];
		}
		if (vm["LRule"].defaulted()==0)  {
			tmp=vm["LRule"].as< std::vector<std::string> >();
			para.learningrule= tmp[0];
		}
		if (vm["noise"].defaulted()==0)  {
			tmp=vm["noise"].as< std::vector<std::string> >();
			para.noisescaling= tmp[0];
		}
		if (vm["loop_test"].defaulted()==0)  {
			tmp=vm["loop_test"].as< std::vector<std::string> >();
			para.loop_test= tmp[0];
		}
		if (vm["activation_function"].defaulted()==0)  {
			tmp=vm["activation_function"].as< std::vector<std::string> >();
			para.activation_function= tmp[0];
		}
		if (vm["interference_signal"].defaulted()==0)  {
			tmp=vm["interference_signal"].as< std::vector<std::string> >();
			para.interference_signal= tmp[0];
		}
		if (vm["tau_mode"].defaulted()==0)  {
			tmp=vm["tau_mode"].as< std::vector<std::string> >();
			para.tau_mode= tmp[0];
		}
		if (vm["Floquet"].defaulted()==0)  {
			tmp=vm["Floquet"].as< std::vector<std::string> >();
			para.floquet= tmp[0];
		}
		if (vm["loop_Floquet_mode"].defaulted()==0)  {
			tmp=vm["loop_Floquet_mode"].as< std::vector<std::string> >();
			para.floquetmode= tmp[0];
		}
		if (vm["ROwrite"].defaulted()==0)  {
			tmp=vm["ROwrite"].as< std::vector<std::string> >();
			para.ROwrite= tmp[0];
		}
		if (vm["FMdynclean"].defaulted()==0)  {
			tmp=vm["FMdynclean"].as< std::vector<std::string> >();
			para.FMdynclean= tmp[0];
		}
	}
	catch(std::exception& e) {
		std::cerr << "error: " << e.what() << std::endl;
		return 1;
	}
	catch(...) {
		std::cerr << "Exception of unknown type." << std::endl;
	}

	// dummy variables
	para.MF_GC_connect="simple";
	para.fcodingmode="global";
	para.gainmode="fixed";
	para.learnmode2="cortex";
	para.Npara2=1;
	para.activation_function_GC="relu";
	para.cortexonlyflag=1;
	para.Ngc=2;

	para.fcodingflag=0;
	para.filterflag=0;

	// target folders
	para.path="sim/";
	para.path2D="sim2D/";
	para.path_connect="connectivity/";
	//~ para.path="sim2/";
	//~ para.path_connect="connectivity2/";

//------------------------------------------------------------------------
// 	set fixed parameters
//------------------------------------------------------------------------

	para.RLSlearn_every=2;					// learn at every nth time step (RLS mode)

	para.duration_FT=3;						// time-window for FT in seconds
	para.duration_trans=15;					// time-window for transient activity in seconds

	para.Ncomp_noise=60;					// number of FT noise components to be extracted
	para.Ncomp_noise_RNN=60;				// number of FT noise components to be extracted

	para.FCcomp_tol=1e-6;					// cutoff for FT component entries

	para.Nwrite= 10; 					// number of components to be written to otput

	// correlation between feedback weights of target and interference signal
	para.rho_Wfb.ones(para.Ncomp_intf+1,para.Ncomp_intf+1);
	para.rho_Wfb.ones(para.Ncomp_intf+1,para.Ncomp_intf+1);

	para.amps.ones(para.Ncomp);
	para.amps_intf.ones(para.Ncomp_intf);
	para.amps_intf(0)=0.0;
	para.phases.zeros(para.Ncomp);
	para.phases_intf.zeros(para.Ncomp_intf);

	//~ para.rho_Wfb={ 	{1, 1, 1}, 				// 3x3 	for N_intf=2
					//~ {1, 1, 1},
					//~ {1, 1, 1} };

	//~ para.rho_Wfb={ 	{1, 1, 1, 1}, 			// 4x4	for N_intf=3
					//~ {1, 1, 1, 1},
					//~ {1, 1, 1, 1},
					//~ {1, 1, 1, 1} };

	//~ para.rho_Wfb={ 	{1, 1, 1, 1, 1}, 		// 5x5
					//~ {1, 1, 1, 1, 1},
					//~ {1, 1, 1, 1, 1},
					//~ {1, 1, 1, 1, 1},
					//~ {1, 1, 1, 1, 1} };
	//~ para.rho_Wfb12=0.5;
	//~ para.rho_Wfb12=para.rho_Wfb(0,2);
	if(para.Ncomp_intf>1){
		para.rho_Wfb(0,2)=para.rho_Wfb12;
		para.rho_Wfb(2,0)=para.rho_Wfb12;
	}
	if(para.Ncomp_intf>2){
		para.rho_Wfb(0,3)=para.rho_Wfb13;
		para.rho_Wfb(3,0)=para.rho_Wfb13;
	}

	para.tau_avrg=0.01;		// average and standard deviation of integration time constant
	para.tau_std=0.01;
	para.tauRNNcutoff=0.003;		// smallest permitted tau

	para.pcon=0.1;				// connection probability between cortical units

	para.phi1=0.0;			// phase of 1st cosine/sine component  	//~ para.phi1=3*M_PI/4.0;
	para.ampmod=0;			// modification of 1st & 2nd cosine components

	para.alpha=1;			// regularisation term of P matrix (for RLS)

	para.sqNsyn=sqrt(para.Nsyn);

//------------------------------------------------------------------------
// set all flags
//------------------------------------------------------------------------
	if 		(para.learnmode == "cortex") 		para.learnflag=1;
    else if (para.learnmode == "none") 			para.learnflag=0;
    else if (para.learnmode == "repeat") 		para.learnflag=99;
    else if (para.learnmode == "scan_sigmaN") 	para.learnflag=100;
    else if (para.learnmode == "scan_Nctx") 	para.learnflag=102;
    else if (para.learnmode == "scan_Ncomp") 	para.learnflag=107;
    else if (para.learnmode == "scan_sigAmp") 	para.learnflag=108;
    else if (para.learnmode == "scan_lambda") 	para.learnflag=109;
    else if (para.learnmode == "scan_Tlearn") 	para.learnflag=110;
    else if (para.learnmode == "scan_g") 		para.learnflag=113;
    else if (para.learnmode == "scan_eta") 		para.learnflag=114;
    else if (para.learnmode == "scan_Phi") 		para.learnflag=118;
    else if (para.learnmode == "scan_Ampmod") 	para.learnflag=119;
    else if (para.learnmode == "scan_sigmaLS") 	para.learnflag=120;
    else if (para.learnmode == "scan_rho_Wfb12") 	para.learnflag=128;
    else if (para.learnmode == "scan_sigAmp_intf") para.learnflag=129;
    else if (para.learnmode == "scan_rho_Wfb13") 	para.learnflag=132;
    else if (para.learnmode == "scan_RNNthreshold") para.learnflag=133;
    else {std::cerr <<"Bad learnmode string. \n";exit(EXIT_FAILURE);}

    if 		(para.learnmode2 == "cortex") 		para.learnflag2=1;
    else if (para.learnmode2 == "none") 		para.learnflag2=0;
    else if (para.learnmode2 == "repeat") 		para.learnflag2=99;
    else if (para.learnmode2 == "scan_sigmaN") 	para.learnflag2=100;
    else if (para.learnmode2 == "scan_Nctx") 	para.learnflag2=102;
    else if (para.learnmode2 == "scan_Ncomp") 	para.learnflag2=107;
    else if (para.learnmode2 == "scan_sigAmp") 	para.learnflag2=108;
    else if (para.learnmode2 == "scan_lambda") 	para.learnflag2=109;
    else if (para.learnmode2 == "scan_Tlearn") 	para.learnflag2=110;
    else if (para.learnmode2 == "scan_g") 		para.learnflag2=113;
    else if (para.learnmode2 == "scan_eta") 	para.learnflag2=114;
    else if (para.learnmode2 == "scan_Phi") 	para.learnflag2=118;
    else if (para.learnmode2 == "scan_Ampmod") 	para.learnflag2=119;
    else if (para.learnmode2 == "scan_sigmaLS") para.learnflag2=120;
    else if (para.learnmode2 == "scan_rho_Wfb12") para.learnflag2=128;
    else if (para.learnmode2 == "scan_sigAmp_intf") para.learnflag2=129;
    else if (para.learnmode2 == "scan_rho_Wfb13") 	para.learnflag2=132;
    else if (para.learnmode2 == "scan_RNNthreshold") para.learnflag2=133;
    else {std::cerr <<"Bad learnmode2 string. \n";exit(EXIT_FAILURE);}

    if 		(para.signaltype == "quadruple_sine") 	para.signalflag=0;
    else if (para.signaltype == "sine") 			para.signalflag=22;
    else if (para.signaltype == "cosine_fixed") 	para.signalflag=23;
    else if (para.signaltype == "cosine_sine") 		para.signalflag=24;
    else if (para.signaltype == "zero_signal") 		para.signalflag=-1;
    else {std::cerr <<"Bad signaltype string. \n";exit(EXIT_FAILURE);}

	if 		(para.noisescaling == "none") 			para.noiseflag=-1;
    else if (para.noisescaling == "constant") 		para.noiseflag=0;
    else if (para.noisescaling == "proportional") 	para.noiseflag=1;
    //~ else if (para.noisescaling == "rate") 			para.noiseflag=2;
    else {std::cerr <<"Bad noisescaling string. \n";exit(EXIT_FAILURE);}

	if 		(para.learningrule == "RLS") 	para.LRflag=0;
    else if (para.learningrule == "DELTA") 	para.LRflag=1;
    else if (para.learningrule == "LS") 	para.LRflag=2;
    else if (para.learningrule == "RR") 	para.LRflag=3;
    else if (para.learningrule == "PINV") 	para.LRflag=4;
    else {std::cerr <<"Bad LRstring string. \n";exit(EXIT_FAILURE);}

    if 		(para.interference_signal == "OFF") para.intrfflag=0;
    else if (para.interference_signal == "ON") 	para.intrfflag=1;
    else {std::cerr <<"Bad interference_signal string. \n";exit(EXIT_FAILURE);}

    if (para.activation_function == "tanh") 		para.actfuncflag=0;
    else if (para.activation_function == "linear") 	para.actfuncflag=1;
	else {std::cerr <<"Bad activation_function string. \n";exit(EXIT_FAILURE);}

	if (para.tau_mode== "OFF") 						para.taudistrbflag=0;
    else if (para.tau_mode == "ON") 				para.taudistrbflag=1;
	else {std::cerr <<"Bad tau_mode string. \n";exit(EXIT_FAILURE);}

	if (para.floquet == "OFF") 						para.floquetflag=0;
    else if (para.floquet == "ON") 					para.floquetflag=1;
	else {std::cerr <<"Bad floquet string. \n";exit(EXIT_FAILURE);}

	if (para.floquetmode == "open") 				para.floquetmodeflag=0;
    else if (para.floquetmode == "closed") 			para.floquetmodeflag=1;
	else {std::cerr <<"Bad floquetmode string. \n";exit(EXIT_FAILURE);}

	if (para.ROwrite == "OFF") 						para.ROwriteflag=0;
    else if (para.ROwrite == "ON") 					para.ROwriteflag=1;
	else {std::cerr <<"Bad ROwrite string. \n";exit(EXIT_FAILURE);}

	if (para.FMdynclean == "OFF") 					para.FMdyncleanflag=0;
    else if (para.FMdynclean == "ON") 				para.FMdyncleanflag=1;
	else {std::cerr <<"Bad FMdynclean string. \n";exit(EXIT_FAILURE);}

	if (para.loop_test == "open") 					para.loop_test_flag=0;
    else if (para.loop_test == "closed") 			para.loop_test_flag=1;
	else {std::cerr <<"Bad loop_test string. \n";exit(EXIT_FAILURE);}

//------------------------------------------------------------------------
// resolve conflicting options
//------------------------------------------------------------------------
	if (para.signalflag==0) para.Ncomp=4;

	if(para.learnflag<99)	{
		para.Nreal=1;
		para.ROwriteflag=1;
	}
    para.Nrep=2;
	if(para.noiseflag==-1){
		para.sigmaNgc=para.sigmaN=0;
	}
	if(para.learnflag==118){
		// if target signal is neither combination of sines nor cosines set to cosine default
		if(para.signalflag==22 || para.signalflag==23 || para.signalflag==24) {
		}
		else{
			para.signaltype="cosine_fixed";
			para.signalflag=23;
		}
	}
	if(para.learnflag==120){
		if(para.LRflag<2){
			para.learningrule="LS";
			para.LRflag=2;
		}
	}

	//~ if( (para.learnflag==128 || para.learnflag2==128 || para.learnflag==129 || para.learnflag2==129) && para.intrfflag==0 ){
		//~ std::cerr << "Learn- and interference-modes incompatible." << std::endl;
		//~ std::cerr << std::endl;
		//~ exit(EXIT_FAILURE);
	//~ }

	// check for incompatible interference signal parameters
	if(para.intrfflag>0 && !(para.signalflag==22 || para.signalflag==23 || para.signalflag==24 || para.signalflag==25) ){
		std::cerr << "Target signal and interference signal incompatible." << std::endl;
		std::cerr << std::endl;
		exit(EXIT_FAILURE);
	}
	if(para.intrfflag>0){
		if(para.rho_Wfb.n_rows==para.Ncomp_intf+1) ;
		else {
			std::cerr << "Size of rho_Wfb should be "<< para.Ncomp_intf+1 <<'x'<<para.Ncomp_intf+1<<" but is " << arma::size(para.rho_Wfb) <<std::endl;
			std::cerr << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	// resolve various flags interdependencies
	if(para.intrfflag==0) para.Ncomp_intf=0;
	if(para.LRflag>=2)	para.loop_learning_flag=0;
	else 				para.loop_learning_flag=1;

	para.multiflag=0;

//------------------------------------------------------------------------
// 	write flags to file
//------------------------------------------------------------------------
	std::ofstream flags_file(para.path+"flags.dat");	if (!flags_file) exit(EXIT_FAILURE);
	nlohmann::json json_flags;
	json_flags["actfuncflag"] = para.actfuncflag;
	json_flags["actfuncflag_GC"] = para.actfuncflag_GC;
	json_flags["cortexonlyflag"] = para.cortexonlyflag;
	json_flags["fcodingflag"] = para.fcodingflag;
	json_flags["filterflag"] = para.filterflag;
	json_flags["floquetflag"] = para.floquetflag;
	json_flags["floquetmodeflag"] = para.floquetmodeflag;
	json_flags["FMdyncleanflag"] = para.FMdyncleanflag;
	json_flags["gainflag"] = para.gainflag;
	json_flags["interferenceflag"] = para.intrfflag;
	json_flags["loopflag"] = para.loop_learning_flag;
	json_flags["loopflag_test"] = para.loop_test_flag;
	json_flags["learnflag"] = para.learnflag;
	json_flags["learnflag2"] = para.learnflag2;
	json_flags["LRflag"] = para.LRflag;
	json_flags["multiflag"] = para.multiflag;
	json_flags["noiseflag"] = para.noiseflag;
	json_flags["signalflag"] = para.signalflag;
	json_flags["simplemapflag"] = para.simplemapflag;
	json_flags["taudistrbflag"] = para.taudistrbflag;
	flags_file << std::setw(4) << json_flags << '\n';
	flags_file.close();

//------------------------------------------------------------------------
// 	generate some output
//------------------------------------------------------------------------
	std::cout << std::endl;
	arma::arma_version ver;
	std::cout << "Version: "<< para.version << std::endl;
	std::cout << "ARMA version: "<< ver.as_string() << std::endl;
	std::cout <<"----------------------------------------------------------------" << std::endl;
	std::cout << std::endl;

	std::cout << "Ncx = " << para.Ncx << "\t g = " << para.g << "\t pcon = " << para.pcon <<std::endl;
	std::cout <<"----------------------------------------------------------------" << std::endl;

	std::cout << "noise scaling: " << para.noisescaling << " with sig_N = " << para.sigmaN <<std::endl;
	std::cout << std::endl;

	if(para.taudistrbflag==1){
		std::cout << "Drawing RNN time-constants from Gamma distribution with avrg " << para.tau_avrg << " and std " << para.tau_std<< std::endl;
	}

	if(para.actfuncflag==0){
		std::cout << "Using tanh activation function for RNN units." << std::endl;
		para.Ptr_ctxIO = &phi_tanh;
	}
	else if(para.actfuncflag==1){
		std::cout << "Using linear activation function for RNN units."<< std::endl;
		para.Ptr_ctxIO = &phi_linear;
	}

	std::cout <<"----------------------------------------------------------------" << std::endl;
	std::cout << "target signal: " << para.signaltype <<  " with Ncomp = " << para.Ncomp << ", lambda = "<< para.lambda <<  ", nf = "<< para.nf << " and "<< std::endl;
	std::cout <<"Amp = " << arma::trans(para.amps.rows(0, para.Ncomp-1));
	std::cout <<"Phi = " << arma::trans(para.phases.rows(0, para.Ncomp-1)) << std::endl;
	if(para.intrfflag>0){
		if(para.intrfflag==1){
			std::cout << "adding interference signal to RNN with Ncomp_interf = " << para.Ncomp_intf << ", sigAmp_interf = " << para.sigAmp_intf << ", rho = " << std::endl;
		}
		std::cout<< para.rho_Wfb << std::endl;
		std::cout <<"Amp_interf = " << arma::trans(para.amps_intf.rows(0, para.Ncomp_intf-1));
		std::cout <<"Phi_interf = " << arma::trans(para.phases_intf.rows(0, para.Ncomp_intf-1));
	}
	std::cout <<"----------------------------------------------------------------" << std::endl;

	if(para.LRflag==2 || para.LRflag==3) std::cout << "Learning rule: " << para.learningrule << " with sigmaLS = " << para.sigmaLS << " with test-phase in " << para.loop_test << " loop." << std::endl;
	if(para.LRflag==4) std::cout << "Learning rule: " << para.learningrule << std::endl;
	if(para.LRflag==0) std::cout << "Learning rule: " << para.learningrule << " with alpha = " << para.alpha  << std::endl;
	if(para.LRflag==1) std::cout << "Learning rule: " << para.learningrule << " with eta_ini = " << para.delta_eta  << std::endl;
	std::cout << "Learning phase: " << para.duration_learn << " seconds." << std::endl;
	std::cout << "Test phase: " << para.duration_test << " seconds." << std::endl;

	std::cout << "Trajectory sampling time step dt = " << para.DT  << std::endl;
	if(para.floquetflag==1) {
		std::cout << "Calculating Floquet Multipliers in " << para.floquetmode << "-loop mode" << std::endl;
	}

	std::cout <<"----------------------------------------------------------------" << std::endl;
	std::cout << "Learning mode: " << para.learnmode << ", with Nreal = " << para.Nreal << " and seed = " << para.seed << std::endl;
	std::cout <<"----------------------------------------------------------------" << std::endl;

	std::cout << "Writing to " << para.path << " and " << para.path_connect <<std::endl;
	if(para.ROwriteflag==1) std::cout << "Writing readout activity to files." <<std::endl;

	return 0;
}
