int set_parameters_flags(netpara &para, int &argc, char* argv[]) {

//------------------------------------------------------------------------
// 	default parameters
//------------------------------------------------------------------------
  	para.DT = 0.001;			// integration time constant
  	para.DT_euler = 0.0001;		// integration time constant for euler test-phase
	para.Ncx=500;				// number of RNN units
	para.Ngc=500;				// number of GCs
	para.Nmf=para.Ncx;			// number of MFs
	para.Nsyn=4;				// number of MF synpases onto each GC
	para.lambda=5;				// signal fundamental frequency
	para.Ncomp=1;				// number of components of combined sine/cosine signal
	para.Ncomp_intf=2;			// number of components of interference signal
	para.g=0.9;
	para.seed=-1;				// default random number seeds
	para.Nreal=10;				// number of network realisations
	para.Npara=10;				// number of parameter steps
	para.N_FMperiods=10;		// number of periods for estimation of Floquet multipliers
	para.duration_learn=10;		// length of learning phase (in periods)
	para.duration_test=20;		// length of test phase (in periods)
	para.sigmaLS=pow(10,-2);	// regularisation parameter for LS or RR learning
	para.sigmaN=5*pow(10,-3);	// noise std added to cortical units
	para.sigmaNgc=para.sigmaN;	// noise std added to GCs
	para.Wgoc=0.75;				// GoC to GC weights
	para.Wgoc_SD=0.0;			// GoC to GC weights standard deviation
	para.rho_Wfb12=1;			//
	para.rho_Wfb13=1;			//
	para.write_every=20;		// read out every nth time step

	para.tau_avrg=0.01;			// average and standard deviation of RNN units integration time constants
	para.tau_std=0.01;

	para.delta_eta=0.0025;		// initial learning rate for DELTA rule
	para.delta_r2_reg=1e-6;		// norm regularisation for DELTA rule

	para.INT_delta_abs=1e-3;	// integrator absolute error
	para.INT_delta_rel=1e-3;	// integrator relative error

	para.GCgain_base=1;
	para.GCgainmax=100;			// maximum allowed GC gain
	para.GCmax=1;				// maximum GC firing rate

	para.sp_alpha=20;			// sharpness parameter of GC softplus function
	para.fc_avrg=0.5;			// GC coding level avrg and std
	para.fc_SD=pow(10,-5);

	para.sigAmp=1.0;			// global amplitude multiplier
	para.sigAmp_intf=1.0;		// global amplitude multiplier (interference signal)

	para.amps={1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};	// all amplitudes of signal / interference signals
	para.amps_intf={0.0,1.0,1.0,1.0,1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

	para.phases={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};		// all phases of signal / interference signals
	para.phases_intf={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

	// target signal: 'quadruple_sine', 'sine', 'cosine_fixed', 'test_signal' or 'zero_signal'
	para.signaltype="cosine_fixed";
	// scale noise with firing rate: 'constant', 'proportional', 'rate' or 'none'
	para.noisescaling="none";
	// learn mode: 'cortex', 'cerebellum', 'none', 'repeat' or 'scan_...'
	para.learnmode="cerebellum";
	// learning rule: 'RLS', 'LS', 'RR', 'PINV' or 'DELTA'
	para.learningrule="LS";
	// testphase in 'closed' or 'open' loop or with open loop with noise exclusively: 'open noise'
	para.loop_test="closed";
	// activation function of reservoir units: 'tanh', 'linear' or 'relu'
	para.activation_function="tanh";
	// activation function of GCs: 'relu', tanh, 'softplus' or 'linear'
	para.activation_function_GC="softplus";
	// MF to GC connectivity: 'mixing' or 'simple'
	para.MF_GC_connect="mixing";
	// adjust GC coding level instead of setting threshold: 'global', 'local', 'OFF'
	para.fcodingmode="local";
	// adjust GC gain: 'global', 'local', 'fixed', 'inverse_norm', 'matched_vectors'
	para.gainmode="fixed";
	// GoC mode: 'ON', 'OFF' or 'MFSD'
	para.GoC_mode="MFSD";
	// add interference signal: 'ON', 'OFF' or 'GCs'
	para.interference_signal="OFF";
	// cortical tau distribution: 'ON' or 'OFF'
	para.tau_mode="OFF";
	// gc tau: 'ON', 'OFF' or HETERO
	para.taugc_mode="OFF";
	// calculate Floquet multipiers: 'ON' or 'OFF'
	para.floquet="OFF";
	// calculate Floquet multipiers: 'closed' or 'open'
	para.floquetmode="open";
	// write readout unit output to file
	para.ROwrite="OFF";

	try {

		po::options_description general("General options", 240);
		general.add_options()
			("help", "Produce help message")
			("seed", po::value<int>(&para.seed)->default_value(para.seed), "random number seed")
			("dt", po::value<double>(&para.DT)->default_value(para.DT), "integration time step")
			("dt_euler", po::value<double>(&para.DT_euler)->default_value(para.DT_euler), "integration time step of euler test-phase")
			("Nreal", po::value<int>(&para.Nreal)->default_value(para.Nreal), "number of network realisations")
			("noise", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.noisescaling), "noise (scaling with firing rate): 'constant', 'proportional' or 'none' ")
			("Floquet", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.floquet), "calculate Floquet multipiers: 'ON' or 'OFF' ")
			("loop_Floquet_mode", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.floquetmode), "Floquet mode: 'closed' or 'open' ")
			("write_every", po::value<int>(&para.write_every)->default_value(para.write_every), "read-out sampling rate (every xth time-step)")
			("ROwrite", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.ROwrite), "write readout to file: 'ON' or 'OFF' ")
			("rho_Wfb12", po::value<double>(&para.rho_Wfb12)->default_value(para.rho_Wfb12), "correlation rho_Wfb12")
			("rho_Wfb13", po::value<double>(&para.rho_Wfb13)->default_value(para.rho_Wfb13), "correlation rho_Wfb13")
			("Npara", po::value<int>(&para.Npara)->default_value(para.Npara), "number of parameter steps")
			("N_FMperiods", po::value<int>(&para.N_FMperiods)->default_value(para.N_FMperiods), "number of periods for calculation of Floquet multiplier")
			("INT_delta_abs", po::value<double>(&para.INT_delta_abs)->default_value(para.INT_delta_abs), "target absolute error for numerical integration")
			("INT_delta_rel", po::value<double>(&para.INT_delta_rel)->default_value(para.INT_delta_rel), "target relative error for numerical integration")
		;

		po::options_description Input("Input/target signal(s) options", 240);
		Input.add_options()
			("signaltype", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.signaltype), "target signal: 'quadruple_sine', 'cosine_fixed'... ")
			("interference_signal", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.interference_signal), "add interference signal: 'ON', 'OFF' or 'GCs' ")
			("lambda", po::value<double>(&para.lambda)->default_value(para.lambda), "signal fundamental frequency")
			("Ncomp", po::value<int>(&para.Ncomp)->default_value(para.Ncomp), "number of signal frequency components")
			("Ncomp_intf", po::value<int>(&para.Ncomp_intf)->default_value(para.Ncomp_intf), "number of interference signal frequency components")
			("sigAmp", po::value<double>(&para.sigAmp)->default_value(para.sigAmp), "global target signal multiplier")
			("sigAmp_intf", po::value<double>(&para.sigAmp_intf)->default_value(para.sigAmp_intf), "global interference signal multiplier")
		;

		po::options_description RNN("RNN options", 240);
		RNN.add_options()
			("Ncx", po::value<int>(&para.Ncx)->default_value(para.Ncx), "number of RNN units")
			("g", po::value<double>(&para.g)->default_value(para.g), "reccurent weight parameter")
			("tau_mode", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.tau_mode), "distribute RNN taus: 'ON' or 'OFF' ")
			("tau_avrg", po::value<double>(&para.tau_avrg)->default_value(para.tau_avrg), "average RNN unit time constant")
			("tau_std", po::value<double>(&para.tau_std)->default_value(para.tau_std), "std of RNN unit time constants")
			("activation_function", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.activation_function), "activation function of RNN units: 'tanh' or 'linear' ")
			("sigmaN", po::value<double>(&para.sigmaN)->default_value(para.sigmaN), "RNN unit noise")
		;

		po::options_description FFL("FFL options", 240);
		FFL.add_options()
			("Ngc", po::value<int>(&para.Ngc)->default_value(para.Ngc), "number of GCs")
			("Nmf", po::value<int>(&para.Nmf)->default_value(para.Nmf), "number of Mfs")
			("sigmaNgc", po::value<double>(&para.sigmaNgc)->default_value(para.sigmaNgc), "GC noise")
			("activation_function_GC", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.activation_function_GC), "activation function of GCs: 'relu', 'tanh', 'softplus' or 'linear' ")
			("MF_GC_connect", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.MF_GC_connect), "MF to GC connectivity: 'mixing' or 'simple' ")
			("taugc_mode", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.taugc_mode), "GC taus: 'ON', 'OFF' or 'HETERO' ")
			("fcodingmode", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.fcodingmode), "GC coding level (global, local or fix threshold): 'global', 'local', 'fixed' ")
			("gainmode", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.gainmode), "GC gain: 'global', 'local', 'fixed', 'inverse_norm', 'matched_vectors', 'local_matched' ")
			("GoC_mode", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.GoC_mode), "Golgi cell: 'ON', 'OFF' or 'MFSD' ")
			("Wgoc", po::value<double>(&para.Wgoc)->default_value(para.Wgoc), "GoC to GC weight average")
			("Wgoc_SD", po::value<double>(&para.Wgoc_SD)->default_value(para.Wgoc_SD), "GoC to GC weight standard deviation")
			("sp_alpha", po::value<double>(&para.sp_alpha)->default_value(para.sp_alpha), "sharpness of GC softplus function")
			("fc_avrg", po::value<double>(&para.fc_avrg)->default_value(para.fc_avrg), "GC coding level average")
			("fc_SD", po::value<double>(&para.fc_SD)->default_value(para.fc_SD), "GC coding level SD")
			("GCgain", po::value<double>(&para.GCgain_base)->default_value(para.GCgain_base), "GC gain")
			("GCgainmax", po::value<double>(&para.GCgainmax)->default_value(para.GCgainmax), "maximum GC gain")
		;

		po::options_description Learning("Learning and testing options", 240);
		Learning.add_options()
			("learnmode", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.learnmode), "learn mode: 'cortex', 'cerebellum', 'none', 'repeat' or 'scan_...' ")
			("LRule", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.learningrule), "learning rule: 'RLS', 'LS', 'RR', 'PINV' or 'DELTA' ")
			("loop_test", po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), para.loop_test), "loop type during test-phase: 'closed' or 'open'")
			("T_learn", po::value<double>(&para.duration_learn)->default_value(para.duration_learn), "length of learning phase (i.e. number of periods)")
			("T_test", po::value<double>(&para.duration_test)->default_value(para.duration_test), "length of test phase  (i.e. number of periods)")
			("sigmaLS", po::value<double>(&para.sigmaLS)->default_value(para.sigmaLS), "regularisation for LS or RR")
			("delta_eta", po::value<double>(&para.delta_eta)->default_value(para.delta_eta), "initial learning rate for DELTA rule")
		;

		// Combine all descriptions
		po::options_description alldesc("Allowed options", 240);
		alldesc.add(general).add(Input).add(Learning).add(RNN).add(FFL);

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, alldesc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			std::cout << alldesc << "\n";
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
		if (vm["activation_function_GC"].defaulted()==0)  {
			tmp=vm["activation_function_GC"].as< std::vector<std::string> >();
			para.activation_function_GC= tmp[0];
		}
		if (vm["MF_GC_connect"].defaulted()==0)  {
			tmp=vm["MF_GC_connect"].as< std::vector<std::string> >();
			para.MF_GC_connect= tmp[0];
		}
		if (vm["fcodingmode"].defaulted()==0)  {
			tmp=vm["fcodingmode"].as< std::vector<std::string> >();
			para.fcodingmode= tmp[0];
		}
		if (vm["gainmode"].defaulted()==0)  {
			tmp=vm["gainmode"].as< std::vector<std::string> >();
			para.gainmode= tmp[0];
		}
		if (vm["GoC_mode"].defaulted()==0)  {
			tmp=vm["GoC_mode"].as< std::vector<std::string> >();
			para.GoC_mode= tmp[0];
		}
		if (vm["interference_signal"].defaulted()==0)  {
			tmp=vm["interference_signal"].as< std::vector<std::string> >();
			para.interference_signal= tmp[0];
		}
		if (vm["tau_mode"].defaulted()==0)  {
			tmp=vm["tau_mode"].as< std::vector<std::string> >();
			para.tau_mode= tmp[0];
		}
		if (vm["taugc_mode"].defaulted()==0)  {
			tmp=vm["taugc_mode"].as< std::vector<std::string> >();
			para.taugc_mode= tmp[0];
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
	}
	catch(std::exception& e) {
		std::cerr << "error: " << e.what() << std::endl;
		return 1;
	}
	catch(...) {
		std::cerr << "Exception of unknown type." << std::endl;
	}

	// target folders
	para.path="simulation_results/";
	para.path_connect="connectivity/";

//------------------------------------------------------------------------
// 	set fixed parameters
//------------------------------------------------------------------------

	para.RLSlearn_every=2;					// learn at every nth time step (RLS mode)

	para.duration_FT=3;						// time-window for FT in seconds
	para.duration_trans=15;					// time-window for transient activity in seconds

	para.Ncomp_noise=60;					// number of FT noise components to be extracted
	para.Ncomp_noise_RNN=60;				// number of FT noise components to be extracted

	para.FCcomp_tol=1e-6;					// cutoff for FT component entries

	para.Nwrite= 10; 						// number of components to be written to otput

	// correlation between feedback weights of target and interference signal
	para.rho_Wfb.ones(2*para.Ncomp+(para.Ncomp_intf-para.Ncomp),2*para.Ncomp+(para.Ncomp_intf-para.Ncomp));

	if(para.Ncomp==1){
		if(para.Ncomp_intf>1){
			para.rho_Wfb(0,2)=para.rho_Wfb12;
			para.rho_Wfb(2,0)=para.rho_Wfb12;
		}
		if(para.Ncomp_intf>2){
			para.rho_Wfb(0,3)=para.rho_Wfb13;
			para.rho_Wfb(3,0)=para.rho_Wfb13;

			if(para.interference_signal == "GCs"){
				para.rho_Wfb(1,2)=0.9;
				para.rho_Wfb(2,1)=0.9;
				para.rho_Wfb(1,3)=0.7;
				para.rho_Wfb(3,1)=0.7;

				para.rho_Wfb(2,3)=0.8;
				para.rho_Wfb(3,2)=0.8;
				para.amps_intf={0.0,1.0,1.0,1.0,1.0};
			}
		}
	}
	if(para.Ncomp==2){
		if(para.Ncomp_intf>3){
			para.rho_Wfb(0,4)=para.rho_Wfb12;
			para.rho_Wfb(4,0)=para.rho_Wfb12;
			para.rho_Wfb(0,5)=para.rho_Wfb13;
			para.rho_Wfb(5,0)=para.rho_Wfb13;

			if(para.interference_signal == "GCs"){
				para.rho_Wfb(1,4)=0.9;
				para.rho_Wfb(4,1)=0.9;
				para.rho_Wfb(1,5)=0.9;
				para.rho_Wfb(5,1)=0.9;
				para.rho_Wfb(2,5)=0.8;
				para.rho_Wfb(5,2)=0.8;
				para.rho_Wfb(2,4)=0.7;
				para.rho_Wfb(4,2)=0.7;
				para.rho_Wfb(3,5)=0.8;
				para.rho_Wfb(5,3)=0.8;
				para.rho_Wfb(3,4)=0.7;
				para.rho_Wfb(4,3)=0.7;

				para.rho_Wfb(5,4)=0.6;
				para.rho_Wfb(4,5)=0.6;

				para.amps_intf={0.0,0.0,1.0,1.0,1.0};
			}
		}
	}

	para.tauRNNcutoff=0.003;		// smallest permitted tau

									// average and standard deviation of GC integration time constant
	para.tauGC_avrg=0.00339159;		// for gamma distribution fitted to data from Schwartz, Rothman et al. 2012
	para.tauGC_std=0.00122581;
	para.tauGCcutoff=0.0014389;		// smallest tauGC encountered in Schwartz, Rothman et al. 2012

	para.pcon=0.1;				// connection probability between cortical units

	para.GCthreshold_base=0;
	para.GCgain=para.GCgain_base;

	para.phi1=0.0;			// phase of 1st cosine/sine component  	//~ para.phi1=3*M_PI/4.0;
	para.ampmod=0;			// modification of 1st & 2nd cosine components

	para.alpha=1;			// regularisation term of P matrix (for RLS)

	para.tol_GCgain=6;		// standard deviations above which GC gain is considered too large

	para.sqNsyn=sqrt(para.Nsyn);

//------------------------------------------------------------------------
// set all flags
//------------------------------------------------------------------------
	if 		(para.learnmode == "cortex") 		para.learnflag=1;
    else if (para.learnmode == "cerebellum") 	para.learnflag=2;
    else if (para.learnmode == "none") 			para.learnflag=0;
    else if (para.learnmode == "repeat") 		para.learnflag=99;
    else if (para.learnmode == "scan_sigmaN") 	para.learnflag=100;
    else if (para.learnmode == "scan_Nctx") 	para.learnflag=102;
    else if (para.learnmode == "scan_Ngc") 		para.learnflag=103;
    else if (para.learnmode == "scan_sigmaNgc") para.learnflag=105;
    else if (para.learnmode == "scan_Wgoc") 	para.learnflag=106;
    else if (para.learnmode == "scan_Ncomp") 	para.learnflag=107;
    else if (para.learnmode == "scan_sigAmp") 	para.learnflag=108;
    else if (para.learnmode == "scan_lambda") 	para.learnflag=109;
    else if (para.learnmode == "scan_Tlearn") 	para.learnflag=110;
    else if (para.learnmode == "scan_g") 		para.learnflag=113;
    else if (para.learnmode == "scan_eta") 		para.learnflag=114;
    else if (para.learnmode == "scan_Phi") 		para.learnflag=118;
    else if (para.learnmode == "scan_Ampmod") 	para.learnflag=119;
    else if (para.learnmode == "scan_sigmaLS") 	para.learnflag=120;
    else if (para.learnmode == "scan_fcoding") 	para.learnflag=123;
    else if (para.learnmode == "scan_sp_alpha") para.learnflag=124;
    else if (para.learnmode == "scan_rho_Wfb12") 	para.learnflag=128;
    else if (para.learnmode == "scan_sigAmp_intf") para.learnflag=129;
    else if (para.learnmode == "scan_gain") 	para.learnflag=131;
    else if (para.learnmode == "scan_rho_Wfb13") 	para.learnflag=132;
    else {std::cerr <<"Bad learnmode string. \n";exit(EXIT_FAILURE);}

    if 		(para.signaltype == "quadruple_sine") 	para.signalflag=0;
    else if (para.signaltype == "sine") 			para.signalflag=22;
    else if (para.signaltype == "cosine_fixed") 	para.signalflag=23;
    else if (para.signaltype == "cosine_sine") 		para.signalflag=24;
    else if (para.signaltype == "zero_signal") 		para.signalflag=-1;
    else {std::cerr <<"Bad signaltype string. \n";exit(EXIT_FAILURE);}

	if 		(para.noisescaling == "none") 			para.noiseflag=-1;
    else if (para.noisescaling == "constant") 		para.noiseflag=0;
    else if (para.noisescaling == "proportional") 	para.noiseflag=1;
    else {std::cerr <<"Bad noisescaling string. \n";exit(EXIT_FAILURE);}

	if 		(para.learningrule == "RLS") 	para.LRflag=0;
    else if (para.learningrule == "DELTA") 	para.LRflag=1;
    else if (para.learningrule == "LS") 	para.LRflag=2;
    else if (para.learningrule == "RR") 	para.LRflag=3;
    else if (para.learningrule == "PINV") 	para.LRflag=4;
    else {std::cerr <<"Bad LRstring string. \n";exit(EXIT_FAILURE);}

    if 		(para.fcodingmode == "OFF") 		para.fcodingflag=0;
    else if (para.fcodingmode == "global") 	para.fcodingflag=1;
    else if (para.fcodingmode == "local") 	para.fcodingflag=2;
	else {std::cerr <<"Bad fcodingmode string. \n";exit(EXIT_FAILURE);}

	if 		(para.gainmode == "global_theo") 	para.gainflag=0;
    else if (para.gainmode == "global") 		para.gainflag=1;
    else if (para.gainmode == "local") 			para.gainflag=2;
    else if (para.gainmode == "fixed") 			para.gainflag=3;
    else if (para.gainmode == "inverse_norm") 	para.gainflag=4;
    else if (para.gainmode == "matched_vectors") para.gainflag=5;
    else if (para.gainmode == "local_matched") 	para.gainflag=6;
    else if (para.gainmode == "maxrate") 		para.gainflag=7;
    else if (para.gainmode == "maxrate_matched") para.gainflag=8;
	else {std::cerr <<"Bad gainmode string. \n";exit(EXIT_FAILURE);}

    if 		(para.interference_signal == "OFF") para.intrfflag=0;
    else if (para.interference_signal == "ON") 	para.intrfflag=1;
    else if (para.interference_signal == "GCs") para.intrfflag=2;
    else {std::cerr <<"Bad interference_signal string. \n";exit(EXIT_FAILURE);}

    if (para.activation_function == "tanh") 		para.actfuncflag=0;
    else if (para.activation_function == "linear") 	para.actfuncflag=1;
	else {std::cerr <<"Bad activation_function string. \n";exit(EXIT_FAILURE);}

    if (para.activation_function_GC == "relu") 			para.actfuncflag_GC=0;
    else if (para.activation_function_GC == "softplus") para.actfuncflag_GC=1;
    else if (para.activation_function_GC == "linear") 	para.actfuncflag_GC=2;
    else if (para.activation_function_GC == "tanh") 	para.actfuncflag_GC=3;
	else {std::cerr <<"Bad activation_function_GC string. \n";exit(EXIT_FAILURE);}

    if (para.MF_GC_connect == "mixing") 		para.simplemapflag=0;
    else if (para.MF_GC_connect == "simple") 	para.simplemapflag=1;
	else {std::cerr <<"Bad MF_GC_connect string. \n";exit(EXIT_FAILURE);}

	if (para.GoC_mode== "OFF") 						para.gocmodeflag=0;
    else if (para.GoC_mode == "ON") 				para.gocmodeflag=1;
    else if (para.GoC_mode == "MFSD") 				para.gocmodeflag=3;
	else {std::cerr <<"Bad GoC_mode string. \n";exit(EXIT_FAILURE);}

	if (para.tau_mode== "OFF") 					para.taudistrbflag=0;
    else if (para.tau_mode == "ON") 				para.taudistrbflag=1;
	else {std::cerr <<"Bad tau_mode string. \n";exit(EXIT_FAILURE);}

	if (para.taugc_mode== "OFF") 					para.taugcmodeflag=0;
    else if (para.taugc_mode == "ON") 				para.taugcmodeflag=1;
    else if (para.taugc_mode == "HETERO") 			para.taugcmodeflag=2;
	else {std::cerr <<"Bad taugc_mode string. \n";exit(EXIT_FAILURE);}

	if (para.floquet == "OFF") 						para.floquetflag=0;
    else if (para.floquet == "ON") 					para.floquetflag=1;
	else {std::cerr <<"Bad floquet string. \n";exit(EXIT_FAILURE);}

	if (para.floquetmode == "open") 				para.floquetmodeflag=0;
    else if (para.floquetmode == "closed") 			para.floquetmodeflag=1;
	else {std::cerr <<"Bad floquetmode string. \n";exit(EXIT_FAILURE);}

	if (para.ROwrite == "OFF") 						para.ROwriteflag=0;
    else if (para.ROwrite == "ON") 					para.ROwriteflag=1;
	else {std::cerr <<"Bad ROwrite string. \n";exit(EXIT_FAILURE);}

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
	if(para.learnflag==103){
		para.simplemapflag=0;
	}
	if(para.learnflag==117){
		para.simplemapflag=0;
		if(para.fcodingflag==0){
			std::cerr << "Learn- and fcoding-modes incompatible." << std::endl;
			std::cerr << std::endl;
			exit(EXIT_FAILURE);
		}
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
	if(para.learnflag==123){
		if(para.fcodingflag==0){
			std::cerr << "Learn- and fcoding-modes incompatible." << std::endl;
			std::cerr << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	if(para.learnflag==124){
		para.actfuncflag_GC=1;
	}
	if(para.gocmodeflag==0){
		para.Wgoc=0;
		para.Wgoc_SD=0;
	}

	if( (para.learnflag==128 || para.learnflag==129) && para.intrfflag==0 ){
		std::cerr << "Learn- and interference-modes incompatible." << std::endl;
		std::cerr << std::endl;
		exit(EXIT_FAILURE);
	}
	if(para.learnflag==131){
		para.gainflag=3;
		para.actfuncflag_GC=2;
		para.gainmode = "fixed";
		para.activation_function_GC = "linear";
	}

	// check for incompatible threshold & gain modes
	if(para.fcodingflag>0 && para.gainflag==0){
		std::cerr << "Threshold- and gain-modes incompatible." << std::endl;
		std::cerr << std::endl;
		exit(EXIT_FAILURE);
	}

	// check for incompatible interference signal parameters
	if(para.intrfflag>0 && !(para.signalflag==22 || para.signalflag==23 || para.signalflag==24 || para.signalflag==25) ){
		std::cerr << "Target signal and interference signal incompatible." << std::endl;
		std::cerr << std::endl;
		exit(EXIT_FAILURE);
	}
	if(para.intrfflag>0){
		//~ if(para.rho_Wfb.n_rows==para.Ncomp_intf+1) ;
		if(para.rho_Wfb.n_rows== (2*para.Ncomp+(para.Ncomp_intf-para.Ncomp))) ;
		else {
			std::cerr << "Size of rho_Wfb should be "<< para.Ncomp_intf+1 <<'x'<<para.Ncomp_intf+1<<" but is " << arma::size(para.rho_Wfb) <<std::endl;
			std::cerr << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	// resolve various flags interdependencies
	if(para.intrfflag==0) para.Ncomp_intf=0;
	else if(para.intrfflag==2) {
		para.actfuncflag_GC=2;
		para.activation_function_GC = "linear";
	}
    if(para.simplemapflag==1)	{
		para.Ngc=para.Nmf=para.Ncx;
		para.Nsyn=1;
		para.sqNsyn=sqrt(para.Nsyn);
	}
	if(para.Nmf>para.Ncx)	{
		para.Nmf=para.Ncx;
	}
	if(para.LRflag>=2)	para.loop_learning_flag=0;		// flag is only for output
	else 				para.loop_learning_flag=1;

//------------------------------------------------------------------------
// 	write flags to file
//------------------------------------------------------------------------
	std::ofstream flags_file(para.path+"flags.dat");	if (!flags_file) exit(EXIT_FAILURE);
	nlohmann::json json_flags;
	json_flags["actfuncflag"] = para.actfuncflag;
	json_flags["actfuncflag_GC"] = para.actfuncflag_GC;
	json_flags["cortexonlyflag"] = para.cortexonlyflag;
	json_flags["fcodingflag"] = para.fcodingflag;
	json_flags["floquetflag"] = para.floquetflag;
	json_flags["floquetmodeflag"] = para.floquetmodeflag;
	json_flags["gainflag"] = para.gainflag;
	json_flags["gocmodeflag"] = para.gocmodeflag;
	json_flags["interferenceflag"] = para.intrfflag;
	json_flags["loopflag"] = para.loop_learning_flag;
	json_flags["loopflag_test"] = para.loop_test_flag;
	json_flags["learnflag"] = para.learnflag;
	json_flags["LRflag"] = para.LRflag;
	json_flags["noiseflag"] = para.noiseflag;
	json_flags["ROwriteflag"] = para.ROwriteflag;
	json_flags["signalflag"] = para.signalflag;
	json_flags["simplemapflag"] = para.simplemapflag;
	json_flags["taudistrbflag"] = para.taudistrbflag;
	json_flags["taugcmodeflag"] = para.taugcmodeflag;
	flags_file << std::setw(4) << json_flags << '\n';
	flags_file.close();

//------------------------------------------------------------------------
// 	generate some output
//------------------------------------------------------------------------
	std::cout << std::endl;
	arma::arma_version ver;
	std::cout << "Version: "<< para.version << std::endl;
	std::cout << "ARMA version: "<< ver.as_string() << std::endl;
	std::cout << "Using Boost version " << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "."<< BOOST_VERSION % 100<< std::endl;

	std::cout <<"----------------------------------------------------------------" << std::endl;
	std::cout << std::endl;

	std::cout << "Ncx = " << para.Ncx << "\t Nmf = " << para.Nmf << "\t Ngc = " << para.Ngc << "\t Nsyn = " << para.Nsyn<<  std::endl;
	std::cout << "g = " << para.g << "\t pcon = " << para.pcon <<std::endl;
	std::cout <<"----------------------------------------------------------------" << std::endl;

	std::cout << "noise scaling: " << para.noisescaling << " with sig_N = " << para.sigmaN << " and sig_N_gc = " << para.sigmaNgc <<std::endl;
	std::cout << std::endl;

	if(para.taudistrbflag==0){
		std::cout << "RNN time-constants set to " << para.tau_avrg << std::endl;
	}
	else if(para.taudistrbflag==1){
		std::cout << "Drawing RNN time-constants from Gamma distribution with avrg " << para.tau_avrg << " and std " << para.tau_std<< std::endl;
	}
	if(para.taugcmodeflag==0){
		std::cout << "Instantaneous GCs"<< std::endl;
	}
	else if(para.taugcmodeflag==1){
		std::cout << "Finite GC time-constants with tauGC = " << para.tauGC_avrg << std::endl;
	}
	else if(para.taugcmodeflag==2){
		std::cout << "Drawing GC time-constants from Gamma distribution with avrg " << para.tauGC_avrg << " and std " << para.tauGC_std<< std::endl;
	}
	std::cout << std::endl;

	if(para.actfuncflag==0){
		std::cout << "Using tanh activation function for RNN units." << std::endl;
		para.Ptr_ctxIO = &phi_tanh;
	}
	else if(para.actfuncflag==1){
		std::cout << "Using linear activation function for RNN units."<< std::endl;
		para.Ptr_ctxIO = &phi_linear;
	}

	if(para.actfuncflag_GC==0){
		std::cout << "Using relu activation function for GCs." << std::endl;
		para.Ptr_gcIO = &phi_relu;
	}
	else if(para.actfuncflag_GC==1){
		std::cout << "Using softplus activation function for GCs with alpha = " << para.sp_alpha << std::endl;
		para.Ptr_gcIO = &phi_softplus;
	}
	else if(para.actfuncflag_GC==2){
		std::cout << "Using linear activation function for GCs."<< std::endl;
		para.Ptr_gcIO = &phi_linear;
	}
	else if(para.actfuncflag_GC==3){
		std::cout << "Using tanh activation function for GCs."<< std::endl;
		para.Ptr_gcIO = &phi_tanh;
	}
	if(para.simplemapflag==1){
		std::cout << std::endl;
		std::cout << "Using one-to-one mapping of ctx activity to GC activity." << std::endl;
	}

	std::cout << std::endl;

	if(para.fcodingflag==1)			std::cout << "Adjusting GC coding levels globally with fc_avrg = " << para.fc_avrg <<std::endl;
	else if(para.fcodingflag==2)	std::cout << "Adjusting GC coding levels individually with fc_avrg = " << para.fc_avrg << " and fc_SD = " <<  para.fc_SD <<std::endl;
	else 							std::cout << "Adjusting thresholds by setting common value for theta." <<std::endl;

	if(para.gainflag==0)		std::cout << "Gains set globally via theory." <<std::endl;
	else if (para.gainflag==1) 	std::cout << "Gains set globally based on population SD with GCgainmax = " << para.GCgainmax <<std::endl;
	else if (para.gainflag==2)	std::cout << "Gains set locally based on individual SD with GCgainmax = " << para.GCgainmax<<std::endl;
	else if (para.gainflag==3)	std::cout << "Gains fixed globally to: " << para.GCgain_base << std::endl;
	else if (para.gainflag==4)	std::cout << "Gains set to inverse norm of first spanning vector." << std::endl;
	else if (para.gainflag==5)	std::cout << "Gains set so as to match norms of first spanning vectors in RNN and GCs with GCgainmax = " << para.GCgainmax << std::endl;
	else if (para.gainflag==6)	std::cout << "Gains set 1) locally and 2) so as to match norms of first spanning vectors in RNN and GCs with GCgainmax = " << para.GCgainmax << std::endl;

	if(para.gocmodeflag==1)			std::cout << "Adding GoC feedback with Wgoc = " << para.Wgoc <<std::endl;
	else if(para.gocmodeflag==3)	std::cout << "Adding FF GoC inhibition based on MF SD with Wgoc = " << para.Wgoc <<std::endl;

	std::cout <<"----------------------------------------------------------------" << std::endl;
	std::cout << "target signal: " << para.signaltype <<  " with Ncomp = " << para.Ncomp << ", lambda = "<< para.lambda << std::endl;
	//~ std::cout <<"Amp = " << arma::trans(para.amps.rows(0, para.Ncomp-1));
	//~ std::cout <<"Phi = " << arma::trans(para.phases.rows(0, para.Ncomp-1)) << std::endl;
	if(para.intrfflag>0){
		if(para.intrfflag==1){
			std::cout << "adding interference signal to RNN with Ncomp_interf = " << para.Ncomp_intf << ", sigAmp_interf = " << para.sigAmp_intf << ", rho = " << std::endl;
		}
		else{
			std::cout << "adding interference signal to GCs with Ncomp_interf = " << para.Ncomp_intf << ", sigAmp_interf = " << para.sigAmp_intf << ", rho = " << std::endl;
		}
		//~ std::cout<< para.rho_Wfb << std::endl;
		//~ std::cout <<"Amp_interf = " << arma::trans(para.amps_intf.rows(0, para.Ncomp_intf-1));
		//~ std::cout <<"Phi_interf = " << arma::trans(para.phases_intf.rows(0, para.Ncomp_intf-1));
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
