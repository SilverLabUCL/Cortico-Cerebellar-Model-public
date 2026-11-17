///////////////////////////////////
// target signal function
///////////////////////////////////
double func_target(double t, netpara &para){
	double f=0;
	for(int j= 0; j < para.Ncomp; j++) {
		f += para.amps_func[j]*cos(para.freqs_func[j]*t + para.phases_func[j]);
	}
	return  f;
}
///////////////////////////////////
// interference signal function
///////////////////////////////////
double func_target_intf(double t, int idx_intf, netpara &para){
	return para.amps_intf_func[idx_intf]*cos(para.freqs_intf_func[idx_intf]*t + para.phases_intf_func[idx_intf]);
}

void target_signal(netpara &para, nlohmann::json &json_parameters,  std::ofstream &log_file) {

	int Tmax_learn=para.Tmax_learn;
	int Tmax_test=para.Tmax_test;
	int Tmax_FT=para.Tmax_FT;
	int Tmax_trans=para.Tmax_trans;
	int Tmax_FM=para.Tmax_FM;

	int Ncomp=para.Ncomp;
	int Ncomp_intf=para.Ncomp_intf;

	para.signalavrg=0;
	para.signalvar=0;

	para.freqs=arma::zeros<arma::vec>(para.Ncomp);
	para.freqs_intf=arma::zeros<arma::vec>(para.Ncomp_intf);

	para.freqs_func=arma::zeros<arma::vec>(para.Ncomp);
	para.freqs_intf_func=arma::zeros<arma::vec>(para.Ncomp_intf);
	para.amps_func=arma::zeros<arma::vec>(para.Ncomp);
	para.amps_intf_func=arma::zeros<arma::vec>(para.Ncomp_intf);
	para.phases_func=arma::zeros<arma::vec>(para.Ncomp);
	para.phases_intf_func=arma::zeros<arma::vec>(para.Ncomp_intf);

	if(para.signalflag==23)	{
		log_file << "\t generate cosine target signal(s)"<< std::endl;
	}
	else{
		log_file << "\t generate target signal(s)"<< std::endl;
	}

	// Mix of 4 sine waves from Sussillo & Abbott 2009, Fig 2d; the original signal is obtained for lambda=0.8;
	if (para.signalflag==0){
		double amp = para.sigAmp*1.6/1.5;

		para.freqs(0) =para.lambda;
		para.freqs(1) =2*para.lambda;
		para.freqs(2) =3*para.lambda;
		para.freqs(3) =4*para.lambda;

		para.amps_func(0) =amp;
		para.amps_func(1) =amp/2.0;
		para.amps_func(2) =amp/6.0;
		para.amps_func(3) =amp/3.0;

		for(int j= 0; j < Ncomp; j++) {
			para.freqs_func(j) =2.0*M_PI*para.freqs(j);
			para.phases_func(j) =M_PI/2.0;
		}

		// write output to file
		nlohmann::json json_freq_tmp;
		nlohmann::json json_amp_tmp;
		nlohmann::json json_phase_tmp;
		for(int j = 0; j < Ncomp; j++) {
			json_freq_tmp.push_back(para.freqs(j));
			json_amp_tmp.push_back(para.amps_func(j));
			json_phase_tmp.push_back(para.phases_func(j));
		}
		json_parameters["frequency"].push_back(json_freq_tmp);
		json_parameters["amp"].push_back(json_amp_tmp);
		json_parameters["phase"].push_back(json_phase_tmp);

	}
	// cos (or sin) waves with Ncomp components with fixed parameters
	else if(para.signalflag==23 || para.signalflag==22){
		arma::vec phase(Ncomp);
		arma::vec phase_intf(Ncomp_intf);
		arma::vec freq(Ncomp);
		arma::vec freq_intf(Ncomp_intf);
		arma::vec amp(Ncomp);
		arma::vec amp_intf(Ncomp_intf);
		double w;

		// phases
		for(int j= 0; j < Ncomp; j++) {
			phase(j)=para.phases(j);
			if(para.signalflag==22){
				phase(j)=para.phases(j) +M_PI/2.0;
			}
		}
		phase(0)=para.phi1;
		if(para.signalflag==22){
			phase(0)=para.phi1 +M_PI/2.0;
		}
		for(int j= 0; j < Ncomp_intf; j++) {
			phase_intf(j)=para.phases_intf(j);
			if(para.signalflag==22){
				phase_intf(j)=para.phases_intf(j) +M_PI/2.0;
			}
		}
		// frequencies
		for(int j= 0; j < Ncomp; j++) {
			freq(j)=double(j+1)*para.lambda;		// overtones = integer multiple of base frequency
		}
		// interference frequencies
		for(int j= 0; j < Ncomp_intf; j++) {
			freq_intf(j)=double(j+1)*para.lambda;
		}
		// norm of amplitudes
		int Ncomp_eff;
		if(Ncomp_intf>0)	Ncomp_eff=accu(para.amps.rows(0, Ncomp-1) != 0)+accu(para.amps_intf.rows(0, Ncomp_intf-1) != 0);
		else 				Ncomp_eff=accu(para.amps.rows(0, Ncomp-1) != 0);
		double ampnorm=sqrt(double(Ncomp_eff));

		// amplitudes
		for(int j= 0; j < Ncomp; j++) {
			amp(j)=para.sigAmp*para.amps(j)/ampnorm;
		}
		for(int j= 0; j < Ncomp_intf; j++) {
			amp_intf(j)=para.sigAmp_intf*para.amps_intf(j)/ampnorm;
		}

		for(int j= 0; j < Ncomp; j++) {
			para.freqs(j) =freq(j);
			para.freqs_func(j) =2.0*M_PI*freq(j);
			para.amps_func(j) =amp(j);
			para.phases_func(j) =phase(j);
		}

		// write output to file
		nlohmann::json json_freq_tmp;
		nlohmann::json json_amp_tmp;
		nlohmann::json json_phase_tmp;
		for(int j = 0; j < Ncomp; j++) {
			json_freq_tmp.push_back(freq(j));
			json_amp_tmp.push_back(amp(j));
			json_phase_tmp.push_back(phase(j));
		}
		json_parameters["frequency"].push_back(json_freq_tmp);
		json_parameters["amp"].push_back(json_amp_tmp);
		json_parameters["phase"].push_back(json_phase_tmp);

		for(int j= 0; j < Ncomp_intf; j++) {
			para.freqs_intf(j) =freq_intf(j);
			para.freqs_intf_func(j) =2.0*M_PI*freq_intf(j);
			para.amps_intf_func(j) =amp_intf(j);
			para.phases_intf_func(j) =phase_intf(j);
		}

		json_freq_tmp.clear();
		json_amp_tmp.clear();
		json_phase_tmp.clear();
		for(int j = 0; j < Ncomp_intf; j++) {
			json_freq_tmp.push_back(freq_intf(j));
			json_amp_tmp.push_back(amp_intf(j));
			json_phase_tmp.push_back(phase_intf(j));
		}
		json_parameters["frequency_inft"].push_back(json_freq_tmp);
		json_parameters["amp_inft"].push_back(json_amp_tmp);
		json_parameters["phase_inft"].push_back(json_phase_tmp);
	}
	// zero signal
	else if(para.signalflag==-1){
		// do nothing
	}

	// calculate signal mean and variance
	para.signalavrg=0;
	para.signalvar=0;
	for(int t= 0; t < Tmax_test; t++) {
		para.signalavrg+=func_target(t*para.dt, para);
		para.signalvar+=func_target(t*para.dt, para)*func_target(t*para.dt, para);
	}
	para.signalavrg=para.signalavrg/double(Tmax_test);
	para.signalvar=para.signalvar/double(Tmax_test) - para.signalavrg*para.signalavrg;

	// write targets to files
	nlohmann::json json_traces;
	nlohmann::json json_tmp;
	for(int t= 0; t < para.Tmax_learn; t++){
		if ((t % para.write_every)==0) json_tmp.push_back(func_target(t*para.dt, para));
	}
	json_traces["target_learn"].push_back(json_tmp);
	json_tmp.clear();

	for(int t= 0; t < para.Tmax_test; t++){
		if ((t % para.write_every)==0) json_tmp.push_back(func_target(t*para.dt, para));
	}
	json_traces["target_test"].push_back(json_tmp);
	json_tmp.clear();

	if (para.learnflag<99) {
		for(int t= 0; t < para.Tmax_FT; t++){
			if ((t % para.write_every)==0) json_tmp.push_back(func_target(t*para.dt, para));
		}
		json_traces["target_FT"].push_back(json_tmp);
		json_tmp.clear();
	}
	if (para.learnflag<99 && para.intrfflag==1) {
		nlohmann::json json_tmp;
		for(int j= 0; j < para.Ncomp_intf; j++) {
			for(int t= 0; t < para.Tmax_FT; t++){
				if ((t % para.write_every)==0) json_tmp.push_back(func_target_intf(t*para.dt, j, para));
			}
			json_traces["target_FT_interf"].push_back(json_tmp);
			json_tmp.clear();
		}
	}

	// write time to files
	for(int t= 0; t < para.Tmax_learn; t++){
		if ((t % para.write_every)==0) json_traces["time_learn"].push_back(t*para.dt);
	}
	for(int t= 0; t < para.Tmax_test; t++){
		if ((t % para.write_every)==0) json_traces["time_test"].push_back(t*para.dt);
	}
	for(int t= 0; t < para.Tmax_FT; t++){
		if ((t % para.write_every)==0) json_traces["time_FT"].push_back(t*para.dt);
	}

	std::ofstream traces_file(para.folder+"targets.json");		if (!traces_file) exit(EXIT_FAILURE);
	traces_file << std::setw(4) << json_traces << '\n';
	traces_file.close();

	return;
}

