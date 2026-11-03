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

void target_signal(netpara &para, std::ofstream &freq_file, std::ofstream &amp_file, std::ofstream &phase_file, std::ofstream &freq_intf_file, std::ofstream &amp_intf_file, std::ofstream &phase_intf_file, std::ofstream &alg_steps_file) {

	int Tmax_learn=para.Tmax_learn;
	int Tmax_test=para.Tmax_test;
	int Tmax_FT=para.Tmax_FT;
	int Tmax_trans=para.Tmax_trans;
	int Tmax_FM=para.Tmax_FM;

	int Ncomp=para.Ncomp;
	int Ncomp_intf=para.Ncomp_intf;

	para.target=arma::zeros<arma::vec>(Tmax_learn);
	para.test=arma::zeros<arma::vec>(Tmax_test);
	para.fourier=arma::zeros<arma::vec>(Tmax_FT);
	para.trans=arma::zeros<arma::vec>(Tmax_trans);
	para.FM=arma::zeros<arma::vec>(Tmax_FM);

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
		alg_steps_file << "\t generate cosine target signal(s)"<< std::endl;
	}
	else{
		alg_steps_file << "\t generate target signal(s)"<< std::endl;
	}

	// Mix of 4 sine waves from Sussillo & Abbott 2009, Fig 2d; the original signal is obtained for lambda=0.8;
	if (para.signalflag==0){
		double amp = 1.6;
		for(int t= 0; t < Tmax_learn; t++) {
			para.target[t] = 	(amp/1.0)*sin(2.0*M_PI*para.lambda*t*para.dt) + (amp/2.0)*sin(2.0*M_PI*2*para.lambda*t*para.dt)
								+ (amp/6.0)*sin(2.0*M_PI*3*para.lambda*t*para.dt) + (amp/3.0)*sin(2.0*M_PI*4*para.lambda*t*para.dt);
			para.target[t] = para.target[t]/1.5;
		}
		for(int t= Tmax_learn; t < (Tmax_learn+Tmax_test); t++) {
			para.test(t-Tmax_learn) = (amp/1.0)*sin(2.0*M_PI*para.lambda*t*para.dt) + (amp/2.0)*sin(2.0*M_PI*2*para.lambda*t*para.dt)
								+ (amp/6.0)*sin(2.0*M_PI*3*para.lambda*t*para.dt) + (amp/3.0)*sin(2.0*M_PI*4*para.lambda*t*para.dt);
			para.test(t-Tmax_learn) = para.test(t-Tmax_learn)/1.5;
		}
		for(int t= 0; t < Tmax_trans; t++) {
			para.trans[t] = 	(amp/1.0)*sin(2.0*M_PI*para.lambda*t*para.dt) + (amp/2.0)*sin(2.0*M_PI*2*para.lambda*t*para.dt)
								+ (amp/6.0)*sin(2.0*M_PI*3*para.lambda*t*para.dt) + (amp/3.0)*sin(2.0*M_PI*4*para.lambda*t*para.dt);
			para.trans[t] = para.trans[t]/1.5;
		}
		for(int t= Tmax_trans; t < (Tmax_trans+Tmax_FT); t++) {
			para.fourier(t-Tmax_trans) = 	(amp/1.0)*sin(2.0*M_PI*para.lambda*t*para.dt) + (amp/2.0)*sin(2.0*M_PI*2*para.lambda*t*para.dt)
								+ (amp/6.0)*sin(2.0*M_PI*3*para.lambda*t*para.dt) + (amp/3.0)*sin(2.0*M_PI*4*para.lambda*t*para.dt);
			para.fourier(t-Tmax_trans) = para.fourier(t-Tmax_trans)/1.5;
		}
		for(int t= Tmax_trans; t < (Tmax_trans+Tmax_FM); t++) {
			para.FM(t-Tmax_trans) = 	(amp/1.0)*sin(2.0*M_PI*para.lambda*t*para.dt) + (amp/2.0)*sin(2.0*M_PI*2*para.lambda*t*para.dt)
								+ (amp/6.0)*sin(2.0*M_PI*3*para.lambda*t*para.dt) + (amp/3.0)*sin(2.0*M_PI*4*para.lambda*t*para.dt);
			para.FM(t-Tmax_trans) = para.FM(t-Tmax_trans)/1.5;
		}

		para.freqs(0) =para.lambda;
		para.freqs(1) =2*para.lambda;
		para.freqs(2) =3*para.lambda;
		para.freqs(3) =4*para.lambda;

		// write to files
		freq_file << para.freqs(0) << "\t" << para.freqs(1) << "\t" << para.freqs(2) << "\t" << para.freqs(3) << "\t";
		amp_file << amp/1.5 << "\t" << amp/(2.0*1.5) << "\t" << amp/(6.0*1.5) << "\t" << amp/(3.0*1.5) << "\t";
		phase_file << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t";
		freq_file << std::endl;
		amp_file << std::endl;
		phase_file << std::endl;

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

		// generate cosine waves
		for(int j= 0; j < Ncomp; j++) {
			w=2.0*M_PI*freq(j);
			for(int t= Tmax_learn; t < (Tmax_learn+Tmax_test); t++) {
				para.test(t-Tmax_learn) += amp(j)*cos(w*t*para.dt  + phase(j));
			}
		}

		// write to files
		for(int j= 0; j < Ncomp; j++) {
			freq_file << freq(j) << "\t";
			amp_file << amp(j) << "\t";
			phase_file << phase(j) << "\t";
			para.freqs(j) =freq(j);
			para.freqs_func(j) =2.0*M_PI*freq(j);
			para.amps_func(j) =amp(j);
			para.phases_func(j) =phase(j);
		}
		freq_file << std::endl;
		amp_file << std::endl;
		phase_file << std::endl;

		for(int j= 0; j < Ncomp_intf; j++) {
			freq_intf_file << freq_intf(j) << "\t";
			amp_intf_file << amp_intf(j) << "\t";
			phase_intf_file << para.phases_intf(j) << "\t";
			para.freqs_intf(j) =freq_intf(j);
			para.freqs_intf_func(j) =2.0*M_PI*freq_intf(j);
			para.amps_intf_func(j) =amp_intf(j);
			para.phases_intf_func(j) =phase_intf(j);
		}
		freq_intf_file << std::endl;
		amp_intf_file << std::endl;
		phase_intf_file << std::endl;
	}
	// zero signal
	else if(para.signalflag==-1){
		// do nothing
	}

	// calculate signal mean and variance
	para.signalavrg=0;
	para.signalvar=0;
	for(int t= Tmax_learn; t < (Tmax_learn+Tmax_test); t++) {
		para.signalavrg+=para.test(t-Tmax_learn);
		para.signalvar+=para.test(t-Tmax_learn)*para.test(t-Tmax_learn);
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

	std::ofstream traces_file(para.folder+"traces.dat");		if (!traces_file) exit(EXIT_FAILURE);
	traces_file << std::setw(4) << json_traces << '\n';
	traces_file.close();

	return;
}

