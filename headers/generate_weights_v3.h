void generate_weights(netpara &para) {

	int Ncx=para.Ncx;
	int Ngc=para.Ngc;
	int Nmf=para.Nmf;
	double sigmaJ=para.g/sqrt(para.pcon*para.Ncx);
	std::vector <int> dummyADJ;
	std::vector <double> dummyW;
	//*******************************
	// recurrent cortical connections
	//*******************************
	arma::mat Adjacency= arma::zeros<arma::mat>(Ncx,Ncx);
	arma::mat Weights= arma::zeros<arma::mat>(Ncx,Ncx);
	for(int i= 0; i < Ncx; i++) {
		for (int j=0; j<Ncx; j++) {
			if(i==j) continue;				// no self interaction
			if(gsl_rng_uniform (para.r)<para.pcon){
				Adjacency.at(i,j)=1;
				Weights.at(i,j)=gsl_ran_gaussian_ziggurat (para.r, sigmaJ);
			}
		}
	}

	dummyADJ.clear(); dummyW.clear();
	para.Wcx.clear();
	para.Wcx.reserve(Ncx);
	for(int i= 0; i < Ncx; i++) {
		int k=0;
		for (int j=0; j<Ncx; j++) {
			if(Adjacency(i,j)==1){
				k++;
				dummyADJ.push_back (j);
				dummyW.push_back (Weights.at(i,j));
			}
		}
		arma::uvec tmp_idx=arma::zeros<arma::uvec>(k);
		arma::vec tmp_weight=arma::zeros<arma::vec>(k);
		for (int j=0; j<k; j++) {
			tmp_idx[j]=dummyADJ[j];
			tmp_weight[j]=dummyW[j];
		}
		para.Wcx.emplace_back(tmp_idx, tmp_weight);
		dummyADJ.clear(); dummyW.clear();
	}

	//*******************************
	// modifiable weights
	//*******************************
	// cortex to readout connections
	para.Wout=arma::zeros<arma::vec>(Ncx);
	if(para.cortexonlyflag==0){
		// GC to PC connections
		para.Wout_PC=arma::zeros<arma::vec>(Ngc);
	}

	//*******************************
	// readout to cortex connections (feedback)
	//*******************************
	para.Wfb=arma::zeros<arma::vec>(Ncx);
	para.Wfb_interf=arma::zeros<arma::mat>(Ncx,para.Ncomp_intf);
	// connect different readout units/PCs to different recurrent units
	if(para.intrfflag==0){
		for(int i= 0; i < Ncx; i++) {
			para.Wfb[i]= gsl_ran_gaussian_ziggurat (para.r, 1.0);
		}
	}
	// 	generate correlation between feedback weights
	else{
		gsl_vector * res= gsl_vector_alloc (para.Ncomp_intf+1);
		gsl_vector * mu = gsl_vector_alloc (para.Ncomp_intf+1);
		gsl_vector_set_zero(mu);
		gsl_matrix *L = gsl_matrix_alloc(para.Ncomp_intf+1, para.Ncomp_intf+1);

		arma::mat RHO;
		arma::mat Lchol;
		RHO.eye(para.Ncomp_intf+1, para.Ncomp_intf+1);

		for(int j= 0; j < para.Ncomp_intf+1; j++) {
			for(int k= 0; k < j; k++) {
				RHO(j,k)=para.rho_Wfb(j,k);
				RHO(k,j)=para.rho_Wfb(j,k);
			}
		}

		bool ok = arma::chol(Lchol, RHO);
		//~ std::cout << ok << std::endl;

		// deal with non positive semi-definite matrices
		if(ok==0){
			nearestSPD(RHO);
			ok = arma::chol(Lchol, RHO);
			//~ std::cout << ok << std::endl;
		}

		for(int j= 0; j < para.Ncomp_intf+1; j++) {
			for(int k= 0; k < para.Ncomp_intf+1; k++) {
				gsl_matrix_set(L, j, k, Lchol(k,j));
				gsl_matrix_set(L, k, j, Lchol(j,k));
			}
		}

		// for cortical interference signal   // FINISH THIS, i.e. declare size of Wfb_interf properly
		for(int i= 0; i < Ncx; i++) {
			gsl_ran_multivariate_gaussian(para.r, mu, L, res);
			for(int k= 0; k < para.Ncomp_intf+1; k++) {
				//~ std::cout <<  gsl_vector_get(res,k) << '\t' ;
				if(k==0)	para.Wfb[i]=gsl_vector_get(res,k);
				else 		para.Wfb_interf.at(i,k-1)=gsl_vector_get(res,k);
			}
			//~ std::cout << std::endl;
		}

		gsl_vector_free (mu);
		gsl_matrix_free (L);
		gsl_vector_free (res);
	}

	//*******************************
	// MF identity & MF-GC connections
	//*******************************
	para.Wmf.clear();
	para.Wmf.reserve(Ngc);
	para.Ncx_to_MF_idx=arma::zeros<arma::uvec>(Nmf);
	if(para.cortexonlyflag==0){
		arma::uvec tmp_idx=arma::zeros<arma::uvec>(para.Nsyn);
		arma::vec tmp_weight=arma::zeros<arma::vec>(para.Nsyn);
		if(para.simplemapflag==0){
			// MF identity: every MF represents a distinct cortical unit
			std::vector<int>shffl_arr_Ncx(Ncx);
			for (int j=0;j<Ncx;j++) shffl_arr_Ncx[j]=j;
			gsl_ran_shuffle (para.r, shffl_arr_Ncx.data(), Ncx, sizeof (int));
			for (int j=0;j<Nmf;j++) para.Ncx_to_MF_idx[j]=shffl_arr_Ncx[j];

			// MF to GC connections: draw exactly Nsyn distinct random connections
			std::vector<int>shffl_arr(Nmf);
			for (int j=0;j<Nmf;j++) shffl_arr[j]=j;
			for (int i=0;i<Ngc;i++){
				gsl_ran_shuffle (para.r, shffl_arr.data(), Nmf, sizeof (int));
				for (int j=0;j<para.Nsyn;j++) {
					tmp_weight[j]=1.0;
					tmp_idx[j]=shffl_arr[j];
				}
				para.Wmf.emplace_back(tmp_idx, tmp_weight);
			}
		}
		// one-to-one ctx to GC mapping (i.e. no mixing, requires Nctx=Ngc=Nmf and Nsyn=1)
		else if(para.simplemapflag==1){
			for (int j=0;j<Nmf;j++) para.Ncx_to_MF_idx[j]=j;
			for (int i=0;i<Ngc;i++){
				arma::uvec tmp_idx(1); tmp_idx[0] = i;
				arma::vec tmp_weight(1); tmp_weight[0] = 1.0;
				para.Wmf.emplace_back(tmp_idx, tmp_weight);
			}
		}

		//*******************************
		// target coding levels
		//*******************************
		para.fc_target=arma::zeros<arma::vec>(Ngc);
		if(para.fcodingflag==2){
			// check whether f_target distribution is well defined
			if(para.fc_avrg <= 0.0 || para.fc_avrg >= 1.0){
				std::cerr << "fc_avrg must be in (0,1)." << std::endl;
				exit(EXIT_FAILURE);
			}
			if( (para.fc_SD*para.fc_SD) >= para.fc_avrg*(1-para.fc_avrg) ){
				std::cerr << "f_target variance too large." << std::endl;
				std::cerr << std::endl;
				exit(EXIT_FAILURE);
			}
			double a=(1.0-para.fc_avrg)*para.fc_avrg*para.fc_avrg/(para.fc_SD*para.fc_SD) - para.fc_avrg;
			double b=a*( 1.0/para.fc_avrg -1);

			for (int i=0;i<Ngc;i++){
				para.fc_target[i]=gsl_ran_beta(para.r, a, b);
			}
		}
		//*******************************
		// Golgi-cell to GC weights
		//*******************************
		para.Wgoc_distrb=arma::zeros<arma::vec>(para.Ngc);
		if(para.gocmodeflag==2){
			double kk=para.Wgoc*para.Wgoc/(para.Wgoc_SD*para.Wgoc_SD);
			double theta_kk=para.Wgoc_SD*para.Wgoc_SD/para.Wgoc;
			for (int i=0;i<Ngc;i++){
				para.Wgoc_distrb[i]=gsl_ran_gamma(para.r, kk, theta_kk);
			}
		}
		else if(para.gocmodeflag==1 || para.gocmodeflag==3){
			for (int i=0;i<Ngc;i++){
				para.Wgoc_distrb[i]=para.Wgoc;
			}
		}
		para.Wgoc_distrb_scaled = para.Wgoc_distrb / double(para.Ncx);
	}

	//*******************************
	// distributed taus
	//*******************************
	para.tau_distrb=arma::ones<arma::vec>(Ncx)*para.tau_avrg;
	if(para.taudistrbflag==1){
		para.tau_distrb=arma::zeros<arma::vec>(para.Ncx);
		double kk=para.tau_avrg*para.tau_avrg/(para.tau_std*para.tau_std);
		double theta_kk=para.tau_std*para.tau_std/para.tau_avrg;
		for (int i=0;i<Ncx;i++){
			while(para.tau_distrb[i]<para.tauRNNcutoff){
				para.tau_distrb[i]=gsl_ran_gamma(para.r, kk, theta_kk);
			}
		}
	}
	para.sqt=sqrt(para.dt);
	para.dttau_distrb=para.dt/para.tau_distrb;
	para.IIdt_distb=(1.0-para.dttau_distrb);
	para.sqDT=para.sqt/sqrt(para.tau_distrb);

	//*******************************
	// distributed GC taus
	//*******************************
	para.tauGC_distrb=arma::ones<arma::vec>(Ngc)*para.tauGC_avrg;
	if(para.taugcmodeflag==2){
		para.tauGC_distrb=arma::zeros<arma::vec>(Ngc);
		double kk=para.tauGC_avrg*para.tauGC_avrg/(para.tauGC_std*para.tauGC_std);
		double theta_kk=para.tauGC_std*para.tauGC_std/para.tauGC_avrg;
		for (int i=0;i<Ngc;i++){
			while(para.tauGC_distrb[i]<para.tauGCcutoff){
				para.tauGC_distrb[i]=gsl_ran_gamma(para.r, kk, theta_kk);
			}
		}
	}
	para.dttauGC_distb=para.dt/para.tauGC_distrb;
	para.IIdtGC_distb=(1.0-para.dttauGC_distb);
	para.sqDTGC=para.sqt/sqrt(para.tauGC_distrb);

	return;
}
