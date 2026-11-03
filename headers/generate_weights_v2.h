void generate_weights(presyn Wcx[], presyn Wfb[], presyn Wfb_interf[], presyn Wmf[],  presyn &Wout, presyn &Wout_PC, int Ncx_to_MF_idx[], arma::colvec& ftarget, arma::colvec& zeta, netpara &para, gsl_rng * r) {

	int Ncx=para.Ncx;
	int Ngc=para.Ngc;
	int Nmf=para.Nmf;
	int Nsyn=para.Nsyn;
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
			// check for reciprocal connection
			//~ if (gsl_rng_uniform (r)<para.precip) {
				//~ if(gsl_rng_uniform (r)<para.pcon){
					//~ Adjacency.at(i,j)=1;
					//~ Adjacency.at(j,i)=1;
					//~ Weights.at(i,j)=gsl_ran_gaussian_ziggurat (r, sigmaJ);
					//~ Weights.at(j,i)=arma::sign(Weights.at(i,j))*std::fabs(gsl_ran_gaussian_ziggurat (r, sigmaJ));	// sgn(W(j,i)) = sgn(W(i,j))
				//~ }
				//~ else{
					//~ Adjacency.at(i,j)=0;
					//~ Adjacency.at(j,i)=0;
				//~ }
			//~ }
			//~ else{
				if(gsl_rng_uniform (r)<para.pcon){
					Adjacency.at(i,j)=1;
					Weights.at(i,j)=gsl_ran_gaussian_ziggurat (r, sigmaJ);
				}
				else{
					Adjacency.at(i,j)=0;
				}
			//~ }
		}
	}

	dummyADJ.clear(); dummyW.clear();
	for(int i= 0; i < Ncx; i++) {
		int k=0;
		for (int j=0; j<Ncx; j++) {
			if(Adjacency(i,j)==1){
				k++;
				dummyADJ.push_back (j);
				dummyW.push_back (Weights.at(i,j));
			}
		}
		Wcx[i].Nsyn=k;
		Wcx[i].idx=arma::zeros<arma::uvec>(k);
		Wcx[i].weight=arma::zeros<arma::vec>(k);
		for (int j=0; j<k; j++) {
			Wcx[i].idx[j]=dummyADJ[j];
			//~ Wcx[i].weight[j]=gsl_ran_gaussian_ziggurat (r, sigmaJ);
			Wcx[i].weight[j]=dummyW[j];
			// Wcx[i].weight[j]=a*heavi(a);
			// if(i==j)	Wcx[i].weight[j]=0;
		}
		dummyADJ.clear(); dummyW.clear();
	}
	//~ for(int i= 0; i < 3; i++) {
		//~ std::cout << Wcx[i].Nsyn << "\t" << Wcx[i].idx.n_elem << "\t" << Wcx[i].weight.n_elem << std::endl;
	//~ }

	//*******************************
	// modifiable weights
	//*******************************
	// cortex to readout connections
	Wout.weight=arma::zeros<arma::vec>(Ncx);
	if(para.cortexonlyflag==0){
		// GC to PC connections
		Wout_PC.weight=arma::zeros<arma::vec>(Ngc);
	}

	//*******************************
	// readout to cortex connections (feedback)
	//*******************************
	double temp;
	int Nfbgroups=round(Ncx/para.Nfb);
	for(int i= 0; i < Ncx; i++) {
		Wfb[i].weight=arma::zeros<arma::vec>(1);
		Wfb_interf[i].weight=arma::zeros<arma::vec>(para.Ncomp_intf);
	}
	//~ for(i= 0; i < Ncx; i++) {
		//~ if(i%Nfbgroups==0) 	{
			//~ temp= gsl_ran_flat(r, -1, 1);
		//~ }
		//~ // Wfb[i].weight[k]= gsl_ran_flat(r, -1, 1);
		//~ Wfb[i].weight[0]= temp;

		//~ // scale weights by number of feedback units
		//~ // Wfb[i].weight[0]= gsl_ran_flat(r, -1, 1);

		//~ // non-zero weight with probability 1/Npc
		//~ // if(gsl_ran_flat(r, 0, 1)<(1.0/double(Npc))) {
			//~ // Wfb[i].weight[k]= gsl_ran_flat(r, -1, 1);
		//~ // }
		//~ // else{
			//~ // Wfb[i].weight[0]= 0.0;
		//~ // }
		//~ // if(i%2==0) 	Wfb[i].weight[0]= 0.5;
		//~ // else  		Wfb[i].weight[0]= -0.5;
	//~ }
	// connect different readout units/PCs to different recurrent units
	if(para.intrfflag==0){
		for(int i= 0; i < Ncx; i++) {
			if(i%Nfbgroups==0) 	{
				temp = gsl_ran_gaussian_ziggurat (r, 1.0);
			}
			else temp=0;
			Wfb[i].weight[0]= temp;
			// set interference feedback weights identical to regular feedback weights
			//~ Wfb_interf[i].weight[0]= Wfb[i].weight[0];
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
		//~ if(para.intrfflag==1){
			for(int i= 0; i < Ncx; i++) {
				gsl_ran_multivariate_gaussian(r, mu, L, res);
				for(int k= 0; k < para.Ncomp_intf+1; k++) {
					//~ std::cout <<  gsl_vector_get(res,k) << '\t' ;
					if(k==0)	Wfb[i].weight[k]=gsl_vector_get(res,k);
					else 		Wfb_interf[i].weight[k-1]=gsl_vector_get(res,k);
				}
				//~ std::cout << std::endl;
			}
		//~ }
		// for cerebellar interference signal
		//~ else if (para.intrfflag==2){
			//~ for(int i= 0; i < Ngc; i++) {
				//~ gsl_ran_multivariate_gaussian(r, mu, L, res);
				//~ for(int k= 0; k < para.Ncomp_intf+1; k++) {
					//~ std::cout <<  gsl_vector_get(res,k) << '\t' ;
					//~ if(k==0)	Wfb[i].weight[k]=gsl_vector_get(res,k);
					//~ else 		Wfb_interf[i].weight[k-1]=gsl_vector_get(res,k);
				//~ }
				//~ std::cout << std::endl;
			//~ }
		//~ }

		gsl_vector_free (mu);
		gsl_matrix_free (L);
		gsl_vector_free (res);
	}

	//*******************************
	// MF identity & MF-GC connections
	//*******************************
	if(para.cortexonlyflag==0){
		if(para.simplemapflag==0){
			// MF identity
			int *shffl_arr_Ncx=new int[Ncx];
			for (int j=0;j<Ncx;j++) shffl_arr_Ncx[j]=j;
			gsl_ran_shuffle (r, shffl_arr_Ncx, Ncx, sizeof (int));
			for (int j=0;j<Nmf;j++) Ncx_to_MF_idx[j]=shffl_arr_Ncx[j];
			delete [] shffl_arr_Ncx;

			// MF to GC connections: draw exactly Nsyn distinct random connections
			int *shffl_arr=new int[Nmf];
			for (int j=0;j<Nmf;j++) shffl_arr[j]=j;
			for (int i=0;i<Ngc;i++){
				gsl_ran_shuffle (r, shffl_arr, Nmf, sizeof (int));
				Wmf[i].idx=arma::zeros<arma::uvec>(Nsyn);
				Wmf[i].weight=arma::zeros<arma::vec>(Nsyn);
				for (int j=0;j<Nsyn;j++) {
					//~ Wmf[i].weight[j]=1.0;
					Wmf[i].weight[j]=gsl_ran_gaussian_ziggurat (r, sigmaJ);
					Wmf[i].idx[j]=shffl_arr[j];
					Wmf[i].Nsyn=Nsyn;
				}
			}
			delete [] shffl_arr;
		}
		// one to one ctx GC mapping (i.e. no mixing, requires Nctx=Ngc=Nmf and Nsyn=1)
		else if(para.simplemapflag==1){
			for (int j=0;j<Nmf;j++) Ncx_to_MF_idx[j]=j;
			for (int i=0;i<Ngc;i++){
				Wmf[i].idx=i;
				Wmf[i].weight=1.0;
				Wmf[i].Nsyn=1;
			}
		}
		//*******************************
		// target coding levels
		//*******************************
		if(para.fcodingflag==2){
			// check whether f_target distribution is well defined
			if( (para.fc_SD*para.fc_SD) >= para.fc_avrg*(1-para.fc_avrg) ){
				std::cerr << "f_target variance too large." << std::endl;
				std::cerr << std::endl;
				exit(EXIT_FAILURE);
			}
			double a=(1.0-para.fc_avrg)*para.fc_avrg*para.fc_avrg/(para.fc_SD*para.fc_SD) - para.fc_avrg;
			double b=a*( 1.0/para.fc_avrg -1);

			for (int i=0;i<Ngc;i++){
				ftarget[i]=gsl_ran_beta(r, a, b);
			}
		}
		//*******************************
		// Golgi-cell to GC weights
		//*******************************
		para.Wgoc_distrb=arma::zeros<arma::vec>(para.Ngc);
		if(para.gocmodeflag==2){
			double kk=para.Wgoc*para.Wgoc/(para.Wgoc_SD*para.Wgoc_SD);
			double theta=para.Wgoc_SD*para.Wgoc_SD/para.Wgoc;
			for (int i=0;i<Ngc;i++){
				para.Wgoc_distrb[i]=gsl_ran_gamma(r, kk, theta);
			}
		}
		else if(para.gocmodeflag==1 || para.gocmodeflag==3){
			for (int i=0;i<Ngc;i++){
				para.Wgoc_distrb[i]=para.Wgoc;
			}
		}
	}

	//*******************************
	// distributed taus
	//*******************************
	para.tau_distb=arma::ones<arma::vec>(Ncx)*para.tau_avrg;
	if(para.taudistrbflag==1){
		double taucutoff=0.003;
		para.tau_distb=arma::zeros<arma::vec>(para.Ncx);
		double kk=para.tau_avrg*para.tau_avrg/(para.tau_std*para.tau_std);
		double theta=para.tau_std*para.tau_std/para.tau_avrg;
		for (int i=0;i<Ncx;i++){
			while(para.tau_distb[i]<taucutoff){
				para.tau_distb[i]=gsl_ran_gamma(r, kk, theta);
			}
		}
	}
	para.sqt=sqrt(para.dt);
	para.dttau_distb=para.dt/para.tau_distb;
	para.IIdt_distb=(1.0-para.dttau_distb);
	para.sqDT=para.sqt/sqrt(para.tau_distb);

	//*******************************
	// distributed GC taus
	//*******************************
	para.tauGC_distb=arma::ones<arma::vec>(Ngc)*para.tauGC_avrg;
	if(para.taugcmodeflag==2){
		// gamma distribution fitted to data from Schwartz, Rothman et al. 2012
		double taucutoff=0.0014389;		// smallest tauGC encountered in Schwartz, Rothman et al. 2012
		para.tauGC_distb=arma::zeros<arma::vec>(Ngc);
		double kk=para.tauGC_avrg*para.tauGC_avrg/(para.tauGC_std*para.tauGC_std);
		double theta=para.tauGC_std*para.tauGC_std/para.tauGC_avrg;
		for (int i=0;i<Ngc;i++){
			while(para.tauGC_distb[i]<taucutoff){
				para.tauGC_distb[i]=gsl_ran_gamma(r, kk, theta);
			}
		}
	}
	else{
		para.tauGC_distb=arma::zeros<arma::vec>(Ngc);
	}
	para.dttauGC_distb=para.dt/para.tauGC_distb;
	para.IIdtGC_distb=(1.0-para.dttauGC_distb);
	para.sqtauGC_distb=sqrt(para.tauGC_distb);
	para.sqDTGC=para.sqt/para.sqtauGC_distb;

	//*******************************
	// random Gaussian vector zeta
	//*******************************
	if(para.spanmodflag==5){
		for (int i=0;i<Ncx;i++){
			zeta.at(i)=gsl_ran_gaussian_ziggurat (r, 1.0);
		}
		zeta=zeta/arma::norm(zeta);
	}

	return;
}
