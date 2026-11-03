void LS_solution(arma::mat RATE, arma::vec &WOUT, double &RMSE_train, arma::uword& rnk, netpara &para, double FACsigmaLS) {

	int N=RATE.n_rows;
	int Nt=RATE.n_cols;

	arma::vec targetrate(Nt);

	arma::uvec Nidx(N);
	arma::ucolvec status(N);

	// subsample target signal
	for(int t= 0; t < Nt; t++) {
		targetrate[t]= func_target(t*para.dt, para);
		//~ targetrate[t]= func_target(t, para); + gsl_ran_gaussian_ziggurat(para.r,0.01); 	// add noise to target
	}

	//********************************************************
	// Learning
	//********************************************************
	arma::vec Wls;
	bool learn_success;
	if(para.LRflag==2){
		if(para.sigmaLS>0){
			double regnoise=para.sigmaLS*FACsigmaLS;
			// add regularisation noise to rates
			for(int i= 0; i < N; i++) {
				for(int t= 0; t < Nt; t++) {
					RATE.at(i,t)=RATE.at(i,t) + gsl_ran_gaussian_ziggurat(para.r,regnoise);
				}
			}
		}

		// compute weights with least-squares
		learn_success=arma::solve(Wls, RATE.t(), targetrate);
		if(learn_success==0) {
			Wls.zeros(N,1);
		}

	}
	else if (para.LRflag==3){
		// compute weights with ridge regression
		Wls= arma::inv( RATE*RATE.t() + para.sigmaLS*para.sigmaLS*arma::eye(N,N) ) * (RATE *targetrate);
	}
	else if (para.LRflag==4){
		// compute weights with pseudo-inverse
		Wls= arma::pinv( RATE*RATE.t() ) * (RATE *targetrate);
	}

	// rank of XX^t
	bool rnk_success = arma::rank(rnk, RATE*RATE.t() );

	// compute output rate
	arma::vec outrate=arma::trans( Wls.t() * RATE);

	// update the error of the linear readout unit
	RMSE_train=0;
	for(int t= 0; t < Nt; t++) {
		RMSE_train+=(outrate[t]-targetrate[t])*(outrate[t]-targetrate[t]);
	}
	RMSE_train=sqrt(RMSE_train/double(Nt))/sqrt(para.signalvar);

	// update the output weights
	WOUT = Wls*N;

	return;
}
