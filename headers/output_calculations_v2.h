////////////////////////////////////////////////////////////////////////////////
double output_rate(arma::vec& WOUT, arma::colvec& INPUT) {

	double OUTPUT=0;
	for (size_t i=0; i<INPUT.n_rows; i++) 	OUTPUT+=INPUT[i]*WOUT[i];
	OUTPUT=OUTPUT/double(INPUT.n_rows);

	return OUTPUT;
}
////////////////////////////////////////////////////////////////////////////////
//~ double calc_corr(arma::mat& RATE) {
	//~ //	*************************
	//~ double avrgcorr;

	//~ // calculate correlation matrix and vectorise
	//~ arma::vec C=arma::vectorise(arma::cor(RATE.t()));

	//~ // remove NaN entries (due to zero neurons)
	//~ arma::uvec nanidx=find_nonfinite(C);
	//~ C.shed_rows(nanidx);
	// std::cout  << C.n_rows << "\t" << C.n_cols << std::endl;
	//~ // return average correlation
	//~ if(C.n_rows==0)		avrgcorr=NAN;
	//~ else  				avrgcorr=arma::mean(arma::abs(C));
	//~ return avrgcorr;
//~ }
////////////////////////////////////////////////////////////////////////////////
double calc_dim(arma::mat& RATE, std::vector<double>& lambdaT, int SVD_tsteps) {
	//	*************************
	// calculate dimensionality
	// 	do temporal PCA
	int N=RATE.n_rows;
	int Nlamb=std::min(SVD_tsteps,N);
	// compute average rate for each unit and centre activity matrix
	double AVRG[N];
	for (int i = 0; i < N; i++) {
		AVRG[i]=0;
		for (int j = 0; j < SVD_tsteps; j++) {
			AVRG[i]+=RATE.at(i, j);
		}
		AVRG[i]=AVRG[i]/double(SVD_tsteps);
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < SVD_tsteps; j++) {
			RATE.at(i, j)=RATE.at(i, j)-AVRG[i];
		}
	}
	// compute singular values
	arma::vec ST = arma::svd(RATE);
	//~ std::cout << Nlamb <<"\t" << ST.n_rows <<"\t" << ST.n_cols << std::endl;
	// compute eigenvalues
	for (int k = 0; k < Nlamb; k++) {
		lambdaT[k]=ST[k]*ST[k]/(double(SVD_tsteps)-1);
	}
	// compute dimensionality
	double a,b;
	a=b=0;
	for (int i = 0; i < Nlamb; i++) {
		a+=lambdaT[i];
		b+=lambdaT[i] *lambdaT[i];
	}
	// return dimensionality
	return a*a/b;
}
////////////////////////////////////////////////////////////////////////////////
double calc_dim_noEV(arma::mat& RATE, int SVD_tsteps) {
	//	*************************
	// calculate dimensionality
	double dim;
	// 	do temporal PCA
	int N=RATE.n_rows;
	int Nlamb=std::min(SVD_tsteps,N);
	// compute average rate for each unit and centre activity matrix
	double AVRG[N];
	for (int i = 0; i < N; i++) {
		AVRG[i]=0;
		for (int j = 0; j < SVD_tsteps; j++) {
			AVRG[i]+=RATE.at(i, j);
		}
		AVRG[i]=AVRG[i]/double(SVD_tsteps);
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < SVD_tsteps; j++) {
			RATE.at(i, j)=RATE.at(i, j)-AVRG[i];
		}
	}
	// compute singular values
	arma::vec ST;
	bool success = arma::svd(ST,RATE);
	//~ std::cout << Nlamb <<"\t" << ST.n_rows <<"\t" << ST.n_cols << std::endl;
	if(success) {
		// compute eigenvalues
		double lambdaT[Nlamb];
		for (int k = 0; k < Nlamb; k++) {
			lambdaT[k]=ST[k]*ST[k]/(double(SVD_tsteps)-1);
		}
		// compute dimensionality
		double a,b;
		a=b=0;
		for (int i = 0; i < Nlamb; i++) {
			a+=lambdaT[i];
			b+=lambdaT[i] *lambdaT[i];
		}
		// return dimensionality
		dim=a*a/b;
	}
	else dim=NAN;
	return dim;
}
////////////////////////////////////////////////////////////////////////////////
inline double calc_wnorm(arma::vec& Wout) {

	return arma::norm(Wout);
}
