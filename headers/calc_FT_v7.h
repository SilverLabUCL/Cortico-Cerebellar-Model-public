// centers spanning vectors
void center_vectors(arma::mat &Space, arma::mat &Space_center){
	double avrg;
	int M=Space.n_cols;
	Space_center=Space;
	for (int i = 0; i < M; i++) {
		avrg=arma::mean(Space_center.col(i));
		Space_center.col(i)=Space_center.col(i)-avrg;
	}
	return;
}

// calculates dimensionality from spanning vectors
double calc_FT_dim(arma::mat &Space){
	double numerator=0;
	double denominator=0;
	double dij,Snorm;
	int M=Space.n_cols;
	for (int i = 0; i < M; i++) {
		Snorm=arma::norm(Space.col(i));
		numerator+=Snorm*Snorm;
		denominator+=Snorm*Snorm*Snorm*Snorm;
		for (int j = 0; j < i; j++) {
			dij=arma::dot(Space.col(i),Space.col(j));
			denominator+= 2.0*dij*dij;
		}
	}
	return numerator*numerator/denominator;
}

// calculates angles between spanning vectors
void calc_FT_angles(arma::mat &Angles, arma::mat &Space){
	int M=Space.n_cols;
	Angles=arma::zeros<arma::mat>(M,M);
	arma::vec norms(M);
	for (int i = 0; i < M; i++) {
		norms(i)=arma::norm(Space.col(i));
		for (int j = 0; j < i; j++) {
			Angles.at(i,j)=std::acos(arma::dot(Space.col(i)/norms.at(i) , Space.col(j)/norms.at(j) ) ) * 180.0 /M_PI;
			Angles.at(j,i)=Angles.at(i,j);
		}
	}
	return;
}

// calculates condition number from non-zero spanning vectors
double calc_FT_cond(arma::mat &Space, double tol){

	double Snorm;
	arma::mat Space_clean;
	for (size_t i = 0; i < Space.n_cols; i++) {
		Snorm=arma::norm(Space.col(i));
		if(Snorm>tol){
			Space_clean=arma::join_horiz(Space_clean, Space.col(i));
		}
	}
	int M=Space_clean.n_cols;

	arma::mat CR(M,M);
	for (int i = 0; i < M; i++) {
		Snorm=arma::norm(Space_clean.col(i));
		CR.at(i,i)=Snorm*Snorm;
		for (int j = 0; j < i; j++) {
			CR.at(i,j)=arma::dot(Space_clean.col(i),Space_clean.col(j));
			CR.at(j,i)=CR.at(i,j);
		}
	}
	return arma::cond(CR);
}

// calculates spanning vetcor components via FT
void FTcomponents(FTcomp &Cmp, arma::mat &RATE, arma::vec &cleantol,netpara &para, int FT_tsteps, int Ncomp_total, int cleanflag, int writeflag, std::string modeflag){

	// baseline frequency
	double w=para.lambda;

	int N=RATE.n_rows;
	int Nthalf=FT_tsteps/2;

	double fs = 1.0/para.dt;		// max frequency

	arma::vec freqaxis(Nthalf);
	for (int i = 0; i < Nthalf; i++) {
		freqaxis(i) = i*fs/double(FT_tsteps);		// fs/double(FT_tsteps) is the frequency step
	}

	//~ std::cout << "# of time steps FT: " << FT_tsteps <<"\t maximum frequency: " << fs <<"\t \t base frequency: " << w <<"\t \t frequency step: " << fs/double(FT_tsteps) << std::endl;
	//~ std::cout << "number of components: " << Ncomp_total << "\t maximum searched frequency: " << Ncomp_total*w << std::endl;
	//~ std::cout << "freqaxis: " << freqaxis(0) <<"\t" << freqaxis(1) <<"\t" << freqaxis(2) <<"\t" << freqaxis(3) << std::endl;

	arma::mat pks_real(N,1);
	arma::mat pks_imag(N,1);

	arma::cx_mat FT =  arma::fft(RATE.t())/double(Nthalf);
	FT.shed_rows( Nthalf, FT_tsteps-1 );

	// clean up Fourier components
	if(cleanflag==1) {
		arma::vec P;
		for (int i = 0; i < N; i++) {
			arma::cx_vec tmp=FT.col(i);
			P=arma::abs(tmp) % arma::abs(tmp);
			for (int f = 0; f < Nthalf; f++) {
				if(P(f)<cleantol(i)) FT.at(f,i)=1e-6;
			}
		}

	}
	// calculate properties of error spectrum
	else{
		for (int i = 0; i < N; i++) {
			arma::cx_vec tmp=FT.col(i);
			arma::vec P=arma::abs(tmp.subvec( Nthalf/2-1, Nthalf-1 ));
			P=P % P;
			//~ if(i==0) std::cout << arma::size(FT) <<"\t" << Nthalf <<"\t" << arma::size(P) << std::endl;
			double tailmean=arma::mean(P);
			double tailstd=arma::stddev(P);
			//~ cleantol(i)=tailmean+4*tailstd;		// remove < 4 times STD
			cleantol(i)=tailmean+3*tailstd;		// remove < 3 times STD
			//~ cleantol(i)=tailmean+1*tailstd;		// remove < 1 times STD
		}
	}

	arma:: mat rr=arma::real(FT.t());
	arma:: mat ii=arma::imag(FT.t());

	if(para.learnflag<99 && writeflag==1) {
		std::ofstream FT_real_file(para.folder+"FT_real_" + modeflag +".txt");		if (!FT_real_file) exit(EXIT_FAILURE);
		rr.save(FT_real_file,arma::raw_ascii);
		FT_real_file.close();
		std::ofstream FT_imag_file(para.folder+"FT_imag_" + modeflag +".txt");		if (!FT_imag_file) exit(EXIT_FAILURE);
		ii.save(FT_imag_file,arma::raw_ascii);
		FT_imag_file.close();
	}

	Cmp.a0=rr.col(0)/2.0;
	// find peaks at integer multiples of the baseline frequency
	double eps_base=0.01;
	double eps;
	int emptyflag;
	arma::uvec tmp;
	for (int i = 0; i < Ncomp_total; i++) {
		int j=i+1;
		arma::uvec idx_loc(Nthalf);
		emptyflag=1;
		eps=eps_base;
		//~ std::cout << j <<"\t" << fs <<"\t" << fs/double(FT_tsteps) <<"\t" << j*w-eps <<"\t" << j*w <<"\t" << j*w+eps << std::endl;
		while(emptyflag==1){
			for (int k = 0; k < Nthalf; k++) {
				idx_loc(k)=( freqaxis(k)>=(j*w-eps) ) & ( freqaxis(k)<=(j*w+eps) );
			}
			tmp=arma::find(idx_loc>0);
			if(tmp.is_empty()){
				emptyflag=1;
				eps=eps*1.01;
				// std::cout << "is_empty" << std::endl;
			}
			else{
				emptyflag=0;
			}
		}
		idx_loc=arma::find(idx_loc>0);

		pks_real=rr.cols( idx_loc );
		pks_imag=ii.cols( idx_loc );
		Cmp.an.col(i)=pks_real;
		Cmp.bn.col(i)=pks_imag;
	}

	return;
}

void calc_FT(arma::mat RATE, FTcomp &Cmp, netpara &para, int FT_tsteps, std::string modeflag) {

	int Ncomp_noise;

	if (modeflag== "ctx") {
		// linear RNN without artificial overtones
		if(para.actfuncflag==1 && para.intrfflag==0)		Ncomp_noise=0;
		// linear RNN with artificial overtones
		else if(para.actfuncflag==1 && para.intrfflag==1)	Ncomp_noise=para.nf*(para.Ncomp_intf-1);
		// non-linear tanh RNN
		//~ else if(para.actfuncflag==0)	Ncomp_noise=6;
		else if(para.actfuncflag==0)	Ncomp_noise=para.Ncomp_noise_RNN;
		//~ Ncomp_noise=std::max(para.Ncomp_intf-1,0);
	}
    else if (modeflag == "cereb"){
		// linear GCs without artificial overtones
		if(para.actfuncflag_GC==2 && para.intrfflag==0)			Ncomp_noise=1;
		// linear GCs with artificial overtones
		else if(para.actfuncflag_GC==2 && para.intrfflag==2)	Ncomp_noise=para.Ncomp_intf-1;
		// non-linear GCs
		else Ncomp_noise=para.Ncomp_noise;
	}
	else {std::cerr <<"Bad modeflag string in 'calc_FT'. \n";exit(EXIT_FAILURE);}

	int Ncomp=para.Ncomp;
	int Ncomp_total=Ncomp_noise+Ncomp;

	int N=RATE.n_rows;
	double sqN=sqrt(N);

	Cmp.a0=arma::zeros<arma::vec>(N);
	Cmp.an=arma::zeros<arma::mat>(N,Ncomp_total);
	Cmp.bn=arma::zeros<arma::mat>(N,Ncomp_total);

	Cmp.S_signal=arma::zeros<arma::mat>(N,2*Ncomp);
	Cmp.S_noise=arma::zeros<arma::mat>(N,2*Ncomp_noise+1);

	//	*************************
	// 	Fourier transforms
	int writeflag;
	int cleanflag;
	arma::vec cleantol=arma::zeros<arma::vec>(N);

	if(para.FMdyncleanflag==1) {
		FTcomp Cmp_noise;
		Cmp_noise.a0=arma::zeros<arma::vec>(N);
		Cmp_noise.an=arma::zeros<arma::mat>(N,Ncomp_total);
		Cmp_noise.bn=arma::zeros<arma::mat>(N,Ncomp_total);
		arma::mat TMP=arma::zeros<arma::mat>(arma::size(RATE));
		TMP.randn(arma::size(RATE));
		 //~ std::cout << arma::stddev(arma::vectorise(TMP*para.sigmaLS)) << std::endl;
		arma::mat NOISE=RATE+TMP*para.sigmaLS;//*para.FACsigmaLS_ctx; // ADJUST THIS

		writeflag=1; cleanflag=0;
		FTcomponents(Cmp_noise, NOISE, cleantol, para, FT_tsteps, Ncomp_total,cleanflag,writeflag,modeflag+"_noise");

		writeflag=1; cleanflag=1;
		FTcomponents(Cmp, RATE, cleantol, para, FT_tsteps, Ncomp_total,cleanflag,writeflag,modeflag+"_cleaned");
	}
	else {
		writeflag=1; cleanflag=0;
		FTcomponents(Cmp, RATE, cleantol, para, FT_tsteps, Ncomp_total,cleanflag,writeflag,modeflag+"_cleaned");
	}

	//	*************************
	// clean up fourier components
	//~ double tol=1e-4;
	double tol=para.FCcomp_tol;
	Cmp.a0.clean( tol);
	Cmp.an.clean( tol);
	Cmp.bn.clean( tol);
	// set spanning vectors to zero if norm negligible
	//~ double tol2=1e-5;
	//~ double tol2=5e-4;
	//~ double tol2=1e-3;
	//~ double tol2=1e-12;
	//~ double na1=arma::norm(Cmp.an.col(0));
	// double nb1=arma::norm(Cmp.bn.col(0));
	//~ for (int i = 0; i < Ncomp_total; i++) {
		//~ if( arma::norm(Cmp.an.col(i))/na1 <tol2 ) {//|| arma::norm(Cmp.an.col(i)) < para.sigmaLS || arma::norm(Cmp.an.col(i)) < para.sigmaN){
			//~ Cmp.an.col(i)=arma::zeros(arma::size(Cmp.an.col(i)));
			//~ std::cout << "a_"<< i << std::endl;
		//~ }
		//~ // if( arma::norm(Cmp.bn.col(i))/nb1 <tol2){
		//~ if( arma::norm(Cmp.bn.col(i))/na1 <tol2 ) {//|| arma::norm(Cmp.bn.col(i)) < para.sigmaLS || arma::norm(Cmp.an.col(i)) < para.sigmaN){
			//~ Cmp.bn.col(i)=arma::zeros(arma::size(Cmp.bn.col(i)));
			//~ std::cout << "b_"<< i << std::endl;
		//~ }
	//~ }

	//	*************************
	// assemble signal and noise subspaces
	for (int i = 0; i < Ncomp; i++) {
		int k=2*i;
		Cmp.S_signal.col(k)=Cmp.an.col(i)/sqN;
		Cmp.S_signal.col(k+1)=Cmp.bn.col(i)/sqN;
	}
	Cmp.S_noise.col(0)=Cmp.a0/sqN;
	for (int i = Ncomp; i < Ncomp_total; i++) {
		int k=2*(i-Ncomp)+1;
		Cmp.S_noise.col(k)=Cmp.an.col(i)/sqN;
		Cmp.S_noise.col(k+1)=Cmp.bn.col(i)/sqN;
	}
	if (Ncomp_noise>0){
		Cmp.S_noise_noa0=arma::zeros<arma::mat>(N,2*Ncomp_noise);
		for (int i = Ncomp; i < Ncomp_total; i++) {
			int k=2*(i-Ncomp);
			Cmp.S_noise_noa0.col(k)=Cmp.an.col(i)/sqN;
			Cmp.S_noise_noa0.col(k+1)=Cmp.bn.col(i)/sqN;
		}
	}
	Cmp.S_total=arma::join_horiz(Cmp.a0/sqN, Cmp.S_signal);
	if (Ncomp_noise>0){
		Cmp.S_total=arma::join_horiz(Cmp.S_total, Cmp.S_noise_noa0);
	}

	// non-normalised (raw) vectors
	Cmp.S_total_sqN=Cmp.S_total*sqN;
	Cmp.S_signal_sqN=Cmp.S_signal*sqN;
	Cmp.S_noise_sqN=Cmp.S_noise*sqN;

	// **************************
	// center vectors
	center_vectors(Cmp.S_total,Cmp.S_total_centered);
	center_vectors(Cmp.S_signal,Cmp.S_signal_centered);
	center_vectors(Cmp.S_noise,Cmp.S_noise_centered);
	center_vectors(Cmp.S_noise_noa0,Cmp.S_noise_noa0_centered);
	//~ Cmp.S_total_centered=Cmp.S_total;
	//~ Cmp.S_signal_centered=Cmp.S_signal;
	//~ Cmp.S_noise_centered=Cmp.S_noise;
	//~ Cmp.S_noise_noa0_centered=Cmp.S_noise_noa0;

	int Ms=Cmp.S_signal.n_cols;
	int Mn=Cmp.S_noise.n_cols;
	Cmp.Ms=Ms;
	Cmp.Mn=Mn;

	//	*************************
	// norms of Fourier components
	Cmp.nn_signal = arma::zeros<arma::vec>(Ms);
	Cmp.nn_noise = arma::zeros<arma::vec>(Mn);
	Cmp.nn_total = arma::zeros<arma::vec>(Ms+Mn);

	for (int i = 0; i < Ms; i++) 		Cmp.nn_signal(i)=arma::norm(Cmp.S_signal_centered.col(i));
	for (int i = 0; i < Mn; i++) 		Cmp.nn_noise(i)=arma::norm(Cmp.S_noise_centered.col(i));
	for (int i = 0; i < (Ms+Mn); i++) 	Cmp.nn_total(i)=arma::norm(Cmp.S_total_centered.col(i));

	//	*************************
	// angles between individual spanning vectors
	calc_FT_angles(Cmp.Angles_all, Cmp.S_total_centered);

	//	*************************
	// signal and noise-spaces based on cosine-components
	arma::mat S_signal_cos;
	arma::mat S_noise_cos;
	for (int i = 0; i < Ms/2; i++) {
		S_signal_cos=arma::join_horiz(S_signal_cos, Cmp.S_signal_centered.col(2*i));
		S_noise_cos=arma::join_horiz(S_noise_cos, Cmp.S_signal_centered.col(2*i+1));
	}
	if (Ncomp_noise>0){
		S_noise_cos=arma::join_horiz(S_noise_cos, Cmp.S_noise_centered);
	}

	//	*************************
	// compute dimensionality from Fourier components
	// 'signal' dimensionality
	Cmp.dim_signal = calc_FT_dim(Cmp.S_signal_centered);
	//~ Cmp.dim_signal = calc_FT_dim(S_signal_cos);

	//	*************************
	// compute condition number from Fourier components
	//~ Cmp.cond_signal = calc_FT_cond(Cmp.S_signal_centered, 1e-3);
	Cmp.cond_signal = calc_FT_cond(S_signal_cos, 1e-3);

	//	*************************
	// compute pairwise dimensionality and pairwise condition number of signal cos/sin components
	Cmp.dim_signal_pw=arma::zeros<arma::vec>(Ncomp);
	for (int i = 0; i < Ncomp; i++) {
		int j=2*i;
		arma::mat join_sig=arma::join_horiz(Cmp.S_signal_centered.col(j), Cmp.S_signal_centered.col(j+1));
		Cmp.dim_signal_pw.at(i,j)=calc_FT_dim(join_sig);
	}

	if (modeflag == "cereb"){
		int maxindex=std::min(1, Ncomp_noise);
		arma::mat condtmp=arma::join_horiz(Cmp.S_signal_centered, Cmp.S_noise_noa0_centered.cols(0,maxindex-1));
		//~ arma::mat condtmp=arma::join_horiz(Cmp.S_signal_centered, Cmp.S_noise_noa0_centered.cols(0,Ncomp_noise-1));
		Cmp.cond_total=calc_FT_cond(condtmp,1e-3);
	}
	else Cmp.cond_total=NAN;
	//~ Cmp.cond_total=NAN;

	//	*************************
	// subspace angles between cos-part of signal-subspace and the rest
	double threshold=sqrt(2)/2;
	//~ Cmp.subangles=arma::zeros<arma::vec>(std::min(Ms/2,Mn));
	Cmp.subangles=arma::zeros<arma::vec>(Ncomp);

	// scale columns
	const double EPSILON=std::numeric_limits<double>::epsilon();
	for (size_t i = 0; i < S_signal_cos.n_cols; i++){
		double Snorm = arma::norm(S_signal_cos.col(i),"inf");
		if(Snorm > EPSILON)		S_signal_cos.col(i)=S_signal_cos.col(i)/Snorm;
	}
	for (size_t i = 0; i < S_noise_cos.n_cols; i++){
		double Snorm = arma::norm(S_noise_cos.col(i),"inf");
		if(Snorm > EPSILON)		S_noise_cos.col(i)=S_noise_cos.col(i)/Snorm;
	}

	arma::mat basis_signal_cos;
	bool bs_succ = arma::orth(basis_signal_cos, S_signal_cos);
	arma::mat basis_noise_cos;
	bool bn_succ = arma::orth(basis_noise_cos, S_noise_cos);

	// if bases are constructed successfully, calculate subspace angles
	if(bs_succ==1 && bn_succ==1 ) {
		arma::vec SV_cos; arma::mat Ys; arma::mat Zs;
		svd_econ(Ys, SV_cos, Zs, basis_signal_cos.t()*basis_noise_cos,"both","std");
		//~ Zs=-Zs;
		arma::vec cos_vals=arma::zeros<arma::vec>(Ncomp);
		size_t nAngles = std::min((size_t)SV_cos.n_elem, (size_t)Ncomp);
		for (size_t i = 0; i < nAngles; i++){
			cos_vals(i)=SV_cos(i);
			Cmp.subangles(i) = std::acos(std::clamp(SV_cos(i), 0.0, 1.0));
		}
		arma::uvec indexsmall = cos_vals > threshold;
		arma::uvec indexsmall_loc=arma::find(indexsmall>0);
		if(arma::any( indexsmall)) {
			arma::mat V = basis_noise_cos*Zs;
			arma::vec SV_cos_small;
			arma::mat Yx; arma::mat Zx;
			arma::mat RG = V.cols(indexsmall_loc);
			svd_econ(Yx, SV_cos_small, Zx, RG-basis_signal_cos*(basis_signal_cos.t()*RG) );
			for (size_t i = 0; i < SV_cos_small.n_elem; i++){
				Cmp.subangles(indexsmall_loc(i)) = std::asin(std::clamp(SV_cos_small(i), 0.0, 1.0));
			}
		}
		Cmp.subangles=arma::sort(Cmp.subangles* 180.0 /M_PI);
	}
	else {
		for (size_t i = 0; i < Cmp.subangles.n_elem; i++){
			Cmp.subangles(i)=arma::datum::nan;
		}
	}

	// calculate sine product of signal angles
	calc_FT_angles(Cmp.Angles_signal, S_signal_cos);
	//~ std::cout << arma::nonzeros(arma::trimatu(Cmp.Angles_signal)) << std::endl;
	Cmp.sineprod_signal=arma::prod(arma::sin(arma::nonzeros(arma::trimatu(Cmp.Angles_signal))*M_PI/180.0));

	return;
}
