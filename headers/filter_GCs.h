arma::colvec filter_GCs(arma::colvec& RATE, filter &butter, netpara &para) {

	int N=RATE.n_rows;			// number of GCs
	int Nt=butter.b.n_elem;		// number of filter coefficients

	arma::vec RATEout=RATE;
	for (int n = 0; n < N; n++) {
		RATEout.at(n)=butter.b.at(0)*RATE.at(n) + butter.w.at(n,0,0);
		for(int i=0; i<Nt-1; i++){
			butter.w.at(n,i,1)=butter.b.at(i+1)*RATE.at(n) + butter.w.at(n,i+1,0) - butter.a.at(i+1)*RATEout.at(n);
		}
		for(int i=0; i<Nt-1; i++){
			butter.w.at(n,i,0)=butter.w.at(n,i,1);
		}
	}

	return RATEout;
}
