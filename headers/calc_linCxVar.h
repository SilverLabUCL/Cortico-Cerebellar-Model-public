double calc_linCxVar(netpara &para) {
	// valid for sum over cosines

	int Ncomp=para.Ncomp;
	// baseline frequency
	double lambda=para.lambda;
	double tau=para.tau;
	double g=para.g;
	double Cxvar=0;
	double w=2*M_PI*lambda*tau;
	for(int k= 1; k <= Ncomp; k++) {
		Cxvar+=0.5/(1+k*k*w*w-g*g);
	}
	//~ std::cout << tau << "\t "<< lambda << "\t "<<Cxvar << std::endl;
	return Cxvar/double(Ncomp*Ncomp);
}
