void RLS_step_FFL(struct RLS &RLS_cereb, double &time, netpara &para) {

	int Ngc=para.Ngc;
	double rPr,Ktemp,c;
	double K[Ngc];

	// update inverse correlation matrix
	for(int i= 0; i < Ngc; i++) {
		K[i]=0;
		for(int j= 0; j < Ngc; j++) {
			K[i]+=RLS_cereb.P.at(i,j)*para.nuGC[j];
		}
	}

	rPr=0;
	for(int i= 0; i < Ngc; i++) {
		rPr+=para.nuGC[i]*K[i];
	}

	c = 1.0/(1.0 + rPr);

	for(int i= 0; i < Ngc; i++) {
		Ktemp=K[i]*c;
		for(int j= 0; j < Ngc; j++) {
			RLS_cereb.P.at(i,j)-= Ktemp*K[j];
		}
	}

	// update the error of the linear readout unit
	para.error_signed_train = para.PC-func_target(time, para);

	// update the output weights
	for (int i=0; i<Ngc; i++){
		RLS_cereb.dwout.at(i) = -para.error_signed_train*K[i]*c;
		para.Wout_PC[i] = para.Wout_PC[i] + RLS_cereb.dwout[i]*double(Ngc);
	}

	return;
}
