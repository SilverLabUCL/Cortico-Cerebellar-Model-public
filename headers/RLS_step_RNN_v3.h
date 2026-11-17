void RLS_step_RNN(struct RLS &RLS_ctx, double &time, netpara &para) {

	int Ncx=para.Ncx;
	double rPr,Ktemp,c;
	double K[Ncx];

	// update inverse correlation matrix
	for(int i= 0; i < Ncx; i++) {
		K[i]=0;
		for(int j= 0; j < Ncx; j++) {
			K[i]+=RLS_ctx.P.at(i,j)*para.nuCx[j];
		}
	}

	rPr=0;
	for(int i= 0; i < Ncx; i++) {
		rPr+=para.nuCx[i]*K[i];
	}

	c = 1.0/(1.0 + rPr);

	for(int i= 0; i < Ncx; i++) {
		Ktemp=K[i]*c;
		for(int j= 0; j < Ncx; j++) {
			RLS_ctx.P.at(i,j)-= Ktemp*K[j];
		}
	}

	// update the error of the linear readout unit
	para.error_signed_train = para.Z-func_target(time, para);

	// update the output weights
	for (int i=0; i<Ncx; i++){
		RLS_ctx.dwout[i] = -para.error_signed_train*K[i]*c;
		para.Wout[i] = para.Wout[i] + RLS_ctx.dwout[i]*double(Ncx);
	}

	return;
}
