void DELTA_step(arma::vec& RATE, arma::vec& Wout, double &time, double &et, netpara &para) {

	int N=RATE.n_elem;

	// output unit
	double ROunit=0;
	for(int i= 0; i < N; i++) {
		ROunit+=RATE[i]*Wout[i];
	}

	// update the error of the linear readout unit
	para.error_signed_train = ROunit-func_target(time, para);

	// eligibility trace
	//~ et = et + 1*para.dt*(-et + para.error_signed_train);

	double eta_error=para.dt*para.error_signed_train*para.delta_eta/para.tau_avrg;
	//~ double eta_error=para.dt*et*para.etaini/para.tau_avrg;
	// update the output weights
	for (int i=0; i<N; i++){
		Wout[i] = Wout[i] -RATE[i]*eta_error;						// global learning rate
	}

	// update the output weights with AdaDelta
	//~ double dwout;
	//~ double gg;
	//~ double eps=1e-6;
	//~ double rhoD=0.999;
	//~ for (int i=0; i<N; i++){
		//~ // compute "gradient"
		//~ gg=para.error_signed_train*RATE[i];
		//~ // accumulate "gradient"
		//~ RLS_data.EG2[i] = rhoD*RLS_data.EG2[i] +(1.0-rhoD)* gg*gg;
		//~ // compute weight update
		//~ dwout = -gg*etaini*sqrt(RLS_data.Edw2[i] + eps)/sqrt(RLS_data.EG2[i] + eps);
		//~ // accumulate weight update
		//~ RLS_data.Edw2[i] = rhoD*RLS_data.Edw2[i] +(1.0-rhoD)* dwout*dwout;
		//~ // apply weight update
		//~ Wout[i] = Wout[i] + dwout;
	//~ }

	// update global learning rate like in Sussillo & Abbott 2009
	//~ RLS_data.LR = RLS_data.LR +para.dt/para.tau_avrg *( -RLS_data.LR*RLS_data.LR + RLS_data.LR*para.error_signed_train*para.error_signed_train );
	//~ if(RLS_data.LR>etaini)	RLS_data.LR=etaini;

	// update global learning rate like in Hoerzer et al. 2012
	//~ RLS_data.LR = etaini/(1+t*para.dt/(20.0*para.lambda));

	return;
}
