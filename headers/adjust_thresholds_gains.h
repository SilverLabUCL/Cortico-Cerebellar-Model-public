void adjust_thresholds_gains(arma::mat& RATE_Fadj_hGC, arma::mat& RATE_Fadj_nuCx, arma::colvec& gain, arma::colvec& theta, arma::colvec& ftarget, int Tmax_Fadj, netpara &para) {

	int FT_Fadj_tsteps=ceil(Tmax_Fadj*1.0/double(para.SVD_FT_every));

	arma::mat RATE_Fadj_GC(para.Ngc, Tmax_Fadj);
	arma::mat RATE_Fadj_GC_FT(para.Ngc, FT_Fadj_tsteps);
	arma::mat RATE_Fadj_nuCx_FT(para.Ncx, FT_Fadj_tsteps);
	////////////////////////////////////////////////////
	// adjust GC thresholds globally
	if(para.fcodingflag==1){
		arma::vec FC={1.0-para.fc_avrg};
		arma::vec Q = arma::quantile(arma::vectorise(RATE_Fadj_hGC),FC);
		para.GCthreshold=0;
		for(int i= 0; i < para.Ngc; i++){
			theta[i] =Q[0];
			para.GCthreshold+=theta[i];
		}
		para.GCthreshold=para.GCthreshold/double(para.Ngc);
	}
	////////////////////////////////////////////////////
	// adjust GC thresholds locally
	else if (para.fcodingflag==2) {
		para.GCthreshold=0;
		for(int i= 0; i < para.Ngc; i++){
			arma::vec FC={1.0-ftarget[i]};
			arma::vec Q = arma::quantile(RATE_Fadj_hGC.row(i),FC);
			theta[i] =Q[0];
			para.GCthreshold+=theta[i];
		}
		para.GCthreshold=para.GCthreshold/double(para.Ngc);
	}
	if(para.actfuncflag_GC==2 || para.actfuncflag_GC==3) {
		para.GCthreshold=0;
		for(int i= 0; i < para.Ngc; i++){
			theta[i] =0;
		}
	}
	////////////////////////////////////////////////////
	// obtain GC activity for gain adjustment
	for(int t= 0; t < Tmax_Fadj; t++){
		for(int i= 0; i < para.Ngc; i++) {
			RATE_Fadj_GC.at(i, t)=(*para.Ptr_gcIO) (RATE_Fadj_hGC.at(i, t)-theta[i],para.sp_alpha);
		}
	}

	////////////////////////////////////////////////////
	// global gain adjustment
	if(para.gainflag==1){
		para.GCgain=arma::stddev(arma::vectorise(RATE_Fadj_hGC))/arma::stddev(arma::vectorise(RATE_Fadj_GC));
		if(para.GCgain>para.GCgainmax)		para.GCgain=para.GCgainmax;
		for(int i= 0; i < para.Ngc; i++) {
			gain[i]=para.GCgain;
		}
	}
	////////////////////////////////////////////////////
	// local gain adjustment (match individual SDs)
	else if(para.gainflag==2){
		for(int i= 0; i < para.Ngc; i++) {
			gain[i]=arma::stddev(RATE_Fadj_hGC.row(i) )/arma::stddev(RATE_Fadj_GC.row(i));
		}
		// set non-finite elements to zero
		gain.elem( arma::find_nonfinite(gain) ).zeros();
		// remove extreme gain values
		double medianGain=arma::median(gain);
		double Gaindeviation=para.tol_GCgain*medianGain;
		for(int i= 0; i < para.Ngc; i++) {
			if(gain[i]>Gaindeviation)		gain[i]=medianGain;
			if(gain[i]>para.GCgainmax)		gain[i]=para.GCgainmax;
		}
	}
	////////////////////////////////////////////////////
	// fixed gains
	else if(para.gainflag==3){
		para.GCgain=para.GCgain_base;
		if(para.GCgain>para.GCgainmax)		para.GCgain=para.GCgainmax;
		for(int i= 0; i < para.Ngc; i++) {
			gain[i]=para.GCgain;
			//~ gain[i]=para.GCgain+0.2*arma::randn();
		}
	}
	////////////////////////////////////////////////////
	// gains = inverse 1st spanning vector norm
	else if(para.gainflag==4){
		// downsample arrays to get data for FT
		int tFT=0;
		for(int t= 0; t < Tmax_Fadj; t++){
			bool bool_time_FT=( (t % para.SVD_FT_every)==0);
			if (bool_time_FT && (tFT< FT_Fadj_tsteps) ){
				for(int i= 0; i < para.Ngc; i++) RATE_Fadj_GC_FT.at(i, tFT)= RATE_Fadj_GC.at(i, t);
				tFT++;
			}
		}

		FTcomp FT_Fadj_comp_gc;
		calc_FT(RATE_Fadj_GC_FT, FT_Fadj_comp_gc, para, FT_Fadj_tsteps,"cereb");
		if(FT_Fadj_comp_gc.nn_signal(0)<0.000001) {
			for(int i= 0; i < para.Ngc; i++) {
				gain[i]=0.0;
			}
		}
		else{
			for(int i= 0; i < para.Ngc; i++) {
				gain[i]=0.42/FT_Fadj_comp_gc.nn_signal(0);
				if(gain[i]>para.GCgainmax)		gain[i]=para.GCgainmax;
			}
		}
	}
	////////////////////////////////////////////////////
	// gain adjusted such as to match 1st spanning vector norms
	else if(para.gainflag==5){
		// downsample arrays to get data for FT
		int tFT=0;
		for(int t= 0; t < Tmax_Fadj; t++){
			bool bool_time_FT=( (t % para.SVD_FT_every)==0);
			if (bool_time_FT && (tFT< FT_Fadj_tsteps) ){
				for(int i= 0; i < para.Ngc; i++) RATE_Fadj_GC_FT.at(i, tFT)= RATE_Fadj_GC.at(i, t);
				for(int i= 0; i < para.Ncx; i++) RATE_Fadj_nuCx_FT.at(i, tFT)= RATE_Fadj_nuCx.at(i, t);
				tFT++;
			}
		}

		FTcomp FT_Fadj_comp_cx;
		calc_FT(RATE_Fadj_nuCx_FT, FT_Fadj_comp_cx, para, FT_Fadj_tsteps,"cereb");
		FTcomp FT_Fadj_comp_gc;
		calc_FT(RATE_Fadj_GC_FT, FT_Fadj_comp_gc, para, FT_Fadj_tsteps,"cereb");
		if(FT_Fadj_comp_gc.nn_signal(0)<0.000001) {
			for(int i= 0; i < para.Ngc; i++) {
				gain[i]=0.0;
			}
		}
		else{
			for(int i= 0; i < para.Ngc; i++) {
				gain[i]=FT_Fadj_comp_cx.nn_signal(0)/FT_Fadj_comp_gc.nn_signal(0);
				if(gain[i]>para.GCgainmax)		gain[i]=para.GCgainmax;
			}
		}
	}
	////////////////////////////////////////////////////
	// first step: local gain adjustment (match every FC to RNN population SD)
	// second step: increase/decrease gain globally such as to match 1st spanning vector norms
	else if(para.gainflag==6){
		// local adjustment (match RNN population SD)
		double popSD=arma::stddev(arma::vectorise(RATE_Fadj_hGC));
		for(int i= 0; i < para.Ngc; i++) {
			gain[i]=popSD/arma::stddev(RATE_Fadj_GC.row(i));
		}
		// set non-finite elements to zero
		gain.elem( arma::find_nonfinite(gain) ).zeros();
		// remove extreme gain values
		double medianGain=arma::median(gain);
		double Gaindeviation=para.tol_GCgain*medianGain;
		for(int i= 0; i < para.Ngc; i++) {
			if(gain[i]>Gaindeviation)		gain[i]=medianGain;
			if(gain[i]>para.GCgainmax)		gain[i]=para.GCgainmax;
		}

		// re-calculate GC activity for 2nd gain adjustment
		for(int t= 0; t < Tmax_Fadj; t++){
			for(int i= 0; i < para.Ngc; i++) {
				RATE_Fadj_GC.at(i, t)=gain[i]*(*para.Ptr_gcIO) (RATE_Fadj_hGC.at(i, t)-theta[i],para.sp_alpha);
			}
		}
		// downsample arrays to get data for FT
		int tFT=0;
		for(int t= 0; t < Tmax_Fadj; t++){
			bool bool_time_FT=( (t % para.SVD_FT_every)==0);
			if (bool_time_FT && (tFT< FT_Fadj_tsteps) ){
				for(int i= 0; i < para.Ngc; i++) RATE_Fadj_GC_FT.at(i, tFT)= RATE_Fadj_GC.at(i, t);
				for(int i= 0; i < para.Ncx; i++) RATE_Fadj_nuCx_FT.at(i, tFT)= RATE_Fadj_nuCx.at(i, t);
				tFT++;
			}
		}

		// global adjustment to match 1st spanning vector
		FTcomp FT_Fadj_comp_cx;
		calc_FT(RATE_Fadj_nuCx_FT, FT_Fadj_comp_cx, para, FT_Fadj_tsteps,"cereb");
		FTcomp FT_Fadj_comp_gc;
		calc_FT(RATE_Fadj_GC_FT, FT_Fadj_comp_gc, para, FT_Fadj_tsteps,"cereb");
		if(FT_Fadj_comp_gc.nn_signal(0)<0.000001) {
			for(int i= 0; i < para.Ngc; i++) {
				gain[i]=0.0;
			}
		}
		else{
			for(int i= 0; i < para.Ngc; i++) {
				gain[i]=gain[i]*FT_Fadj_comp_cx.nn_signal(0)/FT_Fadj_comp_gc.nn_signal(0);
			}
		}

	}

	////////////////////////////////////////////////////
	// set gains so as to set the maximum rate of each GC to GCmax
	else if(para.gainflag==7){
		// find max
		for(int i= 0; i < para.Ngc; i++) {
			double GCmax=arma::max(RATE_Fadj_GC.row(i));
			gain[i]=para.GCmax/GCmax;
		}
		// set non-finite elements to zero
		gain.elem( arma::find_nonfinite(gain) ).zeros();
		// remove extreme gain values
		double medianGain=arma::median(gain);
		double Gaindeviation=para.tol_GCgain*medianGain;
		for(int i= 0; i < para.Ngc; i++) {
			if(gain[i]>Gaindeviation)		gain[i]=medianGain;
			if(gain[i]>para.GCgainmax)		gain[i]=para.GCgainmax;
		}
	}

		////////////////////////////////////////////////////
	// set gains so as to match the maximum rates of GCs and RNN units
	else if(para.gainflag==8){
		// find max of RNN
		arma::vec RNNmax(para.Ncx);
		for(int i= 0; i < para.Ncx; i++) {
			RNNmax.at(i)=arma::max(arma::abs(RATE_Fadj_nuCx.row(i)));
		}
		// find max
		for(int i= 0; i < para.Ngc; i++) {
			double GCmax=arma::max(RATE_Fadj_GC.row(i));
			RNNmax=arma::shuffle(RNNmax);
			gain[i]=RNNmax.at(0)/GCmax;
		}
		// set non-finite elements to zero
		gain.elem( arma::find_nonfinite(gain) ).zeros();
		// remove extreme gain values
		double medianGain=arma::median(gain);
		double Gaindeviation=para.tol_GCgain*medianGain;
		for(int i= 0; i < para.Ngc; i++) {
			if(gain[i]>Gaindeviation)		gain[i]=medianGain;
			if(gain[i]>para.GCgainmax)		gain[i]=para.GCgainmax;
		}

	}

	return;
}
