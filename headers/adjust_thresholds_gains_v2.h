void adjust_thresholds_gains(arma::mat& RATE_Fadj_hGC, arma::mat& RATE_Fadj_nuCx, netpara &para) {

	arma::mat RATE_Fadj_GC(para.Ngc, para.Tmax_Fadj);
	////////////////////////////////////////////////////
	// adjust GC thresholds globally
	if(para.fcodingflag==1){
		arma::vec FC={1.0-para.fc_avrg};
		arma::vec Q = arma::quantile(arma::vectorise(RATE_Fadj_hGC),FC);
		para.GCthreshold=0;
		for(int i= 0; i < para.Ngc; i++){
			para.threshold_vector[i] =Q[0];
			para.GCthreshold+=para.threshold_vector[i];
		}
		para.GCthreshold=para.GCthreshold/double(para.Ngc);
	}
	////////////////////////////////////////////////////
	// adjust GC thresholds locally
	else if (para.fcodingflag==2) {
		para.GCthreshold=0;
		for(int i= 0; i < para.Ngc; i++){
			arma::vec FC={1.0-para.fc_target[i]};
			arma::vec Q = arma::quantile(RATE_Fadj_hGC.row(i),FC);
			para.threshold_vector[i] =Q[0];
			para.GCthreshold+=para.threshold_vector[i];
		}
		para.GCthreshold=para.GCthreshold/double(para.Ngc);
	}
	if(para.actfuncflag_GC==2 || para.actfuncflag_GC==3) {
		para.GCthreshold=0;
		para.threshold_vector.zeros();
	}
	////////////////////////////////////////////////////
	// obtain GC activity for gain adjustment
	for(int t= 0; t < para.Tmax_Fadj; t++){
		for(int i= 0; i < para.Ngc; i++) {
			RATE_Fadj_GC.at(i, t)=(*para.Ptr_gcIO) (RATE_Fadj_hGC.at(i, t)-para.threshold_vector[i],para.sp_alpha);
		}
	}

	////////////////////////////////////////////////////
	// global gain adjustment
	if(para.gainflag==1){
		para.GCgain=arma::stddev(arma::vectorise(RATE_Fadj_hGC))/arma::stddev(arma::vectorise(RATE_Fadj_GC));
		if(para.GCgain>para.GCgainmax)		para.GCgain=para.GCgainmax;
		for(int i= 0; i < para.Ngc; i++) {
			para.gain_vector[i]=para.GCgain;
		}
		//~ std::cout << arma::stddev(arma::vectorise(RATE_Fadj_hGC)) << "\t"<<arma::stddev(arma::vectorise(RATE_Fadj_GC)) << "\t"<<para.GCgain << std::endl;
	}
	////////////////////////////////////////////////////
	// local gain adjustment (match individual SDs)
	else if(para.gainflag==2){
		for(int i= 0; i < para.Ngc; i++) {
			para.gain_vector[i]=arma::stddev(RATE_Fadj_hGC.row(i) )/arma::stddev(RATE_Fadj_GC.row(i));
			//~ para.gain_vector[i]=arma::stddev(RATE_Fadj_nuCx.row(i) )/arma::stddev(RATE_Fadj_GC.row(i));
		}
		// set non-finite elements to zero
		para.gain_vector.elem( arma::find_nonfinite(para.gain_vector) ).zeros();
		// remove extreme gain values
		//~ double meanGain=arma::mean(para.gain_vector);
		//~ double Gaindeviation=meanGain+para.tol_GCgain*arma::stddev(para.gain_vector);
		double medianGain=arma::median(para.gain_vector);
		double Gaindeviation=para.tol_GCgain*medianGain;
		for(int i= 0; i < para.Ngc; i++) {
			if(para.gain_vector[i]>Gaindeviation)		para.gain_vector[i]=medianGain;
			if(para.gain_vector[i]>para.GCgainmax)		para.gain_vector[i]=para.GCgainmax;
		}
	}
	////////////////////////////////////////////////////
	// fixed gains
	else if(para.gainflag==3){
		para.GCgain=para.GCgain_base;
		if(para.GCgain>para.GCgainmax)		para.GCgain=para.GCgainmax;
		for(int i= 0; i < para.Ngc; i++) {
			para.gain_vector[i]=para.GCgain;
			//~ para.gain_vector[i]=para.GCgain+0.2*arma::randn();
		}
	}
	////////////////////////////////////////////////////
	// gains = inverse 1st spanning vector norm
	else if(para.gainflag==4){
		//~ FTcomp FT_Fadj_comp_cx;
		//~ calc_FT(RATE_Fadj_nuCx, FT_Fadj_comp_cx, para, para.Tmax_Fadj,"cereb");
		FTcomp FT_Fadj_comp_gc;
		calc_FT(RATE_Fadj_GC, FT_Fadj_comp_gc, para, para.Tmax_Fadj,"cereb");
		if(FT_Fadj_comp_gc.nn_signal(0)<0.000001) {
			for(int i= 0; i < para.Ngc; i++) {
				para.gain_vector[i]=0.0;
			}
		}
		else{
			for(int i= 0; i < para.Ngc; i++) {
				para.gain_vector[i]=0.42/FT_Fadj_comp_gc.nn_signal(0);
				if(para.gain_vector[i]>para.GCgainmax)		para.gain_vector[i]=para.GCgainmax;
			}
		}
	}
	////////////////////////////////////////////////////
	// gain adjusted such as to match 1st spanning vector norms
	else if(para.gainflag==5){

		FTcomp FT_Fadj_comp_cx;
		calc_FT(RATE_Fadj_nuCx, FT_Fadj_comp_cx, para, para.Tmax_Fadj,"cereb");
		FTcomp FT_Fadj_comp_gc;
		calc_FT(RATE_Fadj_GC, FT_Fadj_comp_gc, para, para.Tmax_Fadj,"cereb");
		if(FT_Fadj_comp_gc.nn_signal(0)<0.000001) {
			for(int i= 0; i < para.Ngc; i++) {
				para.gain_vector[i]=0.0;
			}
		}
		else{
			for(int i= 0; i < para.Ngc; i++) {
				para.gain_vector[i]=FT_Fadj_comp_cx.nn_signal(0)/FT_Fadj_comp_gc.nn_signal(0);
				if(para.gain_vector[i]>para.GCgainmax)		para.gain_vector[i]=para.GCgainmax;
			}
		}
	}
	////////////////////////////////////////////////////
	// first step: local gain adjustment (match individual SDs)
	// second step: increase/decrease gain globally such as to match 1st spanning vector norms
	//~ else if(para.gainflag==6){
		//~ // local adjustment (match individual SDs)
		//~ for(int i= 0; i < para.Ngc; i++) {
			//~ para.gain_vector[i]=arma::stddev(RATE_Fadj_hGC.row(i) )/arma::stddev(RATE_Fadj_GC.row(i));
		//~ }
		//~ // set non-finite elements to zero
		//~ para.gain_vector.elem( arma::find_nonfinite(para.gain_vector) ).zeros();
		//~ // remove extreme gain values
		//~ double medianGain=arma::median(para.gain_vector);
		//~ double Gaindeviation=para.tol_GCgain*medianGain;
		//~ for(int i= 0; i < para.Ngc; i++) {
			//~ if(para.gain_vector[i]>Gaindeviation) para.gain_vector[i]=medianGain;
		//~ }

		//~ // re-calculate GC activity for 2nd gain adjustment
		//~ for(int t= 0; t < para.Tmax_Fadj; t++){
			//~ for(int i= 0; i < para.Ngc; i++) {
				//~ RATE_Fadj_GC.at(i, t)=para.gain_vector[i]*(*para.Ptr_gcIO) (RATE_Fadj_hGC.at(i, t)-para.threshold_vector[i],para.sp_alpha);
			//~ }
		//~ }

		//~ // global adjustment to match 1st spanning vector
		//~ FTcomp FT_Fadj_comp_cx;
		//~ calc_FT(RATE_Fadj_nuCx, FT_Fadj_comp_cx, para, para.Tmax_Fadj);
		//~ FTcomp FT_Fadj_comp_gc;
		//~ calc_FT(RATE_Fadj_GC, FT_Fadj_comp_gc, para, para.Tmax_Fadj);
		//~ if(FT_Fadj_comp_gc.nn_signal(0)<0.000001) {
			//~ for(int i= 0; i < para.Ngc; i++) {
				//~ para.gain_vector[i]=0.0;
			//~ }
		//~ }
		//~ else{
			//~ for(int i= 0; i < para.Ngc; i++) {
				//~ para.gain_vector[i]=para.gain_vector[i]*FT_Fadj_comp_cx.nn_signal(0)/FT_Fadj_comp_gc.nn_signal(0);
			//~ }
		//~ }

	//~ }
	////////////////////////////////////////////////////
	// first step: local gain adjustment (match every FC to RNN population SD)
	// second step: increase/decrease gain globally such as to match 1st spanning vector norms
	else if(para.gainflag==6){
		// local adjustment (match RNN population SD)
		double popSD=arma::stddev(arma::vectorise(RATE_Fadj_hGC));
		for(int i= 0; i < para.Ngc; i++) {
			para.gain_vector[i]=popSD/arma::stddev(RATE_Fadj_GC.row(i));
		}
		// set non-finite elements to zero
		para.gain_vector.elem( arma::find_nonfinite(para.gain_vector) ).zeros();
		// remove extreme gain values
		double medianGain=arma::median(para.gain_vector);
		double Gaindeviation=para.tol_GCgain*medianGain;
		for(int i= 0; i < para.Ngc; i++) {
			if(para.gain_vector[i]>Gaindeviation)		para.gain_vector[i]=medianGain;
			if(para.gain_vector[i]>para.GCgainmax)		para.gain_vector[i]=para.GCgainmax;
		}

		// re-calculate GC activity for 2nd gain adjustment
		for(int t= 0; t < para.Tmax_Fadj; t++){
			for(int i= 0; i < para.Ngc; i++) {
				RATE_Fadj_GC.at(i, t)=para.gain_vector[i]*(*para.Ptr_gcIO) (RATE_Fadj_hGC.at(i, t)-para.threshold_vector[i],para.sp_alpha);
			}
		}

		// global adjustment to match 1st spanning vector
		FTcomp FT_Fadj_comp_cx;
		calc_FT(RATE_Fadj_nuCx, FT_Fadj_comp_cx, para, para.Tmax_Fadj,"cereb");
		FTcomp FT_Fadj_comp_gc;
		calc_FT(RATE_Fadj_GC, FT_Fadj_comp_gc, para, para.Tmax_Fadj,"cereb");
		if(FT_Fadj_comp_gc.nn_signal(0)<0.000001) {
			for(int i= 0; i < para.Ngc; i++) {
				para.gain_vector[i]=0.0;
			}
		}
		else{
			for(int i= 0; i < para.Ngc; i++) {
				para.gain_vector[i]=para.gain_vector[i]*FT_Fadj_comp_cx.nn_signal(0)/FT_Fadj_comp_gc.nn_signal(0);
			}
		}

	}

	////////////////////////////////////////////////////
	// set gains so as to set the maximum rate of each GC to GCmax
	else if(para.gainflag==7){
		// find max
		for(int i= 0; i < para.Ngc; i++) {
			double GCmax=arma::max(RATE_Fadj_GC.row(i));
			para.gain_vector[i]=para.GCmax/GCmax;
			//~ std::cout << "GCgain = " << para.gain_vector[i]<< "\t GCmax = " << GCmax << std::endl;
		}
		// set non-finite elements to zero
		para.gain_vector.elem( arma::find_nonfinite(para.gain_vector) ).zeros();
		// remove extreme gain values
		double medianGain=arma::median(para.gain_vector);
		double Gaindeviation=para.tol_GCgain*medianGain;
		for(int i= 0; i < para.Ngc; i++) {
			if(para.gain_vector[i]>Gaindeviation)		para.gain_vector[i]=medianGain;
			if(para.gain_vector[i]>para.GCgainmax)		para.gain_vector[i]=para.GCgainmax;
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
			para.gain_vector[i]=RNNmax.at(0)/GCmax;
			//~ std::cout << "GCgain = " << para.gain_vector[i]<< "\t GCmax = " << GCmax << std::endl;
		}
		// set non-finite elements to zero
		para.gain_vector.elem( arma::find_nonfinite(para.gain_vector) ).zeros();
		// remove extreme gain values
		double medianGain=arma::median(para.gain_vector);
		double Gaindeviation=para.tol_GCgain*medianGain;
		for(int i= 0; i < para.Ngc; i++) {
			if(para.gain_vector[i]>Gaindeviation)		para.gain_vector[i]=medianGain;
			if(para.gain_vector[i]>para.GCgainmax)		para.gain_vector[i]=para.GCgainmax;
		}

	}

	return;
}
