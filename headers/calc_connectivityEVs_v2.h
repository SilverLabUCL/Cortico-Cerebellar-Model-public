void calc_connectivityEVs(netpara &para, nlohmann::json &json_measures, std::ofstream &alg_steps_file) {

	if( (para.LRflag==2 || para.LRflag==3 || para.LRflag==4) && para.actfuncflag==1){
		alg_steps_file << "\t \t \t Calculate connectivity EVs"<< std::endl;

		arma::mat Weff= arma::zeros<arma::mat>(para.Ncx,para.Ncx);
		for(int i= 0; i < para.Ncx; i++) {
			for (int j=0; j<para.Wcx[i].Nsyn; j++) {
				Weff.at(i,para.Wcx[i].idx[j])=para.Wcx[i].weight[j];
			}
		}
		if(para.actfuncflag_GC==2 && para.simplemapflag==1){
			// if GC activation function is linear and RNN-GC map is simple, use PC readout weights
			for(int i= 0; i < para.Ncx; i++) {
				for (int j=0; j<para.Ncx; j++) {
					Weff.at(i,j)=Weff.at(i,j) + para.Wfb[i]*para.Wout_PC[j];
				}
			}
		}
		else {
			// else compute RNN EVs
			for(int i= 0; i < para.Ncx; i++) {
				for (int j=0; j<para.Ncx; j++) {
					Weff.at(i,j)=Weff.at(i,j) + para.Wfb[i]*para.Wout[j];
					//~ Weff.at(i,j)=gsl_ran_gaussian_ziggurat (r, para.sigmaJ);
				}
			}
		}

		arma::cx_vec Weff_eigval = arma::eig_gen( Weff );
		Weff_eigval=Weff_eigval(arma::sort_index( arma::abs(Weff_eigval), "descend" ));

		arma::cx_vec outliers =Weff_eigval;
		arma::uvec rmidx=arma::find(arma::abs(outliers)<para.g*para.outlier_thresh);

		outliers.shed_rows( rmidx );
		//~ for(int i= 0; i < outliers.n_rows; i++) {
			//~ std::cout << std::real(outliers[i]) <<  "\t" << std::imag(outliers[i]) << std::endl;
		//~ }
		//~ std::cout << std::endl;
		for(size_t i= 0; i < outliers.n_rows; i++) {
			json_measures["EV_outliers_real_ctx"].push_back(std::real(outliers[i]));
			json_measures["EV_outliers_imag_ctx"].push_back(std::imag(outliers[i]));
		}
		for(int i= outliers.n_rows; i < 2*std::max(para.Ncomp,para.Ncomp_intf); i++) {
			json_measures["EV_outliers_real_ctx"].push_back(NAN);
			json_measures["EV_outliers_imag_ctx"].push_back(NAN);
		}
		json_measures["EV_Noutliers_ctx"].push_back(outliers.n_rows);
	}
	else{
		for(int i= 0; i < 2*std::max(para.Ncomp,para.Ncomp_intf); i++) {
			json_measures["EV_outliers_real_ctx"].push_back(NAN);
			json_measures["EV_outliers_imag_ctx"].push_back(NAN);
		}
		json_measures["EV_Noutliers_ctx"].push_back(NAN);
	}

	return;
}
