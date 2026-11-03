void write_parameters(nlohmann::json &json_parameters, netpara &para) {

	int nn= para.Ncomp_intf+1;
	arma::vec rho_Wfb_vec=arma::zeros<arma::vec>((nn*nn - nn)/2 +nn);
	int i=0;
	for(int j= 0; j < nn; j++) {
		for(int k= 0; k <= j; k++) {
			rho_Wfb_vec(i)=para.rho_Wfb(j,k);
			i++;
		}
	}

	json_parameters["Ncx"].push_back(para.Ncx);
	json_parameters["Nmf"].push_back(para.Nmf);
	json_parameters["Ngc"].push_back(para.Ngc);
	json_parameters["Nsyn"].push_back(para.Nsyn);
	json_parameters["Ncomp"].push_back(para.Ncomp);
	json_parameters["Ncomp_intf"].push_back(para.Ncomp_intf);
	json_parameters["Ncomp_noise"].push_back(para.Ncomp_noise);
	json_parameters["Ncomp_noise_RNN"].push_back(para.Ncomp_noise_RNN);
	json_parameters["FCcomp_tol"].push_back(para.FCcomp_tol);
	json_parameters["Nreal"].push_back(para.Nreal);
	json_parameters["N_FMperiods"].push_back(para.N_FMperiods);
	json_parameters["Npara"].push_back(para.Npara);
	json_parameters["RNseed"].push_back(para.seed);
	json_parameters["Ampmod"].push_back(para.ampmod);
	json_parameters["alpha"].push_back(para.alpha);
	json_parameters["amplitudes"].push_back(para.amps);
	json_parameters["amplitudes_intf"].push_back(para.amps_intf);
	json_parameters["BWfilter_order"].push_back(para.norder);
	json_parameters["BWfilter_cutoff_factor"].push_back(para.cutoff);
	json_parameters["DT"].push_back(para.DT);
	json_parameters["DT_euler"].push_back(para.DT_euler);
	json_parameters["dt"].push_back(para.dt);
	json_parameters["dt_euler"].push_back(para.dt_euler);
	json_parameters["delta_eta"].push_back(para.delta_eta);
	json_parameters["delta_r2_reg"].push_back(para.delta_r2_reg);
	json_parameters["FACsigmaLS_ctx"].push_back(para.FACsigmaLS_ctx);
	json_parameters["FACsigmaLS_gc"].push_back(para.FACsigmaLS_gc);
	json_parameters["fc_avrg"].push_back(para.fc_avrg);
	json_parameters["fc_sd"].push_back(para.fc_SD);
	json_parameters["g"].push_back(para.g);
	json_parameters["GCgain"].push_back(para.GCgain);
	json_parameters["GCgainmax"].push_back(para.GCgainmax);
	json_parameters["GCmax"].push_back(para.GCmax);
	json_parameters["tol_GCgain"].push_back(para.tol_GCgain);
	json_parameters["GCthreshold_base"].push_back(para.GCthreshold_base);
	json_parameters["GCthreshold"].push_back(para.GCthreshold*para.sqNsyn);
	json_parameters["INT_delta_abs"].push_back(para.INT_delta_abs);
	json_parameters["INT_delta_rel"].push_back(para.INT_delta_rel);
	json_parameters["lambda"].push_back(para.lambda);
	json_parameters["pcon"].push_back(para.pcon);
	json_parameters["phases"].push_back(para.phases);
	json_parameters["phases_intf"].push_back(para.phases_intf);
	json_parameters["Phi1"].push_back(para.phi1);
	json_parameters["rho_Wfb12"].push_back(para.rho_Wfb12);
	json_parameters["rho_Wfb13"].push_back(para.rho_Wfb13);
	json_parameters["rho_Wfb"].push_back(rho_Wfb_vec);
	json_parameters["sigAmp"].push_back(para.sigAmp);
	json_parameters["sigAmp_intf"].push_back(para.sigAmp_intf);
	json_parameters["sigmaLS"].push_back(para.sigmaLS);
	json_parameters["sigmaN"].push_back(para.sigmaN);
	json_parameters["sigmaNgc"].push_back(para.sigmaNgc);
	json_parameters["sp_alpha"].push_back(para.sp_alpha);
	json_parameters["tau_avrg"].push_back(para.tau_avrg);
	json_parameters["tau_std"].push_back(para.tau_std);
	json_parameters["tauRNNcutoff"].push_back(para.tauRNNcutoff);
	json_parameters["tauGC_avrg"].push_back(para.tauGC_avrg);
	json_parameters["tauGC_std"].push_back(para.tauGC_std);
	json_parameters["tauGCcutoff"].push_back(para.tauGCcutoff);
	json_parameters["Tlearn"].push_back(para.duration_learn/para.lambda);
	json_parameters["Ttest"].push_back(para.duration_test/para.lambda);
	json_parameters["Tft"].push_back(para.duration_FT/para.lambda);
	json_parameters["Ttrans"].push_back(para.duration_trans/para.lambda);
	json_parameters["Wgoc"].push_back(para.Wgoc);
	json_parameters["Wgoc_SD"].push_back(para.Wgoc_SD);


	// open file
	std::ofstream param_file(para.path+"parameters.dat");		if (!param_file) exit(EXIT_FAILURE);
	param_file << std::setw(4) << json_parameters << '\n';
	param_file.close();

	return;
}
