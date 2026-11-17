void write_measures_final(nlohmann::json &json_measures, FTcomp &Cmp, double RMSE_train, double RMSE_test, double Wnorm, std::string mode, netpara &para, bool write_FULL) {

	json_measures["RMSE_train_" + mode].push_back(RMSE_train);
	json_measures["RMSE_test_" + mode].push_back(RMSE_test);

	json_measures["Wnorm_" + mode].push_back(Wnorm);

	nlohmann::json json_tmp;

	json_tmp.clear();
	for (size_t i= 0; i < Cmp.subangles.n_elem; i++) {
		json_tmp.push_back(Cmp.subangles.at(i));
	}
	json_measures["Subspaceangles_cos_" + mode].push_back(json_tmp);

	json_measures["Sineprod_signal_" + mode].push_back(Cmp.sineprod_signal);

	//~ std::string istr,jstr,numistr,numjstr;
	if(write_FULL==1){

		json_tmp.clear();
		for (int i= 0; i < (Cmp.Ms+para.Nwrite); i++) {
			for (int j= 0; j <= i; j++) {
				json_tmp.push_back(Cmp.Angles_all.at(i,j));
			}
		}
		json_measures["Angles_" + mode].push_back(json_tmp);

		json_tmp.clear();
		for (int i= 0; i < (Cmp.Ms+para.Nwrite); i++) {
			json_tmp.push_back(Cmp.nn_total.at(i));
		}
		json_measures["Norms_" + mode].push_back(json_tmp);

	}
	else {
		json_tmp.clear();
		for (int i= 0; i < para.Ncomp; i++) {
			json_tmp.push_back(Cmp.nn_total.at(2*i+1));
		}
		json_measures["Norms_signal_" + mode].push_back(json_tmp);
	}

	// open file
	std::ofstream measures_file(para.folder+"measures.json");		if (!measures_file) exit(EXIT_FAILURE);
	measures_file << std::setw(4) << json_measures << '\n';
	measures_file.close();

	return;
}
