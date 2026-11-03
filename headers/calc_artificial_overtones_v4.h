///////////////////////////////////
// The complement function for multidimensional y
///////////////////////////////////
arma::vec complement_multiD(arma::mat y, const arma::vec rho) {
	//	y - Matrix of variables
	//	rho - Desired correlation vector
	//	threshold - Small value to handle numerical stability
	double threshold = 1e-10;
	const arma::uword N = y.n_rows;
	for (size_t i = 0; i < y.n_cols; i++) {
		y.col(i) = y.col(i)/arma::stddev(y.col(i));
	}

	arma::vec x(N, arma::fill::randn);

	// Calculate y.dual using SVD of y
	arma::mat U; arma::mat V;
	arma::vec singular_values;
	arma::svd_econ( U, singular_values, V, y );
	arma::vec dual_diag = singular_values;
	for (size_t i = 0; i < dual_diag.n_elem; i++) {
		if(singular_values(i) <= threshold)		dual_diag(i)=0.0;
		else									dual_diag(i)=1.0 / singular_values(i);
	}
	arma::mat y_dual = double(N - 1) * U * arma::diagmat(dual_diag) * V.t();

	// Residuals: use U from SVD
	arma::vec y_perp = x - U * (U.t() * x);

	// Calculate sigma^2
	double sigma2 = arma::as_scalar(1.0 - rho.t() * arma::cov(y_dual) * rho )/arma::var(y_perp);

	// Return linear combination if sigma^2 >= 0
	arma::vec z(N, arma::fill::zeros);
	if (sigma2 >= 0){
		z = y_dual * rho + sqrt(sigma2) * y_perp;
	}
	else{
		std::cerr << "Warning: Correlations for generating overtones are impossible, sigma2 = " << sigma2 << std::endl;
		z = y_dual * rho ;
	}
	return z;
}

///////////////////////////////////
// Function to compute the nearest PSD correlation matrix with designated fixed entries
///////////////////////////////////
arma::mat nearestCorrPSD_simple(const arma::mat& R_original, const arma::umat& mask, int maxIter = 1000, double tol = 1e-8){
	arma::mat R_adjusted = 0.5 * (R_original + R_original.t());

	arma::mat Delta = arma::zeros(size(R_adjusted));

	arma::uvec mask_idx = find(mask);

	bool converged = false;
	arma::vec eigval;
	arma::mat eigvec;
	for(size_t iter = 0; iter < maxIter; iter++)	{
		// PSD projection with Dykstra correction
		arma::mat Y = R_adjusted - Delta;

		arma::eig_sym(eigval, eigvec, Y);

		// remove negative EVs
		eigval = arma::clamp(eigval, 0.0, arma::datum::inf);
		arma::mat R = eigvec * arma::diagmat(eigval) * eigvec.t();

		Delta = R - Y;

		// enforce fixed entries
		R_adjusted = R;
		R_adjusted.elem(mask_idx) = R_original.elem(mask_idx);

		// symmetrise
		R_adjusted = 0.5 * (R_adjusted + R_adjusted.t());

		// enforce ones on diagonal
		R_adjusted.diag().ones();

		// convergence test
		eigval = arma::eig_sym( R_adjusted );
		if (eigval.min() >= -tol) {
			converged = true;
			break;
		}
	}
	if (!converged)	{
		std::cerr << "Warning: nearestCorrPSD did not converge within " << maxIter <<" iterations." << std::endl;
	}

	return R_adjusted;
}

void calc_artificial_overtones(FTcomp &FT_TRANS_comp_ctx, netpara &para, int np, nlohmann::json &json_measures){

	para.Intf_weights=arma::zeros<arma::mat>(para.Ncx,para.Ncomp_intf);
	// one overtone with relationship between overtone and fundamental sine components free
	if(para.Ncomp_intf==2){

		arma::vec rho_tmp={para.rho_Wfb12};
		arma::mat spann_vec_RNN(para.Ncx,rho_tmp.n_elem);
		for(int i= 0; i < para.Ncx; i++){
			for(size_t j= 0; j < rho_tmp.n_elem; j++){
				spann_vec_RNN(i,j)=FT_TRANS_comp_ctx.S_total_sqN(i,j+1);		// 1st cosine spanning vector
			}
		}

		// generate interference weights
		arma::vec Intf_vec_tmp = complement_multiD(spann_vec_RNN, rho_tmp);

		for(int i= 0; i < para.Ncx; i++){
				para.Intf_weights(i,1)=Intf_vec_tmp[i];
		}
	}

	// one signal frequency, two overtones
	if(para.Ncomp_intf==3 && para.Ncomp==1){
		// extract existing correlations
		arma::mat spann_vec_RNN(para.Ncx,2*para.Ncomp);
		for(int i= 0; i < para.Ncx; i++){
			for(int j= 0; j < 2*para.Ncomp; j++){
				spann_vec_RNN(i,j)=FT_TRANS_comp_ctx.S_total_sqN(i,j+1);		// j=1: 1st cosine spanning vector, j=2: 1st sine spanning vector
			}
		}
		// compose existing corr-matrix
		arma::mat corrRNN=arma::cor( spann_vec_RNN );
		arma::mat R_existing={{1.0,corrRNN(0,1)},{corrRNN(0,1),1.0}};

		// specify cross- and new correlation matrices
		arma::mat R_cross={{para.rho_Wfb12,para.rho_Wfb13},{para.rho_Wfb(1,2),para.rho_Wfb(1,3)}};
		arma::mat R_new={{1.0,para.rho_Wfb(2,3)},{para.rho_Wfb(2,3),1.0}};

		arma::mat AB = arma::join_horiz( R_existing, R_cross);
		arma::mat CD = arma::join_horiz( R_cross.t(), R_new);
		arma::mat R_full=arma::join_vert( AB, CD);

		// create mask: 1 for values that must remain fixed, 0 for values that can be optimized
		int dim = R_full.n_rows;
		arma::umat mask = arma::zeros<arma::umat>(dim,dim);
		for (int i = 0; i < dim; i++) {
			mask(i,i)=1;
		}
		for (size_t i = 0; i < R_existing.n_rows; i++) {
			for (size_t j = 0; j < R_existing.n_rows; j++) {
				mask(i,j)=1;
			}
		}
		mask(0,2)=1;	// fix correlation between cosine and 1st overtone cosine (para.rho_Wfb12)
		mask(2,0)=1;
		mask(0,3)=1;	// fix correlation between cosine and 2st overtone cosine (para.rho_Wfb13)
		mask(3,0)=1;

		// Adjust to nearest PSD correlation matrix
		arma::mat R_adjusted = nearestCorrPSD_simple(R_full, mask);

		json_measures["R_cereb"].push_back(R_adjusted);

		arma::vec rho_tmp={R_adjusted(0,2),R_adjusted(1,2)};

		// generate interference weights
		arma::vec Intf_vec_tmp = complement_multiD(spann_vec_RNN, rho_tmp);

		for(int i= 0; i < para.Ncx; i++){
			for(int j= 1; j < 2; j++){
				para.Intf_weights(i,j)=Intf_vec_tmp[i];
			}
		}

		//~ arma::vec rho_tmp2={para.rho_Wfb12*0.5};
		arma::vec rho_tmp2={R_adjusted(0,3),R_adjusted(1,3),R_adjusted(2,3)};
		arma::mat spann_vec_RNN_intf1(para.Ncx,rho_tmp2.n_elem);
		for(int i= 0; i < para.Ncx; i++){
			for(size_t j= 0; j < rho_tmp2.n_elem-1; j++){
				spann_vec_RNN_intf1(i,j)=FT_TRANS_comp_ctx.S_total_sqN(i,j+1);		// RNN spanning vectors
			}
			spann_vec_RNN_intf1(i,2)=Intf_vec_tmp(i);								// 1st overtone spanning vector
		}
		arma::vec Intf_vec_tmp2 = complement_multiD(spann_vec_RNN_intf1, rho_tmp2);

		for(int i= 0; i < para.Ncx; i++){
			for(int j= 2; j < para.Ncomp_intf; j++){
				para.Intf_weights(i,j)=Intf_vec_tmp2[i];
			}
		}
	}

	// two signal frequencies, two overtones
	if(para.Ncomp_intf==4 && para.Ncomp==2){
		// extract existing correlations
		arma::mat spann_vec_RNN(para.Ncx,2*para.Ncomp);
		for(int i= 0; i < para.Ncx; i++){
			for(int j= 0; j < 2*para.Ncomp; j++){
				spann_vec_RNN(i,j)=FT_TRANS_comp_ctx.S_signal_sqN(i,j);
			}
		}
		// compose existing corr-matrix
		arma::mat corrRNN=arma::cor( spann_vec_RNN );
		arma::mat R_existing=corrRNN;

		// specify cross- and new correlation matrices
		arma::mat R_cross={{para.rho_Wfb12,para.rho_Wfb13},{para.rho_Wfb(1,4),para.rho_Wfb(1,5)},{para.rho_Wfb(2,4),para.rho_Wfb(2,5)},{para.rho_Wfb(3,4),para.rho_Wfb(3,5)}};
		arma::mat R_new={{1.0,para.rho_Wfb(4,5)},{para.rho_Wfb(5,4),1.0}};
		arma::mat AB = arma::join_horiz( R_existing, R_cross);
		arma::mat CD = arma::join_horiz( R_cross.t(), R_new);
		arma::mat R_full=arma::join_vert( AB, CD);

		// create mask: 1 for values that must remain fixed, 0 for values that can be optimized
		int dim = R_full.n_rows;
		arma::umat mask = arma::zeros<arma::umat>(dim,dim);
		for (int i = 0; i < dim; i++) {
			mask(i,i)=1;
		}
		for (size_t i = 0; i < R_existing.n_rows; i++) {
			for (size_t j = 0; j < R_existing.n_rows; j++) {
				mask(i,j)=1;
			}
		}
		mask(0,4)=1;	// fix correlation between cosine and 1st overtone cosine (para.rho_Wfb12)
		mask(4,0)=1;
		mask(0,5)=1;	// fix correlation between cosine and 2st overtone cosine (para.rho_Wfb13)
		mask(5,0)=1;

		// Adjust to nearest PSD correlation matrix
		arma::mat R_adjusted = nearestCorrPSD_simple(R_full, mask);
		arma::vec eigenvalues_adj = arma::eig_sym( R_adjusted );

		json_measures["R_cereb"].push_back(R_adjusted);

		arma::vec rho_tmp={R_adjusted(0,4),R_adjusted(1,4),R_adjusted(2,4),R_adjusted(3,4)};

		// generate interference weights
		arma::vec Intf_vec_tmp = complement_multiD(spann_vec_RNN, rho_tmp);

		for(int i= 0; i < para.Ncx; i++){
			for(int j= 2; j < 3; j++){
				para.Intf_weights(i,j)=Intf_vec_tmp[i];
			}
		}

		arma::vec rho_tmp2={R_adjusted(0,5),R_adjusted(1,5),R_adjusted(2,5),R_adjusted(3,5),R_adjusted(4,5)};
		arma::mat spann_vec_RNN_intf1(para.Ncx,rho_tmp2.n_elem);
		for(int i= 0; i < para.Ncx; i++){
			for(int j= 0; j < 2*para.Ncomp; j++){
				spann_vec_RNN_intf1(i,j)=FT_TRANS_comp_ctx.S_signal_sqN(i,j);
			}
			spann_vec_RNN_intf1(i,4)=Intf_vec_tmp[i];								// 1st overtone spanning vector
		}

		arma::vec Intf_vec_tmp2 = complement_multiD(spann_vec_RNN_intf1, rho_tmp2);

		for(int i= 0; i < para.Ncx; i++){
			for(int j= 3; j < para.Ncomp_intf; j++){
				para.Intf_weights(i,j)=Intf_vec_tmp2[i];
			}
		}
	}

}

