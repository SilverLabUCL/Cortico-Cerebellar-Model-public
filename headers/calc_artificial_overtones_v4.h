///////////////////////////////////
// The complement function
///////////////////////////////////
//~ arma::vec complement(arma::vec& y, double rho) {
	//~ arma::vec x(y.n_elem, arma::fill::randn);

	//~ arma::vec y_perp = x - y * arma::solve(y, x);
	//~ double sd_y_perp = arma::stddev(y_perp);
	//~ double sd_y = arma::stddev(y);

	//~ arma::vec result(y.n_elem);
	//~ for (size_t i = 0; i < y.n_elem; i) {
		//~ result[i] = rho * sd_y_perp * y[i] + y_perp[i] * sd_y * std::sqrt(1 - rho * rho);
	//~ }

	//~ return result;
//~ }
///////////////////////////////////
// The complement function for multidimensional y
///////////////////////////////////
arma::vec complement_multiD(arma::mat y, const arma::vec rho) {
	//	y - Matrix of variables
	//	rho - Desired correlation vector
	//	threshold - Small value to handle numerical stability
	double threshold = 1e-10;
	const arma::uword N = y.n_rows;
	//~ y = arma::normalise(y); 	// Center and scale y
	for (size_t i = 0; i < y.n_cols; i++) {
		//~ y.col(i) = y.col(i)- arma::mean(y.col(i));        // centre
		y.col(i) = y.col(i)/arma::stddev(y.col(i));
	}

	arma::vec x(N, arma::fill::randn);

	// Residuals from regression of x on y
	//~ arma::vec y_perp = x - y * arma::solve(y, x);

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

	// With zero-mean y:
	//~ arma::mat second_moment = (y_dual.t() * y_dual) / double(N - 1);
	//~ double var_ydual_rho    = arma::as_scalar(rho.t() * second_moment * rho);
	//~ double sigma2           = (1.0 - var_ydual_rho) / arma::var(y_perp);

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

/*arma::mat nearestCorrPSD_simple(const arma::mat& R_original, const arma::umat& mask, int maxIter = 5000, double tol = 1e-10){

	const arma::uvec mask_idx = arma::find(mask);

	// Symmetrise input once before the loop
	arma::mat R_adjusted = 0.5 * (R_original + R_original.t());
	R_adjusted.diag().ones();

	// One Dykstra increment per projection:
	//   Delta_S  : increment for the PSD cone projection
	//   Delta_C  : increment for the constraint-set projection (mask + symmetry + diagonal)
	arma::mat Delta_S = arma::zeros(arma::size(R_adjusted));
	arma::mat Delta_C = arma::zeros(arma::size(R_adjusted));

	arma::vec eigval;
	arma::mat eigvec;
	bool converged = false;

	for (int iter = 0; iter < maxIter; iter++) {
		// Projection 1: onto the PSD cone
		arma::mat Y = R_adjusted + Delta_S;			// Dykstra-corrected input
		arma::eig_sym(eigval, eigvec, Y);
		eigval = arma::clamp(eigval, 0.0, arma::datum::inf);
		arma::mat R_adjusted_psd = eigvec * arma::diagmat(eigval) * eigvec.t();
		Delta_S = Y - R_adjusted_psd;				// update increment for PSD set

		// Projection 2: onto the constraint set
		arma::mat Z = R_adjusted_psd + Delta_C;		// Dykstra-corrected input
		arma::mat R_adjusted_new = Z;
		R_adjusted_new.elem(mask_idx) = R_original.elem(mask_idx);	// enforce fixed entries
		R_adjusted_new.diag().ones();								// enforce unit diagonal
		R_adjusted_new = 0.5 * (R_adjusted_new + R_adjusted_new.t());					// enforce symmetry
		Delta_C = Z - R_adjusted_new;				// update increment for constraint set

		R_adjusted = R_adjusted_new;

		if (eigval.min() >= -tol) {
			converged = true;
			break;
		}
	}

	if (!converged)	{
		std::cerr << "Warning: nearestCorrPSD did not converge within " << maxIter <<" iterations." << std::endl;
	}

	return R_adjusted;
}*/

void calc_artificial_overtones(FTcomp &FT_TRANS_comp_ctx, netpara &para, int np, nlohmann::json &json_measures){

	para.Intf_weights=arma::zeros<arma::mat>(para.Ncx,para.Ncomp_intf);
	// one overtone with relationship between overtone and fundamental sine components free
	if(para.Ncomp_intf==2){

		arma::vec rho_tmp={para.rho_Wfb12};
		//~ arma::vec rho_tmp={R_adjusted(0,2),R_adjusted(1,2)};
		arma::mat spann_vec_RNN(para.Ncx,rho_tmp.n_elem);
		for(int i= 0; i < para.Ncx; i++){
			for(size_t j= 0; j < rho_tmp.n_elem; j++){
				spann_vec_RNN(i,j)=FT_TRANS_comp_ctx.S_total_sqN(i,j+1);		// 1st cosine spanning vector
				//~ spann_vec_RNN(i,1)=FT_TRANS_comp_ctx.S_total_sqN(i,2);		// 1st sine spanning vector
			}
		}

		// generate interference weights
		arma::vec Intf_vec_tmp = complement_multiD(spann_vec_RNN, rho_tmp);
		//~ std::cout << "corr SVs: " << std::endl;
		//~ std::cout << arma::cor( arma::join_horiz( spann_vec_RNN, Intf_vec_tmp) ) << std::endl;

		for(int i= 0; i < para.Ncx; i++){
				para.Intf_weights(i,1)=Intf_vec_tmp[i];
		}
	}
	//  one overtone with relationship between overtone and fundamental sine components fixed
	/*if(para.Ncomp_intf==2){
		// extract existing correlations
		arma::mat spann_vec_RNN(para.Ncx,2);
		for(int i= 0; i < para.Ncx; i++){
			for(int j= 0; j < 2; j++){
				spann_vec_RNN(i,j)=FT_TRANS_comp_ctx.S_total_sqN(i,j+1);		// j=1: 1st cosine spanning vector, j=2: 1st sine spanning vector
			}
		}
		// compose existing corr-matrix
		arma::mat corrRNN=arma::cor( spann_vec_RNN );
		//~ std::cout << "corr RNN SVs: " << std::endl;
		//~ std::cout << corrRNN << std::endl;
		arma::mat R_existing={{1.0,corrRNN(0,1)},{corrRNN(0,1),1.0}};

		// specify cross- and new correlation matrices
		arma::mat R_cross={para.rho_Wfb12,0.2};
		arma::mat R_new={1.0};
		//~ std::cout << "R_cross:\n" << R_cross << std::endl
		//~ std::cout << "R_new:\n" << R_new << std::endl

		arma::mat AB = arma::join_horiz( R_existing, R_cross.t());
		arma::mat CD = arma::join_horiz( R_cross, R_new);
		arma::mat R_full=arma::join_vert( AB, CD);

		//~ std::cout << "R_full:\n" << R_full << std::endl

		// create mask: 1 for values that must remain fixed, 0 for values that can be optimized
		int dim = R_full.n_rows;
		arma::mat mask = arma::zeros<arma::mat>(dim,dim);
		for (int i = 0; i < dim; i++) {
			mask(i,i)=1;
		}
		for (int i = 0; i < R_existing.n_rows; i++) {
			for (int j = 0; j < R_existing.n_rows; j++) {
				mask(i,j)=1;
			}
		}
		mask(0,2)=1;	// fix correlation between cosine and 1st overtone cosine (para.rho_Wfb12)
		mask(2,0)=1;
		//~ std::cout << "Mask:\n" << mask << std::endl

		// Adjust to nearest PSD correlation matrix
		arma::mat R_adjusted = nearestCorrPSD(R_full, R_existing, mask);
		//~ std::cout << "R_adjusted:\n" << R_adjusted << std::endl

		//~ arma::vec rho_tmp={para.rho_Wfb12};
		arma::vec rho_tmp={R_adjusted(0,2),R_adjusted(1,2)};
		// generate interference weights
		arma::vec Intf_vec_tmp = complement_multiD(spann_vec_RNN, rho_tmp);
		//~ std::cout << "corr SVs: " << std::endl;
		//~ std::cout << arma::cor( arma::join_horiz( spann_vec_RNN, Intf_vec_tmp) ) << std::endl;

		for(int i= 0; i < para.Ncx; i++){
				para.Intf_weights(i,1)=Intf_vec_tmp[i];
		}
	}*/

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
		//~ std::cout << "corr RNN SVs: " << std::endl;
		//~ std::cout << corrRNN << std::endl;
		arma::mat R_existing={{1.0,corrRNN(0,1)},{corrRNN(0,1),1.0}};

		// specify cross- and new correlation matrices
		arma::mat R_cross={{para.rho_Wfb12,para.rho_Wfb13},{para.rho_Wfb(1,2),para.rho_Wfb(1,3)}};
		arma::mat R_new={{1.0,para.rho_Wfb(2,3)},{para.rho_Wfb(2,3),1.0}};

		arma::mat AB = arma::join_horiz( R_existing, R_cross);
		arma::mat CD = arma::join_horiz( R_cross.t(), R_new);
		arma::mat R_full=arma::join_vert( AB, CD);
		//~ std::cout << "R_full:\n" << R_full << std::endl

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
		//~ std::cout << "Mask:\n" << mask << std::endl

		// Adjust to nearest PSD correlation matrix
		arma::mat R_adjusted = nearestCorrPSD_simple(R_full, mask);
		//~ std::cout << "R_adjusted:\n" << R_adjusted << std::endl

		json_measures["R_cereb"].push_back(R_adjusted);

		//~ arma::vec rho_tmp={para.rho_Wfb12};
		arma::vec rho_tmp={R_adjusted(0,2),R_adjusted(1,2)};

		// generate interference weights
		arma::vec Intf_vec_tmp = complement_multiD(spann_vec_RNN, rho_tmp);
		//~ std::cout << "corr SVs: " << std::endl;
		//~ std::cout << arma::cor( spann_vec_RNN ) << std::endl;
		//~ std::cout << "corr SVs-Overtones: " << std::endl;
		//~ std::cout << arma::cor( spann_vec_RNN, Intf_vec_tmp ) << std::endl;

		for(int i= 0; i < para.Ncx; i++){
			for(int j= 1; j < 2; j++){
				// Intf_vec[i]=Intf_vec_tmp[i];
				//~ std::cout <<"additional SV: "<< j << std::endl;
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
		//~ std::cout << "corr SVs: " << std::endl;
		//~ std::cout << arma::cor( arma::join_horiz( spann_vec_RNN_intf1, Intf_vec_tmp2) ) << std::endl;

		for(int i= 0; i < para.Ncx; i++){
			for(int j= 2; j < para.Ncomp_intf; j++){
				//~ std::cout <<"additional SV: "<< j << std::endl;
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
				//~ spann_vec_RNN(i,j)=FT_TRANS_comp_ctx.S_total_sqN(i,j+1);		// j=1: 1st cosine spanning vector, j=2: 1st sine spanning vector
				spann_vec_RNN(i,j)=FT_TRANS_comp_ctx.S_signal_sqN(i,j);
			}
		}
		// compose existing corr-matrix
		arma::mat corrRNN=arma::cor( spann_vec_RNN );
		//~ std::cout << "corr RNN SVs: " << std::endl;
		//~ std::cout << corrRNN << std::endl;
		arma::mat R_existing=corrRNN;
		//~ std::cout << "R_existing:\n" << R_existing << std::endl;

		// specify cross- and new correlation matrices
		arma::mat R_cross={{para.rho_Wfb12,para.rho_Wfb13},{para.rho_Wfb(1,4),para.rho_Wfb(1,5)},{para.rho_Wfb(2,4),para.rho_Wfb(2,5)},{para.rho_Wfb(3,4),para.rho_Wfb(3,5)}};
		arma::mat R_new={{1.0,para.rho_Wfb(4,5)},{para.rho_Wfb(5,4),1.0}};
		//~ std::cout << "R_cross:\n" << R_cross << std::endl;
		//~ std::cout << "R_new:\n" << R_new << std::endl;
		arma::mat AB = arma::join_horiz( R_existing, R_cross);
		arma::mat CD = arma::join_horiz( R_cross.t(), R_new);
		arma::mat R_full=arma::join_vert( AB, CD);
		//~ std::cout << "R_full:\n" << R_full << std::endl;

		//~ arma::vec eigenvalues_full = arma::eig_sym( R_full );
		//~ std::cout << "R_full EVs = " << eigenvalues_full << std::endl;

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
		//~ std::cout << "Mask:\n" << mask << std::endl

		// Adjust to nearest PSD correlation matrix
		arma::mat R_adjusted = nearestCorrPSD_simple(R_full, mask);
		//~ std::cout << "R_adjusted:\n" << R_adjusted << std::endl;
		arma::vec eigenvalues_adj = arma::eig_sym( R_adjusted );

		//~ std::cout << "R_adjusted = " << std::endl;
		//~ std::cout << R_adjusted << std::endl;
		//~ std::cout << "min(EV) = " << eigenvalues_adj.min() << std::endl;
		//~ std::cout << "nreadj = " << nreadj << std::endl;

		json_measures["R_cereb"].push_back(R_adjusted);

		arma::vec rho_tmp={R_adjusted(0,4),R_adjusted(1,4),R_adjusted(2,4),R_adjusted(3,4)};

		// generate interference weights
		arma::vec Intf_vec_tmp = complement_multiD(spann_vec_RNN, rho_tmp);
		//~ std::cout << "corr SVs: " << std::endl;
		//~ std::cout << arma::cor( spann_vec_RNN ) << std::endl;
		//~ std::cout << arma::cor( arma::join_horiz( spann_vec_RNN, Intf_vec_tmp) ) << std::endl;
		//~ std::cout << "corr SVs-Overtones: " << std::endl;
		//~ std::cout << arma::cor( spann_vec_RNN, Intf_vec_tmp ) << std::endl;

		for(int i= 0; i < para.Ncx; i++){
			for(int j= 2; j < 3; j++){
				para.Intf_weights(i,j)=Intf_vec_tmp[i];
			}
		}

		arma::vec rho_tmp2={R_adjusted(0,5),R_adjusted(1,5),R_adjusted(2,5),R_adjusted(3,5),R_adjusted(4,5)};
		arma::mat spann_vec_RNN_intf1(para.Ncx,rho_tmp2.n_elem);
		for(int i= 0; i < para.Ncx; i++){
			for(int j= 0; j < 2*para.Ncomp; j++){
				//~ spann_vec_RNN_intf1(i,j)=FT_TRANS_comp_ctx.S_total_sqN(i,j+1);		// RNN spanning vectors
				spann_vec_RNN_intf1(i,j)=FT_TRANS_comp_ctx.S_signal_sqN(i,j);
			}
			spann_vec_RNN_intf1(i,4)=Intf_vec_tmp[i];								// 1st overtone spanning vector
		}

		//~ std::cout << arma::size(para.Intf_weights) << "\t" << arma::size(Intf_vec_tmp) << "\t" << arma::size(spann_vec_RNN_intf1) << std::endl
		//~ std::cout << "corr spann_vec_RNN_intf1: " << std::endl;
		//~ std::cout << arma::cor( spann_vec_RNN_intf1 ) << std::endl;

		arma::vec Intf_vec_tmp2 = complement_multiD(spann_vec_RNN_intf1, rho_tmp2);
		//~ std::cout << "corr SVs: " << std::endl;
		//~ std::cout << arma::cor( arma::join_horiz( spann_vec_RNN_intf1, Intf_vec_tmp2) ) << std::endl;
		//~ std::cout << "rho_Wfb:\n" << para.rho_Wfb << std::endl
		//~ std::cout << "R_adjusted:\n" << R_adjusted << std::endl;

		for(int i= 0; i < para.Ncx; i++){
			for(int j= 3; j < para.Ncomp_intf; j++){
				para.Intf_weights(i,j)=Intf_vec_tmp2[i];
			}
		}
	}

}

