///////////////////////////////////
// The complement function
///////////////////////////////////
//~ arma::vec complement(arma::vec& y, double rho) {
	//~ arma::vec x(y.n_elem, arma::fill::randn);

	//~ arma::vec y_perp = x - y * arma::solve(y, x);
	//~ double sd_y_perp = arma::stddev(y_perp);
	//~ double sd_y = arma::stddev(y);

	//~ arma::vec result(y.n_elem);
	//~ for (size_t i = 0; i < y.n_elem; ++i) {
		//~ result[i] = rho * sd_y_perp * y[i] + y_perp[i] * sd_y * std::sqrt(1 - rho * rho);
	//~ }

	//~ return result;
//~ }
///////////////////////////////////
// The complement function for multidimensional y
///////////////////////////////////
arma::vec complement_multiD(arma::mat y, arma::vec rho) {
	//	y - Matrix of variables
	//	rho - Desired correlation vector
	//	threshold - Small value to handle numerical stability
	double threshold = 1e-12;
	int N = y.n_rows;
	//~ int d = y.n_cols;
	//~ y = arma::normalise(y); 	// Center and scale y
	for (size_t i = 0; i < y.n_cols; ++i) {
		y.col(i) = y.col(i)/arma::stddev(y.col(i));
	}
	//~ std::cout << "variance of y " << arma::var( y ) << std::endl;

	arma::vec x(N, arma::fill::randn);

	// Residuals from regression of x on y
	arma::vec y_perp = x - y * arma::solve(y, x);

	// Calculate y.dual using SVD of y
	arma::mat U; arma::mat V;
	arma::vec singular_values;
	arma::svd_econ( U, singular_values, V, y );
	arma::vec dual_diag = singular_values;
	for (size_t i = 0; i < dual_diag.n_elem; ++i) {
		if(singular_values(i) <= threshold)		dual_diag(i)=0.0;
		else									dual_diag(i)=1.0 / singular_values(i);
	}
	arma::mat y_dual = double(N - 1) * U * arma::diagmat(dual_diag) * V.t();

	// Calculate sigma^2
	double sigma2 = arma::as_scalar(1.0 - rho.t() * arma::cov(y_dual) * rho )/arma::var(y_perp);
	//~ std::cout << "dual_diag = " << dual_diag << std::endl;
	//~ std::cout << "singular values = " << singular_values << std::endl;
	//~ std::cout << "cov(y_dual) = " << arma::cov(y_dual) << std::endl;
	//~ std::cout << "sigma2 = " << sigma2 << "\t var_y_perp = " << arma::var(y_perp) << std::endl;

	// Return linear combination if sigma^2 >= 0
	arma::vec z(N, arma::fill::zeros);
	if (sigma2 >= 0){
		z = y_dual * rho + sqrt(sigma2) * y_perp;
	}
	else{
		std::cerr << "Warning: Correlations for generating overtones are impossible, sigma2 = " << sigma2 << std::endl;
		z = y_dual * rho ;
		//~ std::cout << "rho = " << std::endl;
		//~ std::cout << rho << std::endl;
	}
	return z;
}

///////////////////////////////////
// Struct to pass additional data to optimization functions
///////////////////////////////////
struct OptimizationData {
	arma::mat R_full;
	arma::mat mask;
};

///////////////////////////////////
// Compute Frobenius norm
///////////////////////////////////
double objectiveFunction(const std::vector<double>& X, std::vector<double>& grad, void* data) {
	auto* params = static_cast<OptimizationData*>(data);
	arma::mat R_full = params->R_full;
	int dim = R_full.n_rows;

	// Map flattened vector back to a matrix
	arma::vec R_flattened(dim*dim);
	for (int i = 0; i < dim*dim; ++i) R_flattened[i] = X[i];
	arma::mat X_matrix = arma::reshape( R_flattened, dim, dim );

	// Compute Frobenius norm difference
	return arma::norm(X_matrix-R_full, "fro")*arma::norm(X_matrix-R_full, "fro");
	//~ return arma::norm(X_matrix-R_full, "fro");
}
///////////////////////////////////
// Ensure eigenvalues are non-negative (PSD constraint)
///////////////////////////////////
void psdConstraint(unsigned m, double* result, unsigned n, const double* X, double* grad, void* data) {
	auto* params = static_cast<OptimizationData*>(data);
	int dim = params->R_full.n_rows;

	// Map flattened vector back to a matrix
	arma::vec R_flattened(dim*dim);
	for (int i = 0; i < dim*dim; ++i) R_flattened[i] = X[i];
	arma::mat X_matrix = arma::reshape( R_flattened, dim, dim );
	//~ std::cout << "X_matrix = " << X_matrix << std::endl;
	//~ arma::mat X_matrix(const_cast<double*>(X), dim, dim, false, true);


	// Symmetrize
	X_matrix = (X_matrix + X_matrix.t()) / 2.0;

	// Compute eigenvalues
	arma::vec eigenvalues = arma::eig_sym( X_matrix );

	// Set constraints (eigenvalues must be >= 0)
	for (int i = 0; i < dim; ++i) result[i] = -eigenvalues(i);
}
///////////////////////////////////
// Condition to ensure diagonal elements are 1
///////////////////////////////////
//~ void diagConstraint(unsigned m, double* result, unsigned n, const double* X, double* grad, void* data) {
	//~ auto* params = static_cast<OptimizationData*>(data);
	//~ int dim = params->R_full.n_rows;

	//~ // Map flattened vector back to a matrix
	//~ arma::vec R_flattened(dim*dim);
	//~ for (int i = 0; i < dim*dim; ++i) R_flattened[i] = X[i];
	//~ arma::mat X_matrix = arma::reshape( R_flattened, dim, dim );

	//~ for (int i = 0; i < dim; ++i) result[i] = X_matrix(i, i) - 1.0;
//~ }
///////////////////////////////////
// Condition to preserve specified existing correlations
///////////////////////////////////
//~ void preserveExistingConstraint(unsigned m, double* result, unsigned n, const double* X, double* grad, void* data) {
	//~ auto* params = static_cast<OptimizationData*>(data);
	//~ int dim = params->R_full.n_rows;
	//~ arma::mat mask = params->mask;

	//~ // Map flattened vector back to a matrix
	//~ arma::vec R_flattened(dim * dim);
	//~ for (int i = 0; i < dim * dim; ++i) R_flattened[i] = X[i];
	//~ arma::mat X_matrix = arma::reshape( R_flattened, dim, dim );

	//~ // LOOP OVER MASK: WHEN ENTRY=1, ENFORCE CONDITION
	//~ int index = 0;
	//~ for (int i = 0; i < dim; ++i) {
		//~ for (int j = 0; j < dim; ++j) {
			//~ if(mask(i,j)==1)	{
				//~ result[index] = X_matrix(i, j) - params->R_full(i, j);
				//~ index++;
			//~ }
		//~ }
	//~ }
//~ }
///////////////////////////////////
// Condition to simultaneously preserve specified existing correlations and diagonal
///////////////////////////////////
void DiagExistingConstraint(unsigned m, double* result, unsigned n, const double* X, double* grad, void* data) {
	auto* params = static_cast<OptimizationData*>(data);
	int dim = params->R_full.n_rows;
	arma::mat mask = params->mask;

	// Map flattened vector back to a matrix
	arma::vec R_flattened(dim * dim);
	for (int i = 0; i < dim * dim; ++i) R_flattened[i] = X[i];
	arma::mat X_matrix = arma::reshape( R_flattened, dim, dim );

	// LOOP OVER MASK: WHEN ENTRY=1, ENFORCE CONDITION
	int index = 0;
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			if(mask(i,j)==1)	{
				result[index] = X_matrix(i, j) - params->R_full(i, j);
				index++;
			}
		}
	}
	for (int i = 0; i < dim; ++i) {
		result[index] = X_matrix(i, i) - 1.0;
		index++;
	}

}
///////////////////////////////////
// Function to compute the nearest PSD correlation matrix
///////////////////////////////////
arma::mat nearestCorrPSD(arma::mat R_full, arma::mat& R_existing, arma::mat& mask, int nreadj) {
	int dim = R_full.n_rows;

	arma::mat R_result;
	arma::vec eigenvalues = arma::eig_sym( R_full );
	//~ std::cout << "EVs: " <<  eigenvalues << std::endl;

	// Symmetrize the result
	R_full = (R_full + R_full.t()) / 2.0;

	if( arma::min(eigenvalues) > 0){
		R_result=R_full;
	}
	else{
		// Flatten input matrix
		arma::vec R_flattened = arma::vectorise(R_full);
		std::vector<double> X_flattened(R_flattened.size());
		for (size_t i = 0; i < R_flattened.size(); ++i) X_flattened[i] = R_flattened[i];

		// count ones in mask
		arma::vec mask_flattened = arma::vectorise(mask);
		int mones = int(arma::sum(mask_flattened));
		//~ std::cout << mones << std::endl;

		// Optimization setup
		nlopt::opt opt(nlopt::LN_COBYLA, dim * dim);
		OptimizationData optData = {R_full, mask};

		opt.set_min_objective(objectiveFunction, &optData);
		opt.add_inequality_mconstraint(psdConstraint, &optData, std::vector<double>(dim, 0));
		//~ opt.add_equality_mconstraint(diagConstraint, &optData, std::vector<double>(dim, 0));
		//~ opt.add_equality_mconstraint(preserveExistingConstraint, &optData, std::vector<double>(mones, 0));
		opt.add_equality_mconstraint(DiagExistingConstraint, &optData, std::vector<double>(mones+dim, 0));

		// Set bounds and tolerances
		opt.set_xtol_rel(1e-5);
		//~ opt.set_xtol_rel(1e-8);
		//~ opt.set_xtol_abs(1e-5);
		opt.set_ftol_rel(1e-5);
		opt.set_maxeval(700);

		//~ opt.set_initial_step(1e-2); // or a scalar like 1e-2

		// Optimize
		double minf;
		try {
			nlopt::result result = opt.optimize(X_flattened, minf);
		} catch (const nlopt::roundoff_limited& e) {
			std::cerr << "NLopt failed due to roundoff limitations: " << e.what() << std::endl;
			// different initial guess
			arma::mat x(dim,dim, arma::fill::randn);
			arma::mat R_full_mod=R_full+0.01*x;
			R_full_mod = (R_full_mod + R_full_mod.t()) / 2.0;
			R_flattened = arma::vectorise(R_full_mod);
			for (size_t i = 0; i < R_flattened.size(); ++i) X_flattened[i] = R_flattened[i];
			nlopt::result result = opt.optimize(X_flattened, minf);
		}

		// Convert flattened result back to matrix
		for (size_t i = 0; i < R_flattened.size(); ++i) R_flattened[i] = X_flattened[i];
		R_result = arma::reshape( R_flattened, dim, dim);

		// Symmetrize the result
		R_result = (R_result + R_result.t()) / 2.0;

		// Ensure diagonals are exactly 1
		//~ for (int i = 0; i < dim; ++i) R_result(i, i) = 1.0;

		// Ensure fixed part stays the same
		for (size_t i = 0; i < R_existing.n_rows; ++i) {
			for (size_t j = 0; j < R_existing.n_rows; ++j){
				R_result(i, j) = R_existing(i, j);
			}
		}

	}

	return R_result;
}

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
		//~ std::cout << "R_cross:\n" << R_cross << "\n\n";
		//~ std::cout << "R_new:\n" << R_new << "\n\n";

		arma::mat AB = arma::join_horiz( R_existing, R_cross.t());
		arma::mat CD = arma::join_horiz( R_cross, R_new);
		arma::mat R_full=arma::join_vert( AB, CD);

		//~ std::cout << "R_full:\n" << R_full << "\n\n";

		// create mask: 1 for values that must remain fixed, 0 for values that can be optimized
		int dim = R_full.n_rows;
		arma::mat mask = arma::zeros<arma::mat>(dim,dim);
		for (int i = 0; i < dim; ++i) {
			mask(i,i)=1;
		}
		for (int i = 0; i < R_existing.n_rows; ++i) {
			for (int j = 0; j < R_existing.n_rows; ++j) {
				mask(i,j)=1;
			}
		}
		mask(0,2)=1;	// fix correlation between cosine and 1st overtone cosine (para.rho_Wfb12)
		mask(2,0)=1;
		//~ std::cout << "Mask:\n" << mask << "\n\n";

		// Adjust to nearest PSD correlation matrix
		arma::mat R_adjusted = nearestCorrPSD(R_full, R_existing, mask);
		//~ std::cout << "R_adjusted:\n" << R_adjusted << "\n\n";

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
		//~ std::cout << para.Ncomp << "\t" << para.Ncomp_intf << std::endl;
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
		//~ std::cout << "R_full:\n" << R_full << "\n\n";

		// create mask: 1 for values that must remain fixed, 0 for values that can be optimized
		int dim = R_full.n_rows;
		arma::mat mask = arma::zeros<arma::mat>(dim,dim);
		for (int i = 0; i < dim; ++i) {
			mask(i,i)=1;
		}
		for (size_t i = 0; i < R_existing.n_rows; ++i) {
			for (size_t j = 0; j < R_existing.n_rows; ++j) {
				mask(i,j)=1;
			}
		}
		mask(0,2)=1;	// fix correlation between cosine and 1st overtone cosine (para.rho_Wfb12)
		mask(2,0)=1;
		mask(0,3)=1;	// fix correlation between cosine and 2st overtone cosine (para.rho_Wfb13)
		mask(3,0)=1;
		//~ std::cout << "Mask:\n" << mask << "\n\n";

		// Adjust to nearest PSD correlation matrix
		arma::mat R_adjusted = nearestCorrPSD(R_full, R_existing, mask, 0);
		//~ std::cout << "R_adjusted:\n" << R_adjusted << "\n\n";

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

		// save adjusted correlations for next repetition/parameter set
		//~ para.rho_Wfb=R_adjusted;

		for(int i= 0; i < para.Ncx; i++){
			for(int j= 2; j < para.Ncomp_intf; j++){
				//~ std::cout <<"additional SV: "<< j << std::endl;
				para.Intf_weights(i,j)=Intf_vec_tmp2[i];
			}
		}
	}

	// two signal frequencies, two overtones
	if(para.Ncomp_intf==4 && para.Ncomp==2){
		//~ std::cout << para.Ncomp << "\t" << para.Ncomp_intf << std::endl;
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
		arma::mat mask = arma::zeros<arma::mat>(dim,dim);
		for (int i = 0; i < dim; ++i) {
			mask(i,i)=1;
		}
		for (size_t i = 0; i < R_existing.n_rows; ++i) {
			for (size_t j = 0; j < R_existing.n_rows; ++j) {
				mask(i,j)=1;
			}
		}
		mask(0,4)=1;	// fix correlation between cosine and 1st overtone cosine (para.rho_Wfb12)
		mask(4,0)=1;
		mask(0,5)=1;	// fix correlation between cosine and 2st overtone cosine (para.rho_Wfb13)
		mask(5,0)=1;
		//~ std::cout << "Mask:\n" << mask << "\n\n";

		// Adjust to nearest PSD correlation matrix
		arma::mat R_adjusted = nearestCorrPSD(R_full, R_existing, mask, 0);
		//~ std::cout << "R_adjusted:\n" << R_adjusted << std::endl;
		arma::vec eigenvalues_adj = arma::eig_sym( R_adjusted );
		//~ std::cout << "EVs: " <<  eigenvalues << std::endl;
		int nreadj=0;
		if( arma::min(eigenvalues_adj) < 0){
			//~ std::cout << "R_adjusted invalid, min(EV) = " << eigenvalues_adj.min() << std::endl;
			while((arma::min(eigenvalues_adj) < 0) && (nreadj<25)) {
				// re-run with R_adjusted from last iteration
				R_adjusted = nearestCorrPSD(R_adjusted, R_existing, mask, nreadj);
				eigenvalues_adj = arma::eig_sym( R_adjusted );
				nreadj++;
				//~ if(nreadj>15) std::cout << "nreadj step = " << nreadj << std::endl;
			}
			//~ std::cout << "new R_adjusted = " << std::endl;
			//~ std::cout << R_adjusted << std::endl;
			//~ std::cout << "new min(EV) = " << eigenvalues_adj.min() << std::endl;
			//~ std::cout << "nreadj = " << nreadj << std::endl;
		}
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

		//~ std::cout << arma::size(para.Intf_weights) << "\t" << arma::size(Intf_vec_tmp) << "\t" << arma::size(spann_vec_RNN_intf1) << "\n\n";
		//~ std::cout << "corr spann_vec_RNN_intf1: " << std::endl;
		//~ std::cout << arma::cor( spann_vec_RNN_intf1 ) << std::endl;

		arma::vec Intf_vec_tmp2 = complement_multiD(spann_vec_RNN_intf1, rho_tmp2);
		//~ std::cout << "corr SVs: " << std::endl;
		//~ std::cout << arma::cor( arma::join_horiz( spann_vec_RNN_intf1, Intf_vec_tmp2) ) << std::endl;
		//~ std::cout << "rho_Wfb:\n" << para.rho_Wfb << "\n\n";
		//~ std::cout << "R_adjusted:\n" << R_adjusted << std::endl;

		// save adjusted correlations for next repetition/parameter set
		//~ para.rho_Wfb=R_adjusted;

		for(int i= 0; i < para.Ncx; i++){
			for(int j= 3; j < para.Ncomp_intf; j++){
				para.Intf_weights(i,j)=Intf_vec_tmp2[i];
			}
		}
	}

}

