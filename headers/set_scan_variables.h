arma::vec set_scan_variables(int &Npara, int learnflag, netpara &para) {

	arma::vec SCANPARA;
//------------------------------------------------------------------------
//	parameters for scans
//------------------------------------------------------------------------
	//~ arma::vec SIGMAN = arma::join_cols(arma::zeros(1) ,arma::logspace<arma::vec>(-2, 0.1, Npara-1));
	//~ arma::vec SIGMAN = arma::join_cols(arma::zeros(1) ,arma::logspace<arma::vec>(-1.5, 0.3, Npara-1));
	//~ arma::vec SIGMAN = arma::join_cols(arma::zeros(1) ,arma::logspace<arma::vec>(-7, -1, 19));
	arma::vec SIGMAN = arma::logspace<arma::vec>(-5, -1, 10);
	arma::vec THRESH;
	if(para.actfuncflag==0){
		if(para.simplemapflag==1){
			THRESH = {-4, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 0.70,  0.8,  1, 1.2, 1.5, 2, 2.5};	// non-linear network
			//~ THRESH = { -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5};
			THRESH=THRESH/4.0;
		}
		else{
			THRESH = {-4, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 0.70,  0.8,  1, 1.2, 1.5, 2, 2.5};	// non-linear network
			//~ THRESH = { -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5};
		}
	}
	else{
		if(para.simplemapflag==1){
			//~ THRESH = { -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5};
			//~ THRESH = {-2, -1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 0.70,  0.8,  1, 1.2, 1.5};
			//~ THRESH = {-4, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 0.70,  0.8,  1, 1.2, 1.5};
			//~ THRESH = {-5, -2, -1.5, -1.3, -1.1, -1, -0.9, -0.7, -0.5};
			//~ THRESH = arma::join_cols(-5 *arma::ones(1) ,arma::linspace<arma::vec>(-2, -0.25, 18));
			//~ THRESH = arma::join_cols(-5 *arma::ones(1) ,arma::linspace<arma::vec>(-2, 1.5, 40));
			THRESH = arma::linspace<arma::vec>(-1, 1, 21);
			//~ THRESH=THRESH/4.0;
		}
		else{
			//~ THRESH = { -4, -3, -2, -1, 0, 1, 2, 2.5};
			THRESH = arma::linspace<arma::vec>(-1, 1, 21);
		}
	}
	arma::vec GCGAIN = arma::logspace<arma::vec>(-1.699, 2,20);
	arma::vec NGC = { 250, 400, 500, 600, 700, 800, 1000, 1500, 2000, 3000};
	arma::vec NCX = { 200, 400, 500, 600, 700, 800, 1000, 1500, 2000, 3000};
	arma::vec NPC = { 1, 2, 3, 4, 5, 6, 7, 8};
	arma::vec NFB = { 1, 2, 10, 20, 100, 200, 500, 1000};
	//~ arma::vec WGOC = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5};
	//~ arma::vec WGOC = arma::linspace<arma::vec>(0, 2, 20);
	arma::vec WGOC =  arma::join_cols(1.0*arma::zeros(1),  arma::logspace<arma::vec>(-1, 0.301, 20));
	arma::vec SIGAMP =arma::logspace<arma::vec>(-1, 0.7 ,20);
	arma::vec LAMBDA;
	if(para.actfuncflag==0 && para.g > 1.0){
		LAMBDA = {0.4, 0.8, 1, 1.5, 2, 3, 4, 5, 7, 11, 13, 15, 20, 25};	// non-linear + chaotic network
		//~ LAMBDA = {1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 7, 11, 13, 15, 20};	// non-linear + chaotic network
	}
	//~ else if (para.g>0.8){
		// LAMBDA = {0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50, 60, 70};
		//~ LAMBDA = arma::logspace<arma::vec>(-0.30103, 1.8451, 25);
	//~ }
	else {
		LAMBDA = arma::logspace<arma::vec>(0, 1.8451, 25);
	}
	//~ LAMBDA = arma::logspace<arma::vec>(0, 2.5, 25);
	//~ LAMBDA = arma::linspace<arma::vec>(3, 24, 10);
	//~ LAMBDA = arma::linspace<arma::vec>(30, 40, 11);
	arma::vec PHI  = arma::linspace<arma::vec>(0, M_PI, 20);
	arma::vec AMPMOD  = arma::linspace<arma::vec>(-0.7, 0.7, Npara);
	arma::vec SIGMALS  = {0.00001, 0.0001, 0.001, 0.01, 0.1};
	arma::vec TSIM = {2, 4, 8, 10, 12, 15, 20, 25, 30, 40};
	arma::vec TUNSUP = {0, 10, 20, 40, 80, 120};
	arma::vec G = arma::linspace<arma::vec>(0.05,0.95,19);
	//~ arma::vec ETA = {0.0001, 0.0002, 0.0003, 0.0004, 0.0005};
	arma::vec ETA = arma::logspace<arma::vec>(-4, -2, 10);
	arma::vec ALPHAN = {-1, -0.5, -0.25, 0, 1, 2, 3, 4};
	arma::vec BETAN =arma::linspace<arma::vec>(0.0, 1.0, 10);
	arma::vec PZ = {0, 0.1, 0.3, 0.5, 0.7, 0.9};
	arma::vec FCODING =arma::linspace<arma::vec>(0.01, 0.98, 20);
	//~ arma::vec FCODING =arma::linspace<arma::vec>(0.02, 0.98, 20);
	//~ arma::vec FCODING =arma::linspace<arma::vec>(0.80, 0.98, 4);
	arma::vec SPALPHA =arma::logspace<arma::vec>(0, 2, 10);
	arma::vec PERTDELTA =arma::linspace<arma::vec>(0.0, 2, 10);
	arma::vec PERTDURATION =arma::linspace<arma::vec>(0.0, 1.0, 10);
	arma::vec PERTSTART =arma::linspace<arma::vec>(0.0, 1.0, 10);
	//~ arma::vec RHOWFB =arma::linspace<arma::vec>(0.0, 1.0, 11);
	//~ arma::vec RHOWFB = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.97, 0.98, 0.99, 0.999, 1};
	//~ arma::vec RHOWFB = arma::join_cols( arma::tanh(arma::linspace(0,3,20)), 1.0*arma::ones(1));
	arma::vec RHOWFB = arma::join_cols( arma::tanh(arma::linspace(0,4,20)), 1.0*arma::ones(1));
	arma::vec SIGAMPINTF =arma::logspace<arma::vec>(-2, 2 ,20);
	arma::vec GAINMAX =arma::logspace<arma::vec>(0,2,20);

	para.algswitch=1;
	if(learnflag==1){
		para.algswitch=1;
		para.Nrep=1;
		Npara=1;
		std::cout << "Learning with cortical readout unit." << std::endl;
	}
	else if(learnflag==2){
		para.algswitch=2;
		para.Nrep=1;
		Npara=1;
		std::cout << "Learning with cerebellar cortex." << std::endl;
	}
	else if(learnflag==99){
		Npara=1;
		std::cout << "Repeated learning over " << para.Nreal << " initial conditions." << std::endl;
	}
	else if(learnflag==100) {
		Npara=SIGMAN.n_elem;
		SCANPARA=SIGMAN;
	}
	else if(learnflag==101){
		Npara=THRESH.n_elem;
		SCANPARA=THRESH;
	}
	else if(learnflag==102){
		Npara=NCX.n_elem;
		SCANPARA=NCX;
	}
	else if(learnflag==103){
		Npara=NGC.n_elem;
		SCANPARA=NGC;
	}
	//~ else if(learnflag==104){
		//~ Npara=NPC.n_elem;
		//~ para.Npcmax=NPC(Npara-1);
		//~ SCANPARA=NPC;
	//~ }
	else if(learnflag==105){
		Npara=SIGMAN.n_elem;
		SCANPARA=SIGMAN;
	}
	else if(learnflag==106){
		Npara=WGOC.n_elem;
		SCANPARA=WGOC;
	}
	//~ else if(learnflag==107){
		//~ Npara=NPC.n_elem;
		//~ para.Ncomp=NPC(Npara-1);
		//~ SCANPARA=NPC;
	//~ }
	else if(learnflag==108){
		Npara=SIGAMP.n_elem;
		SCANPARA=SIGAMP;
	}
	else if(learnflag==109){
		Npara=LAMBDA.n_elem;
		SCANPARA=LAMBDA;
	}
	else if(learnflag==110){
		Npara=TSIM.n_elem;
		SCANPARA=TSIM;
	}
	else if(learnflag==111){
		Npara=NFB.n_elem;
		SCANPARA=NCX;
	}
	else if(learnflag==113){
		Npara=G.n_elem;
		SCANPARA=G;
	}
	else if(learnflag==114){
		Npara=ETA.n_elem;
		SCANPARA=ETA;
	}
	else if(learnflag==115){
		Npara=TUNSUP.n_elem;
		SCANPARA=TUNSUP;
	}
	else if(learnflag==116){
		std::cout << "Varying parameters of target signal. Scanning over " << Npara << " signals with " << para.Nreal << " initial conditions per signal." << std::endl;
	}
	else if(learnflag==117){
		Npara=NGC.n_elem;
		SCANPARA=NGC;
		std::cout << "Varying Ngc while fixing number of active GCs. Scanning over " << Npara << " values with " << para.Nreal << " initial conditions each." << std::endl;
		std::cout << "Number of active GCs = " << para.Ngc_act << std::endl;
	}
	else if(learnflag==118){
		Npara=PHI.n_elem;
		SCANPARA=PHI;
	}
	else if(learnflag==119){
		Npara=AMPMOD.n_elem;
		SCANPARA=AMPMOD;
	}
	else if(learnflag==120){
		Npara=SIGMALS.n_elem;
		SCANPARA=SIGMALS;
	}
	else if(learnflag==121){
		Npara=ALPHAN.n_elem;
		SCANPARA=ALPHAN;
	}
	else if(learnflag==122){
		Npara=PZ.n_elem;
		SCANPARA=PZ;
	}
	else if(learnflag==123){
		Npara=FCODING.n_elem;
		SCANPARA=FCODING;
	}
	else if(learnflag==124){
		Npara=SPALPHA.n_elem;
		SCANPARA=SPALPHA;
	}
	else if(learnflag==125){
		Npara=PERTDELTA.n_elem;
		SCANPARA=PERTDELTA;
	}
	else if(learnflag==126){
		Npara=PERTDURATION.n_elem;
		SCANPARA=PERTDURATION;
	}
	else if(learnflag==127){
		Npara=PERTSTART.n_elem;
		SCANPARA=PERTSTART;
	}
	else if(learnflag==128 || learnflag==132){
		Npara=RHOWFB.n_elem;
		SCANPARA=RHOWFB;
	}
	else if(learnflag==129){
		Npara=SIGAMPINTF.n_elem;
		SCANPARA=SIGAMPINTF;
	}
	else if(learnflag==130){
		Npara=BETAN.n_elem;
		SCANPARA=BETAN;
	}
	else if(learnflag==131){
		Npara=GCGAIN.n_elem;
		SCANPARA=GCGAIN;
	}
	else if(learnflag==133){
		Npara=THRESH.n_elem;
		SCANPARA=THRESH;
	}
	else if(learnflag==134){
		Npara=GAINMAX.n_elem;
		SCANPARA=GAINMAX;
	}
	else {
		para.algswitch=0;
		para.Nrep=1;
		Npara=1;
		std::cout << "No learning." << std::endl;
	}

	return SCANPARA;
}
