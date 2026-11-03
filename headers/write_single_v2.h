//~ void write_single(arma::colvec& nuCx, arma::colvec& GC, std::ofstream &single_ctxunit_files, std::ofstream &single_gc_files, std::ofstream &single_output_files, double &tsingle, double &time, netpara &para) {
void write_single(std::ofstream &single_ctxunit_files, std::ofstream &single_gc_files, std::ofstream &single_output_files, double &time, netpara &para) {

	int idx=0;

	single_ctxunit_files  << time << "\t"<< para.nuCx[idx] << "\t"<< para.nuCx[idx+1] << "\t"<< para.nuCx[idx+2] << "\t"<< para.nuCx[idx+3] << std::endl;
	single_gc_files << time << "\t"<< para.nuGC[idx] << "\t"<< para.nuGC[idx+1] << "\t"<< para.nuGC[idx+2] << "\t"<< para.nuGC[idx+3] << std::endl;
	single_output_files<< time << "\t"<< para.Z << "\t"<< para.PC << std::endl;

	return;
}
