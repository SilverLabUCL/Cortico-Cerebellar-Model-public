void write_single_ctx(std::ofstream &single_ctxunit_files, std::ofstream &single_output_files, double &time, netpara &para) {

	int idx=0;

	single_ctxunit_files  << time << "\t"<< para.nuCx[idx] << "\t"<< para.nuCx[idx+1] << "\t"<< para.nuCx[idx+2] << "\t"<< para.nuCx[idx+3] << std::endl;
	single_output_files<< time << "\t"<< para.Z << std::endl;

	return;
}
