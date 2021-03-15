#include "trans_data.h"


void trans_data::computeRelativePositions(int centromere1_start, int centromere1_stop, int centromere2_start, int centromere2_stop) {
	vrb.title("Compute positions relative to centromeres and telomeres");

	//Compute MIN-MAX
	int phenotype_min_pos1 = 1000000000;
	for (int p = 0 ; p < phenotype_count1 ; p++) if (phenotype_start1[p] < phenotype_min_pos1) phenotype_min_pos1 = phenotype_start1[p];
	int phenotype_min_pos2 = 1000000000;
	for (int p = 0 ; p < phenotype_count2 ; p++) if (phenotype_start2[p] < phenotype_min_pos2) phenotype_min_pos2 = phenotype_start2[p];
	int phenotype_max_pos1 = 0;
	for (int p = 0 ; p < phenotype_count1 ; p++) if (phenotype_end1[p] > phenotype_max_pos1) phenotype_max_pos1 = phenotype_end1[p];
	int phenotype_max_pos2 = 0;
	for (int p = 0 ; p < phenotype_count2 ; p++) if (phenotype_end2[p] > phenotype_max_pos2) phenotype_max_pos2 = phenotype_end2[p];
	//Compute arm lengths
	int arm01 = centromere1_start - phenotype_min_pos1;
	int arm11 = phenotype_max_pos1 - centromere1_stop;
	int arm02 = centromere2_start - phenotype_min_pos2;
	int arm12 = phenotype_max_pos2 - centromere2_stop;
	
	//Verbose
	vrb.bullet("chromosome 1 coordinates = ["+ stb.str(phenotype_min_pos1) + ", " +stb.str(centromere1_start) + ", " + stb.str(centromere1_stop) + ", " +stb.str(phenotype_max_pos1)+ "] and arms length =[" +stb.str(arm01)+ ", " +stb.str(arm11)+ "]");
	vrb.bullet("chromosome 2 coordinates = ["+ stb.str(phenotype_min_pos2) + ", " +stb.str(centromere2_start) + ", " + stb.str(centromere2_stop) + ", " +stb.str(phenotype_max_pos2)+ "] and arms length =[" +stb.str(arm02)+ ", " +stb.str(arm12)+ "]");
	
	//Compute relative positions
	phenotype_rpos1 = vector < double > (phenotype_count1, 0.0);
	for (int p = 0 ; p < phenotype_count1 ; p++) {
		int pos1 = (phenotype_start1[p] + phenotype_end1[p]) / 2;
		if (pos1 > centromere1_start && pos1 < centromere1_stop) phenotype_rpos1[p] = 0;
		else if (pos1 < centromere1_start) phenotype_rpos1[p] = (centromere1_start - pos1) * 1.0 / arm01;
		else phenotype_rpos1[p] = (pos1 - centromere1_stop) * 1.0 / arm11;
	}

	phenotype_rpos2 = vector < double > (phenotype_count2, 0.0);
	for (int p = 0 ; p < phenotype_count2 ; p++) {
		int pos2 = (phenotype_start2[p] + phenotype_end2[p]) / 2;
		if (pos2 > centromere2_start && pos2 < centromere2_stop) phenotype_rpos2[p] = 0;
		else if (pos2 < centromere2_start) phenotype_rpos2[p] = (centromere2_start - pos2) * 1.0 / arm02;
		else phenotype_rpos2[p] = (pos2 - centromere2_stop) * 1.0 / arm12;
	}
	vrb.bullet("done");
}

void trans_data::binRelativePositions(int n_bins, string fout) {
	vrb.title("Bin positions relative to centromeres and telomeres");
	vector < vector < int > > counts = vector < vector < int > > (n_bins, vector < int > (n_bins, 0));

	for (int i = 0 ; i < phenotype_count1 ; i ++) {
		for (int j = 0 ; j < phenotype_count2 ; j ++) {
			int idx1 = (int)floor(phenotype_rpos1[i] * n_bins);
			int idx2 = (int)floor(phenotype_rpos2[j] * n_bins);
			if (idx1 < 0 || idx1 >= n_bins || idx2 < 0 || idx2 >= n_bins) cout << idx1 << " " << idx2 << " " << phenotype_rpos1[i] << " " << phenotype_rpos2[j] << endl;
			counts[idx1][idx2]++;
		}
	}

	output_file fdo (fout);
	for (int b1 = 0 ; b1 < n_bins ; b1 ++) for (int b2 = 0 ; b2 < n_bins ; b2 ++) fdo << b1 << " " << b2 << " " << counts[b1][b2] << endl;
	fdo.close();
}
