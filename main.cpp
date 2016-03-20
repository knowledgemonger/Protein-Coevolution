#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

using namespace std;

class MSA {
private:
	int number_of_sequences;
	int alignment_length;
	vector<string> sequences;
	vector<vector<int> > aa_counts;
	vector<vector<float> > aa_frequencies;
	vector<vector<double> > MI_with_gaps;
	vector<vector<double> > MI;
	vector<vector<double> > Pearson;
	vector<vector<double> > Pearson_with_gaps;
	vector<vector<double> > Dunn;
	vector<vector<double> > Dunn_with_gaps;
	vector<vector<double> > GOF;
	vector<vector<double> > OMES;
	vector<vector<double> > McBASC;

	vector<vector<double> > subset_Hx;
	vector<double> Hx_with_gaps;
	vector<vector<double> > Hxy;
	vector<vector<double> > Hxy_with_gaps;
	vector<vector<double> > MI_residual_with_gaps;
	vector<vector<double> > MI_residual;
	vector<vector<int> > subset_size;
	vector<int> MSA_to_reference;
	vector<int> reference_to_MSA;
	vector<vector <int> > helices;
	vector<int> sheets;
	vector<int> turns;
	vector<vector<vector <float> > > atomic_coordinates;
	vector<int> C_beta_index;
	
	vector<int> number_residue_atoms;
	
	vector<vector<float> > minimum_distances;
	vector<vector<float> > residue_distances;
	vector<vector<float> > C_alpha_distances;
	vector<vector<float> > C_beta_distances;
	
	
	vector <double>	mean_MI_with_i;
	vector <double>	std_MI_with_i;
	
	vector<vector<float> > cumulative_PAM1;
	vector<vector<int> > permutation_PAM1;
	
	vector<vector<float> > scoring_matrix;

	string Amino_Acids;
	
	
	int aa_index (char aa) {								// Converts an amino acid char to an index for use in arrays
		const int char_index [78] = {
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,			// -,.,/,0,1,2,3,4,5,6,7,8,9,:,;,<,=,>,?,@
			1,0,2,3,4,5,6,7,8,0,9,10,11,						// A,(B),C,D,E,F,G,H,I,(J),K,L,M
			12,0,13,14,15,16,17,0,18,19,0,20,0,				// N,(O),P,Q,R,S,T,(U),V,W,(X),Y,(Z)
			0,0,0,0,0,0,										// [,\,],^,_,`
			1,0,2,3,4,5,6,7,8,0,9,10,11,						// a,(b),c,d,e,f,g,h,i,(j),k,l,m
			12,0,13,14,15,16,17,0,18,19,0,20,0				// n,(o),p,q,r,s,t,(u),v,w,(x),y,(z)
		};
		return char_index[aa - '-'];
	}
	
	vector<vector<int> > count_joints (int i, int j){
		int i_index, j_index;
		vector<vector<int> > joint_counts (21,vector<int> (21,0));
		
		for (int n=0;n<number_of_sequences;n++) {
			i_index=aa_index(sequences.at(n)[i]);
			j_index=aa_index(sequences.at(n)[j]);
			joint_counts.at(i_index).at(j_index)++;
		}
		return joint_counts;
	}
	
	
	void clear_information () {
		MI_with_gaps.clear();
		MI.clear();
		subset_Hx.clear();
		Hxy.clear();
		Hx_with_gaps.clear();
		Hxy_with_gaps.clear();
		MI_residual_with_gaps.clear();
		MI_residual.clear();
		Dunn_with_gaps.clear();
		Dunn.clear();
		Pearson_with_gaps.clear();
		Pearson.clear();
		GOF.clear();
		OMES.clear();
		McBASC.clear();
		
		subset_size.clear();
	}
	
	void clear_tabulation () {
		aa_counts.clear();
		aa_frequencies.clear();
	}
	
	
public:
		void clear () {
			sequences.clear();
			aa_counts.clear();
			aa_frequencies.clear();
			MI_with_gaps.clear();
			MI.clear();
			subset_Hx.clear();
			Hxy.clear();
			Hx_with_gaps.clear();
			Hxy_with_gaps.clear();
			MI_residual_with_gaps.clear();
			MI_residual.clear();
			Dunn_with_gaps.clear();
			Dunn.clear();
			Pearson_with_gaps.clear();
			Pearson.clear();
			GOF.clear();
			OMES.clear();
			McBASC.clear();

			
			subset_size.clear();
			number_of_sequences = 0;
			alignment_length = 0;
			MSA_to_reference.clear();
			reference_to_MSA.clear();
			helices.clear();
			sheets.clear();
			turns.clear();
			atomic_coordinates.clear();
			C_beta_index.clear();
			number_residue_atoms.clear();
			
			minimum_distances.clear();
			residue_distances.clear();
			C_alpha_distances.clear();
			C_beta_distances.clear();
			
			
			mean_MI_with_i.clear();
			std_MI_with_i.clear();
						
		}
	
	
	MSA() {
		Amino_Acids = "ACDEFGHIKLMNPQRSTVWY";
		
		cumulative_PAM1 = vector < vector <float> > (20, vector <float> (20,0.0));
		permutation_PAM1 = vector < vector <int> > (20, vector <int> (20,0));
		
		scoring_matrix = vector < vector <float> > (21, vector <float> (21,0.0));
		
		cumulative_PAM1.at(0).at(0)=1.000000; cumulative_PAM1.at(0).at(1)=0.013300; cumulative_PAM1.at(0).at(2)=0.013200; cumulative_PAM1.at(0).at(3)=0.012600; cumulative_PAM1.at(0).at(4)=0.011600; cumulative_PAM1.at(0).at(5)=0.011500; cumulative_PAM1.at(0).at(6)=0.009400; cumulative_PAM1.at(0).at(7)=0.009300; cumulative_PAM1.at(0).at(8)=0.009100; cumulative_PAM1.at(0).at(9)=0.008900; cumulative_PAM1.at(0).at(10)=0.008600; cumulative_PAM1.at(0).at(11)=0.008500; cumulative_PAM1.at(0).at(12)=0.008100; cumulative_PAM1.at(0).at(13)=0.006800; cumulative_PAM1.at(0).at(14)=0.006500; cumulative_PAM1.at(0).at(15)=0.006400; cumulative_PAM1.at(0).at(16)=0.003600; cumulative_PAM1.at(0).at(17)=0.001400; cumulative_PAM1.at(0).at(18)=0.000100; cumulative_PAM1.at(0).at(19)=0.000100;
		cumulative_PAM1.at(1).at(0)=1.000000; cumulative_PAM1.at(1).at(1)=0.002700; cumulative_PAM1.at(1).at(2)=0.002400; cumulative_PAM1.at(1).at(3)=0.002400; cumulative_PAM1.at(1).at(4)=0.002400; cumulative_PAM1.at(1).at(5)=0.002400; cumulative_PAM1.at(1).at(6)=0.002300; cumulative_PAM1.at(1).at(7)=0.002200; cumulative_PAM1.at(1).at(8)=0.002000; cumulative_PAM1.at(1).at(9)=0.002000; cumulative_PAM1.at(1).at(10)=0.002000; cumulative_PAM1.at(1).at(11)=0.002000; cumulative_PAM1.at(1).at(12)=0.002000; cumulative_PAM1.at(1).at(13)=0.001900; cumulative_PAM1.at(1).at(14)=0.001900; cumulative_PAM1.at(1).at(15)=0.001800; cumulative_PAM1.at(1).at(16)=0.000700; cumulative_PAM1.at(1).at(17)=0.000600; cumulative_PAM1.at(1).at(18)=0.000300; cumulative_PAM1.at(1).at(19)=0.000300;
		cumulative_PAM1.at(2).at(0)=1.000000; cumulative_PAM1.at(2).at(1)=0.014100; cumulative_PAM1.at(2).at(2)=0.013100; cumulative_PAM1.at(2).at(3)=0.013100; cumulative_PAM1.at(2).at(4)=0.007500; cumulative_PAM1.at(2).at(5)=0.007500; cumulative_PAM1.at(2).at(6)=0.006400; cumulative_PAM1.at(2).at(7)=0.006100; cumulative_PAM1.at(2).at(8)=0.006000; cumulative_PAM1.at(2).at(9)=0.005400; cumulative_PAM1.at(2).at(10)=0.005400; cumulative_PAM1.at(2).at(11)=0.005400; cumulative_PAM1.at(2).at(12)=0.001800; cumulative_PAM1.at(2).at(13)=0.001700; cumulative_PAM1.at(2).at(14)=0.001200; cumulative_PAM1.at(2).at(15)=0.001200; cumulative_PAM1.at(2).at(16)=0.000500; cumulative_PAM1.at(2).at(17)=0.000100; cumulative_PAM1.at(2).at(18)=0.000000; cumulative_PAM1.at(2).at(19)=0.000000;
		cumulative_PAM1.at(3).at(0)=1.000000; cumulative_PAM1.at(3).at(1)=0.013500; cumulative_PAM1.at(3).at(2)=0.011800; cumulative_PAM1.at(3).at(3)=0.011800; cumulative_PAM1.at(3).at(4)=0.006500; cumulative_PAM1.at(3).at(5)=0.006500; cumulative_PAM1.at(3).at(6)=0.005800; cumulative_PAM1.at(3).at(7)=0.005700; cumulative_PAM1.at(3).at(8)=0.005500; cumulative_PAM1.at(3).at(9)=0.004800; cumulative_PAM1.at(3).at(10)=0.004700; cumulative_PAM1.at(3).at(11)=0.004700; cumulative_PAM1.at(3).at(12)=0.004100; cumulative_PAM1.at(3).at(13)=0.003800; cumulative_PAM1.at(3).at(14)=0.001100; cumulative_PAM1.at(3).at(15)=0.001100; cumulative_PAM1.at(3).at(16)=0.000500; cumulative_PAM1.at(3).at(17)=0.000300; cumulative_PAM1.at(3).at(18)=0.000100; cumulative_PAM1.at(3).at(19)=0.000100;
		cumulative_PAM1.at(4).at(0)=1.000200; cumulative_PAM1.at(4).at(1)=0.005600; cumulative_PAM1.at(4).at(2)=0.005400; cumulative_PAM1.at(4).at(3)=0.005400; cumulative_PAM1.at(4).at(4)=0.005400; cumulative_PAM1.at(4).at(5)=0.005400; cumulative_PAM1.at(4).at(6)=0.005300; cumulative_PAM1.at(4).at(7)=0.005100; cumulative_PAM1.at(4).at(8)=0.004400; cumulative_PAM1.at(4).at(9)=0.004400; cumulative_PAM1.at(4).at(10)=0.003100; cumulative_PAM1.at(4).at(11)=0.003000; cumulative_PAM1.at(4).at(12)=0.002900; cumulative_PAM1.at(4).at(13)=0.002800; cumulative_PAM1.at(4).at(14)=0.002800; cumulative_PAM1.at(4).at(15)=0.002700; cumulative_PAM1.at(4).at(16)=0.002400; cumulative_PAM1.at(4).at(17)=0.002300; cumulative_PAM1.at(4).at(18)=0.002200; cumulative_PAM1.at(4).at(19)=0.002100;
		cumulative_PAM1.at(5).at(0)=1.000000; cumulative_PAM1.at(5).at(1)=0.006500; cumulative_PAM1.at(5).at(2)=0.004400; cumulative_PAM1.at(5).at(3)=0.004400; cumulative_PAM1.at(5).at(4)=0.003800; cumulative_PAM1.at(5).at(5)=0.003400; cumulative_PAM1.at(5).at(6)=0.003300; cumulative_PAM1.at(5).at(7)=0.003300; cumulative_PAM1.at(5).at(8)=0.003300; cumulative_PAM1.at(5).at(9)=0.003100; cumulative_PAM1.at(5).at(10)=0.003000; cumulative_PAM1.at(5).at(11)=0.003000; cumulative_PAM1.at(5).at(12)=0.002400; cumulative_PAM1.at(5).at(13)=0.002200; cumulative_PAM1.at(5).at(14)=0.002100; cumulative_PAM1.at(5).at(15)=0.002100; cumulative_PAM1.at(5).at(16)=0.000500; cumulative_PAM1.at(5).at(17)=0.000300; cumulative_PAM1.at(5).at(18)=0.000000; cumulative_PAM1.at(5).at(19)=0.000000;
		cumulative_PAM1.at(6).at(0)=0.999900; cumulative_PAM1.at(6).at(1)=0.008700; cumulative_PAM1.at(6).at(2)=0.008500; cumulative_PAM1.at(6).at(3)=0.008400; cumulative_PAM1.at(6).at(4)=0.008000; cumulative_PAM1.at(6).at(5)=0.007800; cumulative_PAM1.at(6).at(6)=0.007600; cumulative_PAM1.at(6).at(7)=0.007500; cumulative_PAM1.at(6).at(8)=0.007500; cumulative_PAM1.at(6).at(9)=0.007300; cumulative_PAM1.at(6).at(10)=0.006900; cumulative_PAM1.at(6).at(11)=0.006900; cumulative_PAM1.at(6).at(12)=0.004800; cumulative_PAM1.at(6).at(13)=0.004300; cumulative_PAM1.at(6).at(14)=0.002000; cumulative_PAM1.at(6).at(15)=0.001000; cumulative_PAM1.at(6).at(16)=0.000800; cumulative_PAM1.at(6).at(17)=0.000700; cumulative_PAM1.at(6).at(18)=0.000400; cumulative_PAM1.at(6).at(19)=0.000400;
		cumulative_PAM1.at(7).at(0)=1.000100; cumulative_PAM1.at(7).at(1)=0.012900; cumulative_PAM1.at(7).at(2)=0.012300; cumulative_PAM1.at(7).at(3)=0.012200; cumulative_PAM1.at(7).at(4)=0.012100; cumulative_PAM1.at(7).at(5)=0.011800; cumulative_PAM1.at(7).at(6)=0.011000; cumulative_PAM1.at(7).at(7)=0.011000; cumulative_PAM1.at(7).at(8)=0.011000; cumulative_PAM1.at(7).at(9)=0.010600; cumulative_PAM1.at(7).at(10)=0.008400; cumulative_PAM1.at(7).at(11)=0.007900; cumulative_PAM1.at(7).at(12)=0.007600; cumulative_PAM1.at(7).at(13)=0.007500; cumulative_PAM1.at(7).at(14)=0.007400; cumulative_PAM1.at(7).at(15)=0.007100; cumulative_PAM1.at(7).at(16)=0.006900; cumulative_PAM1.at(7).at(17)=0.005800; cumulative_PAM1.at(7).at(18)=0.000100; cumulative_PAM1.at(7).at(19)=0.000100;
		cumulative_PAM1.at(8).at(0)=1.000200; cumulative_PAM1.at(8).at(1)=0.007600; cumulative_PAM1.at(8).at(2)=0.007400; cumulative_PAM1.at(8).at(3)=0.007400; cumulative_PAM1.at(8).at(4)=0.007100; cumulative_PAM1.at(8).at(5)=0.006700; cumulative_PAM1.at(8).at(6)=0.006700; cumulative_PAM1.at(8).at(7)=0.006500; cumulative_PAM1.at(8).at(8)=0.006400; cumulative_PAM1.at(8).at(9)=0.006200; cumulative_PAM1.at(8).at(10)=0.006000; cumulative_PAM1.at(8).at(11)=0.005600; cumulative_PAM1.at(8).at(12)=0.004300; cumulative_PAM1.at(8).at(13)=0.004100; cumulative_PAM1.at(8).at(14)=0.003500; cumulative_PAM1.at(8).at(15)=0.001600; cumulative_PAM1.at(8).at(16)=0.000900; cumulative_PAM1.at(8).at(17)=0.000100; cumulative_PAM1.at(8).at(18)=0.000000; cumulative_PAM1.at(8).at(19)=0.000000;
		cumulative_PAM1.at(9).at(0)=1.000000; cumulative_PAM1.at(9).at(1)=0.005300; cumulative_PAM1.at(9).at(2)=0.004900; cumulative_PAM1.at(9).at(3)=0.004900; cumulative_PAM1.at(9).at(4)=0.004900; cumulative_PAM1.at(9).at(5)=0.004800; cumulative_PAM1.at(9).at(6)=0.004200; cumulative_PAM1.at(9).at(7)=0.004100; cumulative_PAM1.at(9).at(8)=0.004000; cumulative_PAM1.at(9).at(9)=0.003100; cumulative_PAM1.at(9).at(10)=0.003000; cumulative_PAM1.at(9).at(11)=0.002200; cumulative_PAM1.at(9).at(12)=0.002100; cumulative_PAM1.at(9).at(13)=0.001900; cumulative_PAM1.at(9).at(14)=0.001600; cumulative_PAM1.at(9).at(15)=0.001500; cumulative_PAM1.at(9).at(16)=0.001400; cumulative_PAM1.at(9).at(17)=0.001200; cumulative_PAM1.at(9).at(18)=0.000100; cumulative_PAM1.at(9).at(19)=0.000100;
		cumulative_PAM1.at(10).at(0)=0.999900; cumulative_PAM1.at(10).at(1)=0.012500; cumulative_PAM1.at(10).at(2)=0.011900; cumulative_PAM1.at(10).at(3)=0.011900; cumulative_PAM1.at(10).at(4)=0.011900; cumulative_PAM1.at(10).at(5)=0.011800; cumulative_PAM1.at(10).at(6)=0.011400; cumulative_PAM1.at(10).at(7)=0.011300; cumulative_PAM1.at(10).at(8)=0.011300; cumulative_PAM1.at(10).at(9)=0.010100; cumulative_PAM1.at(10).at(10)=0.008100; cumulative_PAM1.at(10).at(11)=0.003600; cumulative_PAM1.at(10).at(12)=0.003600; cumulative_PAM1.at(10).at(13)=0.003500; cumulative_PAM1.at(10).at(14)=0.003100; cumulative_PAM1.at(10).at(15)=0.002700; cumulative_PAM1.at(10).at(16)=0.002300; cumulative_PAM1.at(10).at(17)=0.001700; cumulative_PAM1.at(10).at(18)=0.000000; cumulative_PAM1.at(10).at(19)=0.000000;
		cumulative_PAM1.at(11).at(0)=1.000000; cumulative_PAM1.at(11).at(1)=0.017800; cumulative_PAM1.at(11).at(2)=0.016900; cumulative_PAM1.at(11).at(3)=0.016900; cumulative_PAM1.at(11).at(4)=0.012700; cumulative_PAM1.at(11).at(5)=0.012000; cumulative_PAM1.at(11).at(6)=0.011900; cumulative_PAM1.at(11).at(7)=0.010700; cumulative_PAM1.at(11).at(8)=0.008900; cumulative_PAM1.at(11).at(9)=0.008600; cumulative_PAM1.at(11).at(10)=0.006100; cumulative_PAM1.at(11).at(11)=0.005800; cumulative_PAM1.at(11).at(12)=0.005800; cumulative_PAM1.at(11).at(13)=0.005600; cumulative_PAM1.at(11).at(14)=0.005200; cumulative_PAM1.at(11).at(15)=0.005100; cumulative_PAM1.at(11).at(16)=0.001700; cumulative_PAM1.at(11).at(17)=0.000400; cumulative_PAM1.at(11).at(18)=0.000300; cumulative_PAM1.at(11).at(19)=0.000300;
		cumulative_PAM1.at(12).at(0)=1.000200; cumulative_PAM1.at(12).at(1)=0.007600; cumulative_PAM1.at(12).at(2)=0.005400; cumulative_PAM1.at(12).at(3)=0.005300; cumulative_PAM1.at(12).at(4)=0.005200; cumulative_PAM1.at(12).at(5)=0.004900; cumulative_PAM1.at(12).at(6)=0.004900; cumulative_PAM1.at(12).at(7)=0.004600; cumulative_PAM1.at(12).at(8)=0.004300; cumulative_PAM1.at(12).at(9)=0.004300; cumulative_PAM1.at(12).at(10)=0.004000; cumulative_PAM1.at(12).at(11)=0.003700; cumulative_PAM1.at(12).at(12)=0.003700; cumulative_PAM1.at(12).at(13)=0.003500; cumulative_PAM1.at(12).at(14)=0.002900; cumulative_PAM1.at(12).at(15)=0.002500; cumulative_PAM1.at(12).at(16)=0.000800; cumulative_PAM1.at(12).at(17)=0.000300; cumulative_PAM1.at(12).at(18)=0.000000; cumulative_PAM1.at(12).at(19)=0.000000;
		cumulative_PAM1.at(13).at(0)=1.000000; cumulative_PAM1.at(13).at(1)=0.012400; cumulative_PAM1.at(13).at(2)=0.011600; cumulative_PAM1.at(13).at(3)=0.011600; cumulative_PAM1.at(13).at(4)=0.011000; cumulative_PAM1.at(13).at(5)=0.007500; cumulative_PAM1.at(13).at(6)=0.007500; cumulative_PAM1.at(13).at(7)=0.007200; cumulative_PAM1.at(13).at(8)=0.005200; cumulative_PAM1.at(13).at(9)=0.005100; cumulative_PAM1.at(13).at(10)=0.003900; cumulative_PAM1.at(13).at(11)=0.003300; cumulative_PAM1.at(13).at(12)=0.003100; cumulative_PAM1.at(13).at(13)=0.002700; cumulative_PAM1.at(13).at(14)=0.001900; cumulative_PAM1.at(13).at(15)=0.000900; cumulative_PAM1.at(13).at(16)=0.000500; cumulative_PAM1.at(13).at(17)=0.000200; cumulative_PAM1.at(13).at(18)=0.000000; cumulative_PAM1.at(13).at(19)=0.000000;
		cumulative_PAM1.at(14).at(0)=0.999900; cumulative_PAM1.at(14).at(1)=0.008600; cumulative_PAM1.at(14).at(2)=0.008400; cumulative_PAM1.at(14).at(3)=0.008300; cumulative_PAM1.at(14).at(4)=0.008300; cumulative_PAM1.at(14).at(5)=0.008300; cumulative_PAM1.at(14).at(6)=0.008200; cumulative_PAM1.at(14).at(7)=0.008100; cumulative_PAM1.at(14).at(8)=0.007300; cumulative_PAM1.at(14).at(9)=0.007100; cumulative_PAM1.at(14).at(10)=0.003400; cumulative_PAM1.at(14).at(11)=0.003300; cumulative_PAM1.at(14).at(12)=0.003200; cumulative_PAM1.at(14).at(13)=0.003100; cumulative_PAM1.at(14).at(14)=0.002600; cumulative_PAM1.at(14).at(15)=0.001700; cumulative_PAM1.at(14).at(16)=0.000600; cumulative_PAM1.at(14).at(17)=0.000400; cumulative_PAM1.at(14).at(18)=0.000200; cumulative_PAM1.at(14).at(19)=0.000000;
		cumulative_PAM1.at(15).at(0)=1.000000; cumulative_PAM1.at(15).at(1)=0.016000; cumulative_PAM1.at(15).at(2)=0.012500; cumulative_PAM1.at(15).at(3)=0.012000; cumulative_PAM1.at(15).at(4)=0.011500; cumulative_PAM1.at(15).at(5)=0.011100; cumulative_PAM1.at(15).at(6)=0.010900; cumulative_PAM1.at(15).at(7)=0.008800; cumulative_PAM1.at(15).at(8)=0.008700; cumulative_PAM1.at(15).at(9)=0.008600; cumulative_PAM1.at(15).at(10)=0.007800; cumulative_PAM1.at(15).at(11)=0.007700; cumulative_PAM1.at(15).at(12)=0.007600; cumulative_PAM1.at(15).at(13)=0.005600; cumulative_PAM1.at(15).at(14)=0.004400; cumulative_PAM1.at(15).at(15)=0.004200; cumulative_PAM1.at(15).at(16)=0.003600; cumulative_PAM1.at(15).at(17)=0.000400; cumulative_PAM1.at(15).at(18)=0.000200; cumulative_PAM1.at(15).at(19)=0.000100;
		cumulative_PAM1.at(16).at(0)=1.000200; cumulative_PAM1.at(16).at(1)=0.013100; cumulative_PAM1.at(16).at(2)=0.009900; cumulative_PAM1.at(16).at(3)=0.009800; cumulative_PAM1.at(16).at(4)=0.009500; cumulative_PAM1.at(16).at(5)=0.009300; cumulative_PAM1.at(16).at(6)=0.009200; cumulative_PAM1.at(16).at(7)=0.008900; cumulative_PAM1.at(16).at(8)=0.008800; cumulative_PAM1.at(16).at(9)=0.008100; cumulative_PAM1.at(16).at(10)=0.007000; cumulative_PAM1.at(16).at(11)=0.006700; cumulative_PAM1.at(16).at(12)=0.006500; cumulative_PAM1.at(16).at(13)=0.005600; cumulative_PAM1.at(16).at(14)=0.005200; cumulative_PAM1.at(16).at(15)=0.005000; cumulative_PAM1.at(16).at(16)=0.004900; cumulative_PAM1.at(16).at(17)=0.001100; cumulative_PAM1.at(16).at(18)=0.000100; cumulative_PAM1.at(16).at(19)=0.000100;
		cumulative_PAM1.at(17).at(0)=1.000000; cumulative_PAM1.at(17).at(1)=0.009900; cumulative_PAM1.at(17).at(2)=0.008100; cumulative_PAM1.at(17).at(3)=0.007900; cumulative_PAM1.at(17).at(4)=0.007800; cumulative_PAM1.at(17).at(5)=0.007600; cumulative_PAM1.at(17).at(6)=0.007600; cumulative_PAM1.at(17).at(7)=0.007100; cumulative_PAM1.at(17).at(8)=0.007000; cumulative_PAM1.at(17).at(9)=0.003700; cumulative_PAM1.at(17).at(10)=0.003600; cumulative_PAM1.at(17).at(11)=0.002100; cumulative_PAM1.at(17).at(12)=0.001700; cumulative_PAM1.at(17).at(13)=0.001600; cumulative_PAM1.at(17).at(14)=0.001400; cumulative_PAM1.at(17).at(15)=0.001300; cumulative_PAM1.at(17).at(16)=0.001200; cumulative_PAM1.at(17).at(17)=0.001000; cumulative_PAM1.at(17).at(18)=0.000100; cumulative_PAM1.at(17).at(19)=0.000100;
		cumulative_PAM1.at(18).at(0)=1.000000; cumulative_PAM1.at(18).at(1)=0.002400; cumulative_PAM1.at(18).at(2)=0.002400; cumulative_PAM1.at(18).at(3)=0.002400; cumulative_PAM1.at(18).at(4)=0.002400; cumulative_PAM1.at(18).at(5)=0.002400; cumulative_PAM1.at(18).at(6)=0.002100; cumulative_PAM1.at(18).at(7)=0.002100; cumulative_PAM1.at(18).at(8)=0.002000; cumulative_PAM1.at(18).at(9)=0.002000; cumulative_PAM1.at(18).at(10)=0.002000; cumulative_PAM1.at(18).at(11)=0.001600; cumulative_PAM1.at(18).at(12)=0.001600; cumulative_PAM1.at(18).at(13)=0.001500; cumulative_PAM1.at(18).at(14)=0.001500; cumulative_PAM1.at(18).at(15)=0.001500; cumulative_PAM1.at(18).at(16)=0.000700; cumulative_PAM1.at(18).at(17)=0.000200; cumulative_PAM1.at(18).at(18)=0.000200; cumulative_PAM1.at(18).at(19)=0.000200;
		cumulative_PAM1.at(19).at(0)=0.999800; cumulative_PAM1.at(19).at(1)=0.005300; cumulative_PAM1.at(19).at(2)=0.005100; cumulative_PAM1.at(19).at(3)=0.004800; cumulative_PAM1.at(19).at(4)=0.004800; cumulative_PAM1.at(19).at(5)=0.004700; cumulative_PAM1.at(19).at(6)=0.001900; cumulative_PAM1.at(19).at(7)=0.001900; cumulative_PAM1.at(19).at(8)=0.001500; cumulative_PAM1.at(19).at(9)=0.001400; cumulative_PAM1.at(19).at(10)=0.001300; cumulative_PAM1.at(19).at(11)=0.001100; cumulative_PAM1.at(19).at(12)=0.001100; cumulative_PAM1.at(19).at(13)=0.000700; cumulative_PAM1.at(19).at(14)=0.000700; cumulative_PAM1.at(19).at(15)=0.000700; cumulative_PAM1.at(19).at(16)=0.000700; cumulative_PAM1.at(19).at(17)=0.000500; cumulative_PAM1.at(19).at(18)=0.000300; cumulative_PAM1.at(19).at(19)=0.000100;
		
		permutation_PAM1.at(0).at(0)=0; permutation_PAM1.at(0).at(1)=1; permutation_PAM1.at(0).at(2)=2; permutation_PAM1.at(0).at(3)=3; permutation_PAM1.at(0).at(4)=4; permutation_PAM1.at(0).at(5)=5; permutation_PAM1.at(0).at(6)=6; permutation_PAM1.at(0).at(7)=7; permutation_PAM1.at(0).at(8)=8; permutation_PAM1.at(0).at(9)=9; permutation_PAM1.at(0).at(10)=10; permutation_PAM1.at(0).at(11)=11; permutation_PAM1.at(0).at(12)=12; permutation_PAM1.at(0).at(13)=13; permutation_PAM1.at(0).at(14)=14; permutation_PAM1.at(0).at(15)=15; permutation_PAM1.at(0).at(16)=16; permutation_PAM1.at(0).at(17)=17; permutation_PAM1.at(0).at(18)=18; permutation_PAM1.at(0).at(19)=19;
		permutation_PAM1.at(1).at(0)=1; permutation_PAM1.at(1).at(1)=0; permutation_PAM1.at(1).at(2)=2; permutation_PAM1.at(1).at(3)=3; permutation_PAM1.at(1).at(4)=4; permutation_PAM1.at(1).at(5)=5; permutation_PAM1.at(1).at(6)=6; permutation_PAM1.at(1).at(7)=7; permutation_PAM1.at(1).at(8)=8; permutation_PAM1.at(1).at(9)=9; permutation_PAM1.at(1).at(10)=10; permutation_PAM1.at(1).at(11)=11; permutation_PAM1.at(1).at(12)=12; permutation_PAM1.at(1).at(13)=13; permutation_PAM1.at(1).at(14)=14; permutation_PAM1.at(1).at(15)=15; permutation_PAM1.at(1).at(16)=16; permutation_PAM1.at(1).at(17)=17; permutation_PAM1.at(1).at(18)=18; permutation_PAM1.at(1).at(19)=19;
		permutation_PAM1.at(2).at(0)=2; permutation_PAM1.at(2).at(1)=0; permutation_PAM1.at(2).at(2)=1; permutation_PAM1.at(2).at(3)=3; permutation_PAM1.at(2).at(4)=4; permutation_PAM1.at(2).at(5)=5; permutation_PAM1.at(2).at(6)=6; permutation_PAM1.at(2).at(7)=7; permutation_PAM1.at(2).at(8)=8; permutation_PAM1.at(2).at(9)=9; permutation_PAM1.at(2).at(10)=10; permutation_PAM1.at(2).at(11)=11; permutation_PAM1.at(2).at(12)=12; permutation_PAM1.at(2).at(13)=13; permutation_PAM1.at(2).at(14)=14; permutation_PAM1.at(2).at(15)=15; permutation_PAM1.at(2).at(16)=16; permutation_PAM1.at(2).at(17)=17; permutation_PAM1.at(2).at(18)=18; permutation_PAM1.at(2).at(19)=19;
		permutation_PAM1.at(3).at(0)=3; permutation_PAM1.at(3).at(1)=0; permutation_PAM1.at(3).at(2)=1; permutation_PAM1.at(3).at(3)=2; permutation_PAM1.at(3).at(4)=4; permutation_PAM1.at(3).at(5)=5; permutation_PAM1.at(3).at(6)=6; permutation_PAM1.at(3).at(7)=7; permutation_PAM1.at(3).at(8)=8; permutation_PAM1.at(3).at(9)=9; permutation_PAM1.at(3).at(10)=10; permutation_PAM1.at(3).at(11)=11; permutation_PAM1.at(3).at(12)=12; permutation_PAM1.at(3).at(13)=13; permutation_PAM1.at(3).at(14)=14; permutation_PAM1.at(3).at(15)=15; permutation_PAM1.at(3).at(16)=16; permutation_PAM1.at(3).at(17)=17; permutation_PAM1.at(3).at(18)=18; permutation_PAM1.at(3).at(19)=19;
		permutation_PAM1.at(4).at(0)=4; permutation_PAM1.at(4).at(1)=0; permutation_PAM1.at(4).at(2)=1; permutation_PAM1.at(4).at(3)=2; permutation_PAM1.at(4).at(4)=3; permutation_PAM1.at(4).at(5)=5; permutation_PAM1.at(4).at(6)=6; permutation_PAM1.at(4).at(7)=7; permutation_PAM1.at(4).at(8)=8; permutation_PAM1.at(4).at(9)=9; permutation_PAM1.at(4).at(10)=10; permutation_PAM1.at(4).at(11)=11; permutation_PAM1.at(4).at(12)=12; permutation_PAM1.at(4).at(13)=13; permutation_PAM1.at(4).at(14)=14; permutation_PAM1.at(4).at(15)=15; permutation_PAM1.at(4).at(16)=16; permutation_PAM1.at(4).at(17)=17; permutation_PAM1.at(4).at(18)=18; permutation_PAM1.at(4).at(19)=19;
		permutation_PAM1.at(5).at(0)=5; permutation_PAM1.at(5).at(1)=0; permutation_PAM1.at(5).at(2)=1; permutation_PAM1.at(5).at(3)=2; permutation_PAM1.at(5).at(4)=3; permutation_PAM1.at(5).at(5)=4; permutation_PAM1.at(5).at(6)=6; permutation_PAM1.at(5).at(7)=7; permutation_PAM1.at(5).at(8)=8; permutation_PAM1.at(5).at(9)=9; permutation_PAM1.at(5).at(10)=10; permutation_PAM1.at(5).at(11)=11; permutation_PAM1.at(5).at(12)=12; permutation_PAM1.at(5).at(13)=13; permutation_PAM1.at(5).at(14)=14; permutation_PAM1.at(5).at(15)=15; permutation_PAM1.at(5).at(16)=16; permutation_PAM1.at(5).at(17)=17; permutation_PAM1.at(5).at(18)=18; permutation_PAM1.at(5).at(19)=19;
		permutation_PAM1.at(6).at(0)=6; permutation_PAM1.at(6).at(1)=0; permutation_PAM1.at(6).at(2)=1; permutation_PAM1.at(6).at(3)=2; permutation_PAM1.at(6).at(4)=3; permutation_PAM1.at(6).at(5)=4; permutation_PAM1.at(6).at(6)=5; permutation_PAM1.at(6).at(7)=7; permutation_PAM1.at(6).at(8)=8; permutation_PAM1.at(6).at(9)=9; permutation_PAM1.at(6).at(10)=10; permutation_PAM1.at(6).at(11)=11; permutation_PAM1.at(6).at(12)=12; permutation_PAM1.at(6).at(13)=13; permutation_PAM1.at(6).at(14)=14; permutation_PAM1.at(6).at(15)=15; permutation_PAM1.at(6).at(16)=16; permutation_PAM1.at(6).at(17)=17; permutation_PAM1.at(6).at(18)=18; permutation_PAM1.at(6).at(19)=19;
		permutation_PAM1.at(7).at(0)=7; permutation_PAM1.at(7).at(1)=0; permutation_PAM1.at(7).at(2)=1; permutation_PAM1.at(7).at(3)=2; permutation_PAM1.at(7).at(4)=3; permutation_PAM1.at(7).at(5)=4; permutation_PAM1.at(7).at(6)=5; permutation_PAM1.at(7).at(7)=6; permutation_PAM1.at(7).at(8)=8; permutation_PAM1.at(7).at(9)=9; permutation_PAM1.at(7).at(10)=10; permutation_PAM1.at(7).at(11)=11; permutation_PAM1.at(7).at(12)=12; permutation_PAM1.at(7).at(13)=13; permutation_PAM1.at(7).at(14)=14; permutation_PAM1.at(7).at(15)=15; permutation_PAM1.at(7).at(16)=16; permutation_PAM1.at(7).at(17)=17; permutation_PAM1.at(7).at(18)=18; permutation_PAM1.at(7).at(19)=19;
		permutation_PAM1.at(8).at(0)=8; permutation_PAM1.at(8).at(1)=0; permutation_PAM1.at(8).at(2)=1; permutation_PAM1.at(8).at(3)=2; permutation_PAM1.at(8).at(4)=3; permutation_PAM1.at(8).at(5)=4; permutation_PAM1.at(8).at(6)=5; permutation_PAM1.at(8).at(7)=6; permutation_PAM1.at(8).at(8)=7; permutation_PAM1.at(8).at(9)=9; permutation_PAM1.at(8).at(10)=10; permutation_PAM1.at(8).at(11)=11; permutation_PAM1.at(8).at(12)=12; permutation_PAM1.at(8).at(13)=13; permutation_PAM1.at(8).at(14)=14; permutation_PAM1.at(8).at(15)=15; permutation_PAM1.at(8).at(16)=16; permutation_PAM1.at(8).at(17)=17; permutation_PAM1.at(8).at(18)=18; permutation_PAM1.at(8).at(19)=19;
		permutation_PAM1.at(9).at(0)=9; permutation_PAM1.at(9).at(1)=0; permutation_PAM1.at(9).at(2)=1; permutation_PAM1.at(9).at(3)=2; permutation_PAM1.at(9).at(4)=3; permutation_PAM1.at(9).at(5)=4; permutation_PAM1.at(9).at(6)=5; permutation_PAM1.at(9).at(7)=6; permutation_PAM1.at(9).at(8)=7; permutation_PAM1.at(9).at(9)=8; permutation_PAM1.at(9).at(10)=10; permutation_PAM1.at(9).at(11)=11; permutation_PAM1.at(9).at(12)=12; permutation_PAM1.at(9).at(13)=13; permutation_PAM1.at(9).at(14)=14; permutation_PAM1.at(9).at(15)=15; permutation_PAM1.at(9).at(16)=16; permutation_PAM1.at(9).at(17)=17; permutation_PAM1.at(9).at(18)=18; permutation_PAM1.at(9).at(19)=19;
		permutation_PAM1.at(10).at(0)=10; permutation_PAM1.at(10).at(1)=0; permutation_PAM1.at(10).at(2)=1; permutation_PAM1.at(10).at(3)=2; permutation_PAM1.at(10).at(4)=3; permutation_PAM1.at(10).at(5)=4; permutation_PAM1.at(10).at(6)=5; permutation_PAM1.at(10).at(7)=6; permutation_PAM1.at(10).at(8)=7; permutation_PAM1.at(10).at(9)=8; permutation_PAM1.at(10).at(10)=9; permutation_PAM1.at(10).at(11)=11; permutation_PAM1.at(10).at(12)=12; permutation_PAM1.at(10).at(13)=13; permutation_PAM1.at(10).at(14)=14; permutation_PAM1.at(10).at(15)=15; permutation_PAM1.at(10).at(16)=16; permutation_PAM1.at(10).at(17)=17; permutation_PAM1.at(10).at(18)=18; permutation_PAM1.at(10).at(19)=19;
		permutation_PAM1.at(11).at(0)=11; permutation_PAM1.at(11).at(1)=0; permutation_PAM1.at(11).at(2)=1; permutation_PAM1.at(11).at(3)=2; permutation_PAM1.at(11).at(4)=3; permutation_PAM1.at(11).at(5)=4; permutation_PAM1.at(11).at(6)=5; permutation_PAM1.at(11).at(7)=6; permutation_PAM1.at(11).at(8)=7; permutation_PAM1.at(11).at(9)=8; permutation_PAM1.at(11).at(10)=9; permutation_PAM1.at(11).at(11)=10; permutation_PAM1.at(11).at(12)=12; permutation_PAM1.at(11).at(13)=13; permutation_PAM1.at(11).at(14)=14; permutation_PAM1.at(11).at(15)=15; permutation_PAM1.at(11).at(16)=16; permutation_PAM1.at(11).at(17)=17; permutation_PAM1.at(11).at(18)=18; permutation_PAM1.at(11).at(19)=19;
		permutation_PAM1.at(12).at(0)=12; permutation_PAM1.at(12).at(1)=0; permutation_PAM1.at(12).at(2)=1; permutation_PAM1.at(12).at(3)=2; permutation_PAM1.at(12).at(4)=3; permutation_PAM1.at(12).at(5)=4; permutation_PAM1.at(12).at(6)=5; permutation_PAM1.at(12).at(7)=6; permutation_PAM1.at(12).at(8)=7; permutation_PAM1.at(12).at(9)=8; permutation_PAM1.at(12).at(10)=9; permutation_PAM1.at(12).at(11)=10; permutation_PAM1.at(12).at(12)=11; permutation_PAM1.at(12).at(13)=13; permutation_PAM1.at(12).at(14)=14; permutation_PAM1.at(12).at(15)=15; permutation_PAM1.at(12).at(16)=16; permutation_PAM1.at(12).at(17)=17; permutation_PAM1.at(12).at(18)=18; permutation_PAM1.at(12).at(19)=19;
		permutation_PAM1.at(13).at(0)=13; permutation_PAM1.at(13).at(1)=0; permutation_PAM1.at(13).at(2)=1; permutation_PAM1.at(13).at(3)=2; permutation_PAM1.at(13).at(4)=3; permutation_PAM1.at(13).at(5)=4; permutation_PAM1.at(13).at(6)=5; permutation_PAM1.at(13).at(7)=6; permutation_PAM1.at(13).at(8)=7; permutation_PAM1.at(13).at(9)=8; permutation_PAM1.at(13).at(10)=9; permutation_PAM1.at(13).at(11)=10; permutation_PAM1.at(13).at(12)=11; permutation_PAM1.at(13).at(13)=12; permutation_PAM1.at(13).at(14)=14; permutation_PAM1.at(13).at(15)=15; permutation_PAM1.at(13).at(16)=16; permutation_PAM1.at(13).at(17)=17; permutation_PAM1.at(13).at(18)=18; permutation_PAM1.at(13).at(19)=19;
		permutation_PAM1.at(14).at(0)=14; permutation_PAM1.at(14).at(1)=0; permutation_PAM1.at(14).at(2)=1; permutation_PAM1.at(14).at(3)=2; permutation_PAM1.at(14).at(4)=3; permutation_PAM1.at(14).at(5)=4; permutation_PAM1.at(14).at(6)=5; permutation_PAM1.at(14).at(7)=6; permutation_PAM1.at(14).at(8)=7; permutation_PAM1.at(14).at(9)=8; permutation_PAM1.at(14).at(10)=9; permutation_PAM1.at(14).at(11)=10; permutation_PAM1.at(14).at(12)=11; permutation_PAM1.at(14).at(13)=12; permutation_PAM1.at(14).at(14)=13; permutation_PAM1.at(14).at(15)=15; permutation_PAM1.at(14).at(16)=16; permutation_PAM1.at(14).at(17)=17; permutation_PAM1.at(14).at(18)=18; permutation_PAM1.at(14).at(19)=19;
		permutation_PAM1.at(15).at(0)=15; permutation_PAM1.at(15).at(1)=0; permutation_PAM1.at(15).at(2)=1; permutation_PAM1.at(15).at(3)=2; permutation_PAM1.at(15).at(4)=3; permutation_PAM1.at(15).at(5)=4; permutation_PAM1.at(15).at(6)=5; permutation_PAM1.at(15).at(7)=6; permutation_PAM1.at(15).at(8)=7; permutation_PAM1.at(15).at(9)=8; permutation_PAM1.at(15).at(10)=9; permutation_PAM1.at(15).at(11)=10; permutation_PAM1.at(15).at(12)=11; permutation_PAM1.at(15).at(13)=12; permutation_PAM1.at(15).at(14)=13; permutation_PAM1.at(15).at(15)=14; permutation_PAM1.at(15).at(16)=16; permutation_PAM1.at(15).at(17)=17; permutation_PAM1.at(15).at(18)=18; permutation_PAM1.at(15).at(19)=19;
		permutation_PAM1.at(16).at(0)=16; permutation_PAM1.at(16).at(1)=0; permutation_PAM1.at(16).at(2)=1; permutation_PAM1.at(16).at(3)=2; permutation_PAM1.at(16).at(4)=3; permutation_PAM1.at(16).at(5)=4; permutation_PAM1.at(16).at(6)=5; permutation_PAM1.at(16).at(7)=6; permutation_PAM1.at(16).at(8)=7; permutation_PAM1.at(16).at(9)=8; permutation_PAM1.at(16).at(10)=9; permutation_PAM1.at(16).at(11)=10; permutation_PAM1.at(16).at(12)=11; permutation_PAM1.at(16).at(13)=12; permutation_PAM1.at(16).at(14)=13; permutation_PAM1.at(16).at(15)=14; permutation_PAM1.at(16).at(16)=15; permutation_PAM1.at(16).at(17)=17; permutation_PAM1.at(16).at(18)=18; permutation_PAM1.at(16).at(19)=19;
		permutation_PAM1.at(17).at(0)=17; permutation_PAM1.at(17).at(1)=0; permutation_PAM1.at(17).at(2)=1; permutation_PAM1.at(17).at(3)=2; permutation_PAM1.at(17).at(4)=3; permutation_PAM1.at(17).at(5)=4; permutation_PAM1.at(17).at(6)=5; permutation_PAM1.at(17).at(7)=6; permutation_PAM1.at(17).at(8)=7; permutation_PAM1.at(17).at(9)=8; permutation_PAM1.at(17).at(10)=9; permutation_PAM1.at(17).at(11)=10; permutation_PAM1.at(17).at(12)=11; permutation_PAM1.at(17).at(13)=12; permutation_PAM1.at(17).at(14)=13; permutation_PAM1.at(17).at(15)=14; permutation_PAM1.at(17).at(16)=15; permutation_PAM1.at(17).at(17)=16; permutation_PAM1.at(17).at(18)=18; permutation_PAM1.at(17).at(19)=19;
		permutation_PAM1.at(18).at(0)=18; permutation_PAM1.at(18).at(1)=0; permutation_PAM1.at(18).at(2)=1; permutation_PAM1.at(18).at(3)=2; permutation_PAM1.at(18).at(4)=3; permutation_PAM1.at(18).at(5)=4; permutation_PAM1.at(18).at(6)=5; permutation_PAM1.at(18).at(7)=6; permutation_PAM1.at(18).at(8)=7; permutation_PAM1.at(18).at(9)=8; permutation_PAM1.at(18).at(10)=9; permutation_PAM1.at(18).at(11)=10; permutation_PAM1.at(18).at(12)=11; permutation_PAM1.at(18).at(13)=12; permutation_PAM1.at(18).at(14)=13; permutation_PAM1.at(18).at(15)=14; permutation_PAM1.at(18).at(16)=15; permutation_PAM1.at(18).at(17)=16; permutation_PAM1.at(18).at(18)=17; permutation_PAM1.at(18).at(19)=19;
		permutation_PAM1.at(19).at(0)=19; permutation_PAM1.at(19).at(1)=0; permutation_PAM1.at(19).at(2)=1; permutation_PAM1.at(19).at(3)=2; permutation_PAM1.at(19).at(4)=3; permutation_PAM1.at(19).at(5)=4; permutation_PAM1.at(19).at(6)=5; permutation_PAM1.at(19).at(7)=6; permutation_PAM1.at(19).at(8)=7; permutation_PAM1.at(19).at(9)=8; permutation_PAM1.at(19).at(10)=9; permutation_PAM1.at(19).at(11)=10; permutation_PAM1.at(19).at(12)=11; permutation_PAM1.at(19).at(13)=12; permutation_PAM1.at(19).at(14)=13; permutation_PAM1.at(19).at(15)=14; permutation_PAM1.at(19).at(16)=15; permutation_PAM1.at(19).at(17)=16; permutation_PAM1.at(19).at(18)=17; permutation_PAM1.at(19).at(19)=18;

		scoring_matrix.at(1).at(1)=5.00; scoring_matrix.at(1).at(2)=-2.00; scoring_matrix.at(1).at(3)=0.00; scoring_matrix.at(1).at(4)=1.00; scoring_matrix.at(1).at(5)=-2.00; scoring_matrix.at(1).at(6)=0.00; scoring_matrix.at(1).at(7)=0.00; scoring_matrix.at(1).at(8)=-1.00; scoring_matrix.at(1).at(9)=0.00; scoring_matrix.at(1).at(10)=-1.00; scoring_matrix.at(1).at(11)=0.00; scoring_matrix.at(1).at(12)=0.00; scoring_matrix.at(1).at(13)=1.00; scoring_matrix.at(1).at(14)=0.00; scoring_matrix.at(1).at(15)=-1.00; scoring_matrix.at(1).at(16)=1.00; scoring_matrix.at(1).at(17)=0.00; scoring_matrix.at(1).at(18)=0.00; scoring_matrix.at(1).at(19)=-2.00; scoring_matrix.at(1).at(20)=-2.00;
		scoring_matrix.at(2).at(1)=-2.00; scoring_matrix.at(2).at(2)=6.00; scoring_matrix.at(2).at(3)=-2.00; scoring_matrix.at(2).at(4)=-3.00; scoring_matrix.at(2).at(5)=-3.00; scoring_matrix.at(2).at(6)=-2.00; scoring_matrix.at(2).at(7)=0.00; scoring_matrix.at(2).at(8)=-2.00; scoring_matrix.at(2).at(9)=-3.00; scoring_matrix.at(2).at(10)=-3.00; scoring_matrix.at(2).at(11)=0.00; scoring_matrix.at(2).at(12)=-2.00; scoring_matrix.at(2).at(13)=-3.00; scoring_matrix.at(2).at(14)=-3.00; scoring_matrix.at(2).at(15)=-2.00; scoring_matrix.at(2).at(16)=-1.00; scoring_matrix.at(2).at(17)=-1.00; scoring_matrix.at(2).at(18)=-2.00; scoring_matrix.at(2).at(19)=-1.00; scoring_matrix.at(2).at(20)=-2.00;
		scoring_matrix.at(3).at(1)=0.00; scoring_matrix.at(3).at(2)=-2.00; scoring_matrix.at(3).at(3)=5.00; scoring_matrix.at(3).at(4)=2.00; scoring_matrix.at(3).at(5)=-2.00; scoring_matrix.at(3).at(6)=0.00; scoring_matrix.at(3).at(7)=1.00; scoring_matrix.at(3).at(8)=-3.00; scoring_matrix.at(3).at(9)=0.00; scoring_matrix.at(3).at(10)=-2.00; scoring_matrix.at(3).at(11)=-1.00; scoring_matrix.at(3).at(12)=2.00; scoring_matrix.at(3).at(13)=0.00; scoring_matrix.at(3).at(14)=1.00; scoring_matrix.at(3).at(15)=-2.00; scoring_matrix.at(3).at(16)=0.00; scoring_matrix.at(3).at(17)=0.00; scoring_matrix.at(3).at(18)=-2.00; scoring_matrix.at(3).at(19)=-3.00; scoring_matrix.at(3).at(20)=-2.00;
		scoring_matrix.at(4).at(1)=1.00; scoring_matrix.at(4).at(2)=-3.00; scoring_matrix.at(4).at(3)=2.00; scoring_matrix.at(4).at(4)=5.00; scoring_matrix.at(4).at(5)=-3.00; scoring_matrix.at(4).at(6)=0.00; scoring_matrix.at(4).at(7)=-1.00; scoring_matrix.at(4).at(8)=-2.00; scoring_matrix.at(4).at(9)=1.00; scoring_matrix.at(4).at(10)=-2.00; scoring_matrix.at(4).at(11)=-2.00; scoring_matrix.at(4).at(12)=1.00; scoring_matrix.at(4).at(13)=1.00; scoring_matrix.at(4).at(14)=2.00; scoring_matrix.at(4).at(15)=0.00; scoring_matrix.at(4).at(16)=1.00; scoring_matrix.at(4).at(17)=1.00; scoring_matrix.at(4).at(18)=-1.00; scoring_matrix.at(4).at(19)=-2.00; scoring_matrix.at(4).at(20)=-1.00;
		scoring_matrix.at(5).at(1)=-2.00; scoring_matrix.at(5).at(2)=-3.00; scoring_matrix.at(5).at(3)=-2.00; scoring_matrix.at(5).at(4)=-3.00; scoring_matrix.at(5).at(5)=6.00; scoring_matrix.at(5).at(6)=-3.00; scoring_matrix.at(5).at(7)=1.00; scoring_matrix.at(5).at(8)=0.00; scoring_matrix.at(5).at(9)=-3.00; scoring_matrix.at(5).at(10)=2.00; scoring_matrix.at(5).at(11)=2.00; scoring_matrix.at(5).at(12)=-3.00; scoring_matrix.at(5).at(13)=-2.00; scoring_matrix.at(5).at(14)=-3.00; scoring_matrix.at(5).at(15)=-2.00; scoring_matrix.at(5).at(16)=-1.00; scoring_matrix.at(5).at(17)=-2.00; scoring_matrix.at(5).at(18)=0.00; scoring_matrix.at(5).at(19)=3.00; scoring_matrix.at(5).at(20)=3.00;
		scoring_matrix.at(6).at(1)=0.00; scoring_matrix.at(6).at(2)=-2.00; scoring_matrix.at(6).at(3)=0.00; scoring_matrix.at(6).at(4)=0.00; scoring_matrix.at(6).at(5)=-3.00; scoring_matrix.at(6).at(6)=5.00; scoring_matrix.at(6).at(7)=-1.00; scoring_matrix.at(6).at(8)=-2.00; scoring_matrix.at(6).at(9)=0.00; scoring_matrix.at(6).at(10)=-2.00; scoring_matrix.at(6).at(11)=-2.00; scoring_matrix.at(6).at(12)=0.00; scoring_matrix.at(6).at(13)=0.00; scoring_matrix.at(6).at(14)=-1.00; scoring_matrix.at(6).at(15)=0.00; scoring_matrix.at(6).at(16)=0.00; scoring_matrix.at(6).at(17)=-1.00; scoring_matrix.at(6).at(18)=-1.00; scoring_matrix.at(6).at(19)=-2.00; scoring_matrix.at(6).at(20)=-3.00;
		scoring_matrix.at(7).at(1)=0.00; scoring_matrix.at(7).at(2)=0.00; scoring_matrix.at(7).at(3)=1.00; scoring_matrix.at(7).at(4)=-1.00; scoring_matrix.at(7).at(5)=1.00; scoring_matrix.at(7).at(6)=-1.00; scoring_matrix.at(7).at(7)=5.00; scoring_matrix.at(7).at(8)=-1.00; scoring_matrix.at(7).at(9)=1.00; scoring_matrix.at(7).at(10)=-1.00; scoring_matrix.at(7).at(11)=0.00; scoring_matrix.at(7).at(12)=1.00; scoring_matrix.at(7).at(13)=0.00; scoring_matrix.at(7).at(14)=1.00; scoring_matrix.at(7).at(15)=2.00; scoring_matrix.at(7).at(16)=0.00; scoring_matrix.at(7).at(17)=1.00; scoring_matrix.at(7).at(18)=-1.00; scoring_matrix.at(7).at(19)=0.00; scoring_matrix.at(7).at(20)=1.00;
		scoring_matrix.at(8).at(1)=-1.00; scoring_matrix.at(8).at(2)=-2.00; scoring_matrix.at(8).at(3)=-3.00; scoring_matrix.at(8).at(4)=-2.00; scoring_matrix.at(8).at(5)=0.00; scoring_matrix.at(8).at(6)=-2.00; scoring_matrix.at(8).at(7)=-1.00; scoring_matrix.at(8).at(8)=5.00; scoring_matrix.at(8).at(9)=-2.00; scoring_matrix.at(8).at(10)=2.00; scoring_matrix.at(8).at(11)=2.00; scoring_matrix.at(8).at(12)=-2.00; scoring_matrix.at(8).at(13)=-2.00; scoring_matrix.at(8).at(14)=-3.00; scoring_matrix.at(8).at(15)=-2.00; scoring_matrix.at(8).at(16)=-1.00; scoring_matrix.at(8).at(17)=0.00; scoring_matrix.at(8).at(18)=2.00; scoring_matrix.at(8).at(19)=0.00; scoring_matrix.at(8).at(20)=0.00;
		scoring_matrix.at(9).at(1)=0.00; scoring_matrix.at(9).at(2)=-3.00; scoring_matrix.at(9).at(3)=0.00; scoring_matrix.at(9).at(4)=1.00; scoring_matrix.at(9).at(5)=-3.00; scoring_matrix.at(9).at(6)=0.00; scoring_matrix.at(9).at(7)=1.00; scoring_matrix.at(9).at(8)=-2.00; scoring_matrix.at(9).at(9)=5.00; scoring_matrix.at(9).at(10)=-1.00; scoring_matrix.at(9).at(11)=-2.00; scoring_matrix.at(9).at(12)=1.00; scoring_matrix.at(9).at(13)=0.00; scoring_matrix.at(9).at(14)=1.00; scoring_matrix.at(9).at(15)=2.00; scoring_matrix.at(9).at(16)=0.00; scoring_matrix.at(9).at(17)=0.00; scoring_matrix.at(9).at(18)=-1.00; scoring_matrix.at(9).at(19)=-2.00; scoring_matrix.at(9).at(20)=-2.00;
		scoring_matrix.at(10).at(1)=-1.00; scoring_matrix.at(10).at(2)=-3.00; scoring_matrix.at(10).at(3)=-2.00; scoring_matrix.at(10).at(4)=-2.00; scoring_matrix.at(10).at(5)=2.00; scoring_matrix.at(10).at(6)=-2.00; scoring_matrix.at(10).at(7)=-1.00; scoring_matrix.at(10).at(8)=2.00; scoring_matrix.at(10).at(9)=-1.00; scoring_matrix.at(10).at(10)=5.00; scoring_matrix.at(10).at(11)=3.00; scoring_matrix.at(10).at(12)=-2.00; scoring_matrix.at(10).at(13)=-2.00; scoring_matrix.at(10).at(14)=0.00; scoring_matrix.at(10).at(15)=-1.00; scoring_matrix.at(10).at(16)=-1.00; scoring_matrix.at(10).at(17)=0.00; scoring_matrix.at(10).at(18)=2.00; scoring_matrix.at(10).at(19)=0.00; scoring_matrix.at(10).at(20)=0.00;
		scoring_matrix.at(11).at(1)=0.00; scoring_matrix.at(11).at(2)=0.00; scoring_matrix.at(11).at(3)=-1.00; scoring_matrix.at(11).at(4)=-2.00; scoring_matrix.at(11).at(5)=2.00; scoring_matrix.at(11).at(6)=-2.00; scoring_matrix.at(11).at(7)=0.00; scoring_matrix.at(11).at(8)=2.00; scoring_matrix.at(11).at(9)=-2.00; scoring_matrix.at(11).at(10)=3.00; scoring_matrix.at(11).at(11)=5.00; scoring_matrix.at(11).at(12)=-1.00; scoring_matrix.at(11).at(13)=-2.00; scoring_matrix.at(11).at(14)=0.00; scoring_matrix.at(11).at(15)=-2.00; scoring_matrix.at(11).at(16)=-1.00; scoring_matrix.at(11).at(17)=0.00; scoring_matrix.at(11).at(18)=1.00; scoring_matrix.at(11).at(19)=-2.00; scoring_matrix.at(11).at(20)=-1.00;
		scoring_matrix.at(12).at(1)=0.00; scoring_matrix.at(12).at(2)=-2.00; scoring_matrix.at(12).at(3)=2.00; scoring_matrix.at(12).at(4)=1.00; scoring_matrix.at(12).at(5)=-3.00; scoring_matrix.at(12).at(6)=0.00; scoring_matrix.at(12).at(7)=1.00; scoring_matrix.at(12).at(8)=-2.00; scoring_matrix.at(12).at(9)=1.00; scoring_matrix.at(12).at(10)=-2.00; scoring_matrix.at(12).at(11)=-1.00; scoring_matrix.at(12).at(12)=5.00; scoring_matrix.at(12).at(13)=-2.00; scoring_matrix.at(12).at(14)=1.00; scoring_matrix.at(12).at(15)=0.00; scoring_matrix.at(12).at(16)=2.00; scoring_matrix.at(12).at(17)=0.00; scoring_matrix.at(12).at(18)=-2.00; scoring_matrix.at(12).at(19)=-3.00; scoring_matrix.at(12).at(20)=-1.00;
		scoring_matrix.at(13).at(1)=1.00; scoring_matrix.at(13).at(2)=-3.00; scoring_matrix.at(13).at(3)=0.00; scoring_matrix.at(13).at(4)=1.00; scoring_matrix.at(13).at(5)=-2.00; scoring_matrix.at(13).at(6)=0.00; scoring_matrix.at(13).at(7)=0.00; scoring_matrix.at(13).at(8)=-2.00; scoring_matrix.at(13).at(9)=0.00; scoring_matrix.at(13).at(10)=-2.00; scoring_matrix.at(13).at(11)=-2.00; scoring_matrix.at(13).at(12)=-2.00; scoring_matrix.at(13).at(13)=5.00; scoring_matrix.at(13).at(14)=0.00; scoring_matrix.at(13).at(15)=0.00; scoring_matrix.at(13).at(16)=0.00; scoring_matrix.at(13).at(17)=0.00; scoring_matrix.at(13).at(18)=-1.00; scoring_matrix.at(13).at(19)=-3.00; scoring_matrix.at(13).at(20)=-3.00;
		scoring_matrix.at(14).at(1)=0.00; scoring_matrix.at(14).at(2)=-3.00; scoring_matrix.at(14).at(3)=1.00; scoring_matrix.at(14).at(4)=2.00; scoring_matrix.at(14).at(5)=-3.00; scoring_matrix.at(14).at(6)=-1.00; scoring_matrix.at(14).at(7)=1.00; scoring_matrix.at(14).at(8)=-3.00; scoring_matrix.at(14).at(9)=1.00; scoring_matrix.at(14).at(10)=0.00; scoring_matrix.at(14).at(11)=0.00; scoring_matrix.at(14).at(12)=1.00; scoring_matrix.at(14).at(13)=0.00; scoring_matrix.at(14).at(14)=5.00; scoring_matrix.at(14).at(15)=2.00; scoring_matrix.at(14).at(16)=1.00; scoring_matrix.at(14).at(17)=0.00; scoring_matrix.at(14).at(18)=-1.00; scoring_matrix.at(14).at(19)=-1.00; scoring_matrix.at(14).at(20)=-2.00;
		scoring_matrix.at(15).at(1)=-1.00; scoring_matrix.at(15).at(2)=-2.00; scoring_matrix.at(15).at(3)=-2.00; scoring_matrix.at(15).at(4)=0.00; scoring_matrix.at(15).at(5)=-2.00; scoring_matrix.at(15).at(6)=0.00; scoring_matrix.at(15).at(7)=2.00; scoring_matrix.at(15).at(8)=-2.00; scoring_matrix.at(15).at(9)=2.00; scoring_matrix.at(15).at(10)=-1.00; scoring_matrix.at(15).at(11)=-2.00; scoring_matrix.at(15).at(12)=0.00; scoring_matrix.at(15).at(13)=0.00; scoring_matrix.at(15).at(14)=2.00; scoring_matrix.at(15).at(15)=5.00; scoring_matrix.at(15).at(16)=1.00; scoring_matrix.at(15).at(17)=0.00; scoring_matrix.at(15).at(18)=-1.00; scoring_matrix.at(15).at(19)=0.00; scoring_matrix.at(15).at(20)=-1.00;
		scoring_matrix.at(16).at(1)=1.00; scoring_matrix.at(16).at(2)=-1.00; scoring_matrix.at(16).at(3)=0.00; scoring_matrix.at(16).at(4)=1.00; scoring_matrix.at(16).at(5)=-1.00; scoring_matrix.at(16).at(6)=0.00; scoring_matrix.at(16).at(7)=0.00; scoring_matrix.at(16).at(8)=-1.00; scoring_matrix.at(16).at(9)=0.00; scoring_matrix.at(16).at(10)=-1.00; scoring_matrix.at(16).at(11)=-1.00; scoring_matrix.at(16).at(12)=2.00; scoring_matrix.at(16).at(13)=0.00; scoring_matrix.at(16).at(14)=1.00; scoring_matrix.at(16).at(15)=1.00; scoring_matrix.at(16).at(16)=5.00; scoring_matrix.at(16).at(17)=2.00; scoring_matrix.at(16).at(18)=-1.00; scoring_matrix.at(16).at(19)=0.00; scoring_matrix.at(16).at(20)=0.00;
		scoring_matrix.at(17).at(1)=0.00; scoring_matrix.at(17).at(2)=-1.00; scoring_matrix.at(17).at(3)=0.00; scoring_matrix.at(17).at(4)=1.00; scoring_matrix.at(17).at(5)=-2.00; scoring_matrix.at(17).at(6)=-1.00; scoring_matrix.at(17).at(7)=1.00; scoring_matrix.at(17).at(8)=0.00; scoring_matrix.at(17).at(9)=0.00; scoring_matrix.at(17).at(10)=0.00; scoring_matrix.at(17).at(11)=0.00; scoring_matrix.at(17).at(12)=0.00; scoring_matrix.at(17).at(13)=0.00; scoring_matrix.at(17).at(14)=0.00; scoring_matrix.at(17).at(15)=0.00; scoring_matrix.at(17).at(16)=2.00; scoring_matrix.at(17).at(17)=5.00; scoring_matrix.at(17).at(18)=0.00; scoring_matrix.at(17).at(19)=-1.00; scoring_matrix.at(17).at(20)=-2.00;
		scoring_matrix.at(18).at(1)=0.00; scoring_matrix.at(18).at(2)=-2.00; scoring_matrix.at(18).at(3)=-2.00; scoring_matrix.at(18).at(4)=-1.00; scoring_matrix.at(18).at(5)=0.00; scoring_matrix.at(18).at(6)=-1.00; scoring_matrix.at(18).at(7)=-1.00; scoring_matrix.at(18).at(8)=2.00; scoring_matrix.at(18).at(9)=-1.00; scoring_matrix.at(18).at(10)=2.00; scoring_matrix.at(18).at(11)=1.00; scoring_matrix.at(18).at(12)=-2.00; scoring_matrix.at(18).at(13)=-1.00; scoring_matrix.at(18).at(14)=-1.00; scoring_matrix.at(18).at(15)=-1.00; scoring_matrix.at(18).at(16)=-1.00; scoring_matrix.at(18).at(17)=0.00; scoring_matrix.at(18).at(18)=5.00; scoring_matrix.at(18).at(19)=-1.00; scoring_matrix.at(18).at(20)=0.00;
		scoring_matrix.at(19).at(1)=-2.00; scoring_matrix.at(19).at(2)=-1.00; scoring_matrix.at(19).at(3)=-3.00; scoring_matrix.at(19).at(4)=-2.00; scoring_matrix.at(19).at(5)=3.00; scoring_matrix.at(19).at(6)=-2.00; scoring_matrix.at(19).at(7)=0.00; scoring_matrix.at(19).at(8)=0.00; scoring_matrix.at(19).at(9)=-2.00; scoring_matrix.at(19).at(10)=0.00; scoring_matrix.at(19).at(11)=-2.00; scoring_matrix.at(19).at(12)=-3.00; scoring_matrix.at(19).at(13)=-3.00; scoring_matrix.at(19).at(14)=-1.00; scoring_matrix.at(19).at(15)=0.00; scoring_matrix.at(19).at(16)=0.00; scoring_matrix.at(19).at(17)=-1.00; scoring_matrix.at(19).at(18)=-1.00; scoring_matrix.at(19).at(19)=6.00; scoring_matrix.at(19).at(20)=3.00;
		scoring_matrix.at(20).at(1)=-2.00; scoring_matrix.at(20).at(2)=-2.00; scoring_matrix.at(20).at(3)=-2.00; scoring_matrix.at(20).at(4)=-1.00; scoring_matrix.at(20).at(5)=3.00; scoring_matrix.at(20).at(6)=-3.00; scoring_matrix.at(20).at(7)=1.00; scoring_matrix.at(20).at(8)=0.00; scoring_matrix.at(20).at(9)=-2.00; scoring_matrix.at(20).at(10)=0.00; scoring_matrix.at(20).at(11)=-1.00; scoring_matrix.at(20).at(12)=-1.00; scoring_matrix.at(20).at(13)=-3.00; scoring_matrix.at(20).at(14)=-2.00; scoring_matrix.at(20).at(15)=-1.00; scoring_matrix.at(20).at(16)=0.00; scoring_matrix.at(20).at(17)=-2.00; scoring_matrix.at(20).at(18)=0.00; scoring_matrix.at(20).at(19)=3.00; scoring_matrix.at(20).at(20)=6.00;
																						
		clear();
		srand ( time(NULL) );
	}
	
	
	void get_from_file (string& MSA_file_name) {
		string line;
		ifstream MSA_file;
		int line_count = 0;
		int line_length;
		char replace;
		
		MSA_file.open (MSA_file_name.c_str());
		/*if (!sequences.empty()) {
			cout << "Are you sure you wish to replace the currently loaded MSA? [Y/N, default N] ";
		cin >> replace;
		if (replace != 'N' || replace != 'n')
			return;
		}*/
		
		if (MSA_file.is_open()) {
			if (! MSA_file.eof()) {
				getline (MSA_file, line,';');
				line_count++;
				line_length = line.length();
				getline (MSA_file, line,';');
				while (! MSA_file.eof()) {
					if (line.length() == line_length)
						line_count++;
					else {
						cout << "Operation Aborted: All sequences must be of the same length!\n";
						return;
					}
					getline (MSA_file, line,';');  //reads next sequence from file
				}
				MSA_file.close();
				MSA_file.open(MSA_file_name.c_str());
				MSA_file.seekg(0,ios::beg);
				clear();
				number_of_sequences = line_count;
				alignment_length = line_length;
				getline(MSA_file, line,';');
				while (! MSA_file.eof()) {
					if (!line.length()<4){
						sequences.push_back(line);
					}
					getline(MSA_file, line,';');
					
				}
			}
			else
				cout << "Operation Aborted: File \"" << MSA_file_name << "\" is empty!\n";
		}
		else
			cout << "Operation Aborted: Unable to open file \"" << MSA_file_name << "\"!\n";
		return;
	}
	
	void perspective (int reference_sequence_number) {
		int reference_sequence_length=0;
		
		if (aa_counts.empty()) {
			tabulate ();
		}
		MSA_to_reference.resize(alignment_length,0);
		reference_to_MSA.resize(alignment_length,0);
		for (int i=0; i<alignment_length; i++) {
			if (sequences[reference_sequence_number - 1][i] != '-') {
				reference_sequence_length = reference_sequence_length + 1;
				MSA_to_reference.at(i)=reference_sequence_length - 1;
				reference_to_MSA.at(reference_sequence_length - 1)=i;
			} else {
				MSA_to_reference.at(i)= -1;
			}
		}
		reference_to_MSA.resize(reference_sequence_length);
		return;
	}
	
	void read_secondary_structure (string& PDB_file_name, string& target_chain, int offset) {
		ifstream PDB_file;
		string line, record_type, chain_id, post_structure_id;
		
		istringstream substring (stringstream::in);
		int beginning, ending, helix_class;
		int number_helices=0;
		int number_sheets=0;
		int number_turns=0;
		int start_position=0;
		int end_position=0;
		
		PDB_file.open (PDB_file_name.c_str());
		helices.resize(10, vector <int> (reference_to_MSA.size(), 0));
		sheets.resize(reference_to_MSA.size(), 0);
		turns.resize(reference_to_MSA.size(), 0);
		while ((record_type!="HELIX ")&&(record_type!="SHEET ")&&(record_type!="TURN  ")) {
			getline (PDB_file, line);
			record_type.assign(line,0,6);
			cout << record_type;
		}
		while ((record_type=="HELIX ")||(record_type=="SHEET ")||(record_type=="TURN  ")) {
			if (record_type=="HELIX ") {
				chain_id.assign(line,19,1);
			} if (record_type=="SHEET ") {
				chain_id.assign(line,21,1);
			} if (record_type=="TURN ") {
				chain_id.assign(line,19,1);
			}
			if ((chain_id==target_chain)||(target_chain==" "&&((chain_id=="A")||(chain_id=="a")))) {
				if (record_type=="HELIX ") {
					substring.clear();
					substring.str(line.substr(21,4));
					substring >> skipws >> start_position;
					if (start_position - offset >= 1) {
						beginning = start_position - offset - 1;
					} else {
						beginning = 0;
					}
					substring.clear();
					substring.str(line.substr(33,4));
					substring >> skipws >> end_position;
					if (end_position - offset <= reference_to_MSA.size()) {
						ending = end_position - offset - 1;
					} else {
						ending = reference_to_MSA.size() - 1;
					}
					if ((beginning != 0)||(ending!=reference_to_MSA.size() - 1)||((start_position - offset == 1)&&(end_position - offset == reference_to_MSA.size()))) {
						number_helices=number_helices + 1;
						if (line.substr(38,2)=="  "){
							helix_class=1;
						} else {
							substring.clear();
							substring.str(line.substr(38,2));
							substring >> skipws >> helix_class;
						}
						for (int i=beginning; i <= ending; i++) {
							helices.at(helix_class - 1).at(i)=number_helices;
						}
					}
				} else if (record_type=="SHEET ") {
					substring.clear();
					substring.str(line.substr(22,4));
					substring >> skipws >> start_position;
					if (start_position - offset >= 1) {
						beginning = start_position - offset - 1;
					} else {
						beginning = 0;
					}
					substring.clear();
					substring.str(line.substr(33,4));
					substring >> skipws >> end_position;
					if (end_position - offset <= reference_to_MSA.size()) {
						ending = end_position - offset - 1;
					} else {
						ending = reference_to_MSA.size() - 1;
					}
					if ((beginning != 0)||(ending!=reference_to_MSA.size() - 1)||((start_position - offset == 1)&&(end_position - offset == reference_to_MSA.size()))) {
						number_sheets=number_sheets + 1;
						for (int i=beginning; i <= ending; i++) {
							sheets.at(i)=number_sheets;
						}
					}
				} else if (record_type=="TURN  ") {
					substring.clear();
					substring.str(line.substr(20,4));
					substring >> skipws >> start_position;
					if (start_position - offset >= 1) {
						beginning = start_position - offset - 1;
					} else {
						beginning = 0;
					}
					substring.clear();
					substring.str(line.substr(31,4));
					substring >> skipws >> end_position;
					if (end_position - offset <= reference_to_MSA.size()) {
						ending = end_position - offset - 1;
					} else {
						ending = reference_to_MSA.size() - 1;
					}
					if ((beginning != 0)||(ending!=reference_to_MSA.size() - 1)||((start_position - offset == 1)&&(end_position - offset == reference_to_MSA.size()))) {
						number_turns=number_turns + 1;
						for (int i=beginning; i <= ending; i++) {
							turns.at(i)=number_turns;
						}
					}
				}
			}
			getline (PDB_file, line);
			record_type.assign(line,0,6);
		}
	}			
	
	void read_PDB_file (string& PDB_file_name, int offset) {
		ifstream PDB_file;
		string line,record_type, atom,residue, variance;
		float x=0;
		float y=0;
		float z=0;
		float x2=0;
		float y2=0;
		float z2=0;
		istringstream substring (stringstream::in);
		int read_position=0;
		int read_model=1;
		int residue_position=0;
		int has_variants=0;
		
		
		atomic_coordinates.resize(reference_to_MSA.size(),vector<vector<float> >(100,vector<float>(3,0.0/0.0)));
		C_beta_index.resize(reference_to_MSA.size(),0);
		number_residue_atoms.resize(reference_to_MSA.size(),-1);
		minimum_distances.resize(reference_to_MSA.size(),vector<float> (reference_to_MSA.size(),0.0/0.0));
		residue_distances.resize(reference_to_MSA.size(),vector<float> (reference_to_MSA.size(),0.0/0.0));
		C_alpha_distances.resize(reference_to_MSA.size(),vector<float> (reference_to_MSA.size(),0.0/0.0));
		C_beta_distances.resize(reference_to_MSA.size(),vector<float> (reference_to_MSA.size(),0.0/0.0));
		
		
		PDB_file.open (PDB_file_name.c_str());
		
		while (record_type!="ATOM  ") {
			getline (PDB_file, line);
			record_type.assign(line,0,6);
		}
		while (record_type!="TER   ") {
			if (record_type=="ATOM  ") {
				substring.clear();
				substring.str(line.substr(22,4));
				substring >> read_position;
				if ((read_position - offset >= 1) && (read_position - offset <= reference_to_MSA.size())) {
					substring.clear();
					substring.str(line.substr(30,8));
					substring >> skipws >> x;
					substring.clear();
					substring.str(line.substr(38,8));
					substring >> skipws >> y;
					substring.clear();
					substring.str(line.substr(46,8));
					substring >> skipws >> z;
					atom.assign(line,12,4);
					residue.assign(line,19,3);
					variance.assign(line,16,1);
					if (variance != " ") {
						has_variants=1;
					}
					
					
					//cout << read_position;
					
					if (read_position!=residue_position) {
						residue_position=read_position;
						if (atom == " N  ") {
							number_residue_atoms.at(residue_position - offset - 1) = 0;
							atomic_coordinates.at(residue_position - offset - 1).at(0).at(0)=x;
							atomic_coordinates.at(residue_position - offset - 1).at(0).at(1)=y;
							atomic_coordinates.at(residue_position - offset - 1).at(0).at(2)=z;	
							for (int i=0; i<residue_position - offset - 1; i++) {
								if (number_residue_atoms.at(i)!=-1) {
									x2=atomic_coordinates.at(i).at(0).at(0);
									y2=atomic_coordinates.at(i).at(0).at(1);
									z2=atomic_coordinates.at(i).at(0).at(2);
									minimum_distances.at(i).at(residue_position - offset - 1) = float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0)));
									for (int n=1; n<=3 + number_residue_atoms.at(i); n++) {
										x2=atomic_coordinates.at(i).at(n).at(0);
										y2=atomic_coordinates.at(i).at(n).at(1);
										z2=atomic_coordinates.at(i).at(n).at(2);
										minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
									}
									minimum_distances.at(residue_position - offset - 1).at(i)=minimum_distances.at(i).at(residue_position - offset - 1);
								}
							}
						} else {
							cout << "WARNING! This is a bad PDB file! Improper atom ordering!\n";
						}
					} else {
						if (atom == "CA ") {
							atomic_coordinates.at(residue_position - offset - 1).at(1).at(0)=x;
							atomic_coordinates.at(residue_position - offset - 1).at(1).at(1)=y;
							atomic_coordinates.at(residue_position - offset - 1).at(1).at(2)=z;
							for (int i=0; i<residue_position - offset - 1; i++) {
								if (number_residue_atoms.at(i)!=-1) {
									x2=atomic_coordinates.at(i).at(0).at(0);
									y2=atomic_coordinates.at(i).at(0).at(1);
									z2=atomic_coordinates.at(i).at(0).at(2);
									minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
									
									x2=atomic_coordinates.at(i).at(1).at(0);
									y2=atomic_coordinates.at(i).at(1).at(1);
									z2=atomic_coordinates.at(i).at(1).at(2);
									minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
									C_alpha_distances.at(i).at(residue_position - offset - 1) = float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0)));										
									C_alpha_distances.at(residue_position - offset - 1).at(i) = C_alpha_distances.at(i).at(residue_position - offset - 1);
									
									for (int n=2; n<=3 + number_residue_atoms.at(i); n++) {
										x2=atomic_coordinates.at(i).at(n).at(0);
										y2=atomic_coordinates.at(i).at(n).at(1);
										z2=atomic_coordinates.at(i).at(n).at(2);
										minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
									}
									minimum_distances.at(residue_position - offset - 1).at(i)=minimum_distances.at(i).at(residue_position - offset - 1);
								}
							}
							if (residue=="GLY") {
								C_beta_index.at(residue_position - offset - 1)=1;
								for (int i=0; i<residue_position - offset - 1; i++) {
									if (number_residue_atoms.at(i)!=-1) {
										x2=atomic_coordinates.at(i).at(C_beta_index.at(i)).at(0);
										y2=atomic_coordinates.at(i).at(C_beta_index.at(i)).at(1);
										z2=atomic_coordinates.at(i).at(C_beta_index.at(i)).at(2);
										C_beta_distances.at(i).at(residue_position - offset - 1) = float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y-y2),2.0) + pow(fabs(z-z2),2.0)));
										C_beta_distances.at(residue_position - offset - 1).at(i) = C_beta_distances.at(i).at(residue_position - offset - 1);
									}
								}									
							}
							
						} else if (atom == " C  ") {
							
							atomic_coordinates.at(residue_position - offset - 1).at(2).at(0)=x;
							atomic_coordinates.at(residue_position - offset - 1).at(2).at(1)=y;
							atomic_coordinates.at(residue_position - offset - 1).at(2).at(2)=z;
							for (int i=0; i<residue_position - offset - 1; i++) {
								if (number_residue_atoms.at(i)!=-1) {
									x2=atomic_coordinates.at(i).at(0).at(0);
									y2=atomic_coordinates.at(i).at(0).at(1);
									z2=atomic_coordinates.at(i).at(0).at(2);
									minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
									
									for (int n=1; n<=3 + number_residue_atoms.at(i); n++) {
										x2=atomic_coordinates.at(i).at(n).at(0);
										y2=atomic_coordinates.at(i).at(n).at(1);
										z2=atomic_coordinates.at(i).at(n).at(2);
										minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
									}
									minimum_distances.at(residue_position - offset - 1).at(i)=minimum_distances.at(i).at(residue_position - offset - 1);
								}
							}
						} else if (atom == " O  ") {
							
							atomic_coordinates.at(residue_position - offset - 1).at(3).at(0)=x;
							atomic_coordinates.at(residue_position - offset - 1).at(3).at(1)=y;
							atomic_coordinates.at(residue_position - offset - 1).at(3).at(2)=z;
							for (int i=0; i<residue_position - offset - 1; i++) {
								if (number_residue_atoms.at(i)!=-1) {
									x2=atomic_coordinates.at(i).at(0).at(0);
									y2=atomic_coordinates.at(i).at(0).at(1);
									z2=atomic_coordinates.at(i).at(0).at(2);
									minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
									
									for (int n=1; n<=3 + number_residue_atoms.at(i); n++) {
										x2=atomic_coordinates.at(i).at(n).at(0);
										y2=atomic_coordinates.at(i).at(n).at(1);
										z2=atomic_coordinates.at(i).at(n).at(2);
										minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
									}
									minimum_distances.at(residue_position - offset - 1).at(i)=minimum_distances.at(i).at(residue_position - offset - 1);
								}
							}
						} else if (atom == " CB ") {
							
							number_residue_atoms.at(residue_position - offset - 1) = number_residue_atoms.at(residue_position - offset - 1) + 1;
							C_beta_index.at(residue_position - offset - 1)=number_residue_atoms.at(residue_position - offset - 1) + 3;
							atomic_coordinates.at(residue_position - offset - 1).at(number_residue_atoms.at(residue_position - offset - 1) + 3).at(0)=x;
							atomic_coordinates.at(residue_position - offset - 1).at(number_residue_atoms.at(residue_position - offset - 1) + 3).at(1)=y;
							atomic_coordinates.at(residue_position - offset - 1).at(number_residue_atoms.at(residue_position - offset - 1) + 3).at(2)=z;
							if (number_residue_atoms.at(residue_position - offset - 1) == 1) {
								for (int i=0; i<residue_position - offset - 1; i++) {
									if (number_residue_atoms.at(i)!=-1) {
										x2=atomic_coordinates.at(i).at(C_beta_index.at(i)).at(0);
										y2=atomic_coordinates.at(i).at(C_beta_index.at(i)).at(1);
										z2=atomic_coordinates.at(i).at(C_beta_index.at(i)).at(2);
										C_beta_distances.at(i).at(residue_position - offset - 1) = float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y-y2),2.0) + pow(fabs(z-z2),2.0)));
										C_beta_distances.at(residue_position - offset - 1).at(i) = C_beta_distances.at(i).at(residue_position - offset - 1);
										for (int n=0; n<=3; n++) {
											x2=atomic_coordinates.at(i).at(n).at(0);
											y2=atomic_coordinates.at(i).at(n).at(1);
											z2=atomic_coordinates.at(i).at(n).at(2);
											minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
										}
										if (number_residue_atoms.at(i)>=1) {
											x2=atomic_coordinates.at(i).at(4).at(0);
											y2=atomic_coordinates.at(i).at(4).at(1);
											z2=atomic_coordinates.at(i).at(4).at(2);
											residue_distances.at(i).at(residue_position - offset - 1) = float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0)));
											minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
											
											for (int n=5; n<=3 + number_residue_atoms.at(i); n++) {
												x2=atomic_coordinates.at(i).at(n).at(0);
												y2=atomic_coordinates.at(i).at(n).at(1);
												z2=atomic_coordinates.at(i).at(n).at(2);
												residue_distances.at(i).at(residue_position - offset - 1) = min(residue_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
												minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
											}
										}
										residue_distances.at(residue_position - offset - 1).at(i) = residue_distances.at(i).at(residue_position - offset - 1);
										minimum_distances.at(residue_position - offset - 1).at(i) = minimum_distances.at(i).at(residue_position - offset - 1);
									}
								}
							} else {
								for (int i=0; i<residue_position - offset - 1; i++) {
									if (number_residue_atoms.at(i)!=-1) {
										x2=atomic_coordinates.at(i).at(C_beta_index.at(i)).at(0);
										y2=atomic_coordinates.at(i).at(C_beta_index.at(i)).at(1);
										z2=atomic_coordinates.at(i).at(C_beta_index.at(i)).at(2);
										C_beta_distances.at(i).at(residue_position - offset - 1) = float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y-y2),2.0) + pow(fabs(z-z2),2.0)));
										C_beta_distances.at(residue_position - offset - 1).at(i) = C_beta_distances.at(i).at(residue_position - offset - 1);
										for (int n=0; n<=3; n++) {
											x2=atomic_coordinates.at(i).at(n).at(0);
											y2=atomic_coordinates.at(i).at(n).at(1);
											z2=atomic_coordinates.at(i).at(n).at(2);
											minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
										}
										if (number_residue_atoms.at(i)>=1) {
											for (int n=4; n<=3 + number_residue_atoms.at(i); n++) {
												x2=atomic_coordinates.at(i).at(n).at(0);
												y2=atomic_coordinates.at(i).at(n).at(1);
												z2=atomic_coordinates.at(i).at(n).at(2);
												residue_distances.at(i).at(residue_position - offset - 1) = min(residue_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
												minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
											}
										}
										residue_distances.at(residue_position - offset - 1).at(i) = residue_distances.at(i).at(residue_position - offset - 1);
										minimum_distances.at(residue_position - offset - 1).at(i) = minimum_distances.at(i).at(residue_position - offset - 1);
									}
								}
							}
						} else if (atom != " OXT" & atom.at(1) != 'H') {
							number_residue_atoms.at(residue_position - offset - 1) = number_residue_atoms.at(residue_position - offset - 1) + 1;
							atomic_coordinates.at(residue_position - offset - 1).at(number_residue_atoms.at(residue_position - offset - 1) + 3).at(0)=x;
							atomic_coordinates.at(residue_position - offset - 1).at(number_residue_atoms.at(residue_position - offset - 1) + 3).at(1)=y;
							atomic_coordinates.at(residue_position - offset - 1).at(number_residue_atoms.at(residue_position - offset - 1) + 3).at(2)=z;
							if (number_residue_atoms.at(residue_position - offset - 1) == 1) {
								for (int i=0; i<residue_position - offset - 1; i++) {
									if (number_residue_atoms.at(i)!=-1) {
										for (int n=0; n<=3; n++) {
											x2=atomic_coordinates.at(i).at(n).at(0);
											y2=atomic_coordinates.at(i).at(n).at(1);
											z2=atomic_coordinates.at(i).at(n).at(2);
											minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
										}
										if (number_residue_atoms.at(i)>=1) {
											x2=atomic_coordinates.at(i).at(4).at(0);
											y2=atomic_coordinates.at(i).at(4).at(1);
											z2=atomic_coordinates.at(i).at(4).at(2);
											residue_distances.at(i).at(residue_position - offset - 1) = float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0)));
											minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
											
											for (int n=5; n<=3 + number_residue_atoms.at(i); n++) {
												x2=atomic_coordinates.at(i).at(n).at(0);
												y2=atomic_coordinates.at(i).at(n).at(1);
												z2=atomic_coordinates.at(i).at(n).at(2);
												residue_distances.at(i).at(residue_position - offset - 1) = min(residue_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
												minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
											}
										}
										residue_distances.at(residue_position - offset - 1).at(i) = residue_distances.at(i).at(residue_position - offset - 1);
										minimum_distances.at(residue_position - offset - 1).at(i) = minimum_distances.at(i).at(residue_position - offset - 1);
									}
								}
							} else {
								for (int i=0; i<residue_position - offset - 1; i++) {
									if (number_residue_atoms.at(i)!=-1) {
										for (int n=0; n<=3; n++) {
											x2=atomic_coordinates.at(i).at(n).at(0);
											y2=atomic_coordinates.at(i).at(n).at(1);
											z2=atomic_coordinates.at(i).at(n).at(2);
											minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
										}
										if (number_residue_atoms.at(i)>=1) {
											for (int n=4; n<=3 + number_residue_atoms.at(i); n++) {
												x2=atomic_coordinates.at(i).at(n).at(0);
												y2=atomic_coordinates.at(i).at(n).at(1);
												z2=atomic_coordinates.at(i).at(n).at(2);
												residue_distances.at(i).at(residue_position - offset - 1) = min(residue_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
												minimum_distances.at(i).at(residue_position - offset - 1) = min(minimum_distances.at(i).at(residue_position - offset - 1),float(sqrt(pow(fabs(x - x2),2.0) + pow(fabs(y - y2),2.0) + pow(fabs(z - z2),2.0))));
											}
										}
										residue_distances.at(residue_position - offset - 1).at(i) = residue_distances.at(i).at(residue_position - offset - 1);
										minimum_distances.at(residue_position - offset - 1).at(i) = minimum_distances.at(i).at(residue_position - offset - 1);
									}
								}
							}
						}
					}
				}
			}
			getline (PDB_file, line);
			record_type.assign(line,0,6);
		}
		if (has_variants == 1) {
			cout << "WARNING! This PDB file has variants!";
		}
		return;
	}
	
	void print_distances_to_directory (string& directory_name) {
		ofstream minimum_distance_file, Ca_distance_file, Cb_distance_file, residue_distance_file;
		string minimum_distance_file_name, Ca_distance_file_name, Cb_distance_file_name, residue_distance_file_name;
		
		minimum_distance_file_name=directory_name+"/min_distance_recalculated.txt";
		minimum_distance_file.open(minimum_distance_file_name.c_str());
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				minimum_distance_file << minimum_distances.at(i).at(j) << " ";
			}
			minimum_distance_file << "\n";
		}
		minimum_distance_file.close();
		
		Ca_distance_file_name=directory_name+"/Ca_distance_recalculated.txt";
		Ca_distance_file.open(Ca_distance_file_name.c_str());
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				Ca_distance_file << C_alpha_distances.at(i).at(j) << " ";
			}
			Ca_distance_file << "\n";
		}
		Ca_distance_file.close();
		
		Cb_distance_file_name=directory_name+"/Cb_distance_recalculated.txt";
		Cb_distance_file.open(Cb_distance_file_name.c_str());
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				Cb_distance_file << C_beta_distances.at(i).at(j) << " ";
			}
			Cb_distance_file << "\n";
		}
		Cb_distance_file.close();
		
		residue_distance_file_name=directory_name+"/residue_distance_recalculated.txt";
		residue_distance_file.open(residue_distance_file_name.c_str());
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				residue_distance_file << residue_distances.at(i).at(j) << " ";
			}
			residue_distance_file << "\n";
		}
		residue_distance_file.close();		
		
		return;
	}
	
	void print_secondary_structure_to_directory (string& directory_name) {
		ofstream helix_file, sheet_file, turn_file;
		string helix_file_name, sheet_file_name, turn_file_name;
		
		helix_file_name=directory_name+"/helices.txt";
		helix_file.open(helix_file_name.c_str());
		for (int i=0; i<10; i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				helix_file << helices.at(i).at(j) << " ";
			}
			helix_file << "\n";
		}
		helix_file.close();
		
		sheet_file_name=directory_name+"/sheets.txt";
		sheet_file.open(sheet_file_name.c_str());
		for (int i=0; i<reference_to_MSA.size(); i++) {
			sheet_file << sheets.at(i) << " ";
		}
		sheet_file << "\n";
		sheet_file.close();
		
		turn_file_name=directory_name+"/turns.txt";
		turn_file.open(turn_file_name.c_str());
		for (int i=0; i<reference_to_MSA.size(); i++) {
			turn_file << turns.at(i) << " ";
		}
		turn_file << "\n";
		turn_file.close();
		
		return;
	}
	
	void display_distances () {
		cout << "\n\n\n\n\n\n\nMinimum Distance\n";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << minimum_distances.at(i).at(j) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n\n\n\n\n\nC alpha distance";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << C_alpha_distances.at(i).at(j) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n\n\n\n\n\nC beta distance";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << C_beta_distances.at(i).at(j) << " ";
			}
			cout << "\n";
		}
		
		cout << "\n";
		
		return;
	}
	
	void my_evolve (string ancestral_sequence, const int target_number_of_sequences, const float duplication_rate, int number_of_interactions, const float interaction_strength) {
		int number_of_new_sequences = 0;
		vector<string> new_generation;
		int size_of_last_generation = 0;
		int remainder = 0;
		float random_number;
		vector<vector<int> > interaction_boundaries (number_of_interactions, vector <int> (2,0));
		int original_aa, interacting_aa, replacement_aa;
		int cumulative_interaction;
		
		clear ();
		alignment_length = ancestral_sequence.length();
		
		if (number_of_interactions > alignment_length/2) {
			number_of_interactions = alignment_length/2;
		}
		for (int i=0; i < number_of_interactions; i++) {
			interaction_boundaries.at(i).at(0) = rand() % 19;
			interaction_boundaries.at(i).at(1) = rand() % 19;
		}
		
		sequences.push_back(ancestral_sequence);
		number_of_sequences = 1;
		
		while (number_of_sequences < target_number_of_sequences) {
			size_of_last_generation = number_of_sequences;
			new_generation.clear();
			number_of_new_sequences=0;
			for (int n=0;n<size_of_last_generation;n++){
				if (float(rand())/float(RAND_MAX) < duplication_rate) {
					new_generation.push_back(sequences.at(n));
					number_of_new_sequences++;
				}				
				for (int i=0; i<alignment_length; i++){
					random_number = float(rand())/(float(RAND_MAX));
					original_aa = aa_index(sequences.at(n).at(i)) - 1;
					replacement_aa = 1;
					if (i > (2 * number_of_interactions) - 1) {
						while ((replacement_aa != 20) && (random_number < cumulative_PAM1.at(original_aa).at(replacement_aa))) {
							replacement_aa = replacement_aa + 1;
						}
						sequences.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));						
					} else if (i < number_of_interactions) {
						interacting_aa = aa_index(sequences.at(n).at(i+number_of_interactions)) - 1;
						if ((interacting_aa <= interaction_boundaries.at(i).at(1)) && (original_aa > interaction_boundaries.at(i).at(0))) {
							cumulative_interaction = interaction_boundaries.at(i).at(0) + 1 - int (original_aa <= interaction_boundaries.at(i).at(0));
							while ((replacement_aa != 20) && (random_number < (cumulative_PAM1.at(original_aa).at(replacement_aa) * (1.0 - interaction_strength) + (float (cumulative_interaction) * interaction_strength / (interaction_boundaries.at(i).at(0) + 1))))) {
								replacement_aa = replacement_aa + 1;
								if (cumulative_interaction !=0) {
									cumulative_interaction = cumulative_interaction - 1;
								}
							}
							sequences.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));
						} else if ((interacting_aa > interaction_boundaries.at(i).at(1)) && (original_aa <= interaction_boundaries.at(i).at(0))) {
							cumulative_interaction = 19 - interaction_boundaries.at(i).at(0) - int (original_aa > interaction_boundaries.at(i).at(0));
							while ((replacement_aa != 20) && (random_number < (cumulative_PAM1.at(original_aa).at(replacement_aa) * (1 - interaction_strength) + (float (cumulative_interaction) * interaction_strength / (interaction_boundaries.at(i).at(0) + 1))))) {
								replacement_aa = replacement_aa + 1;
								if (cumulative_interaction > 20 - replacement_aa) {
									cumulative_interaction = cumulative_interaction - 1;
								}
							}
							sequences.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));													
						} else {
							while ((replacement_aa != 20) && (random_number < cumulative_PAM1.at(original_aa).at(replacement_aa))) {
								replacement_aa = replacement_aa + 1;
							}
							sequences.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));													
						}
					} else { //could change to indexing on bool in interaction_boundaries.at(i).at(0/1)
						interacting_aa = aa_index(sequences.at(n).at(i-number_of_interactions)) - 1;
						if ((interacting_aa <= interaction_boundaries.at(i-number_of_interactions).at(0)) && (original_aa > interaction_boundaries.at(i-number_of_interactions).at(1))) {
							cumulative_interaction = interaction_boundaries.at(i-number_of_interactions).at(1) + 1 - int (original_aa <= interaction_boundaries.at(i-number_of_interactions).at(1));
							while ((replacement_aa != 20) && (random_number < (cumulative_PAM1.at(original_aa).at(replacement_aa) * (1 - interaction_strength) + (float (cumulative_interaction) * interaction_strength / (interaction_boundaries.at(i-number_of_interactions).at(1) + 1))))) {
								replacement_aa = replacement_aa + 1;
								if (cumulative_interaction !=0) {
									cumulative_interaction = cumulative_interaction - 1;
								}
							}
							sequences.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));						
						} else if ((interacting_aa > interaction_boundaries.at(i-number_of_interactions).at(0)) && (original_aa <= interaction_boundaries.at(i-number_of_interactions).at(1))) {
							cumulative_interaction = 19 - interaction_boundaries.at(i-number_of_interactions).at(1) - int (original_aa > interaction_boundaries.at(i-number_of_interactions).at(1));
							while ((replacement_aa != 20) && (random_number < (cumulative_PAM1.at(original_aa).at(replacement_aa) * (1 - interaction_strength) + (float (cumulative_interaction) * interaction_strength / (interaction_boundaries.at(i-number_of_interactions).at(1) + 1))))) {
								replacement_aa = replacement_aa + 1;
								if (cumulative_interaction > 20 - replacement_aa) {
									cumulative_interaction = cumulative_interaction - 1;
								}
							}
							sequences.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));													
						} else {
							while ((replacement_aa != 20) && (random_number < cumulative_PAM1.at(original_aa).at(replacement_aa))) {
								replacement_aa = replacement_aa + 1;
							}
							sequences.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));																				
						}
						
					}
				}
			}
			/*if (float(rand())/float(RAND_MAX) < mutation_rate) {
				new_generation.at(n).at(i)=Amino_Acids.at(rand() % 20);
			if (i <= number_of_interactions - 1) {
				new_generation.at(n).at(i+number_of_interactions)=Amino_Acids.at(rand() % 20);
			}
			else if (i <= 2*number_of_interactions -1) {
				new_generation.at(n).at(i-number_of_interactions)=Amino_Acids.at(rand() % 20);
			}
			}*/
			
			remainder = target_number_of_sequences - number_of_sequences;
			if (number_of_new_sequences <= remainder) {
				for (int n=0;n<number_of_new_sequences;n++) {
					for (int i=0; i<alignment_length; i++){
						random_number = float(rand())/(float(RAND_MAX));
						original_aa = aa_index(new_generation.at(n).at(i)) - 1;
						replacement_aa = 1;
						if (i > (2 * number_of_interactions) - 1) {
							while ((replacement_aa != 20) && (random_number < cumulative_PAM1.at(original_aa).at(replacement_aa))) {
								replacement_aa = replacement_aa + 1;
							}
							new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));						
						} else if (i < number_of_interactions) {
							interacting_aa = aa_index(new_generation.at(n).at(i+number_of_interactions)) - 1;
							if ((interacting_aa <= interaction_boundaries.at(i).at(1)) && (original_aa > interaction_boundaries.at(i).at(0))) {
								cumulative_interaction = interaction_boundaries.at(i).at(0) + 1 - int (original_aa <= interaction_boundaries.at(i).at(0));
								while ((replacement_aa != 20) && (random_number < (cumulative_PAM1.at(original_aa).at(replacement_aa) * (1.0 - interaction_strength) + (float (cumulative_interaction) * interaction_strength / (interaction_boundaries.at(i).at(0) + 1))))) {
									replacement_aa = replacement_aa + 1;
									if (cumulative_interaction !=0) {
										cumulative_interaction = cumulative_interaction - 1;
									}
								}
								new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));
							} else if ((interacting_aa > interaction_boundaries.at(i).at(1)) && (original_aa <= interaction_boundaries.at(i).at(0))) {
								cumulative_interaction = 19 - interaction_boundaries.at(i).at(0) - int (original_aa > interaction_boundaries.at(i).at(0));
								while ((replacement_aa != 20) && (random_number < (cumulative_PAM1.at(original_aa).at(replacement_aa) * (1 - interaction_strength) + (float (cumulative_interaction) * interaction_strength / (interaction_boundaries.at(i).at(0) + 1))))) {
									replacement_aa = replacement_aa + 1;
									if (cumulative_interaction > 20 - replacement_aa) {
										cumulative_interaction = cumulative_interaction - 1;
									}
								}
								new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));													
							} else {
								while ((replacement_aa != 20) && (random_number < cumulative_PAM1.at(original_aa).at(replacement_aa))) {
									replacement_aa = replacement_aa + 1;
								}
								new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));													
							}
						} else { //could change to indexing on bool in interaction_boundaries.at(i).at(0/1)
							interacting_aa = aa_index(new_generation.at(n).at(i-number_of_interactions)) - 1;
							if ((interacting_aa <= interaction_boundaries.at(i-number_of_interactions).at(0)) && (original_aa > interaction_boundaries.at(i-number_of_interactions).at(1))) {
								cumulative_interaction = interaction_boundaries.at(i-number_of_interactions).at(1) + 1 - int (original_aa <= interaction_boundaries.at(i-number_of_interactions).at(1));
								while ((replacement_aa != 20) && (random_number < (cumulative_PAM1.at(original_aa).at(replacement_aa) * (1 - interaction_strength) + (float (cumulative_interaction) * interaction_strength / (interaction_boundaries.at(i-number_of_interactions).at(1) + 1))))) {
									replacement_aa = replacement_aa + 1;
									if (cumulative_interaction !=0) {
										cumulative_interaction = cumulative_interaction - 1;
									}
								}
								new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));						
							} else if ((interacting_aa > interaction_boundaries.at(i-number_of_interactions).at(0)) && (original_aa <= interaction_boundaries.at(i-number_of_interactions).at(1))) {
								cumulative_interaction = 19 - interaction_boundaries.at(i-number_of_interactions).at(1) - int (original_aa > interaction_boundaries.at(i-number_of_interactions).at(1));
								while ((replacement_aa != 20) && (random_number < (cumulative_PAM1.at(original_aa).at(replacement_aa) * (1 - interaction_strength) + (float (cumulative_interaction) * interaction_strength / (interaction_boundaries.at(i-number_of_interactions).at(1) + 1))))) {
									replacement_aa = replacement_aa + 1;
									if (cumulative_interaction > 20 - replacement_aa) {
										cumulative_interaction = cumulative_interaction - 1;
									}
								}
								new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));													
							} else {
								while ((replacement_aa != 20) && (random_number < cumulative_PAM1.at(original_aa).at(replacement_aa))) {
									replacement_aa = replacement_aa + 1;
								}
								new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));																				
							}
							
						}
					}					
					sequences.push_back(new_generation.at(n));
				}
			} else {
				for (int n=0;n<number_of_new_sequences;n++) {
					if (rand() % (number_of_new_sequences - n) < remainder) {
						for (int i=0; i<alignment_length; i++) {
							random_number = float(rand())/(float(RAND_MAX));
							original_aa = aa_index(new_generation.at(n).at(i)) - 1;
							replacement_aa = 1;
							if (i > (2 * number_of_interactions) - 1) {
								while ((replacement_aa != 20) && (random_number < cumulative_PAM1.at(original_aa).at(replacement_aa))) {
									replacement_aa = replacement_aa + 1;
								}
								new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));						
							} else if (i < number_of_interactions) {
								interacting_aa = aa_index(new_generation.at(n).at(i+number_of_interactions)) - 1;
								if ((interacting_aa <= interaction_boundaries.at(i).at(1)) && (original_aa > interaction_boundaries.at(i).at(0))) {
									cumulative_interaction = interaction_boundaries.at(i).at(0) + 1 - int (original_aa <= interaction_boundaries.at(i).at(0));
									while ((replacement_aa != 20) && (random_number < (cumulative_PAM1.at(original_aa).at(replacement_aa) * (1.0 - interaction_strength) + (float (cumulative_interaction) * interaction_strength / (interaction_boundaries.at(i).at(0) + 1))))) {
										replacement_aa = replacement_aa + 1;
										if (cumulative_interaction !=0) {
											cumulative_interaction = cumulative_interaction - 1;
										}
									}
									new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));
								} else if ((interacting_aa > interaction_boundaries.at(i).at(1)) && (original_aa <= interaction_boundaries.at(i).at(0))) {
									cumulative_interaction = 19 - interaction_boundaries.at(i).at(0) - int (original_aa > interaction_boundaries.at(i).at(0));
									while ((replacement_aa != 20) && (random_number < (cumulative_PAM1.at(original_aa).at(replacement_aa) * (1 - interaction_strength) + (float (cumulative_interaction) * interaction_strength / (interaction_boundaries.at(i).at(0) + 1))))) {
										replacement_aa = replacement_aa + 1;
										if (cumulative_interaction > 20 - replacement_aa) {
											cumulative_interaction = cumulative_interaction - 1;
										}
									}
									new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));													
								} else {
									while ((replacement_aa != 20) && (random_number < cumulative_PAM1.at(original_aa).at(replacement_aa))) {
										replacement_aa = replacement_aa + 1;
									}
									new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));													
								}
							} else { //could change to indexing on bool in interaction_boundaries.at(i).at(0/1)
								interacting_aa = aa_index(new_generation.at(n).at(i-number_of_interactions)) - 1;
								if ((interacting_aa <= interaction_boundaries.at(i-number_of_interactions).at(0)) && (original_aa > interaction_boundaries.at(i-number_of_interactions).at(1))) {
									cumulative_interaction = interaction_boundaries.at(i-number_of_interactions).at(1) + 1 - int (original_aa <= interaction_boundaries.at(i-number_of_interactions).at(1));
									while ((replacement_aa != 20) && (random_number < (cumulative_PAM1.at(original_aa).at(replacement_aa) * (1 - interaction_strength) + (float (cumulative_interaction) * interaction_strength / (interaction_boundaries.at(i-number_of_interactions).at(1) + 1))))) {
										replacement_aa = replacement_aa + 1;
										if (cumulative_interaction !=0) {
											cumulative_interaction = cumulative_interaction - 1;
										}
									}
									new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));						
								} else if ((interacting_aa > interaction_boundaries.at(i-number_of_interactions).at(0)) && (original_aa <= interaction_boundaries.at(i-number_of_interactions).at(1))) {
									cumulative_interaction = 19 - interaction_boundaries.at(i-number_of_interactions).at(1) - int (original_aa > interaction_boundaries.at(i-number_of_interactions).at(1));
									while ((replacement_aa != 20) && (random_number < (cumulative_PAM1.at(original_aa).at(replacement_aa) * (1 - interaction_strength) + (float (cumulative_interaction) * interaction_strength / (interaction_boundaries.at(i-number_of_interactions).at(1) + 1))))) {
										replacement_aa = replacement_aa + 1;
										if (cumulative_interaction > 20 - replacement_aa) {
											cumulative_interaction = cumulative_interaction - 1;
										}
									}
									new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));													
								} else {
									while ((replacement_aa != 20) && (random_number < cumulative_PAM1.at(original_aa).at(replacement_aa))) {
										replacement_aa = replacement_aa + 1;
									}
									new_generation.at(n).at(i)=Amino_Acids.at(permutation_PAM1.at(original_aa).at(replacement_aa-1));																				
								}
								
							}
						}
						sequences.push_back(new_generation.at(n));
						remainder--;
					}
				}
			}
			number_of_sequences = sequences.size();
		}
		return;
	} 
	
	void my_evolve (const int target_alignment_length, const int target_number_of_sequences, const float duplication_rate, const int number_of_interactions, const float interaction_strength) {
		string ancestral_sequence;
		vector<float> GAAF (21,0.05);
		
		ancestral_sequence.resize(target_alignment_length);
		for (int i=0;i<target_alignment_length;i++) {
			//ancestral_sequence.at(i) = Amino_Acids.at(rand() % 20);
			ancestral_sequence.at(i) = Amino_Acids.at(0);
		}
		
		my_evolve (ancestral_sequence, target_number_of_sequences,duplication_rate,number_of_interactions, interaction_strength);
		return;
	} 	
	
	
	void evolve (string ancestral_sequence, const int target_number_of_sequences, float duplication_rate, float mutation_rate, int number_of_interactions, vector <float> GAAF, float fraction_unconserved, bool time_stamp, bool Dunn) {
		int number_of_new_sequences = 0;
		vector<string> new_generation;
		int size_of_last_generation = 0;
		int remainder = 0;
		
		clear();
		alignment_length = ancestral_sequence.length();
		if (number_of_interactions > alignment_length/2) {
			number_of_interactions = alignment_length/2;
		}
		sequences.push_back(ancestral_sequence);
		number_of_sequences = 1;
		size_of_last_generation = 1;
		while (number_of_sequences < target_number_of_sequences) {
			size_of_last_generation = number_of_sequences;
			new_generation.clear();
			number_of_new_sequences=0;
			for (int n=0;n<size_of_last_generation;n++){
				if (float(rand())/float(RAND_MAX) < duplication_rate) {
					new_generation.push_back(sequences.at(n));
					number_of_new_sequences++;
				}				
				for (int i=0; i<alignment_length; i++){
					if (float(rand())/(float(RAND_MAX)) < mutation_rate) {
						if (i > (2 * number_of_interactions) - 1 ) {
							sequences.at(n).at(i)=Amino_Acids.at(rand() % 20);
						} else if ((float(rand())/(float(RAND_MAX)) < 1.0)) { //0.5
							if (i < number_of_interactions) {
								sequences.at(n).at(i)=Amino_Acids.at(rand() % 20);
								sequences.at(n).at(i+number_of_interactions)=Amino_Acids.at(rand() % 20);
							} else if (i < (2 * number_of_interactions)) {
								sequences.at(n).at(i)=Amino_Acids.at(rand() % 20);
								sequences.at(n).at(i-number_of_interactions)=Amino_Acids.at(rand() % 20);
							}
						}
					}
				}
			}
			/*if (float(rand())/float(RAND_MAX) < mutation_rate) {
				new_generation.at(n).at(i)=Amino_Acids.at(rand() % 20);
			if (i <= number_of_interactions - 1) {
				new_generation.at(n).at(i+number_of_interactions)=Amino_Acids.at(rand() % 20);
			}
			else if (i <= 2*number_of_interactions -1) {
				new_generation.at(n).at(i-number_of_interactions)=Amino_Acids.at(rand() % 20);
			}
			}*/
			
			remainder = target_number_of_sequences - number_of_sequences;
			if (number_of_new_sequences <= remainder) {
				for (int n=0;n<number_of_new_sequences;n++) {
					for (int i=0; i<alignment_length; i++){
						if (float(rand())/float(RAND_MAX) < mutation_rate) {
							if (i > (2 * number_of_interactions) - 1 ) {
								new_generation.at(n).at(i)=Amino_Acids.at(rand() % 20);
							} else if ((float(rand())/(float(RAND_MAX)) < 1.0)) {//0.5
								if (i < number_of_interactions) {
									new_generation.at(n).at(i)=Amino_Acids.at(rand() % 20);
									new_generation.at(n).at(i+number_of_interactions)=Amino_Acids.at(rand() % 20);
								} else if (i < (2 * number_of_interactions)) {
									new_generation.at(n).at(i)=Amino_Acids.at(rand() % 20);
									new_generation.at(n).at(i-number_of_interactions)=Amino_Acids.at(rand() % 20);
								}
							}
						}
					}					
					sequences.push_back(new_generation.at(n));
				}
			} else {
				for (int n=0;n<number_of_new_sequences;n++) {
					if (rand() % (number_of_new_sequences - n) < remainder) {
						for (int i=0; i<alignment_length; i++) {
							if (float(rand())/(RAND_MAX)<mutation_rate) {
								if (i > (2 * number_of_interactions) - 1 ) {
									new_generation.at(n).at(i)=Amino_Acids.at(rand() % 20);
								} else if ((float(rand())/(float(RAND_MAX)) < 1.0)) {//0.5
									if (i < number_of_interactions) {
										new_generation.at(n).at(i)=Amino_Acids.at(rand() % 20);
										new_generation.at(n).at(i+number_of_interactions)=Amino_Acids.at(rand() % 20);
									} else if (i < (2 * number_of_interactions)) {
										new_generation.at(n).at(i)=Amino_Acids.at(rand() % 20);
										new_generation.at(n).at(i-number_of_interactions)=Amino_Acids.at(rand() % 20);
									}
								}
							}
						}
						sequences.push_back(new_generation.at(n));
						remainder--;
					}
				}
			}
			number_of_sequences = sequences.size();
		}
		return;
	} 
	
	void evolve (int target_alignment_length, int target_number_of_sequences, float duplication_rate, float mutation_rate, int number_of_interactions) {
		string ancestral_sequence;
		vector<float> GAAF (21,0.05);
		
		ancestral_sequence.resize(target_alignment_length);
		for (int i=0;i<target_alignment_length;i++) {
			//ancestral_sequence.at(i) = Amino_Acids.at(rand() % 20);
			ancestral_sequence.at(i) = Amino_Acids.at(0);
		}
		evolve (ancestral_sequence, target_number_of_sequences,duplication_rate,mutation_rate,number_of_interactions, GAAF, 1.0, 0, 1);
		return;
	} 
	
	void display () {
		cout << "\n";
		for (int n=0; n<number_of_sequences; n++) {
			cout <<sequences.at(n) << "\n";
		}
		return;
	}
	
	void display_sequence (int sequence_number) {
		cout << sequences.at(sequence_number) <<"\n";
		return;
	}
	
	void tabulate () {	
		string reformat = "-ACDEFGHIKLMNPQRSTVWY";
		if 	((!aa_counts.empty())||(!aa_frequencies.empty())) {
			cout << "Operation Aborted: This MSA has already been tabulated!\n";
			return;
		}
		clear_tabulation ();
		aa_counts.resize(21,vector<int>(alignment_length,0));
		aa_frequencies.resize(21,vector<float>(alignment_length,0));
		for (int i=0; i<alignment_length; i++) {
			for (int n=0; n<number_of_sequences; n++) {
				aa_counts.at(aa_index(sequences[n][i])).at(i)++;
				sequences[n][i] = reformat.at(aa_index(sequences[n][i]));
			}
			aa_frequencies[0][i]=aa_counts[0][i]/number_of_sequences;
			for (int aa=1; aa<=20; aa++) {
				aa_frequencies.at(aa).at(i)=float(aa_counts.at(aa).at(i))/(number_of_sequences-aa_counts.at(0).at(i));
			}
		}
		return;
	}
	
	
	void display_frequencies (int position) {
		for (int aa=1; aa<=20; aa++) {
			cout << aa_frequencies.at(aa).at(position) << " ";
		}
		cout << "\n";
		return;
	}
	
	void calculate_information () {
		int i_index, j_index;
		vector<vector<int> > joint_counts (21,vector<int> (21,0));
		vector<int> amino_acid_counts_at_i (21,0);
		vector<int> amino_acid_counts_at_j (21,0);
		vector <double> sum_MI_with_i_with_gaps;
		vector <double>	sum_MI_with_i;
		vector <double> mean_MI_with_i_with_gaps;
		//vector <double>	mean_MI_with_i; defined above
		vector <double>	prestd_MI_with_i_with_gaps;
		vector <double>	std_MI_with_i_with_gaps;
		vector <double>	prestd_MI_with_i;
		//vector <double>	std_MI_with_i;
		vector<int> number_of_subsets_with_i_above_cutoff (0);
		double sum_meanimeanj_with_gaps (0);
		double sum_meanimeanj (0);
		double sum_MI_with_gaps (0);
		double sum_MI (0);
		int total_number_of_subsets_above_cutoff (0);		
		double mean_meanimeanj_with_gaps (0);
		double mean_meanimeanj (0);
		double mean_MI_with_gaps (0);
		double mean_MI (0);
		double Sx_with_gaps (0);
		double Sy_with_gaps (0);
		double Sxy_with_gaps (0);
		double Sx (0);
		double Sy (0);
		double Sxy (0);
		double slope_with_gaps (0);
		double intercept_with_gaps (0);
		double slope (0);
		double intercept (0);
		
		
		if 	(!MI_with_gaps.empty()) {
			cout << "Operation Aborted: This MSA has already been mutually informated!\n";
			return;
		}
		clear_information();
		MI_with_gaps.resize(alignment_length,vector<double>(alignment_length,0.0/0.0));
		MI.resize(alignment_length,vector<double>(alignment_length,0.0/0.0));
		subset_Hx.resize(alignment_length,vector<double>(alignment_length,0));
		Hxy.resize(alignment_length,vector<double>(alignment_length,0));
		MI_residual_with_gaps.resize(alignment_length,vector<double>(alignment_length,0.0/0.0));
		MI_residual.resize(alignment_length,vector<double>(alignment_length,0.0/0.0));
		Dunn_with_gaps.resize(alignment_length,vector<double>(alignment_length,0.0/0.0));
		Dunn.resize(alignment_length,vector<double>(alignment_length,0.0/0.0));
		Pearson_with_gaps.resize(alignment_length,vector<double>(alignment_length,0.0/0.0));
		Pearson.resize(alignment_length,vector<double>(alignment_length,0.0/0.0));	
		Hx_with_gaps.resize(alignment_length,0.0);
		Hxy_with_gaps.resize(alignment_length,vector<double>(alignment_length,0.0));
		subset_size.resize(alignment_length,vector<int>(alignment_length,0));
		// calculate single distribution entropies with gaps as the 21st amino acid
		for (int i=0;i<alignment_length;i++) {
			for (int aa=0;aa<=20;aa++) {
				if (aa_counts.at(aa).at(i)!=0) {
					Hx_with_gaps.at(i)=Hx_with_gaps.at(i)-(double(aa_counts.at(aa).at(i))/double(number_of_sequences))*  log2(double(aa_counts.at(aa).at(i))/double(number_of_sequences));
				}
			}
		}
		sum_MI_with_i_with_gaps.resize(alignment_length,0.0);
		sum_MI_with_i.resize(alignment_length,0.0);
		mean_MI_with_i_with_gaps.resize(alignment_length,0.0);
		mean_MI_with_i.resize(alignment_length,0.0);
		prestd_MI_with_i_with_gaps.resize(alignment_length,0.0);
		prestd_MI_with_i.resize(alignment_length,0.0);
		std_MI_with_i_with_gaps.resize(alignment_length,0.0);
		std_MI_with_i.resize(alignment_length,0.0);
		
		number_of_subsets_with_i_above_cutoff.resize(alignment_length,0);
		for (int i=0;i<alignment_length-1;i++) {
			for (int j=i+1;j<alignment_length;j++) {
				joint_counts=vector<vector<int> > (21,vector<int> (21,0));
				amino_acid_counts_at_i=vector<int>(21,0);
				amino_acid_counts_at_j=vector<int>(21,0);
				// count the joint occurances of amino acids at i and j
				for (int n=0;n<number_of_sequences;n++) {
					i_index=aa_index(sequences.at(n)[i]);
					j_index=aa_index(sequences.at(n)[j]);
					joint_counts.at(i_index).at(j_index)++;
					if (i_index!=0 && j_index!=0) {
						amino_acid_counts_at_i.at(i_index)++;
						amino_acid_counts_at_j.at(j_index)++;
						subset_size.at(i).at(j)++;
					}
				}
				// calculate single distribution entropies at x igoring gaps with y
				for (int aa=0;aa<=20;aa++) {
					if (amino_acid_counts_at_i.at(aa)!=0) {
						subset_Hx.at(i).at(j)=subset_Hx.at(i).at(j)-(double(amino_acid_counts_at_i.at(aa))/double(subset_size.at(i).at(j))) * log2(double(amino_acid_counts_at_i.at(aa))/double(subset_size.at(i).at(j)));
					}
					if (amino_acid_counts_at_j.at(aa)!=0) {
						subset_Hx.at(j).at(i)=subset_Hx.at(j).at(i)-(double(amino_acid_counts_at_j.at(aa))/double(subset_size.at(i).at(j))) * log2(double(amino_acid_counts_at_j.at(aa))/double(subset_size.at(i).at(j)));
					}
				}
				//calculate joint distribution entropies
				for (int i_aa=0;i_aa<=20;i_aa++) {
					for (int j_aa=0;j_aa<=20;j_aa++) {
						if (joint_counts.at(i_aa).at(j_aa)!=0) {
							Hxy_with_gaps.at(i).at(j)=Hxy_with_gaps.at(i).at(j)-(double(joint_counts.at(i_aa).at(j_aa))/double(number_of_sequences))*  log2(double(joint_counts.at(i_aa).at(j_aa))/double(number_of_sequences));
							if (i_aa!=0 && j_aa!=0) {
								Hxy.at(i).at(j)=Hxy.at(i).at(j)-(double(joint_counts.at(i_aa).at(j_aa))/double(subset_size.at(i).at(j))) * log2(double(joint_counts.at(i_aa).at(j_aa))/double(subset_size.at(i).at(j)));
							}
						}
					}
				}
				subset_size.at(j).at(i)=subset_size.at(i).at(j);
				Hxy_with_gaps.at(j).at(i) = Hxy_with_gaps.at(i).at(j);
				Hxy.at(j).at(i) = Hxy.at(i).at(j);
				//MI = Hx + Hy - Hxy
				MI_with_gaps.at(i).at(j)=Hx_with_gaps.at(i)+Hx_with_gaps.at(j)-Hxy_with_gaps.at(i).at(j);
				MI_with_gaps.at(j).at(i)=MI_with_gaps.at(i).at(j);
				Dunn_with_gaps.at(i).at(j)=(Hx_with_gaps.at(i)+Hx_with_gaps.at(j)-Hxy_with_gaps.at(i).at(j))/Hxy_with_gaps.at(i).at(j);
				Dunn_with_gaps.at(j).at(i)=Dunn_with_gaps.at(i).at(j);
				if (subset_size.at(i).at(j) >= int(0.8 * number_of_sequences)){
					MI.at(i).at(j)=subset_Hx.at(i).at(j)+subset_Hx.at(j).at(i)-Hxy.at(i).at(j);
					MI.at(j).at(i)=MI.at(i).at(j);
					Dunn.at(i).at(j)=(subset_Hx.at(i).at(j)+subset_Hx.at(j).at(i)-Hxy.at(i).at(j))/Hxy.at(i).at(j);
					Dunn.at(j).at(i)=Dunn.at(i).at(j);
				}  else {
					MI.at(i).at(j)=0.0/0.0;
					MI.at(j).at(i)=MI.at(i).at(j);
					Dunn.at(i).at(j)=0.0/0.0;
					Dunn.at(j).at(i)=Dunn.at(i).at(j);
				}
				//sum MI at i across all j
				sum_MI_with_i_with_gaps.at(i)=sum_MI_with_i_with_gaps.at(i)+MI_with_gaps.at(i).at(j);
				sum_MI_with_i_with_gaps.at(j)=sum_MI_with_i_with_gaps.at(j)+MI_with_gaps.at(i).at(j);
				if (subset_size.at(i).at(j) >= int(0.8 * number_of_sequences)) {
					sum_MI_with_i.at(i) = sum_MI_with_i.at(i) + MI.at(i).at(j);
					sum_MI_with_i.at(j) = sum_MI_with_i.at(j) + MI.at(i).at(j);
					number_of_subsets_with_i_above_cutoff.at(i) = number_of_subsets_with_i_above_cutoff.at(i) + 1;
					number_of_subsets_with_i_above_cutoff.at(j) = number_of_subsets_with_i_above_cutoff.at(j) + 1;
				}
			}
			//average MI at i across all j
			mean_MI_with_i_with_gaps.at(i) = sum_MI_with_i_with_gaps.at(i) / (alignment_length - 1);
			mean_MI_with_i.at(i) = sum_MI_with_i.at(i) / number_of_subsets_with_i_above_cutoff.at(i);
		}
		mean_MI_with_i_with_gaps.at(alignment_length-1) = sum_MI_with_i_with_gaps.at(alignment_length-1) / (alignment_length - 1);
		mean_MI_with_i.at(alignment_length-1) = sum_MI_with_i.at(alignment_length-1) / number_of_subsets_with_i_above_cutoff.at(alignment_length-1);
		
		for (int i=0;i<alignment_length-1;i++) {
			for (int j=i+1;j<alignment_length;j++) {
				prestd_MI_with_i_with_gaps.at(i) = prestd_MI_with_i_with_gaps.at(i) + pow(MI.at(i).at(j) - mean_MI_with_i_with_gaps.at(i),2.0);
				prestd_MI_with_i_with_gaps.at(j) = prestd_MI_with_i_with_gaps.at(j) + pow(MI.at(i).at(j) - mean_MI_with_i_with_gaps.at(j),2.0);
				if (subset_size.at(i).at(j) >= int(0.8 * number_of_sequences)) {
					prestd_MI_with_i.at(i) = prestd_MI_with_i.at(i) + pow(MI.at(i).at(j) - mean_MI_with_i.at(i),2.0);
					prestd_MI_with_i.at(j) = prestd_MI_with_i.at(j) + pow(MI.at(i).at(j) - mean_MI_with_i.at(j),2.0);
				}
			}
			std_MI_with_i_with_gaps.at(i) = pow(prestd_MI_with_i_with_gaps.at(i) / (alignment_length - 2.0),0.5);
			std_MI_with_i.at(i) = pow(prestd_MI_with_i.at(i) / (number_of_subsets_with_i_above_cutoff.at(i) - 1.0),0.5);
		}
		std_MI_with_i_with_gaps.at(alignment_length-1) = pow(prestd_MI_with_i_with_gaps.at(alignment_length-1) / (alignment_length - 2.0),0.5);
		std_MI_with_i.at(alignment_length-1) = pow(prestd_MI_with_i.at(alignment_length-1) / (number_of_subsets_with_i_above_cutoff.at(alignment_length-1) - 1.0),0.5);
		
		//Now we will calculate the regression of MI(i,j) onto avg(MI(i))*avg(MI(j)) and the Pearson measure
		for (int i=0;i<alignment_length-1;i++) {
			for (int j=i+1;j<alignment_length;j++) {
				if (MI_with_gaps.at(i).at(j) >= mean_MI_with_i_with_gaps.at(i) && MI_with_gaps.at(i).at(j) >= mean_MI_with_i_with_gaps.at(j)) {
					Pearson_with_gaps.at(i).at(j)=(MI_with_gaps.at(i).at(j) - mean_MI_with_i_with_gaps.at(i))*(MI_with_gaps.at(i).at(j) - mean_MI_with_i_with_gaps.at(j))/(std_MI_with_i_with_gaps.at(i)*std_MI_with_i_with_gaps.at(j));
					Pearson_with_gaps.at(j).at(i)=Pearson_with_gaps.at(i).at(j);
				} else {
					Pearson_with_gaps.at(i).at(j)=0.0/0.0;
					Pearson_with_gaps.at(j).at(i)=Pearson_with_gaps.at(i).at(j);
				}
				if (MI.at(i).at(j) >= mean_MI_with_i.at(i) && MI.at(i).at(j) >= mean_MI_with_i.at(j)) {
					Pearson.at(i).at(j)=(MI.at(i).at(j) - mean_MI_with_i.at(i))*(MI.at(i).at(j) - mean_MI_with_i.at(j))/(std_MI_with_i.at(i)*std_MI_with_i.at(j));
					Pearson.at(j).at(i)=Pearson.at(i).at(j);
				} else {
					Pearson.at(i).at(j)=0.0/0.0;
					Pearson.at(j).at(i)=Pearson.at(i).at(j);
				}			
				sum_meanimeanj_with_gaps = sum_meanimeanj_with_gaps + (mean_MI_with_i_with_gaps.at(i) * mean_MI_with_i_with_gaps.at(j));
				sum_MI_with_gaps = sum_MI_with_gaps + MI_with_gaps.at(i).at(j);
				if (subset_size.at(i).at(j)  >= int(0.8 * number_of_sequences)) {
					sum_meanimeanj = sum_meanimeanj + (mean_MI_with_i.at(i) * mean_MI_with_i.at(j));
					sum_MI = sum_MI + MI.at(i).at(j);
					total_number_of_subsets_above_cutoff++;
				}
			}
		}
		mean_meanimeanj_with_gaps = sum_meanimeanj_with_gaps / double(alignment_length * (alignment_length - 1) / 2);
		mean_MI_with_gaps = sum_MI_with_gaps / double(alignment_length * (alignment_length - 1) / 2);
		mean_meanimeanj = sum_meanimeanj / total_number_of_subsets_above_cutoff;
		mean_MI = sum_MI / total_number_of_subsets_above_cutoff;
		for (int i=0;i<alignment_length-1;i++) {
			for (int j=i+1;j<alignment_length;j++) {
				Sx_with_gaps = Sx_with_gaps + pow(mean_meanimeanj_with_gaps - (mean_MI_with_i_with_gaps.at(i) * mean_MI_with_i_with_gaps.at(j)),double (2));
				Sy_with_gaps = Sy_with_gaps + pow(mean_MI_with_gaps - MI_with_gaps.at(i).at(j),double (2));
				Sxy_with_gaps = Sxy_with_gaps + ((mean_MI_with_gaps - MI_with_gaps.at(i).at(j)) * (mean_meanimeanj_with_gaps - (mean_MI_with_i_with_gaps.at(i) * mean_MI_with_i_with_gaps.at(j))));
				if (subset_size.at(i).at(j)  >= int(0.8 * number_of_sequences)) {
					Sx = Sx + pow(mean_meanimeanj - (mean_MI_with_i.at(i) * mean_MI_with_i.at(j)),double (2));
					Sy = Sy + pow(mean_MI - MI.at(i).at(j),double (2));
					Sxy = Sxy + ((mean_MI - MI.at(i).at(j)) * (mean_meanimeanj - (mean_MI_with_i.at(i) * mean_MI_with_i.at(j))));
				}
			}
		}
		slope_with_gaps = Sxy_with_gaps / Sx_with_gaps;
		intercept_with_gaps = mean_MI_with_gaps - (slope_with_gaps * mean_meanimeanj_with_gaps);
		slope = Sxy / Sx;
		intercept = mean_MI - (slope * mean_meanimeanj);
		for (int i=0;i<alignment_length-1;i++) {
			for (int j=i+1;j<alignment_length;j++) {
				MI_residual_with_gaps.at(i).at(j) = MI_with_gaps.at(i).at(j) - (slope_with_gaps * mean_MI_with_i_with_gaps.at(i) * mean_MI_with_i_with_gaps.at(j)) - intercept_with_gaps;
				if (subset_size.at(i).at(j)  >= int(0.8 * number_of_sequences) ) {
					MI_residual.at(i).at(j) = MI.at(i).at(j) - (slope * mean_MI_with_i.at(i) * mean_MI_with_i.at(j)) - intercept;
					MI_residual.at(j).at(i) = MI_residual.at(i).at(j);
				} else {
					MI_residual.at(i).at(j) = 0.0/0.0;
					MI_residual.at(j).at(i) = MI_residual.at(i).at(j);
				}
			}
		}		
		return;
	}
	
	void print_information_to_directory (string& directory_name) {
		ofstream MI_file, Residual_file, Dunn_file, Pearson_file, Entropy_file;
		string MI_file_name, Residual_file_name, Dunn_file_name, Pearson_file_name, Entropy_file_name;
		
		MI_file_name=directory_name+"/MI.txt";
		MI_file.open(MI_file_name.c_str());
		for (int i=0; i<alignment_length; i++) {
			for (int j=0; j<alignment_length; j++) {
				MI_file << MI.at(i).at(j) << " ";
			}
			MI_file << "\n";
		}
		MI_file.close();
		
		Residual_file_name=directory_name+"/Residual.txt";
		Residual_file.open(Residual_file_name.c_str());
		for (int i=0; i<alignment_length; i++) {
			for (int j=0; j<alignment_length; j++) {
				Residual_file << MI_residual.at(i).at(j) << " ";
			}
			Residual_file << "\n";
		}
		Residual_file.close();
		
		Dunn_file_name=directory_name+"/Dunn.txt";
		Dunn_file.open(Dunn_file_name.c_str());
		for (int i=0; i<alignment_length; i++) {
			for (int j=0; j<alignment_length; j++) {
				Dunn_file << Dunn.at(i).at(j) << " ";
			}
			Dunn_file << "\n";
		}
		Dunn_file.close();
		
		Pearson_file_name=directory_name+"/Pearson.txt";
		Pearson_file.open(Pearson_file_name.c_str());
		for (int i=0; i<alignment_length; i++) {
			for (int j=0; j<alignment_length; j++) {
				Pearson_file << Pearson.at(i).at(j) << " ";
			}
			Pearson_file << "\n";
		}
		Pearson_file.close();
		
		Entropy_file_name=directory_name+"/Entropy.txt";
		Entropy_file.open(Entropy_file_name.c_str());
		for (int i=0; i<alignment_length; i++) {
			for (int j=0; j<alignment_length; j++) {
				Entropy_file << subset_Hx.at(i).at(j) << " ";
			}
			Entropy_file << "\n";
		}
		Entropy_file.close();
		
	}
	
	
	void calculate_other_measures () {
		int i_index, j_index, n_index, m_index, i_index2, j_index2;
		vector<vector<int> > joint_counts (21,vector<int> (21,0));
		vector<int> amino_acid_counts_at_i (21,0);
		vector<int> amino_acid_counts_at_j (21,0);
		vector<double> mean_NxM;
		vector<double> std_NxM;
		double NEX, sum_NxM, sum_sqr_NxM, pre_McBASC, n_NxM, number_under_5, n_NEX;
		
		if (double(alignment_length)*alignment_length*number_of_sequences*number_of_sequences<=6000000000.0){
			return;
		}
		
		mean_NxM.resize(alignment_length,0.0/0.0);
		std_NxM.resize(alignment_length,0.0/0.0);
		subset_size.resize(alignment_length,vector<int>(alignment_length,0));
		McBASC.resize(alignment_length,vector<double>(alignment_length,0.0/0.0));
		GOF.resize(alignment_length,vector<double>(alignment_length,0.0/0.0));
		OMES.resize(alignment_length,vector<double>(alignment_length,0.0/0.0));
		
		for (int i=0;i<alignment_length;i++) {
			n_NxM=0.0;
			sum_NxM=0.0;
			sum_sqr_NxM=0.0;
			for (int n=0;n<number_of_sequences;n++) {
				for (int m=0;m<number_of_sequences;m++) {
					n_index=aa_index(sequences.at(n)[i]);
					m_index=aa_index(sequences.at(m)[i]);
					if (n_index!=0 && m_index!=0) {
						n_NxM++;
						sum_NxM=sum_NxM+ scoring_matrix.at(n_index).at(m_index);
						sum_sqr_NxM=sum_sqr_NxM+pow(scoring_matrix.at(n_index).at(m_index),double(2));
					}
				}
			}	
			mean_NxM.at(i)=sum_NxM/n_NxM;
			std_NxM.at(i)=pow((sum_sqr_NxM-(sum_NxM*mean_NxM.at(i)))/(n_NxM - 1.0),0.5);
		}

		for (int i=0;i<alignment_length-1;i++) {
			for (int j=i+1;j<alignment_length;j++) {
				joint_counts=vector<vector<int> > (21,vector<int> (21,0));
				amino_acid_counts_at_i=vector<int>(21,0);
				amino_acid_counts_at_j=vector<int>(21,0);
				pre_McBASC=0.0;
				// count the joint occurances of amino acids at i and j
				for (int n=0;n<number_of_sequences;n++) {
					i_index=aa_index(sequences.at(n)[i]);
					j_index=aa_index(sequences.at(n)[j]);
					joint_counts.at(i_index).at(j_index)++;
					if (i_index!=0 && j_index!=0) {
						amino_acid_counts_at_i.at(i_index)++;
						amino_acid_counts_at_j.at(j_index)++;
						subset_size.at(i).at(j)++;
						for (int m=0;m<number_of_sequences;m++) {
							i_index2=aa_index(sequences.at(m)[i]);
							j_index2=aa_index(sequences.at(m)[j]);
							if (i_index2!=0 && j_index2!=0) {
								pre_McBASC=pre_McBASC + ((scoring_matrix.at(i_index).at(i_index2) - mean_NxM.at(i))*(scoring_matrix.at(j_index).at(j_index2) - mean_NxM.at(j)));
							}
						}
					}
				}
				McBASC.at(i).at(j)=pre_McBASC/(double(number_of_sequences)*double(number_of_sequences)*std_NxM.at(i)*std_NxM.at(j));
				
				subset_size.at(j).at(i)=subset_size.at(i).at(j);
				
				number_under_5=0.0;
				GOF.at(i).at(j)=0.0;
				OMES.at(i).at(j)=0.0;
				n_NEX=0;
				for (int i_aa=1;i_aa<=20;i_aa++) {
					for (int j_aa=1;j_aa<=20;j_aa++) {
						if ((amino_acid_counts_at_i.at(i_aa)!=0)&(amino_acid_counts_at_j.at(j_aa)!=0)) {
							n_NEX++;
							NEX=double(amino_acid_counts_at_i.at(i_aa))*double(amino_acid_counts_at_j.at(j_aa))/double(subset_size.at(i).at(j));
							OMES.at(i).at(j)=OMES.at(i).at(j)+ ((double(joint_counts.at(i_aa).at(j_aa)) - NEX)*(double(joint_counts.at(i_aa).at(j_aa)) - NEX)/double(subset_size.at(i).at(j)));
							if (NEX<1.0) {
								GOF.at(i).at(j)=0.0/0.0;
							} else {
								GOF.at(i).at(j)=GOF.at(i).at(j)+ ((double(joint_counts.at(i_aa).at(j_aa)) - NEX)*(double(joint_counts.at(i_aa).at(j_aa)) - NEX)/NEX);
								if (NEX<5.0) {
									number_under_5++;
								}
							}
						}
					}
				}
				if (number_under_5>(0.2 * double(n_NEX))){
					GOF.at(i).at(j)=0.0/0.0;
				}
				McBASC.at(j).at(i)=McBASC.at(i).at(j);
				GOF.at(j).at(i)=GOF.at(i).at(j);
				OMES.at(j).at(i)=OMES.at(i).at(j);
			}
		}
		return;
	}

	void print_other_to_directory (string& directory_name) {
		ofstream GOF_file, OMES_file, McBASC_file;
		string GOF_file_name, OMES_file_name, McBASC_file_name;
		
				if (double(alignment_length)*alignment_length*number_of_sequences*number_of_sequences<=6000000000.0){
					return;
				}
				
		
		GOF_file_name=directory_name+"/GOF.txt";
		GOF_file.open(GOF_file_name.c_str());
		for (int i=0; i<alignment_length; i++) {
			for (int j=0; j<alignment_length; j++) {
				GOF_file << GOF.at(i).at(j) << " ";
			}
			GOF_file << "\n";
		}
		GOF_file.close();
		
				OMES_file_name=directory_name+"/OMES.txt";
				OMES_file.open(OMES_file_name.c_str());
				for (int i=0; i<alignment_length; i++) {
					for (int j=0; j<alignment_length; j++) {
						OMES_file << OMES.at(i).at(j) << " ";
					}
					OMES_file << "\n";
				}
				OMES_file.close();				
		
		McBASC_file_name=directory_name+"/McBSAC.txt";
		McBASC_file.open(McBASC_file_name.c_str());
		for (int i=0; i<alignment_length; i++) {
			for (int j=0; j<alignment_length; j++) {
				McBASC_file << McBASC.at(i).at(j) << " ";
			}
			McBASC_file << "\n";
		}
		McBASC_file.close();
		
	}
	
	
	void print_perspective_information_to_directory (string& directory_name, string& PDB_id) {
		ofstream MI_file, Residual_file, Dunn_file, Pearson_file, Entropy_file;
		string MI_file_name, Residual_file_name, Dunn_file_name, Pearson_file_name, Entropy_file_name;
		
		MI_file_name=directory_name+"/MI_" + PDB_id + ".txt";
		MI_file.open(MI_file_name.c_str());
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				MI_file << MI.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			MI_file << "\n";
		}
		MI_file.close();
		
		Residual_file_name=directory_name+"/Residual_" + PDB_id + ".txt";
		Residual_file.open(Residual_file_name.c_str());
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				Residual_file << MI_residual.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			Residual_file << "\n";
		}
		Residual_file.close();
		
		Dunn_file_name=directory_name+"/Dunn_" + PDB_id + ".txt";
		Dunn_file.open(Dunn_file_name.c_str());
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				Dunn_file << Dunn.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			Dunn_file << "\n";
		}
		Dunn_file.close();
		
		Pearson_file_name=directory_name+"/Pearson_" + PDB_id + ".txt";
		Pearson_file.open(Pearson_file_name.c_str());
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				Pearson_file << Pearson.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			Pearson_file << "\n";
		}
		Pearson_file.close();
		
		Entropy_file_name=directory_name+"/Entropy_" + PDB_id + ".txt";
		Entropy_file.open(Entropy_file_name.c_str());
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				Entropy_file << subset_Hx.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			Entropy_file << "\n";
		}
		Entropy_file.close();
		
	}
	
	
	void display_information () {
		cout << "\n";
		cout << "\nMI\n";
		for (int i=0; i<alignment_length; i++) {
			for (int j=0; j<alignment_length; j++) {
				cout << MI.at(i).at(j) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n\n\n\n\n\n\nResiduals\n";
		for (int i=0; i<alignment_length; i++) {
			for (int j=0; j<alignment_length; j++) {
				cout << MI_residual.at(i).at(j) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n\n\n\n\n\n\nPearson\n";
		for (int i=0; i<alignment_length; i++) {
			for (int j=0; j<alignment_length; j++) {
				cout << Pearson.at(i).at(j) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n\n\n\n\n\n\nMI/Hxy\n";
		for (int i=0; i<alignment_length; i++) {
			for (int j=0; j<alignment_length; j++) {
				cout << Dunn.at(i).at(j) << " ";
			}
			cout << "\n";
		}
		
		cout << "\n\n\n\n\n";
		
		return;
	}
	
	
	void display_perspective_information () {
		cout << "\n";
		cout << "\nMI\n";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << MI.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n\n\n\n\n\n\nResiduals\n";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << MI_residual.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n\n\n\n\n\n\nPearson\n";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << Pearson.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n\n\n\n\n\n\nMI/Hxy\n";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << Dunn.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			cout << "\n";
		}
		
		cout << "\n\n\n\n\n";
		/*
		 cout << "\n";
		 for (int i=0; i<reference_to_MSA.size(); i++) {
			 cout << mean_MI_with_i.at(reference_to_MSA.at(i)) << " ";
		 }
		 cout << "\n\n\n\n\n";	
		 cout << "\n";
		 for (int i=0; i<reference_to_MSA.size(); i++) {
			 cout << std_MI_with_i.at(reference_to_MSA.at(i)) << " ";
		 }*/
		cout << "\n\n\n\n\n";	
		
		cout << "Single Entropy\n";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << subset_Hx.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n\n\n";
		cout << "Joint Entropy\n";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << Hxy.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			cout << "\n";
		}
		
		return;
	}
	
	void display_gapped_perspective_information () {
		cout << "\n";
		cout << "\nGapped MI\n";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << MI_with_gaps.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n\n\n\n\n\n\nGapped Residuals\n";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << MI_residual_with_gaps.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n\n\n\n\n\n\nGapped Pearson\n";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << Pearson_with_gaps.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n\n\n\n\n\n\nGapped MI/Hxy\n";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << Dunn_with_gaps.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			cout << "\n";
		}
		
		cout << "\n\n\n\n\n";
		
		return;
	}
	
	void display_entropies_perspective () {
		cout << "\n";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << subset_Hx.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			cout << "\n";
		}
		return;
	}
	
	
	
	void display_perspective_degapped_sizes () {
		cout << "\n";
		for (int i=0; i<reference_to_MSA.size(); i++) {
			for (int j=0; j<reference_to_MSA.size(); j++) {
				cout << subset_size.at(reference_to_MSA.at(i)).at(reference_to_MSA.at(j)) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n\n\n";
		
		return;
	}
	
	void display_information (int position) {
		for (int i=0; i<position-1; i++) {
			cout << MI_residual.at(i).at(position-1) << " ";
		}
		for (int i=position; i<alignment_length; i++) {
			cout << MI_residual.at(position -1).at(i) << " ";
		}		
		cout << "\n\n";
		
		for (int i=0; i<position-1; i++) {
			cout << MI.at(i).at(position-1) / Hxy.at(i).at(position-1) << " ";
		}
		for (int i=position; i<alignment_length; i++) {
			cout << MI.at(position -1).at(i) / Hxy.at(position -1).at(i) << " ";
		}		
		cout << "\n\n";
		
		for (int i=0; i<position-1; i++) {
			cout << MI.at(i).at(position-1) << " ";
		}
		for (int i=position; i<alignment_length; i++) {
			cout << MI.at(position -1).at(i) << " ";
		}		
		cout << "\n";
		
		return;
	}
	
	void test_normalization (int target_alignment_length, int target_number_of_sequences, float duplication_rate, float mutation_rate, int number_of_interactions, float& residual_sensitivity, float& residual_specificity, float& Dunn_sensitivity, float& Dunn_specificity, float& mi_sensitivity, float& mi_specificity) {
		int residual_true_positives (0);
		int residual_false_positives (0);
		int Dunn_true_positives (0);
		int Dunn_false_positives (0);
		int mi_true_positives (0);
		int mi_false_positives (0);
		double mu_residual (0);
		double sigma_residual (0);
		double mu_Dunn (0);
		double sigma_Dunn (0);
		double mu_mi (0);
		double sigma_mi (0);
		
		residual_true_positives=0;
		residual_false_positives=0;
		Dunn_true_positives=0;
		Dunn_false_positives=0;
		evolve (target_alignment_length, target_number_of_sequences, duplication_rate, number_of_interactions, mutation_rate);
		tabulate();
		calculate_information();
		for (int i=0;i<alignment_length-1;i++) {
			for (int j=i+1;j<alignment_length;j++) {
				mu_residual = mu_residual + MI_residual_with_gaps.at(i).at(j);
				if (Hxy_with_gaps.at(i).at(j) != 0) {
					mu_Dunn = mu_Dunn + (MI_with_gaps.at(i).at(j) / Hxy_with_gaps.at(i).at(j));
				}
				mu_mi = mu_mi + MI_with_gaps.at(i).at(j);
			}
		}
		mu_residual = mu_residual / double(alignment_length * (alignment_length - 1)/ 2);
		mu_Dunn = mu_Dunn / double(alignment_length * (alignment_length - 1) / 2);
		mu_mi = mu_mi / double(alignment_length * (alignment_length - 1) / 2);
		for (int i=0;i<alignment_length-1;i++) {
			for (int j=i+1;j<alignment_length;j++) {
				sigma_residual = sigma_residual + pow(MI_residual_with_gaps.at(i).at(j) - mu_residual,double(2));
				if (Hxy_with_gaps.at(i).at(j) != 0) {
					sigma_Dunn = sigma_Dunn + pow(MI_with_gaps.at(i).at(j) / Hxy_with_gaps.at(i).at(j) - mu_Dunn,double(2));
				} else {
					sigma_Dunn = sigma_Dunn + pow(0.0 - mu_Dunn,double(2));
					
				}
				sigma_mi = sigma_mi + pow(MI_with_gaps.at(i).at(j) - mu_mi,double(2));
			}
		}
		sigma_residual = sqrt(sigma_residual / double(alignment_length * (alignment_length - 1) / 2 - 1));
		sigma_Dunn = sqrt(sigma_Dunn / double(alignment_length * (alignment_length - 1) / 2 - 1));
		sigma_mi = sqrt(sigma_mi / double(alignment_length * (alignment_length - 1) / 2 - 1));
		for (int i=0;i<alignment_length-1;i++) {
			for (int j=i+1;j<alignment_length;j++) {
				if (MI_residual_with_gaps.at(i).at(j) > mu_residual + 4 * sigma_residual) {
					if ((j-i == number_of_interactions)&&(j < 2 * number_of_interactions)) {
						residual_true_positives++;
					} else {
						residual_false_positives++;
					}
				}
				if (Hxy_with_gaps.at(i).at(j) != 0) {
					if ((MI_with_gaps.at(i).at(j) / Hxy_with_gaps.at(i).at(j)) > mu_Dunn + 4 * sigma_Dunn) {
						if ((j-i == number_of_interactions)&&(j < 2 * number_of_interactions)) {
							Dunn_true_positives++;
						} else {
							Dunn_false_positives++;
						}
					}
				} else {
					if (0.0 > mu_Dunn + 4 * sigma_Dunn) {
						if ((j-i == number_of_interactions)&&(j < 2 * number_of_interactions)) {
							Dunn_true_positives++;
						} else {
							Dunn_false_positives++;
						}
					}					
				}
				if (MI_with_gaps.at(i).at(j) > mu_mi + 4 * sigma_mi) {
					if ((j-i == number_of_interactions)&&(j < 2 * number_of_interactions)) {
						mi_true_positives++;
					} else {
						mi_false_positives++;
					}
				}
			}				
		}
		
		residual_sensitivity = float(residual_true_positives) / float(number_of_interactions);
		Dunn_sensitivity = float(Dunn_true_positives) / float(number_of_interactions);
		mi_sensitivity = float(mi_true_positives) / float(number_of_interactions);
		
		if (residual_true_positives + residual_false_positives != 0) {
			residual_specificity = float(residual_true_positives) / float(residual_true_positives + residual_false_positives);
		} else {
			residual_specificity = 0;
		}
		
		if (Dunn_true_positives + Dunn_false_positives != 0) {
			Dunn_specificity = float(Dunn_true_positives) / float(Dunn_true_positives + Dunn_false_positives);
		} else {
			Dunn_specificity = 0;
		}
		
		if (mi_true_positives + mi_false_positives != 0) {
			mi_specificity = float(mi_true_positives) / float(mi_true_positives + mi_false_positives);
		} else {
			mi_specificity = 0;
		}		
		
		return;
	}
	
	void my_test (int target_alignment_length, int target_number_of_sequences, float duplication_rate, int number_of_interactions, float interaction_strength, float& residual_sensitivity, float& residual_specificity, float& Dunn_sensitivity, float& Dunn_specificity, float& mi_sensitivity, float& mi_specificity) {
		int residual_true_positives (0);
		int residual_false_positives (0);
		int Dunn_true_positives (0);
		int Dunn_false_positives (0);
		int mi_true_positives (0);
		int mi_false_positives (0);
		double mu_residual (0);
		double sigma_residual (0);
		double mu_Dunn (0);
		double sigma_Dunn (0);
		double mu_mi (0);
		double sigma_mi (0);
		
		residual_true_positives=0;
		residual_false_positives=0;
		Dunn_true_positives=0;
		Dunn_false_positives=0;
		my_evolve (target_alignment_length, target_number_of_sequences, duplication_rate, number_of_interactions, interaction_strength);
		tabulate();
		calculate_information();
		for (int i=0;i<alignment_length-1;i++) {
			for (int j=i+1;j<alignment_length;j++) {
				mu_residual = mu_residual + MI_residual_with_gaps.at(i).at(j);
				if (Hxy_with_gaps.at(i).at(j) != 0) {
					mu_Dunn = mu_Dunn + (MI_with_gaps.at(i).at(j) / Hxy_with_gaps.at(i).at(j));
				}
				mu_mi = mu_mi + MI_with_gaps.at(i).at(j);
			}
		}
		mu_residual = mu_residual / double(alignment_length * (alignment_length - 1)/ 2);
		mu_Dunn = mu_Dunn / double(alignment_length * (alignment_length - 1) / 2);
		mu_mi = mu_mi / double(alignment_length * (alignment_length - 1) / 2);
		for (int i=0;i<alignment_length-1;i++) {
			for (int j=i+1;j<alignment_length;j++) {
				sigma_residual = sigma_residual + pow(MI_residual_with_gaps.at(i).at(j) - mu_residual,double(2));
				if (Hxy_with_gaps.at(i).at(j) != 0) {
					sigma_Dunn = sigma_Dunn + pow(MI_with_gaps.at(i).at(j) / Hxy_with_gaps.at(i).at(j) - mu_Dunn,double(2));
				} else {
					sigma_Dunn = sigma_Dunn + pow(0.0 - mu_Dunn,double(2));
					
				}
				sigma_mi = sigma_mi + pow(MI_with_gaps.at(i).at(j) - mu_mi,double(2));
			}
		}
		sigma_residual = sqrt(sigma_residual / double(alignment_length * (alignment_length - 1) / 2 - 1));
		sigma_Dunn = sqrt(sigma_Dunn / double(alignment_length * (alignment_length - 1) / 2 - 1));
		sigma_mi = sqrt(sigma_mi / double(alignment_length * (alignment_length - 1) / 2 - 1));
		
		for (int i=0;i<alignment_length-1;i++) {
			for (int j=i+1;j<alignment_length;j++) {
				if (MI_residual_with_gaps.at(i).at(j) > mu_residual + 4 * sigma_residual) {
					if ((j-i == number_of_interactions)&&(j < 2 * number_of_interactions)) {
						residual_true_positives++;
					} else {
						residual_false_positives++;
					}
				}
				if (Hxy_with_gaps.at(i).at(j) != 0) {
					if ((MI_with_gaps.at(i).at(j) / Hxy_with_gaps.at(i).at(j)) > mu_Dunn + 4 * sigma_Dunn) {
						if ((j-i == number_of_interactions)&&(j < 2 * number_of_interactions)) {
							Dunn_true_positives++;
						} else {
							Dunn_false_positives++;
						}
					}
				} else {
					if (0.0 > mu_Dunn + 4 * sigma_Dunn) {
						if ((j-i == number_of_interactions)&&(j < 2 * number_of_interactions)) {
							Dunn_true_positives++;
						} else {
							Dunn_false_positives++;
						}
					}					
				}
				if (MI_with_gaps.at(i).at(j) > mu_mi + 4 * sigma_mi) {
					if ((j-i == number_of_interactions)&&(j < 2 * number_of_interactions)) {
						mi_true_positives++;
					} else {
						mi_false_positives++;
					}
				}
			}				
		}
		
		residual_sensitivity = float(residual_true_positives) / float(number_of_interactions);
		Dunn_sensitivity = float(Dunn_true_positives) / float(number_of_interactions);
		mi_sensitivity = float(mi_true_positives) / float(number_of_interactions);
		
		if (residual_true_positives + residual_false_positives != 0) {
			residual_specificity = float(residual_true_positives) / float(residual_true_positives + residual_false_positives);
		} else {
			residual_specificity = 0;
		}
		
		if (Dunn_true_positives + Dunn_false_positives != 0) {
			Dunn_specificity = float(Dunn_true_positives) / float(Dunn_true_positives + Dunn_false_positives);
		} else {
			Dunn_specificity = 0;
		}
		
		if (mi_true_positives + mi_false_positives != 0) {
			mi_specificity = float(mi_true_positives) / float(mi_true_positives + mi_false_positives);
		} else {
			mi_specificity = 0;
		}		
		return;
	}
	
	
};

int main (int argc, char * const argv[]) {
	string directory_path, PF_file_name, MSA_file_name, subdirectory, PDB_id, PDB_file_name, Stuff, spacer, temp;
	ifstream Stuff_file;
	MSA my_MSA;
	vector<float> residual_specificity (51,0);
	vector<float> residual_sensitivity (51,0);
	vector<float> Dunn_specificity (51,0);
	vector<float> Dunn_sensitivity (51,0);
	vector<float> mi_specificity (51,0);
	vector<float> mi_sensitivity (51,0);
	char have_structure;
	int sequence_number,structure_offset;
	float rsen,rspe,dsen,dspe,msen,mspe,trsen,trspe,tdsen,tdspe,tmsen,tmspe;
	
	/*cout << "Which file?\n";
	cin >> MSA_file_name;
	my_MSA.get_from_file (MSA_file_name);
	my_MSA.tabulate();
	my_MSA.calculate_information();
	//subdirectory=directory_path+"/"+PF_file_name;
	//my_MSA.perspective(581);
	//cout << "Which file?\n";
	//cin >> PDB_file_name;
	//my_MSA.read_PDB_file(PDB_file_name,5);
	cout << "Which directory?\n";
	cin >> subdirectory;
	//PDB_id="1q1l";
	//my_MSA.print_perspective_information_to_directory(subdirectory,PDB_id);
	//my_MSA.print_distances_to_directory(subdirectory);	
	my_MSA.print_information_to_directory(subdirectory);
	//my_MSA.clear();
	
	/**/
	//my_MSA.my_evolve(4,100,0.5,1,1);
	cout << "Which directory?\n";
	cin >> directory_path;
	//directory_path="/Users/knowledgemonger/Desktop/Pfam/Xcode";
	//cout << "Which file?\n";
	//cin >> Stuff;
	//Stuff="/Users/knowledgemonger/Desktop/Pfam/Xcode/Stuctures.txt";
	Stuff="/Users/knowledgemonger/Desktop/CoevolutionProject/Pfam_Full_copy/Bad_Pfam_Input.txt";
	Stuff_file.open(Stuff.c_str());
	//cout << "Which PF number?\n";

	Stuff_file >> PF_file_name;
	cout << PF_file_name << "\n";
	while (PF_file_name!="END"){
		MSA_file_name=directory_path+"/"+PF_file_name+"/MSA.txt";
		my_MSA.get_from_file (MSA_file_name);
		my_MSA.tabulate();
		//my_MSA.calculate_information();
		my_MSA.calculate_other_measures();
		subdirectory=directory_path+"/"+PF_file_name;
		//cout << "Is there a Structure?[y/n]\n";

		Stuff_file >> have_structure;
		//cout << have_structure << "\n";

		if (have_structure=='y' | have_structure=='Y') {
			//cout << "which sequence has a structure?\n";
			Stuff_file >> sequence_number;
			//cout << sequence_number << "\n";
			//my_MSA.perspective(sequence_number);
			//cout << "Which PDB file?\n";
			Stuff_file >> PDB_id;
			//cout << PDB_id << "\n";
			//cout << "What is the structure offset?\n";
			Stuff_file >> structure_offset;
			//cout << structure_offset << "\n";
			for (int q=1; q <= structure_offset; q++){
				Stuff_file >> temp;
				Stuff_file >> temp;
			}
			PDB_file_name = subdirectory + "/" + PDB_id + ".pdb";
			spacer=" ";
			//my_MSA.read_secondary_structure(PDB_file_name, spacer,structure_offset);
			//my_MSA.print_secondary_structure_to_directory(subdirectory);
			//my_MSA.read_PDB_file(PDB_file_name,structure_offset);
			//my_MSA.print_perspective_information_to_directory(subdirectory,PDB_id);
			//my_MSA.print_distances_to_directory(subdirectory);
		}		
		//my_MSA.print_information_to_directory(subdirectory);
		my_MSA.print_other_to_directory(subdirectory);
		my_MSA.clear();
		//cout << "Which PF number?\n";
		Stuff_file >> PF_file_name;
		cout << PF_file_name << "\n";

	}
	//my_MSA.display_information();
	//my_MSA.perspective(1);
	//my_MSA.display_perspective_information();
	//my_MSA.display_gapped_perspective_information();
	
	/*cout << "Is there a Structure?[y/n]\n";
	cin >> have_structure;
	if (have_structure=='y' | have_structure=='Y') {
		cout << "which sequence has a structure?\n";
		cin >> sequence_number;
		my_MSA.perspective(sequence_number);
		cout << "Which PDB file?\n";
		cin >> PDB_file_name;
		cout << "What is the structure offset?\n";
		cin >> structure_offset;
		my_MSA.read_PDB_file(PDB_file_name,structure_offset);
		my_MSA.print_distances_to_directory(subdirectory);
	}
	
	/*
	 my_MSA.my_test(100,500,0.05,25,0.5,rsen,rspe,dsen,dspe,msen,mspe);
	 my_MSA.display();
	 my_MSA.display_information(1);
	 cout << "\n" << rsen << " " << rspe << " " << dsen << " " << dspe << " " << msen << " " << mspe << "\n";
	 /*
	  for (int mu = 0; mu <= 50; mu = mu + 1) {
		  trsen = 0.0;
		  trspe = 0.0;
		  tdsen = 0.0;
		  tdspe = 0.0;
		  tmsen = 0.0;
		  tmspe = 0.0;
		  for (int n=0; n<50 ; n++) {
			  my_MSA.test_normalization(100,100,0.01,pow(float(10),float(-mu)/10),25,rsen,rspe,dsen,dspe,msen,mspe);
			  trsen = trsen + rsen;
			  trspe = trspe + rspe;
			  tdsen = tdsen + dsen;
			  tdspe = tdspe + dspe;
			  tmsen = tmsen + msen;
			  tmspe = tmspe + mspe;
			  
		  }
		  residual_sensitivity.at(mu)= trsen/float(50);
		  residual_specificity.at(mu)= trspe/float(50);
		  Dunn_sensitivity.at(mu)= tdsen/float(50);
		  Dunn_specificity.at(mu)= tdspe/float(50);
		  mi_sensitivity.at(mu)= tmsen/float(50);
		  mi_specificity.at(mu)= tmspe/float(50);
	  }
	  cout << "\nresidual sensitivity\n";
	  for (int mu = 0; mu <= 50; mu = mu + 1) {
		  cout << residual_sensitivity.at(mu) <<" ";
	  }
	  cout << "\nDunn sensitivity\n";
	  for (int mu = 0; mu <= 50; mu = mu + 1) {
		  cout << Dunn_sensitivity.at(mu) << " ";
	  }
	  cout << "\nMI sensitivity\n";
	  for (int mu = 0; mu <= 50; mu = mu + 1) {
		  cout << mi_sensitivity.at(mu) << " ";
	  }	
	  cout << "\n\nresidual specificity\n";
	  for (int mu = 0; mu <= 50; mu = mu + 1) {
		  cout << residual_specificity.at(mu) << " ";
	  }
	  cout << "\nDunn specificity\n";
	  for (int mu = 0; mu <= 50; mu = mu + 1) {
		  cout << Dunn_specificity.at(mu) << " ";
	  }
	  cout << "\nMI specificity\n";
	  for (int mu = 0; mu <= 50; mu = mu + 1) {
		  cout << mi_specificity.at(mu) << " ";
	  }*/
	 return 0;
}
