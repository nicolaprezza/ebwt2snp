/*
 * bwt_merger.hpp
 *
 *  Created on: Dec 5, 2018
 *      Author: nico
 *
 * Merges two compressed BWTs of length n1 and n2 of two string collections. Always computes the document array DA (n1+n2 bits).
 *
 * Assumption: terminators are all represented with special character TERM=#
 *
 * If just DA is computed, uses only n1+n2 bits of space on top of the input (compressed) BWTs
 * (plus a small stack of O(sigma*log n) words during merge)
 *
 * Optionally can compute also any subset of these three arrays:
 *
 * - LCP array of the merged BWT. Requires (n1+n2)*size_int bytes, where size_int is the size in
 *   Bytes of an integer (template parameter lcp_int_t)
 *
 * - the LCP-minima array: a bit-array MIN such that MIN[i] = 1 iff LCP[i] is a local minima. Requires n1+n2 additional bits.
 *
 * - the LCP-threshold array: a bit-array T such that T[i] = 1 iff LCP[i]>=lcp_threshold. Requires n1+n2 additional bits.
 *
 * Based on the SA-nodes enumeration algorithm described in
 *
 * "Linear-time string indexing and analysis in small space", Djamal Belazzougui, Fabio Cunial, Juha Kärkkäinen, and Veli Mäkinen
 *
 */

#ifndef INTERNAL_BWT_MERGER_HPP_
#define INTERNAL_BWT_MERGER_HPP_

#include "dna_bwt.hpp"
#include "bwt.hpp"
#include "sdsl/int_vector.hpp"
#include <stack>

using namespace std;
using namespace sdsl;

template<class bwt_t1, class bwt_t2, typename lcp_int_t>
class bwt_merger{

public:

	/*
	 * merge bwt1 and bwt2. Parameters:
	 *
	 * bwt1, bwt2: the two BWTs.
	 * compute_lcp = if true, compute the LCP array of the merged BWT.
	 * out_da = store to file document array.
	 * compute_minima: if true, compute array with 1 at the local minima of the LCP array
	 * lcp_threshold: if >0, compute a boolean array T[i] = 1 iff LCP[i] >= lcp_threshold.
	 *
	 */
	bwt_merger(bwt_t1 * bwt1, bwt_t2 * bwt2, bool compute_lcp = false, bool out_da = false, bool compute_minima = false, uint64_t lcp_threshold = 0){

		this->out_da = out_da;
		this->bwt1 = bwt1;
		this->bwt2 = bwt2;
		this->lcp_threshold = lcp_threshold;

		n = bwt1->size() + bwt2->size();

		DA = bit_vector(n);

		if(compute_lcp){

			LCP = vector<lcp_int_t>(n);
			LCP[0] = 0;
		}

		if(lcp_threshold>0) T = vector<bool>(n);
		if(compute_minima) MIN = vector<bool>(n);

		/*
		 * FIRST PASS: NAVIGATE LEAVES AND MERGE BWTs (i.e. build DA). If enabled, compute LCP inside leaves.
		 */

		cout << "\nNow navigating suffix tree leaves to merge BWTs and computing internal LCP values (if LCP enabled)." << endl;

		uint64_t m = 0;//number of entries filled in DA
		uint64_t leaves = 0;//number of visited nodes
		uint64_t max_stack = 0;
		uint64_t lcp_values = 1;//number of filled LCP values

		{

			stack<pair<sa_leaf, sa_leaf> > S;
			S.push({bwt1->first_leaf(), bwt2->first_leaf()});

			int last_perc = -1;
			int perc = 0;

			while(not S.empty()){

				perc = (100*m)/n;

				if(perc > last_perc){

					cout << perc << "% done" << endl;
					last_perc = perc;

				}

				pair<sa_leaf, sa_leaf> L = S.top();
				S.pop();
				leaves++;

				max_stack = S.size() > max_stack ? S.size() : max_stack;

				sa_leaf L1 = L.first;
				sa_leaf L2 = L.second;

				uint64_t start1 = L1.rn.first + L2.rn.first;//start position of first interval in merged intervals
				uint64_t start2 = L2.rn.first + L1.rn.second;//start position of second interval in merged intervals
				uint64_t end = L1.rn.second + L2.rn.second;//end position of merged intervals

				assert(leaf_size(L)>0);
				assert(end>start1);

				for(uint64_t i = start1; i<start2; ++i){
					DA[i] = 0;
					m++;
				}

				for(uint64_t i = start2; i<end; ++i){
					DA[i] = 1;
					m++;
				}

				assert(L1.depth==L2.depth);

				if(compute_lcp){

					for(uint64_t i = start1+1; i<end; ++i){

						LCP[i] = L1.depth;
						lcp_values++;

					}

				}

				assert(m<=n);

				p_range ext_1 = bwt1->parallel_LF(L1.rn);
				p_range ext_2 = bwt2->parallel_LF(L2.rn);

				if(range_length(ext_1.A)>0 or range_length(ext_2.A)>0) S.push({{ext_1.A, L1.depth+1},{ext_2.A, L2.depth+1}});
				if(range_length(ext_1.C)>0 or range_length(ext_2.C)>0) S.push({{ext_1.C, L1.depth+1},{ext_2.C, L2.depth+1}});
				if(range_length(ext_1.G)>0 or range_length(ext_2.G)>0) S.push({{ext_1.G, L1.depth+1},{ext_2.G, L2.depth+1}});
				if(range_length(ext_1.T)>0 or range_length(ext_2.T)>0) S.push({{ext_1.T, L1.depth+1},{ext_2.T, L2.depth+1}});

			}
		}

		cout << "Computed " << m << "/" << n << " DA entries." << endl;

		if(compute_lcp)
		cout << "Computed " << lcp_values << "/" << n << " LCP values." << endl;

		cout << "Max stack depth = " << max_stack << endl;
		cout << "Processed " << leaves << " suffix-tree leaves." << endl;

		assert(m==n);

		//add rank support to DA for random access to merged BWT
		rank1 = bit_vector::rank_1_type(&DA);

		if(compute_lcp){

			cout << "\nNow navigating suffix tree nodes to compute remaining LCP values." << endl;

			uint64_t nodes = 0;//visited ST nodes
			max_stack = 0;

			stack<pair<sa_node, sa_node> > S;
			S.push({bwt1->root(), bwt2->root()});

			int last_perc = -1;
			int perc = 0;

			while(not S.empty()){

				perc = (100*lcp_values)/n;

				if(perc > last_perc){

					cout << "Computed LCP values: " << perc << "% " << endl;
					last_perc = perc;

				}

				max_stack = S.size() > max_stack ? S.size() : max_stack;

				pair<sa_node, sa_node> N = S.top();
				S.pop();
				nodes++;

				sa_node N1 = N.first;
				sa_node N2 = N.second;

				sa_node merged = merge_nodes(N1, N2);

				if(merged.first_A-merged.first_TERM > 0){
					LCP[merged.first_A] = merged.depth;
					lcp_values++;
				}
				if(merged.first_C-merged.first_A > 0){
					LCP[merged.first_C] = merged.depth;
					lcp_values++;
				}
				if(merged.first_G-merged.first_C > 0){
					LCP[merged.first_G] = merged.depth;
					lcp_values++;
				}
				if(merged.first_T-merged.first_G > 0){
					LCP[merged.first_T] = merged.depth;
					lcp_values++;
				}

				p_node left_exts1 = bwt1->weiner(N1);
				p_node left_exts2 = bwt2->weiner(N2);

				/*cout << "weiners:" << endl;
				print_node(N1);
				print_node(N2);
				cout << endl;
				print_nodes(left_exts1);
				cout << endl;
				print_nodes(left_exts2);
				cout << endl;*/

				if(number_of_children(left_exts1.A, left_exts2.A)>1) S.push({left_exts1.A, left_exts2.A});
				if(number_of_children(left_exts1.C, left_exts2.C)>1) S.push({left_exts1.C, left_exts2.C});
				if(number_of_children(left_exts1.G, left_exts2.G)>1) S.push({left_exts1.G, left_exts2.G});
				if(number_of_children(left_exts1.T, left_exts2.T)>1) S.push({left_exts1.T, left_exts2.T});

			}

			cout << "Computed " << lcp_values << "/" << n << " LCP values." << endl;
			cout << "Max stack depth = " << max_stack << endl;
			cout << "Processed " << nodes << " suffix-tree nodes." << endl;

		}

	}

	/*
	 * access merged BWT at position i
	 */
	uint8_t merged_bwt_at(uint64_t i){

		return DA[i]==0 ? bwt1->operator[](i-rank1(i)) : bwt2->operator[](rank1(i));

	}

	/*
	 * store to file the components that have been built (BWT, DA for sure; optionally, LCP). Adds extensions .merged.bwt, .merged.da, .merged.lcp
	 */
	void save_to_file(string base_path){

		string bwt_path = base_path;
		bwt_path.append(".merged.ebwt");

		string da_path = base_path;
		da_path.append(".merged.da");

		string lcp_path = base_path;
		lcp_path.append(".merged.lcp");

		{
			std::ofstream out(bwt_path);

			uint64_t rank1 = 0;
			for(uint64_t i=0;i<n;++i){

				uint8_t c = DA[i]==0 ? bwt1->operator[](i-rank1) : bwt2->operator[](rank1);
				out.write((char*)&c,sizeof(c));

				rank1 += DA[i];

			}
			out.close();
		}

		if(out_da){
			std::ofstream out(da_path);
			for(uint64_t i=0;i<n;++i){

				uint8_t c = DA[i] ? '1' : '0';
				out.write((char*)&c,sizeof(c));

			}
			out.close();
		}

		if(LCP.size()>0){

			std::ofstream out(lcp_path);

			out.write((char*)LCP.data(),LCP.size() * sizeof(lcp_int_t));

			out.close();

		}

	}

	uint64_t lcp_threshold = 0;

	bit_vector DA; //document array
	//rank support on DA
	bit_vector::rank_1_type rank1;

	vector<bool> T;   //LCP-threshold array: 1 where LCP > lcp_threshold
	vector<bool> MIN; //LCP-minima array: MIN[i] = 1 iff LCP[i] is a local minima
	vector<lcp_int_t> LCP;

private:

	uint64_t n = 0;//total size

	bool out_da = false;

	bwt_t1 * bwt1 = NULL;
	bwt_t2 * bwt2 = NULL;

};



#endif /* INTERNAL_BWT_MERGER_HPP_ */
