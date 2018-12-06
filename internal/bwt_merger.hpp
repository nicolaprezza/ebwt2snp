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
	 * compute_minima: if true, compute array with 1 at the local minima of the LCP array
	 * lcp_threshold: if >0, compute a boolean array T[i] = 1 iff LCP[i] >= lcp_threshold.
	 *
	 */
	bwt_merger(bwt_t1 * bwt1, bwt_t2 * bwt2, bool compute_lcp = false, bool compute_minima = false, uint64_t lcp_threshold = 0){

		this->bwt1 = bwt1;
		this->bwt2 = bwt2;
		this->lcp_threshold = lcp_threshold;

		n = bwt1->size() + bwt2->size();

		DA = bit_vector(n);

		if(compute_lcp) LCP = vector<lcp_int_t>(n);
		if(lcp_threshold>0) T = vector<bool>(n);
		if(compute_minima) MIN = vector<bool>(n);

		/*
		 * FIRST PASS: NAVIGATE LEAVES AND MERGE BWTs (i.e. build DA). If enabled, compute LCP inside leaves.
		 */

		uint64_t m = 0;//number of entries filled in DA
		uint64_t nodes = 0;//number of visited nodes
		uint64_t max_stack = 0;

		{

			stack<pair<sa_leaf, sa_leaf> > S;
			S.push({bwt1->first_leaf(), bwt2->first_leaf()});

			while(not S.empty()){

				auto L = S.top();
				S.pop();
				nodes++;

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

				assert(m<=n);
				assert(L1.depth==L2.depth);

				//leaves with extension A
				sa_leaf l1A = { bwt1->LF(L1.rn,'A'), L1.depth+1 };
				sa_leaf l2A = { bwt2->LF(L2.rn,'A'), L2.depth+1 };
				//leaves with extension C
				sa_leaf l1C = { bwt1->LF(L1.rn,'C'), L1.depth+1 };
				sa_leaf l2C = { bwt2->LF(L2.rn,'C'), L2.depth+1 };
				//leaves with extension G
				sa_leaf l1G = { bwt1->LF(L1.rn,'G'), L1.depth+1 };
				sa_leaf l2G = { bwt2->LF(L2.rn,'G'), L2.depth+1 };
				//leaves with extension T
				sa_leaf l1T;
				sa_leaf l2T;

				if(		leaf_size(l1A)+leaf_size(l1C)+leaf_size(l1G) < leaf_size(L1) or
						leaf_size(l2A)+leaf_size(l2C)+leaf_size(l2G) < leaf_size(L2)){

					l1T = { bwt1->LF(L1.rn,'T'), L1.depth+1 };
					l2T = { bwt2->LF(L2.rn,'T'), L2.depth+1 };

					S.push({l1T,l2T});

				}

				if(leaf_size(l1A)>0 or leaf_size(l2A)>0) S.push({l1A,l2A});
				if(leaf_size(l1C)>0 or leaf_size(l2C)>0) S.push({l1C,l2C});
				if(leaf_size(l1G)>0 or leaf_size(l2G)>0) S.push({l1G,l2G});

			}
		}

		cout << m << "/" << n << endl;
		cout << "max stack depth = " << max_stack << endl;

		//assert(m==n);

		//add rank support to DA for random access to merged BWT
		rank1 = bit_vector::rank_1_type(&DA);

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
		bwt_path.append(".merged.bwt");

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

		{
			std::ofstream out(da_path);
			for(uint64_t i=0;i<n;++i){

				uint8_t c = DA[i];
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

	bwt_t1 * bwt1 = NULL;
	bwt_t2 * bwt2 = NULL;

};



#endif /* INTERNAL_BWT_MERGER_HPP_ */
