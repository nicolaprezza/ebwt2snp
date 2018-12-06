/*
 * bwt_merger.hpp
 *
 *  Created on: Dec 5, 2018
 *      Author: nico
 *
 * Merges two compressed BWTs of length n1 and n2 of two string collections. Always computes the document array DA (n1+n2 bits).
 *
 * Assumption: terminators are all represented with special character '$'
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

		{

			stack<pair<sa_leaf, sa_leaf> > S;
			S.push({bwt1->first_leaf(), bwt2->first_leaf()});

			while(not S.empty()){

				auto L = S.top();
				S.pop();

				sa_leaf L1 = L.first;
				sa_leaf L2 = L.second;

				uint64_t start1 = L1.rn.first + L2.rn.first;//start position of first interval in merged intervals
				uint64_t start2 = L2.rn.first + L1.rn.second;//start position of second interval in merged intervals
				uint64_t end = L1.rn.second + L2.rn.second;//end position of merged intervals

				for(uint64_t i = start1; i<start2; ++i) DA[i] = 0;
				for(uint64_t i = start2; i<end; ++i) DA[i] = 1;

				auto leaves1 = bwt1->next_leaves(L1);
				auto leaves2 = bwt2->next_leaves(L2);

				vector<pair<sa_leaf, sa_leaf> > non_empty_leaves;

				/*
				 * Keep pair of leaves {l1,l2} such that at least one of the two is non-empty: in the merged BWT,
				 * these will generate a non-empty leave
				 */
				for(int c=0;c<4;++c){

					if(leaf_size(leaves1[c])>0 or leaf_size(leaves2[c])>0)
						non_empty_leaves.push_back({leaves1[c],leaves2[c]});

					/*
					 * sort non-empty leaves in ascending order of combined range size.
					 */
					std::sort( non_empty_leaves.begin( ), non_empty_leaves.end( ), [ ]( const pair<sa_leaf, sa_leaf>& lhs, const pair<sa_leaf, sa_leaf>& rhs )
					{
					   return range_length(lhs.first.rn)+range_length(lhs.second.rn) < range_length(rhs.first.rn)+range_length(rhs.second.rn);
					});

					/*
					 * push non-empty leaves in descending order of combined range size.
					 */

					for(int i=non_empty_leaves.size()-1;i>=0;--i) S.push(non_empty_leaves[i]);

				}



			}

		}

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
			for(uint64_t i=0;i<n;++i){

				uint8_t c = merged_bwt_at(i);
				out.write((char*)&c,sizeof(c));

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
