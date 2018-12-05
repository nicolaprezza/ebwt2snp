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

using namespace std;
using namespace sdsl;

template<class bwt_t, typename lcp_int_t>
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
	bwt_merger(bwt_t * bwt1, bwt_t * bwt2, bool compute_lcp = false, bool compute_minima = false, uint64_t lcp_threshold = 0){

		this->bwt1 = bwt1;
		this->bwt2 = bwt2;
		this->lcp_threshold = lcp_threshold;

		n = bwt1.size() + bwt2.size();

		DA = bit_vector(n);





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
	 * store to file the specified components. Adds extensions .merged.bwt, .merged.da, .merged.lcp
	 */
	void save_to_file(string base_path, bool store_bwt = false, bool store_da = false, bool store_lcp = false){

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

	bwt_t * bwt1 = NULL;
	bwt_t * bwt2 = NULL;

};



#endif /* INTERNAL_BWT_MERGER_HPP_ */
