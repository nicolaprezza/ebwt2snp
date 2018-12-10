/*
 * dna_bwt.hpp
 *
 *  Created on: Dec 3, 2018
 *      Author: nico
 *
 *  BWT for big files: built in chunks. The class is a template on the string type (could be run-length, Huffman...)
 *
 */

#include "include.hpp"
#include "rle_string.hpp"
#include "huff_string.hpp"
#include "dna_string.hpp"

#ifndef INTERNAL_DNA_BWT_HPP_
#define INTERNAL_DNA_BWT_HPP_

using namespace sdsl;

template<class str_type>
class dna_bwt{

public:

	dna_bwt(){};

	/*
	 * constructor path of a BWT file containing the BWT in ASCII format
	 */
	dna_bwt(string path){

		//build F column
		F = vector<uint64_t>(256,0);

		n = uint64_t(filesize(path));

		BWT = str_type(path);

		//build F column
		for(uint64_t i=0;i<n;++i){

			assert(BWT[i] < 256);
			F[BWT[i]]++;

		}

		for(uint64_t i=255;i>0;--i)
			F[i] = F[i-1];

		F[0] = 0;

		for(uint64_t i=1;i<256;++i)
			F[i] += F[i-1];

	}

	/*
	 * get full BWT range
	 */
	range_t full_range(){

		//right-exclusive range
		return {0,size()};

	}

	/*
	 * left-extend range by all 4 nucleotides
	 */
	p_range parallel_LF(range_t rn){

		assert(rn.second >= rn.first);

		//number of A,C,T,G before start of interval
		p_rank start = BWT.parallel_rank(rn.first);

		//number of A,C,T,G before end of interval (last position of interval included)
		p_rank end;

		if(rn.second>rn.first)
			end	= BWT.parallel_rank(rn.second);
		else
			end = start;

		assert(start <= end);

		p_rank f = {F['A'],F['C'],F['G'],F['T']};
		p_rank l = f + start;
		p_rank r = f + end;

		assert(r.A < l.C);
		assert(r.C < l.G);
		assert(r.G < l.T);

		return fold_ranks(l,r);

	}

	//backward navigation of the BWT
	uint64_t LF(uint64_t  i){

		auto c = operator[](i);
		return F[c] + rank(i,c);

	}

	/*
	 * access column F at position i
	 */
	uint8_t F_at(uint64_t i){

		uint64_t c = (upper_bound(F.begin(),F.end(),i) - F.begin()) - 1;
		assert(c<256);
		assert(i>=F[c]);

		return uint8_t(c);

	}

	/*
	 * Return BWT range of character c
	 */
	range_t get_char_range(uint8_t c){

		return {F[c], F[c+1]};

	}

	/*
	 * Return BWT range of pattern P
	 */
	range_t count(string &P){

		auto range = full_range();
		uint64_t m = P.size();

		for(uint64_t i=0;i<m and range.second>range.first;++i)
			range = LF(range,P[m-i-1]);

		return range;

	}

	/*
	 * Return number of occurrences of P in the text
	 */
	uint64_t occ(string &P){

		auto rn = count(P);

		return rn.second>=rn.first ? (rn.second-rn.first)+1 : 0;

	}


	uint8_t operator[](uint8_t i){

		return BWT[i];

	}

	/*
	 * number of c before position i excluded
	 */
	uint64_t rank(uint64_t i, uint8_t c){

		assert(i<=n);

		return BWT.rank(i,c);

	}

	/*
	 * return number of occurrences of A,C,T,G in the prefix of length i of the text. At most 1 cache miss!
	 */
	p_rank parallel_rank(uint64_t i){

		return BWT.parallel_rank(i);

	}

	uint64_t size(){

		assert(n == BWT.size());
		return n;

	}

	uint64_t serialize(std::ostream& out){

		uint64_t w_bytes = 0;

		out.write((char*)&n,sizeof(n));
		out.write((char*)F.data(),256*sizeof(uint64_t));

		w_bytes += sizeof(n) + sizeof(uint64_t)*256;

		w_bytes += BWT.serialize(out);

		return w_bytes;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&n,sizeof(n));

		F = vector<uint64_t>(256, 0);
		in.read((char*)F.data(),256*sizeof(uint64_t));

		BWT.load(in);

	}

	void save_to_file(string path){

		std::ofstream out(path);
		serialize(out);
		out.close();

	}

	/*
	 * path = path of an index file
	 */
	void load_from_file(string path){

		std::ifstream in(path);
		load(in);
		in.close();

	}


	/*
	 * functions for suffix tree navigation
	 */

	sa_node root(){

		return {
			F[TERM],
			F['A'],
			F['C'],
			F['G'],
			F['T'],
			F['T'+1],
			0
		};

	}

	/*
	 * depth = LCP inside the leaf.
	 */
	sa_leaf first_leaf(){

		return {{F[TERM], F[uint8_t(TERM)+1]}, 0};

	}

	/*
	 * Input: suffix array node N representing right-maximal string W
	 * Output: vector of sa nodes A.W,...,T.W reached
	 * following Weiner links from N. Nodes could be implicit!
	 */
	p_node weiner(sa_node N){

		p_rank before_TERM;
		p_rank before_A;
		p_rank before_C;
		p_rank before_G;
		p_rank before_T;
		p_rank before_end;

		before_TERM = parallel_rank(N.first_TERM);

		if(N.first_A == N.first_TERM) before_A = before_TERM;
		else before_A = parallel_rank(N.first_A);

		if(N.first_C == N.first_A) before_C = before_A;
		else before_C = parallel_rank(N.first_C);

		if(N.first_G == N.first_C) before_G = before_C;
		else before_G = parallel_rank(N.first_G);

		if(N.first_T == N.first_G) before_T = before_G;
		else before_T = parallel_rank(N.first_T);

		if(N.last == N.first_T) before_end = before_T;
		else before_end = parallel_rank(N.last);

		return {
			{before_TERM.A, before_A.A, before_C.A, before_G.A, before_T.A, before_end.A, N.depth+1},
			{before_TERM.C, before_A.C, before_C.C, before_G.C, before_T.C, before_end.C, N.depth+1},
			{before_TERM.G, before_A.G, before_C.G, before_G.G, before_T.G, before_end.G, N.depth+1},
			{before_TERM.T, before_A.T, before_C.T, before_G.T, before_T.T, before_end.T, N.depth+1}
		};

	}

private:

	void count_chars(string & s){

		for(auto c:s) F[c]++;

	}

	uint64_t n = 0;//BWT length
	vector<uint64_t> F;
	str_type BWT;

};

typedef dna_bwt<dna_string> dna_bwt_t;

#endif /* INTERNAL_DNA_BWT_HPP_ */
