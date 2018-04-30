/*
 * read64_t.hpp
 *
 *  Created on: Apr 23, 2018
 *      Author: nico
 *
 * Implements a packed read of at most 64 bases, with fast (constant-time) operations
 * such as Hamming distance, access, LCP, substring
 */

//#include "include.hpp"

#ifndef READ64_T_HPP_
#define READ64_T_HPP_

using namespace std;

class read64_t{

public:

	read64_t(){}

	read64_t(string s){

		assert(s.length()<=64);

		int i=0;

		while(i<s.length()){

			R = R << 2;

			R += base_to_int(s[i]);

			++i;

		}

		while(i<64){
			R = R << 2;
			++i;
		}

		len = s.length();

	}

	read64_t(uint128_t bits, int l){

		uint128_t M = ~((uint128_t(1) << (128-(2*l))) - 1);
		bits &= M;

		R = bits;
		len = l;

	}

	char operator[](int i){

		assert(i<len);
		return int_to_base((R >> (2*(64-i-1))) & uint128_t(3));

	}

	string to_string(){

		string s;

		for(int i=0;i<len;++i) s += operator[](i);

		return s;

	}

	int length(){
		return len;
	}

	/*
	 * Hamming distance until the end of the shortest of the two reads
	 */
	int dH(read64_t r2){

		auto ir1 = R;
		auto ir2 = r2.R;

		int D = 0;

		int minlen = len < r2.length() ? len : r2.length();

		uint128_t M = ~((uint128_t(1) << (128-(2*minlen))) - 1);

		ir1 &= M;
		ir2 &= M;

		uint128_t diff1 = (ir1^ir2) & MASK;
		uint128_t diff2 = ((ir1^ir2)>>1) & MASK;

		return popcount128(diff1 | diff2);

	}

	/*
	 * substring
	 */
	read64_t substr(int start, int l){

		assert(start < len);

		l = start+l <= len ? l : len - start;

		uint128_t sub = R;

		sub = sub << (start*2);

		return {sub,l};

	}

	/*
	 * substring
	 */
	read64_t substr(int start){

		assert(start < len);

		int l = len - start;

		uint128_t sub = R;

		sub << (start*2);

		return {sub,l};

	}

	/*
	 * LCP
	 */
	int LCP(read64_t r2){

		int minlen = len < r2.length() ? len : r2.length();
		int lcp = (clz_u128(R^r2.R)/2);

		lcp = minlen < lcp ? minlen : lcp;

		return lcp;

	}

	/*
	 * is this read strictly less than r2?
	 */
	bool less_than(read64_t r2){

		int minlen = len < r2.length() ? len : r2.length();
		auto lcp = LCP(r2);

		return lcp < minlen ? R<r2.R : len < r2.length();

	}

private:

	//the packed bases
	uint128_t R = 0;

	//length
	int len = 0;

	static constexpr uint128_t MASK = (uint128_t(0x5555555555555555)<<64) + uint128_t(0x5555555555555555);

};


#endif /* READ64_T_HPP_ */
