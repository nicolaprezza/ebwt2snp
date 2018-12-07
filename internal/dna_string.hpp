/*
 * dna_string.hpp
 *
 *  Created on: Dec 6, 2018
 *      Author: nico
 *
 *  Optimized string with rank on DNA alphabet: {A,C,G,T,TERM,N}, where the terminator TERM is defined in include.hpp
 *
 *  One access or a parallel rank for all 4 letters A,C,G,T causes only 1 cache miss in the worst case
 *
 *  Max string length: 2^48
 *
 *  Supports very efficient (1 cache miss) parallel rank for (A,C,G,T), and (1 cache miss) single rank for other characters combined ($ + N)
 *  Note: if both N's and $'s are present, then can only compute combined rank of N + $!
 *
 *
 *  Data is stored and cache-aligned in blocks of 512 bits (64 bytes)
 *
 *  Each block contains:
 *
 *  - in the first 24 bytes = 3 * 64-bit words = 192 bits, 4 integers of 48 bits each storing the partial ranks
 *    before the block for A, C, G, T
 *  - in the next 40 bytes = 5 * 64-bit words = 320 bits, 105 nucleotides, 'N', or terminator character,
 *    encoded using 3-bits per letter. These characters are stored in blocks of 21 per 64-bit word, wasting only 1 bit per word.
 *
 *    Size of the string: < 4,9 bits / base
 *
 */

#ifndef INTERNAL_DNA_STRING_HPP_
#define INTERNAL_DNA_STRING_HPP_

#define BLOCK_SIZE 105 //characters per block of 512 bits
#define SUB_BLOCK_SIZE 21
#define WORDS_PER_BLOCK 8

#include "include.hpp"

class p_rank{

public:

	uint64_t A;
	uint64_t C;
	uint64_t G;
	uint64_t T;

	p_rank operator+(const p_rank& a) const{

		return {
			a.A + A,
			a.C + C,
			a.G + G,
			a.T + T
		};

	}

};

class dna_string{

public:

	/*
	 * constructor from file
	 */
	dna_string(string path){

		n = uint64_t(filesize(path));

		n_blocks_512 = n/BLOCK_SIZE + (n%BLOCK_SIZE != 0);
		n_blocks_64 = n_blocks_512 * 8;

		data = new uint64_t[n_blocks_64];

		for(uint64_t i=0;i<n_blocks_512;++i) data[i] = 0;

		ifstream ifs(path);

		uint8_t i = 0;
		while(not ifs.eof()){

			uint8_t c;
			ifs.read((char*)&c, sizeof(uint8_t));

			if(c != 'A' and c != 'C' and c != 'G' and c != 'T' and c != 'N' and c != TERM){

				cout << "Error while loading file " << path << ": forbidden character '" << c << "'" << endl;
				exit(1);

			}

			set(i++,c);

		}

	}

	/*
	 * constructor: allocates a string of length size filled with 'A's
	 *
	 * after allocation, fill the string using set(), and then build rank support.
	 *
	 */
	dna_string(uint64_t size){

		n = size;

		n_blocks_512 = size/BLOCK_SIZE + (size%BLOCK_SIZE != 0);
		n_blocks_64 = n_blocks_512 * 8;

		data = new uint64_t[n_blocks_64];

		for(uint64_t i=0;i<n_blocks_512;++i) data[i] = 0;

	}

	/*
	 * set the i-th character to c.
	 */
	void set(uint64_t i, uint8_t c){

		//A is 0x0
		uint64_t b = 	(c == 'C')*  uint64_t(1) +
						(c == 'G')*  uint64_t(2) +
						(c == 'T')*  uint64_t(3) +
						(c == TERM)* uint64_t(4) +
						(c == 'N')*  uint64_t(5);

		uint64_t block_number = i / BLOCK_SIZE;
		uint64_t off_in_block = i % BLOCK_SIZE;
		uint8_t sub_block_number = 	off_in_block / SUB_BLOCK_SIZE;
		uint8_t off_in_sub_block = 	off_in_block % SUB_BLOCK_SIZE;

		data[block_number * WORDS_PER_BLOCK + 3 + sub_block_number] += (b << ((SUB_BLOCK_SIZE - off_in_sub_block -1)*3));

	}

	uint8_t operator[](uint64_t i){

		assert(i<n);

		uint64_t block_number = i / BLOCK_SIZE;
		uint64_t off_in_block = i % BLOCK_SIZE;
		uint8_t sub_block_number = 	off_in_block / SUB_BLOCK_SIZE;
		uint8_t off_in_sub_block = 	off_in_block % SUB_BLOCK_SIZE;

		uint64_t b = (data[block_number * WORDS_PER_BLOCK + 3 + sub_block_number] >> ((SUB_BLOCK_SIZE - off_in_sub_block -1)*3)) & 0x7;

		return 	(b == 0)*'A' +
				(b == 1)*'C' +
				(b == 2)*'G' +
				(b == 3)*'T' +
				(b == 4)*TERM +
				(b == 5)*'N';

	}

	void build_rank_support(){

		p_rank r;

		for(uint64_t i = 1; i<n_blocks_512;++i){

			r = r + block_rank(i-1,BLOCK_SIZE);
			set_counters(i,r);

		}

	}

	/*
	 * return number of occurrences of A,C,T,G in the prefix of length i of the text. At most 1 cache miss!
	 */
	p_rank parallel_rank(uint64_t i){

		assert(i<=n);
		uint64_t bl = i/BLOCK_SIZE;
		return get_counters(bl) + block_rank(bl,i%BLOCK_SIZE);

	}

	/*
	 * standard rank. c can be A,C,G,T, or TERM (TERM = any non-dna symbol is counted, which could include also N if present)
	 */
	uint64_t rank(uint64_t i, uint8_t c){

		p_rank pr = parallel_rank(i);

		switch(c){
			case 'A' : return pr.A; break;
			case 'C' : return pr.C; break;
			case 'G' : return pr.G; break;
			case 'T' : return pr.T; break;
			case TERM : return rank_non_dna(i); break;
		}

		return 0;

	}

	/*
	 * return number of non-dna symbols in the prefix of length i of the text. At most 1 cache miss!
	 */
	uint64_t rank_non_dna(uint64_t i){

		assert(i<=n);
		auto r = parallel_rank(i);

		assert(r.A + r.C + r.G + r.T <= i);

		return i - (r.A + r.C + r.G + r.T);

	}

	~dna_string(){

		delete data;

	}

private:

	/*
	 * A 64-bits word represents a block of SUB_BLOCK_SIZE = 21 characters (leftmost bit wasted).
	 * This function counts in parallel the number of A,C,G,T in the prefix of length i of the block.
	 *
	 */
	inline p_rank word_rank(uint64_t x, uint8_t i){

		if(i==0) return {};

		//First, replace with '111' (unused code) the last SUB_BLOCK_SIZE - i characters

		x |= ( (uint64_t(1) << (3*(SUB_BLOCK_SIZE - i)))-1 );

		//x shifted to the left 1 and 2 positions
		uint64_t x1 = x>>1;
		uint64_t x2 = x1>>1;

		return {
			uint64_t(__builtin_popcountll(((~x2) & MASK) & ((~x1) & MASK) & ((~x) & MASK))),
			uint64_t(__builtin_popcountll(((~x2) & MASK) & ((~x1) & MASK) & ((x) & MASK))),
			uint64_t(__builtin_popcountll(((~x2) & MASK) & ((x1) & MASK) & ((~x) & MASK))),
			uint64_t(__builtin_popcountll(((~x2) & MASK) & ((x1) & MASK) & ((x) & MASK)))
		};

	}

	/*
	 * Parallel rank of (A,C,T,G) for the length-j prefix of the i-th block of size BLOCK_SIZE of the text
	 *
	 */
	inline p_rank block_rank(uint64_t i, uint8_t j){

		if(j==0) return {};

		uint64_t start = i*WORDS_PER_BLOCK + 3;//first word containing the block's characters

		uint8_t bl = j/SUB_BLOCK_SIZE; // block where position j lies.
		uint8_t off = j%SUB_BLOCK_SIZE; // prefix length in the bl-th block

		return	word_rank( data[start],   bl == 0 ? off : SUB_BLOCK_SIZE ) +
				(bl < 1 ? p_rank() : word_rank( data[start+1], bl == 1 ? off : SUB_BLOCK_SIZE)) +
				(bl < 2 ? p_rank() : word_rank( data[start+2], bl == 2 ? off : SUB_BLOCK_SIZE)) +
				(bl < 3 ? p_rank() : word_rank( data[start+3], bl == 3 ? off : SUB_BLOCK_SIZE)) +
				(bl < 4 ? p_rank() : word_rank( data[start+4], bl == 4 ? off : SUB_BLOCK_SIZE));

	}

	/*
	 * set counters in the i-th block to r
	 */
	void set_counters(uint64_t i, p_rank r){

		uint64_t start = i*WORDS_PER_BLOCK;//start of region where to write the four 48-bits numbers

		data[start]   |= (r.A << 16);//first 48 bits of data[start] contain r.A

		data[start]   |= (r.C >> 32);//last 16 bits of data[start] contain first 16 bits of r.C
		data[start+1] |= ((r.C & 0x00000000FFFFFFFF) << 32);//first 32 bits of data[start+1] contain last 32 bits of r.C

		data[start+1] |= (r.G >> 16) ; //last 32 bits of data[start+1] contain first 32 bits of r.G
		data[start+2] |= ((r.G & 0x000000000000FFFF) << 48) ; //first 16 bits of data[start+2] contain last 16 bits of r.G

		data[start+2] |= r.T; //last 48 bits of data[start+2] contain the 48 bits of r.G

	}

	/*
	 * get counters of the i-th block
	 */
	inline p_rank get_counters(uint64_t i){

		uint64_t start = i*WORDS_PER_BLOCK;//start of region where to read the four 48-bits numbers

		return {	data[start]>>16,
					((data[start]&0x000000000000FFFF)<<32) | (data[start+1]>>32),
					((data[start+1]&0x00000000FFFFFFFF)<<16) | (data[start+2]>>48),
					data[start+2]&0x0000FFFFFFFFFFFF
		};

	}

	/*
	 * MASK =
	 *    0001 0010 0100 1001 0010 0100 1001 0010 0100 1001 0010 0100 1001 0010 0100 1001 =
	 * 0x 1    2    4    9    2    4    9    2    4    9    2    4    9    2    4    9
	 */
	static const uint64_t MASK = 0x1249249249249249;

	//data aligned with blocks of 64 bytes = 512 bits
	alignas(64) uint64_t * data = NULL;

	uint64_t n_blocks_512 = 0;
	uint64_t n_blocks_64 = 0;
	uint64_t n = 0;

};


#endif /* INTERNAL_DNA_STRING_HPP_ */
