/*
 * dna_string.hpp
 *
 *  Created on: Dec 6, 2018
 *      Author: nico
 *
 *  Optimized string with rank on DNA alphabet: {A,C,G,T,TERM}, where the terminator TERM is defined in include.hpp
 *
 *  One access or a block of rank for all 4 letters A,C,G,T causes only 1 cache miss in the worst case
 *
 *  Max string length: 2^48
 *
 *  Supports very efficient (1 cache miss) combined rank for (A,C,G,T), and (1 cache miss) single rank for other characters ($ + N)
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
 *    Size of the string: 4,876190476 bits / base
 *
 */

#ifndef INTERNAL_DNA_STRING_HPP_
#define INTERNAL_DNA_STRING_HPP_

#define BLOCK_SIZE 105 //characters per block of 512 bits
#define SUB_BLOCK_SIZE 21
#define WORDS_PER_BLOCK 8

#include "include.hpp"

struct p_rank{

	uint64_t A;
	uint64_t C;
	uint64_t G;
	uint64_t T;

};

class dna_string{

public:

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

	}

	p_rank parallel_rank(uint64_t i){

		assert(i<=n);

		return {};

	}

	~dna_string(){

		delete data;

	}

private:

	//data aligned with blocks of 64 bytes = 512 bits
	alignas(64) uint64_t * data;

	uint64_t n_blocks_512 = 0;
	uint64_t n_blocks_64 = 0;
	uint64_t n = 0;

};


#endif /* INTERNAL_DNA_STRING_HPP_ */
