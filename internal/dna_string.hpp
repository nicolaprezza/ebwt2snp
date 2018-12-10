/*
 * dna_string.hpp
 *
 *  Created on: Dec 6, 2018
 *      Author: nico
 *
 *  Optimized string with rank on DNA alphabet: {A,C,G,T,TERM}, where the terminator TERM is defined in include.hpp
 *
 *  One access or a parallel rank for all 4 letters A,C,G,T causes only 1 cache miss in the worst case
 *
 *  Max string length: 2^64
 *
 *  Supports very efficient (1 cache miss) parallel rank for (A,C,G,T), and (1 cache miss) single rank for TERM
 *
 *  Data is stored and cache-aligned in blocks of 512 bits (64 bytes)
 *
 *  Size of the string: 4n bits, where n = string length
 *
 *
 */

#ifndef INTERNAL_DNA_STRING_HPP_
#define INTERNAL_DNA_STRING_HPP_

#define SUPERBLOCK_SIZE 0x100000000 	//number of characters in a superblock = 2^32 characters
#define BLOCKS_PER_SUPERBLOCK 33554432	//blocks in a superblock
#define BYTES_PER_SUPERBLOCK 2147483648	//bytes in a superblock
#define BLOCK_SIZE 128 					//number of characters inside a block
#define BYTES_PER_BLOCK 64				//bytes in a block of 512 bits
#define ALN 64

#include "include.hpp"

class dna_string{

public:

	dna_string(){}

	/*
	 * constructor from ASCII file
	 */
	dna_string(string path){

		n = uint64_t(filesize(path));

		n_superblocks = (n+1)/SUPERBLOCK_SIZE + ((n+1)%SUPERBLOCK_SIZE != 0);
		n_blocks = (n+1)/128 + ((n+1)%128 != 0);
		nbytes = (n_blocks * BYTES_PER_BLOCK);//number of bytes effectively filled with data

		superblock_ranks = vector<p_rank>(n_superblocks);

		/*
		 * this block of code ensures that data is aligned by 64 bytes = 512 bits
		 */
		memory = vector<uint8_t>(nbytes+ALN,0);
		data = memory.data();
		while(uint64_t(data) % ALN != 0) data++;

		cout << "alignment of data: " << (void*)data << endl;

		{

			ifstream ifs(path);

			for(uint64_t i = 0; i<n; ++i){

				uint8_t c;
				ifs.read((char*)&c, sizeof(uint8_t));

				if(c != 'A' and c != 'C' and c != 'G' and c != 'T' and c != TERM){

					cout << "Error while loading file " << path << ": forbidden character '" << c << "'" << endl;
					exit(1);

				}

				set(i,c);
				assert(operator[](i) == c);

			}

		}

		assert(check_content(path));
		build_rank_support();

	}

	/*
	 * constructor: allocates a string of length size filled with 'A's
	 *
	 * after allocation, fill the string using set(), and then build rank support.
	 *
	 */
	dna_string(uint64_t n){

		this->n = n;

		n_superblocks = (n+1)/SUPERBLOCK_SIZE + ((n+1)%SUPERBLOCK_SIZE != 0);
		n_blocks = (n+1)/128 + ((n+1)%128 != 0);
		nbytes = (n_blocks * BYTES_PER_BLOCK);//number of bytes effectively filled with data

		superblock_ranks = vector<p_rank>(n_superblocks);

		/*
		 * this block of code ensures that data is aligned by 64 bytes = 512 bits
		 */
		memory = vector<uint8_t>(nbytes+ALN,0);
		data = memory.data();
		while(uint64_t(data) % ALN != 0) data++;

		cout << "alignment of data: " << (void*)data << endl;

	}

	/*
	 * set the i-th character to c.
	 */
	void set(uint64_t i, uint8_t c){

		assert(i<n);

		//A is 0x0
		uint64_t b = 	(c == 'C')*  0x1 +
						(c == 'G')*  0x2 +
						(c == 'T')*  0x3 +
						(c == TERM)* 0x4 +
						(c == 'N')*  0x5;

		uint64_t superblock_number = i / SUPERBLOCK_SIZE;
		uint64_t superblock_off = i % SUPERBLOCK_SIZE;
		uint64_t block_number = superblock_off / BLOCK_SIZE;
		uint64_t block_off = superblock_off % BLOCK_SIZE;

		//chars[0..3] contains 1st, 2nd, 3rd bits of the 128 characters
		__uint128_t* chars = (__uint128_t*)(data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK);

		//from rightmost to leftmost bit
		chars[0] |= ((__uint128_t(b&0x1))<<(128-(block_off+1)));
		chars[1] |= ((__uint128_t((b&0x2)>>1))<<(128-(block_off+1)));
		chars[2] |= ((__uint128_t((b&0x4)>>2))<<(128-(block_off+1)));

		assert(operator[](i)==c);

	}

	//return i-th character
	uint8_t operator[](uint64_t i){

		assert(i<n);

		uint64_t superblock_number = i / SUPERBLOCK_SIZE;
		uint64_t superblock_off = i % SUPERBLOCK_SIZE;
		uint64_t block_number = superblock_off / BLOCK_SIZE;
		uint64_t block_off = superblock_off % BLOCK_SIZE;

		//chars[0..3] contains 1st, 2nd, 3rd least significant bits of the 128 characters
		__uint128_t* chars = (__uint128_t*)(data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK);

		uint64_t b =	((chars[0]>>(128-(block_off+1)))&0x1) +
						(((chars[1]>>(128-(block_off+1)))&0x1)<<1) +
						(((chars[2]>>(128-(block_off+1)))&0x1)<<2);

		return 	(b == 0)*'A' +
				(b == 1)*'C' +
				(b == 2)*'G' +
				(b == 3)*'T' +
				(b == 4)*TERM;

	}

	void build_rank_support(){

		p_rank superblock_r = {};
		p_rank block_r = {};

		for(uint64_t bl = 0; bl < n_blocks-1; ++bl){

			uint64_t superblock_number = bl/BLOCKS_PER_SUPERBLOCK;
			uint64_t block_number = bl%BLOCKS_PER_SUPERBLOCK;

			if(block_number == 0){

				superblock_ranks[superblock_number]=superblock_r;
				block_r = {};

			}

			set_counters(superblock_number, block_number,block_r);

			p_rank local_rank = block_rank(superblock_number, block_number);

			block_r = block_r + local_rank;
			superblock_r = superblock_r + local_rank;

		}

		uint64_t superblock_number = (n_blocks-1)/BLOCKS_PER_SUPERBLOCK;
		uint64_t block_number = (n_blocks-1)%BLOCKS_PER_SUPERBLOCK;

		if(block_number == 0){

			superblock_ranks[superblock_number]=superblock_r;
			block_r = {};

		}

		set_counters(superblock_number, block_number,block_r);

		assert(check_rank());

	}

	/*
	 * Parallel rank of (A,C,T,G) for the length-j prefix of the block of 512 bytes (128 characters) starting in start.
	 * The block belongs to the superblock_nr-th superblock
	 *
	 */
	p_rank parallel_rank(uint64_t i){

		uint64_t superblock_number = i / SUPERBLOCK_SIZE;
		uint64_t superblock_off = i % SUPERBLOCK_SIZE;
		uint64_t block_number = superblock_off / BLOCK_SIZE;
		uint64_t block_off = superblock_off % BLOCK_SIZE;

		p_rank superblock_r = superblock_ranks[superblock_number];
		p_rank block_r = get_counters(superblock_number,block_number);

		return superblock_r + block_r + block_rank(superblock_number, block_number, block_off);

	}

	/*
	 * standard rank. c can be A,C,G,T, or TERM
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

	uint64_t serialize(std::ostream& out){

		uint64_t w_bytes = 0;

		out.write((char*)&n,sizeof(n));
		out.write((char*)&nbytes,sizeof(nbytes));
		out.write((char*)&n_superblocks,sizeof(n_superblocks));
		out.write((char*)&n_blocks,sizeof(n_blocks));

		w_bytes += sizeof(n) + sizeof(nbytes) + sizeof(n_superblocks) + sizeof(n_blocks);

		out.write((char*)superblock_ranks.data(),n_superblocks*sizeof(p_rank));
		w_bytes += n_superblocks*sizeof(p_rank);

		out.write((char*)data,nbytes*sizeof(uint8_t));
		w_bytes += nbytes*sizeof(uint8_t);

		return w_bytes;

	}

	void load(std::istream& in) {

		in.read((char*)&n,sizeof(n));
		in.read((char*)&nbytes,sizeof(nbytes));
		in.read((char*)&n_superblocks,sizeof(n_superblocks));
		in.read((char*)&n_blocks,sizeof(n_blocks));

		superblock_ranks = vector<p_rank>(n_superblocks);
		in.read((char*)superblock_ranks.data(),n_superblocks*sizeof(p_rank));

		memory = vector<uint8_t>(nbytes+ALN,0);
		data = memory.data();
		while(uint64_t(data) % ALN != 0) data++;
		in.read((char*)data,nbytes*sizeof(uint8_t));

		assert(check_rank());

	}

	uint64_t size(){
		return n;
	}

private:

	/*
	 * rank in block given as coordinates: superblock, block, offset in block
	 */
	inline p_rank block_rank(uint64_t superblock_number, uint64_t block_number, uint64_t block_off){

		assert(block_off<128);

		//starting address of the block
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK;

		//chars[0..3] contains 1st, 2nd, 3rd bits of the 128 characters
		__uint128_t* chars = (__uint128_t*)(start);

		//partial ranks
		uint32_t * block_ranks = (uint32_t*)(start+48);

		__uint128_t PAD = ((~__uint128_t(0))>>block_off);

		__uint128_t b2 = ~(chars[2] | PAD); //most significant bit, padded and negated

		p_rank res = {

			popcount128(b2 & (~chars[1]) & (~chars[0])),
			popcount128(b2 & (~chars[1]) & (chars[0])),
			popcount128(b2 & (chars[1]) & (~chars[0])),
			popcount128(b2 & (chars[1]) & (chars[0]))

		};

		//assert(check_rank_local(res, superblock_number*SUPERBLOCK_SIZE + block_number*BLOCK_SIZE, block_off));

		return res;

	}

	/*
	 * rank in whole block given as coordinates: superblock, block
	 */
	inline p_rank block_rank(uint64_t superblock_number, uint64_t block_number){

		//starting address of the block
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK;

		//chars[0..3] contains 1st, 2nd, 3rd bits of the 128 characters
		__uint128_t* chars = (__uint128_t*)(start);

		//partial ranks
		uint32_t * block_ranks = (uint32_t*)(start+48);

		__uint128_t b2 = ~chars[2]; //most significant bit

		p_rank res = {

			popcount128(b2 & (~chars[1]) & (~chars[0])),
			popcount128(b2 & (~chars[1]) & (chars[0])),
			popcount128(b2 & (chars[1]) & (~chars[0])),
			popcount128(b2 & (chars[1]) & (chars[0]))

		};

		//assert(check_rank_local(res, superblock_number*SUPERBLOCK_SIZE + block_number*BLOCK_SIZE, block_off));

		return res;

	}

	bool check_rank_local(p_rank p, uint64_t i, uint64_t len){

		p_rank r = {};

		for(uint64_t j = i; j < i+len; ++j){

			if(operator[](j)=='A') r.A++;
			if(operator[](j)=='C') r.C++;
			if(operator[](j)=='G') r.G++;
			if(operator[](j)=='T') r.T++;

		}

		/*if(r != p){
			cout << "Error in local rank "
					<< p.A << "/" << r.A << " "
					<< p.C << "/" << r.C << " "
					<< p.G << "/" << r.G << " "
					<< p.T << "/" << r.T << endl;
		}*/

		return r == p;

	}

	bool check_rank(){

		p_rank p = {};

		bool res = true;

		for(int i=0;i<size();++i){

			auto r = parallel_rank(i);

			if(p != r){

				res = false;
				/*cout << "Error in local rank at position " << n << " "
				<< p.A << "/" << r.A << " "
				<< p.C << "/" << r.C << " "
				<< p.G << "/" << r.G << " "
				<< p.T << "/" << r.T << endl;*/

			}

			p.A += (operator[](i)=='A');
			p.C += (operator[](i)=='C');
			p.G += (operator[](i)=='G');
			p.T += (operator[](i)=='T');

		}

		auto r = parallel_rank(n);

		if(p != r){

			res = false;
			/*cout << "Error in local rank at position " << n << " "
			<< p.A << "/" << r.A << " "
			<< p.C << "/" << r.C << " "
			<< p.G << "/" << r.G << " "
			<< p.T << "/" << r.T << endl;*/

		}

		if(res){

			cout << "rank is correct" << endl;

		}else{

			cout << "rank is not correct" << endl;

		}

		return res;

	}

	/*
	 * check that the string contains exactly the same characters as the file in path
	 */
	bool check_content(string path){

		ifstream ifs(path);

		bool res = true;

		for(uint64_t i=0;i<n;++i){

			uint8_t c;
			ifs.read((char*)&c, sizeof(uint8_t));

			if(operator[](i) != c) res = false;

		}

		if(res){

			cout << "string content is valid" << endl;

		}else{

			cout << "string content is not valid" << endl;

		}

		return res;

	}

	/*
	 * set counters in the i-th block to r
	 */
	void set_counters(uint64_t superblock_number, uint64_t block_number, p_rank r){

		//block start
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK;
		uint32_t * block_ranks = (uint32_t*)(start+48);

		block_ranks[0] = r.A;
		block_ranks[1] = r.C;
		block_ranks[2] = r.G;
		block_ranks[3] = r.T;

		assert(get_counters(superblock_number,block_number) == r);

	}

	/*
	 * get counters of the i-th block
	 */
	inline p_rank get_counters(uint64_t superblock_number, uint64_t superblock_off){

		//block start
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + superblock_off*BYTES_PER_BLOCK;
		uint32_t * block_ranks = (uint32_t*)(start+48);

		return {
			block_ranks[0],
			block_ranks[1],
			block_ranks[2],
			block_ranks[3]
		};

	}

	uint64_t n_superblocks = 0;
	uint64_t n_blocks = 0;

	vector<uint8_t> memory; //allocated memory

	//data aligned with blocks of 64 bytes = 512 bits
	uint8_t * data = NULL;

	vector<p_rank> superblock_ranks;

	uint64_t nbytes = 0; //bytes used in data
	uint64_t n = 0;

};


#endif /* INTERNAL_DNA_STRING_HPP_ */
