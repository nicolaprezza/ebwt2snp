/*
 * bwt.hpp
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

#ifndef INTERNAL_BWT_HPP_
#define INTERNAL_BWT_HPP_

using namespace sdsl;

template<class str_type = rrr_str>
class bwt{

public:

	/*
	 * constructor: the path of a BWT file and number of blocks n_blocks.
	 *
	 * File will be read in chunks of block_size. We build one rle_string per block.
	 *
	 * Default block size: 2^29, i.e. approx. 537 MB
	 *
	 */
	bwt(string path, bool random_n = true, uint64_t block_size = (uint64_t(1)<<29) ){

		//build F column
		F = vector<uint64_t>(256,0);

		bsize = block_size;
		n = uint64_t(filesize(path));

		nblocks = (n/block_size) + ((n%block_size)!=0);

		partial_rank = vector<vector<uint64_t> >(5, vector<uint64_t>(nblocks,0));
		blocks = vector<str_type>(nblocks);

		ifstream input(path);

		uint64_t block_idx = 0;

		for(uint64_t bl = 0;bl<nblocks;++bl) {

			string buffer (bl == nblocks-1 ? n%block_size : block_size,0); //init buffer

			input.read((char*)buffer.data(), buffer.size());

			convert(buffer, random_n);
			count_chars(buffer);

			blocks[block_idx++] = str_type(buffer);

		}

		//build F column
		for(uint64_t i=255;i>0;--i)
			F[i] = F[i-1];

		F[0] = 0;

		for(uint64_t i=1;i<256;++i)
			F[i] += F[i-1];

		for(uint8_t c = 0;c<5;++c){

			for(uint64_t i=1;i<nblocks;++i){

				partial_rank[c][i] = partial_rank[c][i-1] + blocks[i-1].rank(blocks[i-1].size(), int_to_base(c));

			}

		}

	}

	/*
	 * get full BWT range
	 */
	range_t full_range(){

		//inclusive range
		return {0,size()-1};

	}

	/*
	 * \param r inclusive range of a string w
	 * \param c character
	 * \return inclusive range of cw
	 */
	range_t LF(range_t rn, uint8_t c){

		//if character does not appear in the text, return empty pair
		if((c==255 and F[c]==size()) || F[c]>=F[c+1])
			return {1,0};

		//number of c before the interval
		uint64_t c_before = rank(rn.first,c);

		//number of c inside the interval rn
		uint64_t c_inside = rank(rn.second+1,c) - c_before;

		//if there are no c in the interval, return empty range
		if(c_inside==0) return {1,0};

		uint64_t l = F[c] + c_before;

		return {l,l+c_inside-1};

	}

	//backward navigation of the BWT
	uint64_t LF(uint64_t  i){

		auto c = operator[](i);
		return F[c] + rank(i,c);

	}

	//forward navigation of the BWT
	uint64_t FL(uint64_t  i){

		//i-th character in first BWT column
		auto c = F_at(i);

		//this c is the j-th (counting from 0)
		uint64_t j = i - F[c];

		return select(j,uint8_t(c));

	}

	//forward navigation of the BWT, where for efficiency we give c=F[i] as input
	uint64_t FL(uint64_t  i, uint8_t c){

		//i-th character in first BWT column
		assert(c == F_at(i));

		//this c is the j-th (counting from 0)
		uint64_t j = i - F[c];

		return select(j,uint8_t(c));

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

		//if character does not appear in the text, return empty pair
		if((c==255 and F[c]==size()) || F[c]>=F[c+1])
			return {1,0};

		uint64_t l = F[c];
		uint64_t r = size()-1;

		if(c<255)
			r = F[c+1]-1;

		return {l,r};

	}

	/*
	 * Return BWT range of pattern P
	 */
	range_t count(string &P){

		auto range = full_range();
		uint64_t m = P.size();

		for(uint64_t i=0;i<m and range.second>=range.first;++i)
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

		return blocks[i/bsize][i%bsize];

	}

	/*
	 * number of c before position i excluded
	 */
	uint64_t rank(uint64_t i, uint8_t c){

		assert(i<=n);

		uint64_t b = i/bsize;
		uint64_t off = i%bsize;

		return partial_rank[base_to_int(c)][b] + blocks[b].rank(off,c);

	}

	uint64_t size(){
		return n;
	}

	/*
	 * position of i-th character c. i starts from 0!
	 */
	uint64_t select(uint64_t i, uint8_t c){

		//block containing i-th c
		uint64_t idx = std::distance( 	partial_rank[base_to_int(c)].begin(),
										std::upper_bound(	partial_rank[base_to_int(c)].begin(),
															partial_rank[base_to_int(c)].end(),
															i)
		) - 1;

		return blocks[idx].select( i - partial_rank[base_to_int(c)][idx], c ) + idx*bsize;

	}

	uint64_t serialize(std::ostream& out){

		uint64_t w_bytes = 0;

		out.write((char*)&n,sizeof(n));
		out.write((char*)&bsize,sizeof(bsize));
		out.write((char*)&nblocks,sizeof(nblocks));
		out.write((char*)&F,256*sizeof(uint64_t));

		w_bytes += sizeof(n) + sizeof(bsize) + sizeof(nblocks);

		if(n==0) return w_bytes;

		for(uint64_t i=0;i<nblocks;++i){

			uint64_t x = blocks[i].serialize(out);

			w_bytes += x;

		}

		for(uint8_t c=0;c<5;++c){
			for(uint64_t i=0;i<nblocks;++i){

				out.write((char*)partial_rank[c].data(),sizeof(uint64_t)*nblocks);
				w_bytes += sizeof(uint64_t)*nblocks;

			}
		}

		return w_bytes;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&n,sizeof(n));
		in.read((char*)&bsize,sizeof(bsize));
		in.read((char*)&nblocks,sizeof(nblocks));

		F = vector<uint64_t>(256, 0);
		in.read((char*)&F,256*sizeof(uint64_t));

		if(n==0) return;

		blocks = vector<str_type>(nblocks);
		partial_rank = vector<vector<uint64_t> >(5, vector<uint64_t>(nblocks,0));

		for(uint64_t i=0;i<nblocks;++i){

			blocks[i].load(in);

		}

		for(uint8_t c=0;c<5;++c){
			for(uint64_t i=0;i<nblocks;++i){

				in.read((char*)partial_rank[c].data(),sizeof(uint64_t)*nblocks);

			}
		}

	}

	void save_to_file(string path){

		std::ofstream out(path);
		serialize(out);
		out.close();

	}

	void load_from_file(string path){

		std::ifstream in(path);
		load(in);
		in.close();

	}

private:

	void count_chars(string & s){

		for(auto c:s) F[c]++;

	}

	void convert(string & s, bool random_n){

		for(uint64_t i=0;i<s.size();++i)
			s[i] = (not random_n) and (s[i]=='N' or s[i]=='n') ? 'N' : int_to_base(base_to_int(s[i]));

	}

	uint64_t n = 0;//BWT length
	vector<str_type> blocks;

	//partial rank info for each block and character
	vector<vector<uint64_t> > partial_rank;

	uint64_t bsize = 0;
	uint64_t nblocks = 0;

	vector<uint64_t> F;

};

typedef bwt<rle_str> rle_bwt;
typedef bwt<suc_str> suc_bwt;
typedef bwt<rrr_str> rrr_bwt;

#endif /* INTERNAL_BWT_HPP_ */
