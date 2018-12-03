/*
 * rlbwt.hpp
 *
 *  Created on: Dec 3, 2018
 *      Author: nico
 *
 *  run-length BWT for large files.
 *
 */

#include "include.hpp"
#include "rle_string.hpp"

#ifndef INTERNAL_RLBWT_HPP_
#define INTERNAL_RLBWT_HPP_

class rlbwt{

public:

	/*
	 * constructor: the path of a BWT file and number of blocks n_blocks.
	 *
	 * File will be read in chunks of block_size. We build one rle_string per block.
	 *
	 * Default block size: 2^29, i.e. approx. 537 MB
	 *
	 */
	rlbwt(string path, uint64_t block_size = (uint64_t(1)<<29) ){

		bsize = block_size;
		n = uint64_t(filesize(path));

		cout << "file len = " << n << endl;

		nblocks = (n/block_size) + ((n%block_size)!=0);

		cout << "nblocks = " << nblocks << endl;

		partial_rank = vector<vector<uint64_t> >(5, vector<uint64_t>(nblocks,0));
		blocks = vector<rle_string<> >(nblocks);

		cout << "1" << endl;

		ifstream input(path);

		uint64_t block_idx = 0;

		//more efficient for small strings
		if(nblocks==1) block_size = n;

		cout << "2" << endl;

		for(uint64_t bl = 0;bl<nblocks;++bl) {

			string buffer (bl == nblocks-1 ? n%block_size : block_size,0); //init buffer

			input.read((char*)buffer.data(), buffer.size());

			cout << "3" << endl;

			blocks[block_idx++] = rle_string<>(buffer);

			cout << "5" << endl;

		}

		cout << "6" << endl;

		for(uint8_t c = 0;c<5;++c){

			for(uint64_t i=1;i<nblocks;++i){

				partial_rank[c][i] = partial_rank[c][i-1] + blocks[i-1].rank(blocks[i-1].size(), int_to_base(c));

			}

		}

		cout << "7" << endl;


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

		w_bytes += sizeof(n) + sizeof(bsize) + sizeof(nblocks);

		if(n==0) return w_bytes;

		for(uint64_t i=0;i<nblocks;++i){

			w_bytes += blocks[i].serialize(out);

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

		if(n==0) return;

		blocks = vector<rle_string<> >(nblocks);
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

	uint64_t n = 0;//BWT length
	vector<rle_string<> > blocks;

	//partial rank info for each block and character
	vector<vector<uint64_t> > partial_rank;

	uint64_t bsize = 0;
	uint64_t nblocks = 0;

};

#endif /* INTERNAL_RLBWT_HPP_ */
