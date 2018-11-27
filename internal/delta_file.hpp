/*
 * delta_file.hpp
 *
 *  Created on: Nov 27, 2018
 *      Author: nico
 */

#include <iostream>
#include <fstream>
#include <sdsl/int_vector.hpp>
#include <sdsl/enc_vector.hpp>

#ifndef INTERNAL_DELTA_FILE_HPP_
#define INTERNAL_DELTA_FILE_HPP_

using namespace std;
using namespace sdsl;

class delta_file{

public:

	/*
	 * init empty delta file. Specify if the file is open for read/write operations
	 */
	delta_file(string path, bool write){

		this->write = write;

		if(write){

			out = new std::ofstream(path);
			buf_w = vector<uint64_t>();

		}else{

			in = new std::ifstream(path);
			buf_r.load(*in);
			i = 0;

		}

	}

	~delta_file(){
		if(write){

			if(buf_w.size()>0) flush();

			assert(out != NULL);
			delete out;

		}else{

			assert(in!=NULL);
			delete in;

		}
	}

	void push_back(uint64_t x){

		assert(write);

		buf_w.push_back(x);
		if(buf_w.size() == MAX_BUF_LEN) flush();

	}

	uint64_t get(){

		assert(not write);
		assert(i<buf_r.size());

		uint64_t x = buf_r[i++];

		if(i>=buf_r.size() and (not eof())){

			buf_r.load(*in);
			i=0;

		}


		return 0;

	}

	bool eof(){

		return in->eof();

	}

private:

	/*
	 * write to file buffer content
	 */
	void flush(){

		enc_vector<> deltas(buf_w);
		deltas.serialize(*out);

		buf_w = vector<uint64_t>();

	}

	bool write = true;

	istream * in = NULL;
	ostream * out = NULL;

	const int MAX_BUF_LEN = 10000;
	vector<uint64_t> buf_w;

	enc_vector<> buf_r;
	//position in buf_r during read
	int i = 0;

};


#endif /* INTERNAL_DELTA_FILE_HPP_ */
