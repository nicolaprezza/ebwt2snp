/*
 * merge_bwt.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: nico
 */


#include <iostream>
#include <fstream>
#include "internal/bwt.hpp"
#include "internal/bwt_merger.hpp"

using namespace std;

string input_bwt1;
string input_bwt2;
string output_file;

bool lcp = false;

void help(){

	cout << "merge_bwt [options]" << endl <<
	"Options:" << endl <<
	"-h          Print this help" << endl <<
	"-1 <arg>    Input BWT index 1 (REQUIRED)" << endl <<
	"-2 <arg>    Input BWT index 2 (REQUIRED)" << endl <<
	"-o <arg>    Output prefix (REQUIRED)" << endl <<
	"-l          Compute LCP of merged BWT" << endl;
	exit(0);
}

template<class bwt_1_t, class bwt_2_t, typename lcp_t>
void merge(ifstream &in1, ifstream &in2){

	bwt_1_t BWT1;
	bwt_2_t BWT2;

	BWT1.load(in1);
	BWT2.load(in2);

	bwt_merger<bwt_1_t,bwt_2_t, uint8_t> M(&BWT1, &BWT2, lcp);

	M.save_to_file(output_file);

}

int main(int argc, char** argv){

	if(argc < 4) help();

	int opt;
	while ((opt = getopt(argc, argv, "h1:2:o:l")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case '1':
				input_bwt1 = string(optarg);
			break;
			case '2':
				input_bwt2 = string(optarg);
			break;
			case 'o':
				output_file = string(optarg);
			break;
			case 'l':
				lcp = true;
			break;
			default:
				help();
			return -1;
		}
	}

	if(input_bwt1.size()==0) help();
	if(input_bwt2.size()==0) help();
	if(output_file.size()==0) help();

	cout << "Input bwt index file 1: " << input_bwt1 << endl;
	cout << "Input bwt index file 2: " << input_bwt2 << endl;
	cout << "Output prefix: " << output_file << endl;

	std::ifstream in1(input_bwt1);
	std::ifstream in2(input_bwt2);

	uint8_t type1 = 0;
	uint8_t type2 = 0;

	in1.read((char*)&type1,sizeof(type1));
	in2.read((char*)&type2,sizeof(type2));

	switch(type1){

	case 's': switch(type2){

		case 's': merge<suc_bwt,suc_bwt,uint8_t>(in1,in2); break;
		case 'e': merge<suc_bwt,rrr_bwt,uint8_t>(in1,in2); break;
		case 'r': merge<suc_bwt,rle_bwt,uint8_t>(in1,in2); break;
		default: break;

	}; break;

	case 'e': switch(type2){

		case 's': merge<rrr_bwt,suc_bwt,uint8_t>(in1,in2); break;
		case 'e': merge<rrr_bwt,rrr_bwt,uint8_t>(in1,in2); break;
		case 'r': merge<rrr_bwt,rle_bwt,uint8_t>(in1,in2); break;
		default: break;

	}; break;

	case 'r':switch(type2){

		case 's': merge<rle_bwt,suc_bwt,uint8_t>(in1,in2); break;
		case 'e': merge<rle_bwt,rrr_bwt,uint8_t>(in1,in2); break;
		case 'r': merge<rle_bwt,rle_bwt,uint8_t>(in1,in2); break;
		default: break;

	}; break;

	default:break;

	}


}

