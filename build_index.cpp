/*
 * build_index.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: nico
 */


#include <iostream>
#include <fstream>
#include "internal/bwt.hpp"
#include "internal/bwt_merger.hpp"

using namespace std;

string input_file;
string output_file;
bool convert_n = true;
bool rle = false;
bool suc = true;
bool rrr = false;

void help(){

	cout << "build_index [options]" << endl <<
	"Options:" << endl <<
	"-h          Print this help" << endl <<
	"-i <arg>    Input BWT file (REQUIRED)" << endl <<
	"-o <arg>    Output index file (REQUIRED)" << endl <<
	//"-n          Do not convert N's to random bases (DEFAULT: N's are converted)" << endl <<
	"-r          Run-length compressed BWT" << endl <<
	"-s          Huffman + succinct bitvectors (DEFAULT)" << endl <<
	"-e          Huffman + entropy-compressed bitvectors (RRR)" << endl << endl <<
 	"Note: string terminators must be stored as '$'." << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 3) help();

	int opt;
	while ((opt = getopt(argc, argv, "hi:o:rse")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'i':
				input_file = string(optarg);
			break;
			case 'o':
				output_file = string(optarg);
			break;
			case 'n':
				convert_n = false;
			break;
			case 'r':
				rle=true;suc=false;rrr=false;
			break;
			case 's':
				rle=false;suc=true;rrr=false;
			break;
			case 'e':
				rle=false;suc=false;rrr=true;
			break;
			default:
				help();
			return -1;
		}
	}

	if(input_file.size()==0) help();
	if(output_file.size()==0) help();

	srand(time(NULL));

	cout << "Input bwt file: " << input_file << endl;
	cout << "Output index file: " << output_file << endl;

	if(rle){

		cout << "Building run-length compressed BWT" << endl;
		auto BWT = rle_bwt(input_file,convert_n);

		std::ofstream out(output_file);
		uint8_t type = 'r';
		out.write((char*)&type,sizeof(type));
		BWT.serialize(out);
		out.close();

	}

	if(suc){

		cout << "Building succinct+Huffman compressed BWT" << endl;
		auto BWT = suc_bwt(input_file,convert_n);

		std::ofstream out(output_file);
		uint8_t type = 's';
		out.write((char*)&type,sizeof(type));
		BWT.serialize(out);
		out.close();

	}

	if(rrr){

		cout << "Building RRR+Huffman compressed BWT" << endl;
		auto BWT = rrr_bwt(input_file,convert_n);

		std::ofstream out(output_file);
		uint8_t type = 'e';
		out.write((char*)&type,sizeof(type));
		BWT.serialize(out);
		out.close();

	}


}

