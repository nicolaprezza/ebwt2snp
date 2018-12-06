/*
 * build_index.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: nico
 */


#include <iostream>
#include <fstream>
#include "internal/bwt.hpp"
#include "internal/dna_string.hpp"

using namespace std;

string input_file;
string output_file;
bool rle = false;
bool suc = true;
bool rrr = false;

void help(){

	cout << "build_index [options]" << endl <<
	"Options:" << endl <<
	"-h          Print this help" << endl <<
	"-i <arg>    Input BWT file (REQUIRED)" << endl <<
	"-o <arg>    Output index file (REQUIRED)" << endl <<
	"-r          Run-length compressed BWT" << endl <<
	"-s          Huffman + succinct bitvectors (DEFAULT)" << endl <<
	"-e          Huffman + entropy-compressed bitvectors (RRR)" << endl << endl <<
 	"Note: string terminators must be stored as '$'." << endl;
	exit(0);
}

int main(int argc, char** argv){

	/*srand(time(NULL));

	int n = 100000;

	string s(n,'A');
	dna_string dna(n);

	for(int i=0;i<n;++i){

		int x = rand()%6;

		switch(x){

		case 0: s[i]='A'; dna.set(i,'A');break;
		case 1: s[i]='C'; dna.set(i,'C');break;
		case 2: s[i]='G'; dna.set(i,'G');break;
		case 3: s[i]='T'; dna.set(i,'T');break;
		case 4: s[i]=TERM; dna.set(i,TERM);break;
		case 5: s[i]='N'; dna.set(i,'N');break;

		}

	}

	cout << endl;
	for(int i=0;i<n;++i) cout << dna[i];
	cout << endl;
	for(int i=0;i<n;++i) cout << s[i];
	cout << endl;


	for(int i=0;i<n;++i){

		if(dna[i] != s[i]){
			cout << "err: " << i << " " << dna[i] << " " << s[i] << endl;
			exit(0);
		}

	}

	cout << "success" << endl;exit(0);*/

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
		auto BWT = rle_bwt(input_file);

		std::ofstream out(output_file);
		uint8_t type = 'r';
		out.write((char*)&type,sizeof(type));
		BWT.serialize(out);
		out.close();

	}

	if(suc){

		cout << "Building succinct+Huffman compressed BWT" << endl;
		auto BWT = suc_bwt(input_file);

		std::ofstream out(output_file);
		uint8_t type = 's';
		out.write((char*)&type,sizeof(type));
		BWT.serialize(out);
		out.close();

	}

	if(rrr){

		cout << "Building RRR+Huffman compressed BWT" << endl;
		auto BWT = rrr_bwt(input_file);

		std::ofstream out(output_file);
		uint8_t type = 'e';
		out.write((char*)&type,sizeof(type));
		BWT.serialize(out);
		out.close();

	}


}

