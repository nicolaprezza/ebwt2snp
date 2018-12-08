/*
 * build_index.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: nico
 */


#include <iostream>
#include <fstream>
#include "internal/dna_bwt.hpp"
#include "internal/bwt.hpp"
#include "internal/dna_string.hpp"
#include "include.hpp"

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
 	"Note: string terminators must be stored as '#'." << endl;
	exit(0);
}

struct block{
	uint64_t x_1;
	uint64_t x_2;
	uint64_t x_3;
	uint64_t x_4;
	uint64_t x_5;
	uint64_t x_6;
	uint64_t x_7;
	uint64_t x_8;
};

int bin(uint64_t x){

	bool bit=1;

	int cnt=0;

	do{

		bit = (x&uint64_t(1));
		x = x>>1;

		if(not bit) cnt++;

	}while(bit==0);

	return cnt;
}

void printbin(uint64_t x){

	bool bit=1;

	int cnt=0;

	while(x>0){

		bit = (x&uint64_t(1));
		x = x>>1;

		cout << uint64_t(bit);

	};

	cout << endl;
}

int main(int argc, char** argv){


/*	srand(time(NULL));

	int n = 100000000;

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

	/*cout << endl;
	for(int i=0;i<n;++i) cout << dna[i];
	cout << endl;
	for(int i=0;i<n;++i) cout << s[i];
	cout << endl;*/

	/*dna.build_rank_support();






	for(int i=0;i<n;++i){

		if(dna[i] != s[i]){
			cout << "err: " << i << " " << dna[i] << " " << s[i] << endl;
			exit(0);
		}

	}

	p_rank r = {};

	for(int i=0;i<n;++i){

		auto x = dna.parallel_rank(i);

		if(r != x){
			cout << "err in rank at position " << i << endl;
			cout << r.A << " " << x.A << endl;
			cout << r.C << " " << x.C << endl;
			cout << r.G << " " << x.G << endl;
			cout << r.T << " " << x.T << endl;

			exit(0);
		}

		r.A += (s[i]=='A');
		r.C += (s[i]=='C');
		r.G += (s[i]=='G');
		r.T += (s[i]=='T');

	}

	cout << "success" << endl;exit(0);





*/





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

	cout << "Building DNA-optimized BWT ..." << endl;
	auto BWT = dna_bwt_t(input_file);

	cout << "Done. Storing to file ..." << endl;
	BWT.save_to_file(output_file);
	cout << "Done. " << endl;

	/*if(rle){

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

	}*/


}

