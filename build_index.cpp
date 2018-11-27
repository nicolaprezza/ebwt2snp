/*
 * build_index.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: nico
 */


#include <iostream>
#include <fstream>
#include "internal/rle_string.hpp"

using namespace std;

void help(){

	cout << "build_index input.bwt" << endl << endl <<
	"Builds a run-length string with rank/select support on the input BWT." << endl <<
	"Assumption: string terminators are all stored as '$'. Output is stored in input.rlbwt" << endl;

	exit(0);
}


int main(int argc, char** argv){

	if(argc!=2) help();

	string input_bwt = argv[1];
	string output_index = input_bwt.substr(0,input_bwt.rfind(".bwt"));
	output_index.append(".rlbwt");

	cout << "Input bwt file: " << input_bwt << endl;
	cout << "Output index file: " << output_index << endl;

	ifstream is(input_bwt);

	cout << "loading input ... " << flush;
	std::string bwt((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
	cout << "done." << endl;

	cout << "Indexing BWT ... " << flush;
	rle_string_sd rlbwt(bwt);
	cout << "done." << endl;

	cout << "Storing index to file ... " << flush;
	rlbwt.save_to_file(output_index);
	cout << "done." << endl;

}

