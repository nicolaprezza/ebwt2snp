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
	string output_index = input_bwt;
	input_bwt = input_bwt.substr(0,input_bwt.rfind(".bwt"));
	input_bwt.append(".rlbwt");

	ifstream is(input_bwt);

	std::string bwt((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());

	rle_string_sd rlbwt(bwt);

	rlbwt.save_to_file(output_index);

}

