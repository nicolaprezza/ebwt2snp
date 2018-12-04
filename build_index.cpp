/*
 * build_index.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: nico
 */


#include <iostream>
#include <fstream>
#include "internal/rlbwt.hpp"
#include "internal/wt_bwt.hpp"

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

	cout << "Indexing BWT ... " << flush;
	wt_bwt BWT(input_bwt);
	cout << "done." << endl;

	cout << "The BWT has " << BWT.size() << " characters." << endl;

	cout << "Storing index to file ... " << flush;
	BWT.save_to_file(output_index);
	cout << "done." << endl;

}

