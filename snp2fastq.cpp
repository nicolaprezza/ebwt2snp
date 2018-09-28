// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <unistd.h>
#include <sstream>
#include <set>
#include <cstring>

using namespace std;


void help(){

	cout << "snp2fastq calls.snp" << endl << endl <<
	"Converts clust2snp's calls 'calls.snp' into a fastq  file 'calls.snp.fastq'. The  output contains  one " << endl <<
	"read per call, where we put the second individual's DNA in the read's name, and the first individual's " << endl <<
	"DNA in the read DNA. Base qualities are fake (all maximum)." << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc != 2) help();

	string infile = argv[1];

	string outfile = infile;
	outfile.append(".fastq");

	ifstream is(infile);
	ofstream of(outfile);

	string str;
	unsigned int idx=0;

	string header;
	string dna;

	while(getline(is, str)){

		if(idx%4==0){//first line of call

			str = str.substr(1);//remove '>'

			char *p = strtok((char*)str.c_str(), "_");

			header = string(p);

			p = strtok(NULL,"_");
			p = strtok(NULL,"_");
			p = strtok(NULL,"_");//this contains event number

			//pick number before '|'
			string ev(p);
			p = strtok((char*)ev.c_str(), "|");
			ev = string(p);

			header.append("_");
			header.append(ev);
			header.append("_");

		}

		if(idx%4==1){//DNA of first individual

			dna = str;

		}

		//if(idx%4==2) ignore: header of second individual

		if(idx%4==3){//DNA of second individual

			header.append(str);

			//now output fastq entry
			of 	<< "@" << header << endl <<
				dna << endl <<
				"+" << endl <<
				string(dna.length(),'I') << endl;

		}

		idx++;

	}

	is.close();
	of.close();

}
