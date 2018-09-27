// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include "include.hpp"
#include <unistd.h>

using namespace std;

/*
 * parameters
 *
 */
int min_def = 2; //Discard clusters smaller than this value
int K_def = 16; //require an LCP of at least k inside clusters
int k = 0;
string input;

int min_len=0;

//we require that inside a cluster, LCP[i+1] >= LCP[i] + delta
bool bcr = false;

void help(){

	cout << "eBWTclust [options]" << endl <<
	"Options:" << endl <<
	"-h          Print this help" << endl <<
	"-i <arg>    Input fasta file (REQUIRED)" << endl <<
	"-k <arg>    Minimum LCP required in clusters (default: " << K_def << ")" << endl <<
	"-m <arg>    Discard clusters smaller than this value (default: " << min_def << ")" << endl <<
//	"-b          Input index files are in BCR format (default: GESA format)" << endl << endl <<

	"\nTo run eBWTclust, you must first build the Enhanced Generalized Suffix Array of the input" << endl <<
	"sequences. The EGSA must be stored in the input file's folder adding extension .egsa to " << endl <<
	"the name of the input file  (github.com/felipelouza/egsa). Output  is stored in " << endl <<
	"reads.fasta.clusters." << endl;
	 exit(0);
}

void append_entry(ofstream & out, uint64_t start, uint16_t length){

	if(length >=min_len){

		out.write((char*)&start, sizeof(uint64_t));
		out.write((char*)&length, sizeof(uint16_t));

	}

}

/*
 * clusters = regions between local LCP minima (excluding tails where LCP < k)
 */
void cluster_lm(ifstream & egsa,ofstream & out){

	uint64_t null = ~uint64_t(0);

	uint64_t start = null;		//start position of cluster
	uint64_t i = 0;			//current BWT position
	uint64_t off = 0;		//current offset in cluster

	uint64_t prev_lcp = null;//previous LCP value

	//find local minima in the LCP
	t_GSA e1 = read_el(egsa, bcr);
	t_GSA e2 = read_el(egsa, bcr);
	t_GSA e3;

	start = e1.lcp >= k ? 0 :
			e2.lcp >= k ? 1 : null;

	i = 1;//index of e2

	while(not egsa.eof()){

		//read next value
		e3 = read_el(egsa, bcr);

		//cout << e3.text << " " << e3.suff << " " << e3.lcp << " " << char(e3.bwt) << endl;

		//e2 is start of a (possibly flat) local minima
		if(		start != null and
					(		(e1.lcp > e2.lcp and e2.lcp <= e3.lcp) or
							e3.lcp < k
					)
					){

			uint16_t length = (i - start) + 1; //this cluster ends in e2
			append_entry(out, start, length);

			start = null;

		}

		e1 = e2;
		e2 = e3;
		++i;

		if(start==null and e2.lcp >= k){

			start = i;

		}

	}

	/*
	 * check if last cluster was closed. If not, output it
	 */
	if(start != null){

		uint16_t length = (i - start) + 1; //this cluster ends in e2
		append_entry(out, start, length);

		start = null;

	}

}

int main(int argc, char** argv){

	if(argc < 2) help();

	int opt;
	while ((opt = getopt(argc, argv, "hk:i:m:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'b':
				bcr = true;
			break;
			case 'k':
				k = atoi(optarg);
				//cout << "k = " << k << "\n";
			break;
			case 'm':
				min_len = atoi(optarg);
				//cout << "k = " << k << "\n";
			break;
			case 'i':
				input = string(optarg);
				//cout << "input = " << input << "\n";
			break;
			default:
				help();
			return -1;
		}
	}

	k = k==0?K_def:k;
	min_len = min_len==0?min_def:min_len;

	if(input.compare("")==0) help();

	string egsa_path = input;
	egsa_path.append(".gesa");

	cout << "This is eBWTclust. Input index file: " << egsa_path << endl;

	ifstream egsa;

	egsa.open(egsa_path, ios::in | ios::binary);

	if(not egsa.is_open()){
		cout << "Error: missing file " <<egsa_path << endl << endl;
		help();
	}

	string filename_out = input;
	filename_out.append(".clusters");
	ofstream out;
	out.open(filename_out, ios::out | ios::binary);

	cluster_lm(egsa,out);

	egsa.close();
	out.close();

}
