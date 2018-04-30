// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <algorithm>
#include "include.hpp"
#include <unistd.h>
#include <math.h>
#include <iomanip>
#include "read64_t.hpp"

using namespace std;

//int mcov_def = 5;//minimum coverage
int MAX_MEM_def = 2048;
double err_allowed_def = 0.15;
int k_left_def = 20;//extract k_left nucleotides from left of suffix array range, for each entry in the cluster
int k_right_def = 30;//extract k_right nucleotides from right of suffix array range, only for entry with max LCP

//discard clusters whose sizes fall in this area of the tails of Poisson distribution
double pval_def = 0.05;
int max_snvs_def = 3;//maximum number of SNVs allowed in left contexts
int mcov_out_def = 4;//minimum coverage required in the output events

int mcov_out = 0;//minimum coverage required in the output events
int max_snvs = 0;//maximum number of SNVs allowed in left contexts
string input;
uint64_t nr_reads1 = 0;
int mcov1 = 0, mcov2 = 0;//minimum coverage required for each allele; computed automatically
int max_cov1;//maximum coverage allowed; computed automatically
int max_cov2;//maximum coverage allowed; computed automatically
double pval = 0;
double err_allowed = 0;//
double lcpfreq = 0;//output read's length is the max LCP in common of this fraction of contiguous suffixes in the cluster
bool bcr = false;
bool diploid = false;
uint64_t MAX_MEM = 0;//max RAM allowed for events buffer, in MB
double read_len_def = 100;//average read length

//letter with at least this frequency goes in consensus
double min_consensus_freq = 0.65;

uint64_t id_nr = 1;

int k_left = 0;//extract k_left nucleotides from left of suffix array range, for each entry in the cluster
int k_right = 0;//extract k_right nucleotides from right of suffix array range, for each entry in the cluster

double read_len = 0;//average read length
double min_k = 16;//minimum LCP we look at
double err_rate = 0.01;//sequencing error rate
//fraction of average coverage that we can expect to be in the clusters (computed in main)
double frac_cov = 0;

bool discoSNP=true;

/*
 * the two Poisson probability and cumulative functions
 */
vector<double> Poi1;
vector<double> Poi2;
vector<double> Cumul1;
vector<double> Cumul2;

void help(){

	cout << "clust2snp [options]" << endl <<
	"Options:" << endl <<
	"-h          Print this help." << endl <<
	"-i <arg>    Input fasta file (REQUIRED)." << endl <<
	"-n <arg>    Number of reads in the first sample (REQUIRED)." << endl <<
	"-L <arg>    Length of left-context, SNP included (default: " << k_left_def << ")." << endl <<
	"-R <arg>    Maximum length of right context, SNP excluded (default: " << k_right_def << ")." << endl <<
	"-v <arg>    Maximum number of SNPs allowed in a non-isolated group (default: " << max_snvs_def << ")."<< endl <<
	"-m <arg>    Minimum coverage per sample per event (default: " << mcov_out_def << ")."<< endl <<
	"-p <arg>    p-value of discarded cluster sizes (default: " << pval_def << ")."<< endl <<
	"-e <arg>    Tolerated error rate in clusters (default: " << err_allowed_def << ", range: [0, 0.5))." << endl <<
//	"-D          Samples are diploid (experimental! default: haploid)." << endl <<
//	"-b          Input index files are in BCR format (default: GESA format)." << endl << endl <<
//	"-d          Output fasta in discoSNP format (default: off)" << endl << endl<<

	"To run clust2snp, you must first build (1)  the Enhanced Generalized Suffix Array of the input" << endl <<
	"sequences, stored in a file with extension .0.egsa and with the same name of the input file" << endl <<
	"(github.com/felipelouza/egsa), and  (2) the cluster file  built with cluster-bwt. Output is" << endl <<
	"stored in reads.snp.fasta, where reads.fasta is the input fasta file." << endl << endl <<

	"Output: SNPs are output in KisSNP format as a fasta file. " << endl;

	exit(0);
}

/*
 * set of reads containing an interesting event
 *
 * read_prefixes indicates the lengths of the reads'
 * prefixes to be aligned
 *
 */
struct event{

	vector<uint64_t> read_ranks;
	vector<uint16_t> read_prefixes;
	vector<bool> indiv;
	uint64_t good_lcp_t;
	uint16_t good_lcp_s;
	uint16_t good_lcp_len;
	string hapl;
	string snp;

};

/*
 * looks at the frequencies in counts and decides which of these events is more likely:
 *
 * - repetitive region
 * - non-repetitive region, heterozigote/homozigote alleles
 *
 * returns a string with the 4 alleles. E.g.
 *
 * - AACC means that region is not repetitive and Indiv 1 and 2 have
 *   homozigote alleles A/C
 * - AAAC means that region is not repetitive, Indiv 1 has hom. allele
 *   A and indiv. 2 has het. alleles A/C
 * - NNNN means repetitive region in either indiv1/2 (cluster too big)
 *   or cluster not sufficiently covered in one of the 2 indiv, or
 *   allele sets are equal (e.g. ACAC)
 */
string event_type(vector<vector<unsigned int> > counts){

	string out = "NNNN";

	for(int indiv=0;indiv<2;++indiv){

		int mcov = indiv==0?mcov1:mcov2;
		int max_cov = indiv==0?max_cov1:max_cov2;

		int cluster_size=0;

		for(auto x:counts[indiv]) cluster_size+=x;

		/*
		 * first filter: too covered or too less covered event
		 */
		if(	cluster_size < mcov or cluster_size > max_cov){

			return "NNNN";

		}

		/*
		 * second filter: sum of two most frequent nucleotides must
		 * not be too covered or too less covered
		 */

		int max1=std::distance(counts[indiv].begin(),std::max_element(counts[indiv].begin(),counts[indiv].end()));;
		int max2;

		int M = -1;
		for(int b=0;b<4;++b){

			if(b!=max1 and int(counts[indiv][b])>M){

				M = counts[indiv][b];
				max2 = b;

			}

		}

		/*
		 * probability of homozigosis:
		 *
		 * P(allele1 good coverage) * P(allele2 good coverage) = P(allele good coverage)^2
		 *
		 */
		double p_h = indiv==0 ? Poi1[counts[indiv][max1]] : Poi2[counts[indiv][max1]];
		p_h *= p_h;

		/*
		 * probability of heterozigosis
		 *
		 * P(allele1 good coverage) * P(allele2 good coverage)
		 *
		 */
		double p_e = indiv==0 ? Poi1[2*counts[indiv][max1]] : Poi2[2*counts[indiv][max1]];
		p_e *= indiv==0 ? Poi1[2*counts[indiv][max2]] : Poi2[2*counts[indiv][max2]];

		/*
		 * homozigosis more likely
		 */
		if(p_h>p_e){

			if(		counts[indiv][max1] < mcov or
					counts[indiv][max1] > max_cov or
					double(counts[indiv][max1]/double(cluster_size)<(1-err_allowed)))
				return "NNNN";

			out[indiv*2] = int_to_base(max1);
			out[indiv*2+1] = int_to_base(max1);

		}

		/*
		 * heterozigosis more likely
		 */
		if(p_e>=p_h){

			/*
			 * check coverage/frequencies
			 */
			if(	(not diploid) or
				counts[indiv][max1]+counts[indiv][max2] < mcov or
				counts[indiv][max1]+counts[indiv][max2] > max_cov or
				double(counts[indiv][max1]+counts[indiv][max2])/double(cluster_size)<(1-err_allowed)
			)
				return "NNNN";

			out[indiv*2] = int_to_base(max1);
			out[indiv*2+1] = int_to_base(max2);

		}

	}

	return out;

}

/*
 * check if the haplotipes are a valid variant
 */
bool variant(string hapl){

	//if haplotipe calling failed
	if(hapl[0]=='N') return false;

	//if same haplotype SET: e.g. {A,C} {C,A}
	if((hapl[0]==hapl[2] and hapl[1]==hapl[3]) or (hapl[0]==hapl[3] and hapl[1]==hapl[2])) return false;

	//if number of distinct characters is not 2
	string copy = hapl;
	std::sort(hapl.begin(), hapl.end());
	copy.erase(unique(copy.begin(), copy.end()), copy.end());
	if(copy.size()!=2) return false;

	return true;

}

/*
 * get a set of reads from file given their rank.
 *
 * reads must be sorted!
 *
 * output: reads and their IDs
 */
void get_reads(string fasta_path, vector<uint64_t> & read_ranks, vector<string> & out_DNA){

	ifstream fasta;
	fasta.open(fasta_path);

	uint64_t j = 0;//current read in file

	string ID;
	string DNA;
	string line;

	//get read ID
	getline(fasta,ID);

	//read DNA
	getline(fasta,line);
	DNA = line;

	//skip newlines, append DNA segments
	while(not fasta.eof() && line[0] != '>'){

		getline(fasta,line);
		if(line[0] != '>') DNA.append(line);

	}

	for(uint64_t i = 0; i<read_ranks.size();++i){

		while(j < read_ranks[i]){

			ID = line;

			//read DNA
			getline(fasta,line);
			DNA = line;
			while(not fasta.eof() && line[0] != '>'){

				getline(fasta,line);
				if(line[0] != '>') DNA.append(line);

			}

			++j;

		}

		out_DNA.push_back(DNA);

	}

	fasta.close();

}

/*
 * input: range of SA corresponding to an event, plus its associated reads and IDs
 *
 * print reads aligned in the events
 */
void print_event(event & e, vector<string> & r, string & right_context){

	int max_left = 0;

	cout << right_context << endl;

	for(auto s : e.read_prefixes) max_left = s>max_left ? s : max_left;

	for(int i = 0;i<r.size();++i){

		int shift = max_left - e.read_prefixes[i];

		for(int j=0;j<shift;++j) cout << " ";

		for(int j=0;j<e.read_prefixes[i]-1;++j) cout << r[i][j];

		cout << "|" << r[i][e.read_prefixes[i]-1] << " " << int(e.indiv[i]?1:0) << "|";

		for(int j=e.read_prefixes[i];j<r[i].size();++j) cout << r[i][j];

		cout << endl;
	}

	cout << endl;

}


//find consensus in a set of strings, aligning them on the right
//returns empty string if consensus cannot be found
string consensus(vector<string> & S){

	if(S.size()==0) return "";

	int len = S[0].length();

	//consensus of reverse strings
	auto C = vector<vector<int>>(len, vector<int>(4,0));

	for(auto s:S)
		for(int i=0;i<len;++i)
			C[i][base_to_int(s[i])]++;

	string cons;
	bool found_consensus = true;

	for(int i=0;i<len;++i){

		bool base_found = false;

		for(int base=0;base<4;++base){

			double f = double(C[i][base])/double(S.size());

			if(f >= min_consensus_freq ){

				cons += int_to_base(base);
				base_found = true;

			}

		}

		found_consensus = base_found and found_consensus;

	}

	if(found_consensus) return cons;

	return "";

}

/*
 * Hamming distance on strings of same length
 */
int dH(string & a, string & b){

	assert(a.size()==b.size());

	int d=0;

	for(int i=0;i<a.size();++i) d += a[i]!=b[i];

	return d;

}

/*
 *
 */
int distance(string & a, string & b){
	return dH(a,b);
}


/*
 *
 * aligns the reads and calls the event
 */
void call_event(event & e, vector<string> & left_contexts, string & right_context, string fasta_path, ofstream & out_file){

	//cluster left parts (BWT symbol included) by Hamming distance (exact match)
	//cluster representative is longest string.
	vector<string> I1;
	vector<string> I2;

	for(int i = 0;i<left_contexts.size();++i){

		if(e.indiv[i]){//individual 2

			I2.push_back(left_contexts[i]);

		}else{//individual 1

			I1.push_back(left_contexts[i]);

		}

	}

	string c1 = consensus(I1);
	string c2 = consensus(I2);

	if(k_left==1){

		c1 = string();
		c2 = string();

		c1 += e.snp[0];
		c2 += e.snp[1];

	}

	if(c1.length()==0 or c2.length()==0) return;

	int d = distance(c1, c2);

	/*
	 * Last filters: distance between left contexts and minimum coverage
	 */

	if(k_left == 1 or (d <= max_snvs and I1.size()>= mcov_out and I2.size()>= mcov_out) ){

		if(discoSNP){

			/*
			 * sample 1
			 */

			string ID = ">SNP_higher_path_";
			ID.append(to_string(id_nr));
			ID.append("|P_1:");
			ID.append(to_string(right_context.size()));
			ID.append("_");
			ID += e.snp[0];
			ID.append("/");
			ID += e.snp[1];
			ID.append("|");
			ID.append("high");//TODO high/low events
			ID.append("|nb_pol_1");

			out_file << ID << endl;

			string DNA = c1;
			DNA.append(right_context);
			out_file << DNA << endl;

			/*
			 * sample 2
			 */

			ID = ">SNP_lower_path_";
			ID.append(to_string(id_nr));
			ID.append("|P_1:");
			ID.append(to_string(right_context.size()));
			ID.append("_");
			ID += e.snp[0];
			ID.append("/");
			ID += e.snp[1];
			ID.append("|");
			ID.append("high");
			ID.append("|nb_pol_1");

			out_file << ID << endl;

			DNA = c2;
			DNA.append(right_context);
			out_file << DNA << endl;

			id_nr++;

		}else{

			/*
			 * sample 1
			 */

			string ID = ">id:";
			ID.append(to_string(id_nr));
			ID.append(" ");
			ID.append(c1);
			ID.append(" sample:1 ");
			ID.append(e.hapl[0]!=e.hapl[1]?"0|1":"1|1");
			ID.append(" cov:");
			ID.append(to_string(I1.size()));
			out_file << ID << endl;

			string DNA = c1;
			DNA.append(right_context);

			out_file << DNA << endl;

			/*
			 * sample 2
			 */

			ID = ">id:";
			ID.append(to_string(id_nr));
			ID.append(" ");
			ID.append(c2);
			ID.append(" sample:2 ");
			ID.append(e.hapl[2]!=e.hapl[3]?"0|1":"1|1");
			ID.append(" cov:");
			ID.append(to_string(I2.size()));
			out_file << ID << endl;

			DNA = c2;
			DNA.append(right_context);

			out_file << DNA << endl;


			id_nr++;

		}

	}

}

void call_events(vector<event> & event_buffer, string fasta_path, ofstream & out_file){

	cout << "Extracting reads from fasta file ... " << flush;
	vector<uint64_t> read_ranks;

	for(event e : event_buffer){

		read_ranks.push_back(e.good_lcp_t);

		for(auto r : e.read_ranks)
			read_ranks.push_back(r);

	}

	//get read ranks, sort and remove duplicates
	std::sort( read_ranks.begin(), read_ranks.end() );
	auto last = std::unique( read_ranks.begin(), read_ranks.end() );
	read_ranks.erase(last, read_ranks.end());

	vector<string> reads;

	get_reads(fasta_path, read_ranks, reads);

	cout << "Done." << endl;

	int last_perc=0;
	int perc=0;

	cout << "Processing events ... " << endl;

	for(int j=0;j<event_buffer.size();++j){

		auto e = event_buffer[j];

		/*
		 * get right context
		 */
		auto c = e.good_lcp_t;
		auto it = std::find(read_ranks.begin(), read_ranks.end(), c);
		assert(it != read_ranks.end());
		auto idx = std::distance(read_ranks.begin(), it);

		string right_context = reads[idx].substr(e.good_lcp_s, e.good_lcp_len);

		/*
		 * get left contexts
		 */
		//reads of the event
		vector<string> left_contexts;

		for(int k=0;k< e.read_ranks.size(); ++k){

			auto c = e.read_ranks[k];
			auto it = std::find(read_ranks.begin(), read_ranks.end(), c);
			assert(it != read_ranks.end());
			auto idx = std::distance(read_ranks.begin(), it);

			int start = e.read_prefixes[k]-k_left;
			int len = k_left;

			left_contexts.push_back(reads[idx].substr(start, len));

			e.read_prefixes[k] = k_left;

		}

		call_event(e, left_contexts, right_context, fasta_path, out_file);
		//print_event(e, left_contexts, right_context);

		perc = (100*j)/event_buffer.size();

		if(perc>last_perc+4){

			cout << " " << perc << "%" << endl;
			last_perc = perc;

		}

	}

}

/*
 * scans EGSA, clusters and finds interesting clusters. In chunks, extracts the reads
 * from interesting clusters and aligns them.
 */
void find_events(string & egsa_path, string & clusters_path, string fasta_path, string out_path){

	ifstream egsa;
	egsa.open(egsa_path, ios::in | ios::binary);

	ifstream clusters;
	clusters.open(clusters_path, ios::in | ios::binary);

	ofstream out_file;
	out_file.open(out_path, ios::out);

	uint64_t tot_mem = 0;//memory occupied by buffer
	vector<event> event_buffer;

	uint64_t i = 0;//position on suffix array

	//read first egsa entry
	t_GSA e = read_el(egsa, bcr);

	while(not clusters.eof()){

		//read entry in clusters

		uint64_t start;
		uint16_t length;

		clusters.read((char*)&start, sizeof(uint64_t));
		clusters.read((char*)&length, sizeof(uint16_t));

		while(i < start){

			e = read_el(egsa, bcr);
			++i;

		}

		//now e is the start-th entry in egsa.

		auto counts = vector<vector<unsigned int> >(2,vector<unsigned int>(4,0));

		//extract a segment of the EGSA
		vector<uint64_t> read_rank_segment;
		vector<uint16_t> read_pref_segment;
		string bwt_segment;
		vector<bool> sample_segment;

		uint16_t good_lcp_s = 0;
		uint16_t good_lcp_len = 0;
		uint64_t good_lcp_t = 0;

		for(uint64_t l = 0;l<length; ++l){

			/*
			 * find maximum LCP
			 */
			if(e.lcp > good_lcp_len){

				good_lcp_len = e.lcp;
				good_lcp_t = e.text;
				good_lcp_s = e.suff;

			}

			bool sample = e.text < nr_reads1 ? 0 : 1;

			if(e.suff>=k_left and k_left>1){
				bwt_segment += e.bwt;
				read_rank_segment.push_back(e.text);
				read_pref_segment.push_back(e.suff);
				sample_segment.push_back(sample);
			}

			counts[sample?1:0][base_to_int(e.bwt)]++;

			e = read_el(egsa, bcr);
			++i;

		}

		string type = event_type(counts);

		if(variant(type)){

			good_lcp_len = good_lcp_len > k_right ? k_right : good_lcp_len;

			string snp = "NN";

			if(type[0]!=type[2]){
				snp[0]=type[0];
				snp[1]=type[2];
			}else if (type[0]!=type[3]){
				snp[0]=type[0];
				snp[1]=type[3];
			}else if (type[1]!=type[2]){
				snp[0]=type[1];
				snp[1]=type[2];
			}else if (type[1]!=type[3]){
				snp[0]=type[1];
				snp[1]=type[3];
			}

			vector<uint64_t> read_rank_segment_filt;
			vector<uint16_t> read_pref_segment_filt;
			vector<bool> sample_segment_filt;

			for(uint64_t l = 0;l<read_rank_segment.size(); ++l){

				if((sample_segment[l]==0 and bwt_segment[l]==snp[0]) or (sample_segment[l]==1 and bwt_segment[l]==snp[1])){
					read_rank_segment_filt.push_back(read_rank_segment[l]);
					read_pref_segment_filt.push_back(read_pref_segment[l]);
					sample_segment_filt.push_back(sample_segment[l]);
				}

			}

			event_buffer.push_back({read_rank_segment_filt,read_pref_segment_filt,sample_segment_filt, good_lcp_t, good_lcp_s, good_lcp_len, type, snp});

			//max memory that we will use now and in later steps TODO
			/*tot_mem +=  2*read_rank_segment.size()*sizeof(uint64_t) +
						read_pref_segment.size()*sizeof(uint16_t) +
						2*read_rank_segment.size()*k_left;*/

		}

		/*if(tot_mem > MAX_MEM*1048576){

			cout << "Emptying event buffer. " << event_buffer.size() << " events detected" <<endl;

			//extract reads, perform alignments, call events, and store them to file
			call_events(event_buffer, fasta_path, out_file);

			//clear buffer
			event_buffer = vector<event>();
			tot_mem = 0;

		}*/

	}

	//if there are events left to process
	if(event_buffer.size()!=0){

		cout << "Emptying event buffer. " << event_buffer.size() << " events detected" <<endl;

		//extract reads, perform alignments, call events, and store them to file
		call_events(event_buffer, fasta_path, out_file);

	}


	out_file.close();
	clusters.close();
	egsa.close();

}

/*
 * compute coverage statistics
 */
void statistics(string & egsa_path, string & clusters_path){

	ifstream egsa;
	egsa.open(egsa_path, ios::in | ios::binary);

	ifstream clusters;
	clusters.open(clusters_path, ios::in | ios::binary);

	uint64_t i = 0;//position on suffix array

	vector<int> cov_1 = vector<int>(200, 0);
	vector<int> cov_2 = vector<int>(200, 0);
	vector<int> cov = vector<int>(200, 0);

	//read first egsa entry
	t_GSA e = read_el(egsa, bcr);

	uint64_t n_bases = 0;//number of nucleotides

	n_bases += 	e.bwt=='A' or e.bwt=='a' or
				e.bwt=='C' or e.bwt=='c' or
				e.bwt=='G' or e.bwt=='g' or
				e.bwt=='T' or e.bwt=='t';

	uint64_t tot1=0;
	uint64_t tot2=0;


	while(not clusters.eof()){

		//read entry in clusters

		uint64_t start;
		uint16_t length;

		clusters.read((char*)&start, sizeof(uint64_t));
		clusters.read((char*)&length, sizeof(uint16_t));

		while(i < start){

			e = read_el(egsa, bcr);

			n_bases += 	e.bwt=='A' or e.bwt=='a' or
						e.bwt=='C' or e.bwt=='c' or
						e.bwt=='G' or e.bwt=='g' or
						e.bwt=='T' or e.bwt=='t';
			++i;

		}

		int cov1_cnt = 0;
		int cov2_cnt = 0;


		for(uint64_t l = 0;l<length; ++l){

			bool sample = e.text < nr_reads1 ? 0 : 1;

			if(sample) cov2_cnt++;
			else cov1_cnt++;

			e = read_el(egsa, bcr);

			n_bases += 	e.bwt=='A' or e.bwt=='a' or
						e.bwt=='C' or e.bwt=='c' or
						e.bwt=='G' or e.bwt=='g' or
						e.bwt=='T' or e.bwt=='t';

			++i;

		}

		/*
		 * cut clusters that contain less than 3 elements per sample (noise)
		 */
		if(cov1_cnt+cov2_cnt<cov.size()) cov[cov1_cnt+cov2_cnt]++;
		if(cov1_cnt<cov_1.size()) cov_1[cov1_cnt]++;
		if(cov2_cnt<cov_2.size()) cov_2[cov2_cnt]++;

	}

	int max=0;
	for(auto x:cov_1) max=x>max?x:max;
	for(auto x:cov_2) max=x>max?x:max;

	int scale = max/100;

	cout << "\nSample 1, distribution of cluster size: "<< endl;
	cout << "\ncov\tf\tcumul" << endl;
	for(int i=0;i<cov_1.size() and cov_1[i] > max/200 ;++i){

		if(i%1==0){
			cout << i << "\t" << flush;
			for(int j=0;j<cov_1[i]/scale;++j) cout << "-" << flush;
			cout << "   " << cov_1[i] << endl;
		}

	}

	cout << "\nSample 2, distribution of cluster size: "<< endl;
	cout << "\ncov\tf\tcumul" << endl;
	for(int i=0;i<cov_2.size()  and cov_2[i] > max/200;++i){

		if(i%1==0){
			cout << i  << "\t" << flush;
			for(int j=0;j<cov_2[i]/scale;++j) cout << "-" << flush;
			cout << "   " << cov_2[i] << endl;
		}

	}

	/*
	 * now cut the noise
	 */

	/*
	 * Check if there is a peak in the first 2 bases (noise)
	 *

	bool peak_noise1=cov_1[0]>cov_1[1];
	bool peak_noise2=cov_2[0]>cov_2[1];

	for(int i=1;i<=2 and i<cov_1.size();++i)
		if(cov_1[i] >= cov_1[i-1] and cov_1[i] > cov_1[i+1]){
			cout << i << endl;
			peak_noise1=true;
		}

	for(int i=1;i<=2 and i<cov_2.size();++i)
		if(cov_2[i] >= cov_2[i-1] and cov_2[i] > cov_2[i+1])
			peak_noise2=true;

	int local_min1=0;
	int local_min2=0;

	if(not peak_noise1){
		local_min1 = mcov_out;
	}else{

		bool found = false;
		i=1;
		while(not found and i < cov_1.size()-1){
			if(cov_1[i] <= cov_1[i-1] and cov_1[i] < cov_1[i+1]){
				found=true;
				local_min1=i;
			}
			++i;
		}

	}

	if(not peak_noise2){
		local_min2 = mcov_out;
	}else{
		i=1;
		bool found=false;
		while(not found and i < cov_2.size()-1){
			if(cov_2[i] <= cov_2[i-1] and cov_2[i] < cov_2[i+1]){
				found=true;
				local_min2=i;
			}
			++i;
		}
	}

	local_min1 = local_min1>4?local_min1-1:local_min1;
	local_min2 = local_min2>4?local_min2-1:local_min2;
	 */


	int local_min1 = mcov_out;
	int local_min2 = mcov_out;

	for(i=0;i<local_min1;++i) cov_1[i]=0;
	for(i=0;i<local_min2;++i) cov_2[i]=0;

	uint64_t tot=0;

	Poi1 = vector<double>(cov_1.size());
	Poi2 = vector<double>(cov_2.size());
	Cumul1 = vector<double>(cov_1.size());
	Cumul2 = vector<double>(cov_2.size());

	for(auto x:cov_1) tot1+=x;
	for(auto x:cov_2) tot2+=x;
	for(auto x:cov) tot+=x;

	for(int i=1;i<cov_1.size();++i){
		Poi1[i] = double(cov_1[i])/double(tot1);
		Cumul1[i] = Cumul1[i-1] + Poi1[i];
	}
	for(int i=1;i<cov_2.size();++i){
		Poi2[i] = double(cov_2[i])/double(tot2);
		Cumul2[i] = Cumul2[i-1] + Poi2[i];
	}

	max_cov1=0;
	max_cov2=0;
	mcov1=0;
	mcov2=0;

	int mean1 = std::distance(cov_1.begin(),std::max_element(cov_1.begin(),cov_1.end()));
	int mean2 = std::distance(cov_2.begin(),std::max_element(cov_2.begin(),cov_2.end()));

	/*
	 * compute tails of Poisson distributions
	 */

	max_cov1 = cov_1.size()-1;
	for(int i=cov_1.size()-1; i>=0 ;--i)
		if(1-Cumul1[i] <= pval/2) max_cov1 = i;

	max_cov2 = cov_2.size()-1;
	for(int i=cov_2.size()-1; i>=0 ;--i)
		if(1-Cumul2[i] <= pval/2) max_cov2 = i;

	for(int i=0; i<cov_1.size() ;++i)
		if(Cumul1[i] <= pval/2) mcov1 = i;

	for(int i=0; i<cov_2.size() ;++i)
		if(Cumul2[i] <= pval/2) mcov2 = i;

	mcov1 = mcov1<local_min1?local_min1:mcov1;
	mcov2 = mcov2<local_min2?local_min2:mcov2;

	cout << "\nCluster sizes allowed in the two samples: [" << mcov1 << "," << max_cov1 << "] and [" << mcov2 << "," << max_cov2 << "]" << endl;

	clusters.close();
	egsa.close();

}


int main(int argc, char** argv){


	srand(time(NULL));

	if(argc < 3) help();

	int opt;
	while ((opt = getopt(argc, argv, "hi:n:e:p:v:L:R:m:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'b':
				bcr = true;
			break;
			case 'D':
				diploid = true;
			break;
			case 'i':
				input = string(optarg);
				//cout << "input = " << input << "\n";
			break;
			case 'n':
				nr_reads1 = atoi(optarg);
				//cout << "k = " << k << "\n";
			break;
			case 'm':
				mcov_out = atoi(optarg);
				//cout << "m = " << m << "\n";
			break;
			case 'L':
				k_left = atoi(optarg);
				//cout << "k = " << k << "\n";
			break;
			case 'R':
				k_right = atoi(optarg);
				//cout << "m = " << m << "\n";
			break;
			case 'e':
				err_allowed = atof(optarg);
				//cout << "k = " << M << "\n";
			break;
			case 'p':
				pval = atof(optarg);
				//cout << "k = " << M << "\n";
			break;
			/*case 'r':
				read_len = atof(optarg);
				//cout << "k = " << M << "\n";
			break;*/
			/*case 'd':
				discoSNP=true;
				//cout << "k = " << M << "\n";
			break;*/
			case 'v':
				max_snvs = atoi(optarg);
			break;
			default:
				help();
			return -1;
		}
	}

	MAX_MEM = MAX_MEM==0?MAX_MEM_def:MAX_MEM;
	err_allowed = (err_allowed==0?err_allowed_def:err_allowed);
	k_left = k_left==0?k_left_def:k_left;
	k_right = k_right==0?k_right_def:k_right;
	pval = pval==0?pval_def:pval;
	read_len = read_len==0?read_len_def:read_len;
	max_snvs = max_snvs==0?max_snvs_def:max_snvs;
	mcov_out = mcov_out==0?mcov_out_def:mcov_out;

	if(pval <= 0 or pval > 1){

		cout << "Error: argument of -p must be in (0,1]" << endl;
		help();

	}

	if(err_allowed < 0 or err_allowed >= 0.5){

		cout << "Error: argument of -f must be in [0,0.5)" << endl;
		help();

	}

	if(input.compare("")==0 or nr_reads1 == 0) help();

	string egsa_path = input;
	egsa_path.append(".0.gesa");

	{

		ifstream ifs(egsa_path);

		if(not ifs.good()){

			cout << "\nERROR: Could not find EGSA file \"" << egsa_path << "\"" << endl << endl;
			help();


		}

		ifs.close();

	}

	cout << "This is clust2snp." << endl <<
			"Input index file: " << egsa_path << endl <<
			"Tolerated error rate : " << err_allowed << endl <<
			"p-value : " << pval << endl <<
			"Left-extending GSA ranges by " << k_left << " bases." << endl <<
			"Right context length: at most " << k_right << " bases." << endl;

	if(diploid)
		cout << "Diploid samples" << endl;
	else
		cout << "Haploid samples" << endl;

	frac_cov = (1-min_k/read_len)*pow(1-err_rate,min_k);

	string clusters_path = input;
	clusters_path.append(".clusters");

	{

		ifstream ifs(clusters_path);

		if(not ifs.good()){

			cout << "\nERROR: Could not find BWT clusters file \"" << clusters_path << "\"" << endl << endl;
			help();

		}

		ifs.close();

	}

	string filename_out = input;
	filename_out = filename_out.substr(0,filename_out.rfind(".fast"));
	filename_out.append(".snp.fasta");

	cout << "Output events will be stored in " << filename_out << endl;

	statistics(egsa_path, clusters_path);
	find_events(egsa_path, clusters_path, input, filename_out);

	cout << "Done. " <<endl;

}
