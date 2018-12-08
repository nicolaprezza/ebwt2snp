#ifndef INCLUDE_HPP_
#define INCLUDE_HPP_

#include <fstream>
#include <vector>
#include <cassert>

using namespace std;

typedef uint32_t int_text;
typedef uint32_t int_suff;
typedef uint32_t int_lcp;
typedef uint8_t int8;
typedef uint8_t dataTypelenSeq;	//length of the sequences (in biologic case 100)
typedef uint32_t dataTypeNSeq;	//number of sequences
typedef __int128 uint128_t;

typedef pair<uint64_t,uint32_t> coordinate;//suffix array coordinate (text, suff)
typedef pair<uint64_t,uint64_t> range_t;

const uint8_t TERM = '#';

/*
 * EGSA
 */
typedef struct{

	int_text	text; //read nr
	int_suff	suff; //starting position of the suffix in the read
	int_lcp 	lcp;
	int8		bwt;

} t_GSA;

std::ifstream::pos_type filesize(string filename){
    std::ifstream in(filename.c_str(), std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

/*
 * this class abstracts the EGSA type and allows reading from different formats (EGSA/BCR)
 */
class egsa_stream{

public:

	/*
	 * input: reads.fasta
	 *
	 * automatically detects the index files and format
	 *
	 */
	egsa_stream(string & input_path){

		string egsa_path = input_path;
		egsa_path.append(".gesa");

		EGSA.open(egsa_path, ios::in | ios::binary);

		if(EGSA.is_open()){

			egsa = true;

		}else{//else try BCR

			string LCP_path = input_path;
			LCP_path.append(".out.lcp");

			string BWT_path = input_path;
			BWT_path.append(".out");

			string GSA_path = input_path;
			GSA_path.append(".out.pairSA");

			LCP.open(LCP_path, ios::in | ios::binary);
			BWT.open(BWT_path, ios::in | ios::binary);
			GSA.open(GSA_path, ios::in | ios::binary);

			//note: we require LCP and BWT to exist. Missing GSA files
			//will generate 0-fields in the returned elements
			if(LCP.is_open() and BWT.is_open()){

				bcr = true;

			}else{

				cout << "Error: missing index files." << endl;
				exit(1);

			}

		}

	}

	/*
	 * returns true iff index files exist in the input folder
	 */
	bool index_exists(){

		return egsa or bcr;

	}

	bool eof(){

		if(egsa){

			return EGSA.eof();

		}else if (bcr){

			return LCP.eof();

		}

		cout << "Error: missing index files." << endl;
		exit(1);

	}

	/*
	 * overwrite default type byte-sizes (LCP, document array, suffix in read)
	 */
	void set_bytesizes(int lcp_size, int da_size, int suff_size){

		this->lcp_size = lcp_size;
		this->da_size = da_size;
		this->suff_size = suff_size;

	}

	t_GSA read_el(){

		t_GSA e;

		if(egsa){

			switch(da_size){

				case 1 : uint8_t x8; EGSA.read((char*)&x8, 1); e.text = x8;  break;
				case 2 : uint16_t x16; EGSA.read((char*)&x16, 2); e.text = x16; break;
				case 4 : uint32_t x32; EGSA.read((char*)&x32, 4); e.text = x32; break;
				case 8 : uint64_t x64; EGSA.read((char*)&x64, 8); e.text = x64; break;

			}

			switch(suff_size){

				case 1 : uint8_t x8; EGSA.read((char*)&x8, 1); e.suff = x8;  break;
				case 2 : uint16_t x16; EGSA.read((char*)&x16, 2); e.suff = x16; break;
				case 4 : uint32_t x32; EGSA.read((char*)&x32, 4); e.suff = x32; break;
				case 8 : uint64_t x64; EGSA.read((char*)&x64, 8); e.suff = x64; break;

			}

			switch(lcp_size){

				case 1 : uint8_t x8; EGSA.read((char*)&x8, 1); e.lcp = x8;  break;
				case 2 : uint16_t x16; EGSA.read((char*)&x16, 2); e.lcp = x16; break;
				case 4 : uint32_t x32; EGSA.read((char*)&x32, 4); e.lcp = x32; break;
				case 8 : uint64_t x64; EGSA.read((char*)&x64, 8); e.lcp = x64; break;

			}

			uint8_t x8;
			EGSA.read((char*)&x8, 1);
			e.bwt = x8;

		}else if(bcr){


			if(GSA.is_open()){

				switch(suff_size){

					case 1 : uint8_t x8; GSA.read((char*)&x8, 1); e.suff = x8;  break;
					case 2 : uint16_t x16; GSA.read((char*)&x16, 2); e.suff = x16; break;
					case 4 : uint32_t x32; GSA.read((char*)&x32, 4); e.suff = x32; break;
					case 8 : uint64_t x64; GSA.read((char*)&x64, 8); e.suff = x64; break;

				}

				switch(da_size){

					case 1 : uint8_t x8; GSA.read((char*)&x8, 1); e.text = x8;  break;
					case 2 : uint16_t x16; GSA.read((char*)&x16, 2); e.text = x16; break;
					case 4 : uint32_t x32; GSA.read((char*)&x32, 4); e.text = x32; break;
					case 8 : uint64_t x64; GSA.read((char*)&x64, 8); e.text = x64; break;

				}

			}else{

				e.suff = 0;
				e.text = 0;

			}

			switch(lcp_size){

				case 1 : uint8_t x8; LCP.read((char*)&x8, 1); e.lcp = x8;  break;
				case 2 : uint16_t x16; LCP.read((char*)&x16, 2); e.lcp = x16; break;
				case 4 : uint32_t x32; LCP.read((char*)&x32, 4); e.lcp = x32; break;
				case 8 : uint64_t x64; LCP.read((char*)&x64, 8); e.lcp = x64; break;

			}

			uint8_t x8;
			BWT.read((char*)&x8, 1);
			e.bwt = x8;

		}else{

			cout << "Error: missing index files." << endl;
			exit(1);

		}

		return e;

	}

private:

	bool egsa = false;
	bool bcr = false;

	//byte size of components
	int lcp_size = 1; //LCP values
	int da_size = 4; //document array (read number)
	int suff_size = 1; //position inside read

	//the EGSA index
	ifstream EGSA;

	//the BCR index
	ifstream LCP;
	ifstream BWT;
	ifstream GSA;//pairs

};

unsigned char int_to_base(int i){

	switch(i){

		case 0: return 'A'; break;
		case 1: return 'C'; break;
		case 2: return 'G'; break;
		case 3: return 'T'; break;
		case 4: return TERM; break;

	}

	return 'A';

}

int base_to_int(unsigned char c){

	switch(c){

		case 'A': case 'a': return 0; break;
		case 'C': case 'c': return 1; break;
		case 'G': case 'g': return 2; break;
		case 'T': case 't': return 3; break;
		case TERM: return 4; break;
		default: return rand()%4; break;

	}

	return 0;

}

unsigned char RC(unsigned char c){

	switch(c){

		case 'A': case 'a': return 'T'; break;
		case 'C': case 'c': return 'G'; break;
		case 'G': case 'g': return 'C'; break;
		case 'T': case 't': return 'A'; break;
		case TERM: return TERM; break;
		default: break;

	}

	return 'N';

}

string RC(string & s){

	if(s.length()==0) return string();

	string rc = s;

	for(int i=0;i<s.length();++i) rc[rc.length()-i-1] = RC(s[i]);

	return rc;

}


string rev(string & a){

	string reversed(a.rbegin(), a.rend());

	return reversed;

}


inline uint16_t clz_u128 (uint128_t u) {

	uint64_t hi = u>>64;
	uint64_t lo = u;

	return hi!= 0 ? __builtin_clzll(hi) : 64 + __builtin_clzll(lo);

}

inline int popcount128(uint128_t u){

	return __builtin_popcountl(uint64_t(u)) + __builtin_popcountl(uint64_t(u>>64));

}

class cons{

public:

	cons(int size){

		counts = vector<vector<int>>(size,vector<int>(4,0));
		C = string(size,'A');

	}

	unsigned char operator[](int i){
		return C[i];
	}

	void increment(int i, unsigned char b){

		int b_i = base_to_int(b);

		counts[i][b_i]++;

		if( counts[i][b_i] > counts[i][ base_to_int(C[i]) ])
			C[i] = b;

	}

	string to_string(){

		return C;

	}

private:

	string C;//the current consensus
	vector<vector<int>> counts;//base counts

};

uint8_t MASK_TERM = uint8_t(1);
uint8_t MASK_A = uint8_t(2);
uint8_t MASK_C = uint8_t(4);
uint8_t MASK_G = uint8_t(8);
uint8_t MASK_T = uint8_t(16);


/*
 * representation of a right-maximal substring (SA node) as a list of BWT intervals
 */
struct sa_node{

	//right-maximal substring: string W such that Wa_1, ..., Wa_k occur in the text for
	//at least k>=2 characters a_1, ..., a_k

	//Length k+1. Inclusive bwt range of W.chars[i] is <first[i], first[i+1]-1>
	//Range of W is <first[0], first[k]-1>
	//Equivalently, number of suffixes smaller than W.chars[i] is first[i]
	//vector<uint64_t> first;

	uint64_t first_TERM;
	uint64_t first_A;
	uint64_t first_C;
	uint64_t first_G;
	uint64_t first_T;
	uint64_t last;

	//depth = |W|
	uint64_t depth;

	uint8_t h;//head = first character of suffix = W[0]. Equal to 0 if W = empty string.

};

/*
 * suffix array leaf = BWT range (inclusive) of W.TERM, for some string W.
 *
 */
struct sa_leaf{

	//rn.first = first position of range. Equivalently, number of suffixes smaller than W.TERM (valid also if W.TERM does not occur)
	//rn.second = last position (excluded) of interval.  Equivalently, number of suffixes smaller than W.TERM + number of occurrences of W.TERM
	//if last == first, then W.TERM does not occur (however, 'first' is in any case number of suffixes smaller than W.TERM)
	range_t rn;

	//depth = |W.TERM|
	uint64_t depth;

};

uint64_t range_length(range_t r){
	assert(r.second >= r.first);
	return r.second - r.first;
}

uint64_t leaf_size(sa_leaf L){
	return range_length(L.rn);
}

uint64_t leaf_size(pair<sa_leaf, sa_leaf> P){
	return leaf_size(P.first) + leaf_size(P.second);
}


struct p_range{

	range_t A;
	range_t C;
	range_t G;
	range_t T;

};

struct p_rank{

public:

	uint64_t A;
	uint64_t C;
	uint64_t G;
	uint64_t T;

	p_rank operator+(const p_rank& a) const{

		return {
			a.A + A,
			a.C + C,
			a.G + G,
			a.T + T
		};

	}

	/*p_rank operator-(const p_rank& a) const{

		return {
			A-a.A,
			C-a.C,
			G-a.G,
			T-a.T
		};

	}*/

	bool operator==(const p_rank& a) const{

		return a.A == A and a.C == C and a.G == G and a.T == T;

	}

	bool operator!=(const p_rank& a) const{

		return a.A != A or a.C != C or a.G != G or a.T != T;

	}

	bool operator<=(const p_rank& a) const{

		return A <= a.A and C <= a.C and G <= a.G and T <= a.T;

	}

};

p_range fold_ranks(p_rank &a, p_rank &b){

	return {{a.A, b.A},{a.C, b.C},{a.G, b.G},{a.T, b.T}};

}

uint64_t popcount128(__uint128_t x){

	return __builtin_popcountll(uint64_t(x>>64)) + __builtin_popcountll( x & 0xFFFFFFFFFFFFFFFF );

}


#endif /* INCLUDE_HPP_ */

