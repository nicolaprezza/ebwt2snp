#include <fstream>

using namespace std;

typedef uint32_t int_text;
typedef uint32_t int_suff;
typedef uint32_t int_lcp;
typedef uint8_t int8;
typedef uint8_t dataTypelenSeq;	//length of the sequences (in biologic case 100)
typedef uint32_t dataTypeNSeq;	//number of sequences
typedef __int128 uint128_t;

typedef pair<uint64_t,uint32_t> coordinate;//suffix array coordinate (text, suff)

/*
 * EGSA
 */
typedef struct{

	int_text	text; //read nr
	int_suff	suff; //starting position of the suffix in the read
	int_lcp 	lcp;
	int8		bwt;

} t_GSA;

t_GSA read_el(ifstream & egsa, bool bcr){

	t_GSA e;

	if(bcr){

        dataTypeNSeq txt;
        uint8_t suf;
        uint8_t lcp;

        egsa.read((char*)&txt, sizeof(dataTypeNSeq));
        egsa.read((char*)&suf, sizeof(dataTypelenSeq));
        egsa.read((char*)&lcp, sizeof(dataTypelenSeq));
        egsa.read((char*)&e.bwt, sizeof(uint8_t));

        e.suff = suf;
        e.text = txt;
        e.lcp = lcp;

	}else{

		egsa.read((char*)&e, sizeof(int_text)+sizeof(int_suff)+sizeof(int_lcp)+sizeof(int8));

	}

	return e;

}

unsigned char int_to_base(int i){

	switch(i){

		case 0: return 'A'; break;
		case 1: return 'C'; break;
		case 2: return 'G'; break;
		case 3: return 'T'; break;

	}

	return 'A';

}

int base_to_int(unsigned char c){

	switch(c){

		case 'A': case 'a': return 0; break;
		case 'C': case 'c': return 1; break;
		case 'G': case 'g': return 2; break;
		case 'T': case 't': return 3; break;
		case 'N': case 'n': return rand()%4; break;

	}

	return 0;

}

unsigned char RC(unsigned char c){

	switch(c){

		case 'A': case 'a': return 'T'; break;
		case 'C': case 'c': return 'G'; break;
		case 'G': case 'g': return 'C'; break;
		case 'T': case 't': return 'A'; break;
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
