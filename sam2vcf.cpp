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
#include "include.hpp"
#include <algorithm>

using namespace std;

int indel_deduplicate_def = 5;//if 2 indels are within this number of bases, keep only one of the two.
int indel_deduplicate = 0;

bool non_isolated = true;

void help(){

	cout << "sam2vcf [OPTIONS]" << endl << endl <<
	"Converts the aligned calls (with bwa-mem) 'calls.sam' of clust2snp into a vcf file 'calls.sam.vcf'." << endl <<
	"Options:" << endl <<
		"-h          Print this help." << endl <<
		"-x          Disable non-isolated SNPs (default: enabled)." << endl <<
		"-s <arg>    Input SAM file. REQUIRED" << endl <<
		"-d <arg>    Keep only one indel in pairs within <arg> bases. Default: " <<  indel_deduplicate_def << ")" << endl;
	exit(0);
}

struct vcf_entry{

	string chr;
	uint64_t pos;

	string REF;
	string ALT;

	bool indel;

	bool operator<(const vcf_entry & a) const{

		if(chr.compare(a.chr) < 0 ) return true;
		else if(chr.compare(a.chr) == 0 ){

			return pos < a.pos;

		}else return false;

	}

	bool operator==(const vcf_entry & a) const{

		return (chr.compare(a.chr)== 0) and pos == a.pos and REF.compare(a.REF) == 0 and ALT.compare(a.ALT) == 0;

	}


};

int main(int argc, char** argv){

	if(argc < 2) help();

	string infile;

	int opt;
	while ((opt = getopt(argc, argv, "xhd:s:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'x':
				non_isolated = false;
			break;
			case 'd':
				indel_deduplicate = atoi(optarg);
			break;
			case 's':
				infile = string(optarg);
			break;
			default:
				help();
			return -1;
		}
	}

	indel_deduplicate = indel_deduplicate==0 ? indel_deduplicate_def : indel_deduplicate;

	if(infile.compare("")==0) help();

	string outfile = infile;
	outfile.append(".vcf");

	ifstream is(infile);
	ofstream of(outfile);

	string str;

	of << "#CHROM\tPOS\tID\tREF\tALT\tINFO" << endl;
	//cout << "#CHROM\tPOS\tID\tREF\tALT\tINFO" << endl;

	vector<vcf_entry> VCF;

	while(getline(is, str)){

		if(str[0]!='@'){//skip header

			string type;//INDEL or SNP
			string REF;//reference allele
			string ALT;//alternative allele
			string REF_dna;//reference dna
			string ALT_dna;//alternative dna
			uint64_t pos;//alignment position

			string str_copy = str;

			char *p = strtok((char*)str_copy.c_str(), "_");
			type = string(p);

			str_copy = str;
			p = strtok((char*)str_copy.c_str(), "|");
			p = strtok(NULL, "|");
			string tmp1 = string(p);//P_1:30_T/

			p = strtok(NULL, "|");
			p = strtok(NULL, "|");
			string tmp2 = string(p);//nb_pol_1_DNA

			p = strtok((char*)tmp1.c_str(), "_");
			p = strtok(NULL, "_");//1:30

			string tmp3 = string(p);

			p = strtok(NULL, "_");
			string REFALT = string(p);

			if(REFALT[0]=='/'){

				ALT = REFALT.substr(1);
				REF = string();

			}else{

				p = strtok((char*)REFALT.c_str(), "/");
				REF = string(p);
				p = strtok(NULL, "/");
				if(p!=NULL)
					ALT = string(p);
				else
					ALT = string();

			}

			//cout << type << " '" << REF << "' -> '" << ALT << "'" << flush;

			p = strtok((char*)tmp2.c_str(), "_");
			p = strtok(NULL, "_");
			p = strtok(NULL, "_");
			p = strtok(NULL, "_");
			tmp2 = string(p);

			p = strtok((char*)tmp3.c_str(), ":");
			p = strtok(NULL, "_");

			/*
			 * snp_pos will contain the starting position of the snp/indel on REF on the forward strand.
			 * Note: in case of indel, snp_pos contains the position before the beginning of the indel, since
			 * we want to report at least 1 base in both REF/ALT in the VCF
			 */
			int snp_pos = atoi(p);

			std::istringstream iss(tmp2);
			std::string token;
			getline(iss, token, '\t'); //ALT DNA
			ALT_dna = token;

			getline(iss, token, '\t'); //flag
			unsigned int f = atoi(token.c_str());

			bool reversed = (f & (unsigned int)16) != 0;
			bool indel = type.compare("INDEL")==0;

			if(reversed){//then ALT_dna is on RC

				ALT_dna = RC(ALT_dna);
				REF = RC(REF);
				ALT = RC(ALT);

				if(indel){

					snp_pos--;

				}

			}

			getline(iss, token, '\t'); //chr

			string chr = token;

			getline(iss, token, '\t'); //pos

			pos = atoi(token.c_str());

			getline(iss, token, '\t');
			getline(iss, token, '\t');
			getline(iss, token, '\t');
			getline(iss, token, '\t');
			getline(iss, token, '\t');
			getline(iss, token, '\t');//REF_dna

			REF_dna = token;

			//adjust snp_pos in the case we are on FW strand
			if(not reversed){

				if(indel){

					if(REF.length()>0){

						//insert in REF
						int indel_len = REF.length();
						snp_pos = REF_dna.length() - snp_pos - indel_len -1;

					}else{

						//insert in ALT
						snp_pos = REF_dna.length() - snp_pos -1;

					}

				}else{

					snp_pos = (REF_dna.length() - snp_pos)-1;

				}

			}

			if((f==0 or f==16) and snp_pos >= 0 and chr.compare("*") != 0){

				vcf_entry v;

				if(indel){

					if(REF.length()>0){

						v = {
										chr,
										pos + snp_pos,
										REF_dna.substr(snp_pos,REF.length()+1),
										ALT_dna.substr(snp_pos,1),
										indel
						};

					}else{

						v = {
										chr,
										pos + snp_pos,
										REF_dna.substr(snp_pos,1),
										ALT_dna.substr(snp_pos,ALT.length()+1),
										indel
						};

					}


				}else{

					v = {
									chr,
									pos + snp_pos,
									REF_dna.substr(snp_pos,1),
									ALT_dna.substr(snp_pos,1),
									indel
					};

				}

				VCF.push_back(v);

				/*
				 * find non-isolated SNPs
				 */

				if(non_isolated){

					if(indel){

						if(reversed){

							//non-isolated SNPs are on the right of the end position of indel

							int indel_length = REF.length() > 0 ? REF.length() : ALT.length();

							int L = std::max(REF_dna.length(), ALT_dna.length());//length of fragment containing the insert

							int len_right = L - (snp_pos+indel_length) -1; //length of right part in common (with potential SNPs)

							int snp_pos_ref = REF.length() > 0 ? snp_pos + indel_length +1 : snp_pos +1;
							int snp_pos_alt = REF.length() > 0 ? snp_pos +1 : snp_pos + indel_length +1;

							for(int i=0;i<len_right;++i){

								if(REF_dna[snp_pos_ref+i] != ALT_dna[snp_pos_alt+i]){

									vcf_entry v  = {
														chr,
														pos + i,
														REF_dna.substr(snp_pos_ref+i,1),
														ALT_dna.substr(snp_pos_alt+i,1),
														false
										};

									VCF.push_back(v);

								}

							}

						}else{

							//non-isolated SNPs are on the left of the start position of indel

							for(int i=0;i<snp_pos;++i){

								if(REF_dna[i] != ALT_dna[i]){

									vcf_entry v  = {
														chr,
														pos + i,
														REF_dna.substr(i,1),
														ALT_dna.substr(i,1),
														false
										};

									VCF.push_back(v);

								}

							}

						}

					}else{

						if(reversed){

							//non-isolated SNPs are on the right

							for(int i=snp_pos+1;i<REF_dna.length();++i){

								if(REF_dna[i] != ALT_dna[i]){

									vcf_entry v  = {
														chr,
														pos + i,
														REF_dna.substr(i,1),
														ALT_dna.substr(i,1),
														false
										};

									VCF.push_back(v);

								}

							}

						}else{

							//non-isolated SNPs are on the left

							for(int i=0;i<snp_pos;++i){

								if(REF_dna[i] != ALT_dna[i]){

									vcf_entry v  = {
														chr,
														pos + i,
														REF_dna.substr(i,1),
														ALT_dna.substr(i,1),
														false
										};

									VCF.push_back(v);

								}

							}

						}

					}

				}

			}

		}

	}

	//remove duplicates

	uint64_t old_size = VCF.size();

	std::sort(VCF.begin(),VCF.end());

	/*for(int i=0;i<VCF.size()-1;++i){

		if(VCF[i]==VCF[i+1]){

			cout << VCF[i].chr << "\t" << VCF[i].pos << "\t" << ".\t" << VCF[i].REF << "\t" << VCF[i].ALT << "\t" << (VCF[i].indel?"INDEL":"SNP") << endl;
			cout << VCF[i+1].chr << "\t" << VCF[i+1].pos << "\t" << ".\t" << VCF[i+1].REF << "\t" << VCF[i+1].ALT << "\t" << (VCF[i+1].indel?"INDEL":"SNP") << endl;

			cout << endl;

		}

	}*/

	VCF.erase( unique( VCF.begin(), VCF.end() ), VCF.end() );

	//remove close indels

	vector<vcf_entry> VCF_filt;
	for(int i=0;i<VCF.size()-1;++i){

		if(		VCF[i].indel and
				VCF[i+1].indel and
				VCF[i].chr.compare(VCF[i+1].chr)==0 and
				std::abs( int(VCF[i].pos)-int(VCF[i+1].pos) ) <= indel_deduplicate  ){


			VCF_filt.push_back(VCF[i]);//keep first indel
			i++;//skip second indel

		}else{

			VCF_filt.push_back(VCF[i]);

		}

	}

	VCF_filt.push_back(VCF[VCF.size()-1]);

	uint64_t new_size = VCF_filt.size();

	cout << (old_size-new_size) << " duplicates found. Saving remaining " << new_size << " unique SNPS/indels." << endl;

	uint64_t n_indels = 0;
	uint64_t n_snps = 0;

	for(auto v:VCF_filt){

		//cout << v.chr << "\t" << v.pos << "\t" << ".\t" << v.REF << "\t" << v.ALT << "\t" << (indel?"INDEL":"SNP") << endl;
		of << v.chr << "\t" << v.pos << "\t" << ".\t" << v.REF << "\t" << v.ALT << "\t" << (v.indel?"INDEL":"SNP") << endl;

		n_indels += v.indel;
		n_snps += (not v.indel);

	}

	cout << "Number of SNPs found: " << n_snps << endl;
	cout << "Number of indels found: " << n_indels << endl;

	is.close();
	of.close();

}
