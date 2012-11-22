/*
  fragrep c++ rebuild */

// std and iostream headers
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <list>
#include <stdio.h>
#include <map>
#include <string>
#include <sstream>

// fragrep's own includes
#include "omit_pattern_matching.h"
#include <sequtil/position_weight_matrix.h>

#include <sequtil/fasta.h>
#include <sequtil/bbq_util.h>

#define CHUNK_SIZE 4096


/** \mainpage

welcome to fragrep.

*/

using namespace fragrep;
using namespace bbq;

typedef omit_pattern_matching::match_list match_list;
typedef omit_pattern_matching::result_list result_list;
typedef omit_pattern_matching::complete_result complete_result;
typedef omit_pattern_matching::complete_result_pos_map complete_result_pos_map;
typedef omit_pattern_matching::int_vector int_vector;
typedef omit_pattern_matching::float_vector float_vector;
typedef omit_pattern_matching::match_list_vector match_list_vector;
typedef omit_pattern_matching::string_vector string_vector;
typedef omit_pattern_matching::int_vector_vector int_vector_vector;

std::ifstream ifs;
std::list<std::pair<char*,int> > chunk_list;
int last_chunk_length = 0;
int total_length = 0;
char* chunk_block;
bool final_chunk;
char fasta_header[CHUNK_SIZE];

/// \c genome sequence
std::string T;

std::string T_name;

/// \c k sequence fragments
string_vector s;


/// n[j] is the length of s[j]
int_vector n;

/// N is the length of T
int N;

/// m[j] number of mismatches allowed for occurences of s[j]
int_vector m;

/// d[j] specifies the maximum number of deletions for occurences of s[j]
int_vector d;

/// delta specifies the maximum number of blocks that can be deleted
int delta;

/// t[j] threshold for pwm match score of s[j]
float_vector t;

PWM_vector PWMs;

/// u[j] is an upper bound for the gap length between s[j-1] and s[j].
int_vector u;

/// l[j] is a lower bound for the gap length between s[j-1] and s[j].
int_vector l;

/// number of binding site fragments
int k;

void fill_reverse_segments();
void print_help_message();

/// Are the parameters passed to the program legal?
bool params_good;

/// number of command-line options
int num_of_opts;

/// consider reverse complemented segments also?
bool reverse;

/// fragments are defined by position weight matrices rather than
/// consensus strings if this flag is true.
bool use_pwms;

/// print simplifed output?
bool simple_output;

/// allow deletion of leading blocks at zero cost, suitable for
/// 'cut off' occurences in scaffolds of pre-assemblies
bool scaffold_cutoff_scoring;

/// false if -u option (unweighted) is given.
bool weighted;

/// specify extra output to stderr
bool verbose;

/// report only maximal match in case there is more than one match at each position
bool unique;

/// num of characters of extra output before a matching block
int pref_length;

/// num of characters of extra output after a matching block
int suff_length;


bool read_parameters(int argc, char** argv)
{
	if (argc<3)
	{
		print_help_message();		
		return false;
	}
	
	weighted = false;
	params_good = true;
	reverse = true;
	verbose = false;
	unique  = false;
	simple_output = false;
	scaffold_cutoff_scoring = false;
	num_of_opts = 0;
	delta = 0;
	int i;
	for (i=1; i<argc && argv[i][0]=='-'; i++)
	{
		num_of_opts++;
		switch (argv[i][1])
		{
			case 's' : {
				if (strlen(argv[i])>2)
				{
					std::cerr<<"Unknown option: "<<argv[i]<<"\n";
					return false;
				}
				if (argc<i+3)
				{
					std::cerr<<"Option -s takes two integer parameters.\n";
					return false;
				}
				if (!is_int(argv[i+1]) || !is_int(argv[i+2]))
				{
					std::cerr<<"Option -s takes two integer parameters.\n";
					return false;
				}
				pref_length = atoi(argv[i+1]);
				suff_length = atoi(argv[i+2]);				
				num_of_opts += 2;
				break;
			}
				// case 'd' : {
				// 				if (strlen(argv[i])>2)
				// 				{
				// 					std::cerr<<"Unknown option: "<<argv[i]<<"\n";
				// 					return false;
				// 				}
				// 				if (argc<i+2)
				// 				{
				// 					std::cerr<<"Option -d takes an integer parameter.\n";
				// 					return false;
				// 				}
				// 				if (!is_int(argv[i+1]))
				// 				{
				// 					std::cerr<<"Option -d takes an integer parameter.\n";
				// 					return false;
				// 				}
				// 				delta = atoi(argv[i+1]);
				// 				num_of_opts++;
				// 				break;
				// 			}
			case 'u' : {
				weighted = false;
				if (strlen(argv[i])>2)
				{
					std::cerr<<"Unknown option: "<<argv[i]<<"\n";
					return false;
				}
				break;
			}
			case 'v' : {
				verbose = true;
				if (strlen(argv[i])>2)
				{
					std::cerr<<"Unknown option: "<<argv[i]<<"\n";
					return false;
				}
				break;
			}
			case 'c' : {
				scaffold_cutoff_scoring = true;
				if (strlen(argv[i])>2)
				{
					std::cerr<<"Unknown option: "<<argv[i]<<"\n";
					return false;
				}
				break;
			}
			case 'q' : {
				unique = true;
				if (strlen(argv[i])>2)
				{
					std::cerr<<"Unknown option: "<<argv[i]<<"\n";
					return false;
				}
				break;
			}
			case 'r' : {
				reverse = true;
				if (strlen(argv[i])>2)
				{
					std::cerr<<"Unknown option: "<<argv[i]<<"\n";
					return false;
				}
				break;
			}
			case 'f' : {
				reverse = false;
				if (strlen(argv[i])>2)
				{
					std::cerr<<"Unknown option: "<<argv[i]<<"\n";
					return false;
				}
				break;
			}
			case 'S' : {
				simple_output = true;
				if (strlen(argv[i])>2)
				{
					std::cerr<<"Unknown option: "<<argv[i]<<"\n";
					return false;
				}
				break;
			} case 'H' : {
				print_help_message();
				break;
			} case 'w' : {
				weighted = true;
				if (strlen(argv[i])>3)
				{
					std::cerr<<"Unknown option: "<<argv[i]<<"\n";
					return false;
				}
				if (strlen(argv[i])>2)
					switch (argv[i][2])
					{
						case 'm' :
							break;
						case 'p':
							break;
						default:
							std::cerr<<"ignoring unknown weight type "<<argv[i][2]<<std::endl;
							break;
					}
				
				break;
			}
			default: {
				std::cerr<<"unknown option: "<<argv[i][1]<<std::endl;
				return false;
			}
		}
	}
	return params_good;
}

void print_help_message()
{
	std::cerr<<"Use fragrep as follows:\n\n";
	std::cerr<<"fragrep [options] <motif-file> [<genome-file>],\n\n";
	std::cerr<<"where <motif-file> is a filename as well as <genome-file>.\n";
	std::cerr<<"if <genome-file> is not specified, stdin is read instead, which is expected to be in fasta format.\n";
	std::cerr<<"valid options are:\n";
	std::cerr<<"-u (unweighted) to switch off weighting and optimize for intersection cardinality;\n";
	std::cerr<<"-w select weighting based on based on p-values\n";
	std::cerr<<"-v verbous output\n";
	std::cerr<<"-c scaffold cutoff scoring\n";
	std::cerr<<"-r include reverse complemented matches (default for fragrep-0.4 and later)\n";
	std::cerr<<"-f do not include reverse complemented matches\n";
	std::cerr<<"-s <prefixlength> <suffixlength> to extend extraction of result sequences by a specified number of nucleotides\n";
	std::cerr<<"-S use simplified output (leave out aligned query sequences and print reverse sequences for reverse matches)\n";
	std::cerr<<"   output produced with -S can be plugged into alignment programs such as clustalw immediately, but is less legible.\n";
	std::cerr<<"\n";
}

void get_next_line(std::ifstream& file, std::string& line)
{
	bool comment=false;
	do {
		getline(file,line);
		if (file.good() && line.size()>0)
			if (line[0]=='#')
				comment=true;
			else
				comment=false;
	} while (file.good() && comment);
}

bool read_data(int argc, char** argv)
{
	
	int j;
	params_good = true;

	if (argc-num_of_opts<2)
	{
		print_help_message();
		params_good = false;
	}
	
	// 0. READ NUMBER OF FRAGMENTS
	std::ifstream frag_file(argv[num_of_opts+1]);
	if (!frag_file.good())
	{
		std::cerr<<"could not open file "<<argv[num_of_opts+1]<<"\n";
		return false;
	}
	std::string read_buf;
	get_next_line(frag_file,read_buf);
	std::istringstream tokens(read_buf);
	tokens>>k;
	if (tokens.fail() || k<0)
	{
		std::cerr<<"expected a positive integer in pattern file.\n";
 		exit(-1);
	}
	// check whether conventional consensus patterns or position weight
	// matrices are used
	tokens>>read_buf;
	if (!tokens.fail() && read_buf=="matrices")
		use_pwms = true;

	// read delta (max. number of block deletions), if available
	tokens>>read_buf;
	if (!tokens.fail() && read_buf.size()>0)
	{
		if (read_buf=="-")
		{
			tokens>>delta;
			if (tokens.fail() || delta<0)
			{
				std::cerr<<"expected a positive integer in pattern file.\n";
				exit(-1);
			}
		} else if (read_buf[0]=='-') {
			std::istringstream numstr(read_buf.substr(1));
			numstr>>delta;
			if (numstr.fail() || delta<0)
			{
				std::cerr<<"expected a positive integer in pattern file.\n";
				exit(-1);
			}		
		}
	}
	
	s.resize(k);
	m.resize(k);
	d.resize(k);
	n.resize(k);
	u.resize(k);
	l.resize(k);
	t.resize(k);
	PWMs.resize(k);
	
	for (j=0; j<k && frag_file.good(); j++)
	{
		// read lower length bound
		get_next_line(frag_file,read_buf);
		if (read_buf.empty() || read_buf[0]=='#')
		{
			j--;
			continue;
		}
		std::istringstream line_tokens(read_buf);
		line_tokens>>l[j];
		if (line_tokens.fail())
		{
			std::cerr<<"expected an integer in pattern file.\n";
			exit(-1);
		}
		// read upper length bound
		line_tokens>>u[j];
		if (line_tokens.fail())
		{
			std::cerr<<"expected an integer in pattern file.\n";
			exit(-1);
		}		
		// read fragment sequence, currently delimited to 8192 characters
		line_tokens>>s[j];
		//cout<<"> "<<x.substr(0,x.find('$'))<<"\n";
		//cout<<"> "<<x.substr(x.find('$')+1,x.size())<<"\n";

		//line_tokens>>read_buf;
		//int s_len = read_buf.size();
		//s[j].resize(s_len+1);
		//strcpy(s[j],read_buf.c_str());
		n[j] = s[j].size();
		if (!use_pwms)
		{
			// read maximum edit distance for the fragment
			line_tokens>>m[j];
			if (line_tokens.fail())
			{
				std::cerr<<"expected an integer in pattern file.\n";
				exit(-1);
			}
		} else {
			// read maximum edit distance for the fragment
			line_tokens>>t[j];
			if (line_tokens.fail())
			{
				std::cerr<<"expected a float in pattern file.\n";
				exit(-1);
			}
			if (t[j]<=0. || t[j]>1.)
			{
				std::cerr<<"WARNING: fishy threshold for matrix matching.\n";
			}
		}
		line_tokens>>d[j];
		if (line_tokens.fail())
			d[j]=0;
	}
	if (j<k)
	{
		std::cerr<<"less patterns than specified.\n";
		exit(-1);
	}
		
	if (use_pwms)
	{
		for (j=0; j<k && frag_file.good(); j++)
		{
			get_next_line(frag_file,read_buf);
			if (read_buf!=s[j])
			{
				std::cerr<<"wrong matrix ID\n";
				exit(-1);
			}
			frag_file>>PWMs[j];
			PWMs[j].set_RNA();
			n[j] = PWMs[j].length;
			//std::cout<<"matrix "<<j<<":\n"<<PWMs[j]<<"\n\n";
		}
	}

	return true;
	
}

void revert_segments()
{

	if (!reverse)
		return;
	
	int i;
	
	std::string compl_table;
	compl_table.resize(256);

	for(i=0; i<256; i++) compl_table[i]=(int)'-';
	
	compl_table[(int)'A']='T';
	compl_table[(int)'C']='G';
	compl_table[(int)'G']='C';
	compl_table[(int)'T']='A';
	compl_table[(int)'R']='Y';
	compl_table[(int)'Y']='R';
	compl_table[(int)'M']='K';
	compl_table[(int)'K']='M';
	compl_table[(int)'S']='S';
	compl_table[(int)'W']='W';
	compl_table[(int)'H']='D';
	compl_table[(int)'B']='V';
	compl_table[(int)'V']='B';
	compl_table[(int)'D']='H';
	compl_table[(int)'N']='N';
	compl_table[(int)'a']='t';
	compl_table[(int)'c']='g';
	compl_table[(int)'g']='c';
	compl_table[(int)'t']='a';
	compl_table[(int)'r']='y';
	compl_table[(int)'y']='r';
	compl_table[(int)'m']='k';
	compl_table[(int)'k']='m';
	compl_table[(int)'s']='s';
	compl_table[(int)'w']='w';
	compl_table[(int)'h']='d';
	compl_table[(int)'b']='v';
	compl_table[(int)'v']='b';
	compl_table[(int)'d']='h';
	compl_table[(int)'n']='n';
	
	int j, iota;
}


void free_all()
{
}

int main(int argc, char** argv)
{
	
	T.clear();
	s.clear();
	n.clear();
	l.clear();
	u.clear();
	m.clear();
	d.clear();

	if (!read_parameters(argc,argv))
	{
		exit(0);
		return 0;
	}		
	if (!params_good)
	{
		print_help_message();
		exit(0);
		return 0;
	}
	if (!read_data(argc,argv))
	{
		exit(0);
		return 0;
	}
	if (use_pwms)
	{
		std::cerr<<k<<" matrices read.\n";
	}
	int jj;

	omit_pattern_matching pm;
	if (!use_pwms)
		pm.init(k, s, l, u, m);
	else
		pm.init(k, PWMs, l, u, t);
		
	pm.set_output_region(pref_length,suff_length);
	pm.set_simple_output(simple_output);
	pm.set_scaffold_scoring(scaffold_cutoff_scoring);
	pm.set_verbose(verbose);
	pm.set_unique(unique);
	pm.set_deletions(d);
	pm.set_delta(delta);
	// by default, fastafile reads from stdin.

	fasta_reader fr;
	if (argc-num_of_opts>2)
	{
		fr.open(argv[num_of_opts+2]);
		std::string xxh,xxs;
	}
	
	std::string header,seq;
	if (verbose)
		std::cerr<<"checking for next sequence...\n";

	while (fr.get_seq (header, seq)) {
		std::string T_name = header;
		if (verbose)
			std::cerr<<"searching sequence "<<header<<"\n";
		pm.set_genome(seq,header);
		pm.run(delta);
		if (reverse)
		{
			if (verbose)
				std::cerr<<"searching reverse complement...\n";
			pm.reverse_compl();
			pm.run(delta);
			// in case this was not the last genome sequence, we need to
			// restore the original status...
			pm.reverse_compl();
		}
		if (verbose)
			std::cerr<<"checking for next sequence...\n";
	} 
	if (verbose)
		std::cerr<<"done."<<header<<"\n";

	pm.clear();
	free_all();
	
	return 0;
	
}



