#include <sequtil/alignment_reader.hpp>
#include "pattern_matching.h"
#include <sequtil/bbq_util.h>
#include <sequtil/position_weight_matrix.h>
#include <iostream>
#include <cstdio>
#include <list>
#include <vector>
#include <fstream>
#include <iomanip>
#include "math.h"
#include "eps_output.h"

typedef std::pair<int,int> fragrep_block;
typedef std::list<fragrep_block> fragrep_pattern;
typedef std::pair<int,int> int_pair;
typedef std::list<int_pair> int_pair_list;
typedef std::list<int> int_list;
typedef std::list<float> float_list;
typedef std::string::iterator str_iterator;
typedef std::string::const_iterator str_c_iterator;
typedef std::vector<float> float_vector;
typedef std::pair<float,char> float_char;
typedef std::vector<float_char> float_char_vector;


typedef bbq::position_weight_matrix PWM;
typedef std::list<PWM> PWM_list;

/// object for importing clustalW formatted alignments
alignment_reader AR;

/// object handling some core functionality of fragrep
fragrep::pattern_matching PM;


std::string block_line;

/// description of the positions of the conserved alignments within
/// the annotated mutliple alignment
fragrep_pattern all_blocks;

float sensitivity;

float left_offset;

/// counter for reading command line options
int num_of_opts;

/// flag for checking validty of command line options
bool params_good;

/// flag for producing matrix patterns (true) or IUPAC consensus
/// patterns (false)
bool make_matrices;

/// formatting information for eps output of matrix patterns
float rewind_length;

/// formatting information for eps output of matrix patterns
float max_line_length;

float GC_ratio;

/// determine maximal number of mismatches for one block of an IUPAC
/// consensus sequences
// ==============================================================
int compute_mismatches(const std::string& pattern, const std::string& text)
// ==============================================================
{

	int m_cost = 0;

	str_c_iterator p_it = pattern.begin();
	str_c_iterator p_end = pattern.end();
	str_c_iterator t_it = text.begin();
	str_c_iterator t_end = text.end();
	
	for (; p_it!=p_end && t_it!=t_end; ++t_it, ++p_it)
	{
		if (PM.match_table[PM.index_table[*p_it]][PM.index_table[*t_it]]==0)
			m_cost++;
	}
	return m_cost;
}


/// determine maximal number of deletions for one block of an IUPAC
/// consensus sequences
// ==============================================================
int compute_deletions(const std::string& pattern)
// ==============================================================
{

	int m_cost = 0;

	str_c_iterator p_it = pattern.begin();
	str_c_iterator p_end = pattern.end();
	for (; p_it!=p_end; ++p_it)
		if (*p_it=='-' || *p_it=='.')
			m_cost++;
	return m_cost;
}


/// helper for parsing the command line options
// ==============================================================
void read_long_parameter(int argc, char** argv, int& i)
// ==============================================================
{
	if (!strcmp(argv[i],"--sensitivity"))
	{ 
		if (argc<=i+1 || !is_flt(argv[i+1]))
		{
			throw bbq_exception("expected float specifying specificty.\n");
		}
		if (i+1<argc && argv[i+1][0]!='-')
		{
			// read num_of_group_characters
			sensitivity = atof(argv[i+1]);
			i++;
			num_of_opts++;
		}
	} else if (!strcmp(argv[i],"--matrices")) {
		make_matrices = true;
	} else { 
		throw bbq_exception(std::string("unknown option: ")+std::string(argv[i]));
	}
}

// ==============================================================
void print_help_message()
// ==============================================================
{
	std::cout<<"\nusage: aln2pattern [options] [filename]\n\n";
	std::cout<<"Creates an input pattern for fragrep from an annotated alignment.\n \
Reads alignments in ClustalW format annotated by \n \
an extra alignment sequence termed FRAGREP which consists of\n \
* (consider column as conserved) and - (non-conserved) characters only.\n \
aln2pattern writes to stdout and creates the file aln2pattern.eps \n \
containing a visualization of the pattern, if option -m is used.\n\n";
	std::cout<<"If filename is not specified, aln2pattern\nreads the alignment form stdin.\n\n";
	std::cout<<"options:\n";
	std::cout<<" -m  create matrix pattern\n";
	std::cout<<" -s <threshold>  specify threshold for consensus character determination\n    (ignored if -m is specified).\n";
	std::cout<<" -L <max_length>  specify maximum line length for eps output\n    (default 100, ignored if -m is not specified).\n\n";
}

/// parse command line options
// ==============================================================
void read_parameters(int argc, char** argv)
// ==============================================================
{
	params_good = true;
	make_matrices = false;
	sensitivity = .75;
	num_of_opts = 0;
	max_line_length=50.;
	int i;
	for (i=1; i<argc && argv[i][0]=='-'; i++)
	{
		num_of_opts++;
		switch (argv[i][1])
		{
			case '-' : {
				read_long_parameter(argc,argv,i);
				break;
			}
			case 's' : {
				if (strlen(argv[i])>3)
				{
					throw bbq_exception(std::string("Unknown option: ")+std::string(argv[i])+std::string("\n"));
				}
				if (i+1<argc)
				{
					if (argc<=i+1 || !is_flt(argv[i+1]))
					{
						throw bbq_exception("expected float specifying sensitivity.\n");
					}
					sensitivity = atof(argv[i+1]);
					i++;
					num_of_opts++;
				} else {
					params_good = false;
				}
				break;
			}
			case 'L' : {
				if (strlen(argv[i])>3)
				{
					throw bbq_exception(std::string("Unknown option: ")+std::string(argv[i])+std::string("\n"));
				}
				if (i+1<argc)
				{
					if (argc<=i+1 || !is_flt(argv[i+1]))
					{
						throw bbq_exception("expected float specifying line length scale factor\n (in percent, i.e., default is -L 100).\n");
					}
					max_line_length = atof(argv[i+1])/2.;
					i++;
					num_of_opts++;
				} else {
					params_good = false;
				}
				break;
			}
			case 'm' : {
				if (strlen(argv[i])>3)
				{
					throw bbq_exception(std::string("Unknown option: ")+std::string(argv[i])+std::string("\n"));
				}
				make_matrices = true;
				break;
			}
			case 'h' : {
				print_help_message();
				exit(0);
			}
		}
	}
}

/// parameter driven heuristic for constructing the IUPAC code of a
/// given distribution over the nucleotide alphabet
// ==============================================================
char consensus_char(int f_A, int f_C, int f_G, int f_U)
// ==============================================================
{
	float sum = (float)(f_A+f_C+f_G+f_U);
	float r_A = ((float)f_A)/sum;
	float r_C = ((float)f_C)/sum;
	float r_G = ((float)f_G)/sum;
	float r_U = ((float)f_U)/sum;
	if (r_A>=sensitivity)
		return 'A';
	if (r_C>=sensitivity)
		return 'C';
	if (r_G>=sensitivity)
		return 'G';
	if (r_U>=sensitivity)
		return 'U';
	if (r_A+r_C>=sensitivity)
		return 'M';
	if (r_A+r_G>=sensitivity)
		return 'R';
	if (r_A+r_U>=sensitivity)
		return 'W';
	if (r_C+r_G>=sensitivity)
		return 'S';
	if (r_C+r_U>=sensitivity)
		return 'Y';
	if (r_G+r_U>=sensitivity)
		return 'K';
	if (r_A+r_C+r_G>=sensitivity)
		return 'V';
	if (r_A+r_C+r_U>=sensitivity)
		return 'H';
	if (r_A+r_G+r_U>=sensitivity)
		return 'D';
	if (r_C+r_G+r_U>=sensitivity)
		return 'B';
	return 'N';
}

/// compute the IUPAC coded consensus for a conserved block in the
/// annotated alignment
// ==============================================================
std::string get_consensus_seq(int gap_start, int gap_end)
// ==============================================================
{
	std::string result;
	result.clear();
	int freq_A, freq_C, freq_G, freq_U, freq_gap;
	int gap_len = gap_end-gap_start;
	int max_gaps = 0;
	int min_gaps = gap_len;
	int i;
	for (i=gap_start; i<gap_end; i++)
	{
		freq_A = 0;
		freq_C = 0;
		freq_G = 0;
		freq_U = 0;
		freq_gap = 0;
		alignment_reader::iterator a_it = AR.begin();
		alignment_reader::iterator a_end = AR.end();
		for (; a_it!=a_end; ++a_it)
		{
			switch (a_it->second[i])
			{
				case 'A': freq_A++; break;
				case 'C': freq_C++; break;
				case 'G': freq_G++; break;
				case 'U': freq_U++; break;
				case 'T': freq_U++; break;
				case '#': break;
				case '-': freq_gap++;
				case '.': freq_gap++;
			}
		}
		result.push_back(consensus_char(freq_A,freq_C,freq_G,freq_U));
	}
	return result;
}


/// determine mismatch parameters for a pattern in the IUPAC world
// ==============================================================
int get_mismatches(const std::string& consensus_pattern, int gap_start, int gap_end)
// ==============================================================
{
	int max_mm = 0;
	int len = gap_end-gap_start;
	int i;
	for (i=gap_start; i<gap_end; i++)
	{
		alignment_reader::iterator a_it = AR.begin();
		alignment_reader::iterator a_end = AR.end();
		for (; a_it!=a_end; ++a_it)
		{
			int this_mm = compute_mismatches(consensus_pattern,a_it->second.substr(gap_start,len));
			if (this_mm>max_mm)
				max_mm=this_mm;
		}
		
	}
	return max_mm;	
}

// ==============================================================
int get_deletions(int gap_start, int gap_end)
// ==============================================================
{
	int max_del = 0;
	int len = gap_end-gap_start;
	int i;
	for (i=gap_start; i<gap_end; i++)
	{
		alignment_reader::iterator a_it = AR.begin();
		alignment_reader::iterator a_end = AR.end();
		for (; a_it!=a_end; ++a_it)
		{
			int this_del = compute_deletions(a_it->second.substr(gap_start,len));
			if (this_del>max_del)
				max_del=this_del;
		}
		
	}
	return max_del;	
}

/// compute upper and lower bounds for gap lengths in between two
/// conserved blocks, the first block ending at gap_end and the second
/// one ending at gap_end
// ==============================================================
int_pair get_bounds(int gap_start, int gap_end)
// ==============================================================
{
	int i;
	int gap_len = gap_end-gap_start;
	int max_gaps = 0;
	int min_gaps = gap_len;
	alignment_reader::iterator a_it = AR.begin();
	alignment_reader::iterator a_end = AR.end();
	for (; a_it!=a_end; ++a_it)
	{
		int gaps_i = 0;
		for (i=gap_start; i<gap_end; i++)
		{
			if (a_it->second[i]=='-' || a_it->second[i]=='.')
				gaps_i++;
		}
		if (gaps_i>max_gaps)
			max_gaps = gaps_i;
		if (gaps_i<min_gaps)
			min_gaps = gaps_i;
	}
	return int_pair(gap_len-max_gaps,gap_len-min_gaps);
}

/// compute background GC ratio
// ==============================================================
float get_GC_ratio()
// ==============================================================
{
	int i,j;
	int AT_count=0;
	int GC_count=0;

	alignment_reader::iterator a_it = AR.begin();
	alignment_reader::iterator a_end = AR.end();
	int len = AR.length();
	for (; a_it!=a_end; ++a_it)
	{
		std::string curr_line = a_it->second;
		int clen = curr_line.size();
		for (j=0; j<clen; j++)
			if (curr_line[j]=='A' || curr_line[j]=='T')
				AT_count++;
			else if (curr_line[j]=='C' || curr_line[j]=='G')
				GC_count++;
	}
	if (AT_count+GC_count==0)
		return 0.;
	return ((float)GC_count) / (float)(AT_count+GC_count);

}

// ==============================================================
int get_indel_bounds(int gap_start, int gap_end)
// ==============================================================
{
	int_pair bounds = get_bounds(gap_start,gap_end);
	int gap_len = gap_end-gap_start;
	return -(bounds.first-gap_len);
}

///convert position within the alignment to a position within one of
///the aligned sequences by 'gap counting'
// ==============================================================
int aln_to_seq_pos(int seq, int aln_pos)
// ==============================================================
{
	int seq_pos = aln_pos;
	int len = AR.length();
	if (aln_pos>=len)
		return len-1;
	int j;
	for (j=0; j<len; j++)
	{
		if (AR(seq,j)=='-' || AR(seq,j)=='.')
			seq_pos--;
	}
	return seq_pos;
}

/// compute anchors that can be used for a dialign alignment
// ==============================================================
void write_anchors(const fragrep_block& B, std::ofstream& anchor_file)
// ==============================================================
{
	int hght = AR.height();
	int ref_beg = aln_to_seq_pos(0,B.first);
	int ref_end = aln_to_seq_pos(0,B.second);
	int ref_len = B.second-B.first;
	int i;
	for (i=1; i<hght; i++)
	{
		int beg=aln_to_seq_pos(i,B.first);
		int end=aln_to_seq_pos(i,B.second);
		anchor_file<<"1 "<<i+1<<" "<<ref_beg<<" "<<beg<<" "<<ref_len<<" 1.0\n";
	}
}

// ==============================================================
void remove_gaps(std::string& str)
// ==============================================================
{
	int len = str.length();
	int i,ins_pos=0;
	int dels;
	for (i=0; i<len; ++i)
	{
		if (str[i]!='-' && str[i]!='.') {
			str[ins_pos]=str[i];
			ins_pos++;
		}
	}
	str.resize(ins_pos);
}

/// write a IUPAC coded consensus search pattern for fragrep to stdout
// ==============================================================
void write_pattern()
// ==============================================================
{
	int hght = AR.height();
	if (all_blocks.empty())
		return;
	std::cout<<all_blocks.size()<<"\n";
	int freq_A, freq_C, freq_G, freq_U;
	fragrep_pattern::iterator f_it = all_blocks.begin();
	fragrep_pattern::iterator f_end = all_blocks.end();
	int prev_block_end = f_it->first;
	for (; f_it!=f_end; ++f_it)
	{
		int_pair bounds = get_bounds(prev_block_end,f_it->first);
		std::string consensus = get_consensus_seq(f_it->first,f_it->second);
		int mismatches = get_mismatches(consensus,f_it->first,f_it->second);
		int deletions = get_deletions(f_it->first,f_it->second);	
		int indel_bounds = get_indel_bounds(prev_block_end,f_it->first);
		std::cout
			<<bounds.first<<" "<<bounds.second<<" "<<consensus
			<<" "<<mismatches<<" "<<deletions<<"\n";
		prev_block_end = f_it->second;
	}
}

/// compute a position frequency matrix form the alignment positions
/// between gap_start and gap_end, returning a [0,1]-valued
/// match_score and a maximum number of deletions in dels.
// ==============================================================
PWM get_PWM(int gap_start, int gap_end, float& match_score, int& dels)
// ==============================================================
{
	int freq_A, freq_C, freq_G, freq_U, freq_gap;
	float_vector A,C,G,T;
	int gap_len = gap_end-gap_start;
	int max_gaps = 0;
	int min_gaps = gap_len;
	dels = get_deletions(gap_start,gap_end);
	int i;
	for (i=gap_start; i<gap_end; i++)
	{
		freq_A = 0;
		freq_C = 0;
		freq_G = 0;
		freq_U = 0;
		freq_gap = 0;
		alignment_reader::iterator a_it = AR.begin();
		alignment_reader::iterator a_end = AR.end();
		int dels=0;
		for (; a_it!=a_end; ++a_it)
		{
			char ch=a_it->second[i];
			switch (a_it->second[i])
			{
				case 'A': freq_A++; break;
				case 'C': freq_C++; break;
				case 'G': freq_G++; break;
				case 'U': freq_U++; break;
				case 'T': freq_U++; break;
				case '#': break;
				case '-' : break;
				case '.' : break;
				default : std::cerr<<"UNKNOWN CHAR: "<<a_it->second[i]<<"\n";
			}
		}
		A.push_back((float)freq_A);
		C.push_back((float)freq_C);
		G.push_back((float)freq_G);
		T.push_back((float)freq_U);
		
		//result.push_back(consensus_char(freq_A,freq_C,freq_G,freq_U));
	}
	PWM result(A,C,G,T);
	match_score=1.;
	alignment_reader::iterator a_it = AR.begin();
	alignment_reader::iterator a_end = AR.end();
	for (; a_it!=a_end; ++a_it)
	{
		typedef std::pair<float,float> float_pair;
		std::string s_cstr = a_it->second.substr(gap_start,gap_end-gap_start);
		remove_gaps(s_cstr);
		result.set_indel_costs(1.,1.,dels);
		float_pair m_score = result.get_frac_score(s_cstr.c_str());
		if (m_score.first<match_score)
			match_score=m_score.first;
	}
	match_score = floor(100.*match_score)/100.;
	return result;
}

/// 
// ==============================================================
void write_matrix_pattern()
// ==============================================================
{
	std::string anchor_filename = "aln2pattern.anchor";
	std::cerr<<"writing anchor file...\n";
	std::ofstream anchor_file(anchor_filename.c_str(),std::ios::out);
	AR.write_fasta("aln2pattern.mfa");
	std::cerr<<"writing fasta file...\n";

	PWM_list PWMs;
	int_pair_list BOUNDS;
	int_list DELS;
	float_list SCORES;
	
	int i;
	float m_score;
	int deletions;
	int hght = AR.height();
	if (all_blocks.empty())
		return;
	std::cout<<all_blocks.size()<<" matrices\n";
	fragrep_pattern::iterator f_it = all_blocks.begin();
	fragrep_pattern::iterator f_end = all_blocks.end();
	int prev_block_end = f_it->first;
	for (i=0; f_it!=f_end; ++i, ++f_it)
	{
		write_anchors(*f_it,anchor_file);
		int_pair bounds = get_bounds(prev_block_end,f_it->first);
		PWMs.push_back(get_PWM(f_it->first,f_it->second,m_score,deletions));
		BOUNDS.push_back(bounds);
		DELS.push_back(deletions);
		SCORES.push_back(m_score);
		//int deletions = get_deletions(f_it->first,f_it->second);	
		//int mismatches = get_mismatches(consensus,f_it->first,f_it->second);
		//int indel_bounds = get_indel_bounds(prev_block_end,f_it->first);
		PWMs.back().set_RNA();
		std::cout
			<<std::setw(5)<<bounds.first<<" "
			<<std::setw(5)<<bounds.second<<" M"<<i
			<<":"<<(PWMs.back()).most_informative_string(GC_ratio)
			<<" "<<m_score<<" "<<deletions<<"\n";
		prev_block_end = f_it->second;
	}
	
	fragrep::write_eps_file(PWMs,BOUNDS,DELS,SCORES,max_line_length,"aln2pattern.eps");
	
	PWM_list::iterator m_it = PWMs.begin();
	PWM_list::iterator m_end = PWMs.end();
	for (i=0; m_it!=m_end; ++i, ++m_it)
	{
		std::cout
			<<"M"<<i
			<<":"<<m_it->most_informative_string(GC_ratio)
			<<"\n"<<*m_it;
		int len=m_it->length;
		int j;
		std::cout<<"#  ";
		for (j=0; j<len; j++)
		{
			std::cout<<m_it->most_informative_character(j,GC_ratio)<<"    ";
		}
		std::cout<<"\n";
	}
	anchor_file.close();
	
 	m_it = PWMs.begin();
	m_end = PWMs.end();
	int_pair_list::iterator B_it = BOUNDS.begin();
	int_pair_list::iterator B_end = BOUNDS.end();
	int_list::iterator D_it = DELS.begin();
	int_list::iterator D_end = DELS.end();
	float_list::iterator S_it = SCORES.begin();
	float_list::iterator S_end = SCORES.end();

	m_it = PWMs.begin();
	m_end = PWMs.end();
	int num_of_lines = 1;
	float line_length=5.5;
	float longest_line = 0.;
	bool line_break=false;
	for (i=0; m_it!=m_end; ++i, ++m_it)
	{
		if (line_break)
		{
			num_of_lines++;
			line_break = false;
		}
 		line_length += 1.5*(float)m_it->length+6.;			
		if (line_length>longest_line)
		{
			longest_line=line_length;
		}
		if (line_length>max_line_length)
		{
			line_length=5.5;
			line_break=true;
		}
	}
	
	m_it = PWMs.begin();
	m_end = PWMs.end();
	B_it = BOUNDS.begin();
	B_end = BOUNDS.end();
	D_it = DELS.begin();
	D_end = DELS.end();
	S_it = SCORES.begin();
	S_end = SCORES.end();	

}


// ==============================================================
void scan_block_line()
// ==============================================================
{
	std::string::iterator b_it = block_line.begin();
	std::string::iterator b_end = block_line.end();
	int block_start=0, block_end=0;
	char curr_block_id = 0;
	fragrep_block curr_block;
	
	for (; b_it!=b_end; ++b_it, ++block_end)
	{
		int new_block_id;
		if (*b_it==' ' || *b_it=='-' || *b_it=='.')
			new_block_id = 0;
		else
			new_block_id = *b_it;
		if (new_block_id!=curr_block_id)
		{
			if (curr_block_id!=0)
			{
				all_blocks.push_back(fragrep_block(block_start,block_end));
			}
			block_start = block_end;
			curr_block_id = new_block_id;
		}
	}
	if (curr_block_id!=0)
		all_blocks.push_back(fragrep_block(block_start,block_end));
}

// ==============================================================
void clean_sequences()
// ==============================================================
{
	int len = AR.length();
	int hght = AR.height();
	int i,j;
	for (i=0; i<hght; i++)
		for (j=0; j<len; j++)
		{
			switch (AR(i,j))
			{
				case 'A': case 'C': case 'G': case 'U': case 'T': case '-': case '.': break;
				case 'a': AR(i,j)='A'; break;
				case 'c': AR(i,j)='C'; break;
				case 'g': AR(i,j)='G'; break;
				case 't': AR(i,j)='U'; break;
				case 'u': AR(i,j)='U'; break;
				default: std::cerr<<"replacing "<<AR(i,j)<<"\n"; AR(i,j)='#'; 
			}
		}
}

// ==============================================================
int main(int argc, char** argv)
// ==============================================================
{
	
	read_parameters(argc,argv);
	if (argc-1>num_of_opts)
		AR = alignment_reader(argv[num_of_opts+1]);
	else
		AR = alignment_reader(stdin);

	alignment_reader::iterator a_it = AR.begin();
	alignment_reader::iterator a_fragrep = AR.end();
	alignment_reader::iterator a_end = AR.end();
	for (; a_it!=a_end; ++a_it)
	{
		if (a_it->first=="FRAGREP" || a_it->first=="#=FRAGREP")
		{
			block_line = a_it->second;
			a_fragrep = a_it;
		}
	}
	if (a_fragrep==a_end)
	{
		std::cout<<"No block defining sequence found.\n";
		return 0;
	}

	block_line = a_fragrep->second;
	AR.erase(a_fragrep);

	clean_sequences();
	GC_ratio = get_GC_ratio();
	scan_block_line();
	if (!make_matrices)
		write_pattern();
	else {
		write_matrix_pattern();
		//write_eps_output();
	}
	
}

