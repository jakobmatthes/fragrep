/*
  fragrep c++ rebuild */

// std and iostream headers
#include <cstdlib>
#include <cstring>
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
#include "eps_output.h"
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


float max_line_length;

bool use_pwms;

string_vector s;

/// n[j] is the length of s[j]
int_vector n;

/// N is the length of T
int N;

/// m[j] number of mismatches allowed for occurences of s[j]
int_vector m;

/// d[j] specifies the maximum number of deletions for occurences of s[j]
int_vector d;

/// t[j] threshold for pwm match score of s[j]
float_vector t;

PWM_vector PWMs;

PWM_list PWMs_l;

int_pair_list BOUNDS;

int_list DELS;

float_list SCORES;

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

void print_help_message()
{
	std::cerr<<"usage:\n\n";
	std::cerr<<" pattern2eps [-L <line_length>] <pattern-file>\n\n";
	std::cerr<<"where <pattern-file> is the name of a fragrep pattern.\n";
	std::cerr<<"<line_length> secifies the width of the output, the default being 100.\n";
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

	if (argc<2)
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
	use_pwms = false;
	if (!tokens.fail() && read_buf=="matrices")
		use_pwms = true;
	if (!use_pwms)
	{
		std::cerr<<"pattern2eps works on matrix patterns only\n";
		return false;
	}

	int delta;
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

/// parse command line options
// ==============================================================
bool read_parameters(int argc, char** argv)
// ==============================================================
{

	params_good = true;
	num_of_opts = 0;
	max_line_length=50.;
	
	int i;
	for (i=1; i<argc && argv[i][0]=='-'; i++)
	{
		num_of_opts++;
		switch (argv[i][1])
		{
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
			case 'h' : {
				print_help_message();
				exit(0);
			}
		}
	}
	return params_good;
}




int main(int argc, char** argv)
{
	
	s.clear();
	n.clear();
	l.clear();
	u.clear();
	m.clear();
	d.clear();

	max_line_length = 50.;
	
	if (!read_parameters(argc,argv))
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

	PWMs_l.clear();
	BOUNDS.clear();
	DELS.clear();
	SCORES.clear();
	
	int ii;
	for (ii=0; ii<k; ii++)
	{
		PWMs_l.push_back(PWMs[ii]);
		BOUNDS.push_back(int_pair(l[ii],u[ii]));
		DELS.push_back(d[ii]);
		SCORES.push_back(t[ii]);
	}
	write_eps_file(PWMs_l, BOUNDS, DELS, SCORES, max_line_length,"pattern2eps.eps");
	std::cerr<<"output written to pattern2eps.eps\n";
	
	return 0;
	
}



