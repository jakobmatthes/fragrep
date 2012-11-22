// bbq - Copyright (C) 2005, Axel Mosig, University of Leipzig. See
// main.cpp for details.

#ifndef PATTERN_MATCHING_H
#define PATTERN_MATCHING_H

#include <set>
#include <string>
#include <map>
#include <list>
#include <vector>

#include <sequtil/position_weight_matrix.h>

#include <sequtil/triple.h>


//typedef uint_bitset bit_set;
//typedef uint64_bitset bit_set;
//typedef weighted_variable_bitset w_bit_set;



/** \c fragrep is the namespace that contains all classes except for those
	 used for reading parameters etc. in main.cpp */
namespace fragrep {

	typedef bbq::position_weight_matrix PWM;
	typedef std::vector<PWM> PWM_vector;

class result_item
{

public:

	int pos;

	double p_val;

	double E_val;

	result_item(int p, double pv, double ev);

};


	
class pattern_matching
{
	
public:

	typedef std::pair<int,bool> match_pair;
	typedef std::list<match_pair> match_list;
	typedef std::list<result_item> result_list;
	typedef std::pair<float, result_list> complete_result;
	typedef std::map<int, complete_result> complete_result_pos_map;

	pattern_matching();
	
	void init(int k_, char** s_, int* n_, int* l_, int* u_, int* m_);

	void init(int k_, const PWM_vector& PWMs_, int* n_, int* l_, int* u_, float* t_);

	void init(char* T_, const std::string& T_name_, int N_, int k_, char** s_, int* n_, int* l_, int* u_, int* m_);

	~pattern_matching();
	
	void init();

	void run();

	void clear();

	void set_genome(const char* T_, const std::string& T_name_, int N_);

	void set_output_region(int pref_len, int suff_len);

	void set_reverse(bool rev);

	void set_unique(bool unq);

	void set_deletions(int* d_);

	void set_verbose(bool verb);

	void set_simple_output(bool so);

	/// compute the reverse complement of the patterns, so that in effect
	/// the genome is searched in backward direction for occurences.
	void reverse_compl();

	void reverse_compl_old();

	void reverse_compl(char* S, int len_S);

	/// for matching wild card characters
	unsigned char** match_table;

	/// for converting C,G,T,A to 0,1,2,3
	int* index_table;
	
private:

	/// \c genome sequence
	char* T;

	/// name associated with T
	std::string T_name;
	
	/// \c k sequence fragments
	char** s;

	PWM_vector PWMs;
	
	/// index of the most informative pattern, i.e., the one with the
	/// lowest E-value.
	int j_inf;
		
	/// n[j] is the length of s[j]
	int* n;
	
	/// N is the length of T
	int N;
	
	/// m[j] number of mismatches allowed for occurences of s[j]
	int* m;
 
	/// d[j] specifies the maximum number of deletions for occurences of s[j]
	int* d;
	
	/// t[j] is the match threshold for PWM[j]
	float* t;

	/// u[j] is an upper bound for the gap length between s[j-1] and s[j].
	int* u;
	
	/// l[j] is a lower bound for the gap length between s[j-1] and s[j].
	int* l;
	
	/// number of binding site fragments
	int k;

	/// done[j] denotes the pos. up to which T has been searched for s[j].
	int* done;

	int* l_bound;

	int* u_bound;

	/// table for computing reverse complemements
	char* compl_table;
		
	/// match_lists[j] contains all curently "active" matches of s[j]
	/// (where only such matches are active that are relevant for
	/// occurences of the current match of s[j_inf].)
	match_list* match_lists;
	
	/// frequencies[i][x] contains the absolute frequency of nucleotide
	/// x in T[i], for x=C,G,T,A.
	int* frequencies;

	/// frequencies[i][x][y] contains the absolute frequency of the
	/// dinucleotide sequence xy in T[i], for x,y=C,G,T,A.
	int** di_frequencies;

	int number_of_matches;
	
	/// search in both directions?
	bool match_reverse;

	/// report only the best match for each position, in case there is
	/// more than one.
	bool report_unique_matches;
	
	/// fragments already reverse-complemented?
	bool reverted;

	/// set type of output to stderr
	bool verbose;
	
	/// leave out printing the aligned query sequences below the matches?
	bool simple_output;

	/// trigger use of position weight matrix based matching rather
	/// than consensus pattern based procedure.
	bool use_pwms;
	
	/// genome sequence prefix length to be appended before each matching region
	int prefix_length;

	/// genome sequence suffix length to be appended after each matching region	
	int suffix_length;

	/// in case uniqe results option is set, this removes all but the
	/// best match at each position from the list.
	void print_unique_results();

	/// collection of all (unique) results along with their E-value
	complete_result_pos_map all_results;
	
	/// in case uniqe results option is set, this method adds result to
	/// the collection of all results in case there is not a better solution
	/// at the same position as result.
	void add_unique_result(result_list& result);
	
	/// prints the subsequence of T corresponding to the fragment
	/// occs. contained in the list result, headed by a fasta header
	/// and with the originally specified fragments aligned below.
	void print_result(result_list& result);
	
	/// search new occurences in T of each s[j], limited to the regions
	/// that correspond to the occurence of s[j_inf] at pos. Returns
	/// true if for each s[j], at least one match was found in the
	/// admissible region, false otherwise.
	bool update_match_lists(int pos);

	/// perform dynamic programming in the match lists produced by the
	/// most recent call of update_match_lists(). Returns true if at
	/// least one match can be reconstructed that needs to be printed,
	/// false otherwise.
	bool check_match_lists();

	/// update the index bounds indicating the regions where each s[j]
	/// needs to be searched. This procedure is called whenever a new
	/// occurence of s[j_inf] was found.
	void update_bounds(int pos);

	/// determine whether there is a match of PWMs[j] at T[ii]
	bool match_PWM_fragment(int j, int ii);

	/// determine whether there is a match of s[j] at T[ii]
	bool match_consensus_fragment(int j, int ii);

	/// determine whether there is a match of fragment j at T[ii]	
	bool match_fragment(int j, int ii);

	/// perform the complete matching task after all initialization has
	/// been performed in run().
	void match();

	/// get the position of next occurence if s[j] in T by searching
	/// the region starting with T[done[j]], but not exceeding
	/// T[up_to].
	int get_next_occurence(int j, int up_to);

	/// recursive output function. pushes (if available) one occurence
	/// of s[j] obtained from the current match_lists in front of
	/// result and proceeds recursively with j-1. If j=-1, the list
	/// result is printed to stdout.
	void output_matches(int j, result_list& result);

	/// output all matches corresponding to one occurence of the most
	/// informative fragment s[y_inf].
	void output_matches();	

	/// obtain the index j_inf of the most informative fragment
	/// s[j_inf], i.e., the fragment with the smalles E-value.
	int get_j_inf();

	/// determine E-value of s[j] occuring at position pos in T; if
	/// pos=-1, determine the E-value of an exact occurence of s[j] in
	/// T.
	double E_value(int j, int pos);

	/// determine p-value of s[j] occuring at position pos in T; if
	/// pos=-1, determine the p-value of an exact occurence of s[j] in
	/// T.
	double p_value(int j, int pos);

	/// compute all nucleotide- and dinucleotide-frequencies and fill
	/// in the corresponding tablkes (needed for p- and E-value
	/// computation).
	void get_freq();

	/// Help procedure for extending dinucleotide freqs. to ambiguity
	/// codes.
	void sum_up_frq(char O, char P, char Q);
	
	/// initialization of the amiguity-code-tables needed for string
	/// matching.
	void compute_index_table();
	
	/// initialization of the amiguity-code-tables needed for string
	/// matching.
	void compute_match_table();

};



} // end namespace fragrep


#endif


