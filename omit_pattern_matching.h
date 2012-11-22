// bbq - Copyright (C) 2005, Axel Mosig, University of Leipzig. See
// main.cpp for details.

#ifndef OMIT_PATTERN_MATCHING_H
#define OMIT_PATTERN_MATCHING_H

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

class weighted_match;
typedef std::list<weighted_match> w_match_list;
typedef w_match_list::iterator w_match_ref;
typedef std::set<w_match_ref> w_match_ref_set;
typedef std::pair<w_match_ref,w_match_ref> w_match_interval;
typedef std::vector<w_match_interval> w_match_interval_vector;
typedef std::vector<w_match_list> w_match_list_vector;
	

bool operator<(const w_match_ref& x, const w_match_ref& y);
	
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
class result_item
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
{

public:

	int pos;

	int j;

	double p_val;

	double E_val;

	std::string match_seq;

	result_item(int p, int j_, double pv, double ev, const std::string& seq);

};
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// A \c weighted_match represents one occurence of a particular
/// fragment (typically represented by a consensus sequence or a
/// position weight matrix) in a (genome) sequence; Since these
/// occurences are partially ordered, a \c weighted_match needs to
/// keep track of predecessors and successors in this partial order.
class weighted_match
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
{

public:
	
	weighted_match(int,int,float,const std::string&);

	weighted_match();

	/// fragment index
	int j;

	/// position in the (genome) sequence
	int pos;

	/// match score or -log-p-value of the PWM match of matrix j at
	/// position pos
	float weight;
	
	/// counter for dynamic programming 
	int deletions;

	/// keep track of summed up weights along the optimal DP path
	float score;

	/// matching string, evtl. after some deletions...
	std::string match_seq;
	
	/// reference to optimal predecessor in the dynamic programming
	/// recurrence
	w_match_ref DP_pred;

	/// reference to optimal successor in the dynamic programming
	/// recurrence
	w_match_ref DP_succ;

	/// flag for dynamic programming to avoid multiple reports of
	/// overlapping matches
	bool visited;

	/// set of all predecessors in the partial order on all
	/// weighted_matches in an instance
	w_match_ref_set pred;

	/// set of all successors in the partial order on all
	/// weighted_matches in an instance
	w_match_ref_set succ;

	/// lexicographic order w.r.t. j and pos
	bool operator<(const weighted_match&);
	
	/// two weighted_matches are equal iff their indices j and
	/// their position are equal
	bool operator==(const weighted_match&);
	
};
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
class match_pair
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
{
public:
	
	/// position of th match
	int first;

	/// flag for DP
	bool second;

	/// matching sequence, eventually involving gaps
	std::string match_seq;

	match_pair(int f, const std::string& seq);
	
};
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
class omit_pattern_matching
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
{
	
public:


	typedef std::list<match_pair> match_list;
	typedef std::list<result_item> result_list;
	typedef std::pair<float, result_list> complete_result;
	typedef std::map<int, complete_result> complete_result_pos_map;
	typedef std::vector<int> int_vector;
	typedef std::vector<int_vector> int_array;
	typedef std::vector<float> float_vector;
	typedef std::vector<match_list> match_list_vector;
	typedef std::vector<std::string> string_vector;
	typedef std::vector<int_vector> int_vector_vector;
	typedef std::list<int> int_list;
	typedef std::pair<float,int> float_int;
	typedef std::vector<float_int> float_int_vector;

	omit_pattern_matching();
	
	void init(int k_, const string_vector& s_, int_vector l_, int_vector u_, int_vector m_);

	void init(int k_, const PWM_vector& PWMs_, int_vector l_, int_vector u_, float_vector t_);

	void init(const std::string& T_, const std::string& T_name_, int k_, const string_vector& s_, int_vector l_, int_vector u_, int_vector m_);

	~omit_pattern_matching();
	
	void init();
	
	/// run the matching process with no omissions allowed.
	void run();

	/// run the matching process with up to delta_ omissions allowed.
	void run(int delta_);

	void clear();

	/// initialize omit_u_bound and omit_l_bound
	void compute_omit_bounds();

	/// initialize omit_u, and omit_l
	void compute_omit_ul();

	void set_genome(const std::string& T_, const std::string& T_name_);

	void set_output_region(int pref_len, int suff_len);

	void set_reverse(bool rev);

	void set_unique(bool unq);

	void set_delta(int delta_);

	/// set costs for fractional programming when matching PWMs
	/// with insertions and deletions
	void set_indel_costs(float ins, float del);
	
	void set_deletions(int_vector d_);

	void set_verbose(bool verb);

	void set_simple_output(bool so);

	void set_scaffold_scoring(bool sc);
	
	/// compute the reverse complement of the patterns, so that in effect
	/// the genome is searched in backward direction for occurences.
	void reverse_compl();

	void reverse_compl_old();

	void reverse_compl(char* S, int len_S);

	/// for matching wild card characters
	string_vector match_table;

	/// for converting C,G,T,A to 0,1,2,3
	int_vector index_table;
	
private:

	/// \c genome sequence
	std::string T;

	/// name associated with T
	std::string T_name;
	
	/// \c k sequence fragments
	string_vector s;

	PWM_vector PWMs;
	
	/// index of the most informative pattern, i.e., the one with the
	/// lowest E-value.
	int j_inf;

	/// for matches with omissions: pointer to the position up to where
	/// occurences of the most informative patterns have been searched
	/// in \c get_next_j_inf_occurence()
	int j_inf_start;

	/// in case delta many blocks can be deleted, we need to keep track
	/// of the (delta+1) most informative blocks, which we index in
	/// \c j_inf_list.
	int_list j_inf_list;
		
	/// n[j] is the length of s[j]
	int_vector n;
	
	/// N is the length of T
	int N;
	
	/// m[j] number of mismatches allowed for occurences of s[j]
	int_vector m;
 
	/// d[j] specifies the maximum number of deletions for occurences of s[j]
	int_vector d;

	/// deletion_penalties[j] specifies the cost for deleting block s[j]
	float_vector deletion_penalties;
	
	/// parameter for fractional programming
	float INS_COST;

	/// parameter for fractional programming
	float DEL_COST;
	
	/// t[j] is the match threshold for PWM[j]
	float_vector t;

	/// u[j] is an upper bound for the gap length between s[j-1] and s[j].
	int_vector u;
	
	/// l[j] is a lower bound for the gap length between s[j-1] and s[j].
	int_vector l;
	
	/// omit_u[b][j] is an upper bound for the gap length between s[b] and s[j].
	int_array omit_u;
	
	/// omit_l[b][j] is a lower bound for the gap length between s[b] and s[j].
	int_array omit_l;

	/// max width of a match
	int max_width;
	
	/// number of binding site fragments
	int k;

	/// done[j] denotes the pos. up to which T has been searched for s[j].
	int_vector done;

	/// pre-compute how far to look for patterns relative to an
	/// occurence of the most significant pattern
	int_vector l_bound;
	int_vector u_bound;

	/// pre-compute how far to look for patterns relative to an
	/// occurence of one of the relevant the most significant patterns
	int_vector_vector omit_l_bound;
	int_vector_vector omit_u_bound;
	
	/// table for computing reverse complemements
	std::string compl_table;
		
	/// match_lists[j] contains all curently "active" matches of s[j]
	/// (where only such matches are active that are relevant for
	/// occurences of the current match of s[j_inf].)
	match_list_vector match_lists;

	/// represents all matches that are currently under consideration
	/// in \c check_omit_match_lists().
	w_match_list all_w_matches;

	/// w_match_lists[j] represents the interval in all_w_matches
	/// that (exclusively) contains all w_matches X with X.j=j.
	w_match_list_vector w_match_lists;
	
	/// frequencies[i][x] contains the absolute frequency of nucleotide
	/// x in T[i], for x=C,G,T,A.
	int_vector frequencies;

	/// frequencies[i][x][y] contains the absolute frequency of the
	/// dinucleotide sequence xy in T[i], for x,y=C,G,T,A.
	int_vector_vector di_frequencies;

	/// a representation of the last match of a fragment s[j] or PWMs[j]
	string_vector last_match_sequences;
	
	int number_of_matches;
	
	/// search in both directions?
	bool match_reverse;

	/// report only the best match for each position, in case there is
	/// more than one.
	bool report_unique_matches;

	/// allow deletion of leading blocks at zero cost, suitable for
	/// 'cut off' occurences in scaffolds of pre-assemblies
	bool scaffold_cutoff_scoring;
	
	/// fragments already reverse-complemented?
	bool reverted;

	/// maximum number of blocks that can be omitted
	int delta;
	
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

	/// memo pointer for unique matches
	int curr_first_pos;

	/// memo pointer for unique matches
	int curr_last_pos;

	/// memo list for unique matches
	result_list curr_first_match;

	/// memo list for unique matches
	result_list curr_last_match;

	/// memo for avoiding redundant match enumeration
	int last_reported_match_pos;
	
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

	bool update_omit_match_lists(int pos, int curr_j_inf, const weighted_match& j_inf_match);

	int update_lower_match_lists(int pos, int curr_j_inf);

	int update_upper_match_lists(int pos, int curr_j_inf);

	void update_omit_bounds(int pos, int curr_j_inf);

	/// update partial order graph
	void add_to_partial_order(const w_match_list::iterator& match_ref);

	/// update partial order graph
	void remove_from_partial_order(const w_match_list::iterator& match_ref);
	
	/// perform dynamic programming in the match lists produced by the
	/// most recent call of update_match_lists(). Returns true if at
	/// least one match can be reconstructed that needs to be printed,
	/// false otherwise.
	bool check_match_lists();

	/// version of check_match_lists for block omission version
	bool check_omit_match_lists();
		
	/// update the index bounds indicating the regions where each s[j]
	/// needs to be searched. This procedure is called whenever a new
	/// occurence of s[j_inf] was found.
	void update_bounds(int pos);

	/// determine whether there is a match of s[j] at T[ii]
	bool match_consensus_fragment(int j, int ii);

	/// determine whether there is a match of fragment j at T[ii]	
	bool match_fragment(int j, int ii);

	/// determine whether there is a match of fragment j at T[ii] for
	/// the version with omissions
	bool match_fragment(int j, int ii, weighted_match& w_match);
	
	/// perform the complete matching task after all initialization has
	/// been performed in run().
	void match();

	void omit_match();

	void output_omit_matches();

	bool trace_back(w_match_ref tail, result_list& result);
	
	/// get the position of next occurence if s[j] in T by searching
	/// the region starting with T[done[j]], but not exceeding
	/// T[up_to].
	int get_next_occurence(int j, int up_to);

	/// version of next occurence finding for block omission version
	int get_next_occurence(int j, int up_to, weighted_match& w_match);

	int get_j_inf_occurence(int& j_inf, weighted_match& w_match);
	
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

	void get_j_inf_list();
	
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
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++



} // end namespace fragrep


#endif


