#ifndef POSITION_WEIGHT_MATRIX_H
#define POSITION_WEIGHT_MATRIX_H


#include <fstream>
#include <list>
#include <vector>
#include <iostream>
#include <utility>
#include "triple.h"

namespace bbq {

typedef std::pair<float,float> float_pair;
	
class position_weight_matrix
	{
		
	public:

		/// first indicates the position, second determines the match
		/// score, while third indicates the direction of the match
		/// (false iff reverse match)

		typedef std::triple<int,float,bool> match;
		typedef std::list< match > match_list;		
		typedef match_list::iterator match_iterator;

		typedef std::vector<int> int_vector;
		typedef std::vector<char> char_vector;
		typedef std::vector<float> float_vector;
		typedef std::vector<double> double_vector;
		
		typedef std::vector<float_vector> float_array;
		typedef std::vector<float_array> float_table;
		
		typedef std::pair<char,char> char_pair;
		typedef std::pair<int,int> int_pair;
		typedef std::list<char_pair> alignment_sequence;
		typedef std::list<int_pair> alignment_path;
		
		static int match_mode;

		/// true iff both forward matches and reverse complement matches
		/// are to be computed
		bool match_reverse_complement;

		/// if truem, \c most_informative_character() outputs 'T',
		/// otherwise 'U' for RNA sequences
		bool DNA;
		
		float_vector A;
		float_vector C;
		float_vector G;
		float_vector T;

		float_vector min_freq;
		float_vector max_freq;
		float_vector information_content;
		float_array freq;

		float Max;
		float Min;
		float_vector Rest_Bound;

		float_vector rev_min_freq;
		float_vector rev_max_freq;
		float_vector rev_information_content;
		float_array rev_freq;
		
		char_vector compl_table;
		
		int length;

		float max_score;
		float max_score_rev;
		
		int max_pos;

		/// insertion costs for fractional programming 
		float INS_COST;

		/// deletion costs for fractional programming 
		float DEL_COST;

		/// maximum number of deletions for fractional programming
		int max_dels;
		
		/// number of iterations for fractional programming
		int iterations;

		/// alignment information for fractional programming matches w/
		/// insertions and deletions
		alignment_sequence curr_aln_seq;

		/// aligned (sub-)sequence of the matrix following the last DP
		/// trace path
		std::string curr_match_seq;

		/// pre-computed most informative sequence w/o any deletions
		std::string mis_sequence;

		/// pre-computed reverse complemented most informative sequence w/o any deletions		
		std::string rev_mis_sequence;
		
		/// this feature is not yet supported...
		alignment_path curr_aln_path;

		position_weight_matrix();
		
		position_weight_matrix(const float_vector& A_, const float_vector& C_, const float_vector& G_, const float_vector& T_);
		
		//position_weight_matrix(const position_weight_matrix&);

		//position_weight_matrix& operator=(const position_weight_matrix&);

		friend std::ifstream& operator>>(std::ifstream& ifs, position_weight_matrix& M);

		match_list matches;
		
		bool read(const char* fname);
		
		bool get_plain_match_score(const char* seq, float theta);

		float_pair get_plain_match_score(const char* seq);
		
		float_pair get_match_score(const char* seq);
		
		float get_fwd_match_score(const char* seq);
		
		float get_rev_match_score(const char* seq);

		const std::string& get_match_sequence() const; 
		
		float get_fb_match_score(
			const char* seq,
			const float_array& fb_freq,
			const float_vector& fb_ic,
			const float_vector& fb_min_freq,
			const float_vector& fb_max_freq
			);
		
		float_pair get_frac_score(const char* seq);

		bool match_frac_score(const char* seq, float threshold);

		float get_fwd_frac_score(const char* seq);

		float get_rev_frac_score(const char* seq);

		std::string most_informative_string(float GC_ratio) const;
		
		char most_informative_character(int i, float GC_ratio) const;
		
		std::string rev_most_informative_string(float GC_ratio) const;
		
		char rev_most_informative_character(int i, float GC_ratio) const;
		
		bool fractional_programming_score(
			float theta,
			const char* seq,
			const float_array& fb_freq,
			const float_vector& fb_ic,
			const float_vector& fb_min_freq,
			const float_vector& fb_max_freq,
			bool trace_back
			);

		void trace_back(
			float theta,
			const char* seq,
			const float_array& fb_freq,
			const float_vector& fb_ic,
			const float_vector& fb_min_freq,
			const float_vector& fb_max_freq,
			int best_d,
			const float_table& FTABLE
			);

		float get_information_content();

		float get_p_value(double observed_score, const double_vector& bg_freq);
		
		void set_matches(const char* seq, float threshold);

		void set_MATCH_matches(const char* seq, float threshold);

		void set_naive_matches(const char* seq, float threshold);

		/// make \c most_informative_character() output ACGT
		void set_DNA();

		/// make \c most_informative_character() output ACGU
		void set_RNA();

		/// set cost parameters for fractional programming.
		/// Insertions and deletions are termed
		/// 'with respect to the sequence' (not 'with respect to the matrix')
		void set_indel_costs(float ins, float del, int m_d);
		
		/// if true is passed to this function, matches will be reported
		/// in BOTH, fwd and bwd direction. If false is passed, matches
		/// will be reported in fwd direction only.
		void set_reverse_complement(bool);
		
		void compute_information_content();

		void clear();

		void init();

		void reinit();

		match_list get_match_list();

		match_iterator begin();
		match_iterator end();
		
	};


std::ifstream& operator>>(std::ifstream& ifs, position_weight_matrix& M);

std::ostream& operator<<(std::ostream& ofs, const position_weight_matrix& M);

std::ostream& operator<<(std::ostream& ofs, const position_weight_matrix::match& M);


}


#endif

