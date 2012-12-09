// bbq - Copyright (C) 2005, Axel Mosig, University of Leipzig. See
// main.cpp for details.

#include <cstring>

#include <limits.h>
#include <float.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <ios>
#include <math.h>

#include <list>
#include <vector>
#include <string>
#include <algorithm>
#include <sequtil/triple.h>

#include "omit_pattern_matching.h"

#include <sequtil/position_weight_matrix.h>

using namespace bbq;

namespace fragrep {

	
/// A cluster triple is a triplet (p,e,j), where p is the position, e
/// the edit distance and j the fragment index of an occurence in some
/// T_i.
typedef std::triple< int,int,int > cluster_triple;
typedef std::set<cluster_triple> cluster_set;
typedef std::multimap<int,cluster_triple> cluster_map;
typedef std::vector<int> int_vector;
typedef std::vector<double> double_vector;
typedef std::vector<int_vector> int_int_vector;

#define max(x,y) ( ( (x) > (y) ) ? (x) : (y) )

typedef std::list< std::pair<int,float> > open_interval_list;


// ==============================================================	
bool operator<(const w_match_ref& x, const w_match_ref& y)
// ==============================================================	
{
	return (x->j<y->j || (x->j==y->j && x->pos<y->pos));
}
	
// ==============================================================	
result_item::result_item(int p, int j_, double pv, double ev, const std::string& seq)
// ==============================================================	
{
	pos=p; p_val=pv; E_val=ev;
	j=j_; match_seq=seq;
}

// ==============================================================	
match_pair::match_pair(int f, const std::string& seq)
// ==============================================================	
{
	first=f;
	second=false;
	match_seq=seq;
}

// ==============================================================	
weighted_match::weighted_match(int j_, int pos_, float weight_, const std::string& seq)
// ==============================================================	
{
	j=j_; pos=pos_; weight=weight_; score = 0.; match_seq = seq;
}
	
// ==============================================================	
weighted_match::weighted_match()
// ==============================================================	
{
	j=-1; pos=-1; weight=0; score = 0.;
}
	
// ==============================================================
omit_pattern_matching::omit_pattern_matching()
// ==============================================================
{

	init();

	use_pwms = false;
	
	reverted  = false;
	simple_output = false;
	
	k = 0;
	s.clear();
	n.clear();
	m.clear();
	l.clear();
	u.clear();
	T = "";
	T_name = "";

	INS_COST = 1.;
	DEL_COST = 1.;

	compute_index_table();
	compute_match_table();

}


// ==============================================================
void omit_pattern_matching::init(const std::string& T_, const std::string& T_name_, int k_, const string_vector& s_, int_vector l_, int_vector u_, int_vector m_)
// ==============================================================
{

	//init();
	
	reverted  = false;
	simple_output = false;
	k = k_;
	s = s_;
	T = T_;
	N = T.length();
	T_name = T_name_;
	n.resize(k);
	m = m_;	
	l = l_;
	u = u_;
	all_w_matches.clear();
	w_match_lists.resize(k);
	last_match_sequences.resize(k);

	int i;
	for (i=0; i<k; ++i)
	{
		n[i] = s[i].size();
		deletion_penalties[i] = .5;
		w_match_lists[i].clear();
		//w_match_lists[i].first = all_w_matches.end();
		//w_match_lists[i].second = all_w_matches.end();
	}

	match_reverse = false;
	reverted = false;
	prefix_length = 0;
	suffix_length = 0;

	compute_index_table();
	compute_match_table();
}


// ==============================================================
void omit_pattern_matching::init(int k_, const string_vector& s_, int_vector l_, int_vector u_, int_vector m_)
// ==============================================================
{
	//init();
	
	reverted  = false;
	simple_output = false;
	k = k_;
	s = s_;
	m = m_;	
	l = l_;
	u = u_;
	T = "";
	T_name = "";

	n.resize(k);
	deletion_penalties.resize(k);
	all_w_matches.clear();
	w_match_lists.resize(k);
	last_match_sequences.resize(k);
		
	int i;
	for (i=0; i<k; ++i)
	{
		n[i] = s[i].size();
		deletion_penalties[i] = .5;
		w_match_lists[i].clear();
		//w_match_lists[i].first = all_w_matches.end();
		//w_match_lists[i].second = all_w_matches.end();
	}
	
	match_reverse = false;
	reverted = false;
	prefix_length = 0;
	suffix_length = 0;

	compute_index_table();
	compute_match_table();
}

// ==============================================================
void omit_pattern_matching::init(int k_, const PWM_vector& PWMs_, int_vector l_, int_vector u_, float_vector t_)
// ==============================================================
{

	//init();
	use_pwms = true;
	
	reverted  = false;
	simple_output = false;
	k = k_;
	PWMs = PWMs_;
	all_w_matches.clear();

	deletion_penalties.resize(k);
	w_match_lists.resize(k);
	last_match_sequences.resize(k);
		
	int i;
	n.resize(k);
	for (i=0; i<k; ++i)
	{
		n[i] = PWMs[i].length;
		deletion_penalties[i] = .5;
		w_match_lists.clear();
		//w_match_lists[i].first = all_w_matches.end();
		//w_match_lists[i].second = all_w_matches.end();
	}
	
	t = t_;	
	l = l_;
	u = u_;
	T = "";
	s.clear();
	m.clear();
	T_name = "";

	match_reverse = false;
	reverted = false;
	prefix_length = 0;
	suffix_length = 0;

	compute_index_table();
	compute_match_table();
}

// ==============================================================
void omit_pattern_matching::init()
// ==============================================================
{
	compl_table.resize(256);

	int i;
	for(i=0; i<256; i++) compl_table[i]=(int)'*';
	
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
}

// ==============================================================
void omit_pattern_matching::set_output_region(int pref_len, int suff_len)
// ==============================================================
{
	prefix_length = pref_len;
	suffix_length = suff_len;
}

// ==============================================================
void omit_pattern_matching::set_reverse(bool rev)
// ==============================================================
{
	match_reverse = rev;
}

// ==============================================================
void omit_pattern_matching::set_unique(bool unq)
// ==============================================================
{
	report_unique_matches = unq;
}

// ==============================================================
	void omit_pattern_matching::set_delta(int delta_)
// ==============================================================
{
	delta = delta_;
}

// ==============================================================
void omit_pattern_matching::set_indel_costs(float ins, float del)
// ==============================================================
{
	INS_COST = ins;
	DEL_COST = del;
}

// ==============================================================
void omit_pattern_matching::set_deletions(int_vector d_)
// ==============================================================
{
	d = d_;
}

// ==============================================================
void omit_pattern_matching::set_verbose(bool verb)
// ==============================================================
{
	verbose = verb;
}
	
// ==============================================================
void omit_pattern_matching::set_simple_output(bool so)
// ==============================================================
{
	simple_output = so;
}
	
// ==============================================================
void omit_pattern_matching::set_scaffold_scoring(bool sc)
// ==============================================================
{
	scaffold_cutoff_scoring = sc;
}
	

// ==============================================================
	void omit_pattern_matching::set_genome(const std::string& T_, const std::string& T_name_)
// ==============================================================
{
	T = T_;
	T_name = T_name_;
	N = T.length();
}

// ==============================================================	
bool omit_pattern_matching::match_consensus_fragment(int j, int ii)
// ==============================================================
{
	int jj, m_cost, d_cost, i_cost;
	int m_j;
	std::string s_j;

	m_j = m[j];
	s_j = s[j];

	int n_j = n[j];
	int best_m_cost = 0;
	m_cost = 0;
	if (d[j]==0)
	{
		for (jj=0; (jj<n_j); jj++)
		{
			//std::cout<<"("<<s_j[jj]<<T[ii+jj]<<")";
			if (match_table[index_table[(int)s_j[jj]]][index_table[T[ii+jj]]]==0)
			//if (s_j[jj]!=T[ii+jj])
				m_cost++;
		}
		if (m_cost<=m_j) {
			//std::cout<<"*\n";
			done[j] = ii+1;
			return true;
		}
		//std::cout<<"\n";
	} else {
		// perform "smith waterman" style dynamic programming on the
		// partial order graph
		int d_j = d[j];
		int NN = n_j;
		
		// start w/ first column filled w/ all 0 for counting deletions
		int_int_vector column_1(d_j+1,int_vector(NN));
		int_int_vector column_2(d_j+1,int_vector(NN));
		for (int s=0;s<=d_j;s++)
		{
			for (int i=0;i<NN;i++)
			{
				column_1[s][i]=-1; column_2[s][i]=-1;
			}
		}
		//column_1[0][0]=0; column_2[0][0]=0;
		int_int_vector* prev_column = &column_1;
		int_int_vector* curr_column = &column_2;
		int_int_vector* tmp_column;
		
		for (jj=0; (jj<n_j); jj++)
		{
			int kk,ss;
			int min_kk=jj-d_j-1;
			if (min_kk<0)
				min_kk=0;
			min_kk=0;
			//int max_kk=jj;
			//if (max_kk>NN)
			int max_kk=NN;
			for (kk=min_kk; (kk<max_kk); kk++)
			{
				for (ss=0; (ss<=d_j); ss++)
				{
					if (kk>0)
					{
						m_cost=(*prev_column)[ss][kk-1];
					} else {
						if (jj>0 || ss>0)
							m_cost=-1;
						else
							m_cost=0;
					}
					// compute costs for (mis)match and deletion

					if (m_cost>=0 && match_table[index_table[(int)s_j[jj]]][index_table[T[ii+kk]]]==0)
						m_cost++;
					if (ss>0)
						d_cost=(*prev_column)[ss-1][kk];
					else
						d_cost=-1;
					if (m_cost<=d_cost && m_cost!=-1 || d_cost==-1) {
						(*curr_column)[ss][kk] = m_cost;
					} else {
						(*curr_column)[ss][kk] = d_cost;
					}
				}

			}
			//std::cout<<"\ndel column "<<jj<<": \n";
			/*for (int yy=0; yy<NN; yy++)
			{
				std::cout<<"(";
				for (int s=0; s<=d_j; s++) std::cout<<(*prev_column)[s][yy]<<" ";
				std::cout<<")";
			}
			std::cout<<"\n";
			std::cout.flush();*/
			// swap columns
			tmp_column = curr_column;
			curr_column = prev_column;
			prev_column = tmp_column;
		}
		/*for (int yy=0; yy<NN; yy++)
		{
			std::cout<<"(";
			for (int s=0; s<=d_j; s++) std::cout<<(*prev_column)[s][yy]<<" ";
			std::cout<<")";
		}
		std::cout<<"\n";*/
		int min_mm=INT_MAX;
		int min_kk = n_j-d_j-1;
		if (min_kk<0)
			min_kk=0;
		int max_kk = n_j;
		if (max_kk>NN)
			max_kk=NN;
		int min_d=0;
		for (int kk=min_kk; kk<max_kk; kk++)
		{
			int xxx=0;
			int this_d = n_j-kk-1;
			int this_mm;
			if (this_d>=0 && this_d<=d_j) {
				this_mm=(*prev_column)[this_d][kk];
			} else {
				this_mm = -1;
			}
			if (this_mm!=-1 && this_mm<min_mm) {
				min_mm=this_mm;
				min_d=this_d;
			}
			if (this_mm>=0 && this_mm<=m_j)
			{
				done[j]=ii+1;
				// create a dummy string to represent a match whose length
				// at least corresponds to the length of the actual match;
				// sooner or later this should be accompanied by tracing
				// back the aligned sequence throught the DP matrix.
				last_match_sequences[j]=std::string(n_j-this_d,'*');
				return true;
			}
		}
	}
	return false;
}


// ==============================================================	
bool omit_pattern_matching::match_fragment(int j, int ii)
// ==============================================================
{	
	if (!use_pwms)
	{
		return match_consensus_fragment(j,ii);
	} else {
		int n_j=n[j]+d[j];
		//int n_j=n[j];
		if (ii+n_j>N)
			n_j=N-ii;
		char T_ii[n_j+1];
		strncpy(T_ii,&(T[ii]),n_j);
		T_ii[n_j]='\0';
		//float_pair fp = PWMs[j].get_frac_score(T_ii);
		//if (fp.first>=t[j]) {
		if (PWMs[j].match_frac_score(T_ii,t[j])) {
			//last_match_sequences[j]=std::string(PWMs[j].length,'*');
			last_match_sequences[j]=PWMs[j].get_match_sequence();
			done[j] = ii+1;
			return true;
		}
	}
	return false;
}
	
// ==============================================================	
bool omit_pattern_matching::match_fragment(int j, int ii, weighted_match& w_match)
// ==============================================================
{	
	if (!use_pwms)
	{
		return match_consensus_fragment(j,ii);
		w_match = weighted_match(j,ii,1.,last_match_sequences[j]);
	} else {
		int n_j=n[j];
		char T_ii[n_j+1];
		strncpy(T_ii,&(T[ii]),n_j);
		T_ii[n_j]='\0';
		//float_pair fp = PWMs[j].get_frac_score(T_ii);
		//if (fp.first>=t[j]) {
		if (PWMs[j].match_frac_score(T_ii,t[j])) {
			//done[j] = ii+1;
			//last_match_sequences[j]=std::string(PWMs[j].length,'*');
			//w_match = weighted_match(j,ii,1.,std::string(PWMs[j].length,'*'));
			last_match_sequences[j] = PWMs[j].get_match_sequence();
			w_match = weighted_match(j,ii,1.,PWMs[j].get_match_sequence());
			return true;
		}
	}
	return false;
}
	

// ==============================================================
	int omit_pattern_matching::get_next_occurence(int j, int up_to)
// ==============================================================
{

	if (up_to<done[j])
		return -1;
	if (up_to>N-n[j])
		up_to = N-n[j];
	
	int ii = done[j];
	
	for (; ii<=up_to; ii++)
	{
		if (match_fragment(j,ii))
			return ii;
		if (verbose)
			if (ii%100000==0)
				std::cerr<<'*';
	}
	done[j] = up_to+1;
	return -1;
	
}

// =============================================================
int omit_pattern_matching::get_next_occurence(int j, int up_to, weighted_match& w_match)
// ==============================================================
{

	if (up_to<done[j])
		return -1;
	if (up_to>N-n[j])
		up_to = N-n[j];
	
	int ii = done[j];
	
	for (; ii<=up_to; ii++)
	{
		if (match_fragment(j,ii,w_match))
		{
			w_match.DP_pred = all_w_matches.end();
			w_match.DP_succ = all_w_matches.end();
			done[j]=ii+1;
			return ii;
		}
		if (verbose)
			if (ii%100000==0)
				std::cerr<<'*';
	}
	done[j] = up_to+1;
	return -1;
	
}

// =============================================================
int omit_pattern_matching::get_j_inf_occurence(int& j_inf, weighted_match& w_match)
// ==============================================================
{

	int_list::iterator j_it = j_inf_list.begin();
	int_list::iterator j_end = j_inf_list.end();
	int_vector up_to(k);
	int ii=j_inf_start;
	int j=0;
	for (; j<k; j++)
	{
		up_to[j] = N-n[j];
	}
	
	if (ii>N)
		return -1;
	
	for (; ii<=N; ii++)
	{
		j_it = j_inf_list.begin();
		j_end = j_inf_list.end();
		for (; j_it!=j_end; ++j_it)
		{
			j_inf = *j_it;
			if (ii<=up_to[j_inf])
			{
				if (match_fragment(j_inf,ii,w_match))
				{
					// store the match position as the starting point for
					// the next call of this method
					j_inf_start=ii+1;
					// return the position of the match
					w_match.DP_pred = all_w_matches.end();
					w_match.DP_succ = all_w_matches.end();
					return ii;
				}
			}
		}
		if (verbose)
			if (ii%100000==0)
				std::cerr<<'*';
	}
	return -1;
}

/*
// ==============================================================
int omit_pattern_matching::get_next_occurence(int_list j_list, int& curr_j_inf, int up_to)
// ==============================================================
{

	curr_j_inf=0;
	int_list::iterator j_it = j_list.begin();
	int_list::iterator j_end = j_list.end();
	if (j_it==j_end)
		return -1;
	
	int min_ii = *j_it;
	bool all_done = true;
	for (; j_it!=j_end; ++j_it)
	{
		int curr_j = *j_it;
		if (up_to>=done[curr_j])
			all_done=false;
		if (up_to>N-n[curr_j])
			up_to = N-n[curr_j];
		if (done[curr_j]<min_ii)
			min_ii=done[curr_j];
	}
	int ii = min_ii;
	if (all_done)
		return -1;
	
	for (; ii<=up_to; ii++)
	{
		if (verbose && ii%100000==0)
			std::cerr<<'*';
		j_it = j_list.begin();
		for (; j_it!=j_end; ++j_it)
		{
			if (ii>=done[*j_it])
			{
				//int curr_j=*j_it;
				//done[curr_j]=*j_it;
				if (match_fragment(*j_it,ii))
				{
					curr_j_inf = *j_it;
					return ii;
				}
			}
		}
	}

	j_it = j_list.begin();
	for (; j_it!=j_end; ++j_it)
		done[*j_it] = up_to+1;
	return -1;
	
	}*/


// ==============================================================
bool omit_pattern_matching::check_match_lists()
// ==============================================================
{

	int j;
	bool hope_left;

	// if there is no hope left that we can find any match for the
	// current occurence of s[j_inf], hope_left is set to false, which
	// gives a speedup because the for loop in this procedure can be
	// terminated early. Initially, we have not lost all hopes for
	// success, so hope_left is set to true :)
	hope_left = true;
	
	match_list_vector::iterator match_list_j = match_lists.begin();
	match_list::iterator m_it = match_lists[0].begin();
	match_list::iterator m_end = match_lists[0].end();
	while (m_it != m_end)
	{
		m_it->second = true;
		m_it++;
	}
	for (j=1; j<k && hope_left; j++)
	{
		match_list_j = match_lists.begin()+j;
		match_list_vector::iterator match_list_prev = match_lists.begin()+(j-1);

		m_it = match_list_j->begin();
		m_end = match_list_j->end();
		hope_left = false;
		while (m_it != m_end)
		{
			m_it->second = false;
			match_list::iterator p_it = match_list_prev->begin();
			match_list::iterator p_end = match_list_prev->end();
			while (p_it != p_end)
			{
				if (p_it->second)
				{
					int m_pos = m_it->first;
					int p_pos = p_it->first;
					int dist = m_pos-p_pos-n[j-1];
					if (dist>=l[j] && dist<=u[j])
					{
						m_it->second = true;
						hope_left = true;
					}	
				}
				p_it++;
			}
			m_it++;
		}
	}
	return hope_left;
	
}

// ==============================================================
bool omit_pattern_matching::check_omit_match_lists()
// ==============================================================
{
	int j;

	// if there is no hope left that we can find any match for the
	// current occurence of s[j_inf], lost_columns is set to false, which
	// gives a speedup because the for loop in this procedure can be
	// terminated early. 
	int lost_columns = -1;
	
	w_match_list_vector::iterator match_list_j = w_match_lists.begin();
	w_match_list_vector::iterator match_list_end = w_match_lists.end();
	for (; match_list_j!=match_list_end && lost_columns<=delta; ++match_list_j)
		//for (; match_list_j!=match_list_end; ++match_list_j)
	{
		w_match_list::iterator m_it = match_list_j->begin();
		w_match_list::iterator m_end = match_list_j->end();
		
		// perform dynamic programming
		int curr_j=0;
		float curr_deletion_penalty=0.;

		// flag for early termination
		bool pred_found=false;
		for (; m_it != m_end; ++m_it)
		{
			m_it->DP_pred = all_w_matches.end();
			//
			int curr_j = m_it->j;
			int min_del = curr_j;
			float best_score = m_it->score;
			int curr_pos = m_it->pos;
			w_match_ref best_pred = all_w_matches.end();
			
			w_match_ref_set::iterator p_it = m_it->pred.begin();
			w_match_ref_set::iterator p_end = m_it->pred.end();
			for (; p_it!=p_end; ++p_it)
			{
				weighted_match p_match = **p_it;
				int curr_del = curr_j - p_match.j - 1 + p_match.deletions;
				if ( curr_del < min_del)
				{
					best_pred = *p_it;
					best_score = m_it->score + p_match.score;
					min_del = curr_del;
					pred_found=true;
				}
			}
			
			m_it->deletions = min_del;
			m_it->score = best_score;
			m_it->DP_pred = best_pred;
			if (best_pred!=all_w_matches.end())
			{
				best_pred->DP_succ = m_it;
			}
		}
		if (!pred_found)
			lost_columns++;
					
	}
	if (lost_columns>delta)
		return false;
	// check whether dynamic programming was successful by browsing
	// through the deletion entries of the partial order members in
	// w_match_lists[k-1]
	match_list_j = w_match_lists.begin()+k-1-delta;
	match_list_end = w_match_lists.end();
	for (; match_list_j!=match_list_end; ++match_list_j)
	{
		w_match_list::reverse_iterator r_it = match_list_j->rbegin();
		w_match_list::reverse_iterator r_end = match_list_j->rend();
		for (; r_it != r_end; ++r_it)
		{
			if (r_it->deletions+(k-1-r_it->j) <= delta)
			{
				return true;
			}
			if (r_it->j < k-1-delta)
				break;
		}	
	}
	return false;
}

// ==============================================================
bool omit_pattern_matching::update_match_lists(int pos)
// ==============================================================
{
	// if there is no hope left that we can find any match for the
	// current occurence of s[j_inf], hope_left is set to false, which
	// gives a sppedup because the for loop in this procedure can be
	// terminated early. Initially, we have not lost all hopes for
	// success, so hope_left is set to true :)
	bool hope_left = true;
	int j;
	for (j=0; j<k && hope_left; j++)
	{
		// remove old matches from match_lists[j]
		int curr_l_bound_j = pos+l_bound[j];
		match_list_vector::iterator match_list_j = match_lists.begin()+j;
		if (j!=j_inf)
		{
			match_list::iterator m_it = match_list_j->begin();
			match_list::iterator m_end = match_list_j->end();
			bool abort = false;
			if (match_list_j->empty())
				abort = true; // GODFUCKING GCC
			if (abort==false)
			{
				while ((m_it!=m_end) && !abort)
				{
					if (m_it->first < curr_l_bound_j)
					{
						m_it++;
						match_list_j->pop_front();
					} else {
						abort = true;
					}
				}
			}
			while (m_it!=m_end)
			{
				m_it->second = false;
				m_it++;
			}
			
			// add new matches to match_lists[j]
			int curr_u_bound_j = pos+u_bound[j];
			int j_pos;
			while ( (j_pos=get_next_occurence(j,curr_u_bound_j)) >= 0 )
			{
				match_list_j->push_back(match_pair(j_pos,last_match_sequences[j]));
			}
		} else {
			match_list_j->clear();
			match_list_j->push_back(match_pair(pos,last_match_sequences[j_inf]));
		}
		hope_left = !match_list_j->empty();
	}
	return hope_left;
}

/// for debug purposes only
void print_w_matches(w_match_list& L, w_match_ref L_it, w_match_ref L_end)
{
	for (; L_it!=L_end; ++L_it)
	{
		std::cerr<<"("<<L_it->j<<","<<L_it->pos<<"):";
	}
	std::cerr<<"\n";
}

// ==============================================================
bool omit_pattern_matching::update_omit_match_lists(
		int pos, int curr_j_inf, const weighted_match& j_inf_match)
// ==============================================================
{

	int empty_count = 0;
	int j;
	int lower_bound = pos-max_width;
	int upper_bound = pos+max_width;

	for (j=0; j<k && empty_count<=delta; j++)
	{
		// remove old matches from w_match_lists[j] note that we only need
		// to look at the l_bounds to decide which elements we may kick
		// out
		w_match_list_vector::iterator match_list_j = w_match_lists.begin()+j;
		bool abort = false;
		w_match_list::iterator m_it = match_list_j->begin();
		w_match_list::iterator m_end = match_list_j->end();
		while ((m_it!=m_end) && !abort)
		{
			if (m_it->pos < lower_bound)
			{
				// ERASE ALL REFERENCES TO *m_it IN THE PARTIAL ORDER GRAPH
				remove_from_partial_order(m_it);
				m_it++;
			} else {
				abort = true;
			}
		}
		match_list_j->erase(match_list_j->begin(),m_it);

		// mark all remaining matches as unvisited by the dynamic
		// programming to follow
		while (m_it!=m_end)
		{
			m_it->visited = false;
			m_it++;
		}
		
		// add new matches to w_match_lists[j]
		int upper_bound = pos+max_width;
		int j_pos;
		weighted_match j_match;
		while ( (j_pos=get_next_occurence(j,upper_bound,j_match)) >= 0 )
		{
			w_match_ref new_match = match_list_j->insert(match_list_j->end(),j_match);
			add_to_partial_order(new_match);
		}
		//print_w_matches(all_w_matches, all_w_matches.begin(),all_w_matches.end());
		if (match_list_j->empty())
			empty_count++;
		//std::cerr<<"LIST_"<<j<<":\n";
		//print_w_matches(*match_list_j, match_list_j->begin(),match_list_j->end());
	}
	return (empty_count<=delta);
	
// 	int lower_bound = pos-max_width;
// 	int upper_bound = pos+max_width;
	
// 	// insert j_inf_match into w_match_lists[curr_j_inf]
// 	w_match_interval_vector::iterator j_inf_list = w_match_lists.begin()+curr_j_inf;
// 	all_w_matches.insert(j_inf_list->second,j_inf_match);
// 	// if j_inf_list was empty before, the list starts at the newly
// 	// inserted element
// 	if (j_inf_list->second==j_inf_list->first)
// 	{
// 		j_inf_list->first--;
// 	}
	
// 	// if there is no hope left that we can find any match for the
// 	// current occurence of s[curr_j_inf], hope_left is set to false, which
// 	// gives a speedup since the for-loop in this procedure can be
// 	// terminated early. Initially, we have not lost all hopes for
// 	// success, indicated by a zero empty_count.

// 	int empty_count = update_lower_match_lists(pos,curr_j_inf);
// 	empty_count += update_upper_match_lists(pos,curr_j_inf);
	
}

// ==============================================================
void omit_pattern_matching::add_to_partial_order(const w_match_list::iterator& match_ref)
// ==============================================================
{
	w_match_list::iterator w_it;
	w_match_list::iterator w_end;
	int j = match_ref->j;
	int pos = match_ref->pos;
	match_ref->deletions = j;
	match_ref->visited = false;

	int jj = j-delta-1;
	if (jj<0)
		jj=0;
	
	for (; jj<j; ++jj)
	{
		w_it  = w_match_lists[jj].begin();
		w_end = w_match_lists[jj].end();
		for (; w_it!=w_end; ++w_it)
		{
			int w_j = w_it->j;
			int w_pos = w_it->pos;
			int pos_diff = pos-w_pos-n[jj];
			if (pos_diff>=omit_l[w_j][j] && pos_diff<=omit_u[w_j][j])
			{
				w_it->succ.insert(match_ref);
				match_ref->pred.insert(w_it);
			}
		}
	}

	if (j>=k-1)
		//no successors in the partial order to be updated
		return;
	for (jj=j+1; jj<k && jj<=j+delta+1; ++jj)
	{
		w_it = w_match_lists[jj].begin();
		w_end = w_match_lists[jj].end();
		for (; w_it!=w_end; ++w_it)
		{
			int w_j = w_it->j;
			int w_pos = w_it->pos;
			int pos_diff = w_pos-pos-n[j];
			if (pos_diff>=omit_l[j][w_j] && pos_diff<=omit_u[j][w_j])
			{
				w_it->pred.insert(match_ref);
				match_ref->succ.insert(w_it);
			}
		}
	}
}
	
// ==============================================================
void omit_pattern_matching::remove_from_partial_order(const w_match_list::iterator& match_ref)
// ==============================================================
{
	w_match_list::iterator w_it;
	w_match_list::iterator w_end;
	int j = match_ref->j;
	int pos = match_ref->pos;

	int jj = j-delta-1;
	if (jj<0)
		jj=0;
	
	for (; jj<k && jj<=j+delta+1; ++jj)
	{
		w_it = w_match_lists[jj].begin();
		w_end = w_match_lists[jj].end();
		for (; w_it!=w_end; ++w_it)
		{
			w_it->succ.erase(match_ref);
			w_it->pred.erase(match_ref);
		}
	}	
}


// ==============================================================
void omit_pattern_matching::update_bounds(int pos)
// ==============================================================
{
	int j;
	int_vector::iterator l_bound_j = l_bound.begin();
	for (j=0; j<k; j++)
	{
		int new_b = pos + (*l_bound_j);
		if ( (new_b>=0) && (new_b>done[j]) )
			done[j] = new_b;
		l_bound_j++;
	}
}

// ==============================================================
void omit_pattern_matching::update_omit_bounds(int pos, int curr_j_inf)
// ==============================================================
{
	int j;
	for (j=0; j<k; j++)
	{
		int new_b = pos-max_width;
		if (new_b>done[j])
			done[j]=new_b;
	}
// 	int j;
// 	int_vector::iterator l_bound_j = omit_l_bound[curr_j_inf].begin();
// 	for (j=0; j<k; j++)
// 	{
// 		int new_b = pos + (*l_bound_j);
// 		if ( (new_b>=0) && (new_b>done[j]) )
// 			done[j] = new_b;
// 		l_bound_j++;
// 	}
}

// ==============================================================
void omit_pattern_matching::print_unique_results()
// ==============================================================
{
	typedef std::map<int,result_item> pos_result_map;
	pos_result_map r_map;
	complete_result_pos_map::iterator r_it = all_results.begin();
	complete_result_pos_map::iterator r_end = all_results.end();
	for (; r_it!=r_end ; ++r_it)
	{
		print_result(r_it->second.second);
	}

}
	
// ==============================================================
void omit_pattern_matching::add_unique_result(result_list& result)
// ==============================================================
{
	int jj = 0;
	char dir_str[6];
	double E_val=0., p_val=1.;
	result_list::iterator r_it = result.begin();
	result_list::iterator r_end = result.end();
	if (r_it==r_end)
		return;
	int pos = r_it->pos+1;
	for (jj=0 ; r_it!=r_end ; ++r_it, jj++ )
	{
		p_val *= ((double)(u[jj]-l[jj]+1)/(double)N) * r_it->p_val;			
	}
	E_val = -log(p_val);
	complete_result_pos_map::iterator cr_pos = all_results.find(pos);
	if (cr_pos!=all_results.end())
	{
		float ex_E_val = cr_pos->second.first;
		if (ex_E_val>=E_val)
			return;
		all_results.erase(cr_pos);
	} else {
		// indicate that a match was found
		if (verbose)
			std::cerr<<'+';
	}
	all_results[pos] = complete_result(E_val,result);
		
}
	
// ==============================================================
void omit_pattern_matching::print_result(result_list& result)
// ==============================================================
{
	
	number_of_matches++;	
	/// PRINT HEADER
	int jj = 0;
	char dir_str[6];
	double E_val=0., p_val=1.;
	result_list::iterator r_it = result.begin();
	result_list::iterator r_end = result.end();
	if (r_it==r_end)
		return;
	for (jj=0 ; r_it!=r_end ; ++r_it, jj++ )
	{
		p_val *= ((double)(u[jj]-l[jj]+1)/(double)N) * r_it->p_val;			
	}
	E_val = -log(p_val);
	
	if (!simple_output)
	{
		if (reverted)
			strcpy(dir_str,"-bwd");
		else
			strcpy(dir_str,"-fwd");
	} else {
		if (reverted)
			strcpy(dir_str,"-");
		else
			strcpy(dir_str,"+");
	}
	int dir_pos;
	if (!reverted)
		dir_pos = result.begin()->pos+1;
	else
		dir_pos = N-(result.rbegin()->pos+n[k-1])+1;
	if (!simple_output)
		std::cout<<T_name<<dir_str<<":pos"<<dir_pos;
	else
		std::cout<<T_name<<dir_str<<dir_pos;
	if (!simple_output)
		std::cout<<" weight="<<E_val<<" p-value="<<p_val;
	std::cout<<"\n";

	/// PRINT SEQUENCE
	int seq_beg = result.front().pos-prefix_length;
	if (seq_beg<0)
		seq_beg = 0;
	int seq_end = result.back().pos+n[k-1]+suffix_length;
	if (seq_end > N)
		seq_end = N;
	int i;
	std::cout<<" ";
	if (reverted && simple_output) {
		for (i=seq_end-1; i>=seq_beg; i--)
			std::cout<<compl_table[T[i]];
	} else {
		for (i=seq_beg; i<seq_end; i++)
			std::cout<<T[i];
	}
	std::cout<<"\n";

	/// PRINT FRAGMENTS BELOW SEQUENCE
	if (!simple_output)
	{
		std::cout<<">matchsequence\n";
		int suff_len = suffix_length;
		if (result.begin()->pos<suff_len)
			suff_len = result.begin()->pos;
		std::cout<<" ";
		for (i=0; i<suff_len; i++)
			std::cout<<"-";
		r_it = result.begin();
		std::cout<<r_it->match_seq;
		int prev_pos = r_it->pos;
		int prev_len = r_it->match_seq.size();
		r_it++;
		r_end = result.end();
		for (jj=1 ; r_it!=r_end ; ++r_it, jj++ )
		{
			if (!use_pwms)
			{
				int gap_len = r_it->pos-prev_pos-prev_len;
				prev_pos = r_it->pos;
				prev_len = r_it->match_seq.size();
				for (i=0; i<gap_len; i++)
					std::cout<<"-";
				std::cout<<r_it->match_seq;
			} else {
				int gap_len = r_it->pos-prev_pos-prev_len;
				prev_pos = r_it->pos;
				prev_len = r_it->match_seq.size();
				for (i=0; i<gap_len; i++)
					std::cout<<"-";
				std::cout<<r_it->match_seq;
			}
		}
		int suff_beg = result.back().pos + prev_len;
		int suff_end = result.back().pos + prev_len + suffix_length;
		if (suff_end >= N)
			suff_end = N;
		for (jj=suff_beg; jj<suff_end; jj++)
			std::cout<<"-";
		std::cout<<"\n";
	}
}

// ==============================================================
void omit_pattern_matching::output_matches(int j, result_list& result)
// ==============================================================
{
	if (j>=0) {
		match_list::iterator m_it = match_lists[j].begin();
		match_list::iterator m_end = match_lists[j].end();
		for (; m_it!=m_end; ++m_it)
		{
			if (m_it->second)
			{
				int pos = m_it->first;
				int diff, n_pos;
				bool valid = true;
				if (!result.empty())
				{
					n_pos = result.front().pos;
					diff = n_pos-pos-n[j];
					if (diff<l[j+1] || diff>u[j+1])
						valid = false;
				}
				if (valid)
				{
					float p_val = p_value(j,pos);
					float E_val = -FLT_MAX;
					if (p_val>0)
						E_val=-log(p_val);
						
					result.push_front(
						result_item(pos,j,p_val,E_val,m_it->match_seq)
											);
					output_matches(j-1,result);
					result.pop_front();
				}
			}
		}
	} else {
		if (!report_unique_matches) {
			print_result(result);
		} else {
			if (result.empty())
				return;
			if (result.front().pos<=curr_first_pos)
			{
				curr_first_pos = result.front().pos;
				curr_first_match = result;
			}
			if (result.back().pos>=curr_last_pos)
			{
				curr_last_pos = result.back().pos;
				curr_last_match = result;
			}
		}
	}		
}

// ==============================================================
void omit_pattern_matching::output_matches()
// ==============================================================
{
	result_list result;
	result.clear();
	curr_first_pos=N;
	curr_last_pos=-1;
	output_matches(k-1,result);
	if (report_unique_matches && curr_last_pos>-1 && curr_first_pos>last_reported_match_pos)
	{
		print_result(curr_first_match);
		//print_result(curr_last_match);
		last_reported_match_pos = curr_last_pos;
	}
}

// ==============================================================
bool omit_pattern_matching::trace_back(w_match_ref tail, result_list& result)
// ==============================================================
{
	result.clear();
	bool reported = false;
	while (tail!=all_w_matches.end())
	{
		reported==reported || tail->visited;
		tail->visited = true;
		//std::cerr<<"("<<tail->j<<";"<<tail->pos<<"):";
		int t_j = tail->j;
		float t_sc = tail->score;
		float t_w = tail->weight;
		std::string t_seq = tail->match_seq;
		result.push_front(result_item(tail->pos,tail->j,tail->score,tail->weight,tail->match_seq));
		tail = tail->DP_pred;
	}
	if (!report_unique_matches)
		return reported;
	if (result.front().pos<=curr_first_pos)
	{
		curr_first_pos = result.front().pos;
		curr_first_match = result;
	}
	if (result.back().pos>=curr_last_pos)
	{
		curr_last_pos = result.back().pos;
		curr_last_match = result;
	}
	return reported;
	//std::cerr<<"\n";
}

// ==============================================================
void omit_pattern_matching::output_omit_matches()
// ==============================================================
{
	// check whether dynamic programming was successful by browsing
	// through the deletion entries of the partial order members in
	// w_match_lists[k-1]
	//w_match_list_vector::iterator match_list_j = w_match_lists.begin()+k-1-delta;
	/*w_match_list_vector::iterator match_list_j = w_match_lists.begin();
	w_match_list_vector::iterator match_list_end = w_match_lists.end();
	for (; match_list_j!=match_list_end; ++match_list_j)
	{
		w_match_list::iterator r_it = match_list_j->begin();
		w_match_list::iterator r_end = match_list_j->end();
		for (; r_it != r_end; ++r_it)
		{
			if (r_it->DP_pred!=all_w_matches.end())
			{
				if (r_it->deletions+(k-1-r_it->j) <= delta)
				{
					result_list result;
					trace_back(r_it,result);
				}
			}
		}
		}*/


	
	curr_first_pos=N;
	curr_last_pos=-1;
	result_list result;
	result.clear();
	// check whether dynamic programming was successful by browsing
	// through the deletion entries of the partial order members in
	// the last delta many w_match_lists
	w_match_list_vector::iterator w_list_j = w_match_lists.begin()+k-1-delta;
	w_match_list_vector::iterator w_list_end = w_match_lists.end();
	int jj=0;
	bool reported = false;
	for (; w_list_j!=w_list_end; ++w_list_j, ++jj)
	{
		w_match_ref m_it = w_list_j->begin();
		w_match_ref m_end = w_list_j->end();
		float best_score = -FLT_MAX;
		w_match_ref best_path_end = all_w_matches.end();
		int leftmost_pos = INT_MAX;
		int rightmost_pos = 0;
		for (; m_it != m_end; ++m_it)
		{
			if (m_it->visited==true)
				continue;
			weighted_match w_m = *m_it;
			
			if (m_it->deletions + (k-1-m_it->j) <= delta)
			{
				if (m_it->pos+n[m_it->j]>rightmost_pos)
				{
					rightmost_pos = m_it->pos+n[m_it->j];
				}
				if (m_it->score > best_score)
				{
					best_score = m_it->score;
					best_path_end = m_it;
				}
				reported = reported || trace_back(m_it,result);
				if (result.begin()!=result.end())
				{
					if (result.begin()->pos<leftmost_pos)
						leftmost_pos = result.begin()->pos;
				}
			}
		}
	}
	if (report_unique_matches && curr_last_pos>-1 && curr_first_pos>last_reported_match_pos)
	{
		print_result(curr_first_match);
		//print_result(curr_last_match);
		last_reported_match_pos = curr_last_pos;
	}	else if (!reported) {
		add_unique_result(result);
	}
	
}

// ==============================================================
void omit_pattern_matching::match()
// ==============================================================
{

	int j_inf_pos;
	
	last_reported_match_pos = -1;
	
	while ( (j_inf_pos=get_next_occurence(j_inf,N)) >= 0 )
	{
		update_bounds(j_inf_pos);
		if (update_match_lists(j_inf_pos))
			if (check_match_lists())
				output_matches();
	}
	//if (report_unique_matches)
	//	print_unique_results();
	
}

// ==============================================================
void omit_pattern_matching::omit_match()
// ==============================================================
{

	int j_inf_pos;
	int curr_j_inf;
	weighted_match j_inf_match;

	// initialization of the start position required by
	// get_j_inf_occurence()
	j_inf_start = 0;

	last_reported_match_pos = -1;
	
	while ( (j_inf_pos=get_j_inf_occurence(curr_j_inf,j_inf_match)) >= 0 )
	{
		update_omit_bounds(j_inf_pos, curr_j_inf);
		if (update_omit_match_lists(j_inf_pos, curr_j_inf, j_inf_match))
			if (check_omit_match_lists())
				output_omit_matches();
	}
	print_unique_results();
	
}

// ==============================================================
void omit_pattern_matching::run()
// ==============================================================
{
	
	int i, j, ii, jj, eta;

	number_of_matches = 0;
	
	get_freq();
	//if (!use_pwms)
	j_inf = get_j_inf();
	//else
	//j_inf = 0;

	done.resize(k);
 	l_bound.resize(k);
	u_bound.resize(k);
	match_lists.resize(k);
	w_match_lists.resize(k);
	match_lists[0].clear();
	all_w_matches.clear();
	last_match_sequences.resize(k);
	
	l_bound[j_inf] = 0;
	u_bound[j_inf] = 0;
 	
	for (j=0; j<k; j++)
	{
		done[j] = 0;
		match_lists[j] = match_list();
		match_lists[j].clear();
		PWMs[j].set_indel_costs(INS_COST,DEL_COST,d[j]);
	}
	for (j=j_inf; j>0; j--)
	{
		l_bound[j-1] = l_bound[j]-u[j]-n[j-1];
		u_bound[j-1] = u_bound[j]-l[j]-n[j-1];
	}
	for (j=j_inf+1; j<k; j++)
	{
		l_bound[j] = l_bound[j-1]+l[j]+n[j-1];
		u_bound[j] = u_bound[j-1]+u[j]+n[j-1];
	}

	match();

	if (verbose)
		std::cerr<<number_of_matches<<" match(es) found.\n";
	
}

// ==============================================================
void omit_pattern_matching::compute_omit_bounds()
// ==============================================================
{
	// compute omit_l_bounds and omit_u_bounds;
	// these bounds are used in update_omit_match_lists()
	int curr_j_inf;
	int j;
	for (curr_j_inf=0; curr_j_inf<k; curr_j_inf++)
	{
		omit_l_bound[curr_j_inf][curr_j_inf] = 0;
		omit_u_bound[curr_j_inf][curr_j_inf] = 0;
		for (j=curr_j_inf; j>0; j--)
		{
			omit_l_bound[curr_j_inf][j-1] = omit_l_bound[curr_j_inf][j]-u[j]-n[j-1];
			omit_u_bound[curr_j_inf][j-1] = omit_u_bound[curr_j_inf][j]-l[j]-n[j-1];
		}
		for (j=curr_j_inf+1; j<k; j++)
		{
			omit_l_bound[curr_j_inf][j] = omit_l_bound[curr_j_inf][j-1]+l[j]+n[j-1];
			omit_u_bound[curr_j_inf][j] = omit_u_bound[curr_j_inf][j-1]+u[j]+n[j-1];
		}
	}
}

// ==============================================================
void omit_pattern_matching::compute_omit_ul()
// ==============================================================
{
	// compute omit_l and omit_u;
	// these bounds are used in the dynamic programming in check_omit_bounds();
	// note the slight difference to omit_l_bound and omit_u_bound...
	int b;
	int j;
	max_width=0;
	for (b=0; b<k; b++)
	{
		max_width += (u[b] + n[b]);
		for (j=0; j<b; j++)
		{
			omit_l[b][j] = 0;
			omit_u[b][j] = -1;
		}

		omit_l[b][b] = -n[b];
		omit_u[b][b] = -n[b];
		for (j=b+1; j<k; j++)
		{
			omit_l[b][j] = omit_l[b][j-1]+l[j]+n[j-1];
			omit_u[b][j] = omit_u[b][j-1]+u[j]+n[j-1];
		}
	}
}

// ==============================================================
void omit_pattern_matching::run(int delta_)
// ==============================================================
{
	
	int i, j, ii, jj, eta;

	if (delta_==0)
	{
		run();
		return;
	}
	if (delta>=k)
		delta=k-1;

	number_of_matches = 0;
	
	get_freq();
	//if (!use_pwms)
	get_j_inf_list();
	
	//else
	//j_inf = 0;

	done.resize(k);
 	omit_l_bound.resize(k);
 	omit_u_bound.resize(k);
	l_bound.resize(k);
	u_bound.resize(k);
	match_lists.resize(k);
	match_lists[0].clear();
	omit_u.resize(k);
	omit_l.resize(k);
	w_match_lists.resize(k);
	last_match_sequences.resize(k);
		
	all_w_matches.clear();
 	
	for (j=0; j<k; j++)
	{
		done[j] = 0;
		match_lists[j] = match_list();
		match_lists[j].clear();
		omit_l_bound[j].resize(k);
		omit_u_bound[j].resize(k);
		w_match_lists[j].clear();
		omit_l[j].resize(k);
		omit_u[j].resize(k);
		PWMs[j].set_indel_costs(INS_COST,DEL_COST,d[j]);
	}

	compute_omit_bounds();
	compute_omit_ul();

	omit_match();

	if (verbose)
		std::cerr<<number_of_matches<<" match(es) found.\n";
	
}


/// compute the reverse complement of sequence S
//==============================================================
	void omit_pattern_matching::reverse_compl(char* S, int len_S)
// ==============================================================
{
	int iota;
	int half = len_S/2;
	for (iota=0; iota<half; iota++)
	{
		char tmp = compl_table[S[iota]];
		S[iota] = compl_table[S[len_S-iota-1]];
		S[len_S-iota-1] = tmp;
	}
	if (len_S%2==1)
	{
		S[half] = compl_table[S[half]] ;
	}
}
	
/// compute the reverse complement of the patterns, so that in effect
/// the genome is searched in backward direction for L-occurences.
//==============================================================
void omit_pattern_matching::reverse_compl()
// ==============================================================
{

	int i;
	int j, iota;

	reverted = !reverted;

	// any previous result is invalidated
	all_results.clear();


	int half = N/2;
	for (iota=0; iota<half; iota++)
	{
		char tmp = compl_table[T[iota]];
		T[iota] = compl_table[T[N-iota-1]];
		T[N-iota-1] = tmp;
	}
	if (N%2==1)
	{
		T[half] = compl_table[T[half]] ;
	}

}

/// compute the reverse complement of the patterns, so that in effect
/// the genome is searched in backward direction for L-occurences.
//==============================================================
void omit_pattern_matching::reverse_compl_old()
// ==============================================================
{

	int i;
	int j, iota;

	reverted = !reverted;

	// any previous result is invalidated
	all_results.clear();

	if (!use_pwms)
		for (j=0; j<k; j++)
		{
			// turn s[j] into its reverse complement
			int n_j = n[j];
			int half = n_j/2;
			for (iota=0; iota<half; iota++)
			{
				char tmp = compl_table[s[j][iota]];
				s[j][iota] = compl_table[s[j][n_j-iota-1]];
				s[j][n_j-iota-1] = tmp;
			}
			if (n_j%2==1)
			{
				s[j][half] = compl_table[s[j][half]] ;
			}
		}
	int k_2 = k/2;
	for (j=0; j<k_2; j++)
	{
		if (!use_pwms)
		{
			std::string s_tmp = s[j];
			s[j] = s[k-j-1];
			s[k-j-1] = s_tmp;
			
			int m_tmp = m[j];
			m[j] = m[k-j-1];
			m[k-j-1] = m_tmp;		
		} else {
			PWM PWM_tmp = PWMs[j];
			PWMs[j] = PWMs[k-j-1];
			PWMs[k-j-1] = PWM_tmp;
		}
		int n_tmp = n[j];
		n[j] = n[k-j-1];
		n[k-j-1] = n_tmp;
		
	}

	k_2 = (k-1)/2+1;
	for (j=1; j<k_2; j++)
	{
		int u_tmp = u[j];
		u[j] = u[k-j];
		u[k-j] = u_tmp;
		
		int l_tmp = l[j];
		l[j] = l[k-j];
		l[k-j] = l_tmp;
	}

}


// ==============================================================
void omit_pattern_matching::compute_match_table()
// ==============================================================
{

	int i, j, k, flag_match;

	bool case_sense = false;
	unsigned char asymmetric = 1;

	match_table.resize(31);
	for(i=0; i<31; i++) 
		match_table[i].resize(31);

	unsigned char letters[][8] = {
		//	 A,C,G,T, a,c,g,t
		{1,0,0,0, 0,0,0,0}, // A
	   {0,1,0,0, 0,0,0,0}, // C
	   {0,0,1,0, 0,0,0,0}, // G
	   {0,0,0,1, 0,0,0,0}, // T
	   {1,0,1,0, 0,0,0,0}, // R
	   {0,1,0,1, 0,0,0,0}, // Y
	   {1,1,0,0, 0,0,0,0}, // M
	   {0,0,1,1, 0,0,0,0}, // K
	   {0,1,1,0, 0,0,0,0}, // S
	   {1,0,0,1, 0,0,0,0}, // W
	   {1,1,0,1, 0,0,0,0}, // H
	   {0,1,1,1, 0,0,0,0}, // B
	   {1,1,1,0, 0,0,0,0}, // V
	   {1,0,1,1, 0,0,0,0}, // D
	   {1,1,1,1, 0,0,0,0}, // N
		
		{0,0,0,0, 1,0,0,0}, // a
	   {0,0,0,0, 0,1,0,0}, // c
	   {0,0,0,0, 0,0,1,0}, // g
	   {0,0,0,0, 0,0,0,1}, // t
	   {0,0,0,0, 1,0,1,0}, // r
	   {0,0,0,0, 0,1,0,1}, // y
	   {0,0,0,0, 1,1,0,0}, // m
	   {0,0,0,0, 0,0,1,1}, // k
	   {0,0,0,0, 0,1,1,0}, // s
	   {0,0,0,0, 1,0,0,1}, // w 
	   {0,0,0,0, 1,1,0,1}, // h
	   {0,0,0,0, 0,1,1,1}, // b
	   {0,0,0,0, 1,1,1,0}, // v 
	   {0,0,0,0, 1,0,1,1}, // d
	   {0,0,0,0, 1,1,1,1}, // n
 
		{0,0,0,0, 0,0,0,0}  // all other chars;
	};

	for(i=0; i<31; i++)
	{
		
		for(j=0; j<31; j++)
		{
			// Default: no intersection present; 
			match_table[i][j] = 0;                       
			match_table[j][i] = 0;
			
			flag_match = 0;

			// watch for intersections; 
			for(k=0; k<8; k++) {                              
				
				if(case_sense) {
					// mind capital letters; 
					if(letters[i][k] + letters[j][k] == 2) {    
						flag_match = 1;
						break;
					}
				}
				
				else {
					// do not care about capitals; 
					if(k<4) {
						if(letters[i][k] + letters[j][k] + letters[i][k+4] + letters[j][k+4] >= 2) {
							flag_match = 1;
							break;
						}
					}
					else {
						if (letters[i][k] + letters[j][k] + letters[i][k-4] + letters[j][k-4] >= 2) {
							flag_match = 1;
							break;
						}
					}
				}
			}
			
			if(flag_match==1) {
				// there was a match -> non-empty intersection;    
				match_table[i][j] = 1;                    
				match_table[j][i] = 1;
			}
		}
	}

	if (asymmetric)
	{
		for(i=0; i<31; i++)
		{
			match_table[i][14] = 0;
			match_table[i][29] = 0;
		}
		match_table[14][14] = 1;
		match_table[14][29] = 1;
		match_table[29][14] = 1;
		match_table[29][29] = 1;
	}
}


// ==============================================================
void omit_pattern_matching::clear()
// ==============================================================
{
	int i,j,nu;
}


// ==============================================================
omit_pattern_matching::~omit_pattern_matching()
// ==============================================================
{
	clear();
}


// ==============================================================
int omit_pattern_matching::get_j_inf()
// ==============================================================
{
	int j;
	int my_j_inf = 0;
	double E_inf = E_value(0,-1);

	if (!use_pwms) {
		for (j=1; j<k; j++)
		{
			double E_val = E_value(j,-1);
			if (E_val<E_inf)
			{
				my_j_inf = j;
				E_inf = E_val;
			}
		}
	} else {
		float best_ic=0.;
		for (j=1; j<k; j++)
		{
			// sort the matrix columns ascending by their information content
			PWM::float_vector ic = PWMs[j].information_content;
			int len=ic.size()-d[j];
			sort(ic.begin(),ic.end());
			int kk;
			float ic_j=0;
			// sum up the (length-d[j]) least informative columns, in order
			// to obtain the 'worst case information content' regarding
			// deletion of the most informative columns.
			for (kk=0; kk<len; kk++)
				ic_j += ic[kk];
			if (ic_j>=best_ic)
				my_j_inf=j;
		}
	}

	return my_j_inf;
	
}

// ==============================================================
void omit_pattern_matching::get_j_inf_list()
// ==============================================================
{
	int j;
	
	float_int_vector frag_ranks;
	frag_ranks.resize(k);
	
	if (!use_pwms) {
		for (j=1; j<k; j++)
		{
			double E_val = E_value(j,-1);
			frag_ranks[j] = float_int(E_val,j);
		}
	} else {
		float best_ic=0.;
		for (j=1; j<k; j++)
		{
			// sort the matrix columns ascending by their information content
			PWM::float_vector ic = PWMs[j].information_content;
			int len=ic.size()-d[j];
			sort(ic.begin(),ic.end());
			int kk;
			float ic_j=0;
			// sum up the (length-d[j]) least informative columns, in order
			// to obtain the 'worst case information content' regarding
			// deletion of the most informative columns.
			for (kk=0; kk<len; kk++)
				ic_j += ic[kk];
			frag_ranks[j] = float_int(ic_j,j);
		}
	}

	j_inf_list.clear();
	sort(frag_ranks.begin(),frag_ranks.end());
	int i;
	
	// get the delta+1 highest ranking fragments into j_inf_list
	for (i=k-1; i>=k-delta-1 && i>0; i--)
		j_inf_list.push_back(frag_ranks[i].second);
		
	
}

// ==============================================================
double omit_pattern_matching::E_value(int j, int pos)
// ==============================================================
{
	double p_val = p_value(j,pos);

	/// non-positive p-values don't make sense -- thus we better rank those by an E-value of 0.
	if (p_val<=0.)
		return 0.;
	
	return -log(p_val);
}


// ==============================================================
double omit_pattern_matching::p_value(int j, int pos)
// ==============================================================
{

	if (j>=k)
		return 1.;

	if (use_pwms)
	{
		double_vector bg(4);
		bg[0]=.25; bg[1]=.25; bg[2]=.25; bg[3]=.25;
		float res;
		if (pos>=0)
		{
			int n_j=n[j]+d[j];
			//int n_j=n[j];
			int ii=pos;
			if (ii+n_j>N)
				n_j=N-ii;
			char T_ii[n_j+1];
			strncpy(T_ii,&(T[ii]),n_j);
			T_ii[n_j]='\0';
			//float_pair fp = PWMs[j].get_frac_score(T_ii);
			//if (fp.first>=t[j]) {
			res = PWMs[j].get_fwd_frac_score(T_ii);
		} else {
			res = t[j];
		}
		return PWMs[j].get_p_value(res,bg);
	}
	
	int jj, kk;
	double seq_N, L;
	seq_N = (double)N;	
	if (j==0)
		L = (double)N;
	else
		L = (double)(u[j]-l[j]);
	if (L<1.)
		L = 1.;

	double p = 1.;
	
	std::string txt;
	if (pos>=0)
		txt = T;
	else {
		txt = s[j];
		pos = 0;
	}
	//if (s[j][0]==txt[pos])
	if (match_table[index_table[(int)s[j][0]]][index_table[txt[pos]]]==1)	
		p *= (double)frequencies[index_table[s[j][0]]]/seq_N;
	int kk_end = n[j]-1;
	for (kk=0, jj=pos; kk<kk_end; jj++, kk++)
	{
		//if (s[j][kk]==txt[jj])
		if (match_table[index_table[(int)s[j][kk]]][index_table[txt[jj]]]==1)
			//if (s[j][kk+1]==txt[jj+1])
			if (match_table[index_table[(int)s[j][kk+1]]][index_table[txt[jj+1]]]==1)
			{
				double p_local =
					(double)di_frequencies[index_table[txt[jj]]][index_table[txt[jj+1]]]/seq_N;
				p *= p_local;
			} else {
				double p_local = (double)frequencies[index_table[txt[jj+1]]]/seq_N;					
				p *= p_local;
			}
	}
	if (p<1.)
		return ( 1 - exp( log(1-p)*L ) ); 
	else
		return 1.;
	
}

// ==============================================================
void omit_pattern_matching::get_freq()
// ==============================================================
{
	int a,i,j;
	di_frequencies.clear();
	frequencies.clear();
	di_frequencies.resize(32);
	frequencies.resize(32);

	for (i=0; i<32; i++)
	{
		// we initialize all absolute frequencies with 4 and all
		//	dinucl.-freq with 1, so that we do not encounter 0
		//	probabilities, which may yield unwanted 0 probabilities for
		//	p-values. This corresponds to pretending the presence of a
		//	CCGGTTAACAGCTGATC in the genome sequence) 
		frequencies[i] = 4;
		di_frequencies[i].resize(32);
		for (j=0; j<32; j++)
			di_frequencies[i][j] = 1;
	}
	
	int nn = N;
	for (i=0; i<nn-1; i++)
	{
		if ( (T[i]=='C' || T[i]=='G' || T[i]=='T' || T[i]=='A'
				  || T[i]=='c' || T[i]=='g' || T[i]=='t' || T[i]=='a') && 
			(T[i+1]=='C' || T[i+1]=='G' || T[i+1]=='T' || T[i+1]=='A'
				|| T[i+1]=='c' || T[i+1]=='g' || T[i+1]=='t' || T[i+1]=='a') )
		{
			di_frequencies[index_table[T[i]]][index_table[T[i+1]]]++;
			frequencies[index_table[T[i]]]++;
		}
	}
	
	sum_up_frq('M','A','C');
	sum_up_frq('R','A','G');
	sum_up_frq('W','A','T');
	sum_up_frq('S','C','G');
	sum_up_frq('Y','C','T');
	sum_up_frq('K','G','T');
	sum_up_frq('V','M','G');
	sum_up_frq('H','M','T');
	sum_up_frq('D','R','T');
	sum_up_frq('B','S','T');
	sum_up_frq('N','M','K');
}

// ==============================================================
void omit_pattern_matching::sum_up_frq(char O, char P, char Q)
// ==============================================================	
{
	const char* X = "ACGTMRWSYKVHDBN";
	int i;
	
	frequencies[index_table[O]] = frequencies[index_table[P]] + frequencies[index_table[Q]];

	for (i=0; i<15; i++)
	{
		di_frequencies[index_table[X[i]]][index_table[O]] =
			  di_frequencies[index_table[X[i]]][index_table[P]] 
			+ di_frequencies[index_table[X[i]]][index_table[Q]];
		
		di_frequencies[index_table[O]][index_table[X[i]]] =
			  di_frequencies[index_table[P]][index_table[X[i]]] 
			+ di_frequencies[index_table[Q]][index_table[X[i]]];
	}
}

// ==============================================================
void omit_pattern_matching::compute_index_table()
// ==============================================================	
{

	int i;
	
	/// by default, we want to use asymmetric match tables, i.e.,
	/// wildcard Ns are only wildcards if they occur in a pattern
	/// fragment, but not in the text.
	bool asymmetric = true;

	index_table.resize(256);

	for(i=0; i<256; i++) index_table[i] = 30;
	
	index_table[65] = 0; // A 
	index_table[67] = 1; // C 
	index_table[71] = 2; // G 
	index_table[84] = 3; // T 
	index_table[85] = 3; // U == T 
	index_table[82] = 4; // R 
	index_table[89] = 5; // Y 
	index_table[77] = 6; // M 
	index_table[75] = 7; // K 
	index_table[83] = 8; // S 
	index_table[87] = 9; // W 
	index_table[72] = 10;// H 
	index_table[66] = 11;// B 
	index_table[86] = 12;// V 
	index_table[68] = 13;// D 
	if (asymmetric)
		index_table[78] = 14;  // N 
	
	index_table[97] = 0;  // a 
	index_table[99] = 1;  // c 
	index_table[103] = 2; // g  
	index_table[116] = 3; // t 
	index_table[117] = 3; // u == t 
	index_table[114] = 4; // r 
	index_table[121] = 5; // y 
	index_table[109] = 6; // m 
	index_table[107] = 7; // k 
	index_table[115] = 8; // s 
	index_table[119] = 9; // w 
	index_table[104] = 10; // h 
	index_table[98] = 11;  // b 
	index_table[118] = 12; // v 
	index_table[100] = 13; // d 
	if (asymmetric)
		index_table[110] = 14; // n 
	//index_table[97] = 15;  // a 
	//index_table[99] = 16;  // c 
	//index_table[103] = 17; // g  
	//index_table[116] = 18; // t 
	//index_table[117] = 18; // u == t 
	//index_table[114] = 19; // r 
	//index_table[121] = 20; // y 
	//index_table[109] = 21; // m 
	//index_table[107] = 22; // k 
	//index_table[115] = 23; // s 
	//index_table[119] = 24; // w 
	//index_table[104] = 25; // h 
	//index_table[98] = 26;  // b 
	//index_table[118] = 27; // v 
	//index_table[100] = 28; // d 
	//if (asymmetric)
		//index_table[110] = 29; // n 

}

} // end namespace fragrep
