// bbq - Copyright (C) 2005, Axel Mosig, University of Leipzig. See
// main.cpp for details.

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

#include "pattern_matching.h"

#include "sequtil/position_weight_matrix.h"

using namespace bbq;

namespace fragrep {

	
/// A cluster triple is a triplet (p,e,j), where p is the position, e
/// the edit distance and j the fragment index of an occurence in some
/// T_i.
typedef std::triple< int,int,int > cluster_triple;
typedef std::set<cluster_triple> cluster_set;
typedef std::multimap<int,cluster_triple> cluster_map;
typedef std::vector<int> int_vector;
typedef std::vector<int_vector> int_int_vector;

#define max(x,y) ( ( (x) > (y) ) ? (x) : (y) )

typedef std::list< std::pair<int,float> > open_interval_list;

// ==============================================================	
result_item::result_item(int p, double pv, double ev)
// ==============================================================	
{
	pos=p; p_val=pv; E_val=ev;
}

// ==============================================================
pattern_matching::pattern_matching()
// ==============================================================
{

	init();

	use_pwms = false;
	
	reverted  = false;
	simple_output = false;
	
	k = 0;
	s = NULL;
	n = NULL;
	m = NULL;	
	l = NULL;
	u = NULL;
	T = NULL;
	T_name = "";

	compute_index_table();
	compute_match_table();

}


// ==============================================================
void pattern_matching::init(char* T_, const std::string& T_name_, int N_, int k_, char** s_, int* n_, int* l_, int* u_, int* m_)
// ==============================================================
{

	//init();
	
	reverted  = false;
	simple_output = false;
	k = k_;
	N = N_;
	s = s_;
	T = new char[N_+1000];
	strcpy(T,T_);
	T_name = T_name_;
	n = n_;
	m = m_;	
	l = l_;
	u = u_;

	match_reverse = false;
	reverted = false;
	prefix_length = 0;
	suffix_length = 0;

	compute_index_table();
	compute_match_table();
}


// ==============================================================
void pattern_matching::init(int k_, char** s_, int* n_, int* l_, int* u_, int* m_)
// ==============================================================
{

	//init();
	
	reverted  = false;
	simple_output = false;
	k = k_;
	s = s_;
	n = n_;
	m = m_;	
	l = l_;
	u = u_;
	T = NULL;
	T_name = "";

	match_reverse = false;
	reverted = false;
	prefix_length = 0;
	suffix_length = 0;

	compute_index_table();
	compute_match_table();
}

// ==============================================================
void pattern_matching::init(int k_, const PWM_vector& PWMs_, int* n_, int* l_, int* u_, float* t_)
// ==============================================================
{

	//init();
	use_pwms = true;
	
	reverted  = false;
	simple_output = false;
	k = k_;
	PWMs = PWMs_;
	n = n_;
	t = t_;	
	l = l_;
	u = u_;
	T = NULL;
	s = NULL;
	m = NULL;
	T_name = "";

	match_reverse = false;
	reverted = false;
	prefix_length = 0;
	suffix_length = 0;

	compute_index_table();
	compute_match_table();
}

// ==============================================================
void pattern_matching::init()
// ==============================================================
{
	compl_table = new char[256];

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
void pattern_matching::set_output_region(int pref_len, int suff_len)
// ==============================================================
{
	prefix_length = pref_len;
	suffix_length = suff_len;
}

// ==============================================================
void pattern_matching::set_reverse(bool rev)
// ==============================================================
{
	match_reverse = rev;
}

// ==============================================================
void pattern_matching::set_unique(bool unq)
// ==============================================================
{
	report_unique_matches = unq;
}

// ==============================================================
void pattern_matching::set_deletions(int* d_)
// ==============================================================
{
	d = d_;
}

// ==============================================================
void pattern_matching::set_verbose(bool verb)
// ==============================================================
{
	verbose = verb;
}
	
// ==============================================================
void pattern_matching::set_simple_output(bool so)
// ==============================================================
{
	simple_output = so;
}
	

// ==============================================================
void pattern_matching::set_genome(const char* T_, const std::string& T_name_, int N_)
// ==============================================================
{
	if (T!=NULL)
	{
		delete[] T;
		T = NULL;
	}
	T = new char[N_+100];
	strcpy(T,T_);

	T_name = T_name_;
	N = N_;
}

// ==============================================================	
bool pattern_matching::match_PWM_fragment(int j, int ii)
// ==============================================================
{
}

// ==============================================================	
bool pattern_matching::match_consensus_fragment(int j, int ii)
// ==============================================================
{
	int jj, m_cost, d_cost, i_cost;
	int m_j;
	const char* s_j;

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
		// perform "smith waterman" style dynamic programming
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
					if (m_cost<=d_cost && m_cost!=-1 || d_cost==-1)
						(*curr_column)[ss][kk] = m_cost;
					else 
						(*curr_column)[ss][kk] = d_cost;
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
		for (int kk=min_kk; kk<max_kk; kk++)
		{
			int xxx=0;
			//for (int s=0; s<=d_j; s++)
			//{
				int this_d = n_j-kk-1;
				int this_mm;
				if (this_d>=0 && this_d<=d_j)
					this_mm=(*prev_column)[this_d][kk];
				else
					this_mm = -1;
				if (this_mm!=-1 && this_mm<min_mm)
					min_mm=this_mm;
				//if (this_mm<=m_j && this_mm>=0)
				//	std::cout<<"match w/ len "<<kk+1<<" and "<<this_d<<" deletion(s).\n";
				//}
				if (this_mm>=0 && this_mm<=m_j)
				{
					done[j]=ii+1;
					return true;
				}
		}
	}
	return false;
}


// ==============================================================	
bool pattern_matching::match_fragment(int j, int ii)
// ==============================================================
{	
	if (!use_pwms)
	{
		return match_consensus_fragment(j,ii);
	} else {
		int n_j=n[j];
		char T_ii[n_j+1];
		strncpy(T_ii,&(T[ii]),n_j);
		T_ii[n_j]='\0';
		float_pair fp = PWMs[j].get_frac_score(T_ii);
		//if ((!reverted && fp.first>=t[j]) || (reverted && fp.second>=t[j])) {
		if (fp.first>=t[j]) {
			//std::cerr<<"pwm "<<j<<" match at pos "<<ii<<"\n";
			done[j] = ii+1;
			return true;
		}
	}
	return false;
}
	

// ==============================================================
int pattern_matching::get_next_occurence(int j, int up_to)
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


// ==============================================================
bool pattern_matching::check_match_lists()
// ==============================================================
{

	int j;
	bool hope_left;

	// if there is no hope left that we can find any match for the
	// current occurence of s[j_inf], hope_left is set to false, which
	// gives a sppedup because the for loop in this procedure can be
	// terminated early. Initially, we have not lost all hopes for
	// success, so hope_left is set to true :)
	hope_left = true;
	
	match_list* match_list_j = &(match_lists[0]);
	match_list::iterator m_it = match_lists[0].begin();
	match_list::iterator m_end = match_lists[0].end();
	while (m_it != m_end)
	{
		m_it->second = true;
		m_it++;
	}
	for (j=1; j<k && hope_left; j++)
	{
		match_list_j = &(match_lists[j]);
		match_list* match_list_prev = &(match_lists[j-1]);

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
bool pattern_matching::update_match_lists(int pos)
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
		match_list* match_list_j = &(match_lists[j]);
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
				match_list_j->push_back(match_pair(j_pos,false));
			}
		} else {
			match_list_j->clear();
			match_list_j->push_back(match_pair(pos,false));
		}
		hope_left = !match_list_j->empty();
	}
	return hope_left;
}

// ==============================================================
void pattern_matching::update_bounds(int pos)
// ==============================================================
{
	int j;
	int* l_bound_j = &(l_bound[0]);
	for (j=0; j<k; j++)
	{
		int new_b = pos + (*l_bound_j);
		if ( (new_b>=0) && (new_b>done[j]) )
			done[j] = new_b;
		l_bound_j++;
	}
}

// ==============================================================
void pattern_matching::print_unique_results()
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
void pattern_matching::add_unique_result(result_list& result)
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
void pattern_matching::print_result(result_list& result)
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
		if (!use_pwms)
			std::cout<<s[0];
		else
			std::cout<<std::string(n[0],'*');
		r_it = result.begin();
		int prev_pos = r_it->pos;
		r_it++;
		r_end = result.end();
		for (jj=1 ; r_it!=r_end ; ++r_it, jj++ )
		{
			int gap_len = r_it->pos-prev_pos-n[jj-1];
			prev_pos = r_it->pos;
			for (i=0; i<gap_len; i++)
				std::cout<<"-";
			if (!use_pwms)
				std::cout<<s[jj];
			else
				std::cout<<std::string(n[jj],'*');
		}
		int suff_beg = result.back().pos + n[k-1];
		int suff_end = result.back().pos + n[k-1] + suffix_length;
		if (suff_end >= N)
			suff_end = N;
		for (jj=suff_beg; jj<suff_end; jj++)
			std::cout<<"-";
		std::cout<<"\n";
	}
}

// ==============================================================
void pattern_matching::output_matches(int j, result_list& result)
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
					result.push_front(result_item(pos,p_value(j,pos),E_value(j,pos)));
					output_matches(j-1,result);
					result.pop_front();
				}
			}
		}
	} else {
		if (!report_unique_matches)
			print_result(result);
		else
			add_unique_result(result);

	}		
}

// ==============================================================
void pattern_matching::output_matches()
// ==============================================================
{
	result_list result;
	result.clear();
	output_matches(k-1,result);
}

// ==============================================================
void pattern_matching::match()
// ==============================================================
{

	int j_inf_pos;
	
	while ( (j_inf_pos=get_next_occurence(j_inf,N)) >= 0 )
	{
		update_bounds(j_inf_pos);
		if (update_match_lists(j_inf_pos))
			if (check_match_lists())
				output_matches();
	}
	if (report_unique_matches)
		print_unique_results();
	
}

// ==============================================================
void pattern_matching::run()
// ==============================================================
{
	
	int i, j, ii, jj, eta;

	number_of_matches = 0;
	
	get_freq();
	//if (!use_pwms)
	j_inf = get_j_inf();
	//else
	//j_inf = 0;

	done = new int[k];
 	l_bound = new int[k];
	u_bound = new int[k];
	match_lists = new match_list[k];
	match_lists[0].clear();
	
	l_bound[j_inf] = 0;
	u_bound[j_inf] = 0;
 	
	for (j=0; j<k; j++)
	{
		done[j] = 0;
		match_lists[j] = match_list();
		match_lists[j].clear();
		PWMs[j].set_indel_costs(1.,0.,d[j]);
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


/// compute the reverse complement of sequence S
//==============================================================
	void pattern_matching::reverse_compl(char* S, int len_S)
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
void pattern_matching::reverse_compl()
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
void pattern_matching::reverse_compl_old()
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
			char* s_tmp = s[j];
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
void pattern_matching::compute_match_table()
// ==============================================================
{

	int i, j, k, flag_match;

	bool case_sense = false;
	unsigned char asymmetric = 1;

	match_table = new unsigned char*[31];
	for(i=0; i<31; i++) 
		match_table[i] = new unsigned char[31];

	unsigned char letters[][8] = {
		/*	 A,C,G,T, a,c,g,t */
		
		{1,0,0,0, 0,0,0,0},    /* A */
	   {0,1,0,0, 0,0,0,0},    /* C */
	   {0,0,1,0, 0,0,0,0},    /* G */
	   {0,0,0,1, 0,0,0,0},    /* T */
	   {1,0,1,0, 0,0,0,0},	  /* R */
	   {0,1,0,1, 0,0,0,0},    /* Y */
	   {1,1,0,0, 0,0,0,0},    /* M */
	   {0,0,1,1, 0,0,0,0},    /* K */
	   {0,1,1,0, 0,0,0,0},    /* S */
	   {1,0,0,1, 0,0,0,0},    /* W */
	   {1,1,0,1, 0,0,0,0},    /* H */
	   {0,1,1,1, 0,0,0,0},    /* B */
	   {1,1,1,0, 0,0,0,0},    /* V */
	   {1,0,1,1, 0,0,0,0},    /* D */
	   {1,1,1,1, 0,0,0,0},    /* N */
		
		{0,0,0,0, 1,0,0,0},    /* a */
	   {0,0,0,0, 0,1,0,0},    /* c */
	   {0,0,0,0, 0,0,1,0},    /* g */
	   {0,0,0,0, 0,0,0,1},    /* t */
	   {0,0,0,0, 1,0,1,0},	  /* r */
	   {0,0,0,0, 0,1,0,1},    /* y */
	   {0,0,0,0, 1,1,0,0},    /* m */
	   {0,0,0,0, 0,0,1,1},    /* k */
	   {0,0,0,0, 0,1,1,0},    /* s */
	   {0,0,0,0, 1,0,0,1},    /* w */ 
	   {0,0,0,0, 1,1,0,1},    /* h */
	   {0,0,0,0, 0,1,1,1},    /* b */
	   {0,0,0,0, 1,1,1,0},    /* v */ 
	   {0,0,0,0, 1,0,1,1},    /* d */
	   {0,0,0,0, 1,1,1,1},    /* n */
 
		{0,0,0,0, 0,0,0,0}     /* all other chars; */
	};

	for(i=0; i<31; i++)
	{
		
		for(j=0; j<31; j++)
		{
			
			match_table[i][j] = 0;                       /* Default: no intersection present; */
			match_table[j][i] = 0;
			
			flag_match = 0;
			
			for(k=0; k<8; k++) {                              /* watch for intersections; */
				
				if(case_sense) {
					if(letters[i][k] + letters[j][k] == 2) {    /* mind capital letters; */
						flag_match = 1;
						break;
					}
				}
				
				else {                                         /* do not casre about capitals; */
					if(k<4) {
						if(letters[i][k] + letters[j][k] + letters[i][k+4] + letters[j][k+4] >= 2) {
							flag_match = 1;
							break;
						}
					}
					else {
						if(letters[i][k] + letters[j][k] + letters[i][k-4] + letters[j][k-4] >= 2) {
							flag_match = 1;
							break;
						}
					}
				}
			}
			
			if(flag_match==1) {
				match_table[i][j] = 1;                    /* there was a match -> non-empty intersection;    */
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
void pattern_matching::clear()
// ==============================================================
{
	int i,j,nu;
}


// ==============================================================
pattern_matching::~pattern_matching()
// ==============================================================
{
	clear();
}


// ==============================================================
int pattern_matching::get_j_inf()
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
double pattern_matching::E_value(int j, int pos)
// ==============================================================
{
	double p_val = p_value(j,pos);

	/// non-positive p-values don't make sense -- thus we better rank those by an E-value of 0.
	if (p_val<=0.)
		return 0.;
	
	return -log(p_val);
}


// ==============================================================
double pattern_matching::p_value(int j, int pos)
// ==============================================================
{

	if (j>=k)
		return 1.;

	if (use_pwms)
	{
		float res;
		if (pos>=0)
		{
			float_pair fp = PWMs[j].get_match_score(&(T[pos]));
			res = fp.first;
		} else {
			res = PWMs[j].get_information_content();
		}
			
		return res;
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
	
	const char* txt;
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
				double p_local = (double)di_frequencies[index_table[txt[jj]]][index_table[txt[jj+1]]]/seq_N;
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
void pattern_matching::get_freq()
// ==============================================================
{
	int a,i,j;
	di_frequencies = new int*[32];
	frequencies = new int[32];

	int* frq = frequencies;
	int** di_frq = di_frequencies;
	const char* txt = T;
		
	for (i=0; i<32; i++)
	{
		/* we initialize all absolute frequencies with 4 and all
			dinucl.-freq with 1, so that we do not encounter 0
			probabilities, which may yield unwanted 0 probabilities for
			p-values. This corresponds to pretending the presence of a
			CCGGTTAACAGCTGATC in the genome sequence) */
		frq[i] = 4;
		di_frq[i] = new int[32];
		for (j=0; j<32; j++)
			di_frq[i][j] = 1;
	}
	
	int nn = N;
	for (i=0; i<nn-1; i++)
	{
		if ( (txt[i]=='C' || txt[i]=='G' || txt[i]=='T' || txt[i]=='A'
				  || txt[i]=='c' || txt[i]=='g' || txt[i]=='t' || txt[i]=='a') && 
			(txt[i+1]=='C' || txt[i+1]=='G' || txt[i+1]=='T' || txt[i+1]=='A'
				|| txt[i+1]=='c' || txt[i+1]=='g' || txt[i+1]=='t' || txt[i+1]=='a') )
		{
			di_frq[index_table[txt[i]]][index_table[txt[i+1]]]++;
			frq[index_table[txt[i]]]++;
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
void pattern_matching::sum_up_frq(char O, char P, char Q)
// ==============================================================	
{
	const char* X = "ACGTMRWSYKVHDBN";
	int* frq = frequencies;
	int** di_frq = di_frequencies;
	int i;
	
	frq[index_table[O]] = frq[index_table[P]] + frq[index_table[Q]];

	for (i=0; i<15; i++)
	{
		di_frq[index_table[X[i]]][index_table[O]] =
			  di_frq[index_table[X[i]]][index_table[P]] 
			+ di_frq[index_table[X[i]]][index_table[Q]];
		
		di_frq[index_table[O]][index_table[X[i]]] =
			  di_frq[index_table[P]][index_table[X[i]]] 
			+ di_frq[index_table[Q]][index_table[X[i]]];
	}
}

// ==============================================================
void pattern_matching::compute_index_table()
// ==============================================================	
{

	int i;

	
	/// by default, we want to use asymmetric match tables, i.e.,
	/// wildcard Ns are only wildcards if they occur in a pattern
	/// fragment, but not in the text.
	bool asymmetric = true;

	index_table = new int[256];

	for(i=0; i<256; i++) index_table[i] = 30;
	
	index_table[65] = 0;   /* A */
	index_table[67] = 1;   /* C */
	index_table[71] = 2;   /* G */
	index_table[84] = 3;   /* T */
	index_table[85] = 3;   /* U == T */
	index_table[82] = 4;   /* R */
	index_table[89] = 5;   /* Y */
	index_table[77] = 6;   /* M */
	index_table[75] = 7;   /* K */
	index_table[83] = 8;   /* S */
	index_table[87] = 9;   /* W */
	index_table[72] = 10;  /* H */
	index_table[66] = 11;  /* B */
	index_table[86] = 12;  /* V */
	index_table[68] = 13;  /* D */
	if (asymmetric)
		index_table[78] = 14;  /* N */
	
	index_table[97] = 0;  /* a */
	index_table[99] = 1;  /* c */
	index_table[103] = 2; /* g  */
	index_table[116] = 3; /* t */
	index_table[117] = 3; /* u == t */
	index_table[114] = 4; /* r */
	index_table[121] = 5; /* y */
	index_table[109] = 6; /* m */
	index_table[107] = 7; /* k */
	index_table[115] = 8; /* s */
	index_table[119] = 9; /* w */
	index_table[104] = 10; /* h */
	index_table[98] = 11;  /* b */
	index_table[118] = 12; /* v */
	index_table[100] = 13; /* d */
	if (asymmetric)
		index_table[110] = 14; /* n */
	//index_table[97] = 15;  /* a */
	//index_table[99] = 16;  /* c */
	//index_table[103] = 17; /* g  */
	//index_table[116] = 18; /* t */
	//index_table[117] = 18; /* u == t */
	//index_table[114] = 19; /* r */
	//index_table[121] = 20; /* y */
	//index_table[109] = 21; /* m */
	//index_table[107] = 22; /* k */
	//index_table[115] = 23; /* s */
	//index_table[119] = 24; /* w */
	//index_table[104] = 25; /* h */
	//index_table[98] = 26;  /* b */
	//index_table[118] = 27; /* v */
	//index_table[100] = 28; /* d */
	//if (asymmetric)
		//index_table[110] = 29; /* n */

}

} // end namespace fragrep
