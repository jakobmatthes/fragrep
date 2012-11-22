#include "position_weight_matrix.h"
#include "bbq_tokenizer.h"
#include "bbq_util.h"
#include <list>
#include <math.h>
#include <limits>
#include <float.h>
#include <iomanip>

#ifdef USE_GSL
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#endif

typedef std::list<float> float_list;
typedef std::vector<float> float_vector;
typedef std::vector<float_vector> float_float_vector;

typedef std::vector<bool> bool_vector;
typedef std::vector<bool_vector> bool_bool_vector;


// ==============================================================
float _min(float a, float b, float c, float d)
// ==============================================================
{
	float min = a;
	if (b<=min)
		min = b;
	if (c<=min)
		min = c;
	if (d<=min)
		min = d;
	return min;
}

// ==============================================================
float _max(float a, float b, float c, float d)
// ==============================================================
{
	float max = a;
	if (b>=max)
		max = b;
	if (c>=max)
		max = c;
	if (d>=max)
		max = d;
	return max;
}

// ==============================================================
template <class number>
number max(number x, number y, number z)
// ==============================================================
{
	if (x>y)
		if (x>z)
			return x;
		else // z>=x
			if (z>y)
				return z;
			else // y>=z>=y
				return y;
	else // y>=x
		if (y>z)
			return y;
		else // z>=y>=x
			return z;
}


namespace bbq {

int position_weight_matrix::match_mode = 2;
	
// ==============================================================
position_weight_matrix::position_weight_matrix()
// ==============================================================
{
	std::cerr<<'P';std::cerr.flush();
	match_reverse_complement = false;
	length = 0;
	DNA = true;

	INS_COST = 1.;
	DEL_COST = 0.;
	max_dels = 0;
	iterations = 10;
	
	A.clear();
	C.clear();
	G.clear();
	T.clear();
	information_content.clear();
	freq.clear();

	curr_aln_seq.clear();
	curr_match_seq.clear();
	mis_sequence.clear();
	rev_mis_sequence.clear();
	curr_aln_path.clear();
	std::cerr<<':';std::cerr.flush();	
}


// ==============================================================
position_weight_matrix::position_weight_matrix(const float_vector& A_, const float_vector& C_, const float_vector& G_, const float_vector& T_)
// ==============================================================
{
	A = A_;
	C = C_;
	G = G_;
	T = T_;

	INS_COST = 1.;
	DEL_COST = 0.;
	max_dels = 0;
	DNA = true;
	iterations = 10;
	
	length = A.size();
	if (C.size()!=length || G.size()!=length || T.size()!=length)
		throw bbq_exception("incompatible vector length for PWM initalization.\n");
	init();
}

	/*position_weight_matrix::position_weight_matrix(const position_weight_matrix& M)
{
	//clear();
	match_reverse_complement = M.match_reverse_complement;
	length = M.length;
	A = M.A;
	C.reserve(length);
	G.reserve(length);
	T.reserve(length);

	INS_COST = M.INS_COST;
	DEL_COST = M.DEL_COST;
	max_dels = M.max_dels;
	iterations = M.iterations;	

	int i;
	for (i=0; i<length; i++)
	{
		A[i] = M.A[i];
		C[i] = M.C[i];
		G[i] = M.G[i];
		T[i] = M.T[i];
	}
	
	init();
	}*/

	/*position_weight_matrix& position_weight_matrix::operator=(const position_weight_matrix& M)
{
	clear();
	match_reverse_complement = M.match_reverse_complement;
	length = M.length;
	A.reserve(length);
	C.reserve(length);
	G.reserve(length);
	T.reserve(length);

	INS_COST = M.INS_COST;
	DEL_COST = M.DEL_COST;
	max_dels = M.max_dels;
	DNA = M.DNA;
	iterations = M.iterations;	

	int i;
	for (i=0; i<length; i++)
	{
		A[i] = M.A[i];
		C[i] = M.C[i];
		G[i] = M.G[i];
		T[i] = M.T[i];
	}

	init();
	}*/

// ==============================================================
void position_weight_matrix::init()
// ==============================================================
{
	information_content.reserve(length);
	min_freq.reserve(length);
	max_freq.reserve(length);
	freq.clear();
	freq.resize(CHAR_MAX);
	int i;
	for(i=0; i<CHAR_MAX; i++)
		freq[i].clear();
	compl_table.resize(256);
	for(i=0; i<256; i++)
		compl_table[i]=(int)'*';
	
	compl_table[(int)'A']='T';
	compl_table[(int)'C']='G';
	compl_table[(int)'G']='C';
	compl_table[(int)'T']='A';

	reinit();
}

// ==============================================================
void position_weight_matrix::reinit()
// ==============================================================
{

	int i;
	min_freq.resize(length);
	max_freq.resize(length);
	for (i=0; i<length; i++)
	{
		min_freq[i] = _min(A[i],C[i],G[i],T[i]);
		max_freq[i] = _max(A[i],C[i],G[i],T[i]);
	}
	freq.resize(CHAR_MAX);
	freq[(int)'A'] = A;
	freq[(int)'C'] = C;
	freq[(int)'G'] = G;
	freq[(int)'T'] = T;
	freq[(int)'U'] = T;
	freq[(int)'a'] = A;
	freq[(int)'c'] = C;
	freq[(int)'g'] = G;
	freq[(int)'t'] = T;
	freq[(int)'U'] = T;

	rev_freq.resize(CHAR_MAX);
	rev_freq['A'].resize(length);
	rev_freq['C'].resize(length);
	rev_freq['G'].resize(length);
	rev_freq['T'].resize(length);
	rev_freq['U'].resize(length);
	rev_information_content.resize(length);
	rev_min_freq.resize(length);
	rev_max_freq.resize(length);

	compute_information_content();

	for (i=0; i<length; i++)
	{
		rev_freq['A'][i]=freq['T'][length-i-1];
		rev_freq['C'][i]=freq['G'][length-i-1];
		rev_freq['G'][i]=freq['C'][length-i-1];
		rev_freq['T'][i]=freq['A'][length-i-1];
		rev_freq['U'][i]=freq['A'][length-i-1];
		rev_information_content[i]=information_content[length-i-1];
		rev_min_freq[i]=min_freq[length-i-1];
		rev_max_freq[i]=max_freq[length-i-1];
	}

	mis_sequence = most_informative_string(.5);
	rev_mis_sequence = rev_most_informative_string(.5);

	Min=0.;
	Max=0.;
	Rest_Bound.resize(length+1);
	Rest_Bound[length]=0.;
	int j;
	for (j=length-1; j>=0; j--)
	{
		Min += information_content[j]*min_freq[j];
		Max += information_content[j]*max_freq[j];
		Rest_Bound[j] = Max;
	}
	
}

// ==============================================================
std::ifstream& operator>>(std::ifstream& ifs, position_weight_matrix& M)
// ==============================================================
{

	int j=0;
	std::string line;
	std::string filename;
	float_list l_A,l_C,l_G,l_T;
	// read A line

	if (ifs.good())
	{
		getline(ifs,line);
		tokenizer tokens(line.c_str());
		tokenizer::iterator tok_it = tokens.begin();
		while (tok_it != tokens.end())
		{
			std::string tok_str = *tok_it;
			if (!is_flt(tok_str.c_str()))
				throw bbq_exception("float value expected while reading matrix file.\n");
			float f_val = atof(tok_str.c_str());
			l_A.push_back(f_val);
			tok_it++;
		}
		M.length = l_A.size();
	} else {
		throw bbq_exception(std::string("error reading matrix file\n"));
	}
	
	float_list::iterator A_it = l_A.begin();
	float_list::iterator A_end = l_A.end();
	M.A.reserve(M.length);
	for (j=0; j<M.length; A_it++, j++)
	{
		M.A.push_back(*A_it);
	}
	// read C line
	if (ifs.good())
	{
		getline(ifs,line);
		tokenizer tokens(line.c_str());
		tokenizer::iterator tok_it = tokens.begin();
		while (tok_it != tokens.end())
		{
			std::string tok_str = *tok_it;
			if (!is_flt(tok_str.c_str()))
				throw bbq_exception("float value expected while reading matrix file.\n");
			float f_val = atof(tok_str.c_str());
			l_C.push_back(f_val);
			tok_it++;
		}
		if (l_C.size()!=M.length)
				throw bbq_exception("error: different line length in pwm file.\n");
	} else {
		throw bbq_exception(std::string("error reading matrix file\n"));
	}
	float_list::iterator C_it = l_C.begin();
	float_list::iterator C_end = l_C.end();
	M.C.reserve(M.length);
	for (j=0; j<M.length; C_it++, j++)
	{
		M.C.push_back(*C_it);
	}

	// read G line
	if (ifs.good())
	{
		getline(ifs,line);
		tokenizer tokens(line.c_str());
		tokenizer::iterator tok_it = tokens.begin();
		while (tok_it != tokens.end())
		{
			std::string tok_str = *tok_it;
			if (!is_flt(tok_str.c_str()))
				throw bbq_exception("float value expected while reading matrix file.\n");
			float f_val = atof(tok_str.c_str());
			l_G.push_back(f_val);
			tok_it++;
		}
		if (l_G.size()<M.length)
				throw bbq_exception("error: different line length in pwm file.\n");
	} else {
		throw bbq_exception(std::string("error reading matrix file\n"));
	}
	float_list::iterator G_it = l_G.begin();
	float_list::iterator G_end = l_G.end();
	M.G.reserve(M.length);
	for (j=0; j<M.length; G_it++, j++)
	{
		M.G.push_back(*G_it);
	}

	// read T line
	if (ifs.good())
	{
		getline(ifs,line);
		tokenizer tokens(line.c_str());
		tokenizer::iterator tok_it = tokens.begin();
		while (tok_it != tokens.end())
		{
			std::string tok_str = *tok_it;
			if (!is_flt(tok_str.c_str()))
				throw bbq_exception("float value expected while reading matrix file.\n");
			float f_val = atof(tok_str.c_str());
			l_T.push_back(f_val);
			tok_it++;
		}
		if (l_T.size()<M.length)
			throw bbq_exception("error: different line length in pwm file.\n");
	} else {
		throw bbq_exception(std::string("error reading matrix file\n"));
	}
	float_list::iterator T_it = l_T.begin();
	float_list::iterator T_end = l_T.end();
	M.T.reserve(M.length);
	for (j=0; j<M.length; T_it++, j++)
	{
		M.T.push_back(*T_it);
	}

	if (M.match_mode==0 || M.match_mode==1)
	{
		for (j=0; j<M.length; j++)
		{
			
			float sum = M.A[j]+M.C[j]+M.G[j]+M.T[j];
			if (sum<=0.)
				sum = 1.;
			
			M.A[j] /= sum;
			M.C[j] /= sum;
			M.G[j] /= sum;
			M.T[j] /= sum;
		}
	}

	M.init();
	
	return ifs;
}

// ==============================================================
void position_weight_matrix::clear()
// ==============================================================
{
	
	A.clear();
	C.clear();
	G.clear();
	T.clear();
	information_content.clear();
	freq.clear();
	
}

// ==============================================================
void position_weight_matrix::set_DNA()
// ==============================================================
{
	DNA=true;
}
	
// ==============================================================
void position_weight_matrix::set_RNA()
// ==============================================================
{
	DNA=false;
}
	
// ==============================================================
void position_weight_matrix::set_reverse_complement(bool rev_comp)
// ==============================================================
{

	match_reverse_complement = rev_comp;
	// I.: invert order of rows
	/*int l_2 = length / 2;
	int j;
	for (j=0; j<l_2; j++)
	{
		
		float tmp = A[j];
		A[j] = A[length-j-1];
		A[length-j-1] = tmp;
		tmp = C[j];
		C[j] = C[length-j-1];
		C[length-j-1] = tmp;
		tmp = G[j];
		G[j] = G[length-j-1];
		G[length-j-1] = tmp;
		tmp = T[j];
		T[j] = T[length-j-1];
		T[length-j-1] = tmp;
		
	}

	// II.: swap column order
	for (j=0; j<length; j++)
	{
		
		float tmp = A[j];
		A[j] = T[j];
		T[j] = tmp;	
		tmp = C[j];
		C[j] = G[j];
		G[j] = tmp;	
		
	}

	reinit();*/
	
}

// ==============================================================
bool position_weight_matrix::read(const char* fname)
// ==============================================================
{
	
	int j=0;
	std::ifstream ifs(fname);

	std::string line;
	std::string filename;
	float_list l_A,l_C,l_G,l_T;
	// read A line

	if (ifs.good())
	{
		getline(ifs,line);
		tokenizer tokens(line.c_str());
		tokenizer::iterator tok_it = tokens.begin();
		while (tok_it != tokens.end())
		{
			std::string tok_str = *tok_it;
			if (!is_flt(tok_str.c_str()))
				throw bbq_exception("float value expected while reading matrix file.\n");
			float f_val = atof(tok_str.c_str());
			l_A.push_back(f_val);
			tok_it++;
		}
		length = l_A.size();
	} else {
		throw bbq_exception(std::string("error reading matrix file ") + std::string(fname) + std::string(".\n"));
	}
	
	float_list::iterator A_it = l_A.begin();
	float_list::iterator A_end = l_A.end();
	A.reserve(length);
	for (j=0; j<length; A_it++, j++)
	{
		A.push_back(*A_it);
	}
	// read C line
	if (ifs.good())
	{
		getline(ifs,line);
		tokenizer tokens(line.c_str());
		tokenizer::iterator tok_it = tokens.begin();
		while (tok_it != tokens.end())
		{
			std::string tok_str = *tok_it;
			if (!is_flt(tok_str.c_str()))
				throw bbq_exception("float value expected while reading matrix file.\n");
			float f_val = atof(tok_str.c_str());
			l_C.push_back(f_val);
			tok_it++;
		}
		if (l_C.size()!=length)
				throw bbq_exception("error: different line length in pwm file.\n");
	} else {
		throw bbq_exception(std::string("error reading matrix file ") + std::string(fname) + std::string(".\n"));
	}
	float_list::iterator C_it = l_C.begin();
	float_list::iterator C_end = l_C.end();
	C.reserve(length);
	for (j=0; j<length; C_it++, j++)
	{
		C.push_back(*C_it);
	}

	// read G line
	if (ifs.good())
	{
		getline(ifs,line);
		tokenizer tokens(line.c_str());
		tokenizer::iterator tok_it = tokens.begin();
		while (tok_it != tokens.end())
		{
			std::string tok_str = *tok_it;
			if (!is_flt(tok_str.c_str()))
				throw bbq_exception("float value expected while reading matrix file.\n");
			float f_val = atof(tok_str.c_str());
			l_G.push_back(f_val);
			tok_it++;
		}
		if (l_G.size()<length)
				throw bbq_exception("error: different line length in pwm file.\n");
	} else {
		throw bbq_exception(std::string("error reading matrix file ") + std::string(fname) + std::string(".\n"));
	}
	float_list::iterator G_it = l_G.begin();
	float_list::iterator G_end = l_G.end();
	G.reserve(length);
	for (j=0; j<length; G_it++, j++)
	{
		G.push_back(*G_it);
	}

	// read T line
	if (ifs.good())
	{
		getline(ifs,line);
		tokenizer tokens(line.c_str());
		tokenizer::iterator tok_it = tokens.begin();
		while (tok_it != tokens.end())
		{
			std::string tok_str = *tok_it;
			if (!is_flt(tok_str.c_str()))
				throw bbq_exception("float value expected while reading matrix file.\n");
			float f_val = atof(tok_str.c_str());
			l_T.push_back(f_val);
			tok_it++;
		}
		if (l_T.size()<length)
			throw bbq_exception("error: different line length in pwm file.\n");
	} else {
		throw bbq_exception(std::string("error reading matrix file ") + std::string(fname) + std::string(".\n"));
	}
	float_list::iterator T_it = l_T.begin();
	float_list::iterator T_end = l_T.end();
	T.reserve(length);
	for (j=0; j<length; T_it++, j++)
	{
		T.push_back(*T_it);
	}

	if (match_mode==0 || match_mode==1)
	{
		for (j=0; j<length; j++)
		{
			
			float sum = A[j]+C[j]+G[j]+T[j];
			if (sum<=0.)
				sum = 1.;
			
			A[j] /= sum;
			C[j] /= sum;
			G[j] /= sum;
			T[j] /= sum;
		}
	}

	init();
	
	return true;
}

// ==============================================================
void position_weight_matrix::compute_information_content()
// ==============================================================
{
	int i;
	information_content.resize(length);
	for (i=0; i<length; i++)
	{
		float sum = A[i]+C[i]+G[i]+T[i];
		if (sum<=0.)
			sum = 1.;
		float ic_i  = A[i]==0. ? 0. : (A[i]/sum)*log(4.*A[i]/sum);
		ic_i += C[i]==0. ? 0. : (C[i]/sum)*log(4.*C[i]/sum);
		ic_i += G[i]==0. ? 0. : (G[i]/sum)*log(4.*G[i]/sum);
		ic_i += T[i]==0. ? 0. : (T[i]/sum)*log(4.*T[i]/sum);
		information_content[i] = ic_i;		
	}
}

// ==============================================================
float position_weight_matrix::get_information_content()
// ==============================================================
{
	int i;
	float ic = 0.;
	for (i=0; i<length; i++)
	{
		ic += information_content[i];
	}
	return ic/(2.*(float)length);
}

// ==============================================================
float position_weight_matrix::get_p_value(double observed_score, const double_vector& bg_freq)
// ==============================================================
{
#ifdef USE_GSL

	if (bg_freq.size()<4)
		return 0.;
	
	// compute expected value and standard deviation of the matrix
	double mu=0.;
	double sigma_sq=0.;
	int i,j;
	for (i=0; i<length; ++i)
	{
		float sum_i = (A[i]+C[i]+G[i]+T[i]);
		if (sum_i<=0.)
			return 0.;
		float A_i = A[i]/sum_i;
		float C_i = C[i]/sum_i;
		float G_i = G[i]/sum_i;
		float T_i = T[i]/sum_i;
		double mu_sum_i = A_i*bg_freq[0];
		mu_sum_i += C_i*bg_freq[1];
		mu_sum_i += G_i*bg_freq[2];
		mu_sum_i += T_i*bg_freq[3];
		mu += information_content[i]*mu_sum_i;

		double sigma_sum_i = A_i*A_i*bg_freq[0];
		sigma_sum_i += C_i*C_i*bg_freq[1];
		sigma_sum_i += G_i*G_i*bg_freq[2];
		sigma_sum_i += T_i*T_i*bg_freq[3];
		sigma_sq += information_content[i]*information_content[i]*sigma_sum_i;
		
	}
	double sigma = sigma_sq-mu*mu;
	if (sigma<0)
		sigma=-sigma;
	sigma=sqrt(sigma);
	double p_val = gsl_cdf_gaussian_P(mu-observed_score, sigma);
	if (max_dels==0)
		return p_val;
	// compute binomial coeff. by simple dynamic programming
	double length_choose_delta;
	double_vector B(length+1);
	B[0]=1.;
	for (i=0; i<=length; ++i)
	{
		B[0]=1.;
		for (i=0; i<=length; ++i)
			B[i] += B[i-1];
	}
	double p = B[max_dels]*p_val;
	if (p<=1.)
		return p;
	return 1.;
	
#else

	// if libgsl is not available, return the match score as a "bogus p-value"
	return observed_score;
#endif
}

// ==============================================================
void position_weight_matrix::set_matches(const char* seq, float threshold)
// ==============================================================
{
	matches.clear();
	if (match_mode==0 || match_mode==1)
		set_naive_matches(seq, threshold);
	else if (match_mode==2)
		set_MATCH_matches(seq, threshold);
}

// ==============================================================
void position_weight_matrix::set_MATCH_matches(const char* seq, float threshold)
// ==============================================================
{
	matches.clear();
	max_score = 0.;
	max_pos = 0;
	bool curr_rev;
	max_score_rev = false;
	
	int slen = strlen(seq);
	int max_idx = slen-length+1;
	int i,j;
	for (i=0; i<max_idx; i++)
	{
		float score,rev_score,Min,Max,Current,rev_Current;
		score = 0.;
		rev_score = 0.;
		Min = 0.;
		Max = 0.;
		Current = 0.;
		rev_Current = 0.;
		for (j=0; j<length; j++)
		{
			Min += information_content[j]*min_freq[j];
			Max += information_content[j]*max_freq[j];
			if (!freq[(int)seq[i+j]].empty())
			{
				Current += information_content[j]*freq[seq[i+j]][j];
				rev_Current += information_content[length-j-1]*freq[compl_table[seq[i+j]]][length-j-1];
			}
		}
		score = (Current-Min) / (Max-Min);
		rev_score = (rev_Current-Min) / (Max-Min);
		if (score>=threshold)
			// biologists' sequences start with one, so we add 1 to our position...
			matches.push_back(match(i+1,score,false));
		if (rev_score>=threshold && match_reverse_complement)
			matches.push_back(match(i+1,rev_score,true));

		if (score>max_score)
		{
			max_score_rev = false;
			max_score = score;
			max_pos = i;
		}

		if (rev_score>max_score && match_reverse_complement)
		{
			max_score_rev = true;
			max_score = rev_score;
			max_pos = i;
		}
	}
}

// ==============================================================
void position_weight_matrix::set_indel_costs(float ins, float del, int m_d)
// ==============================================================
{
	INS_COST=ins;
	DEL_COST=del;
	max_dels=m_d;
}


// ==============================================================
const std::string& position_weight_matrix::get_match_sequence() const
// ==============================================================
{
	return curr_match_seq;
}		

/// check whether the matrix match score with sequence \c seq exceeds
/// threshold theta. 
// ==============================================================
bool position_weight_matrix::get_plain_match_score(const char* seq, float theta)
// ==============================================================
{
	int j;
	float score,rev_score,Current=0.;
	for (j=0; j<length; j++)
	{
		if (!freq[(int)seq[j]].empty())
		{
			Current += information_content[j]*freq[seq[j]][j];
			// check whether score theta can still be achieved through
			// the upper bound table computed in reinit()
			if (Current+Rest_Bound[j+1]/(Max-Min)<theta)
				return false;
		}
	}
	score=(Current)/(Max-Min);
	return score>=theta;
}

/// check whether the matrix match score with sequence \c seq exceeds
/// threshold theta. 
// ==============================================================
float_pair position_weight_matrix::get_plain_match_score(const char* seq)
// ==============================================================
{
	int j;
	float score,rev_score,Current=0.,rev_Current=0.;
	for (j=0; j<length; j++)
	{
		if (!freq[(int)seq[j]].empty())
		{
			rev_Current += information_content[length-j-1]*freq[compl_table[seq[j]]][length-j-1];
			Current += information_content[j]*freq[seq[j]][j];
		}
	}
	score=((Current-Min) / (Max-Min));
	rev_score=((rev_Current-Min) / (Max-Min));
	return float_pair(score, rev_score);
}

// ==============================================================
float_pair position_weight_matrix::get_frac_score(const char* seq)
// ==============================================================
{
	if (max_dels==0)
		return get_plain_match_score(seq);
	//std::cerr<<"computing fwd score...\n";
	float fwd_score = get_fwd_frac_score(seq);
	//std::cerr<<"computing bwd score...\n";
	float rev_score = get_rev_frac_score(seq);
	return float_pair(fwd_score,rev_score);
}
 
// ==============================================================
bool position_weight_matrix::match_frac_score(const char* seq, float threshold)
// ==============================================================
{
	if (max_dels==0)
	{
		if (get_plain_match_score(seq,threshold))
		{
			//curr_match_seq = most_informative_string(.5);
			curr_match_seq = mis_sequence;
			return true;
		} else {
			return false;
		}
	}

	return fractional_programming_score 
		(threshold,seq,
			freq,information_content,min_freq,max_freq,true);
}

// ==============================================================
float_pair position_weight_matrix::get_match_score(const char* seq)
// ==============================================================
{
	if (max_dels==0)
		return get_plain_match_score(seq);
	float fwd_score = get_fwd_match_score(seq);
	float rev_score = get_rev_match_score(seq);
	return float_pair(fwd_score,rev_score);
}

// ==============================================================
bool position_weight_matrix::fractional_programming_score(
	float theta,
	const char* seq,
	const position_weight_matrix::float_array& fb_freq,
	const position_weight_matrix::float_vector& fb_ic,
	const position_weight_matrix::float_vector& fb_min_freq,
	const position_weight_matrix::float_vector& fb_max_freq,
	bool trace
	)
// ==============================================================
{
	int i,j,d;
	int N = strlen(seq);
	int diff=N-length;
	if (diff>max_dels || diff<-max_dels)
		return false;

	// matrix for fractional programming
	float_table FTABLE(length+1);
	for (i=0; i<=length; i++)
	{
		FTABLE[i].resize(N+1);
		for (j=0; j<=N; j++)
		{
			FTABLE[i][j].resize(max_dels+1);
		}
	}

	float_vector Dyn_Rest_Bound;
	float Bound=0.;
	Dyn_Rest_Bound.resize(length+1);
	Dyn_Rest_Bound[length]=0.;
	for (j=length-1; j>=0; j--)
	{
		Bound += fb_ic[j]*(fb_max_freq[j]-fb_min_freq[j]);
		Dyn_Rest_Bound[j] = Bound;
	}

	FTABLE[0][0][0]=0.;
	//std::cerr<<"["<<FTABLE[0][0][0]<<";";
	for (d=1; d<=max_dels; d++)
	{
		FTABLE[0][0][d] = -FLT_MAX;
		//std::cerr<<FTABLE[0][0][d]<<";";
	}
	//std::cerr<<"]\t";
	// initialize 1st row and first column of FTABLE
	for (i=1; i<=length; i++)
	{
		FTABLE[i][0][0] = -(float)i*INS_COST;
		for (d=1; d<=max_dels; d++)
			FTABLE[i][0][d] = FTABLE[i-1][0][d-1]-fb_ic[i-1]*DEL_COST;
	}
	for (j=1; j<=N; j++)
	{
		FTABLE[0][j][0]=-FLT_MAX;
		for (d=1; d<=max_dels; d++)
			FTABLE[0][j][d] = -FLT_MAX;
		////std::cerr<<FTABLE[0][j]<<"\t";
	}
	//std::cerr<<"\n";
	// Do Dynamic Programming
	for (i=1; i<=length; i++)
	{
		////std::cerr<<FTABLE[i][0]<<"\t";

		// keep track of maximum column scores for eventually terminate
		// DP early according to the dynamic programming bounds
		float col_max = -FLT_MAX;
		int j_start = i-max_dels;
		int j_end = i+max_dels;
		if (j_start<1) j_start=1;
		if (j_end>N) j_end=N;
		if (j_start>0) 
			for (d=0; d<=max_dels; d++)
				FTABLE[i][j_start-1][d] = -FLT_MAX;
		for (j=j_start; j<=j_end; j++)
		{
			//std::cerr<<"[";
			for (d=0; d<=max_dels; d++)
			{
				float m_score =
					FTABLE[i-1][j-1][d]
					+ fb_ic[i-1]*(fb_freq[seq[j-1]][i-1]
				- (1.-theta)*fb_min_freq[i-1]
						- theta*fb_max_freq[i-1]);
				float i_score = FTABLE[i][j-1][d] - INS_COST;
				float d_score;
				if (d>0)
					d_score = FTABLE[i-1][j][d-1] - fb_ic[i-1]*DEL_COST;
				else
					d_score = -FLT_MAX;
				float max_score = max<float>(m_score,i_score,d_score);
				FTABLE[i][j][d] = max_score;
				if (max_score>col_max)
					col_max=max_score;
				//std::cerr<<FTABLE[i][j][d]<<";";
			}
			if (j_end<N) 
				for (d=0; d<=max_dels; d++)
					FTABLE[i][j_end+1][d] = -FLT_MAX;
				
			//std::cerr<<"]\t";
		}
		if (col_max+(1.-theta)*Dyn_Rest_Bound[i]<0.)
			return false;
		//std::cerr<<"\n";
	}
	//std::cerr<<"\n";
	
	float max_total_score = -FLT_MAX;
	int best_d = -1;
	for (d=0; d<=max_dels; d++)
	{
		float score_d = FTABLE[length][N][d];
		if (FTABLE[length][N][d]>=max_total_score)
		{
			max_total_score = score_d;
			best_d = d;
		}
	}
	if (max_total_score<0.)
		return false;
	//std::cerr<<"best_d = "<<best_d<<"\n";
	//std::cerr<<"best_score = "<<max_total_score<<"\n";
	if (trace)
	{
		trace_back(
			theta,seq,
			fb_freq,fb_ic,fb_min_freq,fb_max_freq,
			best_d,FTABLE
			);
		//alignment_path::iterator a_it = curr_aln_path.begin();
		//alignment_path::iterator a_end = curr_aln_path.end();
		alignment_sequence::iterator a_it = curr_aln_seq.begin();
		alignment_sequence::iterator a_end = curr_aln_seq.end();
		for (; a_it!=a_end; ++a_it)
		{
			//std::cerr<<"("<<a_it->first<<","<<a_it->second<<"):";
		}
		//std::cerr<<"\n";
	}
	return true;
}

// ==============================================================
void position_weight_matrix::trace_back(
			float theta,
			const char* seq,
			const float_array& fb_freq,
			const float_vector& fb_ic,
			const float_vector& fb_min_freq,
			const float_vector& fb_max_freq,
			int best_d,
			const float_table& FTABLE
			)
// ==============================================================
{

	int N = strlen(seq);
	int i = length;
	int j = N;
	int d = best_d;
	
	curr_aln_path.clear();
	curr_aln_seq.clear();
	curr_match_seq.clear();
	
	while (i>0 && j>0)
	{
		curr_aln_path.push_front(int_pair(i,j));
		float m_score =
			  FTABLE[i-1][j-1][d]
			+ fb_ic[i-1]*(fb_freq[seq[j-1]][i-1]
			- (1.-theta)*fb_min_freq[i-1]
			- theta*fb_max_freq[i-1]);
		float i_score =  FTABLE[i][j-1][d] - INS_COST;
		float d_score =  FTABLE[i-1][j][d-1] - fb_ic[i-1]*DEL_COST;
		if (FTABLE[i][j][d]==m_score && i>0 && j>0)
		{
			curr_aln_seq.push_front(char_pair(seq[j-1],most_informative_character(i-1,.5)));
			curr_match_seq.append(1,most_informative_character(i-1,.5));
			//std::cerr<<"("<<seq[j]<<","<<most_informative_character(i)<<"):";
			i--; j--;
		} else if (FTABLE[i][j][d]==i_score && i>0) {
			curr_aln_seq.push_front(char_pair('-',most_informative_character(i-1,.5)));
			//std::cerr<<"("<<"-"<<","<<most_informative_character(i)<<"):";
			j--;
		} else if (FTABLE[i][j][d]==d_score && j>0 && d>0) {
			curr_aln_seq.push_front(char_pair(seq[j-1],'-'));
			curr_match_seq.append(1,'-');
			//std::cerr<<"("<<seq[j]<<","<<"-"<<"):";
			i--; d--;
		} else {
			curr_aln_seq.push_front(char_pair('?','?'));
			//std::cerr<<"(?,?):";
			break;
		}
	}
	while (i>0)
	{
		curr_aln_path.push_front(int_pair(i,j));
		curr_match_seq.append(1,most_informative_character(i-1,.5));
		curr_aln_seq.push_front(char_pair('-',most_informative_character(i-1,.5)));
		//std::cerr<<"("<<"-"<<","<<most_informative_character(i)<<"):";
		i--;
	}
	while (j>0)
	{
		curr_aln_path.push_front(int_pair(i,j));
		curr_aln_seq.push_front(char_pair(seq[j-1],'-'));
		curr_match_seq.append(1,'-');
		//std::cerr<<"("<<seq[j]<<","<<"-"<<"):";
		j--;
	}
	//std::cerr<<"\n";*/
	reverse(curr_match_seq.begin(),curr_match_seq.end());
}

/// implementation using fractional programming
// ==============================================================
float position_weight_matrix::get_fwd_frac_score(const char* seq)
// ==============================================================
{
	float theta,offset;
	
	if (fractional_programming_score (
			 1.,seq,freq,information_content,min_freq,max_freq,true
			 ))
		return 1.;
	if (!fractional_programming_score (
			 0.,seq,freq,information_content,min_freq,max_freq,true
			 ))
		return 0.;
	
	int i;
	float best_theta=0.;
	theta=.5;
	offset=.25;
	for (i=0; i<iterations; i++)
	{
		//std::cerr<<"Iteration "<<i<<"; theta="<<theta<<"\n";
		if (fractional_programming_score(
				 theta,seq,freq,information_content,min_freq,max_freq,false
				 ))
		{
			best_theta=theta;
			theta+=offset;
		} else {
			theta-=offset;
		}
		offset /= 2.;
	}
	// do the dynamic programming one more time for the best theta with
	// positive solution in order to reconstruct a valid alignment path
	fractional_programming_score(
		theta,seq,freq,information_content,min_freq,max_freq,true
		);
	return best_theta;
}


// ==============================================================
float position_weight_matrix::get_rev_frac_score(const char* seq)
// ==============================================================
{
	float theta,offset;
	
	if (fractional_programming_score (
				 1.,seq,
				 rev_freq,rev_information_content,rev_min_freq,rev_max_freq,true
			 ))
		return 1.;
	if (!fractional_programming_score (
			 0.,seq,
				 rev_freq,rev_information_content,rev_min_freq,rev_max_freq,true
			 ))
		return 0.;
	
	int i;
	float best_theta=0.;
	theta=.5;
	offset=.25;
	for (i=0; i<iterations; i++)
	{
		//std::cerr<<"Iteration "<<i<<"; theta="<<theta<<"\n";
		if (fractional_programming_score(
				 theta,seq,
				 rev_freq,rev_information_content,rev_min_freq,rev_max_freq,false
				 ))
		{
			best_theta=theta;
			theta+=offset;
		} else {
			theta-=offset;
		}
		offset /= 2.;
	}
	// do the dynamic programming one more time for the best theta with
	// positive solution in order to reconstruct a valid alignment path
	fractional_programming_score(
		theta,seq,
		rev_freq,rev_information_content,rev_min_freq,rev_max_freq,true
		);
	return best_theta;

	//return get_fb_match_score(
	//	seq,max_dels,rev_freq,rev_information_content,rev_min_freq,rev_max_freq
	//	);
}

// ==============================================================	
float position_weight_matrix::get_rev_match_score(const char* seq)
// ==============================================================
{
	return get_fb_match_score(
		seq,rev_freq,rev_information_content,rev_min_freq,rev_max_freq
		);
}

// ==============================================================
float position_weight_matrix::get_fwd_match_score(const char* seq)
// ==============================================================
{
	return get_fb_match_score(
		seq,freq,information_content,min_freq,max_freq
		);
}
	
// ==============================================================
float position_weight_matrix::get_fb_match_score(
	const char* seq,
	const float_array& fb_freq,
	const float_vector& fb_ic,
	const float_vector& fb_min_freq,
	const float_vector& fb_max_freq
	)
// ==============================================================
{
	
	int j,jj;
	int N = strlen(seq);
	if (N>length)
		N = length;
	float m_score, d_score, m_Min, d_Min, m_Max, d_Max, m_Current, d_Current;
	
	float score,rev_score,Min=0.,Max=0.,Current=0.,rev_Current=0.;

	float NEG_INF;
	if (std::numeric_limits<float>::has_infinity)
		NEG_INF=-std::numeric_limits<float>::infinity();
	else
		NEG_INF=-std::numeric_limits<float>::max();
	// start w/ first column filled w/ all 0 for counting deletions
	float_float_vector Min_column_1(max_dels+1,float_vector(N));
	float_float_vector Min_column_2(max_dels+1,float_vector(N));
	float_float_vector Max_column_1(max_dels+1,float_vector(N));
	float_float_vector Max_column_2(max_dels+1,float_vector(N));
	float_float_vector SCurrent_column_1(max_dels+1,float_vector(N));
	float_float_vector SCurrent_column_2(max_dels+1,float_vector(N));
	float_float_vector score_column_1(max_dels+1,float_vector(N));
	float_float_vector score_column_2(max_dels+1,float_vector(N));
	bool_bool_vector gap_column_1(max_dels+1,bool_vector(N));	
	bool_bool_vector gap_column_2(max_dels+1,bool_vector(N));
	
	for (int s=0;s<=max_dels;s++)
	{
		for (int i=0;i<N;i++)
		{
			score_column_1[s][i]=NEG_INF; score_column_2[s][i]=NEG_INF;
			Min_column_1[s][i]=NEG_INF; Min_column_2[s][i]=NEG_INF;
			Max_column_1[s][i]=NEG_INF; Max_column_2[s][i]=NEG_INF;
			SCurrent_column_1[s][i]=NEG_INF; SCurrent_column_2[s][i]=NEG_INF;
			gap_column_1[s][i]=false; gap_column_2[s][i]=false;
		}
	}
	/*score_column_1[0][0]=0.; score_column_2[0][0]=0.;
	Min_column_1[0][0]=0.; Min_column_2[0][0]=0.;
	Max_column_1[0][0]=0.; Max_column_2[0][0]=0.;
	SCurrent_column_1[0][0]=0.; SCurrent_column_2[0][0]=0.;*/
	//column_1[0][0]=0; column_2[0][0]=0;
	float_float_vector* prev_Min_column = &Min_column_1;
	float_float_vector* curr_Min_column = &Min_column_2;
	float_float_vector* prev_Max_column = &Max_column_1;
	float_float_vector* curr_Max_column = &Max_column_2;
	float_float_vector* prev_SCurrent_column = &SCurrent_column_1;
	float_float_vector* curr_SCurrent_column = &SCurrent_column_2;
	float_float_vector* prev_score_column = &score_column_1;
	float_float_vector* curr_score_column = &score_column_2;
	float_float_vector* tmp_column;
	bool_bool_vector* curr_gap_column=&gap_column_1;
	bool_bool_vector* prev_gap_column=&gap_column_2;
	bool_bool_vector* tmp_gap_column;
	
		
	for (jj=0; (jj<length); jj++)
	{
		int kk,ss;
		int min_kk=jj-max_dels-1;
		if (min_kk<0)
			min_kk=0;
		min_kk=0;
		//int max_kk=jj;
		//if (max_kk>N)
		int max_kk=N;
		for (kk=min_kk; (kk<max_kk); kk++)
		{
			int xxx=1;
			for (ss=0; (ss<=max_dels); ss++)
			{
				if (kk>0)
				{
					m_Min = (*prev_Min_column)[ss][kk-1] + fb_ic[jj]*fb_min_freq[jj];
					m_Max = (*prev_Max_column)[ss][kk-1] + fb_ic[jj]*fb_max_freq[jj];
					if (!fb_freq[(int)seq[kk]].empty())
					{
						m_Current = (*prev_SCurrent_column)[ss][kk-1] + fb_ic[jj]*fb_freq[seq[kk]][jj];
					}
					if ((*prev_gap_column)[ss][kk-1] && jj>0)
					{
						m_Min += fb_ic[jj-1]*fb_min_freq[jj-1];
						m_Max += fb_ic[jj-1]*fb_max_freq[jj-1];
						m_Current += fb_ic[jj-1]*fb_freq[seq[kk-1]][jj-1];
					}
					if (m_Max-m_Min>0.)
						m_score = (m_Current-m_Min) / (m_Max-m_Min);
					else
						m_score = 0.;
					//rev_score = (rev_Current-Min) / (Max-Min);
				} else {
					if (jj>0 || ss>0) {
						m_score=NEG_INF;
						m_Min=NEG_INF;
						m_Max=NEG_INF;
					 	m_Current=NEG_INF;
					} else {
						m_Current=fb_ic[jj]*fb_freq[seq[kk]][jj];;
						m_Min=fb_ic[jj]*fb_min_freq[jj];;
						m_Max=fb_ic[jj]*fb_max_freq[jj];;
						if (m_Max-m_Min>0.)
							m_score = (m_Current-m_Min) / (m_Max-m_Min);
						else
							m_score = 0.;
						}
				}
				// compute costs for (mis)match and deletion
				
				//if (m_cost>=0 && match_table[index_table[(int)s_j[jj]]][index_table[T[ii+kk]]]==0)
				//if (m_cost>=0 && s_j[jj]!=T[ii+kk])
				//	m_cost++;
				if (ss>0)
				{
					d_Min = (*prev_Min_column)[ss-1][kk];
					d_Max = (*prev_Max_column)[ss-1][kk];
					d_Current = (*prev_SCurrent_column)[ss-1][kk];
					if (!(*prev_gap_column)[ss-1][kk] && jj>0)
					{
						d_Min -= fb_ic[jj-1]*fb_min_freq[jj-1];
						d_Max -= fb_ic[jj-1]*fb_max_freq[jj-1];
						d_Current -= fb_ic[jj-1]*fb_freq[seq[kk]][jj-1];
					}
					d_score = (*prev_score_column)[ss-1][kk];
				} else {
					d_score=NEG_INF;
					d_Current=NEG_INF;
					d_Min=NEG_INF;
					d_Max=NEG_INF;
				}
				if (m_score>=d_score)
				{
					(*curr_gap_column)[ss][kk]=false;
					(*curr_score_column)[ss][kk] = m_score;
					(*curr_Min_column)[ss][kk] = m_Min;
					(*curr_Max_column)[ss][kk] = m_Max;
					(*curr_SCurrent_column)[ss][kk] = m_Current;
				} else {
					(*curr_gap_column)[ss][kk]=true;
					(*curr_score_column)[ss][kk] = d_score;
					(*curr_Min_column)[ss][kk] = d_Min;
					(*curr_Max_column)[ss][kk] = d_Max;
					(*curr_SCurrent_column)[ss][kk] = d_Current;
				}
			}

		}
		//std::cout<<"\ndel column "<<jj<<": \n";
		/*for (int yy=0; yy<N; yy++)
		{
		  std::cout<<"(";
		  for (int s=0; s<=max_dels; s++) std::cout<<(*prev_score_column)[s][yy]<<" ";
		  std::cout<<")";
		}
		std::cout<<"\n"; std::cout.flush();*/
		// swap columns
		tmp_column = curr_Min_column;
		curr_Min_column = prev_Min_column;
		prev_Min_column = tmp_column;
		tmp_column = curr_Max_column;
		curr_Max_column = prev_Max_column;
		prev_Max_column = tmp_column;
		tmp_column = curr_SCurrent_column;
		curr_SCurrent_column = prev_SCurrent_column;
		prev_SCurrent_column = tmp_column;
		tmp_column = curr_score_column;
		curr_score_column = prev_score_column;
		prev_score_column = tmp_column;
		tmp_gap_column = curr_gap_column;
		curr_gap_column = prev_gap_column;
		prev_gap_column=tmp_gap_column;
	}
	/*for (int yy=0; yy<N; yy++)
	{
		std::cout<<"(";
		for (int s=0; s<=max_dels; s++) std::cout<<(*prev_score_column)[s][yy]<<" ";
		std::cout<<")";
	}
	std::cout<<"\n"; std::cout.flush();*/
	float max_mm=NEG_INF;
	int min_kk = length-max_dels-1;
	if (min_kk<0)
		min_kk=0;
	int max_kk = length;
	if (max_kk>N)
		max_kk=N;
	for (int kk=min_kk; kk<max_kk; kk++)
	{
		int this_d = length-kk-1;
		float this_mm;
		if (this_d>=0 && this_d<=max_dels)
			this_mm=(*prev_score_column)[this_d][kk];
		else
			this_mm = NEG_INF;
		if (this_mm>max_mm)
		{
			max_mm=this_mm;
		}
	}
	if (max_mm==NEG_INF)
		max_mm=0.;
	return max_mm;

}

/// yields the 'most informative sequence' (Freyhult et
/// al 2004)  of the matrix, elements in columns with
/// frequency greater than the background frequency are projected into
/// iupac notation. Columns where gaps are over-represented are in
/// lower case.
// ==============================================================
std::string position_weight_matrix::most_informative_string(float GC_ratio) const
// ==============================================================
{
	std::string result(length,' ');
	int i;
	for (i=0; i<length; i++)
	{
		result[i]=most_informative_character(i,GC_ratio);
	}
	return result;
}

/// yields the 'most informative sequence' (Freyhult et
/// al 2004) for column i of the matrix, elements in columns with
/// frequency greater than the background frequency are projected into
/// iupac notation. Columns where gaps are over-represented are in
/// lower case.
// ==============================================================
char position_weight_matrix::most_informative_character(int i, float GC_ratio) const
// ==============================================================
{
	std::string IUP;
	if (DNA)
		IUP = ".ACMGRSVTWYHKDBN";
	else
		IUP = ".ACMGRSVUWYHKDBN";

	if (i>=length || i<0)
		return '.';
	
	float GC_freq = GC_ratio/2.;
	float AT_freq = .5-GC_freq;
	int code = 0;
	float sum_i = (float)(A[i]+C[i]+G[i]+T[i]);
	
	if (sum_i<=0.)
		return '.';
	
	if (((float)T[i])/sum_i>=AT_freq)
		code++;
	code<<=1;
	if (((float)G[i])/sum_i>=GC_freq)
		code++;
	code<<=1;
	if (((float)C[i])/sum_i>=GC_freq)
		code++;
	code<<=1;
	if (((float)A[i])/sum_i>=AT_freq)
		code++;
	return IUP[code];
}

/// yields the 'most informative sequence' (Freyhult et
/// al 2004)  of the matrix, elements in columns with
/// frequency greater than the background frequency are projected into
/// iupac notation. Columns where gaps are over-represented are in
/// lower case.
// ==============================================================
std::string position_weight_matrix::rev_most_informative_string(float GC_ratio) const
// ==============================================================
{
	std::string result(length,' ');
	int i;
	for (i=0; i<length; i++)
	{
		result[i]=rev_most_informative_character(i,GC_ratio);
	}
	return result;
}

/// yields the 'most informative sequence' (Freyhult et
/// al 2004) for column i of the matrix, elements in columns with
/// frequency greater than the background frequency are projected into
/// iupac notation. Columns where gaps are over-represented are in
/// lower case.
// ==============================================================
char position_weight_matrix::rev_most_informative_character(int i, float GC_ratio) const
// ==============================================================
{
	std::string IUP;
	if (DNA)
		IUP = ".ACMGRSVTWYHKDBN";
	else
		IUP = ".ACMGRSVUWYHKDBN";

	if (i>=length || i<0)
		return '.';

	int ii = length-i-1;
	
	float GC_freq = GC_ratio/2.;
	float AT_freq = .5-GC_freq;
	int code = 0;
	float sum_i = (float)(A[ii]+C[ii]+G[ii]+T[ii]);
	
	if (sum_i<=0.)
		return '.';
	
	if (((float)A[ii])/sum_i>=AT_freq)
		code++;
	code<<=1;
	if (((float)C[ii])/sum_i>=GC_freq)
		code++;
	code<<=1;
	if (((float)G[ii])/sum_i>=GC_freq)
		code++;
	code<<=1;
	if (((float)T[ii])/sum_i>=AT_freq)
		code++;
	return IUP[code];
}


// ==============================================================
void position_weight_matrix::set_naive_matches(const char* seq, float threshold)
// ==============================================================
{

	matches.clear();
	float max_score = 0.;
	int max_pos = 0;
	
	int slen = strlen(seq);
	int max_idx = slen-length;
	int i,j;
	for (i=0; i<max_idx; i++)
	{
		float score;
		if (match_mode==0)
			score = 0.;
		else
			score = 1.;

		for (j=0; j<length; j++)
		{
			float new_score;
			switch (seq[i+j])
			{
				case 'a':
				case 'A': 
					new_score = A[j];
					break;
				case 'c':
				case 'C':
					new_score = C[j];
					break;
				case 'g':
				case 'G':
					new_score = G[j];
					break;
				case 't':
				case 'T':
					new_score = T[j];
					break;
				default:
					new_score = 0.;
					break;
			}
			if (match_mode==0)
				score += new_score;
			else
				score *= new_score;				
		}
		if (score>=threshold)
			matches.push_back(match(i,score,false));
		if (score>=max_score)
		{
			max_score = score;
			max_pos = i;
		}
	}
}

// ==============================================================
position_weight_matrix::match_iterator position_weight_matrix::begin()
// ==============================================================
{
	return matches.begin();
}

// ==============================================================
position_weight_matrix::match_iterator position_weight_matrix::end()
// ==============================================================
{
	return matches.end();
}
	
// ==============================================================
position_weight_matrix::match_list position_weight_matrix::get_match_list()
// ==============================================================
{
	return matches;
}

// ==============================================================
std::ostream& operator<<(std::ostream& ofs, const position_weight_matrix& M)
// ==============================================================
{
	int l = M.length;
	int j;
	
	for (j=0; j<l; j++)
		ofs<<std::setw(4)<<M.A[j]<<" ";
	ofs<<std::endl;
	for (j=0; j<l; j++)
		ofs<<std::setw(4)<<M.C[j]<<" ";
	ofs<<std::endl;
	for (j=0; j<l; j++)
		ofs<<std::setw(4)<<M.G[j]<<" ";
	ofs<<std::setw(4)<<std::endl;
	for (j=0; j<l; j++)
		ofs<<std::setw(4)<<M.T[j]<<" ";
	ofs<<std::endl;
	return ofs;
}

}

