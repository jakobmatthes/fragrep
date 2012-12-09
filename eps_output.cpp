#include "eps_output.h"
#include <sequtil/bbq_util.h>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include "math.h"
#include <algorithm>

namespace fragrep
{

float left_offset;
	
/// formatting information for eps output of matrix patterns
float rewind_length;

// ==============================================================
void print_help_message()
// ==============================================================
{
	std::cout<<"\nusage: pattern2eps [filename]\n\n";
	std::cout<<"Creates for a fragrep pattern, typically obtained\n";
	std::cout<<"using aln2pattern on a multiple sequence alignment.\n";
}


/// write an eps-formatted weblogo representation of matrix,
/// associated as matrix number \c index, annotated with match score
/// \c score and deletion bound ]c dels, into \c eps_file.
// ==============================================================
void write_logo_def(
	int index,
	const PWM& matrix,
	float score,
	int dels,
	std::ofstream& eps_file)
{
	eps_file<<"/M"<<index<<" {\n";
	eps_file<<"/numOfCols "<<(float)matrix.length<<" def\n";
	eps_file<<"/logoWidth "<<1.5*(float)matrix.length+1.<<" cm def\n";
	eps_file<<"("<<score<<" / "<<dels<<")\n";
	eps_file<<"DrawBox\n";
	eps_file<<"StartLogo\n";
	eps_file<<"StartLine % line number 1\n";
	int i=0;
	for (; i<matrix.length; ++i)
	{
		eps_file<<"("<<i<<") startstack\n";
		float_char_vector weights(4);
		float A_i,C_i,G_i,T_i,sum_i;
		A_i = matrix.A[i];
		C_i = matrix.C[i];
		G_i = matrix.G[i];
		T_i = matrix.T[i];
		sum_i = A_i+C_i+G_i+T_i;
		A_i /= sum_i; C_i /= sum_i;
		G_i /= sum_i; T_i /= sum_i;
		weights[0] = A_i>0. ? float_char(A_i,'A') : float_char(0.,'A');
		weights[1] = C_i>0. ? float_char(C_i,'C') : float_char(0.,'C');
		weights[2] = G_i>0. ? float_char(G_i,'G') : float_char(0.,'G');
		weights[3] = T_i>0. ? float_char(T_i,'T') : float_char(0.,'T');
		float ic  = A_i==0. ? 0. : A_i*log(A_i)/log(2.);
		ic += C_i==0. ? 0. : C_i*log(C_i)/log(2.);
		ic += G_i==0. ? 0. : G_i*log(G_i)/log(2.);
		ic += T_i==0. ? 0. : T_i*log(T_i)/log(2.);
		ic = 2.+ic;
		std::sort(weights.begin(),weights.end());
		if (weights[0].first>0.)
			eps_file<<weights[0].first*ic<<"  ("<<weights[0].second<<") numchar\n";
		if (weights[1].first>0.)
			eps_file<<weights[1].first*ic<<"  ("<<weights[1].second<<") numchar\n";
		if (weights[2].first>0.)
			eps_file<<weights[2].first*ic<<"  ("<<weights[2].second<<") numchar\n";
		if (weights[3].first>0.)
			eps_file<<weights[3].first*ic<<"  ("<<weights[3].second<<") numchar\n";
		eps_file<<"endstack\n";
	}
	eps_file<<"EndLine\nEndLogo\n} def\n\n";
}

/// produce an eps formmated representation of a gap in between two
/// conserved blocks, annotated with the length constraints given in
/// \c bounds, written into \c eps_file
// ==============================================================
void write_gap_logo(
	int_pair bounds,
	bool line_break,
	std::ofstream& eps_file,
	float max_line_length)
// ==============================================================
{
	rewind_length+=5.;
	if (rewind_length>max_line_length)
		eps_file<<rewind_length<<" cm\n";
	eps_file<<"("<<bounds.first<<".."<<bounds.second<<")\n";
	if (rewind_length>max_line_length) {
 		eps_file<<"DrawBreakGap\n";
		rewind_length=left_offset;
	} else {
 		eps_file<<"DrawGap\n";
	}
}

/// 
// ==============================================================
void write_eps_file(
	PWM_list& PWMs,	int_pair_list& BOUNDS,
	int_list& DELS, float_list& SCORES,
	float max_line_length,
	std::string eps_filename)
// ==============================================================
{
	std::ofstream eps_file(eps_filename.c_str(),std::ios::out);	
	int i;
	float m_score;
	int deletions;

	PWM_list::iterator m_it;
	PWM_list::iterator m_end;
	
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
	
	float bb_height = (9.+((float)num_of_lines-1.)*9.25)*72./2.54;
	float bb_width = (longest_line+2.5)*72./2.54;

	eps_file.precision(2);
	eps_file<<"%!PS-Adobe-3.0 EPSF-3.0\n";
	eps_file<<"%%Title: Sequence Logo : \n";
	eps_file<<"%%Creator: \n";
	eps_file<<"%%CreationDate: \n";
	eps_file<<"%%BoundingBox:   0  0  "<<bb_width<<" "<<bb_height<<"\n";
	eps_file<<"%%Pages: 0\n";
	eps_file<<"%%DocumentFonts: \n";
	eps_file<<"%%EndComments\n";

	eps_file<<
#include"template.tps"
		;

	
	m_it = PWMs.begin();
	m_end = PWMs.end();
	for (i=0; m_it!=m_end; ++i, ++m_it, ++B_it, ++D_it, ++S_it)
	{
		write_logo_def(i,*m_it,*S_it,*D_it,eps_file);
	}

	m_it = PWMs.begin();
	m_end = PWMs.end();
	B_it = BOUNDS.begin();
	B_end = BOUNDS.end();
	D_it = DELS.begin();
	D_end = DELS.end();
	S_it = SCORES.begin();
	S_end = SCORES.end();
	
	eps_file<<"\n\nGapLength "<<9.*(num_of_lines-1)+2.5<<" cm translate\n";
	left_offset=2.5;
	rewind_length=left_offset;
	for (i=0; m_it!=m_end; ++i, ++m_it, ++B_it, ++D_it, ++S_it)
	{
		if (i>0)
			write_gap_logo(*B_it,0.,eps_file,max_line_length);
		eps_file<<"M"<<i<<"\n";
		rewind_length += 1.5*(float)m_it->length+1.;
		if (i==0)
			eps_file<<"/yaxis false def\n\n";
	}

	eps_file<<"showpage\n\n";
	eps_file<<"%%EOF\n\n";

}


} // end namespace fragrep


