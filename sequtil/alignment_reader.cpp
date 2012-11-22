#include "alignment_reader.hpp"
#include <iostream>
#include <fstream>
#define MAX_ALIGNMENTS 1024

// ==============================================================
alignment_reader::alignment_reader()
// ==============================================================
{
}

// ==============================================================
alignment_reader::alignment_reader(FILE* clustal_file)
// ==============================================================
{
	init(clustal_file);
}

// ==============================================================
void alignment_reader::init(FILE* clustal_file)
// ==============================================================
{
	char** aligned_seqs = new char*[MAX_ALIGNMENTS];
	char** names = new char*[MAX_ALIGNMENTS];
	int i;
	for (i=0; i<MAX_ALIGNMENTS; i++)
	{
		aligned_seqs[i] = NULL;
		names[i] = NULL;
	}
	int hght = read_alignment(clustal_file, aligned_seqs, names);
	fclose(clustal_file);
	for (i=0; i<hght; i++)
	{
		push_back(aligned_sequence(names[i],aligned_seqs[i]));
		//std::cout<<names[i]<<"\n"<<aligned_seqs[i]<<"\n";
  		delete[] names[i];
		delete[] aligned_seqs[i];
	}
	delete[] names;
	delete[] aligned_seqs;	
}

// ==============================================================
alignment_reader::alignment_reader(const std::string& filename)
// ==============================================================
{
	FILE* clustal_file = fopen(filename.c_str(),"r");
	init(clustal_file);
}

// ==============================================================
int alignment_reader::length()
// ==============================================================
{
	if (height()==0)
		return 0;
	return begin()->second.size();
}

// ==============================================================
int alignment_reader::height()
// ==============================================================
{
	return size();
}

// ==============================================================
char& alignment_reader::operator()(int i, int j)
// ==============================================================
{
	return (*this)[i].second[j];
}

// ==============================================================
void alignment_reader::write_fasta(const std::string& filename)
// ==============================================================
{
	std::ofstream ofs(filename.c_str(),std::ios::out);
	iterator s_it = begin();
	iterator s_end = end();
	for (; s_it!=s_end; ++s_it)
	{
		ofs<<">"<<s_it->first<<"\n";
		std::string seq = s_it->second;
		int len = seq.size();
		int i;
		for (i=0; i<len; i++)
			if (seq[i]!='-')
				ofs<<seq[i];
		ofs<<"\n";
	}
	ofs.close();
	
}
