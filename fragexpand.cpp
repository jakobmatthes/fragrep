/*
  fragrep c++ rebuild */

// std and iostream headers
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

// libsequence headers for reading fasta files
//#include <Sequence/Fasta.hpp>
//#include <Sequence/SeqRegexes.hpp>
//#include <Sequence/SeqProperties.hpp>

#include <sequtil/fasta.h>
/** \mainpage

welcome to fragrep.

*/

/// genome sequence
std::string T;

/// fragment sequence
std::string s;


int main(int argc, char** argv)
{

	if (argc<4)
	{
		std::cout<<"\nusage: fragexpand <pref_len> <suff_len> <frag> [<fasta_file>]\n\n";
		std::cout<<"input: nucleoteide sequence <frag> and (multi-) fasta formatted file <fasta_file>.\n";
		std::cout<<"  if <fasta_file> is not specified, input is taken from stdin.\n";
		std::cout<<"output: for each occurence of <frag> in <fasta_file>, <frag> plus a prefix and a suffix of specified length is printed.\n\n";
		return 0;
	}
	int pref_len;
	std::istringstream pref_istr(argv[1]);
	pref_istr>>pref_len;

	int suff_len;
	std::istringstream suff_istr(argv[2]);
	suff_istr>>suff_len;

	s = argv[3];
	
	// by default, fr reads from stdin
	fasta_reader fr;
	if (argc>4)
	{
		fr.open(argv[4]);
		std::string xxh,xxs;
		//fr.get_seq(xxh,xxs);
		//std::cout<<xxh<<"\n"<<xxs<<"\n";
	}

	std::string header,seq;

	while (fr.get_seq (header, seq)) {
		int up_to=0;
		while (true)
		{
			up_to=seq.find(s,up_to);
			if (up_to<0)
				break;
			int len=pref_len+suff_len+s.size();
			int pos=up_to-pref_len;
			if (pos<0) { len+=pos; pos=0;}
			std::cout<<seq.substr(pos,len)<<"\n";
			up_to++;
		}
		// now do the reverse complement
		std::string rev_seq = seq;
		fasta_reader::reverse_complement(rev_seq);
		up_to=0;
		while (true)
		{
			up_to=rev_seq.find(s,up_to);
			if (up_to<0)
				break;
			int len=pref_len+suff_len+s.size();
			int pos=up_to-pref_len;
			if (pos<0) { len+=pos; pos=0;}
			std::cout<<rev_seq.substr(pos,len)<<"\n";
			up_to++;
		}
		
	} 
	return 0;
	
	
}



