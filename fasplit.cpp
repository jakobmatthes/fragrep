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

	int file_size=-1;
	std::string prefix = std::string("fasplit_");
	int i=1;
	while (i<argc && argv[i][0]=='-')
	{
		if (argv[i][1]=='s')
		{
			++i;
			std::istringstream istr(argv[i]);
			if ((istr>>file_size))
			{
				if (!istr.eof())
				{
					std::cerr<<"read DIM\n";
					char DIM;
					if (istr>>DIM)
					{
						std::cerr<<"  > "<<DIM<<"\n";
						if (DIM=='k' || DIM=='K')
							file_size *= 1024;
						else if (DIM=='m' || DIM=='M')
							file_size *= (1024*1024);
						else {
							std::cerr<<"ignoring unknown magnitude '"<<DIM;
							std::cerr<<"';\n valid magnitudes are 'K' and 'M'\n";
						}
					}
				}
			} else {
				std::cerr<<"expected integer file size.\n";
			}
		} else if (argv[i][1]=='p')
		{
			++i;
			prefix=std::string(argv[i]);
		} 
		++i;
	}
	
	if (i>argc)
		std::cerr<<"usage: fasplit [-s <file_size>] [-p <prefix>]\n\nfasta input is read from stdin.\n";
	// by default, fr reads from stdin
	fasta_reader fr;
	std::cerr<<"PRE: "<<prefix<<"\n";
	std::cerr<<"NUM: "<<file_size<<"\n";

	std::string header,seq;
	i=0;
	bool done=false;
	while (!done)
	{
		std::ostringstream of_str;
		of_str<<prefix<<i<<".fa";
		int curr_size=0;
		std::string of_name = of_str.str();
		std::ofstream ofs(of_name.c_str());
		do {
			if (fr.get_seq (header, seq))
			{
				ofs<<header<<"\n"<<seq<<"\n";
				curr_size += seq.size();
			} else {
				done=true;
			}
		} while (curr_size<=file_size && !done);
		ofs.close();
		++i;
	} 
	return 0;
	
	
}



