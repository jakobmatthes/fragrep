#include <sequtil/alignment_reader.hpp>
#include <iostream>
#include <cstdio>

alignment_reader AR;

int main(int argc, char** argv)
{
	if (argc>1)
		AR = alignment_reader(argv[1]);
	else
		AR = alignment_reader(stdin);

	alignment_reader::iterator a_it = AR.begin();
	alignment_reader::iterator a_end = AR.end();
	for (; a_it!=a_end; ++a_it)
	{
		int a_len = a_it->first.size();
		if (a_len>16) {
			std::cout<<a_it->first.substr(0,15)<<" ";
		} else {
			std::cout<<a_it->first<<std::string(16-a_len,' ');
		}
		std::cout<<a_it->second<<"\n";
	}
}

