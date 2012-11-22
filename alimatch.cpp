#include <sequtil/alignment_reader.hpp>
#include <iostream>
#include <cstdio>

alignment_reader AR;

float column_score(int j)
{
	int hght = AR.height()-1;
	int i;
	int freq[256];
	for (i=0; i<255; i++)
		freq[i]=0;
	for (i=0; i<hght; i++)
	{
		char ch_i_j = AR(i,j);
		freq[ch_i_j]++;
	}
	return (float)freq[AR(AR.height()-1,j)]/(float)hght;
}

int main(int argc, char** argv)
{
	if (argc>1)
		AR = alignment_reader(argv[1]);
	else
		AR = alignment_reader(stdin);
	//alignment_reader::iterator a_it = AR.begin();
	//alignment_reader::iterator a_end = AR.end();

	//for (; a_it!=a_end; ++a_it)
	//{
	//	std::cout<<a_it->first<<": "<<a_it->second.length()<<" chrs\n";
	//}

	int j;
	int len = AR.length();
	float total_score = 0;
	for (j=0; j<len; j++)
	{
		float col_sc_j = column_score(j);
		total_score += col_sc_j;
	}
	std::cout<<"\nabsolute score: "<<total_score<<"\n";
	std::cout<<"normalized score: "<<total_score/(float)AR.length()<<"\n";
	
}

