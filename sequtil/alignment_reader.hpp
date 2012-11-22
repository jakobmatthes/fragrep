#ifndef ALIGNMENT_READER__
#define ALIGNMENT_READER__

#include <cstdio>
extern "C" {
#include "aln_util.h"
}

#include <string>
#include <vector>
#include <cstdio>

typedef std::pair<std::string, std::string> aligned_sequence;

typedef std::vector<aligned_sequence> sequence_vector;


class alignment_reader : public sequence_vector
{
	
protected:

	void init(FILE*);
	
public:

	alignment_reader();

	alignment_reader(FILE*);
	
	/// read clustalw format from the file specified in filename
	alignment_reader(const std::string& filename);

	/// return char at pos j of seq i
	char& operator()(int i, int j);
	
	int length();

	int height();

	void write_fasta(const std::string& filename);
	
};


#endif
