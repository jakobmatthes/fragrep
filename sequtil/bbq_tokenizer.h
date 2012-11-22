#ifndef BBQ_TOKENIZER_H
#define BBQ_TOKENIZER_H


#include <string>


namespace bbq {

class tokenizer;
	
class tokenizer_iterator
{

public:

	tokenizer* parent;
	
	int token_pos;

	int token_end;

	std::string curr_token;
	
	tokenizer_iterator(int pos, tokenizer* parent_);

	tokenizer_iterator& operator++();

	tokenizer_iterator& operator++(int);

	bool is_separator(char x);
	
	bool operator!=(const tokenizer_iterator& t);

	bool operator==(const tokenizer_iterator& t);

	std::string operator*();
};

	
class tokenizer
{
	
public:

	std::string str;

	std::string separators;
	
	typedef tokenizer_iterator iterator;
	
	tokenizer(const std::string& str_);

	void set_separators(const std::string& seps);
	
	iterator begin();

	iterator end();
	
};


}

#endif

