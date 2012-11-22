#include "bbq_tokenizer.h"


namespace bbq
{

	bool tokenizer_iterator::is_separator(char x)
	{
		return (parent->separators.find(x)!=std::string::npos);
	}
	
	bool tokenizer_iterator::operator!=(const tokenizer_iterator& t)
	{
		return !(*this==t);
	}
	
	bool tokenizer_iterator::operator==(const tokenizer_iterator& t)
	{
		if (t.parent==this->parent && t.token_pos==this->token_pos)
			return true;
		return false;
	}
	
	tokenizer_iterator::tokenizer_iterator(int pos_, tokenizer* parent_)
	{
		curr_token = std::string("");
		token_pos = pos_;
		parent = parent_;
		token_end = pos_;
		
		while (token_pos<parent->str.size() && is_separator(parent->str[token_pos]))
			token_pos++;
		token_end = token_pos;
		int length;
		while (token_end<parent->str.size() && !is_separator(parent->str[token_end]))
		{
			curr_token.push_back(parent->str[token_end]);
			token_end++;
		}
	}

	std::string tokenizer_iterator::operator*()
	{
		return curr_token;
	}
	
	tokenizer_iterator& tokenizer_iterator::operator++()
	{
		curr_token = std::string("");

		int sz = parent->str.size();
		bool abort;
		abort  = !(token_end<sz);
		if (!abort)
			abort = !is_separator(parent->str[token_end]);
		while (!abort)
		{
			token_end++;
			abort  = !(token_end<parent->str.size());
			if (!abort)
				abort = !is_separator(parent->str[token_end]);
		}
		token_pos = token_end;
		
		abort  = !(token_end<parent->str.size());
		if (!abort)
			abort = is_separator(parent->str[token_end]);
		while (!abort)
		{
			curr_token.push_back(parent->str[token_end]);
			token_end++;
			abort  = !(token_end<parent->str.size());
			if (!abort)
				abort = is_separator(parent->str[token_end]);
		}
	}

	tokenizer_iterator& tokenizer_iterator::operator++(int)
	{
		curr_token = std::string("");

		int sz = parent->str.size();
		bool abort;
		abort  = !(token_end<sz);
		if (!abort)
			abort = !is_separator(parent->str[token_end]);
		while (!abort)
		{
			token_end++;
			abort  = !(token_end<parent->str.size());
			if (!abort)
				abort = !is_separator(parent->str[token_end]);
		}
		token_pos = token_end;
		
		abort  = !(token_end<parent->str.size());
		if (!abort)
			abort = is_separator(parent->str[token_end]);
		while (!abort)
		{
			curr_token.push_back(parent->str[token_end]);
			token_end++;
			abort  = !(token_end<parent->str.size());
			if (!abort)
				abort = is_separator(parent->str[token_end]);
		}
	}
	
	tokenizer::tokenizer(const std::string& str_)
	{
		str = str_;
		separators = " \t";
	}

	void tokenizer::set_separators(const std::string& seps)
	{
		separators = seps;
	}
	
	tokenizer::iterator tokenizer::begin()
	{
		return iterator(0,this);
	}
	
	tokenizer::iterator tokenizer::end()
	{
		return iterator(str.size(),this);
	}


}

