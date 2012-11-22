#ifndef BBQ_UTIL_H
#define BBQ_UTIL_H


#include <set>
#include <map>
#include <list>
#include <string>
#include <fstream>

typedef std::set<std::string> string_set;
typedef std::map<std::string, int> str_int_map;
typedef std::list<std::string> str_list;




bool is_int(const char* s);

bool is_flt(const char* s);

/// read over all useless characters such as spaces, tabs and line
/// breaks
void skip(std::ifstream& ifs);


class bbq_exception
{
	
	std::string MSG;
	
public:
	
	bbq_exception(const std::string&);
	
	bbq_exception();
	
	std::string msg() const; 
	
};

#endif

