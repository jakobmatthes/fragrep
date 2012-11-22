#include "bbq_util.h"


// ==============================================================
bool is_int(const char* s)
// ==============================================================
{
	int i=0;
	while (s[i]!=0)
	{
		if (s[i]<'0' || s[i]>'9')
			return false;
		i++;
	}
	return true;
}

// ==============================================================
bool is_flt(const char* s)
// ==============================================================
{
	int i=0;
	int dot_cnt = 0;
	while (s[i]!=0)
	{
		if (s[i]=='.')
			dot_cnt++;
		if (dot_cnt>1)
		{
			return false;
		}
		if ((s[i]<'0' || s[i]>'9') && s[i]!='.')
		{
			return false;
		}
		i++;
	}
	return true;
}

/// read over all useless characters such as spaces, tabs and line
/// breaks
// ==============================================================
void skip(std::ifstream& ifs)
// ==============================================================
{ 
	if (!ifs.good())
		return;
	char ch = 0;
	ifs.get(ch);
	while (ifs.good() && ch<=32)
		ifs.get(ch);
	if (ifs.good())
		ifs.unget();
}

// ==============================================================
std::string bbq_exception::msg() const
// ==============================================================
{
	return MSG;
}

// ==============================================================
bbq_exception::bbq_exception()
// ==============================================================
{
	MSG=std::string("Unspecified error.\n");
}


// ==============================================================
bbq_exception::bbq_exception(const std::string& m)
// ==============================================================
{
	MSG=m;
}

