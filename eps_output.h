#include <sequtil/position_weight_matrix.h>
#include <list>
#include <vector>
#include <string>

using namespace bbq;

namespace fragrep
{
	
typedef std::pair<int,int> int_pair;
typedef std::list<int_pair> int_pair_list;
typedef std::list<int> int_list;
typedef std::list<float> float_list;
typedef std::pair<float,char> float_char;
typedef std::vector<float_char> float_char_vector;

typedef bbq::position_weight_matrix PWM;
typedef std::list<PWM> PWM_list;

void write_eps_file(
	PWM_list& PWMs,	int_pair_list& BOUNDS,
	int_list& DELS, float_list& SCORES,
	float max_line_length,
	std::string eps_file_name);

} // end namespace fragrep

