#include <sequtil/position_weight_matrix.h>
#include <limits>

using namespace bbq;

typedef std::vector<float> float_vector;

int main(int argc, char** argv)
{
	float_vector A,C,G,T;
	A.push_back(0); A.push_back(0); A.push_back(0); A.push_back(0); A.push_back(1);
	C.push_back(0); C.push_back(0); C.push_back(0); C.push_back(0); C.push_back(0);
	G.push_back(0); G.push_back(1); G.push_back(0); G.push_back(1); G.push_back(0);
	T.push_back(1); T.push_back(0); T.push_back(1); T.push_back(0); T.push_back(0);
	position_weight_matrix M(A,C,G,T);
	M.set_indel_costs(1.,0.,2);
	float_pair fp = M.get_frac_score(argv[1]);
	std::cout<<"fwd score="<<fp.first<<"; bwd score="<<fp.second<<"\n";
	float_pair fp2 = M.get_match_score(argv[1]);
	std::cout<<fp2.first<<"\n";
}
