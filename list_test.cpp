#include <iostream>
#include <list>
#include <sstream>

using namespace std;

list<int> x;

void print()
{
	list<int>::iterator it=x.begin();
	for (;it!=x.end();++it)
		cout<<*it<<"; ";
	cout<<"\n";
}

void r_print()
{
	list<int>::reverse_iterator it=x.rbegin();
	for (;it!=x.rend();++it)
		cout<<*it<<"; ";
	cout<<"\n";
}

int main()
{
	//x.insert(x.end(),5);
	//x.insert(x.end(),7);
	//print();
	x.push_back(9);
	x.push_back(21);
	x.push_back(7);
	list<int>::reverse_iterator it=x.rbegin();
	it++;
	print();
	r_print();
	std::cout<<"i : "<<*(x.insert(x.end(),12))<<"\n";
	x.erase(x.begin(),x.begin());
	print();
	istringstream iss("5-3 matrices -6 - 7");
	while (!iss.fail())
	{
		string s;
		iss>>s;
		std::cerr<<"> "<<s<<"\n";
	}
	return 0;
	
}
