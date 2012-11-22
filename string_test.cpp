#include <iostream>
#include <string>

using namespace std;

int main()
{
	string x("M0$CGTTGC");
	cout<<x<<" ->\n";
	cout<<"> "<<x.substr(0,x.find('$'))<<"\n";
	cout<<"> "<<x.substr(x.find('$')+1,x.size())<<"\n";
	string y("M0");
	cout<<y<<" ->\n";
	cout<<"> "<<y.substr(0,y.find('$'))<<"\n";
	cout<<"> "<<y.substr(y.find('$')+1,y.size())<<"\n";
	string z("$CGTTGC");
	cout<<z<<" ->\n";
	cout<<"> "<<z.substr(0,z.find('$'))<<"\n";
	cout<<"> "<<z.substr(z.find('$')+1,z.size())<<"\n";
}
