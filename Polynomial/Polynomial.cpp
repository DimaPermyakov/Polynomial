#include "Header.h"

int main() {
	std::vector<double> vec;
	vec.push_back(5);
	vec.push_back(4);
	vec.push_back(3);
	vec.push_back(2);
	Polynomial<double> one = vec;
	one.Print();
	one.Degree(3);

	std::cout << std::endl;
	size_t size = one.Size();
	one.Print();

	return 0;
}
