#include <iostream>
#include "Polynomial.hpp"

int main()
{
	Polynomial<double> b;
	b.AddCoef(2, 2);
	b.AddCoef(3, 0);
	Polynomial<double> o;
	o.AddCoef(1, 4);
	o.AddCoef(2, 3);
	o.AddCoef(2, 2);
	o.AddCoef(10, 1);
	o.AddCoef(1, 0);
	Polynomial<double> a;
	a.AddCoef(1, 4);
	a.AddCoef(1, 3);
	a.AddCoef(1, 2);
	a.AddCoef(1, 1);
	a.AddCoef(1, 0);
	o /= b;
//	o %= b;
	std::cout << o << std::endl;

	//for (Polynomial<double>::Iterator it = o.begin(); it != o.end(); ++it)
	//{
	//	std::cout << *it << ' ';
	//}

	//std::cout << std::endl;
//	Polynomial<double> a;
//	a = o * b;
//	b *= a;
//	b.print();
//	double y = o(2);
//	double x = o(1, 2);
//	std::cout << x << std::endl;
//	std::cout << y << std::endl;
//	--o;
//	std::cout << "        " << b << std::endl;
//	std::cout << o << std::endl;
//	b *= o;
//	try {
//	std::cout << o[6] << std::endl;
//	}
//	catch (std::out_of_range & a)
//	{
//		std::cout << a.what();
//	}

//	Polynomial<int> c;
//	std::cin >> c;
//	std::cout << c;
//	o /= b;
//	o.print();


	system("pause");
	return 0;
}