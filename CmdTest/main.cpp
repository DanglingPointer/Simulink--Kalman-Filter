#include"matrix.h"

int main()
{
	IMatrix *p = new Matrix<3, 3>;
	p->at(0, 0) = 1;
	p->at(0, 1) = 2;
	p->at(1, 0) = 3;
	p->at(1, 1) = 4;
	p->at(2, 2) = 1;
	//Print(p);
	std::cout << "Det = " << p->Det() << '\n';
	delete p;


	system("pause");
	return 0;
}