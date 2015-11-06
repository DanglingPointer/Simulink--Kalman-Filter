#include"matrix.h"
using namespace mvkf;
int main()
{
	IMatrix *p = new Matrix<3, 3>;
	p->at(0, 0) = 7;
	p->at(0, 1) = 2;
	p->at(0, 2) = 1;
	p->at(1, 0) = 0;
	p->at(1, 1) = 3;
	p->at(1, 2) = -1;
	p->at(2, 0) = -3;
	p->at(2, 1) = 4;
	p->at(2, 2) = -2;
	Print(p);
	std::cout << "Det = " << p->Det() << '\n';

	IMatrix *pnew = p->Inverse();
	Print(pnew);

	delete p;
	delete pnew;
	system("pause");
	return 0;
}