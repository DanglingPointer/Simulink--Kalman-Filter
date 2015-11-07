#include"filter.h"
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

	IMatrix *pinv = p->Inverse();
	Print(pinv);

	IMatrix *pmult = pinv->Copy();
	Multiply(p, pinv, pmult);
	Print(pmult);

	IMatrix *pmultmult = pinv->Copy();
	Multiply(p, pinv, pmult, pmultmult);
	Print(pmultmult);

	IMatrix* peye = p->Eye();
	Print(peye);

	IMatrix* pone = new Matrix<1, 1>;
	IMatrix* poneeye = pone->Eye();
	Print(pone);
	Print(poneeye);


	delete p;
	delete pinv;
	delete pmult;
	delete pmultmult;
	delete peye;
	delete pone;
	delete poneeye;
	system("pause");
	return 0;
}