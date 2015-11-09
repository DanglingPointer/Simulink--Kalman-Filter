// For debugging purposes only
//#include"kffilter.h"
#include"kfudf.h"
#include"kffilter.h"
using namespace mvkf;
int main()
{
	//Matrix *p = new Matrix(3, 3);
	//p->at(0, 0) = 7;
	//p->at(0, 1) = 2;
	//p->at(0, 2) = 1;
	//p->at(1, 0) = 0;
	//p->at(1, 1) = 3;
	//p->at(1, 2) = -1;
	//p->at(2, 0) = -3;
	//p->at(2, 1) = 4;
	//p->at(2, 2) = -2;
	//Print(*p);
	//Matrix *pcpy = new Matrix(*p);
	//Print(*pcpy);
	//std::cout << "Det = " << p->Det() << '\n';

	//Matrix *pinv = new Matrix(p->Inverse());
	//Print(*pinv);

	//Print((*p)*(*pinv));
	//Print((*p)*(*pinv)*(p->Eye()));

	//delete p;
	//delete pcpy;
	//delete pinv;


	ISysMat *sysmat = new SysMat;
	Filter *pf = new Filter(sysmat);
	
	// Step 1:
	pf->ComputeGain();
	
	// Retrieving inputs:
	Matrix y(1, 1);
	y.at(0, 0) = 0;						// measurement is first input
	Matrix u(1,1);
	u.at(0, 0) = 1;						// reference signal is second input
	
	// Step 2:
	pf->UpdateEstimate(y);
	
	// Step 3:
	pf->ComputeCovariance();
	
	// Step 4:
	pf->ProjectAhead(u);

	delete pf; 
	delete sysmat;
	system("pause");
	return 0;
}


