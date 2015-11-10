// For debugging purposes only
#include"udf.h"
#include"filter.h"
using namespace mvkf;
int main()
{
    Sysmat *sm = new Sysmat;
    Filter<Dim>* f = new Filter<Dim>(sm);

    f->ComputeGain();
    Matrix<1, 1> y, u;
    y(0, 0) = 1;
    u(0, 0) = 1;
    f->UpdateEstimate(y, u);
    f->ComputeCovariance();
    f->ProjectAhead(u);

	delete f; 
	delete sm;
	system("pause");
	return 0;
}


