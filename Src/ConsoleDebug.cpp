// For debugging purposes only
#include"udf.h"
#include"filter.h"
using namespace mvkf;
int main()
{
    Sysmat *sm = new Sysmat;
    Filter<Dim>* f = new Filter<Dim>(sm);

    for (int i = 0; i < 100; ++i)
    {
        std::cout << "\n\nIteration " << i << std::endl;
        std::cout << "\nGain:\n";
        f->ComputeGain();
        Matrix<1, 1> y, u;
        double sign = (rand() % 2 == 1) ? -1 : 1;
        y(0, 0) = 1.0 + sign * (double)rand() / RAND_MAX;
        u(0, 0) = 1.0;
        std::cout << "\nState estimates:\n";
        f->UpdateEstimate(y, u);
        std::cout << "\nCovariance:\n";
        f->ComputeCovariance();
        std::cout << "\nProject ahead estimates:\n";
        f->ProjectAhead(u);
    }

	delete f; 
	delete sm;
	system("pause");
	return 0;
}


