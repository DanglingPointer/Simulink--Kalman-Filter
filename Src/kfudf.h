#pragma once
/*
 * User-defined functionalities
 */
#include"kfmatrix.h"
using namespace mvkf;
// Constant system matrices
class SysMat :public ISysMat
{
public:
	SysMat() :m_time(0.0), m_currentstates(nullptr)
	{
		m_a = new Matrix<5, 5>(true);
		m_a->at(0, 0) = 0.9969;
		m_a->at(0, 1) = 0.09982;
		m_a->at(1, 0) = -0.0611;
		m_a->at(1, 1) = 0.9955;
		m_a->at(2, 3) = 0.09994;
		m_a->at(2, 4) = -9.965e-06;
		m_a->at(3, 3) = 0.9989;
		m_a->at(3, 4) = -0.0001993;

		m_b = new Matrix<5, 1>;
		m_b->at(2, 0) = 9.965e-06;
		m_b->at(3, 0) = 0.0001993;

		m_g = new Matrix<5, 2>;
		m_g->at(0, 0) = 2.847e-06;
		m_g->at(1, 0) = 5.691e-05;
		m_g->at(2, 1) = -3.322e-07;
		m_g->at(3, 1) = -9.965e-06;
		m_g->at(4, 1) = 0.1;

		m_c = new Matrix<1, 5>;
		m_c->at(0, 1) = 1;
		m_c->at(0, 2) = 1;

		m_q = new Matrix<2, 2>;
		m_q->at(0, 0) = 30;
		m_q->at(1, 1) = 1e-06;

		m_r = new Matrix<1, 1>;
		m_r->at(0, 0) = 1.0002e-08;			// = (var(v) / simulation_time)

		m_isv = new Matrix<5, 1>;
		
		m_icov = new Matrix<5, 5>(true);
		m_icov->at(1, 1) = 0.013;
		m_icov->at(2, 2) = 3.14159265359*3.14159265359;
		m_icov->at(4, 4) = 2.5e-04;
	}
	~SysMat()
	{
		delete m_a; delete m_b; delete m_g; delete m_c; delete m_q; 
		delete m_r; delete m_isv; delete m_icov;
		if (m_currentstates != nullptr) delete m_currentstates;
	}
	const IMatrix* get_A() const
	{
		return m_a;
	}
	const IMatrix* get_B() const
	{
		return m_b;
	}
	const IMatrix* get_G() const
	{
		return m_g;
	}
	const IMatrix* get_C() const
	{
		return m_c;
	}
	const IMatrix* get_Q() const
	{
		return m_q;
	}
	const IMatrix* get_R() const
	{
		return m_r;
	}
	const IMatrix* get_InitSV() const
	{
		return m_isv;
	}
	const IMatrix* get_InitCov() const
	{
		return m_icov;
	}
	void set_States(const IMatrix *states)
	{
		if (m_currentstates != nullptr) delete m_currentstates;
		m_currentstates = states->Copy();		// unnecessary?
	}
	void set_Time(double time)
	{
		m_time = time;
	}
private:
	IMatrix *m_a, *m_b, *m_g, *m_c, *m_q, *m_r, *m_isv, *m_icov, *m_currentstates;
	double m_time;
};