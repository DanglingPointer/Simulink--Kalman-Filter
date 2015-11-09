#pragma once
/*
 * User-defined functionalities
 */
#include"kffilter.h"
using namespace mvkf;
// Constant system matrices
class SysMat :public ISysMat
{
public:
	SysMat() :m_a(5,5), m_b(5,1), m_g(5,2), m_c(1,5), m_q(2,2), m_r(1,1), m_isv(5,1), m_icov(5,5), m_pcurrentstates(nullptr), m_time(0.0)
	{
		m_a.at(0, 0) = 0.9969;
		m_a.at(0, 1) = 0.09982;
		m_a.at(1, 0) = -0.0611;
		m_a.at(1, 1) = 0.9955;
		m_a.at(2, 2) = 1.0;
		m_a.at(2, 3) = 0.09994;
		m_a.at(2, 4) = -9.965e-06;
		m_a.at(3, 3) = 0.9989;
		m_a.at(3, 4) = -0.0001993;
		m_a.at(4, 4) = 1.0;

		m_b.at(2, 0) = 9.965e-06;
		m_b.at(3, 0) = 0.0001993;

		m_g.at(0, 0) = 2.847e-06;
		m_g.at(1, 0) = 5.691e-05;
		m_g.at(2, 1) = -3.322e-07;
		m_g.at(3, 1) = -9.965e-06;
		m_g.at(4, 1) = 0.1;

		m_c.at(0, 1) = 1;
		m_c.at(0, 2) = 1;

		m_q.at(0, 0) = 30;
		m_q.at(1, 1) = 1e-06;

		m_r.at(0, 0) = 1.0002e-08;			// = (var(v) / simulation_time)

		m_icov.at(0, 0) = 1.0;
		m_icov.at(1, 1) = 0.013;
		m_icov.at(2, 2) = 3.14159265359*3.14159265359;
		m_icov.at(3, 3) = 1.0;
		m_icov.at(4, 4) = 2.5e-04;
	}
	const Matrix& get_A() const
	{
		return m_a;
	}
	const Matrix& get_B() const
	{
		return m_b;
	}
	const Matrix& get_G() const
	{
		return m_g;
	}
	const Matrix& get_C() const
	{
		return m_c;
	}
	const Matrix& get_Q() const
	{
		return m_q;
	}
	const Matrix& get_R() const
	{
		return m_r;
	}
	const Matrix& get_InitSV() const
	{
		return m_isv;
	}
	const Matrix& get_InitCov() const
	{
		return m_icov;
	}
	void set_States(const Matrix *states)
	{
		m_pcurrentstates = states;
	}
	void set_Time(double time)
	{
		m_time = time;
	}
private:
	Matrix m_a, m_b, m_g, m_c, m_q, m_r, m_isv, m_icov;
	const Matrix *m_pcurrentstates;
	double m_time;
};