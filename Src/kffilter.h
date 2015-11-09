#pragma once
/*
 * Generic Kalman-filter (no use of Simulink API).
 * User has to learn IMatrix class in order to implement state matrix functions
 * as well as initial states vector and initial error covariance matrix.
 * For order of use of the loop functions see Brown&Whang p147.
 *
 * System must be on the form:
 * x[k+1] = A*x[k] + B*u[k] + G*w
 * y[k] = C*x[k] + v
 */
#include"kfsimplemat.h"
namespace mvkf
{
	// User must derive from this class
	__interface ISysMat
	{
		const Matrix& get_A() const;
		const Matrix& get_B() const;
		const Matrix& get_G() const;
		const Matrix& get_C() const;
		const Matrix& get_Q() const;
		const Matrix& get_R() const;
		const Matrix& get_InitSV() const;
		const Matrix& get_InitCov() const;
		void set_States(const Matrix *states);
		void set_Time(double time);
	};
	class Filter
	{
	public:
		Filter(ISysMat *sm) :m_psysmat(sm),
			m_ec(m_psysmat->get_InitCov()), m_eca(m_psysmat->get_InitCov()), m_sv(m_psysmat->get_InitSV()), 
			m_sva(m_psysmat->get_InitSV()), m_gain(m_psysmat->get_C().Transpose())
		{
			m_psysmat->set_States(&m_sv);
		}
		Filter(const Filter&) = delete;
		/* 
		 * Filter loop functions
		 */
		// p.146 eq(4.2.17)
		void ComputeGain()
		{
			Matrix C = m_psysmat->get_C();
			Matrix parenth = ((C * m_eca) * C.Transpose()) + m_psysmat->get_R();
			m_gain = (m_eca * C.Transpose()) * parenth.Inverse();
			Print(m_gain);
		}
		// p.144 eq(4.2.8) expanded
		void UpdateEstimate(const Matrix& y)
		{
			Matrix C = m_psysmat->get_C();
			m_sv = m_sva + (m_gain * (y - (C * m_sva)));
			Print(m_sv);
		}
		// p.146 eq(4.2.18)
		void ComputeCovariance()
		{
			Matrix C = m_psysmat->get_C();
			Matrix R = m_psysmat->get_R();
			Matrix I = m_eca.Eye();
			Matrix parenth = I - (m_gain*C);
			Matrix prod = (m_gain * R) * m_gain.Transpose();
			m_ec = ((parenth * m_eca) * parenth.Transpose()) + prod;
			Print(m_ec);
		}
		// assignment 5 prob.2 f) and g)
		void ProjectAhead(const Matrix& u)
		{
			Matrix A = m_psysmat->get_A();
			Matrix Q = m_psysmat->get_Q();
			Matrix B = m_psysmat->get_B();
			Matrix G = m_psysmat->get_G();
			
			m_sva = (A * m_sv) + (B * u);
			m_eca = ((A * m_ec) * A.Transpose()) + ((G * Q) * G.Transpose());
			Print(m_sva); Print(m_eca);
		}
		/*
		 * Mutator functions
		 */
		// Initial states from arguments
		void InitStates(const Matrix &initsv, const Matrix &initcov)
		{
			m_eca = initcov;
			m_sva = initsv;
		}
		void set_Time(double time)
		{
			m_psysmat->set_Time(time);
		}
		/*
		 * Accessor functions
		 */
		// statenum starts at 0
		double get_State(uint statenum) const
		{
			return m_sv.at(statenum, 0);
		}
		// Doesn't create a new object
		const Matrix& get_Statevec() const
		{
			return m_sv;
		}
	private:
		// Object for retrieving system matrices:
		ISysMat *m_psysmat;
		// Data:
		Matrix m_ec;		// error covariance
		Matrix m_eca;	// error covariance projection ahead
		Matrix m_sv;		// states vector
		Matrix m_sva;	// states vector projection ahead
		Matrix m_gain;	// Kalman gain
	};
} // namespace