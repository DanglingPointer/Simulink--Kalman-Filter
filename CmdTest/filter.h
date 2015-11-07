#pragma once
#include"matrix.h"
/*
 * Generic Kalman-filter (no use of MatLab/Simulink libraries).
 * User has to learn IMatrix class in order to implement state matrix functions
 * as well as initial states vector and initial error covariance matrix.
 * For order of use of the loop functions see Brown&Whang p147.
 */
namespace mvkf
{
	class Filter
	{
	public:
		Filter(IMatrix*(*phi)(unsigned long long), IMatrix*(*h)(unsigned long long), 
			   IMatrix*(*q)(unsigned long long), IMatrix*(*r)(unsigned long long))
			:m_time(0), m_initialized(false), m_phi(phi), m_h(h), m_q(q), m_r(r), 
			m_pec(nullptr), m_peca(nullptr), m_psv(nullptr), m_psva(nullptr), m_pgain(nullptr)
		{ }
		~Filter()
		{
			if (m_pec != nullptr) delete m_pec;
			if (m_peca != nullptr) delete m_peca;
			if (m_psv != nullptr) delete m_psv;
			if (m_psva != nullptr) delete m_psva;
			if (m_pgain != nullptr) delete m_pgain;
		}
		// ---Filter loop functions---
		void ComputeGain()
		{	//p.146 eq(4.2.17)
			IMatrix *H = m_h(m_time);
			IMatrix *H_T = H->Transpose();
			IMatrix *R = m_r(m_time);

			if (m_pgain == nullptr)
				m_pgain = H_T->Copy();	// dim(K) = dim(H^T)

			IMatrix *innerproduct = R->Copy();
			Multiply(H, m_peca, H_T, innerproduct);

			IMatrix *parenth = Add(innerproduct, R);
			IMatrix *parinv = parenth->Inverse();

			Multiply(m_peca, H_T, parinv, m_pgain);

			delete innerproduct; delete R;
			delete parenth; delete parinv;
			delete H; delete H_T;
		}
		void UpdateEstimate(IMatrix *measurement)
		{	// p.144 eq(4.2.8) expanded
			if (m_psv != nullptr)
				delete m_psv;

			IMatrix *H = m_h(m_time);
			
			IMatrix *firstprod = m_psva->Copy();
			Multiply(m_pgain, measurement, firstprod);

			IMatrix *secondprod = m_psva->Copy();
			Multiply(m_pgain, H, m_psva, secondprod);

			IMatrix *sum = Add(m_psva, firstprod);
			m_psv = Substract(sum, secondprod);
			
			delete H; delete firstprod;
			delete secondprod; delete sum;
		}
		void ComputeCovariance()
		{	// p.146 eq(4.2.18)
			if (m_pec != nullptr)
				delete m_pec;
			
			IMatrix *H = m_h(m_time);
			IMatrix *R = m_r(m_time);

			IMatrix *eye = m_peca->Eye();
			IMatrix *firstprod = m_peca->Copy();
			Multiply(m_pgain, H, firstprod); //K*H

			IMatrix *parenth = Substract(eye, firstprod); // (I-K*H)
			IMatrix *parenthT = parenth->Transpose();

			IMatrix *secondprod = m_peca->Copy();
			Multiply(parenth, m_peca, parenthT, secondprod); // (I-KH)(P-)(I-KH)^T

			IMatrix *gainT = m_pgain->Transpose();

			IMatrix *thirdprod = m_peca->Copy();
			Multiply(m_pgain, R, gainT, thirdprod); // K*R*(K^T)

			m_pec = Add(secondprod, thirdprod);

			delete H; delete R; delete eye; delete firstprod;
			delete parenth; delete parenthT; delete secondprod;
			delete gainT; delete thirdprod;
		}
		void ProjectAhead(bool increment_time = true)
		{
			IMatrix *phi = m_phi(m_time);
			IMatrix *phiT = phi->Transpose();
			IMatrix *Q = m_q(m_time);

			if (m_peca != nullptr)
				delete m_peca;

			Multiply(phi, m_psv, m_psva);

			IMatrix *prod = m_pec->Copy();
			Multiply(phi, m_pec, phiT, prod);
			m_peca = Add(prod, Q);

			if (increment_time) ++m_time;

			delete phi; delete phiT; delete Q;
			delete prod;
		}
		// ---Mutators---
		void InitStates(IMatrix *pinitsv, IMatrix *pinitec)
		{
			if (m_peca != nullptr) delete m_peca;
			m_peca = pinitec;
			if (m_psva != nullptr) delete m_psva;
			m_psva = pinitsv;
		}
		void set_Time(unsigned long time)
		{
			m_time = time;
		}
		void Reset()
		{
			if (m_pec != nullptr) delete m_pec;
			m_pec = nullptr;
			if (m_psv != nullptr) delete m_psv;
			m_psv = nullptr;
			if (m_pgain != nullptr) delete m_pgain;
			m_pgain = nullptr;
		}
		// ---Accessors---
		double get_State(uint statenum) const // statenum from 1
		{
			return m_psv->at(statenum - 1, 0);
		}
		IMatrix* get_StateM() const
		{
			return m_psv->Copy();
		}
	private:
		unsigned long long m_time;
		bool m_initialized;
		// System matrices are user-defined functions of time:
		IMatrix*(*m_phi)(unsigned long long);
		IMatrix*(*m_h)(unsigned long long);
		IMatrix*(*m_q)(unsigned long long);
		IMatrix*(*m_r)(unsigned long long);
		// Data:
		IMatrix *m_pec;		// error covariance
		IMatrix *m_peca;	// error covariance projection ahead
		IMatrix *m_psv;		// states vector
		IMatrix *m_psva;	// states vector projection ahead
		IMatrix *m_pgain;	// Kalman gain
	};
} // namespace