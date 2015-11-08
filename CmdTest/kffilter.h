#pragma once
/*
 * Generic Kalman-filter (no use of Simulink API).
 * User has to learn IMatrix class in order to implement state matrix functions
 * as well as initial states vector and initial error covariance matrix.
 * For order of use of the loop functions see Brown&Whang p147.
 */
#include"kfmatrix.h"
namespace mvkf
{
	class Filter
	{
	public:
		Filter(IMatrix*(*phi)(double), IMatrix*(*h)(double), 
			   IMatrix*(*q)(double), IMatrix*(*r)(double))
			:m_time(0.0), m_initialized(false), m_phi(phi), m_h(h), m_q(q), m_r(r),
			m_pec(nullptr), m_peca(nullptr), m_psv(nullptr), m_psva(nullptr), m_pgain(nullptr)
		{ }
		Filter(const Filter&) = delete;
		~Filter()
		{
			if (m_pec != nullptr) delete m_pec;
			if (m_peca != nullptr) delete m_peca;
			if (m_psv != nullptr) delete m_psv;
			if (m_psva != nullptr) delete m_psva;
			if (m_pgain != nullptr) delete m_pgain;
		}
		/* 
		 * Filter loop functions
		 */
		// p.146 eq(4.2.17)
		void ComputeGain()
		{	
			if (!m_initialized)
				throw std::runtime_error("Filter::ComputeGain()");

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
		// p.144 eq(4.2.8) expanded
		void UpdateEstimate(const IMatrix* measurement)
		{	
			if (!m_initialized)
				throw std::runtime_error("Filter::UpdateEstimate()");

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
		// p.146 eq(4.2.18)
		void ComputeCovariance()
		{
			if (!m_initialized)
				throw std::runtime_error("Filter::ComputeCovariance()");

			if (m_pec != nullptr)
				delete m_pec;
			
			IMatrix *H = m_h(m_time);
			IMatrix *R = m_r(m_time);

			IMatrix *eye = m_peca->Eye();
			IMatrix *firstprod = m_peca->Copy();
			Multiply(m_pgain, H, firstprod); // K*H

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
		// p.147 eq(4.2.23) & (4.2.25)
		void ProjectAhead()
		{
			if (!m_initialized)
				throw std::runtime_error("Filter::ProjectAhead()");

			IMatrix *phi = m_phi(m_time);
			IMatrix *phiT = phi->Transpose();
			IMatrix *Q = m_q(m_time);

			if (m_peca != nullptr)
				delete m_peca;

			Multiply(phi, m_psv, m_psva);

			IMatrix *prod = m_pec->Copy();
			Multiply(phi, m_pec, phiT, prod);
			m_peca = Add(prod, Q);

			delete phi; delete phiT; delete Q;
			delete prod;
		}
		/*
		 * Mutator functions
		 */
		void InitStates(const IMatrix *pinitsv, const IMatrix *pinitec)
		{
			if (m_peca != nullptr) delete m_peca;
			m_peca = pinitec->Copy();
			if (m_psva != nullptr) delete m_psva;
			m_psva = pinitsv->Copy();
			m_initialized = true;
		}
		void set_Time(double time)
		{
			m_time = time;
		}
		void Reset(bool full = false)
		{
			if (m_pec != nullptr) delete m_pec;
			m_pec = nullptr;
			if (m_psv != nullptr) delete m_psv;
			m_psv = nullptr;
			if (m_pgain != nullptr) delete m_pgain;
			m_pgain = nullptr;

			if (full)
			{
				if (m_peca != nullptr) delete m_peca;
				m_peca = nullptr;
				if (m_psva != nullptr) delete m_psva;
				m_psva = nullptr;
				m_initialized = false;
			}
		}
		/*
		 * Accessor functions
		 */
		// statenum starts at 0
		double get_State(uint statenum) const
		{
			return m_psv->at(statenum, 0);
		}
		// Creates a new dynamic object
		IMatrix* get_StateM() const
		{
			return m_psv->Copy();
		}
		bool Initialized() const
		{
			return m_initialized;
		}
	private:
		double m_time;
		bool m_initialized;
		// System matrices are user-defined functions of time:
		IMatrix*(*m_phi)(double);
		IMatrix*(*m_h)(double);
		IMatrix*(*m_q)(double);
		IMatrix*(*m_r)(double);
		// Data:
		IMatrix *m_pec;		// error covariance
		IMatrix *m_peca;	// error covariance projection ahead
		IMatrix *m_psv;		// states vector
		IMatrix *m_psva;	// states vector projection ahead
		IMatrix *m_pgain;	// Kalman gain
	};
} // namespace