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
#include"kfmatrix.h"
namespace mvkf
{
	// User must derive from this class
	__interface ISysMat
	{
		const IMatrix* get_A() const;
		const IMatrix* get_B() const;
		const IMatrix* get_G() const;
		const IMatrix* get_C() const;
		const IMatrix* get_Q() const;
		const IMatrix* get_R() const;
		const IMatrix* get_InitSV() const;
		const IMatrix* get_InitCov() const;
		void set_States(const IMatrix *states);
		void set_Time(double time);
	};
	class Filter
	{
	public:
		Filter(ISysMat *sm) :m_initialized(false), m_psysmat(sm),
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

			const IMatrix *C = m_psysmat->get_C();
			IMatrix *C_T = C->Transpose();
			const IMatrix *R = m_psysmat->get_R();

			if (m_pgain == nullptr)
				m_pgain = C_T->Copy();	// dim(K) = dim(C^T)

			IMatrix *innerproduct = R->Copy();
			Multiply(C, m_peca, C_T, innerproduct);

			IMatrix *parenth = Add(innerproduct, R);
			IMatrix *parinv = parenth->Inverse();

			Multiply(m_peca, C_T, parinv, m_pgain);

			delete innerproduct;
			delete parenth; delete parinv;
			delete C_T;
		}
		// p.144 eq(4.2.8) expanded
		void UpdateEstimate(const IMatrix *y)
		{	
			if (!m_initialized)
				throw std::runtime_error("Filter::UpdateEstimate()");

			if (m_psv != nullptr)
				delete m_psv;

			const IMatrix *C = m_psysmat->get_C();
			
			IMatrix *firstprod = m_psva->Copy();
			Multiply(m_pgain, y, firstprod);

			IMatrix *secondprod = m_psva->Copy();
			Multiply(m_pgain, C, m_psva, secondprod);

			IMatrix *sum = Add(m_psva, firstprod);
			m_psv = Substract(sum, secondprod);
			
			delete firstprod; delete secondprod; delete sum;
		}
		// p.146 eq(4.2.18)
		void ComputeCovariance()
		{
			if (!m_initialized)
				throw std::runtime_error("Filter::ComputeCovariance()");

			if (m_pec != nullptr)
				delete m_pec;
			
			const IMatrix *C = m_psysmat->get_C();
			const IMatrix *R =m_psysmat->get_R();

			IMatrix *eye = m_peca->Eye();
			IMatrix *firstprod = m_peca->Copy();
			Multiply(m_pgain, C, firstprod); // K*C

			IMatrix *parenth = Substract(eye, firstprod); // (I-K*C)
			IMatrix *parenthT = parenth->Transpose();

			IMatrix *secondprod = m_peca->Copy();
			Multiply(parenth, m_peca, parenthT, secondprod); // (I-KC)(P-)(I-KC)^T

			IMatrix *gainT = m_pgain->Transpose();

			IMatrix *thirdprod = m_peca->Copy();
			Multiply(m_pgain, R, gainT, thirdprod); // K*R*(K^T)

			m_pec = Add(secondprod, thirdprod);

			delete eye; delete firstprod;
			delete parenth; delete parenthT; delete secondprod;
			delete gainT; delete thirdprod;
		}
		// assignment 5 prob.2 f) and g)
		void ProjectAhead(const IMatrix *u)
		{
			if (!m_initialized)
				throw std::runtime_error("Filter::ProjectAhead()");

			const IMatrix *A = m_psysmat->get_A();
			IMatrix *AT = A->Transpose();
			const IMatrix *Q = m_psysmat->get_Q();
			const IMatrix *B = m_psysmat->get_B();
			const IMatrix *G = m_psysmat->get_G();
			IMatrix *GT = G->Transpose();

			if (m_peca != nullptr) delete m_peca;
			if (m_psva != nullptr) delete m_psva;

			IMatrix *firstprod = m_psv->Copy();
			Multiply(A, m_psv, firstprod);

			IMatrix *secondprod = m_psv->Copy();
			Multiply(B, u, secondprod);
			
			m_psva = Add(firstprod, secondprod);

			delete firstprod; delete secondprod;

			firstprod = m_pec->Copy();
			Multiply(A, m_pec, AT, firstprod);

			secondprod = m_pec->Copy();
			Multiply(G, Q, GT, secondprod);
			m_peca = Add(firstprod, secondprod);

			delete AT; delete GT;
			delete firstprod; delete secondprod;
		}
		/*
		 * Mutator functions
		 */
		// Initial states from ISysMat object
		void InitStates()
		{
			if (m_peca != nullptr) delete m_peca;
			m_peca = m_psysmat->get_InitCov()->Copy();
			if (m_psva != nullptr) delete m_psva;
			m_psva = m_psysmat->get_InitSV()->Copy();
			m_initialized = true;
		}
		// Initial states from parameters
		void InitStates(const IMatrix *initsv, const IMatrix *initcov)
		{
			if (m_peca != nullptr) delete m_peca;
			m_peca = initcov->Copy();
			if (m_psva != nullptr) delete m_psva;
			m_psva = initsv->Copy();
			m_initialized = true;
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
			return m_psv->at(statenum, 0);
		}
		// Doesn't create a new object
		const IMatrix* get_Statevec() const
		{
			return m_psv;
		}
		bool Initialized() const
		{
			return m_initialized;
		}
	private:
		bool m_initialized;
		// Object for retrieving system matrices:
		ISysMat *m_psysmat;
		// Data:
		IMatrix *m_pec;		// error covariance
		IMatrix *m_peca;	// error covariance projection ahead
		IMatrix *m_psv;		// states vector
		IMatrix *m_psva;	// states vector projection ahead
		IMatrix *m_pgain;	// Kalman gain
	};
} // namespace