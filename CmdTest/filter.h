#pragma once
#include"matrix.h"
namespace mvkf
{
	class Filter
	{
	public:
		Filter(IMatrix*(*phi)(unsigned long), IMatrix*(*h)(unsigned long), IMatrix*(*q)(unsigned long), IMatrix*(*r)(unsigned long))
			:m_time(0), m_phi(phi), m_h(h), m_q(q), m_r(r)
		{ }
		~Filter()
		{
			// delete all dynamic matrices
		}
		// Filter loop functions
		void ComputeGain()
		{
			//
		}
		void ProjectAhead()
		{
			//
		}
		void UpdateEstimate()
		{
			//
		}
		void ComputeCovariance()
		{
			//
		}
		// Mutators
		void InitStates(IMatrix *pinitsv, IMatrix *pinitec)
		{
			m_peca = pinitec;
			m_psva = pinitsv;
		}
		void TimeStep() // possibly in ProjectAhead() instead
		{
			++m_time;
		}
		void set_Measurement(uint measnum, double value) // possibly in UpdateEstimate() instead
		{
			//
		}
		// Accessors
		double get_State(uint statenum) const // statenum from 1
		{
			return m_psv->at(statenum - 1, 0);
		}
	private:
		unsigned long m_time;
		// Functions for computation of system matrices at time k
		IMatrix*(*m_phi)(unsigned long);
		IMatrix*(*m_h)(unsigned long);
		IMatrix*(*m_q)(unsigned long);
		IMatrix*(*m_r)(unsigned long);
		// Data
		IMatrix* m_pec;		// error covariance
		IMatrix* m_peca;	// error covariance projection ahead
		IMatrix* m_psv;		// states vector
		IMatrix* m_psva;	// states vector projection ahead
	};
} // namespace