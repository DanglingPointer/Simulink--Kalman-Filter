#pragma once
//
// System:
// x[k] = A*x[x] + B*u[k] + G*w[k]
// y[k] = C*x[k] + D*u[k] + H*v[k]
//
#include"matrix.h"
namespace mvkf
{
    // Matrix typelist
    template<uint x_dim, uint u_dim, uint y_dim, uint q_dim, uint r_dim, uint w_dim, uint v_dim>
    struct Matlist
    {
        typedef Matrix<x_dim, x_dim> MatA;
        typedef Matrix<x_dim, x_dim> MatP;    // ???
        typedef Matrix<x_dim, u_dim> MatB;
        typedef Matrix<x_dim, w_dim> MatG;
        typedef Matrix<y_dim, x_dim> MatC;
        typedef Matrix<x_dim, y_dim> MatK;
        typedef Matrix<y_dim, u_dim> MatD;
        typedef Matrix<y_dim, v_dim> MatH;
        typedef Matrix<q_dim, q_dim> MatQ;
        typedef Matrix<r_dim, r_dim> MatR;
        typedef Matrix<x_dim, 1> VecX;
        typedef Matrix<u_dim, 1> VecU;
        typedef Matrix<y_dim, 1> VecY;
    };
    // To be implemented by user; 'tlist_t' must be Matlist
    template<class tlist_t> class ISysmat
    {
    public:
        typedef typename tlist_t::MatA MatA;
        typedef typename tlist_t::MatB MatB;
        typedef typename tlist_t::MatC MatC;
        typedef typename tlist_t::MatD MatD;
        typedef typename tlist_t::MatG MatG;
        typedef typename tlist_t::MatH MatH;
        typedef typename tlist_t::MatR MatR;
        typedef typename tlist_t::MatQ MatQ;
        typedef typename tlist_t::MatP MatP;
        typedef typename tlist_t::VecX VecX;

        virtual const MatA& get_A() const = 0;
        virtual const MatB& get_B() const = 0;
        virtual const MatG& get_G() const = 0;
        virtual const MatC& get_C() const = 0;
        virtual const MatD& get_D() const = 0;
        virtual const MatH& get_H() const = 0;
        virtual const MatQ& get_Q() const = 0;
        virtual const MatR& get_R() const = 0;
        virtual const VecX& get_InitSV() const = 0;
        virtual const MatP& get_InitCov() const = 0;
        virtual void set_States(const VecX *states)
        { }
        virtual void set_Time(double time)
        { }
    };
    template<class tlist_t> class Filter
    {
    public:
        typedef typename tlist_t::MatA MatA;
        typedef typename tlist_t::MatB MatB;
        typedef typename tlist_t::MatC MatC;
        typedef typename tlist_t::MatD MatD;
        typedef typename tlist_t::MatG MatG;
        typedef typename tlist_t::MatH MatH;
        typedef typename tlist_t::MatR MatR;
        typedef typename tlist_t::MatQ MatQ;
        typedef typename tlist_t::MatP MatP;
        typedef typename tlist_t::MatK MatK;
        typedef typename tlist_t::VecX VecX;
        typedef typename tlist_t::VecY VecY;
        typedef typename tlist_t::VecU VecU;

        Filter(ISysmat<tlist_t> *sysmat) :m_psysmat(sysmat), m_eca(m_psysmat->get_InitCov()), m_sva(m_psysmat->get_InitSV())
        {
            m_psysmat->set_States(&m_sv);
        }
        Filter(const Filter<tlist_t>&) = delete;
        // p.146 eq(4.2.17)
        void ComputeGain()
        {
            MatC C = m_psysmat->get_C();
            MatR R = m_psysmat->get_R();
            m_gain = (m_eca * T(C)) * Inv(((C * m_eca) * T(C)) + R);
            Print(m_gain);
        }
        // p.144 eq(4.2.8)
        void UpdateEstimate(const VecY& y)
        {
            MatC C = m_psysmat->get_C();
            m_sv = m_sva + (m_gain * (y - (C * m_sva)));
            Print(m_sv);
        }
        // p.146 eq(4.2.18)
        void ComputeCovariance()
        {
            MatC C = m_psysmat->get_C();
            MatR R = m_psysmat->get_R();
            MatP I = I(m_eca);
            MatP parenth = I - (m_gain*C);
            MatP prod = (m_gain * R) * T(m_gain);
            m_ec = ((parenth * m_eca) * T(parenth)) + prod;
            Print(m_ec);
        }
        // assignment 5 prob.2 f) and g)
        void ProjectAhead(const VecU& u)
        {
            MatA A = m_psysmat->get_A();
            MatQ Q = m_psysmat->get_Q();
            MatB B = m_psysmat->get_B();
            MatG G = m_psysmat->get_G();

            m_sva = (A * m_sv) + (B * u);
            m_eca = ((A * m_ec) * T(A)) + ((G * Q) * T(G));
            Print(m_sva); Print(m_eca);
        }
        // Initial states from arguments
        void InitStates(const VecX &initsv, const MatP &initcov)
        {
            m_eca = initcov;
            m_sva = initsv;
        }
        void set_Time(double time)
        {
            m_psysmat->set_Time(time);
        }
        // First state is no.0
        double get_State(uint statenum) const
        {
            return m_sv(statenum, 0);
        }
        const VecX& get_Statevec() const
        {
            return m_sv;
        }
    private:
        ISysmat<tlist_t> *m_psysmat;
        MatP   m_ec;       // error covariance
        MatP   m_eca;      // error covariance projection ahead
        VecX   m_sv;       // states vector
        VecX   m_sva;      // states vector projection ahead
        MatK   m_gain;     // Kalman gain
    };
}