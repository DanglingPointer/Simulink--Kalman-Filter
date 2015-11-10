#pragma once
#include"filter.h"
using namespace mvkf;

typedef Matlist<5, 1, 1, 2, 1> Dim;

class Sysmat :public ISysmat<Dim>
{
public:
    Sysmat()
    {
        m_a(0, 0) = 0.9969;
        m_a(0, 1) = 0.09982;
        m_a(1, 0) = -0.0611;
        m_a(1, 1) = 0.9955;
        m_a(2, 2) = 1.0;
        m_a(2, 3) = 0.09994;
        m_a(2, 4) = -9.965e-06;
        m_a(3, 3) = 0.9989;
        m_a(3, 4) = -0.0001993;
        m_a(4, 4) = 1.0;

        m_b(2, 0) = 9.965e-06;
        m_b(3, 0) = 0.0001993;

        m_c(0, 1) = 1;
        m_c(0, 2) = 1;

        m_d(0, 0) = 0.0;

        m_g(0, 0) = 2.847e-06;
        m_g(1, 0) = 5.691e-05;
        m_g(2, 1) = -3.322e-07;
        m_g(3, 1) = -9.965e-06;
        m_g(4, 1) = 0.1;

        m_h(0, 0) = 1.0;
        
        m_p(0, 0) = 1.0;
        m_p(1, 1) = 0.013;
        m_p(2, 2) = 3.14159265359*3.14159265359;
        m_p(3, 3) = 1.0;
        m_p(4, 4) = 2.5e-04;

        m_q(0, 0) = 30;
        m_q(1, 1) = 1e-06;

        m_r(0, 0) = 1.0002e-08;			// = (var(v) / simulation_time)
    }
    virtual const MatA& get_A() const
    {
        return m_a;
    }
    virtual const MatB& get_B() const
    {
        return m_b;
    }
    virtual const MatG& get_G() const
    {
        return m_g;
    }
    virtual const MatC& get_C() const
    {
        return m_c;
    }
    virtual const MatD& get_D() const
    {
        return m_d;
    }
    virtual const MatH& get_H() const
    {
        return m_h;
    }
    virtual const MatQ& get_Q() const
    {
        return m_q;
    }
    virtual const MatR& get_R() const
    {
        return m_r;
    }
    virtual const VecX& get_InitSV() const
    {
        return m_x;
    }
    virtual const MatP& get_InitCov() const
    {
        return m_p;
    }
private:
    MatA m_a;
    MatB m_b;
    MatC m_c;
    MatD m_d;
    MatG m_g;
    MatH m_h;
    MatQ m_q;
    MatR m_r;
    MatP m_p;
    VecX m_x;
};