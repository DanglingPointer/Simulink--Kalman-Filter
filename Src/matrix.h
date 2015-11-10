#pragma once
// Generic matrix implementation using templates
#include<iostream>
#include<cstring>
namespace mvkf
{
    typedef unsigned int uint;
    template<uint nrow, uint ncol> class Matrix
    {
    public:
        typedef Matrix<nrow, ncol> myt;
        enum
        {
            NUMROW = nrow,
            NUMCOL = ncol
        };
        Matrix()
        {
            for (uint i = 0; i < nrow*ncol; ++i)
                m_data[i] = 0;
        }
        Matrix(const myt& rhs)
        {
            std::memcpy(m_data, rhs.m_data, sizeof m_data);
        }
        double& operator()(uint row, uint col)
        {
            return m_data[row*ncol + col];
        }
        const double& operator()(uint row, uint col) const
        {
            return m_data[row*ncol + col];
        }
        template<uint col2> Matrix<nrow, col2> operator*(const Matrix<ncol, col2>& rhs) const
        {   // row2 = col1
            Matrix<nrow, col2> result;
            for (uint row = 0; row < nrow; ++row)
                for (uint col = 0; col < col2; ++col)
                    for (uint k = 0; k < ncol; ++k)
                        result(row, col) += (operator()(row, k) * rhs(k, col));
            return result;
        }
        myt& operator*=(double factor)
        {
            for (uint i = 0; i < nrow*ncol; ++i)
                m_data[i] *= factor;
            return *this;
        }
        myt operator*(double factor) const
        {
            myt temp(*this);
            return temp *= factor;
        }
        myt operator+(const myt& rhs) const
        {
            myt result(*this);
            for (uint i = 0; i < nrow*ncol; ++i)
                result.m_data[i] += rhs.m_data[i];
            return result;
        }
        myt operator-(const myt& rhs) const
        {
            myt result(*this);
            for (uint i = 0; i < nrow*ncol; ++i)
                result.m_data[i] -= rhs.m_data[i];
            return result;
        }
    private:
        double m_data[nrow*ncol];
    };
    template<uint nrow, uint ncol> Matrix<nrow, ncol> operator*(double factor, const Matrix<nrow, ncol>& rhs)
    {
        return rhs*factor;
    }

#define Coeff(j,k) ( ((j+k) % 2 == 0) ? 1:-1 )
#define Det(mat) MatOperation::Determinant(mat)
#define Inv(mat) MatOperation::Inverse(mat)
#define T(mat) MatOperation::Transpose(mat)
#define I(mat) MatOperation::Identity(mat)
#ifdef _DEBUG
#define Print(mat) MatOperation::SendToStream(mat, std::cout)
#else
#define Print(mat)
#endif
    // Matrix operations
    class MatOperation
    {
    public:
        template<uint size>
        static double Determinant(const Matrix<size, size>& mat)
        {
            uint drow = 0;
            double det = 0;
            for (uint dcol = 0; dcol < size; ++dcol)
            {
                Matrix<size - 1, size - 1> min = Minormat(mat, drow, dcol);
                det += (Coeff(drow, dcol) * mat(drow, dcol) * Determinant(min));
            }
            return det;
        }
        static double Determinant(const Matrix<1, 1>& mat)
        {
            return mat(0, 0);
        }
        template<uint size>
        static Matrix<size, size> Inverse(const Matrix<size, size>& mat)
        {
            Matrix<size, size> cf; // cofactor matrix
            for (uint row = 0; row < size; ++row)
                for (uint col = 0; col < size; ++col)
                {
                    Matrix<size - 1, size - 1> min = Minormat(mat, row, col);
                    cf(row, col) = Coeff(row, col) * Determinant(min);
                }
            Matrix<size, size> inv = T(cf);
            double detA = Determinant(mat);
            for (uint row = 0; row < size; ++row)
                for (uint col = 0; col < size; ++col)
                    inv(row, col) *= (1 / detA);
            return inv;
        }
        static Matrix<1, 1> Inverse(const Matrix<1, 1>& mat)
        {
            Matrix<1, 1> temp;
            temp(0, 0) = 1.0 / mat(0, 0);
            return temp;
        }
        template<uint nrow, uint ncol>
        static Matrix<ncol, nrow> Transpose(const Matrix<nrow, ncol>& mat)
        {
            Matrix<ncol, nrow> temp;
            for (uint row = 0; row < nrow; ++row)
                for (uint col = 0; col < ncol; ++col)
                    temp(col, row) = mat(row, col);
            return temp;
        }
        template<uint size>
        static Matrix<size, size> Identity(const Matrix<size, size>& mat)
        {
            Matrix<size, size> temp;
            for (uint i = 0; i < size; ++i)
                temp(i, i) = 1;
            return temp;
        }
        template<uint nrow, uint ncol>
        static void SendToStream(const Matrix<nrow, ncol>& mat, std::ostream& out)
        {
            for (uint row = 0; row < nrow; ++row)
            {
                for (uint col = 0; col < ncol; ++col)
                    out << mat(row, col) << ' ';
                out << '\n';
            }
        }
    private:
        template<uint size>
        static Matrix<size - 1, size - 1> Minormat(const Matrix<size, size>& mat, uint norow, uint nocol)
        {
            Matrix<size - 1, size - 1> temp;
            uint newrow = 0;
            for (uint row = 0; row < size; ++row)
                if (row != norow)
                {
                    uint newcol = 0;
                    for (uint col = 0; col < size; ++col)
                        if (col != nocol)
                        {
                            temp(newrow, newcol) = mat(row, col);
                            ++newcol;
                        }
                    ++newrow;
                }
            return temp;
        }
        static Matrix<1, 1> Minormat(const Matrix<2, 2>& mat, uint norow, uint nocol)
        {
            Matrix<1, 1> temp;
            uint row = (norow == 1) ? 0 : 1;
            uint col = (nocol == 1) ? 0 : 1;
            temp(0, 0) = mat(row, col);
            return temp;
        }
    };
}