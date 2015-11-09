#pragma once
#pragma once
#include<cstring>
#include<stdexcept>
#ifdef _DEBUG
#include<iostream>
#endif
#define Coeff(j,k) ( ((j+k) % 2 == 0) ? 1:-1 )
namespace mvkf
{
	typedef unsigned int uint;
	class Matrix
	{
	public:
		Matrix(uint rows, uint cols) :nrow(rows), ncol(cols), m_pdata(new double[nrow*ncol]{ })
		{ }
		Matrix(const Matrix& rhs) :Matrix(rhs.nrow, rhs.ncol)
		{
			std::memcpy(m_pdata, rhs.m_pdata, nrow*ncol*sizeof(*(rhs.m_pdata)));
		}
		virtual ~Matrix()
		{
			delete[]m_pdata;
		}
		const uint nrow;
		const uint ncol;
		double at(uint row, uint col) const
		{
			if (row >= nrow || col >= ncol)
				throw std::invalid_argument("Matrix::at()");
			return m_pdata[row*ncol + col];
		}
		double& at(uint row, uint col)
		{
			if (row >= nrow || col >= ncol)
				throw std::invalid_argument("Matrix::at()");
			return m_pdata[row*ncol + col];
		}
		Matrix& operator=(const Matrix& rhs)
		{
			if (!(nrow == rhs.nrow && ncol == rhs.ncol))
				throw std::invalid_argument("Matrix::operator=()");
			std::memcpy(m_pdata, rhs.m_pdata, nrow*ncol*sizeof(*(rhs.m_pdata)));
			return *this;
		}
		Matrix Transpose() const
		{
			Matrix temp(ncol, nrow);
			for (uint row = 0; row < temp.nrow; ++row)
				for (uint col = 0; col < temp.ncol; ++col)
					temp.at(row, col) = at(col, row);
			return temp;
		}
		Matrix Eye() const
		{
			if (nrow != ncol)
				throw std::invalid_argument("Matrix::Eye()");
			Matrix temp(nrow, ncol);
			for (uint i = 0; i < nrow; ++i)
				temp.at(i, i) = 1;
			return temp;
		}
		double Det() const
		{
			if (ncol != nrow)
				throw std::invalid_argument("Matrix::Det()");
			if (ncol == 1 && nrow == 1)
				return at(0, 0);

			uint drow = 0;
			double det = 0;
			for (uint dcol = 0; dcol < ncol; ++dcol)
			{
				Matrix min = Minormat(drow, dcol);
				det += (Coeff(drow, dcol) * at(drow, dcol) * (min.Det()));
			}
			return det;
		}
		Matrix Inverse() const
		{
			if (nrow != ncol)
				throw std::invalid_argument("Matrix::Inverse()");


			Matrix cf(*this); // cofactor matrix
			if (nrow == 1 && ncol == 1)
			{
				cf.at(0, 0) = 1.0/cf.at(0, 0);
				return cf;
			}
			for (uint row = 0; row < nrow; ++row)
				for (uint col = 0; col < ncol; ++col)
				{
					Matrix min = Minormat(row, col);
					cf.at(row, col) = Coeff(row, col) * (min.Det());
				}
			Matrix inv = cf.Transpose();
			double detA = this->Det();
			for (uint row = 0; row < nrow; ++row)
				for (uint col = 0; col < ncol; ++col)
					inv.at(row, col) *= (1 / detA);
			return inv;
		}
		Matrix& operator+=(const Matrix& rhs)
		{
			if (!(nrow == rhs.nrow && ncol == rhs.ncol))
				throw std::invalid_argument("Matrix::operator+=()");

			for (uint i = 0; i < nrow*ncol; ++i)
				*(m_pdata + i) += *(rhs.m_pdata + i);
			return *this;
		}
		Matrix& operator-=(const Matrix& rhs)
		{
			if (!(nrow == rhs.nrow && ncol == rhs.ncol))
				throw std::invalid_argument("Matrix::operator+=()");

			for (uint i = 0; i < nrow*ncol; ++i)
				*(m_pdata + i) -= *(rhs.m_pdata + i);
			return *this;
		}
		Matrix operator+(const Matrix& rhs) const
		{
			if (!(nrow == rhs.nrow && ncol == rhs.ncol))
				throw std::invalid_argument("Matrix::operator+()");
			Matrix temp(*this);
			return temp += rhs;
		}
		Matrix operator-(const Matrix& rhs) const
		{
			if (!(nrow == rhs.nrow && ncol == rhs.ncol))
				throw std::invalid_argument("Matrix::operator-()");
			Matrix temp(*this);
			return temp -= rhs;
		}
		Matrix operator*(const Matrix& rhs) const
		{
			if (ncol != rhs.nrow)
				throw std::invalid_argument("Matrix::operator*()");

			Matrix result(nrow, rhs.ncol);
			for (uint row = 0; row < nrow; ++row)
				for (uint col = 0; col < rhs.ncol; ++col)
					for (uint k = 0; k < ncol; ++k)
						result.at(row, col) += (at(row, k) * rhs.at(k, col));
			return result;
		}
	private:
		Matrix Minormat(uint norow, uint nocol) const
		{
			if (ncol != nrow || ncol < 2)
				throw std::invalid_argument("Matrix::Minormat()");
			Matrix temp(nrow - 1, ncol - 1);
			uint newrow = 0;
			for (uint row = 0; row < nrow; ++row)
				if (row != norow)
				{
					uint newcol = 0;
					for (uint col = 0; col < ncol; ++col)
						if (col != nocol)
						{
							temp.at(newrow, newcol) = at(row, col);
							++newcol;
						}
					++newrow;
				}
			return temp;
		}
		double* m_pdata;
	};
	void Print(const Matrix& mat)
	{
	#ifdef _DEBUG
		for (uint row = 0; row < mat.nrow; ++row)
		{
			for (uint col = 0; col < mat.ncol; ++col)
				std::cout << mat.at(row, col) << ' ';
			std::cout << '\n';
		}
	#endif
	}
} // namespace