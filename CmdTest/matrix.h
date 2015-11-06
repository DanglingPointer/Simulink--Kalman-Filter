#pragma once
#include<cstring>
#include<stdexcept>
#ifdef _DEBUG
 #include<iostream>
#endif
#define Coeff(j,k) ( ((j+k) % 2 == 0) ? 1:-1 )

typedef unsigned int uint;
class IMatrix
{
protected:
	IMatrix(uint rows, uint cols):nrow(rows), ncol(cols)
	{ }
	IMatrix(const IMatrix&) = delete;
public:
	virtual ~IMatrix(){ }
	const uint nrow;
	const uint ncol;
	virtual double at(uint row, uint col) const = 0;
	virtual double& at(uint row, uint col) = 0;
	virtual IMatrix* Copy() const = 0;				// creates new dynamic matrix object
	virtual IMatrix* Transpose() const = 0;			// creates new dynamic matrix object
	virtual double Det() const = 0;
	//virtual IMatrix* Invert() const = 0;			// creates new dynamic matrix object
};
#ifdef _DEBUG
void Print(const IMatrix *p)
{
	for (uint row = 0; row < p->nrow; ++row)
	{
		for (uint col = 0; col < p->ncol; ++col)
			std::cout << p->at(row, col) << ' ';
		std::cout << '\n';
	}
} 
#endif
template<uint n_rows, uint n_cols> class Matrix :public IMatrix
{
public:
	typedef Matrix<n_rows, n_cols> myt;
	Matrix(bool is_eye = false) :IMatrix(n_rows, n_cols)
	{
		for (uint i = 0; i < nrow*ncol; ++i)
			m_data[i] = 0;
		if (nrow == ncol && is_eye)
			for (uint i = 0; i < nrow; ++i)
				at(i, i) = 1;

	}
	Matrix(const myt&) = delete;
	double at(uint row, uint col) const
	{
		return m_data[row*ncol + col];
	}
	double& at(uint row, uint col)
	{
		return m_data[row*ncol + col];
	}
	IMatrix* Copy() const
	{
		myt *pnew = new myt;
		std::memcpy(pnew->m_data, m_data, n_rows*n_cols*sizeof(*m_data));
		return static_cast<IMatrix*>(pnew);
	}
	IMatrix* Transpose() const
	{
		Matrix<n_cols, n_rows> *pnew = new Matrix<n_cols, n_rows>;
		for (uint row = 0; row < pnew->nrow; ++row)
			for (uint col = 0; col < pnew->ncol; ++col)
				pnew->at(row, col) = at(col, row);
		return static_cast<IMatrix*>(pnew);
	}
	double Det() const;
private:
	double m_data[n_rows*n_cols];
};
template<uint n_rows, uint n_cols> inline double Matrix<n_rows, n_cols>::Det() const
{
	if (ncol != nrow)
		throw std::invalid_argument("Determinant of a non-square matrix");

	uint drow = 0;
	double det = 0;
	for (uint dcol = 0; dcol < ncol; ++dcol)
	{
		IMatrix *pnew = new Matrix<n_rows - 1, n_cols - 1>;
		uint newrow = 0, newcol = 0;
		for (uint row = 0; row < nrow; ++row)
			if (row != drow)
			{
				for (uint col = 0; col < ncol; ++col)
					if (col != dcol)
					{
						pnew->at(newrow, newcol) = at(row, col);
						++newcol;
					}
				++newrow;
			}
		Print(this);
		det += (Coeff(drow, dcol) * at(drow, dcol) * (pnew->Det()));
		delete pnew;
	}
	return det;
}
template<> inline double Matrix<1, 1>::Det() const
{
	Print(this);
	return at(0, 0);
}