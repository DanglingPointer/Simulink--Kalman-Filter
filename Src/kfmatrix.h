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
	class IMatrix
	{
	protected:
		IMatrix(uint rows, uint cols) :nrow(rows), ncol(cols)
		{ }
		IMatrix(const IMatrix&) = delete;
	public:
		virtual ~IMatrix() { }
		const uint nrow;
		const uint ncol;
		virtual double at(uint row, uint col) const = 0;
		virtual double& at(uint row, uint col) = 0;
		// Creates new dynamic matrix object
		virtual IMatrix* Copy() const = 0;
		// Creates new dynamic matrix object
		virtual IMatrix* Transpose() const = 0;
		// Creates new dynamic matrix object
		virtual IMatrix* Eye() const = 0;
		// Doesn't modify anything
		virtual double Det() const = 0;
		// Creates new dynamic matrix object
		virtual IMatrix* Inverse() const = 0;
		// Modifies calling object
		virtual IMatrix* Times(double factor) = 0;
	};
	void Print(const IMatrix *p)
	{
	#ifdef _DEBUG
		for (uint row = 0; row < p->nrow; ++row)
		{
			for (uint col = 0; col < p->ncol; ++col)
				std::cout << p->at(row, col) << ' ';
			std::cout << '\n';
		}
	#endif
	}
	template<uint n_rows, uint n_cols> class Matrix :public IMatrix
	{
	public:
		typedef Matrix<n_rows, n_cols> myt;
		Matrix(bool is_eye = false) :IMatrix(n_rows, n_cols)
		{
			for (uint i = 0; i < nrow*ncol; ++i)
				m_data[i] = 0;
			if (is_eye && nrow == ncol)
				for (uint i = 0; i < nrow; ++i)
					at(i, i) = 1;
		}
		Matrix(const myt&) = delete;
		double at(uint row, uint col) const
		{
			if (row >= nrow || col >= ncol)
				throw std::invalid_argument("Matrix::at()");
			return m_data[row*ncol + col];
		}
		double& at(uint row, uint col)
		{
			if (row >= nrow || col >= ncol)
				throw std::invalid_argument("Matrix::at()");
			return m_data[row*ncol + col];
		}
		IMatrix* Copy() const
		{
			myt *pnew = new myt;
			std::memcpy(pnew->m_data, m_data, sizeof m_data);
			return static_cast<IMatrix*>(pnew);
		}
		IMatrix* Transpose() const
		{
			IMatrix *pnew = new Matrix<n_cols, n_rows>;
			for (uint row = 0; row < pnew->nrow; ++row)
				for (uint col = 0; col < pnew->ncol; ++col)
					pnew->at(row, col) = at(col, row);
			return pnew;
		}
		IMatrix* Eye() const
		{
			if (nrow != ncol)
				throw std::invalid_argument("Matrix::Eye()");
			IMatrix *pnew = new myt(true);
			return pnew;
		}
		double Det() const;
		IMatrix* Inverse() const
		{
			IMatrix *pcf = new myt; // cofactor matrix
			for (uint row = 0; row < nrow; ++row)
				for (uint col = 0; col < ncol; ++col)
				{
					IMatrix *pminor = Minormat(row, col);
					pcf->at(row, col) = Coeff(row, col) * (pminor->Det());
					delete pminor;
				}
			IMatrix *pinv = pcf->Transpose();
			delete pcf;
			double detA = this->Det();

			for (uint row = 0; row < nrow; ++row)
				for (uint col = 0; col < ncol; ++col)
					pinv->at(row, col) *= (1 / detA);
			return pinv;
		}
		IMatrix* Times(double factor)
		{
			for (uint i = 0; i < nrow*ncol; ++i)
				m_data[i] *= factor;
			return this;
		}
	private:
		IMatrix* Minormat(uint norow, uint nocol) const;
		double m_data[n_rows*n_cols];
	};
	template<uint n_rows, uint n_cols> inline double Matrix<n_rows, n_cols>::Det() const
	{
		if (ncol != nrow)
			throw std::invalid_argument("Matrix::Det()");

		uint drow = 0;
		double det = 0;
		for (uint dcol = 0; dcol < ncol; ++dcol)
		{
			IMatrix *pnew = Minormat(drow, dcol);
			det += (Coeff(drow, dcol) * at(drow, dcol) * (pnew->Det()));
			delete pnew;
		}
		return det;
	}
	template<> inline double Matrix<1, 1>::Det() const
	{
		return at(0, 0);
	}
	// Creates a new (dynamic) matrix without norow and nocol
	template<uint n_rows, uint n_cols> inline IMatrix* Matrix<n_rows, n_cols>::Minormat(uint norow, uint nocol) const
	{
		IMatrix *pnew = new Matrix<n_rows - 1, n_cols - 1>;
		uint newrow = 0;
		for (uint row = 0; row < nrow; ++row)
			if (row != norow)
			{
				uint newcol = 0;
				for (uint col = 0; col < ncol; ++col)
					if (col != nocol)
					{
						pnew->at(newrow, newcol) = at(row, col);
						++newcol;
					}
				++newrow;
			}
		return pnew;
	}
	template<> inline IMatrix* Matrix<1, 1>::Minormat(uint norow, uint nocol) const
	{
		throw std::logic_error("Matrix<1,1>::Minormat()");
		return nullptr;
	}
	// Creates new dynamic matrix
	IMatrix* Add(const IMatrix *first, const IMatrix *second)
	{
		if (!(first->nrow == second->nrow && first->ncol == second->ncol))
			throw std::invalid_argument("Add()");
		IMatrix *p = first->Copy();
		for (uint row = 0; row < first->nrow; ++row)
			for (uint col = 0; col < first->ncol; ++col)
				p->at(row, col) = first->at(row, col) + second->at(row, col);
		return p;
	}
	// Creates new dynamic matrix
	IMatrix* Substract(const IMatrix *first, const IMatrix *second)
	{
		if (!(first->nrow == second->nrow && first->ncol == second->ncol))
			throw std::invalid_argument("Substract()");
		IMatrix *p = first->Copy();
		for (uint row = 0; row < first->nrow; ++row)
			for (uint col = 0; col < first->ncol; ++col)
				p->at(row, col) = first->at(row,col) - second->at(row, col);
		return p;
	}
	// Writes result to *result
	IMatrix* Multiply(const IMatrix *first, const IMatrix *second, IMatrix *result)
	{
		if (first->ncol != second->nrow)
			throw std::invalid_argument("Multiply(2)");
		for (uint row = 0; row < first->nrow; ++row)
			for (uint col = 0; col < second->ncol; ++col)
			{
				result->at(row, col) = 0;
				for (uint k = 0; k < first->ncol; ++k)
					result->at(row, col) += (first->at(row, k) * second->at(k, col));
			}
		return result;
	}
	// Writes result to *result
	IMatrix* Multiply(const IMatrix *first, const IMatrix *second, const IMatrix *third, IMatrix *result)
	{
		if (!(result->nrow == first->nrow && result->ncol == third->ncol))
			throw std::invalid_argument("Multiply(3)");
		class Temp
		{	// Ugly, but inevitable
		public:
			Temp(uint nrow, uint ncol):m_pdata(new double[nrow*ncol]), nrow(nrow), ncol(ncol){ }
			~Temp()	{ delete[] m_pdata; }
			double& at(uint row, uint col) { return *(m_pdata + row*ncol + col); }
			uint nrow, ncol;
		private:
			double *m_pdata;
		};
		if (first->ncol != second->nrow)
			throw std::invalid_argument("Multiply(3)");
		Temp tres(first->nrow, second->ncol);
		for (uint row = 0; row < first->nrow; ++row)
			for (uint col = 0; col < second->ncol; ++col)
			{
				tres.at(row, col) = 0;
				for (uint k = 0; k < first->ncol; ++k)
					tres.at(row, col) += (first->at(row, k) * second->at(k, col));
			}

		if (tres.ncol != third->nrow)
			throw std::invalid_argument("Multiply(3)");
		for (uint row = 0; row < tres.nrow; ++row)
			for (uint col = 0; col < third->ncol; ++col)
			{
				result->at(row, col) = 0;
				for (uint k = 0; k < tres.ncol; ++k)
					result->at(row, col) += (tres.at(row, k) * third->at(k, col));
			}
		return result;
	}
} // namespace