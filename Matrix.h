#pragma once
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

using namespace std;

enum class Operation_Type
{
	ROW,
	COLUMN
};

template <typename T>
class Matrix;

template <typename T>
Matrix<T> elementary_matrix(const int size)
{
	Matrix<T> result(size);
	
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			result(i, j) = T(i == j);
		}
	}
	
	return result;
}

template <typename T>
class Matrix
{
private:
	int rows, columns;
	T* values;

public:
	Matrix(const int rows, const int columns)
	{
		this->rows = rows;
		this->columns = columns;

		this->values = new T[this->rows * this->columns];
	}

	Matrix(const int size) : Matrix(size, size) {}

	Matrix(const Matrix& matrix) : Matrix(matrix.rows, matrix.columns)
	{
		for (int i = 0; i < this->rows * this->columns; i++)
		{
			this->values[i] = matrix.values[i];
		}
	}

	~Matrix()
	{
		delete[] this->values;
	}

	Matrix& operator= (const Matrix& matrix)
	{
		this->~Matrix();

		this->rows = matrix.rows;
		this->columns = matrix.columns;

		this->values = new T[this->rows * this->columns];

		for (int i = 0; i < this->rows * this->columns; i++)
		{
			this->values[i] = matrix.values[i];
		}

		return *this;
	}

	T& operator() (const int row, const int column) const
	{
		return this->values[row * this->columns + column];
	}
	
	bool operator== (const Matrix& matrix) const
	{
		if (this->rows != matrix.rows || this->columns != matrix.columns)
		{
			return false;
		}
		
		for (int i = 0; i < this->rows * this->columns; i++)
		{
			if (this->values[i] != matrix.values[i])
			{
				return false;
			}
		}

		return true;
	}

	bool operator!= (const Matrix& matrix) const
	{
		return !(*this == matrix);
	}

	Matrix operator+ (const Matrix& matrix) const
	{
		if (this->rows != matrix.rows || this->columns != matrix.columns)
		{
			throw std::exception("matrices have not equal number of rows and columns");
		}
		
		Matrix result(this->rows, this->columns);

		for (int i = 0; i < this->rows * this->columns; i++)
		{
			result.values[i] = this->values[i] + matrix.values[i];
		}

		return result;
	}

	Matrix operator- (const Matrix& matrix) const
	{
		if (this->rows != matrix.rows || this->columns != matrix.columns)
		{
			throw std::exception("matrices have not equal number of rows and columns");
		}
		
		Matrix result(this->rows, this->columns);

		for (int i = 0; i < this->rows * this->columns; i++)
		{
			result.values[i] = this->values[i] - matrix.values[i];
		}

		return result;
	}

	Matrix operator* (const Matrix& matrix) const
	{
		if (this->columns != matrix.rows)
		{
			throw std::exception("matrices have incorrect number of rows and columns");
		}
		
		Matrix result(this->rows, matrix.columns);

		for (int i = 0; i < this->rows; i++)
		{
			for (int j = 0; j < matrix.columns; j++)
			{
				result(i, j) = T();

				for (int k = 0; k < this->columns; k++)
				{
					result(i, j) += this->operator()(i, k) * matrix(k, j);
				}
			}
		}

		return result;
	}

	void operator*= (const Matrix& matrix)
	{
		*this = *this * matrix;
	}

	Matrix operator^ (const int extent) const
	{
		if (extent == 0)
		{
			return elementary_matrix<T>(this->rows);
		}

		else if (extent > 0)
		{
			Matrix result = *this;

			for (int i = 1; i < extent; i++)
			{
				result *= *this;
			}

			return result;
		}

		else
		{
			return this->inverted() ^ -extent;
		}
	}

	Matrix operator/ (const Matrix& matrix) const
	{
		return *this * matrix.inverted();
	}

	void operator*= (const T& value)
	{
		for (int i = 0; i < this->rows * this->columns; i++)
		{
			this->values[i] *= value;
		}
	}

	void operator/= (const T& value)
	{
		for (int i = 0; i < this->rows * this->columns; i++)
		{
			this->values[i] /= value;
		}
	}

	void first_elementary(const int number1, const int number2, Operation_Type type)
	{
		if (type == Operation_Type::ROW)
		{
			if (number1 < 0 || number1 >= this->rows || number2 < 0 || number2 >= this->columns)
			{
				throw std::exception("incorrect function arguments");
			}
			
			for (int i = 0; i < this->columns; i++)
			{
				swap(this->operator()(number1, i), this->operator()(number2, i));
			}
		}
		
		if (type == Operation_Type::COLUMN)
		{
			if (number1 < 0 || number1 >= this->columns || number2 < 0 || number2 >= this->columns)
			{
				throw std::exception("incorrect function arguments");
			}
			
			for (int i = 0; i < this->rows; i++)
			{
				swap(this->operator()(i, number1), this->operator()(i, number2));
			}
		}
	}

	void second_elementary(const int number, const T& value, Operation_Type type)
	{
		if (type == Operation_Type::ROW)
		{
			if (number < 0 || number >= this->rows)
			{
				throw std::exception("incorrect function arguments");
			}
			
			for (int i = 0; i < this->columns; i++)
			{
				this->operator()(number, i) *= value;
			}
		}
		
		if (type == Operation_Type::COLUMN)
		{
			if (number < 0 || number >= this->columns)
			{
				throw std::exception("incorrect function arguments");
			}
			
			for (int i = 0; i < this->rows; i++)
			{
				this->operator()(i, number) *= value;
			}
		}
	}

	void third_elementary(const int number1, const int number2, const T& value, Operation_Type type)
	{
		if (type == Operation_Type::ROW)
		{
			if (number1 < 0 || number1 >= this->rows || number2 < 0 || number1 >= this->rows)
			{
				throw std::exception("incorrect function arguments");
			}
			
			for (int i = 0; i < this->columns; i++)
			{
				this->operator()(number1, i) += value * this->operator()(number2, i);
			}
		}
		
		if (type == Operation_Type::COLUMN)
		{
			if (number1 < 0 || number1 >= this->columns || number2 < 0 || number2 >= this->columns)
			{
				throw std::exception("incorrect function arguments");
			}
			
			for (int i = 0; i < this->rows; i++)
			{
				this->operator()(i, number1) += value * this->operator()(i, number2);
			}
		}
	}

	Matrix transposed() const
	{
		Matrix result(this->columns, this->rows);

		for (int i = 0; i < this->rows; i++)
		{
			for (int j = 0; j < this->columns; j++)
			{
				result(j, i) = this->operator()(i, j);
			}
		}

		return result;
	}

	Matrix minor_matrix(const int row, const int column) const
	{
		if (row < 0 || row >= this->rows || column < 0 || column >= this->columns)
		{
			throw std::exception("incorrect function arguments");
		}
		
		Matrix result(this->rows - 1, this->columns - 1);

		int curi = 0, curj = 0;

		for (int i = 0; i < this->rows; i++)
		{
			if (i == row)
			{
				curi--;
			}

			for (int j = 0; j < this->columns; j++)
			{
				if (j == column)
				{
					curj--;
				}

				if (i != row && j != column)
				{
					result(curi, curj) = this->operator()(i, j);
				}

				curj++;
			}

			curj = 0;
			
			curi++;
		}

		return result;
	}

	T minor(const int row, const int column) const
	{
		return this->minor_matrix(row, column).det();
	}

	T algebraic_complement(const int row, const int column) const
	{
		return (((row + column) % 2 == 0) ? 1 : -1) * this->minor(row, column);
	}

	T det() const
	{
		if (this->rows != this->columns)
		{
			throw std::exception("matrix is not square");
		}
		
		if (this->rows == 1)
		{
			return this->operator()(0, 0);
		}

		else
		{
			T det = T();

			for (int i = 0; i < this->rows; i++)
			{
				det += this->operator()(i, 0) * algebraic_complement(i, 0);
			}

			return det;
		}
	}

	Matrix inverted() const
	{
		T det = this->det();

		if (det == T())
		{
			throw std::exception("matrix is not invertible");
		}
		
		Matrix result(this->rows);

		for (int i = 0; i < this->rows; i++)
		{
			for (int j = 0; j < this->columns; j++)
			{
				result(i, j) = this->algebraic_complement(j, i);
			}
		}

		result /= det;

		return result;
	}

	int rank() const
	{
		if (this->rows <= this->columns)
		{
			Matrix result = *this;
			
			int rank = this->rows;
			
			bool check;

			for (int i = 0; i < this->rows; i++)
			{
				for (int j = i + 1; j < this->rows; j++)
				{
					result.third_elementary(j, i, -result(j, i) / result(i, i), Operation_Type::ROW);
				}
			}

			for (int i = 0; i < this->rows; i++)
			{
				check = 0;
				
				for (int j = 0; j < this->columns; j++)
				{
					if (result(i, j) != T())
					{
						check = true;
						break;
					}
				}

				if (!check)
				{
					rank--;
				}
			}

			return rank;
		}

		else
		{
			return this->transposed().rank();
		}
	}

	Matrix to_triangular() const
	{
		if (this->rows != this->columns)
		{
			throw std::exception("matrix is not square");
		}
		
		Matrix result = *this;

		for (int i = 0; i < this->rows; i++)
		{
			for (int j = i + 1; j < this->rows; j++)
			{
				result.third_elementary(j, i, -result(j, i) / result(i, i), Operation_Type::ROW);
			}
		}

		return result;
	}

	Matrix to_diagonal() const
	{
		Matrix result = *this;

		result = result.to_triangular().transposed();

		return result.to_triangular();
	}

	Matrix LU_L() const
	{
		if (this->rows != this->columns)
		{
			throw std::exception("matrix is not square");
		}

		Matrix result = elementary_matrix<T>(this->rows), temp = *this;

		for (int i = 0; i < this->rows; i++)
		{
			for (int j = i + 1; j < this->rows; j++)
			{
				result(j, i) = temp(j, i) / temp(i, i);
				
				temp.third_elementary(j, i, -temp(j, i) / temp(i, i), Operation_Type::ROW);
			}
		}

		return result;
	}

	Matrix LU_U() const
	{
		return this->to_triangular();
	}

	Matrix Cholesky() const
	{
		if (*this != this->transposed())
		{
			throw std::exception("matrix is not symmetric");
		}

		Matrix result = *this;

		for (int i = 0; i < this->rows; i++)
		{
			for (int j = 0; j < this->columns; j++)
			{
				if (i == j)
				{
					for (int k = 0; k < j; k++)
					{
						result(i, j) -= pow(result(i, k), 2);
					}

					result(i, j) = sqrt(result(i, j));
				}

				else if (i > j)
				{
					for (int k = 0; k < j; k++)
					{
						result(i, j) -= result(i, k) * result(j, k);
					}

					result(i, j) /= result(j, j);
				}

				else
				{
					result(i, j) = T();
				}
			}
		}

		return result;
	}

	friend istream& operator>> <> (istream& s, Matrix& matrix);
	friend ostream& operator<< <> (ostream& s, const Matrix& matrix);

	friend ifstream& operator>> <> (ifstream& fs, Matrix& matrix);
	friend ofstream& operator<< <> (ofstream& fs, const Matrix& matrix);
};

template <typename T>
istream& operator>> (istream& s, Matrix<T>& matrix)
{
	for (int i = 0; i < matrix.rows * matrix.columns; i++)
	{
		s >> matrix.values[i];
	}
	
	return s;
}

template <typename T>
ostream& operator<< <> (ostream& s, const Matrix<T>& matrix)
{
	for (int i = 0; i < matrix.rows; i++)
	{
		for (int j = 0; j < matrix.columns; j++)
		{
			s << matrix(i, j) << " ";
		}
		
		s << endl;
	}
	
	return s;
}

template <typename T>
ifstream& operator>> (ifstream& fs, Matrix<T>& matrix)
{
	for (int i = 0; i < matrix.rows * matrix.columns; i++)
	{
		fs >> matrix.values[i];
	}
	
	return fs;
}

template <typename T>
ofstream& operator<< <> (ofstream& fs, const Matrix<T>& matrix)
{
	for (int i = 0; i < matrix.rows; i++)
	{
		for (int j = 0; j < matrix.columns; j++)
		{
			fs << matrix(i, j) << " ";
		}
		
		fs << endl;
	}
	
	return fs;
}
