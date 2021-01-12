#pragma once
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

template <class T>
class Matrix
{
private:
	int rows, cols;
	T* pMatrix;

public:
	Matrix(int Rows, int Cols)
	{
		rows = Rows;
		cols = Cols;

		pMatrix = new T[rows * cols];
	}

	Matrix(int Size)
	{
		rows = Size;
		cols = Size;

		pMatrix = new T[rows * cols];
	}

	Matrix(const Matrix& matrix)
	{
		rows = matrix.rows;
		cols = matrix.cols;

		pMatrix = new T[rows * cols];

		for (int i = 0; i < rows * cols; i++)
		{
			pMatrix[i] = matrix.pMatrix[i];
		}
	}

	~Matrix()
	{
		delete[] pMatrix;
	}

	Matrix& operator= (const Matrix& matrix)
	{
		this->~Matrix();

		rows = matrix.rows;
		cols = matrix.cols;

		pMatrix = new T[rows * cols];

		for (int i = 0; i < rows * cols; i++)
		{
			pMatrix[i] = matrix.pMatrix[i];
		}

		return *this;
	}

	bool operator== (const Matrix& matrix)
	{
		if (rows == matrix.rows && cols == matrix.cols)
		{
			for (int i = 0; i < rows * cols; i++)
			{
				if (pMatrix[i] != matrix.pMatrix[i])
				{
					return false;
				}
			}
			return true;
		}
		else
		{
			return false;
		}
	}

	bool operator!= (const Matrix& matrix)
	{
		return 1 - (*this == matrix);
	}

	Matrix operator+ (const Matrix& matrix)
	{
		if (rows == matrix.rows && cols == matrix.cols)
		{
			Matrix result(rows, cols);

			for (int i = 0; i < rows * cols; i++)
			{
				result.pMatrix[i] = pMatrix[i] + matrix.pMatrix[i];
			}

			return result;
		}
	}

	Matrix operator- (const Matrix& matrix)
	{
		if (rows == matrix.rows && cols == matrix.cols)
		{
			Matrix result(rows, cols);

			for (int i = 0; i < rows * cols; i++)
			{
				result.pMatrix[i] = pMatrix[i] - matrix.pMatrix[i];
			}

			return result;
		}
	}

	Matrix operator* (const Matrix& matrix)
	{
		if (cols == matrix.rows)
		{
			Matrix result(rows, matrix.cols);

			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < matrix.cols; j++)
				{
					result.pMatrix[i * result.cols + j] = 0;
					for (int k = 0; k < cols; k++)
					{
						result.pMatrix[i * result.cols + j] += pMatrix[i * cols + k] * matrix.pMatrix[k * matrix.cols + j];
					}
				}
			}
			return result;
		}
	}

	void operator*= (const Matrix& matrix)
	{
		*this = *this * matrix;
	}

	Matrix operator^ (const int extent)
	{
		if (rows == cols)
		{
			if (extent == 0)
			{
				Matrix result(rows);

				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						result.pMatrix[i * cols + j] = i == j;
					}
				}

				return result;
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
				return this->Inv() ^ -extent;
			}
		}
	}

	void operator*= (const T& value)
	{
		for (int i = 0; i < rows * cols; i++)
		{
			pMatrix[i] *= value;
		}
	}

	void operator/= (const T& value)
	{
		for (int i = 0; i < rows * cols; i++)
		{
			pMatrix[i] /= value;
		}
	}

	void FirstE(int num1, int num2, bool isRow)
	{
		if (isRow)
		{
			if (0 < num1 && num1 <= rows && 0 < num2 && num2 <= rows)
			{
				T swap;
				for (int i = 0; i < cols; i++)
				{
					swap = pMatrix[cols * (num1 - 1) + i];
					pMatrix[cols * (num1 - 1) + i] = pMatrix[cols * (num2 - 1) + i];
					pMatrix[cols * (num2 - 1) + i] = swap;
				}
			}
		}
		else
		{
			if (0 < num1 && num1 <= cols && 0 < num2 && num2 <= cols)
			{
				T swap;
				for (int i = 0; i < rows; i++)
				{
					swap = pMatrix[cols * i + num1 - 1];
					pMatrix[cols * i + num1 - 1] = pMatrix[cols * i + num2 - 1];
					pMatrix[cols * i + num2 - 1] = swap;
				}
			}
		}
	}

	void SecondE(int num, T value, bool isRow)
	{
		if (isRow)
		{
			if (0 < num && num <= rows)
			{
				for (int i = 0; i < cols; i++)
				{
					pMatrix[cols * (num - 1) + i] *= value;
				}
			}
		}
		else
		{
			if (0 < num && num <= cols)
			{
				for (int i = 0; i < rows; i++)
				{
					pMatrix[cols * i + num - 1] *= value;
				}
			}
		}
	}

	void ThirdE(int num1, int num2, T value, bool isRow)
	{
		if (isRow)
		{
			if (0 < num1 && num1 <= rows && 0 < num2 && num2 <= rows)
			{
				for (int i = 0; i < cols; i++)
				{
					pMatrix[cols * (num1 - 1) + i] += value * pMatrix[cols * (num2 - 1) + i];
				}
			}
		}
		else
		{
			if (0 < num1 && num1 <= cols && 0 < num2 && num2 <= cols)
			{
				for (int i = 0; i < rows; i++)
				{
					pMatrix[cols * i + num1 - 1] += value * pMatrix[cols * i + num2 - 1];
				}
			}
		}
	}

	Matrix Transpose()
	{
		Matrix result(cols, rows);

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				result.pMatrix[j * rows + i] = pMatrix[i * cols + j];
			}
		}

		return result;
	}

	Matrix Minor_Matrix(int Row, int Col)
	{
		if (Row <= rows && Col <= cols)
		{
			Matrix result(rows - 1, cols - 1);

			int curi = 0, curj = 0;

			for (int i = 0; i < rows; i++)
			{
				if (i == Row - 1)
					curi--;

				for (int j = 0; j < cols; j++)
				{
					if (j == Col - 1)
						curj--;

					if (i != Row - 1 && j != Col - 1)
						result.pMatrix[curi * (cols - 1) + curj] = pMatrix[i * cols + j];

					curj++;
				}

				curj = 0;
				curi++;
			}

			return result;
		}
	}

	T Minor(int Row, int Col)
	{
		if (rows == cols && Row <= rows && Col <= cols)
		{
			return this->Minor_Matrix(Row, Col).Det();
		}
	}

	T AlgCompl(int Row, int Col)
	{
		if (rows == cols && Row <= rows && Col <= cols)
		{
			return (((Row + Col) % 2 == 0) ? 1 : -1) * this->Minor(Row, Col);
		}
	}

	T Det()
	{
		if (rows == cols)
		{
			if (rows == 1)
			{
				return pMatrix[0];
			}

			else if (rows == 2)
			{
				return pMatrix[0] * pMatrix[3] - pMatrix[1] * pMatrix[2];
			}

			else
			{
				T det = pMatrix[0] * AlgCompl(1, 1);

				for (int i = 2; i <= rows; i++)
				{
					det += pMatrix[(i - 1) * cols] * AlgCompl(i, 1);
				}

				return det;
			}
		}
	}

	Matrix Inv()
	{
		if (rows == cols)
		{
			T det = this->Det();

			if (det != 0)
			{
				Matrix result(rows);

				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						result.pMatrix[i * cols + j] = this->AlgCompl(j + 1, i + 1);
					}
				}

				result /= det;

				return result;
			}
		}
	}

	int Rank()
	{
		if (rows <= cols)
		{
			Matrix result = *this;
			int rank = rows;
			bool check;

			for (int i = 0; i < rows; i++)
			{
				for (int j = i + 1; j < rows; j++)
				{
					result.ThirdE(j + 1, i + 1, -result.pMatrix[j * cols + i] / result.pMatrix[i * (cols + 1)], 1);
				}
			}

			for (int i = 0; i < rows; i++)
			{
				check = 0;
				for (int j = 0; j < cols; j++)
				{
					if (result.pMatrix[i * cols + j] != 0)
					{
						check = 1;
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
			return this->Transpose().Rank();
		}
	}

	Matrix ToTriangular()
	{
		if (rows == cols)
		{
			Matrix result = *this;

			for (int i = 0; i < rows; i++)
			{
				for (int j = i + 1; j < rows; j++)
				{
					result.ThirdE(j + 1, i + 1, -result.pMatrix[j * cols + i] / result.pMatrix[i * (cols + 1)], 1);
				}
			}

			return result;
		}
	}

	Matrix ToDiagonal()
	{
		if (rows == cols)
		{
			Matrix result = *this;

			result = result.ToTriangular().Transpose();

			return result.ToTriangular();
		}
	}

	Matrix LU_L()
	{
		if (rows == cols)
		{
			Matrix result(rows), temp = *this;

			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < cols; j++)
				{
					result.pMatrix[i * cols + j] = i == j;
				}
			}

			for (int i = 0; i < rows; i++)
			{
				for (int j = i + 1; j < rows; j++)
				{
					result.pMatrix[j * cols + i] = temp.pMatrix[j * cols + i] / temp.pMatrix[i * (cols + 1)];
					temp.ThirdE(j + 1, i + 1, -temp.pMatrix[j * cols + i] / temp.pMatrix[i * (cols + 1)], 1);
				}
			}

			return result;
		}
	}

	Matrix LU_U()
	{
		if (rows == cols)
		{
			return this->ToTriangular();
		}
	}

	Matrix Cholesky()
	{
		if (*this == this->Transpose())
		{
			Matrix result = *this;

			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < cols; j++)
				{
					if (i == j)
					{
						for (int k = 0; k < j; k++)
						{
							result.pMatrix[i * cols + j] -= pow(result.pMatrix[i * cols + k], 2);
						}

						result.pMatrix[i * cols + j] = sqrt(result.pMatrix[i * cols + j]);
					}

					else if (i > j)
					{
						for (int k = 0; k < j; k++)
						{
							result.pMatrix[i * cols + j] -= result.pMatrix[i * cols + k] * result.pMatrix[j * cols + k];
						}

						result.pMatrix[i * cols + j] /= result.pMatrix[j * cols + j];
					}

					else
					{
						result.pMatrix[i * cols + j] = 0;
					}
				}
			}

			return result;
		}
	}

	friend istream& operator>> <> (istream& s, Matrix& matrix);
	friend ostream& operator<< <> (ostream& s, const Matrix& matrix);

	friend ifstream& operator>> <> (ifstream& fs, Matrix& matrix);
	friend ofstream& operator<< <> (ofstream& fs, const Matrix& matrix);
};

template <typename T>
istream& operator>> (istream& s, Matrix<T>& matrix)
{
	for (int i = 0; i < matrix.rows * matrix.cols; i++)
	{
		s >> matrix.pMatrix[i];
	}
	return s;
}

template <typename T>
ostream& operator<< <> (ostream& s, const Matrix<T>& matrix)
{
	for (int i = 0; i < matrix.rows; i++)
	{
		for (int j = 0; j < matrix.cols; j++)
		{
			s << matrix.pMatrix[i * matrix.cols + j] << " ";
		}
		s << endl;
	}
	return s;
}

template <typename T>
ifstream& operator>> (ifstream& fs, Matrix<T>& matrix)
{
	for (int i = 0; i < matrix.rows * matrix.cols; i++)
	{
		fs >> matrix.pMatrix[i];
	}
	return fs;
}

template <typename T>
ofstream& operator<< <> (ofstream& fs, const Matrix<T>& matrix)
{
	for (int i = 0; i < matrix.rows; i++)
	{
		for (int j = 0; j < matrix.cols; j++)
		{
			fs << matrix.pMatrix[i * matrix.cols + j] << " ";
		}
		fs << endl;
	}
	return fs;
}
