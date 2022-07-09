# Matrix Class C++

## Create matrix

```c
Matrix<float> A(3), B(3, 3);
Matrix<double> C(2, 1);
Matrix<bool> E = elementary_matrix<bool>(3);
```

## Read from console

```c
std::cin >> A;
```

## Write to console

```c
std::cout << A;
```

## Read from file

```c
std::ifstream fin("input.txt");
fin >> A;
fin.close();
```

## Write to file

```c
std::ofstream fout("output.txt");
fout << A;
fout.close();
```

## Get/Set element by row and column

```c
int value = A(2, 2);
```

## Compare matrices

```c
if (A == B) //do something
if (A != B) //do something
```

## Binary operations

```c
Matrix<float> result = A + B;
Matrix<float> result = A - B;
Matrix<float> result = A * B;
Matrix<float> result = A / B;
```

## Elementary operations

```c
// Operation_Type: Row or Column
first_elementary(number1, number2, type);
second_elementary(number, value, type);
third_elementary(number1, number2, value, type);
```

## Transposition

```c
Matrix<float> A_T = A.transposed();
```

## Determinant

```c
float result = A.det();
```

## Invertion

```c
matrix<float> A_inv = A.inverted();
```

## Rank

```c
int rank = A.rank();
```

## Triangular and Diagonal forms

```c
A.to_triangular();
A.to_diagonal();
```

## Decompompositions

```c
Matrix<float> L = A.LU_L(), U = LU_U();
A.Cholesky();
A.Frobenius();
```
