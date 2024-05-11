#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Tian Putra Perdana
// 21120122130072

// Fungsi LU decomposition
void luDecomposition(const vector<vector<double>> &A, vector<vector<double>> &L, vector<vector<double>> &U, int n)
{
    // Decomposition
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j < i)
                L[j][i] = 0;
            else
            {
                L[j][i] = A[j][i];
                for (int k = 0; k < i; k++)
                    L[j][i] -= L[j][k] * U[k][i];
            }
        }
        for (int j = 0; j < n; j++)
        {
            if (j < i)
                U[i][j] = 0;
            else if (j == i)
                U[i][j] = 1;
            else
            {
                U[i][j] = A[i][j] / L[i][i];
                for (int k = 0; k < i; k++)
                    U[i][j] -= ((L[i][k] * U[k][j]) / L[i][i]);
            }
        }
    }
}

// Fungsi Crout decomposition
void croutDecomposition(const vector<vector<double>> &A, vector<vector<double>> &L, vector<vector<double>> &U, int n)
{
    // Decomposition
    for (int i = 0; i < n; i++)
    {
        U[i][i] = 1;

        // Uupper triangular matrix
        for (int k = i; k < n; k++)
        {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += (L[i][j] * U[j][k]);
            U[i][k] = A[i][k] - sum;
        }

        // Lower triangular matrix
        for (int k = i + 1; k < n; k++)
        {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += (L[k][j] * U[j][i]);
            L[k][i] = (A[k][i] - sum) / U[i][i];
        }
    }
}

// Fungsi untuk menyelesaikan sistem persamaan linear dengan metode LU decomposition
vector<double> solveEquationsLU(const vector<vector<double>> &A, const vector<double> &B)
{
    int n = A.size();

    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));

    luDecomposition(A, L, U, n);

    vector<double> Y(n, 0);
    vector<double> X(n, 0);

    // Solve LY = B
    for (int i = 0; i < n; i++)
    {
        Y[i] = B[i];
        for (int j = 0; j < i; j++)
        {
            Y[i] -= L[i][j] * Y[j];
        }
        Y[i] /= L[i][i];
    }

    // Solve UX = Y
    for (int i = n - 1; i >= 0; i--)
    {
        X[i] = Y[i];
        for (int j = i + 1; j < n; j++)
        {
            X[i] -= U[i][j] * X[j];
        }
        X[i] /= U[i][i];
    }

    return X;
}

// Fungsi untuk menyelesaikan sistem persamaan linear dengan metode Crout decomposition
vector<double> solveEquationsCrout(const vector<vector<double>> &A, const vector<double> &B)
{
    int n = A.size();

    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));

    croutDecomposition(A, L, U, n);

    vector<double> Y(n, 0);
    vector<double> X(n, 0);

    // Solve LY = B
    for (int i = 0; i < n; i++)
    {
        Y[i] = B[i];
        for (int j = 0; j < i; j++)
        {
            Y[i] -= L[i][j] * Y[j];
        }
    }

    // Solve UX = Y
    for (int i = n - 1; i >= 0; i--)
    {
        X[i] = Y[i];
        for (int j = i + 1; j < n; j++)
        {
            X[i] -= U[i][j] * X[j];
        }
        X[i] /= U[i][i];
    }

    return X;
}

// Fungsi untuk mencari invers matrix
vector<vector<double>> inverseMatrix(const vector<vector<double>> &matrix)
{
    int n = matrix.size();

    // Identity matrix
    vector<vector<double>> identity(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i)
        identity[i][i] = 1;

    vector<vector<double>> A = matrix;

    // Forward elimination
    for (int i = 0; i < n; ++i)
    {
        double pivot = A[i][i];
        for (int j = 0; j < n; ++j)
        {
            A[i][j] /= pivot;
            identity[i][j] /= pivot;
        }
        for (int k = i + 1; k < n; ++k)
        {
            double factor = A[k][i];
            for (int j = 0; j < n; ++j)
            {
                A[k][j] -= factor * A[i][j];
                identity[k][j] -= factor * identity[i][j];
            }
        }
    }

    // Backward elimination
    for (int i = n - 1; i > 0; --i)
    {
        for (int k = i - 1; k >= 0; --k)
        {
            double factor = A[k][i];
            for (int j = 0; j < n; ++j)
            {
                A[k][j] -= factor * A[i][j];
                identity[k][j] -= factor * identity[i][j];
            }
        }
    }

    return identity;
}

// Fungsi untuk mengalikan 2 matrix
vector<double> matrixMultiply(const vector<vector<double>> &A, const vector<double> &B)
{
    int n = A.size();
    vector<double> result(n, 0);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result[i] += A[i][j] * B[j];
        }
    }
    return result;
}

// Fungsi untuk menyelesaikan sistem persamaan linear menggunakan metode invers matrix
vector<double> solveLinearSystem(const vector<vector<double>> &A, const vector<double> &B)
{
    vector<vector<double>> invA = inverseMatrix(A);
    return matrixMultiply(invA, B);
}

int main()
{
    vector<vector<double>> A = {{1, -1, 2},
                                {3, 0, 1},
                                {1, 0, 2}};
    vector<double> B = {5, 10, 5};

    cout << "Menggunakan metode LU decomposition:" << endl;
    vector<double> X_LU = solveEquationsLU(A, B);
    for (size_t i = 0; i < X_LU.size(); ++i)
    {
        cout << "X[" << i << "] = " << X_LU[i] << endl;
    }

    cout << "\nMenggunakan metoder Crout decomposition:" << endl;
    vector<double> X_Crout = solveEquationsCrout(A, B);
    for (size_t i = 0; i < X_Crout.size(); ++i)
    {
        cout << "X[" << i << "] = " << X_Crout[i] << endl;
    }

    cout << "\nMenggunakan metode inverse matrix method:" << endl;
    vector<double> solution = solveLinearSystem(A, B);
    for (size_t i = 0; i < solution.size(); ++i)
    {
        cout << "x[" << i << "] = " << solution[i] << endl;
    }

    return 0;
}
