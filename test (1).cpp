#include <iostream>
#include "Matrix.hpp"
using namespace std;

int main() {
    // Создание пустой матрицы
    Matrix<int> emptyMatrix;
    cout << "Empty Matrix:" <<endl;
    cout << emptyMatrix << endl;

    // Создание матрицы заданных размеров
    Matrix<int> sizedMatrix(3, 3);
    cout << "Sized Matrix (3x3):" << endl;
    cout << sizedMatrix << endl;

    // Создание матрицы из вектора векторов
    vector<vector<int>> data = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Matrix<int> vectorMatrix(data);
    cout << "Vector Matrix:" << endl;
    cout << vectorMatrix << endl;

    // Копирование матрицы
    Matrix<int> copiedMatrix(vectorMatrix);
    cout << "Copied Matrix:" << endl;
    cout << copiedMatrix << endl;

    // Тестирование оператора присваивания
    copiedMatrix = sizedMatrix;
    cout << "Copied Matrix (after =):" << endl;
    cout << copiedMatrix << endl;

    // Доступ к элементам матрицы
    cout << "Element at (1, 2) of sizedMatrix: " << sizedMatrix(1, 2) << endl;

    // Получение количества строк и столбцов
    cout << "Number of rows in sizedMatrix: " << sizedMatrix.get_rows() << endl;
    cout << "Number of cols in sizedMatrix: " << sizedMatrix.get_cols() << endl;

    // Сравнение матриц
    cout << "sizedMatrix == vectorMatrix: " << (sizedMatrix == vectorMatrix) << endl;
    cout << "sizedMatrix != copiedMatrix: " << (sizedMatrix != copiedMatrix) << endl;

    // Умножение матрицы на скаляр
    Matrix<int> scaledMatrix = vectorMatrix * 2;
    cout << "Scaled Matrix:" << endl;
    cout << scaledMatrix << endl;

    // Сложение матриц
    Matrix<int> sumMatrix = sizedMatrix + vectorMatrix;
    cout << "Sum Matrix:" << endl;
    cout << sumMatrix << endl;

    // Вычитание матриц
    Matrix<int> diffMatrix = sizedMatrix - vectorMatrix;
    cout << "Difference Matrix:" << endl;
    cout << diffMatrix << endl;

    // Умножение матриц
    Matrix<int> mulMatrix = vectorMatrix * sumMatrix;
    cout << "Matrix multiplication:" << endl;
    cout << mulMatrix << endl;

  
    vector<vector<float>> data3 = {{1, 0, 3}, {0, 5, 6}, {7, 8, 0}};
    Matrix<float> floatMatrix2(data3);

    // Вычисление определителя
    cout << "Determinant of floatMatrix: " << floatMatrix2.determinant() << endl;
    cout << "Determinant of floatMatrix: " << floatMatrix2.deter() << endl;
    // Вычисление обратной матрицы
    Matrix<float> floatMatrix3(data3);
    Matrix<float> inverseMatrix3 = floatMatrix3.inverse2();
    cout << "Inverse Matrix:" << endl;
    cout << inverseMatrix3 << endl;

    // Вычисление ранга матрицы
    vector<vector<float>> data4 = {{2, -1, 3, -2, 4}, {4, -2, 5, 1, 7}, {2, -1, 1, 8, 2}};
    Matrix<float> rankMatrix(data4);
    cout << "Rank Matrix:" << endl;
    cout << rankMatrix.ranking() << endl;

    // Создание единичной матрицы
    Matrix<int> identityMatrix = Matrix<int>::eye(3);
    cout << "Identity Matrix:" << endl;
    cout << identityMatrix << endl;

    // Вычисление обратной матрицы
    vector<vector<float>> data2 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Matrix<float> floatMatrix(data2);
    Matrix<float> inverseMatrix = floatMatrix.inverse2();
    cout << "Inverse Matrix:" << endl;
    cout << inverseMatrix << endl;
    return 0;
}