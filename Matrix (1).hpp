#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
using namespace std;

template<typename T> //T означает, что мы можем использовать любой тип данных
class Matrix {
private:
size_t rows;// количество строк
size_t cols;// количество столбцов
vector<vector<T>> data; // данные матрицы

public:
    // Конструктор пустой матрицы
    Matrix() : rows(0), cols(0) {}
    // Конструктор матрицы заданных размеров
    Matrix(size_t rows, size_t cols) : rows(rows), cols(cols), data(vector<vector<T>>(rows, vector<T>(cols))) {}
    // Конструктор матрицы из вектора векторов
    Matrix(const vector<vector<T>>& v) : rows(v.size()), cols(v[0].size()), data(v) {}
    // Конструктор копирования
    Matrix(const Matrix<T>& other) : rows(other.rows), cols(other.cols), data(other.data) {}

    // Оператор присваивания
    Matrix<T>& operator=(const Matrix<T>& other){
      if(this == &other){
        return *this;
      } else{
        rows = other.rows;
        cols = other.cols;
        data = other.data;
        return *this;
      }
    }

    // Доступ к элементу (i, j)
    T& operator()(size_t i, size_t j){
      if(i >= rows || j >= cols){
        throw out_of_range("Index out of range"); // вызов исключения
      }
      return data[i][j];
    }
    // Доступ к элементу (i, j) для константной матрицы
    const T& operator()(size_t i, size_t j) const{
      if(i >= rows || j >= cols){
        throw out_of_range("Index out of range"); // вызов исключения
      }
      return data[i][j];
    }

    // Получение количества строк
    size_t get_rows() const{
      return rows;
    }
    // Получение количества столбцов
    size_t get_cols() const{
      return cols;
    }

    // Сравнение матриц
    bool operator==(const Matrix<T>& other) const{
      if(rows != other.rows || cols != other.cols){
        return false;
      }
      return data == other.data;
    }
    bool operator!=(const Matrix<T>& other) const{
      return !(*this == other);
    }

    //Арифметические операции
    //Умножение на скаляр
    Matrix<T> operator*(const T& x) const{
      Matrix<T> result(rows, cols);
      size_t i, j;
      for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
          result(i, j) = data[i][j] * x;
        }
      }
      return result;
    }
    // Сложение матриц
    Matrix<T> operator+(const Matrix<T>& other) const{
      if(rows != other.rows || cols != other.cols){
        throw invalid_argument("Matrices have different sizes");
      }
      Matrix<T> result(rows, cols);
      size_t i, j;
      for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
          result(i, j) = data[i][j] + other(i, j);
        }
      }
      return result;
    }
    // Вычитание матриц
    // конст в конце, т.к. мы не изменяем исходный объект
    Matrix<T> operator-(const Matrix<T>& other) const{
      if(rows != other.rows || cols != other.cols){
        throw invalid_argument("Matrices have different sizes");
      }
      Matrix<T> result(rows, cols);
      size_t i, j;
      for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
          result(i, j) = data[i][j] - other(i, j);
        }
      }
      return result;
    }
    // Отрицание матрицы
    Matrix<T> operator-() const{
      Matrix<T> result(rows, cols);
      size_t i, j;
      for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
          result(i, j) = -data[i][j];
        }
      }
      return result;
    }
    // Умножение матриц
    Matrix<T> operator*(const Matrix<T>& other) const{
      if(cols != other.rows){
        throw invalid_argument("Matrices have incompatible sizes");
      }
      Matrix<T> result(rows, other.cols);
      size_t i, j, k;
      for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
          T sum = 0; //любой тип может быть
          for(k = 0; k < cols; k++){
            sum += data[i][k] * other(k, j);
          }
          result(i, j) = sum;
        }
      }
      return result;
    }
    // Транспонирование матрицы
    Matrix<T> transpose() const{
      Matrix<T> result(cols, rows);
      size_t i, j;
      for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
          result(j, i) = data[i][j];
        }
      }
      return result;
    } 
    // Сложение строки к строке.
    Matrix<T> plus_row(size_t i, size_t k){
      if(i >= rows || k >= rows){
        throw out_of_range("Index out of range");
      }
      size_t j;
      for(j = 0; j < cols; j++){
        data[i][j] += data[k][j];
      }
    }
    // Сложение столбца к столбцу.
    Matrix<T> plus_col(size_t i, size_t k){
      if(i >= cols || k >= cols){
        throw out_of_range("Index out of range");
      }
      size_t j;
      for(j = 0; j < cols; j++){
        data[j][i] += data[j][k];
      }
    }

    // Получение подматрицы, начиная с позиции (row, col) и размерами (rows, cols)
    Matrix<T> submatrix(size_t row, size_t col, size_t rows, size_t cols) const{
      if(row + rows > this->rows || col + cols > this->cols){
        throw out_of_range("Index out of range");
      }
      Matrix<T> result(rows, cols);
      size_t i, j;
      for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
          result(i, j) = data[row + i][col + j];
        }
      }
      return result;
    }

    // Создание единичной матрицы
    static Matrix<T> eye(size_t n){
      Matrix<T> result(n, n);
      size_t i, j;
      for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
          if(i == j){
            result(i, j) = 1;
          } else{
            result(i, j) = 0;
          }
        }
      }
      return result;
    }
    // Создание нулевой матрицы
    static Matrix<T> zeros(size_t rows, size_t cols){
      Matrix<T> result(rows, cols);
      size_t i, j;
      for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            result(i, j) = 0;
        }
      }
      return result;
    }

    // Конкатенации
    // Конкатенация матриц по горизонтали
    Matrix<T> horizontal_concatenate(const Matrix<T>& other) const{
      if(rows != other.rows){
        throw invalid_argument("Matrices have different number of rows");
      }
      Matrix<T> result(rows, cols + other.cols);
      size_t i, j;
      for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
          result(i, j) = data[i][j];
        }
        for(j = 0; j < other.cols; j++){
          result(i, cols + j) = other.data[i][j];
        }
      }
      return result;
    }
    // Конкатенация матриц по вертикали
    Matrix<T> vertical_concatenate(const Matrix<T>& other) const{
      if(cols != other.cols){
        throw invalid_argument("Matrices have different number of columns");
      }
      Matrix<T> result(rows + other.rows, cols);
      size_t i, j;
      for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
          result(i, j) = data[i][j];
        }
      }
      for(i = 0; i < other.rows; i++){
        for(j = 0; j < cols; j++){
          result(rows + i, j) = other.data[i][j];
        }
      }
      return result;
    }

    // Поменять две строки местами
    void swap_rows(size_t i, size_t k){
      if(i >= rows || k >= rows){
        throw out_of_range("Index out of range");
      }
      swap(data[i], data[k]);
    }
    // Поменять два столбца местами
    void swap_cols(size_t i, size_t k){
      if(i >= cols || k >= cols){
        throw out_of_range("Index out of range");
      }
      size_t j;
      for(j = 0; j < rows; j++){
        swap(data[j][i], data[j][k]);
      }
    }

    //Чтение матрицы с консоли. Перед этим нужно создать пустую матрицу в main
    void read_from_console(){
      cout << "Enter the number of rows: ";
      cin >> rows;
      cout << "Enter the number of columns: ";
      cin >> cols;
      data.resize(rows, vector<T>(cols));

      cout<< "Enter the elements of the matrix: " << endl;
      size_t i, j;
      for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
          cin >> data[i][j];
        }
      }
    }
    // Определитель матрицы
    T determinant(){
      if(rows != cols){
        throw invalid_argument("Matrix is not square");
      }
      Matrix<T> temp(*this);
      T det = 1;
      size_t i, j, k;
      for(i = 0; i < rows; i++){
        size_t pivot_row = i;
        for(j = i + 1; j < rows; j++){
          if(abs(temp(j, i)) > abs(temp(pivot_row, i))){
            pivot_row = j;
          }
        }
        if(pivot_row != i){
          temp.swap_rows(i, pivot_row);
          det = -det;
        }

        T pivot = temp(i, i);
        if(pivot == 0){
          return 0;
        }

        det *= pivot;
        for(j = i + 1; j < rows; j++){
          T s = temp(j, i) / pivot;
          for(k = i; k < cols; k++){
            temp(j, k) -= s * temp(i, k);
          }
        }
      }
      return det;
    }
    Matrix<T> inverse2() const {
        if (rows != cols) {
            throw logic_error("Matrix is not square");
        }
        size_t i, j, k;
        Matrix<T> augmentedMatrix(rows, cols * 2);
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                augmentedMatrix(i, j) = data[i][j];
                augmentedMatrix(i, j + cols) = (i == j) ? 1 : 0;
            }
        }
        for (i = 0; i < rows; i++) {
            if (augmentedMatrix(i, i) == 0) {
                throw runtime_error("Matrix is singular");
            }
            for (j = 0; j < rows; j++) {
                if (i != j) {
                    T s = augmentedMatrix(j, i) / augmentedMatrix(i, i);
                    for (k = 0; k < 2 * cols; k++) {
                        augmentedMatrix(j, k) -= s * augmentedMatrix(i, k);
                    }
                }
            }
        }
        for( i = 0; i < rows; i++) {
            T s = augmentedMatrix(i, i);
            for (j = 0; j < 2 * cols; j++) {
                augmentedMatrix(i, j) /= s;
            }
        }
        Matrix<T> result(rows, cols);
        for (i = 0; i < rows; ++i) {
            for (j = 0; j < cols; ++j) {
                result(i, j) = augmentedMatrix(i, j + cols);
            }
        }

        return result;
    }
    T deter() {
        if(rows != cols){
          throw invalid_argument("Matrix is not square");
        }
        size_t i, j, z;
        double det = 1;
        for(i = 0; i < rows; i++) {
            T num = data[i][i];
            size_t k = i;
            for(j = i; j < rows; j++)
                if (fabs(data[j][i]) > num){
                    num = data[j][i];
                    k = j;
                }
            if (num == 0) {
                //Матрица вырождена
                det = 0;
                break;
            }
            if(i != k){
                swap_rows(i, k);
                det *= -1;
            }
            det *= data[i][i];
            for(j = 0; j < rows; j++){
                if(i != j){
                    double s = data[j][i] / data[i][i];
                    for(z = 0; z < rows; z++){
                        data[j][z] -= data[i][z] * s;
                    }
                }
            }

        }
        if(fabs(det) < 0.001){
            return 0;
        }
        return det;
    }

    T ranking(){
        size_t i, j, rank, z, d;
        rank = cols;
        int used[rows];
        //выбираем строки(помечаем). изначально 0, т.к. ничего не выбрали
        for (i = 0; i < rows; i++) used[i] = 0;
        for (i = 0; i < cols; i++) {
            //если среди невыбранных строк нет ненулывых, то пропускаем шаг и уменьшаем ранг
            for (j = 0; j < rows; j++) {
                if(!used[j] && fabs(data[j][i]) > 0.000000001){
                    break;
                }
            }
            if(j == rows)
                rank--;
                //если есть ненулевая строка в столбце, то помечаем строчку и просто делаем гаусса.
            else{
                used[j] = 1;
                for(z = i + 1; z < cols; z++)
                    data[j][z] /= data[j][i];
                for(z = 0; z < rows; z++){
                    if(z != j && fabs(data[z][i]) > 0.0000001){
                        for(d = i + 1; d < cols; d++){
                            data[z][d] -= data[j][d] * data[z][i];
                        }
                    }
                }
            }
        }
        return rank;
    }
    // «Красивый» вывод матрицы
    template <typename Tstream> 
    friend ostream &operator<<(ostream &out, const Matrix<Tstream>& m);


};

template <typename Tstream> 
ostream &operator<<(ostream &out, const Matrix<Tstream>& m) {
    const int MAX_NUM_DIGITS = 5;
    out << endl;
    for (int i = 0; i < m.rows; ++i) {
        for (int j = 0; j < m.cols; ++j) {
            out << setw(MAX_NUM_DIGITS) << m(i, j) << " ";
        }
        out << endl;
    }
    return out;
}