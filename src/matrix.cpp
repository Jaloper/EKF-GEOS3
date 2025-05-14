#include "..\include\matrix.hpp"

Matrix::Matrix(){
	this->n_row=0;
	this->n_column=0;
	this->data=nullptr;
}

Matrix::Matrix(const int n_row, const int n_column) {
    if (n_row <= 0 || n_column <= 0) {
		cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = n_row;
	this->n_column = n_column;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n_row; i++) {
		this->data[i] = (double *) malloc(n_column*sizeof(double));
	}
}

Matrix::Matrix(const int n) {
    if (n <= 0) {
		cout << "Matrix create: error in dimensions\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = 1;
	this->n_column = n;
	this->data = (double **) malloc(n*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	this->data[0] = (double *) malloc(n_column*sizeof(double));
}

double& Matrix::operator () (const int row, const int column) {
	if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
		cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[row - 1][column - 1];
}

double Matrix::operator()(int row, int column) const {
    if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
		cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[row - 1][column - 1];
}


double& Matrix::operator () (const int n) {
	if (n <= 0 || n > this->n_column) {
		cout << "Matrix get: error in dimensions \n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[(n - 1)/this->n_column][(n - 1)%this->n_column];
}

Matrix& Matrix::operator + (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sum: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + m(i,j);
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator - (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sub: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - m(i,j);
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator * (Matrix &m) {
    if (this->n_column != m.n_row) {
        cout << "Matrix multiplication: error in dimensions (columns of first matrix must match rows of second)\n";
        exit(EXIT_FAILURE);
    }
    
    Matrix *m_aux = new Matrix(this->n_row, m.n_column);
    
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= m.n_column; j++) {
            double sum = 0.0;
            for(int k = 1; k <= this->n_column; k++) {
                sum =sum + (*this)(i, k) * m(k, j);
            }
            (*m_aux)(i, j) = sum;
        }
    }
    
    return *m_aux;
}

Matrix& Matrix::operator * (const Matrix &m) const{
    if (this->n_column != m.n_row) {
        cout << "Matrix multiplication: error in dimensions (columns of first matrix must match rows of second)\n";
        exit(EXIT_FAILURE);
    }
    
    Matrix *m_aux = new Matrix(this->n_row, m.n_column);
    
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= m.n_column; j++) {
            double sum = 0.0;
            for(int k = 1; k <= this->n_column; k++) {
                sum =sum + (*this)(i, k) * m(k, j);
            }
            (*m_aux)(i, j) = sum;
        }
    }
    
    return *m_aux;
}

Matrix& Matrix::operator / (Matrix &m) {
    if (m.n_row != m.n_column) {
        cout << "Matrix division: matrix must be square\n";
        exit(EXIT_FAILURE);
    }
    
    if (this->n_column != m.n_row) {
        cout << "Matrix division: error in dimensions (columns of first matrix must match rows of second)\n";
        exit(EXIT_FAILURE);
    }
    
    Matrix inv_m = inv(m);
    Matrix *m_aux = new Matrix((*this) * inv_m);
    
    return *m_aux;
}

Matrix& Matrix::operator = (Matrix &m) {
    if (this == &m) return *this;
	
    for (int i = 0; i < this->n_row; i++) free(this->data[i]);
    free(this->data);

    this->n_row = m.n_row;
    this->n_column = m.n_column;
    this->data = (double **)malloc(m.n_row * sizeof(double *));
	
    if (this->data == NULL) {
        cout << "Matrix assignment: memory of rows could not be assigned\n";
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < m.n_row; i++) {
        this->data[i] = (double *)malloc(m.n_column * sizeof(double));
        if (this->data[i] == NULL) {
             cout << "Matrix assignment: memory of columns could not be assigned\n";
        exit(EXIT_FAILURE);
		}
        for (int j = 0; j < m.n_column; j++) {
            this->data[i][j] = m.data[i][j];
        }
    }

    return *this;
}

Matrix& Matrix::operator + (double d) {
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + d;
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator - (double d) {
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - d;
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator * (double d) {
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) * d;
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator / (double d) {
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) / d;
		}
	}
	
	return *m_aux;
}

ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}

Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

Matrix& inv(Matrix &m) {
    if (m.n_row != m.n_column) {
        cout << "Matrix inversion: only square matrices can be inverted\n";
        exit(EXIT_FAILURE);
    }
    
    int n = m.n_row;
    Matrix *result = new Matrix(n, n);
    Matrix m_aux(n, n);
    
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            m_aux(i, j) = m(i, j);
        }
    }
    
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            if(i == j) (*result)(i, j) = 1.0;
			else (*result)(i, j) = 0.0;
        }
    }
    
    for(int col = 1; col <= n; col++) {
        int max_row = col;
        for(int row = col + 1; row <= n; row++) {
            if(fabs(m_aux(row, col)) > fabs(m_aux(max_row, col))) {
                max_row = row;
            }
        }
        
        if(max_row != col) {
            for(int j = 1; j <= n; j++) {
                std::swap(m_aux(col, j), m_aux(max_row, j));
                std::swap((*result)(col, j), (*result)(max_row, j));
            }
        }
        
        if(fabs(m_aux(col, col)) < 1e-10) {
            cout << "Matrix inversion: matrix is singular (non-invertible)\n";
            exit(EXIT_FAILURE);
        }
        
        double pivot = m_aux(col, col);
        for(int j = 1; j <= n; j++) {
            m_aux(col, j) /= pivot;
            (*result)(col, j) /= pivot;
        }
        
        for(int row = 1; row <= n; row++) {
            if(row != col && fabs(m_aux(row, col)) > 1e-10) {
                double factor = m_aux(row, col);
                for(int j = 1; j <= n; j++) {
                    m_aux(row, j) -= factor * m_aux(col, j);
                    (*result)(row, j) -= factor * (*result)(col, j);
                }
            }
        }
    }
    
    return *result;
}

Matrix& eye(int n) {
   	Matrix *m_aux = new Matrix(n, n);
	
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n; j++) {
			if (i == j) (*m_aux)(i,j) = 1;
			else (*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

Matrix& transpose(Matrix &m) {
    Matrix *m_aux = new Matrix(m.n_column, m.n_row);
		for(int i = 1; i <= m.n_row; i++) {
			for(int j = 1; j <= m.n_column; j++) {
				 (*m_aux)(j, i) = m(i, j);
			}
		}
	return (*m_aux);
}

Matrix& zeros(const int n) {
	Matrix *m_aux = new Matrix(n);
	
	for(int i = 1; i <= n; i++) {
		(*m_aux)(i) = 0;
	}
	
	return (*m_aux);
}

double norm(Matrix &m){
    double sum = 0.0;
    for(int i = 1; i <= m.n_row; i++) {
        for(int j = 1; j <= m.n_column; j++) {
            sum += m(i,j) * m(i,j);
        }
    }
    return sqrt(sum);
}

double dot(Matrix& v1, Matrix& v2) {
    if (v1.n_row != 1 || v2.n_row != 1 || v1.n_column != v2.n_column) {
         cout << "Vector dot: dimensions of vectors  are incompatible\n";
		 exit(EXIT_FAILURE);
    }
    double result = 0.0;
		for(int j = 1; j <= v1.n_column; j++) {
			 result += v1(1, j)*v2(1, j);
		}
    return result;
}

Matrix& cross(Matrix& v1, Matrix& v2) {
    if (v1.n_row != 1 || v2.n_row != 1 || v1.n_column != 3 || v2.n_column != 3) {
        cout << "Vector cross: both vectors must have 3 components\n";
        exit(EXIT_FAILURE);
    }

     Matrix *m_aux = new Matrix(1,3);
    (*m_aux)(1, 1) = v1(1, 2)*v2(1, 3) - v1(1, 3)*v2(1, 2);
    (*m_aux)(1, 2) = v1(1, 3)*v2(1, 1) - v1(1, 1)*v2(1, 3);
    (*m_aux)(1, 3) = v1(1, 1)*v2(1, 2) - v1(1, 2)*v2(1, 1); 

    return *m_aux;
}

Matrix& extract_vector(Matrix &v, int init, int fin) {
    if (v.n_row != 1 || init < 1 || fin > v.n_column || init > fin) {
        cout << "Matrix extract_vector: invalid vector extraction\n";
        exit(EXIT_FAILURE);
    }

    int subvector_size = fin - init + 1;
    Matrix* subvector = new Matrix(1,subvector_size);

    for (int j = 1; j <= subvector_size; j++) {
        (*subvector)(1, j) = v( 1,init + j - 1);
    }

    return *subvector;
}

Matrix& union_vector(Matrix &v1, Matrix &v2) {
    if (v1.n_row != 1 || v2.n_row != 1) {
        cerr << "Matrix union_vector: invalid vector union\n";
        exit(EXIT_FAILURE);
    }
    
    Matrix* union_vec = new Matrix(1, v1.n_column + v2.n_column);
    
    for (int j = 1; j <= v1.n_column; ++j) {
        (*union_vec)(1, j) = v1(1, j);
    }
    
    for (int j = 1; j <= v2.n_column; ++j) {
        (*union_vec)(1, v1.n_column + j) = v2(1, j);
    }
    
    return *union_vec;
}

Matrix& extract_row(Matrix &m, int n){
	if (m.n_row < n || n < 1) {
		cout << "Matrix extract_row: row token is out of bounds \n";
		 exit(EXIT_FAILURE);
    }
	
	Matrix *m_aux = new Matrix(1, m.n_column);
	for(int j = 1; j <= m.n_column; j++) {
		(*m_aux)(1, j) = m(n, j);
	}
	
	 return *m_aux;
}

Matrix& extract_column(Matrix &m, int n){
	if (m.n_column < n || n < 1) {
		cout << "Matrix extract_column: column token is out of bounds \n";
		 exit(EXIT_FAILURE);
    }
	
	Matrix *m_aux = new Matrix(m.n_row, 1);
	for(int i = 1; i <= m.n_row; i++) {
		(*m_aux)(i, 1) = m(i, n);
	}
	
	 return *m_aux;
}

Matrix& assign_row(Matrix &m, Matrix &row, int n){
    if (m.n_row < n || n < 1 || row.n_row != 1) {
        cout << "Matrix assign_row: invalid row assignment\n";
        exit(EXIT_FAILURE);
    }
	
    Matrix *m_aux = new Matrix(m.n_row,max(row.n_column, m.n_column));
    *m_aux = m;

    for (int j = 1; j <= m_aux->n_column; j++) {
        (*m_aux)(n, j) = 0.0;
    }
    for (int j = 1; j <= row.n_column; j++) {
        (*m_aux)(n, j) = row(1, j);
    }

    return *m_aux;
}

Matrix &assign_column(Matrix &m, Matrix &col, int n) {
    if (m.n_column < n || n < 1 || col.n_column != 1) {
        cerr << "Matrix assign_column: invalid column assignment\n";
        exit(EXIT_FAILURE);
    }
	int n_columns = max(col.n_row, m.n_row);
    Matrix* m_aux = new Matrix(n_columns, m.n_column);
    
    for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++) {
            (*m_aux)(i, j) = m(i, j);
        }
    }

    for (int i = 1; i <= col.n_row; i++) {
        (*m_aux)(i, n) = col(i, 1);
    }

    for (int i = col.n_row + 1; i <= n_columns; i++) {
        (*m_aux)(i, n) = 0.0;
    }

    return *m_aux; 
}