// $Header$
//--------------------------------------------------------------------------------
// Matrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/02
//
/** @file matrix.hpp
 *  @brief This header file declares the Matrix class for matrix operations.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class Matrix {
public:
    int n_row, n_column;
    double **data;

    // Parameterized constructor
    /** @brief Default constructor that initializes an empty matrix with zero rows and columns and null data pointer. */
    Matrix();

    /** @brief Constructor that creates a matrix with specified rows and columns.
     *  @param n_row Number of rows in the matrix.
     *  @param n_column Number of columns in the matrix.
     */
    Matrix(const int n_row, const int n_column);

    /** @brief Constructor that creates a row vector (1xN matrix).
     *  @param n Number of elements in the row vector.
     */
    Matrix(const int n);
    
    // Member operators
    /** @brief Accesses and modifies the element at the specified row and column (1-based indexing).
     *  @param row Row index (1-based).
     *  @param column Column index (1-based).
     *  @return Reference to the element at (row, column).
     */
    double& operator () (const int row, const int column);

    /** @brief Accesses the element at the specified row and column (1-based indexing) for const objects.
     *  @param row Row index (1-based).
     *  @param column Column index (1-based).
     *  @return Value of the element at (row, column).
     */
    double operator()(int row, int column) const;  

    /** @brief Accesses and modifies the element at the specified linear index (1-based) for vector-like matrices.
     *  @param n Linear index (1-based).
     *  @return Reference to the element at the specified index.
     */
    double& operator () (const int n);

    /** @brief Adds two matrices element-wise, returning a new matrix.
     *  @param m Matrix to add.
     *  @return New matrix containing the element-wise sum.
     */
    Matrix& operator + (Matrix &m);

    /** @brief Adds two matrices element-wise for const objects, returning a new matrix.
     *  @param m Matrix to add.
     *  @return New matrix containing the element-wise sum.
     */
    Matrix& operator + (const Matrix &m) const;

    /** @brief Subtracts one matrix from another element-wise, returning a new matrix.
     *  @param m Matrix to subtract.
     *  @return New matrix containing the element-wise difference.
     */
    Matrix& operator - (Matrix &m);

    /** @brief Multiplies two matrices, returning a new matrix.
     *  @param m Matrix to multiply with.
     *  @return New matrix containing the matrix product.
     */
    Matrix& operator * (Matrix &m);

    /** @brief Multiplies two matrices for const objects, returning a new matrix.
     *  @param m Matrix to multiply with.
     *  @return New matrix containing the matrix product.
     */
    Matrix& operator * (const Matrix &m) const;

    /** @brief Divides the current matrix by another matrix (using inverse), returning a new matrix.
     *  @param m Matrix to divide by (must be square and invertible).
     *  @return New matrix containing the result of division.
     */
    Matrix& operator / (Matrix &m);

    /** @brief Assigns the contents of another matrix to the current matrix.
     *  @param m Matrix to assign from.
     *  @return Reference to the current matrix after assignment.
     */
    Matrix& operator = (Matrix &m);

    /** @brief Assigns the contents of another matrix to the current matrix for const objects.
     *  @param m Matrix to assign from.
     *  @return Reference to the current matrix after assignment.
     */
    Matrix& operator = (const Matrix &m);

    /** @brief Adds a scalar to each element of the matrix, returning a new matrix.
     *  @param d Scalar value to add.
     *  @return New matrix with the scalar added to each element.
     */
    Matrix& operator + (double d);

    /** @brief Subtracts a scalar from each element of the matrix, returning a new matrix.
     *  @param d Scalar value to subtract.
     *  @return New matrix with the scalar subtracted from each element.
     */
    Matrix& operator - (double d);

    /** @brief Multiplies each element of the matrix by a scalar, returning a new matrix.
     *  @param d Scalar value to multiply by.
     *  @return New matrix with each element multiplied by the scalar.
     */
    Matrix& operator * (double d);

    /** @brief Divides each element of the matrix by a scalar, returning a new matrix.
     *  @param d Scalar value to divide by.
     *  @return New matrix with each element divided by the scalar.
     */
    Matrix& operator / (double d);

    /** @brief Destructor that frees the dynamically allocated memory for the matrix data. */
    ~Matrix();

    /** @brief Copy constructor that creates a deep copy of another matrix.
     *  @param other Matrix to copy from.
     */
    Matrix(const Matrix& other);
    
    // Non-member operators
    /** @brief Outputs the matrix to an output stream in a formatted manner.
     *  @param o Output stream to write to.
     *  @param m Matrix to output.
     *  @return Reference to the output stream.
     */
    friend ostream& operator << (ostream &o, Matrix &m);
};

// Methods
/** @brief Creates a matrix of zeros with specified dimensions.
 *  @param n_row Number of rows.
 *  @param n_column Number of columns.
 *  @return New matrix filled with zeros.
 */
Matrix& zeros(const int n_row, const int n_column);

/** @brief Computes the inverse of a square matrix using Gaussian elimination.
 *  @param m Matrix to invert (must be square).
 *  @return New matrix containing the inverse.
 */
Matrix& inv(Matrix &m);

/** @brief Creates an identity matrix of size NxN.
 *  @param n Size of the square matrix.
 *  @return New identity matrix.
 */
Matrix& eye(int n);

/** @brief Computes the transpose of a matrix.
 *  @param m Matrix to transpose.
 *  @return New matrix containing the transpose.
 */
Matrix& transpose(Matrix &m);

/** @brief Creates a row vector of zeros with specified length.
 *  @param n Number of elements in the row vector.
 *  @return New row vector filled with zeros.
 */
Matrix& zeros(const int n);

/** @brief Computes the Euclidean norm (Frobenius norm) of a matrix.
 *  @param m Matrix to compute the norm for.
 *  @return Double value representing the norm.
 */
double norm(Matrix &m);

/** @brief Computes the dot product of two row vectors.
 *  @param v1 First vector (1xN).
 *  @param v2 Second vector (1xN).
 *  @return Double value representing the dot product.
 */
double dot(Matrix& v1, Matrix& v2);

/** @brief Computes the cross product of two 3D row vectors.
 *  @param v1 First vector (1x3).
 *  @param v2 Second vector (1x3).
 *  @return New row vector (1x3) containing the cross product.
 */
Matrix& cross(Matrix &v1, Matrix &v2);

/** @brief Extracts a subvector from a row vector.
 *  @param v Row vector to extract from.
 *  @param init Starting index (1-based).
 *  @param fin Ending index (1-based).
 *  @return New row vector containing the extracted elements.
 */
Matrix& extract_vector(Matrix &v, int init, int fin);

/** @brief Concatenates two row vectors.
 *  @param v1 First row vector.
 *  @param v2 Second row vector.
 *  @return New row vector containing the concatenated elements.
 */
Matrix& union_vector(Matrix &v1, Matrix &v2);

/** @brief Extracts a row from a matrix.
 *  @param m Matrix to extract from.
 *  @param n Row index to extract (1-based).
 *  @return New row vector containing the extracted row.
 */
Matrix& extract_row(Matrix &m, int n);

/** @brief Extracts a column from a matrix.
 *  @param m Matrix to extract from.
 *  @param n Column index to extract (1-based).
 *  @return New column vector containing the extracted column.
 */
Matrix& extract_column(Matrix &m, int n);

/** @brief Assigns a row vector to a specific row in a matrix.
 *  @param m Matrix to modify.
 *  @param row Row vector to assign.
 *  @param n Row index to assign to (1-based).
 *  @return New matrix with the specified row replaced.
 */
Matrix& assign_row(Matrix &m, Matrix &row, int n);

/** @brief Assigns a column vector to a specific column in a matrix.
 *  @param m Matrix to modify.
 *  @param row Column vector to assign.
 *  @param n Column index to assign to (1-based).
 *  @return New matrix with the specified column replaced.
 */
Matrix& assign_column(Matrix &m, Matrix &row, int n);

#endif