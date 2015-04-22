#ifndef EIGENMAT_H
#define EIGENMAT_H

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "config.h"

typedef Eigen::SparseMatrix<VtxType> VtxMat;
typedef Eigen::SparseMatrix<WgtType> WgtMat;
typedef Eigen::MatrixXd DenseMat;
typedef Eigen::VectorXd DenseVec;
typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;

#endif