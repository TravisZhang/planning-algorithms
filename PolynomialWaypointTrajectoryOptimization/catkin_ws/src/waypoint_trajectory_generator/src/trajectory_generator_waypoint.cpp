#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

// ooqp related headers
#include <ooqp/QpGenData.h>
#include <ooqp/QpGenVars.h>
#include <ooqp/QpGenResiduals.h>
#include <ooqp/GondzioSolver.h>
#include <ooqp/QpGenSparseMa27.h>

// ooqpei related headers
// #include "ooqp_eigen_interface/OoqpEigenInterface.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>

// osqp related headers
#include "osqp.h"

using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for(int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}
/*

    STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGeneration function

    variable declaration: input       const int d_order,                    // the order of derivative
                                      const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
                                      const Eigen::MatrixXd &Vel,           // boundary velocity
                                      const Eigen::MatrixXd &Acc,           // boundary acceleration
                                      const Eigen::VectorXd &Time)          // time allocation in each segment
                          output      MatrixXd PolyCoeff(m, 3 * p_num1d);   // position(x,y,z), so we need (3 * p_num1d) coefficients

*/

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment
    int d_num = d_order;  // dimension of one point for close form sol. (>=3)

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);
    ROS_WARN("[TG] d_order: %d", d_order);
    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */
    for (int idx = 0; idx < Path.cols(); ++idx) {
        ROS_WARN("[TG] calculating %d dimension trajectory", idx);
        VectorXd start_pose = VectorXd::Zero(d_num);
        start_pose(0) = Path(0, idx);
        start_pose(1) = Vel(0, idx);
        start_pose(2) = Acc(0, idx);
        VectorXd end_pose = VectorXd::Zero(d_num);
        end_pose(0) = Path(m, idx);
        end_pose(1) = Vel(1, idx);
        end_pose(2) = Acc(1, idx);
        int n_all_poly = m * p_num1d;
        
    /*   Produce the derivatives in X, Y and Z axis directly.  */
        MatrixXd M = MatrixXd::Zero(2*(d_num)*m, p_num1d*m);
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < d_order; ++k) {
                M(j*p_num1d+k, j*p_num1d+k) = Factorial(k);
                for (int i = k; i <= p_order; ++i) {
                    M(j*p_num1d+d_num+k, j*p_num1d+i) = Factorial(i) / Factorial(i-k) * pow(Time(j), i-k);
                }
            }
        }
        ROS_INFO("[TG] cal M done, size: %ld, %ld", M.rows(), M.cols());
        // std::cout << M << std::endl << std::endl;
        MatrixXd Ct = MatrixXd::Zero(2*d_num*m, d_num*(m+1));
        for (int k = 0; k < d_order; ++k) {
            Ct(k, k) = 1;
            Ct(k+d_num+2*d_num*(m-1), k+d_num+m-1) = 1;
        }
        for (int j = 0; j < m-1; ++j) {
            Ct(d_num+2*d_num*j,d_num+j) = 1;
            Ct(2*d_num*(j+1), d_num+j) = 1;
            for (int i = 0; i < d_order-1; ++i) {
                Ct(d_num+1+2*d_num*j+i, 2*d_num+m-1+j*(d_order-1)+i) = 1;
                Ct(1+2*d_num*(j+1)+i, 2*d_num+m-1+j*(d_order-1)+i) = 1;
            }
        }
        ROS_INFO("[TG] cal Ct done, size: %ld, %ld", Ct.rows(), Ct.cols());
        // std::cout << Ct << std::endl << std::endl;
    /*   Produce the Minimum Snap cost function, the Hessian Matrix   */
        MatrixXd Q = MatrixXd::Zero(m*p_num1d, m*p_num1d);
        for (int j = 0; j < m; ++j) {
            for (int k = 4; k <= p_order; ++k) {
                for (int l = 4; l <= p_order; ++l) {
                    Q(k+j*p_num1d, l+j*p_num1d) = k*(k-1)*(k-2)*(k-3)*l*(l-1)*(l-2)*(l-3)*pow(Time(j), k+l-7)/(k+l-7);
                }
            }
        }
        ROS_INFO("[TG] cal Q done, size: %ld, %ld", Q.rows(), Q.cols());
        // std::cout << Q << std::endl << std::endl;
        // compute R
        MatrixXd R = Ct.transpose() * (M.inverse()).transpose() * Q * M.inverse() * Ct;
        ROS_INFO("[TG] cal R done, size: %ld, %ld", R.rows(), R.cols());
        // std::cout << R << std::endl << std::endl;
        // compute df
        VectorXd df = VectorXd::Zero(2*d_num+m-1);
        df.segment(0,d_num) = start_pose;
        df.segment(d_num+m-1,d_num) = end_pose;
        for (int j = 0; j < m-1; ++j) {
            df(j+d_num) = Path(j+1, idx);
        }
        ROS_INFO("[TG] cal df done, size: %ld, %ld", df.rows(), df.cols());
        // std::cout << df << std::endl << std::endl;
        MatrixXd R_pp = R.block(2*d_num+m-1, 2*d_num+m-1, (d_order-1)*(m-1), (d_order-1)*(m-1));
        ROS_INFO("[TG] cal R_pp done, size: %ld, %ld", R_pp.rows(), R_pp.cols());
        MatrixXd R_fp = R.block(0, 2*d_num+m-1, 2*d_num+m-1, (d_order-1)*(m-1));
        ROS_INFO("[TG] cal R_fp done, size: %ld, %ld", R_fp.rows(), R_fp.cols());
        VectorXd dp_best = -R_pp.inverse() * R_fp.transpose() * df;
        ROS_INFO("[TG] cal dp_best done, size: %ld, %ld", dp_best.rows(), dp_best.cols());
        VectorXd dfp;
        dfp.resize(d_num * (m+1));
        dfp << df,
               dp_best;
        ROS_INFO("[TG] cal dfp done, size: %ld, %ld", dfp.rows(), dfp.cols());
        VectorXd coeff = M.inverse() * Ct * dfp;
        ROS_INFO("[TG] cal coeff done, size: %ld, %ld", coeff.rows(), coeff.cols());
        // std::cout << coeff << std::endl << std::endl;
        std::cout << "printing coeffs " << std::endl;
        for (int j = 0; j < m; ++j) {
            PolyCoeff.block(j, idx*p_num1d, 1, p_num1d) = coeff.segment(j*p_num1d, p_num1d).transpose();
            std::cout << "seg [" << j << "]" << std::endl;
            std::cout << PolyCoeff.block(j, idx*p_num1d, 1, p_num1d) << std::endl << std::endl;
        }
    }
    LogData(PolyCoeff, "closed_form_result");
    return PolyCoeff;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::SolvebyOOQP(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time) {
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment
    int d_num = d_order;  // dimension of one point for close form sol. (>=3)

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    ROS_WARN("[TG] d_order: %d", d_order);
    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */
    for (int idx = 0; idx < Path.cols(); ++idx) {   

        int nnzA = GetnnzA(m, d_order);
        ROS_INFO("[TG][OOQP] nnzA: %d", nnzA);
        int *irowA = new (std::nothrow) int[nnzA];
        int *jcolA = new (std::nothrow) int[nnzA];
        double *dA = new (std::nothrow) double[nnzA];
        if (!irowA || !jcolA || !dA) {
            // Handle error
            ROS_ERROR("ERROR: memory allocation failed !!!!");
            return PolyCoeff;
        }
        int k_A = 0;
        int row_idx_A = 0;

        ROS_WARN("[TG] calculating %d dimension trajectory", idx);
        VectorXd start_pose = VectorXd::Zero(d_num);
        start_pose(0) = Path(0, idx);
        start_pose(1) = Vel(0, idx);
        start_pose(2) = Acc(0, idx);
        VectorXd end_pose = VectorXd::Zero(d_num);
        end_pose(0) = Path(m, idx);
        end_pose(1) = Vel(1, idx);
        end_pose(2) = Acc(1, idx);
        int n_all_poly = m * p_num1d;
        
        // calculate Aeq_start
        MatrixXd Aeq_start = MatrixXd::Zero(d_num, n_all_poly);
        VectorXd beq_start = VectorXd::Zero(d_num);
        GenStartConstraint(m, d_order, Time, 0, Aeq_start, &row_idx_A, irowA, jcolA, dA, &k_A);
        beq_start = start_pose;
        ROS_INFO("[TG] cal Aeq_start done, size: %ld, %ld", Aeq_start.rows(), Aeq_start.cols());
        // calculate Aeq_end
        MatrixXd Aeq_end = MatrixXd::Zero(d_num, n_all_poly);
        VectorXd beq_end = VectorXd::Zero(d_num);
        GenEndConstraint(m, d_order, Time, 0, Aeq_end, &row_idx_A, irowA, jcolA, dA, &k_A);
        beq_end = end_pose;
        ROS_INFO("[TG] cal Aeq_end done, size: %ld, %ld", Aeq_end.rows(), Aeq_end.cols());
        // calculate Aeq_wp
        MatrixXd Aeq_wp = MatrixXd::Zero(m-1, n_all_poly);
        VectorXd beq_wp = VectorXd::Zero(m-1);
        GenWPConstraint(m, d_order, Time, 0, Path.block(0, idx, Path.rows(), 1), Aeq_wp, beq_wp, &row_idx_A, irowA, jcolA, dA, &k_A);
        ROS_INFO("[TG] cal Aeq_wp done, size: %ld, %ld", Aeq_wp.rows(), Aeq_wp.cols());
        // calculate pos continuity
        MatrixXd Aeq_con_p = MatrixXd::Zero(m-1, n_all_poly);
        GenContinuityConstraint(m, d_order, Time, 0, Aeq_con_p, &row_idx_A, irowA, jcolA, dA, &k_A);
        ROS_INFO("[TG] cal Aeq_con_p done, size: %ld, %ld", Aeq_con_p.rows(), Aeq_con_p.cols());
        // calculate vel continuity
        MatrixXd Aeq_con_v = MatrixXd::Zero(m-1, n_all_poly);
        GenContinuityConstraint(m, d_order, Time, 1, Aeq_con_v, &row_idx_A, irowA, jcolA, dA, &k_A);
        ROS_INFO("[TG] cal Aeq_con_v done, size: %ld, %ld", Aeq_con_v.rows(), Aeq_con_v.cols());
        // calculate acc continuity
        MatrixXd Aeq_con_a = MatrixXd::Zero(m-1, n_all_poly);
        GenContinuityConstraint(m, d_order, Time, 2, Aeq_con_a, &row_idx_A, irowA, jcolA, dA, &k_A);
        ROS_INFO("[TG] cal Aeq_con_a done, size: %ld, %ld", Aeq_con_a.rows(), Aeq_con_a.cols());
        // calculate jerk continuity
        MatrixXd Aeq_con_j = MatrixXd::Zero(m-1, n_all_poly);
        GenContinuityConstraint(m, d_order, Time, 3, Aeq_con_j, &row_idx_A, irowA, jcolA, dA, &k_A);
        ROS_INFO("[TG] cal Aeq_con_j done, size: %ld, %ld", Aeq_con_j.rows(), Aeq_con_j.cols());
        MatrixXd Aeq;
        VectorXd beq;
        int Aeq_row_num = 2*(d_num) + 5*(m-1);
        Aeq.resize(Aeq_row_num, n_all_poly);
        beq.resize(Aeq_row_num);
        Aeq << Aeq_start,
               Aeq_end,
               Aeq_wp,
               Aeq_con_p,
               Aeq_con_v,
               Aeq_con_a,
               Aeq_con_j;
        beq << beq_start,
               beq_end,
               beq_wp,
               VectorXd::Zero(4*(m-1));
        ROS_WARN("[TG] concatenated Aeq & beq");
        ROS_INFO("[TG][OOQP] row_idx_A: %d", row_idx_A);

        MatrixXd Q = MatrixXd::Zero(m*p_num1d, m*p_num1d);
        for (int j = 0; j < m; ++j) {
            for (int k = 4; k <= p_order; ++k) {
                for (int l = 4; l <= p_order; ++l) {
                    Q(k+j*p_num1d, l+j*p_num1d) = k*(k-1)*(k-2)*(k-3)*l*(l-1)*(l-2)*(l-3)*pow(Time(j), k+l-7)/(k+l-7);
                }
            }
        }
        ROS_INFO("[TG] cal Q done, size: %ld, %ld", Q.rows(), Q.cols());
        // std::cout << Q << std::endl << std::endl;

        int my = beq.rows();
        double *b = new (std::nothrow) double[my];
        if (!b) {
            // Handle error
            ROS_ERROR("ERROR: memory allocation failed !!!!");
            return PolyCoeff;
        }
        for (int i = 0; i < my; ++i) {
            b[i] = beq(i);
        }

        delete [] b;
        delete [] irowA;
        delete [] jcolA;
        delete [] dA;

    }
    
    return PolyCoeff;

}

void TrajectoryGeneratorWaypoint::GenStartConstraint(
        const int n_seq, 
        const int d_order, 
        const Eigen::VectorXd &Time, 
        const int order_deri, 
        Eigen::MatrixXd& Aeq_start,
        int *p_row_idx_A,
        int *irowA,
        int *jcolA,
        double *dA,
        int *p_k_A) {
    int p_order    = 2 * d_order - 1;              // the order of polynomial
    int p_num1d    = p_order + 1;
    int n_all_poly = n_seq * p_num1d;
    int d_num      = d_order;  // dimension of one point for close form sol. (>=3)
    Aeq_start = MatrixXd::Zero(d_num, n_all_poly);        
    for (int k = 0; k < d_num; ++k) {
        Aeq_start(k,k) = Factorial(k);
        AddtodA(Aeq_start(k,k), k, k, *p_row_idx_A, irowA, jcolA, dA, p_k_A);
    }
    *p_row_idx_A += d_order;
    ROS_INFO("[TG][OOQP] row_idx_A: %d, k_A: %d", *p_row_idx_A, *p_k_A);
}

void TrajectoryGeneratorWaypoint::GenEndConstraint(
        const int n_seq, 
        const int d_order, 
        const Eigen::VectorXd &Time, 
        const int order_deri, 
        Eigen::MatrixXd& Aeq_end,
        int *p_row_idx_A,
        int *irowA,
        int *jcolA,
        double *dA,
        int *p_k_A) {
    int p_order    = 2 * d_order - 1;              // the order of polynomial
    int p_num1d    = p_order + 1;
    int n_all_poly = n_seq * p_num1d;
    int d_num      = d_order;  // dimension of one point for close form sol. (>=3)
    Aeq_end = MatrixXd::Zero(d_num, n_all_poly);
    for (int k = 0; k < d_num; ++k) {
        for (int i = k; i <= p_order; ++i) {
            Aeq_end(k, (n_seq-1)*p_num1d+i) = Factorial(i) / Factorial(i-k) * pow(Time(n_seq-1), i-k);
            AddtodA(Aeq_end(k, (n_seq-1)*p_num1d+i), k, (n_seq-1)*p_num1d+i, *p_row_idx_A, irowA, jcolA, dA, p_k_A);
        }
    }
    *p_row_idx_A += d_order;
    ROS_INFO("[TG][OOQP] row_idx_A: %d, k_A: %d", *p_row_idx_A, *p_k_A);
}

void TrajectoryGeneratorWaypoint::GenWPConstraint(
        const int n_seq, 
        const int d_order, 
        const Eigen::VectorXd &Time, 
        const int order_deri, 
        const Eigen::VectorXd& Path,
        Eigen::MatrixXd& Aeq_wp,
        Eigen::VectorXd& beq_wp,
        int *p_row_idx_A,
        int *irowA,
        int *jcolA,
        double *dA,
        int *p_k_A) {
    int p_order    = 2 * d_order - 1;              // the order of polynomial
    int p_num1d    = p_order + 1;
    int n_all_poly = n_seq * p_num1d;
    int d_num      = d_order;  // dimension of one point for close form sol. (>=3)
    Aeq_wp = MatrixXd::Zero(n_seq-1, n_all_poly);
    beq_wp = VectorXd::Zero(n_seq-1);
    for (int j = 0; j < n_seq-1; ++j) {
        Aeq_wp(j, (j+1)*p_num1d) = 1;
        beq_wp(j) = Path(j+1);
        AddtodA(Aeq_wp(j, (j+1)*p_num1d), j, (j+1)*p_num1d, *p_row_idx_A, irowA, jcolA, dA, p_k_A);
    }
    *p_row_idx_A += n_seq-1;
    ROS_INFO("[TG][OOQP] row_idx_A: %d, k_A: %d", *p_row_idx_A, *p_k_A);
} 

void TrajectoryGeneratorWaypoint::GenContinuityConstraint(
        const int n_seq, 
        const int d_order, 
        const Eigen::VectorXd &Time, 
        const int order_deri, 
        Eigen::MatrixXd& Aeq_con,
        int *p_row_idx_A,
        int *irowA,
        int *jcolA,
        double *dA,
        int *p_k_A) {
    int p_order    = 2 * d_order - 1;              // the order of polynomial
    int p_num1d    = p_order + 1;
    int n_all_poly = n_seq * p_num1d;
    int k          = order_deri;
    Aeq_con = MatrixXd::Zero(n_seq-1, n_all_poly);
    for (int j = 0; j < n_seq-1; j++) {
        for (int i = k; i <= p_order; ++i) {
            Aeq_con(j, j*p_num1d+i) = Factorial(i) / Factorial(i-k) * pow(Time(j), i-k);
            AddtodA(Aeq_con(j, j*p_num1d+i), j, j*p_num1d+i, *p_row_idx_A, irowA, jcolA, dA, p_k_A);
        }
        Aeq_con(j,(j+1)*p_num1d+k) = - Factorial(k);
        AddtodA(Aeq_con(j,(j+1)*p_num1d+k), j, (j+1)*p_num1d+k, *p_row_idx_A, irowA, jcolA, dA, p_k_A);
    }
    *p_row_idx_A += n_seq-1;
    ROS_INFO("[TG][OOQP] row_idx_A: %d, k_A: %d", *p_row_idx_A, *p_k_A);
}

void TrajectoryGeneratorWaypoint::AddtodA(
        const double value, 
        const int irow,
        const int jcol,
        const int row_idx_A,
        int *irowA,
        int *jcolA,
        double *dA,
        int *p_k_A) {
    irowA[*p_k_A] = irow + row_idx_A;
    jcolA[*p_k_A] = jcol;
    dA[*p_k_A] = value;
    *p_k_A += 1;
}

int TrajectoryGeneratorWaypoint::GetnnzA(const int n_seq, const int d_order) {
    int p_order    = 2 * d_order - 1;              // the order of polynomial
    int p_num1d    = p_order + 1;
    int n_all_poly = n_seq * p_num1d;
    int d_num      = d_order;  // dimension of one point for close form sol. (>=3)
    int nnzA = 0;
    // start
    for (int k = 0; k < d_num; ++k) {
        nnzA += 1;
    }
    // end
    for (int k = 0; k < d_num; ++k) {
        for (int i = k; i <= p_order; ++i) {
            nnzA += 1;
        }
    }
    // wp
    for (int j = 0; j < n_seq-1; ++j) {
        nnzA += 1;
    }
    // continuity
    for (int k = 0; k < 4; ++k) {
        for (int j = 0; j < n_seq-1; j++) {
            for (int i = k; i <= p_order; ++i) {
                nnzA += 1;
            }
            nnzA += 1;
        }
    }
    return nnzA;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::SolvebyOOQPwithEigen(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // we are using eigen's sparse matrix 
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor> EigenSparse;
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator EigenSparseIter;
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment
    int d_num = d_order;  // dimension of one point for close form sol. (>=3)

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);
    ROS_WARN("[TG] d_order: %d", d_order);
    LogData(Path, "Path");
    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */
    for (int idx = 0; idx < Path.cols(); ++idx) {
        ROS_WARN("[TG] calculating %d dimension trajectory", idx);
        VectorXd start_pose = VectorXd::Zero(d_num);
        start_pose(0) = Path(0, idx);
        start_pose(1) = Vel(0, idx);
        start_pose(2) = Acc(0, idx);
        VectorXd end_pose = VectorXd::Zero(d_num);
        end_pose(0) = Path(m, idx);
        end_pose(1) = Vel(1, idx);
        end_pose(2) = Acc(1, idx);
        int n_all_poly = m * p_num1d;
        ROS_INFO("[TG] n_all_poly: %d", n_all_poly);
        
        // calculate Aeq_start
        EigenSparse Aeq_start(d_num, n_all_poly);
        VectorXd beq_start = VectorXd::Zero(d_num);
        for (int k = 0; k < d_order; ++k) {
            Aeq_start.insert(k,k) = Factorial(k);
        }
        beq_start = start_pose;
        ROS_INFO("[TG] cal Aeq_start done, size: %ld, %ld", Aeq_start.rows(), Aeq_start.cols());
        // calculate Aeq_end
        EigenSparse Aeq_end(d_num, n_all_poly);
        VectorXd beq_end = VectorXd::Zero(d_num);
        for (int k = 0; k < d_order; ++k) {
            for (int i = k; i <= p_order; ++i) {
                Aeq_end.insert(k, (m-1)*p_num1d+i) = Factorial(i) / Factorial(i-k) * pow(Time(m-1), i-k);
            }
        }
        beq_end = end_pose;
        ROS_INFO("[TG] cal Aeq_end done, size: %ld, %ld", Aeq_end.rows(), Aeq_end.cols());
        // calculate Aeq_wp
        EigenSparse Aeq_wp(m-1, n_all_poly);
        VectorXd beq_wp = VectorXd::Zero(m-1);
        for (int j = 0; j < m-1; ++j) {
            Aeq_wp.insert(j, (j+1)*p_num1d) = 1;
            beq_wp(j) = Path(j+1, idx);
        }
        ROS_INFO("[TG] cal Aeq_wp done, size: %ld, %ld", Aeq_wp.rows(), Aeq_wp.cols());
        // calculate pos continuity
        EigenSparse Aeq_con_p(m-1, n_all_poly);
        GenContinuityConstraint(m, d_order, Time, 0, Aeq_con_p);
        ROS_INFO("[TG] cal Aeq_con_p done, size: %ld, %ld", Aeq_con_p.rows(), Aeq_con_p.cols());
        // calculate vel continuity
        EigenSparse Aeq_con_v(m-1, n_all_poly);
        GenContinuityConstraint(m, d_order, Time, 1, Aeq_con_v);
        ROS_INFO("[TG] cal Aeq_con_v done, size: %ld, %ld", Aeq_con_v.rows(), Aeq_con_v.cols());
        // calculate acc continuity
        EigenSparse Aeq_con_a(m-1, n_all_poly);
        GenContinuityConstraint(m, d_order, Time, 2, Aeq_con_a);
        ROS_INFO("[TG] cal Aeq_con_a done, size: %ld, %ld", Aeq_con_a.rows(), Aeq_con_a.cols());
        // calculate jerk continuity
        EigenSparse Aeq_con_j(m-1, n_all_poly);
        GenContinuityConstraint(m, d_order, Time, 3, Aeq_con_j);
        ROS_INFO("[TG] cal Aeq_con_j done, size: %ld, %ld", Aeq_con_j.rows(), Aeq_con_j.cols());
        // concatenate Aeq & beq
        EigenSparse Aeq;
        VectorXd beq;
        int Aeq_row_num = 2*(d_num) + 5*(m-1);
        Aeq.resize(Aeq_row_num, n_all_poly);
        beq.resize(Aeq_row_num);
        std::vector<Eigen::Triplet<double, size_t>> triplets;
        triplets.reserve(Aeq_start.nonZeros() + Aeq_end.nonZeros() + Aeq_wp.nonZeros() + Aeq_con_p.nonZeros() + Aeq_con_v.nonZeros() + Aeq_con_a.nonZeros() + Aeq_con_j.nonZeros());
        for (Index c = 0; c < Aeq.cols(); ++c) {
            // ROS_INFO("[TG][OOQP] concatenating %ld th col of Aeq", c);
            Aeq.startVec(c);
            size_t row_indices = 0;
            for(EigenSparseIter it(Aeq_start, c); it; ++it) {
                Aeq.insertBack(it.row() + row_indices, c) = it.value();
            }
            row_indices += Aeq_start.rows();
            for(EigenSparseIter it(Aeq_end, c); it; ++it) {
                Aeq.insertBack(it.row() + row_indices, c) = it.value();
            }
            row_indices += Aeq_end.rows();
            for(EigenSparseIter it(Aeq_wp, c); it; ++it) {
                Aeq.insertBack(it.row() + row_indices, c) = it.value();
            }
            row_indices += Aeq_wp.rows();
            for(EigenSparseIter it(Aeq_con_p, c); it; ++it) {
                Aeq.insertBack(it.row() + row_indices, c) = it.value();
            }
            row_indices += Aeq_con_p.rows();
            for(EigenSparseIter it(Aeq_con_v, c); it; ++it) {
                Aeq.insertBack(it.row() + row_indices, c) = it.value();
            }
            row_indices += Aeq_con_v.rows();
            for(EigenSparseIter it(Aeq_con_a, c); it; ++it) {
                Aeq.insertBack(it.row() + row_indices, c) = it.value();
            }
            row_indices += Aeq_con_a.rows();
            for(EigenSparseIter it(Aeq_con_j, c); it; ++it) {
                Aeq.insertBack(it.row() + row_indices, c) = it.value();
            }
            row_indices += Aeq_con_j.rows();
            // ROS_INFO("[TG][OOQP] concatenate done !!!, row_indices: %ld", row_indices);
        }
        Aeq.finalize();
        Aeq.makeCompressed();
        LogData(Aeq, "Aeq");
        beq << beq_start,
               beq_end,
               beq_wp,
               VectorXd::Zero(4*(m-1));
        LogData(beq, "beq_" + std::to_string(idx));
        // LogData(beq_start, "beq_start_" + std::to_string(idx));
        // LogData(beq_end, "beq_end_" + std::to_string(idx));
        // LogData(beq_wp, "beq_wp_" + std::to_string(idx));
        ROS_WARN("[TG] concatenated Aeq & beq");
    /*   Produce the Minimum Snap cost function, the Hessian Matrix   */
        EigenSparse Q(m*p_num1d, m*p_num1d);
        for (int j = 0; j < m; ++j) {
            for (int k = 4; k <= p_order; ++k) {
                for (int l = 4; l <= p_order; ++l) {
                    Q.insert(k+j*p_num1d, l+j*p_num1d) = k*(k-1)*(k-2)*(k-3)*l*(l-1)*(l-2)*(l-3)*pow(Time(j), k+l-7)/(k+l-7);
                }
            }
        }
        LogData(Q, "Q");
        ROS_INFO("[TG] cal Q done, size: %ld, %ld", Q.rows(), Q.cols());
        // std::cout << Q << std::endl << std::endl;

        int solver_selection = 1;
        Eigen::VectorXd coeff;
        if (solver_selection == 0) {
            // solve using OOQP
            coeff = OOQPSolver(Q, Aeq, beq);
        } else {
            auto lower_l(beq);
            auto upper_l(beq);
            coeff = OSQPSolver(Q, Aeq, lower_l, upper_l);
        }

        for (int j = 0; j < m; ++j) {
            PolyCoeff.block(j, idx*p_num1d, 1, p_num1d) = coeff.segment(j*p_num1d, p_num1d).transpose();
            std::cout << "seg [" << j << "]" << std::endl;
            std::cout << PolyCoeff.block(j, idx*p_num1d, 1, p_num1d) << std::endl << std::endl;
        }
    }
    LogData(PolyCoeff, "ooqp_result");
    return PolyCoeff;
}

Eigen::VectorXd TrajectoryGeneratorWaypoint::OOQPSolver(
        const Eigen::SparseMatrix<double, Eigen::ColMajor>& Q,
        const Eigen::SparseMatrix<double, Eigen::ColMajor>& Aeq,
        const Eigen::VectorXd& beq) {

    // we are using eigen's sparse matrix 
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor> EigenSparse;
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator EigenSparseIter;
    // Make sure Q is in lower triangular form (Q is symmetric).
    EigenSparse Q_triangular = Q.triangularView<Lower>();
    // Compress sparse Eigen matrices (refer to Eigen Sparse Matrix user manual).
    Q_triangular.makeCompressed();
    LogData(Q_triangular, "Q_tri");
    // make copy
    auto beq_copy(beq);
    
    // ooqp related tasks
    ROS_INFO("[TG][OOQP] data prep");
    // fill in data for ooqp
    int nnzA = Aeq.nonZeros();
    int nnzQ = Q_triangular.nonZeros();
    int nnzC = 0;
    int my = beq.rows();
    int mz = 0;
    int nx = Q.rows(); // n_all_poly
    ROS_WARN("[TG][OOQP] nnzA:%d, nnzQ:%d, nnzC:%d, my:%d, mz:%d, nx:%d",nnzA,nnzQ,nnzC,my,mz,nx);

    double* cp   =  new double[nx]; //&ccopy.coeffRef(0);
    int* irowQ   =  new int[nnzQ]; //Q_triangular.innerIndexPtr(); // row indices for coln major
    int* jcolQ   =  new int[nnzQ]; //Q_triangular.outerIndexPtr(); // coln indices for coln major
    double* dQ   =  new double[nnzQ]; //Q_triangular.valuePtr();
    double* xlow =  new double[nx]; //&lowerLimitForX.coeffRef(0);
    char* ixlow  =  new char[nx]; //&useLowerLimitForX.coeffRef(0);
    double* xupp =  new double[nx]; //&upperLimitForX.coeffRef(0);
    char* ixupp  =  new char[nx]; //&useUpperLimitForX.coeffRef(0);
    int* irowA   =  new int[nnzA]; //Aeq.innerIndexPtr(); // row indices for coln major
    int* jcolA   =  new int[nnzA]; //Aeq.outerIndexPtr(); // coln indices for coln major
    double* dA   =  new double[nnzA]; //Aeq.valuePtr();
    double* bA   = &beq_copy.coeffRef(0);
    int* irowC   =  0; //Ccopy.outerIndexPtr();
    int* jcolC   =  0; //Ccopy.innerIndexPtr();
    double* dC   =  0; //Ccopy.valuePtr();
    double* clow =  0; //&lowerLimitForInequalityConstraints.coeffRef(0);
    char* iclow  =  0; //&useLowerLimitForInequalityConstraints.coeffRef(0);
    double* cupp =  0; //&upperLimitForInequalityConstraints.coeffRef(0);
    char* icupp  =  0; //&useUpperLimitForInequalityConstraints.coeffRef(0);

    // must sort all matrics (e.g. A, Q ...) into row-major forms!!!
    // doubleLexSort( irowA, nnzA, jcolA, dA );

    // fill in indices
    int idx_Aeq = 0;
    EigenSparse Aeq_trans = Aeq.transpose();
    for (int i = 0; i < Aeq_trans.outerSize(); ++i) {
        for (EigenSparseIter it(Aeq_trans, i); it; ++it) {
            irowA[idx_Aeq] = it.col();
            jcolA[idx_Aeq] = it.row();
            dA[idx_Aeq] = it.value();
            idx_Aeq += 1;
        }
    }
    int idx_Q = 0;
    EigenSparse Q_tri_trans = Q_triangular.transpose();
    for (int i = 0; i < Q_tri_trans.outerSize(); ++i) {
        for (EigenSparseIter it(Q_tri_trans, i); it; ++it) {
            irowQ[idx_Q] = it.col();
            jcolQ[idx_Q] = it.row();
            dQ[idx_Q] = it.value();
            idx_Q += 1;
        }
    }
    ROS_INFO("[TG][OOQP] idx_Aeq: %d, idx_Q: %d", idx_Aeq, idx_Q);
    for (int i = 0; i < nx; ++i) {
        cp[i] = 0;
        xlow[i] = 0;
        ixlow[i] = 0;
        xupp[i] = 0;
        ixupp[i] = 0;
    }
    ROS_INFO("[TG][OOQP] cp, xlow, ixlow, xupp, ixupp, all set to 0");

    // test for transpose
    // EigenSparse Aeq_trans = Aeq.transpose();
    // LogData(Aeq_trans, "Aeqt");
    // int idx_Aeqt = 0;
    // int *irowAt = new int[nnzA];
    // int *jcolAt = new int[nnzA];
    // double *dAt = Aeq_trans.valuePtr();
    // for (int i = 0; i < Aeq_trans.outerSize(); ++i) {
    //     for (EigenSparseIter it(Aeq_trans, i); it; ++it) {
    //         irowAt[idx_Aeqt] = it.row();
    //         jcolAt[idx_Aeqt] = it.col();
    //         idx_Aeqt += 1;
    //     }
    // }
    // LogData(dAt, nnzA, irowAt, jcolAt, "dAt");

    // log datas
    LogData(dA, nnzA, irowA, jcolA, "dA");
    int krowb[my] = {0};
    // LogData(bA, my, krowb, "bA_" + std::to_string(idx));
    LogData(dQ, nnzQ, irowQ, jcolQ, "dQ");

    bool monitor_flag = 0;
    int status = 0;
    VectorXd coeff(nx);

    QpGenSparseMa27 * qp = new QpGenSparseMa27(nx, my, mz, nnzQ, nnzA, nnzC);

    // QpGenData * prob = (QpGenData *) qp->makeData(cp, irowQ, jcolQ, dQ,
    //                                         xlow, ixlow, xupp, ixupp,
    //                                         irowA, jcolA, dA, bA,
    //                                         irowC, jcolC, dC,
    //                                         clow, iclow, cupp, icupp);
    QpGenData * prob = (QpGenData *) qp->copyDataFromSparseTriple(
                                            cp, irowQ, nnzQ, jcolQ, dQ,
                                            xlow, ixlow, xupp, ixupp,
                                            irowA, nnzA, jcolA, dA, bA,
                                            irowC, nnzC, jcolC, dC, 
                                            clow, iclow, cupp, icupp);
    // Create object to store problem variables.
    QpGenVars* vars = (QpGenVars*) qp->makeVariables(prob);

    // Create object to store problem residual data.
    QpGenResiduals* resid = (QpGenResiduals*) qp->makeResiduals(prob);

    // Create solver object.
    GondzioSolver* s = new GondzioSolver(qp, prob);

    // monitor itself
    if (monitor_flag == 1) {
        s->monitorSelf();
    }

    // Solve.
    status = s->solve(prob, vars, resid);

    if ((status == 0)) {
        vars->x->copyIntoArray(&coeff.coeffRef(0));
        ROS_WARN("[TG][OOQP] Solved !!!");
        // std::cout << coeff << std::endl << std::endl;
    }
    else if( status == 3)
        std::cout << "The program is probably infeasible, check the formulation.\n";
    else if (status == 4)
        std::cout << "Ther program is very slow in convergence, may have numerical issue.\n";
    else
        std::cout << "Don't know what the fuck it is, should not happen.\n";

    delete [] irowA;
    delete [] jcolA;
    delete [] irowQ;
    delete [] jcolQ;
    delete [] xlow;
    delete [] ixlow;
    delete [] xupp;
    delete [] ixupp;

    delete s;
    delete resid;
    delete vars;
    delete prob;
    delete qp;

    return coeff;
}

Eigen::VectorXd TrajectoryGeneratorWaypoint::OSQPSolver(
        const Eigen::SparseMatrix<double, Eigen::ColMajor>& Q,
        const Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
        const Eigen::VectorXd& lower_l,
        const Eigen::VectorXd& upper_l) {

    // we are using eigen's sparse matrix 
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor> EigenSparse;
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator EigenSparseIter;        
    // Make sure Q is in upper triangular form (Q is symmetric)(note osqp is different with ooqp!!!).
    EigenSparse Q_triangular = Q.triangularView<Upper>();
    // Compress sparse Eigen matrices (refer to Eigen Sparse Matrix user manual).
    Q_triangular.makeCompressed();
    LogData(Q_triangular, "Q_tri");
    // make copy
    auto Acopy(A);
    auto l_copy(lower_l), u_copy(upper_l);
    
    // osqp related tasks
    ROS_INFO("[TG][OSQP] data prep");
    // Load problem data
    c_int   n    = Q.rows(); // n_all_poly
    c_int   m    = A.rows(); // number of all constrains(eq & uneql)
    c_int   P_nnz= Q_triangular.nonZeros(); // maximum number of entries(non zero elements)
    c_float *P_x = Q_triangular.valuePtr(); // entries' values
    c_int   *P_i = new c_int[P_nnz]; // row indices in column major, size nzmax starting from 0
    c_int   *P_p = new c_int[n+1]; // column start indices(ptrs) in P_i, the last element is number of non-zeros
    c_float *q   = new c_float[n]; // dense array for linear part of cost function (size n)
    c_int   A_nnz= A.nonZeros();
    c_float *A_x = Acopy.valuePtr(); // entries' values
    c_int   *A_i = new c_int[A_nnz];
    c_int   *A_p = new c_int[n+1];
    c_float *l   = &l_copy.coeffRef(0); //new c_float[m]; // dense array for lower bound (size m)
    c_float *u   = &u_copy.coeffRef(0); //new c_float[m]; // dense array for lower bound (size m)
    ROS_WARN("[TG][OSQP] A_nnz:%lld, P_nnz:%lld, m:%lld, n:%lld",A_nnz,P_nnz,m,n);

    // fill in indices
    for (c_int i = 0; i < P_nnz; ++i) {
        P_i[i] = Q_triangular.innerIndexPtr()[i];
    }
    for (c_int i = 0; i < n; ++i) {
        P_p[i] = Q_triangular.outerIndexPtr()[i];
    }
    P_p[n] = P_nnz;

    for (c_int i = 0; i < A_nnz; ++i) {
        A_i[i] = A.innerIndexPtr()[i];
    }
    for (c_int i = 0; i < n; ++i) {
        A_p[i] = A.outerIndexPtr()[i];
    }
    A_p[n] = A_nnz;

    for (c_int i = 0; i < n; ++i) {
        q[i] = 0.0;
    }

    // log datas
    int *P_i_int = new int[P_nnz];
    std::copy(P_i, P_i + P_nnz, P_i_int);
    LogData(P_x, P_nnz, P_i_int, "P_x");
    int *P_p_int = new int[n+1];
    std::copy(P_p, P_p + n + 1, P_p_int);
    int P_p_indices[n+1] = {0};
    LogData(P_p_int, n+1, P_p_indices, "P_p");
    int *A_i_int = new int[A_nnz];
    std::copy(A_i, A_i + A_nnz, A_i_int);
    LogData(A_x, A_nnz, A_i_int, "A_x");
    int *A_p_int = new int[n+1];
    std::copy(A_p, A_p + n + 1, A_p_int);
    int A_p_indices[n+1] = {0};
    LogData(A_p_int, n+1, A_p_indices, "A_p");
    delete [] P_i_int;
    delete [] P_p_int;
    delete [] A_i_int;
    delete [] A_p_int;

    /*
    // create triplet csc matrix
    ROS_INFO("[TG][OSQP] Creating Triplets");
    csc *A_triplet = csc_spalloc(m, n, A_nnz, 0, 0);
    c_int idx_A = 0;
    c_int *ATtoC = new c_int[A_nnz];
    A_triplet->i = new c_int[A_nnz];
    A_triplet->p = new c_int[A_nnz];
    A_triplet->x = new c_float[A_nnz];
    EigenSparse A_trans = A.transpose();
    for (c_int i = 0; i < A_trans.outerSize(); ++i) {
        for (EigenSparseIter it(A_trans, i); it; ++it) {
            A_triplet->i[idx_A] = it.row();
            A_triplet->p[idx_A] = it.col();
            A_triplet->x[idx_A] = it.value();
            ATtoC[idx_A] = idx_A;
            idx_A += 1;
            // ROS_INFO("[TG][OSQP] idx_A: %lld", idx_A);
        }
    }
    csc *Q_triplet = csc_spalloc(n, n, P_nnz, 1, 0);
    c_int idx_Q = 0;
    c_int *QTtoC = new c_int[P_nnz];
    Q_triplet->i = new c_int[P_nnz];
    Q_triplet->p = new c_int[P_nnz];
    Q_triplet->x = new c_float[P_nnz];
    EigenSparse Q_trans = Q_triangular.transpose();
    for (c_int i = 0; i < Q_trans.outerSize(); ++i) {
        for (EigenSparseIter it(Q_trans, i); it; ++it) {
            Q_triplet->i[idx_Q] = it.row();
            Q_triplet->p[idx_Q] = it.col();
            Q_triplet->x[idx_Q] = it.value();
            QTtoC[idx_Q] = idx_Q;
            idx_Q += 1;
            // ROS_INFO("[TG][OSQP] idx_Q: %lld", idx_Q);
        }
    }
    ROS_WARN("[TG][OSQP] idx_A: %lld, idx_Q: %lld", idx_A, idx_Q);
    csc *A_csc = triplet_to_csc(A_triplet, ATtoC);
    csc *Q_csc = triplet_to_csc(Q_triplet, QTtoC);
    delete [] ATtoC;
    delete [] QTtoC;
    */

    // Exitflag
    c_int exitflag = 0;

    // Workspace structures
    OSQPWorkspace *work;
    OSQPSettings  *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
    OSQPData      *data     = (OSQPData *)c_malloc(sizeof(OSQPData));

    // Populate data
    if (data) {
        data->n = n; // number of columns
        data->m = m; // number of rows
        data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p); // csc stands for "Compressd Column", csr for "Compressd Row"
        data->q = q;
        data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
        data->l = l;
        data->u = u;
    }

    // Define solver settings as default
    if (settings) osqp_set_default_settings(settings);

    // Setup workspace
    exitflag = osqp_setup(&work, data, settings);
    if (exitflag != 0) {
        ROS_ERROR("[TG][OSQP] Setup Failed !!!!!!");
    } else {
        ROS_INFO("[TG][OSQP] Setup Sucessful");
    }

    // Solve Problem
    osqp_solve(work);

    // copy data from an array to Eigen vector
    c_float* solution = work->solution->x;
    VectorXd coeff = Eigen::Map<Eigen::VectorXd>(solution, work->data->n, 1);

    // Clean workspace
    osqp_cleanup(work);
    if (data) {
        if (data->A) c_free(data->A);
        if (data->P) c_free(data->P);
        c_free(data);
    }
    if (settings)  c_free(settings);

    delete [] P_i;
    delete [] P_p;
    delete [] A_i;
    delete [] A_p;

    return coeff;
}

void TrajectoryGeneratorWaypoint::GenContinuityConstraint(
        const int n_seq, 
        const int d_order, 
        const Eigen::VectorXd &Time, 
        const int order_deri, 
        Eigen::SparseMatrix<double, Eigen::ColMajor>& Aeq_con) {
    int p_order    = 2 * d_order - 1;              // the order of polynomial
    int p_num1d    = p_order + 1;
    int n_all_poly = n_seq * p_num1d;
    int k          = order_deri;
    Aeq_con.resize(n_seq-1, n_all_poly);
    for (int j = 0; j < n_seq-1; j++) {
        for (int i = k; i <= p_order; ++i) {
            Aeq_con.insert(j, j*p_num1d+i) = Factorial(i) / Factorial(i-k) * pow(Time(j), i-k); 
        }
        Aeq_con.insert(j,(j+1)*p_num1d+k) = - Factorial(k);
    }
}

void TrajectoryGeneratorWaypoint::LogData(const Eigen::MatrixXd& result, const std::string file_name) {
    std::string log_file_path = _package_path + "/output/" + file_name + ".csv";
    std::ofstream fp;
    fp.open(log_file_path, std::fstream::in | std::fstream::out | std::fstream::trunc);
    for (int i = 0; i < result.rows(); ++i) {
        for (int j = 0; j < result.cols(); ++j) {
            fp << result(i,j) << ",";
        }
        fp << std::endl;
    }
    fp.close();
}

void TrajectoryGeneratorWaypoint::LogData(const Eigen::VectorXd& result, const std::string file_name) {
    std::string log_file_path = _package_path + "/output/" + file_name + ".csv";
    std::ofstream fp;
    fp.open(log_file_path, std::fstream::in | std::fstream::out | std::fstream::trunc);
    for (int i = 0; i < result.size(); ++i) {
        fp << result(i) << ",";
        fp << std::endl;
    }
    fp.close();
}

void TrajectoryGeneratorWaypoint::LogData(const Eigen::SparseMatrix<double, Eigen::ColMajor>& result, const std::string file_name) {
    std::string log_file_path = _package_path + "/output/" + file_name + ".csv";
    std::ofstream fp;
    fp.open(log_file_path, std::fstream::in | std::fstream::out | std::fstream::trunc);
    MatrixXd result_dense(result);
    for (int i = 0; i < result_dense.rows(); ++i) {
        for (int j = 0; j < result_dense.cols(); ++j) {
            fp << result_dense(i,j) << ",";
        }
        fp << std::endl;
    }
    fp.close();
}

void TrajectoryGeneratorWaypoint::LogData(const double data[], const int size, const int idices[], const std::string file_name) {
    std::string log_file_path = _package_path + "/output/" + file_name + ".csv";
    std::ofstream fp;
    fp.open(log_file_path, std::fstream::in | std::fstream::out | std::fstream::trunc);
    for (int i = 0; i < size; ++i) {
        fp << data[i] << "," << idices[i] << "," << std::endl;
    }
    fp.close();
}

void TrajectoryGeneratorWaypoint::LogData(const double data[], const int size, const int rows[], const int cols[], const std::string file_name) {
    std::string log_file_path = _package_path + "/output/" + file_name + ".csv";
    std::ofstream fp;
    fp.open(log_file_path, std::fstream::in | std::fstream::out | std::fstream::trunc);
    for (int i = 0; i < size; ++i) {
        fp << data[i] << "," << rows[i] << "," << cols[i] << "," << std::endl;
    }
    fp.close();
}

void TrajectoryGeneratorWaypoint::LogData(const int data[], const int size, const int idices[], const std::string file_name) {
    std::string log_file_path = _package_path + "/output/" + file_name + ".csv";
    std::ofstream fp;
    fp.open(log_file_path, std::fstream::in | std::fstream::out | std::fstream::trunc);
    for (int i = 0; i < size; ++i) {
        fp << data[i] << "," << idices[i] << "," << std::endl;
    }
    fp.close();
}

void TrajectoryGeneratorWaypoint::LogData(const int data[], const int size, const int rows[], const int cols[], const std::string file_name) {
    std::string log_file_path = _package_path + "/output/" + file_name + ".csv";
    std::ofstream fp;
    fp.open(log_file_path, std::fstream::in | std::fstream::out | std::fstream::trunc);
    for (int i = 0; i < size; ++i) {
        fp << data[i] << "," << rows[i] << "," << cols[i] << "," << std::endl;
    }
    fp.close();
}