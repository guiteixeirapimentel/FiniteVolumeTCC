#pragma once

#define AMGCL_NO_BOOST
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

#include <amgcl/preconditioner/dummy.hpp>
#include "preconditionerDummy2.hpp"

#include "Timer.h"

typedef amgcl::backend::builtin<double> Backend;

typedef amgcl::make_solver<
    // Use AMG as preconditioner:
    amgcl::amg<
    Backend,
    amgcl::coarsening::smoothed_aggregation,
    amgcl::relaxation::spai0
    >,
    // And BiCGStab as iterative solver:
    amgcl::solver::bicgstab<Backend>
> Solver;

typedef amgcl::make_solver<
	amgcl::preconditioner::dummy2<Backend>,
	// And BiCGStab as iterative solver:
	amgcl::solver::bicgstab<Backend>
> SolverWithoutAMG;

class LSSolver
{
public:
    LSSolver();
    ~LSSolver();

    void SolveSparseCRSWithAMG(std::vector<double>&& x, int nCol, const std::vector<int>& ptr, std::vector<int>& col,
        const std::vector<double>& val, const std::vector<double>& rhs, double& errorOut, int& nItOut, double tol = 1e-9);

	void PrecondtionCRS(int nCol, const std::vector<int>& ptr, std::vector<int>& col,
		const std::vector<double>& val, double tol = 1e-9);

	void SolvePreconditionedCRS(std::vector<double>&& x, int nCol, const std::vector<int>& ptr, std::vector<int>& col,
		const std::vector<double>& val, const std::vector<double>& rhs, double& errorOut, int& nItOut);

	void SolveSparseCRS(std::vector<double>&& x, int nCol, const std::vector<int>& ptr, std::vector<int>& col,
		const std::vector<double>& val, const std::vector<double>& rhs, double& errorOut, int& nItOut, double tol = 1e-9);
	
private:
	Solver* cPSolverPreconditioned;
};