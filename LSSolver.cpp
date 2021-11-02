#include "LSSolver.h"

void LSSolver::SolveSparseCRS(std::vector<double>&& x, int nCol, const std::vector<int>& ptr, std::vector<int>& col,
    const std::vector<double>& val, const std::vector<double>& rhs, double& errorOut, int& nItOut, double tol)
{
	SolverWithoutAMG::params prm;
	prm.solver.tol = tol;

    SolverWithoutAMG solve(std::tie(nCol, ptr, col, val), prm);
    std::tie(nItOut, errorOut) = solve(std::tie(nCol, ptr, col, val), rhs, std::move(x));
}

LSSolver::LSSolver()
	:
	cPSolverPreconditioned(nullptr)
{}
LSSolver::~LSSolver()
{
	if (cPSolverPreconditioned)
	{
		delete cPSolverPreconditioned;
		cPSolverPreconditioned = nullptr;
	}
}

void LSSolver::PrecondtionCRS(int nCol, const std::vector<int>& ptr, std::vector<int>& col,
	const std::vector<double>& val, double tol)
{
	if (cPSolverPreconditioned)
	{
		delete cPSolverPreconditioned;
		cPSolverPreconditioned = nullptr;
	}

	Solver::params prm;
	prm.solver.tol = tol;
	cPSolverPreconditioned = new Solver(std::tie(nCol, ptr, col, val), prm);
}

void LSSolver::SolvePreconditionedCRS(std::vector<double>&& x, int nCol, const std::vector<int>& ptr, std::vector<int>& col,
	const std::vector<double>& val, const std::vector<double>& rhs, double& errorOut, int& nItOut)
{
	std::tie(nItOut, errorOut) = (*cPSolverPreconditioned)(std::tie(nCol, ptr, col, val), rhs, std::move(x));
}

void LSSolver::SolveSparseCRSWithAMG(std::vector<double>&& x, int nCol, const std::vector<int>& ptr, std::vector<int>& col,
	const std::vector<double>& val, const std::vector<double>& rhs, double& errorOut, int& nItOut,double tol)
{
	Solver::params prm;
	prm.solver.tol = tol;

	Solver solve(std::tie(nCol, ptr, col, val), prm);

	std::tie(nItOut, errorOut) = solve(rhs, std::move(x));
}