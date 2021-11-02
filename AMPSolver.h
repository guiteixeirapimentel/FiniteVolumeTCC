#pragma once
#include <vector>

#include <amp.h>
#include <amp_math.h>

#include <iostream>

class AMPSolver
{
public:
	AMPSolver() 
	{
		concurrency::accelerator gpu;
		
		auto descr = gpu.get_description();
		std::wcout << descr << std::endl;
		
		bool suportsDouble = gpu.supports_double_precision;
		std::cout << "suports double " << suportsDouble << std::endl;
	}
	~AMPSolver() { concurrency::amp_uninitialize(); }

	void SolveWithGS(std::vector<double>&& x, int nCol, const std::vector<int>& ptr, std::vector<int>& col,
		const std::vector<double>& val, const std::vector<double>& rhs, double& errorOut, int& nItOut, double tol = 1e-9);
private:
};