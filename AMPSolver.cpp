#include "AMPSolver.h"

void AMPSolver::SolveWithGS(std::vector<double>&& x, const int nCol, const std::vector<int>& ptr, std::vector<int>& col,
	const std::vector<double>& val, const std::vector<double>& rhs, double& errorOut, int& nItOut, const double tol)
  
{
    const int maxIt = 100;

	errorOut = 0.0;
	nItOut = 0;

    concurrency::array_view<const int, 1> ptr_gpu((int)ptr.size(), ptr.data());
    concurrency::array_view<const int, 1> col_gpu((int)col.size(), col.data());
    concurrency::array_view<const double, 1> val_gpu((int)val.size(), val.data());

    concurrency::array_view<const double, 1> rhs_gpu((int)rhs.size(), rhs.data());

    concurrency::array_view<double, 1> x_gpu(nCol, x.data());
    concurrency::array_view<double, 1> residual_gpu(nCol);

    residual_gpu.discard_data();

    
    for (int it = 0; it < maxIt; it++)
    {		
        concurrency::parallel_for_each(concurrency::extent<1>(1),
            [=](concurrency::index<1> idx) restrict(amp)
            {
                const double old0 = x_gpu[0];

                x_gpu[0] = rhs_gpu[0];

                for (int index = ptr_gpu[0] + 1; index < ptr_gpu[1]; index++)
                {
                    x_gpu[0] -= val_gpu[index] * x_gpu[col_gpu[index]];
                }
                x_gpu[0] /= val_gpu[ptr_gpu[0]];

                residual_gpu[0] = concurrency::precise_math::fabs(x_gpu[0] - old0);

                const double oldnm = x_gpu[nCol - 1];

                x_gpu[nCol - 1] = rhs_gpu[nCol - 1];

                for (int index = ptr_gpu[nCol - 1] + 1; index < ptr_gpu[nCol]; index++)
                {
                    x_gpu[nCol - 1] -= val_gpu[index] * x_gpu[col_gpu[index]];
                }
                x_gpu[nCol - 1] /= val_gpu[ptr_gpu[nCol - 1]];

                residual_gpu[nCol - 1] = concurrency::precise_math::fabs(x_gpu[nCol - 1] - oldnm);
            });

        // pares
		concurrency::parallel_for_each(concurrency::extent<1>((nCol/2)),
            [=](concurrency::index<1> idx) restrict(amp)
            {
                {
                    const concurrency::index<1> i = (2 * idx);

                    const double old = x_gpu[i];

                    x_gpu[i] = rhs_gpu[i];

                    for (int index = ptr_gpu[i] + 1; index < ptr_gpu[i + 1]; index++)
                    {
                        x_gpu[i] -= val_gpu[index] * x_gpu[col_gpu[index]];
                    }

                    x_gpu[i] /= val_gpu[ptr_gpu[i]];
                    
                    residual_gpu[i] = concurrency::precise_math::fabs(x_gpu[i] - old);
                }
			});


        //impares
        concurrency::parallel_for_each(concurrency::extent<1>((nCol /2)),
            [=](concurrency::index<1> idx) restrict(amp)
                {
                const concurrency::index<1> i = (2 * idx) + 1;

                const double old = x_gpu[i];

                x_gpu[i] = rhs_gpu[i];

                for (int index = ptr_gpu[i] + 1; index < ptr_gpu[i + 1]; index++)
                {
                    x_gpu[i] -= val_gpu[index] * x_gpu[col_gpu[index]];
                }

                x_gpu[i] /= val_gpu[ptr_gpu[i]];

                residual_gpu[i] = concurrency::precise_math::fabs(x_gpu[i] - old);               
        
                });
				

		for (int shift = nCol / 2; shift > 0; shift /= 2)
		{
			concurrency::parallel_for_each(concurrency::extent<1>(shift), [=](concurrency::index<1> idx) restrict(amp)
				{
					residual_gpu[idx] = concurrency::precise_math::fmax(residual_gpu[idx], residual_gpu[idx + shift]);

					if (shift % 2)
					{ //If odd, each thread includes a shifted entry. One will match the end of the queue
						residual_gpu[idx] = concurrency::precise_math::fmax(residual_gpu[idx], residual_gpu[idx + shift + 1]);
					}
				});
		}

		const double res = residual_gpu[0];

		if (res < tol || it + 1 >= maxIt)
		{
			nItOut = it;
			errorOut = res;
			break;
		}        
    }

    x_gpu.synchronize();

}