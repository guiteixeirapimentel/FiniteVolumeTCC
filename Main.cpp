#include <iostream>
#include <ctime>

//#define AMGCL_DEBUG

#include "FieldRectCylinder.h"
#include "FieldRectCylinder2D.h"
#include "LSSolver.h"
#include "AMPSolver.h"

#include "Timer.h"

#include <cinttypes>

#include <Windows.h>

typedef std::numeric_limits<double> doublelimits;

int SIMPLER();

int REMIXCSIMPROVED();

int REMIXCSIMPROVED2D();

template<class T>
T GetValue(const char* message);

int main()
{
	return REMIXCSIMPROVED2D();
}

int SIMPLER()
{
	Timer timer;
	const int tidSimple = timer.SetTimer("SIMPLE iteration");
	const int tidMomentum = timer.SetTimer("All 3 nonlinear momentum eqts");

	const int tidUMom = timer.SetTimer("U-Momentum iteration");
	const int tidVMom = timer.SetTimer("V-Momentum iteration");
	const int tidWMom = timer.SetTimer("W-Momentum iteration");

	const int tidPress = timer.SetTimer("Pressure LS solver");
	const int tidPressCorr = timer.SetTimer("Pressure Correction LS solver");

	const int tidCorrsResiduals = timer.SetTimer("Corrections and residuals calcs");

	const int tidCreatUMomLS = timer.SetTimer("Creation of UMomLS");
	const int tidCreatVMomLS = timer.SetTimer("Creation of VMomLS");
	const int tidCreatWMomLS = timer.SetTimer("Creation of WMomLS");

	const int tidCreatPress = timer.SetTimer("Creation of pressure LS");
	const int tidCreatPressCorr = timer.SetTimer("Creation of PressCorr LS");

	const int tidCalcHatValues = timer.SetTimer("Calc. of hat values");

	constexpr int NITMAXSIMPLE = 1000;

	constexpr double l = 0.02;

	constexpr double L = 21.0 * l;
	constexpr double H = 18.0 * l;
	constexpr double W = 10.0 * l;

	constexpr double w = W;
	constexpr double h = l;

	constexpr double dd = 0.1 * l;

	constexpr double dt = 0.1 * l;

	constexpr double maxTime = 10.0;
	constexpr int NITMAX = int(maxTime / dt) + 1;

	constexpr double rho = 1.0;
	constexpr double mu = (3.0 / 2.0) * 1e-5;

	constexpr double ufreestream = 0.15;
	constexpr double vfreestream = 0.0;
	constexpr double wfreestream = 0.0;

	constexpr double maxuResidual = 1e-4;
	constexpr double maxvResidual = 1e-4;
	constexpr double maxwResidual = 1e-4;
	constexpr double maxpResidual = 1e-4;

	constexpr double maxContinuityResidual = 1e-13;

	// FOR DEBUGGING PURPOSES

	constexpr int64_t NX = int64_t(L / dd) + 2;
	constexpr int64_t NY = int64_t(W / dd) + 2;
	constexpr int64_t NZ = int64_t(H / dd) + 2;

	constexpr int64_t NN = NX * NY * NZ;

	if (NN > 500000)
	{
		std::cout << "Atencao: o valor de pontos eh de " << NN << " sera necessario aprox. " << 0.7 * double(NN) / 1000000 << "GB de memoria" << std::endl;
		std::cout << "Continuar? (y/n): ";
		char b;
		std::cin >> b;
		if (b != 'y')
			return 0;
		else
			std::cout << "Continuando....\n";
	}

	constexpr double ReL = rho * ufreestream * L / mu;
	constexpr double ReD = rho * ufreestream * l / mu;

	const double delta = 5.0 * L / sqrt(ReL);

	constexpr double Peclet = rho * ufreestream / (mu / dd);

	constexpr double CFL = 3.0 * ufreestream * dt / dd;

	// FOR DEBUGGING PURPOSES

	FieldRectCylinder field(L, H, W, l, dd, rho, mu, ufreestream, vfreestream, wfreestream);

	field.SetBCFlatPlate();

	//field.ReadField("output/field10.010000");

	double tempoPassado = 0.0 + dt;

#ifdef _DEBUG
	field.SaveIKCutFieldToCSVFile(NY / 2, "OUT_DEBUG/initlong");
	field.SaveIJCutFieldToCSVFile(NZ / 2, "OUT_DEBUG/initlat");
#endif

	//FieldRectCylinder fieldnm = field;

	doubleField3D ulastIField = field.cUVelField;
	doubleField3D vlastIField = field.cVVelField;
	doubleField3D wlastIField = field.cWVelField;
	doubleField3D pLastIField = field.cPressureField;

	doubleField3D u0Field = field.cUVelField;
	doubleField3D v0Field = field.cVVelField;
	doubleField3D w0Field = field.cWVelField;

	std::vector<int>   ptr, col;
	std::vector<double> val, rhs;

	LSSolver linearSystemSolver;

	//AMPSolver linearSystemSolverAMP;


	bool showText = true;

	for (int itTime = 0; itTime < NITMAX; itTime++)
	{
		const time_t tSimple = time(nullptr);

		// SET INLET PERTUBATION
		{
			for (int j = 1; j < NY - 1; j++)
			{
				for (int k = 1; k < NZ - 1; k++)
				{
					field.cUVelField[0 + (k * NX) + (j * NX * NZ)] = ufreestream +
						(sin(tempoPassado + (k * dd)) * ufreestream * 0.01);

					field.cUVelField[1 + (k * NX) + (j * NX * NZ)] = ufreestream +
						(sin(tempoPassado + (k * dd)) * ufreestream * 0.01);
				}
			}
		}

		for (int itSimple = 0; itSimple < NITMAXSIMPLE; itSimple++)
		{
			timer.Tick(tidSimple);

			if (GetAsyncKeyState('H') & 0x8000)
			{
				showText = false;
				std::cout << "@@@@@ CHANGED ALL TEXT TO HIDE (HOLD H TO HIDE - S TO SHOW) @@@@@\n";
			}
			else if (GetAsyncKeyState('S') & 0x8000)
			{
				showText = true;
				std::cout << "@@@@@ CHANGED ALL TEXT TO SHOW (HOLD H TO HIDE - S TO SHOW) @@@@@\n";
			}

#ifdef _DEBUG
			double resSimplePosMom = 0.0;
			for (int j = 2; j < field.cNY - 1; j++)
				for (int k = 2; k < field.cNZ - 1; k++)
					for (int i = 2; i < field.cNX - 1; i++)
					{
						resSimplePosMom = fmax(
							fabs(
							(
								field.cUVelField[(j * NZ * NX) + (k * NX) + i] - field.cUVelField[(j * NZ * NX) + (k * NX) + i + 1] +
								field.cVVelField[(j * NZ * NX) + (k * NX) + i] - field.cVVelField[((j + 1) * NZ * NX) + (k * NX) + i] +
								field.cWVelField[(j * NZ * NX) + (k * NX) + i] - field.cWVelField[(j * NZ * NX) + ((k + 1) * NX) + i]
								)
							), resSimplePosMom);
			}

			std::cout << "Residual simple pre momentum eqts. " << resSimplePosMom << std::endl;
#endif

			//ulastIField = field.cUVelField;
			//vlastIField = field.cVVelField;
			//wlastIField = field.cWVelField;
			//pLastIField = field.cPressureField;

			memcpy(ulastIField.data(), field.cUVelField.data(), field.cUVelField.size() * sizeof(double));
			memcpy(vlastIField.data(), field.cVVelField.data(), field.cVVelField.size() * sizeof(double));
			memcpy(wlastIField.data(), field.cWVelField.data(), field.cWVelField.size() * sizeof(double));
			memcpy(pLastIField.data(), field.cPressureField.data(), field.cPressureField.size() * sizeof(double));

			time_t t1 = time(nullptr);

			// Calculate field hat values
			timer.Tick(tidCalcHatValues);

			field.CalcUHatValuesFI(u0Field, dt);
			field.CalcVHatValuesFI(v0Field, dt);
			field.CalcWHatValuesFI(w0Field, dt);

			timer.Tock(tidCalcHatValues);

			// Solve pressure equation
			timer.Tick(tidPress);

			timer.Tick(tidCreatPress);
			int nn = field.CreatePressureLSCSRFI(ptr, col, val, rhs);
			timer.Tock(tidCreatPress);

			double errPressLS = 0.0;
			int nitPressLS = 0;
			// solve LS Pressure eqt

			static int oneTimePrecon = 0;
			if (oneTimePrecon == 0)
			{
				linearSystemSolver.PrecondtionCRS(nn, ptr, col, val);

				oneTimePrecon = 1;
			}

			linearSystemSolver.SolvePreconditionedCRS(std::move(field.cPressureField), nn, ptr, col,
				val, rhs, errPressLS, nitPressLS);


			if (isnan(errPressLS))
			{
				std::cout << "Erro ls Pressure - nao foi possivel resolver ls\n";
				int x = 0;
				std::cin >> x;
			}

			timer.Tock(tidPress);

			if (showText)
			{
				std::cout << "Solved Pressure LS it " << nitPressLS << " residual " << errPressLS << std::endl;
			}

			// Solve Momentum equations		

			timer.Tick(tidMomentum);

			timer.Tick(tidUMom);

			timer.Tick(tidCreatUMomLS);
			nn = field.CreateUMomentumLSCSR(u0Field, dt, ptr, col, val, rhs);
			timer.Tock(tidCreatUMomLS);

			double errULS = 0.0;
			int nitULS = 0;

			// solve LS U
			linearSystemSolver.SolveSparseCRS(std::move(field.cUVelField),
				nn, ptr, col, val, rhs, errULS, nitULS, maxuResidual * 1e-3);

			//linearSystemSolverAMP.SolveWithGS(std::move(field.cUVelField),
			//	nn, ptr, col, val, rhs, errULS, nitULS, maxuResidual * 1e-3);


			if (isnan(errULS))
			{
				std::cout << "Erro ls U - nao foi possivel resolver ls\n";
				int x = 0;
				std::cin >> x;
			}

			timer.Tock(tidUMom);
			if (showText)
			{
				std::cout << "Solved U LS it " << nitULS << " residual " << errULS << std::endl;
			}

			timer.Tick(tidVMom);

			timer.Tick(tidCreatVMomLS);
			nn = field.CreateVMomentumLSCSRParallel(v0Field, dt, ptr, col, val, rhs);
			timer.Tock(tidCreatVMomLS);

			double errVLS = 0.0f;
			int nitVLS = 0;

			linearSystemSolver.SolveSparseCRS(std::move(field.cVVelField),
				nn, ptr, col, val, rhs, errVLS, nitVLS, maxvResidual * 1e-3);

			//linearSystemSolverAMP.SolveWithGS(std::move(field.cVVelField),
			//	nn, ptr, col, val, rhs, errVLS, nitVLS, maxvResidual * 1e-3);

			if (isnan(errVLS))
			{
				std::cout << "Erro ls V - nao foi possivel resolver ls\n";
				int x = 0;
				std::cin >> x;
			}

			timer.Tock(tidVMom);

			if (showText)
			{
				std::cout << "Solved V LS it " << nitVLS << " residual " << errVLS << std::endl;
			}

			timer.Tick(tidWMom);

			timer.Tick(tidCreatWMomLS);
			nn = field.CreateWMomentumLSCSRParallel(w0Field, dt, ptr, col, val, rhs);
			timer.Tock(tidCreatWMomLS);

			double errWLS = 0.0f;
			int nitWLS = 0;

			linearSystemSolver.SolveSparseCRS(std::move(field.cWVelField),
				nn, ptr, col, val, rhs, errWLS, nitWLS, maxwResidual * 1e-3);

			//linearSystemSolverAMP.SolveWithGS(std::move(field.cWVelField),
			//		nn, ptr, col, val, rhs, errWLS, nitWLS, maxwResidual * 1e-3);

			if (isnan(errWLS))
			{
				std::cout << "Erro ls W - nao foi possivel resolver ls\n";
				int x = 0;
				std::cin >> x;
			}

			timer.Tock(tidWMom);
			if (showText)
			{
				std::cout << "Solved W LS it " << nitWLS << " residual " << errWLS << std::endl;
			}

			timer.Tock(tidMomentum);

#ifdef _DEBUG
			field.SaveIKCutFieldToCSVFile(NY / 2, "OUT_DEBUG/posmomlong");
			field.SaveJKCutFieldToCSVFile(NX - 2, "OUT_DEBUG/posmomtrans");
#endif   

			// Solve Pressure Correction Equations
			{
				timer.Tick(tidCreatPressCorr);
				int nn = field.CreatePressureCorrectionLSCSR(ptr, col, val, rhs);
				timer.Tock(tidCreatPressCorr);

				double errPCLS = 0.0;
				int nitPCLS = 0;

				timer.Tick(tidPressCorr);

				if (((itSimple % 10) == 0) || (nitPCLS > 15 && (itSimple % 5) == 0))
				{
					std::cout << "Re-created precondtioner Press. corr LS" << std::endl;

					linearSystemSolver.PrecondtionCRS(nn, ptr, col, val, maxContinuityResidual * 1e-1);
				}


				linearSystemSolver.SolvePreconditionedCRS(std::move(field.cPressureCorrField),
					nn, ptr, col, val, rhs, errPCLS, nitPCLS);

				timer.Tock(tidPressCorr);
				if (isnan(errPCLS))
				{
					std::cout << "Erro ls PC - nao foi possivel resolver ls\n";
					int x = 0;
					std::cin >> x;

				}
				if (showText)
				{
					std::cout << "Solved pressure corr LS it " << nitPCLS << " residual " << errPCLS << std::endl;
				}
			}


			timer.Tick(tidCorrsResiduals);

			// CORRECT VELOCITIES
			// u
			for (int J = 1; J < NY - 1; J++)
			{
				for (int K = 1; K < NZ - 1; K++)
				{
					for (int i = 2; i < NX - 1; i++)
					{
						if (i >= field.ciCyInit && i < field.ciCyEnd + 1 && J >= field.cjCyInit && J < field.cjCyEnd &&
							K >= field.ckCyInit && K < field.ckCyEnd)
						{
							// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO
						}
						else
						{
							field.cUVelField[(J * NZ * NX) + (K * NX) + i] += (dd * dd / field.caijkU[(J * NZ * NX) + (K * NX) + i]) *
								(field.cPressureCorrField[(J * NZ * NX) + (K * NX) + i - 1] - field.cPressureCorrField[(J * NZ * NX) + (K * NX) + i]);
						}
					}
				}
			}
			//v
			for (int J = 2; J < NY - 1; J++)
			{
				for (int K = 1; K < NZ - 1; K++)
				{
					for (int i = 1; i < NX - 1; i++)
					{
						if (i >= field.ciCyInit && i < field.ciCyEnd && J >= field.cjCyInit && J < field.cjCyEnd &&
							K >= field.ckCyInit && K < field.ckCyEnd)
						{
							// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO
						}
						else
						{
							field.cVVelField[(J * NZ * NX) + (K * NX) + i] += (dd * dd / field.caijkV[(J * NZ * NX) + (K * NX) + i]) *
								(field.cPressureCorrField[((J - 1) * NZ * NX) + (K * NX) + i] - field.cPressureCorrField[(J * NZ * NX) + (K * NX) + i]);
						}
					}
				}
			}
			//w
			for (int J = 1; J < NY - 1; J++)
			{
				for (int K = 2; K < NZ - 1; K++)
				{
					for (int i = 1; i < NX - 1; i++)
					{
						if (i >= field.ciCyInit && i < field.ciCyEnd && J >= field.cjCyInit && J < field.cjCyEnd &&
							K >= field.ckCyInit && K < field.ckCyEnd + 1)
						{
							// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO
						}
						else
						{
							field.cWVelField[(J * NZ * NX) + (K * NX) + i] += (dd * dd / field.caijkW[(J * NZ * NX) + (K * NX) + i]) *
								(field.cPressureCorrField[(J * NZ * NX) + ((K - 1) * NX) + i] - field.cPressureCorrField[(J * NZ * NX) + (K * NX) + i]);
						}
					}
				}
			}

			double ures = 0.0;
			double vres = 0.0;
			double wres = 0.0;
			double pres = 0.0;

			// calc residuals velocities
			for (int j = 0; j < NY; j++)
				for (int k = 0; k < NZ; k++)
					for (int i = 0; i < NX; i++)
					{
						const double uResidual = fabs(field.cUVelField[(j * NZ * NX) + (k * NX) + i] - ulastIField[(j * NZ * NX) + (k * NX) + i]);

						if (uResidual > ures)
						{
							ures = uResidual;
						}

						const double vResidual = fabs(field.cVVelField[(j * NZ * NX) + (k * NX) + i] - vlastIField[(j * NZ * NX) + (k * NX) + i]);

						if (vResidual > vres)
						{
							vres = vResidual;
						}

						const double wResidual = fabs(field.cWVelField[(j * NZ * NX) + (k * NX) + i] - wlastIField[(j * NZ * NX) + (k * NX) + i]);

						if (wResidual > wres)
						{
							wres = wResidual;
						}
					}

			int iPres = 0;
			int jPres = 0;
			int kPres = 0;
			double pij = 0.0;
			double presAcum = 0.0;

			for (int j = 2; j < NY - 2; j++)
				for (int k = 2; k < NZ - 2; k++)
					for (int i = 2; i < NX - 2; i++)
					{
						const double residual = fabs((field.cPressureField[(j * NZ * NX) + (k * NX) + i] -
							pLastIField[(j * NZ * NX) + (k * NX) + i]));

						if (residual > pres)
						{
							pres = residual;
							iPres = i;
							jPres = j;
							kPres = k;
							pij = field.cPressureField[(j * NZ * NX) + (k * NX) + i];
						}
					}

			for (int j = 2; j < NY - 2; j++)
				for (int k = 2; k < NZ - 2; k++)
					for (int i = 2; i < NX - 2; i++)
					{
						const double residual = fabs((field.cPressureField[(j * NZ * NX) + (k * NX) + i] -
							pLastIField[(j * NZ * NX) + (k * NX) + i]));

						presAcum += residual;
					}

			double resSimple = 0.0;
			int iResSimple = 0;
			int jResSimple = 0;
			int kResSimple = 0;

			// SET/EXTRAPOLATE/ENFORCES OUTLET B.C.
			{
				for (int64_t j = 1; j < NY - 1; j++)
					for (int64_t k = 1; k < NZ - 1; k++)
					{
						const int64_t i = NX - 2;

						field.cUVelField[(j * NZ * NX) + (k * NX) + i + 1] = field.cUVelField[(j * NZ * NX) + (k * NX) + i] +
							field.cVVelField[(j * NZ * NX) + (k * NX) + i] - field.cVVelField[((j + 1) * NZ * NX) + (k * NX) + i] +
							field.cWVelField[(j * NZ * NX) + (k * NX) + i] - field.cWVelField[(j * NZ * NX) + ((k + 1) * NX) + i];
					}

				field.ExtrapolateBackwardVjK(NX - 1);
				field.ExtrapolateBackwardWJk(NX - 1);
			}

			for (int j = 2; j < field.cNY - 1; j++)
				for (int k = 2; k < field.cNZ - 1; k++)
					for (int i = 2; i < field.cNX - 1; i++)
					{
						const double residual = fabs(
							field.cUVelField[(j * NZ * NX) + (k * NX) + i] - field.cUVelField[(j * NZ * NX) + (k * NX) + i + 1] +
							field.cVVelField[(j * NZ * NX) + (k * NX) + i] - field.cVVelField[((j + 1) * NZ * NX) + (k * NX) + i] +
							field.cWVelField[(j * NZ * NX) + (k * NX) + i] - field.cWVelField[(j * NZ * NX) + ((k + 1) * NX) + i]
						) * field.cRHO * field.cdA;

						if (residual > resSimple)
						{
							resSimple = residual;

							iResSimple = i;
							jResSimple = j;
							kResSimple = k;
						}
					}

			timer.Tock(tidCorrsResiduals);

			if (showText)
			{
				std::cout << "##########################################\n";
				std::cout << "##########################################\n";

				std::cout << "resSimple: " << resSimple << " ures "
					<< ures << " vres " << vres << " wres " << wres << " pres " << pres << std::endl;
				std::cout << "maior mudanca em pres i " << iPres << " j " << jPres << " k " << kPres << " pij " << pij << std::endl;
				std::cout << "pres acumulado " << presAcum << std::endl;
				std::cout << "maior mudanca em res simple i " << iResSimple << " j " << jResSimple << " k " << kResSimple << std::endl;
			}

#ifdef _DEBUG
			field.SaveIKCutFieldToCSVFile(NY / 2, "OUT_DEBUG/poscorrlong");
			field.SaveJKCutFieldToCSVFile(NX - 2, "OUT_DEBUG/poscorrtrans");
#endif

			if ((resSimple < maxContinuityResidual && pres < maxpResidual) && ures < maxuResidual
				&& vres < maxvResidual && wres < maxwResidual)
			{
				std::cout << "Simple converged; time: " << tempoPassado << std::endl;
				std::cout << "Tempo para convergir simple: " << time(nullptr) - tSimple << " \n";

				std::cout << "##########################################\n";
				std::cout << "##########################################\n";

				if ((itTime % 250) == 0)
				{
					field.SaveIKCutFieldToCSVFile(NY / 2, "turbFlatPlate/corteLong" + std::to_string(tempoPassado));
					field.SaveJKCutFieldToCSVFile(NX / 2, "turbFlatPlate/corteTrans" + std::to_string(tempoPassado));

					field.SaveField("output/field" + std::to_string(tempoPassado));
				}

				break;
			}
			if (showText)
			{
				std::cout << "Last simple iteration spent " << time(nullptr) - t1 << " seconds\n";
			}

			timer.Tock(tidSimple);

			if (showText)
			{
				timer.WriteToCoutAllTickTock();
			}

		}

		u0Field = field.cUVelField;
		v0Field = field.cVVelField;
		w0Field = field.cWVelField;

		tempoPassado += dt;
		if (showText)
		{
			std::cout << "Tempo aprox restante para conclusao : " << (NITMAX - itTime) * (time(nullptr) - tSimple) << " segs  ou "
				<< ((NITMAX - itTime) * (time(nullptr) - tSimple)) / 3600 << "horas" << std::endl;
		}

	}

	return 0;
}

int REMIXCSIMPROVED()
{
	Timer timer;
	const int tidSimple = timer.SetTimer("SIMPLE iteration");
	
	const int tidCalcHatValues = timer.SetTimer("Calc. of hat values");
	
	const int tidCreatPress = timer.SetTimer("Creation of pressure LS");
	const int tidPress = timer.SetTimer("Pressure LS solver");

	const int tidCorrsResiduals = timer.SetTimer("Corrections and residuals calcs");
	
	constexpr int NITMAXSIMPLE = 1000;
	
	double dd = 0.001;

	double l = dd * 20;

	double L = 26.0 * l;
	double H = 15.0 * l;
	double W = 10.0 * l;

	double w = W;
	double h = l;
		
	double dt = 0.0005;

	double maxTime = 15.0;
	
	double rho = 1.0;
	double mu = 3.0e-5;

	double ufreestream = 0.15;
	double vfreestream = 0.0;
	double wfreestream = 0.0;

	double maxuResidual = 5e-5;
	double maxvResidual = 5e-5;
	double maxwResidual = 5e-5;
	double maxpResidual = 5e-5;

	double maxContinuityResidual = 1e-9;
	
	const char answer = GetValue<char>("Recomecar? (y/n): ");
	std::string baseFileName = "";
	double initialTime = 0.0;

	if (answer == 'y' || answer == 'Y')
	{
		baseFileName = GetValue<std::string>("Nome base arqs .dat (output/field10.010000): ");
		initialTime = GetValue<double>("Tempo inicial representado pelo arquivo (10.01): ");

		l = GetValue<double>("Diametro cilindro (l): ");

		const double nL = GetValue<double>("Comp. relativo do campo (L/l)[19.0-21.0]: ");
		const double nH = GetValue<double>("Alt. relativa do campo (H/l)[15.0-18.0]: ");
		const double nW = GetValue<double>("Larg. relativa do campo (W/l)[5.0-10.0]: ");

		L = nL * l;
		H = nH * l;
		W = nW * l;

		w = W;
		h = l;

		dd = GetValue<double>("dd (10 % * l): ");

		dt = GetValue<double>("dt (lookout for 0.05<CFL<0.30)");

		maxTime = GetValue<double>("Max time (maxTime): ");
		
		rho = GetValue<double>("Densidade (1.0): ");
		mu = GetValue<double>("Mu [(3/2)*1e-5]: ");

		ufreestream = GetValue<double>("U freestream (0.15): ");
		vfreestream = GetValue<double>("V freestream (0.0): ");
		wfreestream = GetValue<double>("W freestream (0.0): ");
				
		const double CFL = 3.0 * ufreestream * dt / dd;
		const double ReD = rho * ufreestream * l / mu;
		
		std::cout << "CFL " << CFL << std::endl;
		std::cout << "ReD " << ReD << std::endl;
	}


	const int NITMAX = int(maxTime / dt) + 2;

	std::cout << "Max. iteracoes tempo: " << NITMAX << std::endl;


	// FOR DEBUGGING PURPOSES

	const int64_t NX = int64_t(L / dd) + 2;
	const int64_t NY = int64_t(W / dd) + 2;
	const int64_t NZ = int64_t(H / dd) + 2;

	const int64_t NN = NX * NY * NZ;

	const double ReL = rho * ufreestream * L / mu;
	const double ReD = rho * ufreestream * l / mu;

	const double delta = 5.0 * L / sqrt(ReL);

	const double Peclet = rho * ufreestream / (mu / dd);

	const double CFL = 3.0 * ufreestream * dt / dd;

	const double f0 = 0.16 * ufreestream / l;
	const double f1 = 0.2 * ufreestream / l;

	const double T0 = 1.0 / f0;
	const double T1 = 1.0 / f1;

	const double saveEverySec = (T1 / 10.0);

	const int nItToSave = 366;// (saveEverySec / dt);

	std::cout << "Atencao: o valor de pontos eh de " << NN << " sera necessario aprox. " << 0.7 * double(NN) / 1000000 << "GB de memoria" << std::endl;

	if (NN > 500000)
	{
		std::cout << "Continuar? (y/n): ";
		char b;
		std::cin >> b;
		if (b != 'y')
			return 0;
		else
			std::cout << "Continuando....\n";
	}
		
	// FOR DEBUGGING PURPOSES

	FieldRectCylinder field(L, H, W, l, dd, rho, mu, ufreestream, vfreestream, wfreestream);

	field.SetBCFlatPlate();

	if (answer == 'y' || answer == 'Y')
	{
		field.ReadField(baseFileName);

		std::cout << "Leu " << baseFileName << std::endl;
	}
	

	double tempoPassado = initialTime + dt;
	
#ifdef _DEBUG
	//field.SaveIKCutFieldToCSVFile(NY / 2, "OUT_DEBUG/initlong");
	//field.SaveIJCutFieldToCSVFile(NZ / 2, "OUT_DEBUG/initlat");
#endif

	//FieldRectCylinder fieldnm = field;

	doubleField3D ulastIField = field.cUVelField;
	doubleField3D vlastIField = field.cVVelField;
	doubleField3D wlastIField = field.cWVelField;
	doubleField3D pLastIField = field.cPressureField;

	doubleField3D u0Field = field.cUVelField;
	doubleField3D v0Field = field.cVVelField;
	doubleField3D w0Field = field.cWVelField;

	std::vector<int>   ptr, col;
	std::vector<double> val, rhs;

	LSSolver linearSystemSolver;
	   
	bool showText = true;

	int nItPressure = 100;
	const int maxNItPressureToReCondition = 40;

	for (int itTime = (tempoPassado/dt); itTime < NITMAX; itTime++)
	{
		const time_t tSimple = time(nullptr);

		// SET INLET PERTUBATION
		{
			for (int j = 1; j < NY - 1; j++)
			{
				for (int k = 1; k < NZ - 1; k++)
				{
					field.cUVelField[0 + (k * NX) + (j * NX * NZ)] = ufreestream +
						(sin(tempoPassado + (k * dd)) * ufreestream * 0.01);
		
					field.cUVelField[1 + (k * NX) + (j * NX * NZ)] = ufreestream +
						(sin(tempoPassado + (k * dd)) * ufreestream * 0.01);
				}
			}
		}

		for (int itSimple = 0; itSimple < NITMAXSIMPLE; itSimple++)
		{
			timer.Tick(tidSimple);

			if (GetAsyncKeyState('H') & 0x8000)
			{
				showText = false;
				std::cout << "@@@@@ CHANGED ALL TEXT TO HIDE (HOLD H TO HIDE - S TO SHOW) @@@@@\n";
			}
			else if (GetAsyncKeyState('S') & 0x8000)
			{
				showText = true;
				std::cout << "@@@@@ CHANGED ALL TEXT TO SHOW (HOLD H TO HIDE - S TO SHOW) @@@@@\n";
			}

#ifdef _DDEBUG
			double resSimplePosMom = 0.0;
			for (int j = 2; j < field.cNY - 1; j++)
				for (int k = 2; k < field.cNZ - 1; k++)
					for (int i = 2; i < field.cNX - 1; i++)
					{
						resSimplePosMom = fmax(
							fabs(
							(
								field.cUVelField[(j * NZ * NX) + (k * NX) + i] - field.cUVelField[(j * NZ * NX) + (k * NX) + i + 1] +
								field.cVVelField[(j * NZ * NX) + (k * NX) + i] - field.cVVelField[((j + 1) * NZ * NX) + (k * NX) + i] +
								field.cWVelField[(j * NZ * NX) + (k * NX) + i] - field.cWVelField[(j * NZ * NX) + ((k + 1) * NX) + i]
								)
							), resSimplePosMom);
					}

			std::cout << "Residual simple pre momentum eqts. " << resSimplePosMom << std::endl;
#endif

			//ulastIField = field.cUVelField;
			//vlastIField = field.cVVelField;
			//wlastIField = field.cWVelField;
			//pLastIField = field.cPressureField;

			memcpy(ulastIField.data(), field.cUVelField.data(), field.cUVelField.size() * sizeof(double));
			memcpy(vlastIField.data(), field.cVVelField.data(), field.cVVelField.size() * sizeof(double));
			memcpy(wlastIField.data(), field.cWVelField.data(), field.cWVelField.size() * sizeof(double));
			memcpy(pLastIField.data(), field.cPressureField.data(), field.cPressureField.size() * sizeof(double));

			time_t t1 = time(nullptr);

			// Calculate field hat values
			timer.Tick(tidCalcHatValues);

			field.CalcUHatValuesCN(u0Field, dt);
			field.CalcVHatValuesCN(v0Field, dt);
			field.CalcWHatValuesCN(w0Field, dt);

			timer.Tock(tidCalcHatValues);

			// Solve pressure equation
			timer.Tick(tidPress);

			timer.Tick(tidCreatPress);
			int nn = field.CreatePressureLSCSRFI(ptr, col, val, rhs);
			timer.Tock(tidCreatPress);

			double errPressLS = 0.0;
			int nitPressLS = 0;
			// solve LS Pressure eqt

			if (nItPressure > maxNItPressureToReCondition)
			{
				linearSystemSolver.PrecondtionCRS(nn, ptr, col, val, maxContinuityResidual * 1e-3);
			}

			linearSystemSolver.SolvePreconditionedCRS(std::move(field.cPressureField), nn, ptr, col,
				val, rhs, errPressLS, nitPressLS);

			//linearSystemSolver.SolveSparseCRS(std::move(field.cPressureField), nn, ptr, col, val, rhs, errPressLS, nitPressLS);

			nItPressure = nitPressLS;
			

			if (isnan(errPressLS))
			{
				std::cout << "Erro ls Pressure - nao foi possivel resolver ls\n";
				int x = 0;
				std::cin >> x;
			}

			timer.Tock(tidPress);

			if (showText)
			{
				std::cout << "Solved Pressure LS it " << nitPressLS << " residual " << errPressLS << std::endl;
			}
			
#ifdef _DEBUG
			field.SaveIKCutFieldToCSVFile(NY / 2, "OUT_DEBUG/posmomlong");
			field.SaveJKCutFieldToCSVFile(NX - 2, "OUT_DEBUG/posmomtrans");
#endif   

			

			timer.Tick(tidCorrsResiduals);

			// CALCULATE TRUE VELOCITIES
			// u
			for (int J = 1; J < NY - 1; J++)
			{
				for (int K = 1; K < NZ - 1; K++)
				{
					for (int i = 2; i < NX - 1; i++)
					{
						if (i >= field.ciCyInit && i < field.ciCyEnd + 1 && J >= field.cjCyInit && J < field.cjCyEnd &&
							K >= field.ckCyInit && K < field.ckCyEnd)
						{
							// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO
						}
						else
						{
							field.cUVelField[(J * NZ * NX) + (K * NX) + i] = field.cUHatVelField[(J * NZ * NX) + (K * NX) + i]
								+ ((field.cPressureField[(J * NZ * NX) + (K * NX) + i - 1] - field.cPressureField[(J * NZ * NX) + (K * NX) + i])
										*(dd * dd / field.caijkU[(J * NZ * NX) + (K * NX) + i]));
								
								/*(dd * dd / field.caijkU[(J * NZ * NX) + (K * NX) + i]) *
								(field.cPressureCorrField[(J * NZ * NX) + (K * NX) + i - 1] - field.cPressureCorrField[(J * NZ * NX) + (K * NX) + i]);
								*/
						}
					}
				}
			}
			//v
			for (int J = 2; J < NY - 1; J++)
			{
				for (int K = 1; K < NZ - 1; K++)
				{
					for (int i = 1; i < NX - 1; i++)
					{
						if (i >= field.ciCyInit && i < field.ciCyEnd && J >= field.cjCyInit && J < field.cjCyEnd &&
							K >= field.ckCyInit && K < field.ckCyEnd)
						{
							// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO
						}
						else
						{
							//field.cVVelField[(J * NZ * NX) + (K * NX) + i] += (dd * dd / field.caijkV[(J * NZ * NX) + (K * NX) + i]) *
							//	(field.cPressureCorrField[((J - 1) * NZ * NX) + (K * NX) + i] - field.cPressureCorrField[(J * NZ * NX) + (K * NX) + i]);

							field.cVVelField[(J * NZ * NX) + (K * NX) + i] = field.cVHatVelField[(J * NZ * NX) + (K * NX) + i]
								+ ((field.cPressureField[((J-1) * NZ * NX) + (K * NX) + i] - field.cPressureField[(J * NZ * NX) + (K * NX) + i])
								* (dd * dd / field.caijkV[(J * NZ * NX) + (K * NX) + i]));
						}
					}
				}
			}
			//w
			for (int J = 1; J < NY - 1; J++)
			{
				for (int K = 2; K < NZ - 1; K++)
				{
					for (int i = 1; i < NX - 1; i++)
					{
						if (i >= field.ciCyInit && i < field.ciCyEnd && J >= field.cjCyInit && J < field.cjCyEnd &&
							K >= field.ckCyInit && K < field.ckCyEnd + 1)
						{
							// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO
						}
						else
						{
							//field.cWVelField[(J * NZ * NX) + (K * NX) + i] += (dd * dd / field.caijkW[(J * NZ * NX) + (K * NX) + i]) *
							//	(field.cPressureCorrField[(J * NZ * NX) + ((K - 1) * NX) + i] - field.cPressureCorrField[(J * NZ * NX) + (K * NX) + i]);
							
							field.cWVelField[(J * NZ * NX) + (K * NX) + i] = field.cWHatVelField[(J * NZ * NX) + (K * NX) + i]
								+ ((field.cPressureField[(J * NZ * NX) + ((K - 1) * NX) + i] - field.cPressureField[(J * NZ * NX) + (K * NX) + i])
									* (dd * dd / field.caijkW[(J * NZ * NX) + (K * NX) + i]));
						}
					}
				}
			}

			double ures = 0.0;
			double vres = 0.0;
			double wres = 0.0;
			double pres = 0.0;

			double sumCFL = 0.0;

			double CFLmax = 0.0;

			// calc residuals velocities
			for (int j = 0; j < NY; j++)
				for (int k = 0; k < NZ; k++)
					for (int i = 0; i < NX; i++)
					{
						const double uResidual = fabs(field.cUVelField[(j * NZ * NX) + (k * NX) + i] - ulastIField[(j * NZ * NX) + (k * NX) + i]);

						if (uResidual > ures)
						{
							ures = uResidual;
						}

						const double vResidual = fabs(field.cVVelField[(j * NZ * NX) + (k * NX) + i] - vlastIField[(j * NZ * NX) + (k * NX) + i]);

						if (vResidual > vres)
						{
							vres = vResidual;
						}

						const double wResidual = fabs(field.cWVelField[(j * NZ * NX) + (k * NX) + i] - wlastIField[(j * NZ * NX) + (k * NX) + i]);

						if (wResidual > wres)
						{
							wres = wResidual;
						}

						const double cfl = (fabs(field.cUVelField[(j * NZ * NX)]) + fabs(field.cVVelField[(j * NZ * NX) + (k * NX) + i])
							+ fabs(field.cWVelField[(j * NZ * NX) + (k * NX) + i])) * dt / dd;

						sumCFL += cfl;

						if (cfl > CFLmax)
							CFLmax = cfl;
					}

			const double CFLMean = sumCFL / (field.cNX * field.cNY * field.cNZ);

			int iPres = 0;
			int jPres = 0;
			int kPres = 0;
			double pij = 0.0;
			double presAcum = 0.0;

			for (int j = 2; j < NY - 2; j++)
				for (int k = 2; k < NZ - 2; k++)
					for (int i = 2; i < NX - 2; i++)
					{
						const double residual = fabs((field.cPressureField[(j * NZ * NX) + (k * NX) + i] -
							pLastIField[(j * NZ * NX) + (k * NX) + i]));

						if (residual > pres)
						{
							pres = residual;
							iPres = i;
							jPres = j;
							kPres = k;
							pij = field.cPressureField[(j * NZ * NX) + (k * NX) + i];
						}
					}

			for (int j = 2; j < NY - 2; j++)
				for (int k = 2; k < NZ - 2; k++)
					for (int i = 2; i < NX - 2; i++)
					{
						const double residual = fabs((field.cPressureField[(j * NZ * NX) + (k * NX) + i] -
							pLastIField[(j * NZ * NX) + (k * NX) + i]));

						presAcum += residual;
					}

			double resSimple = 0.0;
			int iResSimple = 0;
			int jResSimple = 0;
			int kResSimple = 0;

			// SET/EXTRAPOLATE/ENFORCES OUTLET B.C.
			{
				for (int64_t j = 1; j < NY - 1; j++)
					for (int64_t k = 1; k < NZ - 1; k++)
					{
						const int64_t i = NX - 2;

						field.cUVelField[(j * NZ * NX) + (k * NX) + i + 1] = field.cUVelField[(j * NZ * NX) + (k * NX) + i] +
							field.cVVelField[(j * NZ * NX) + (k * NX) + i] - field.cVVelField[((j + 1) * NZ * NX) + (k * NX) + i] +
							field.cWVelField[(j * NZ * NX) + (k * NX) + i] - field.cWVelField[(j * NZ * NX) + ((k + 1) * NX) + i];
					}

				field.ExtrapolateBackwardVjK(NX - 1);
				field.ExtrapolateBackwardWJk(NX - 1);
			}

			for (int j = 2; j < field.cNY - 1; j++)
				for (int k = 2; k < field.cNZ - 1; k++)
					for (int i = 2; i < field.cNX - 1; i++)
					{
						const double residual = fabs(
							field.cUVelField[(j * NZ * NX) + (k * NX) + i] - field.cUVelField[(j * NZ * NX) + (k * NX) + i + 1] +
							field.cVVelField[(j * NZ * NX) + (k * NX) + i] - field.cVVelField[((j + 1) * NZ * NX) + (k * NX) + i] +
							field.cWVelField[(j * NZ * NX) + (k * NX) + i] - field.cWVelField[(j * NZ * NX) + ((k + 1) * NX) + i]
						) * field.cRHO * field.cdA;

						if (residual > resSimple)
						{
							resSimple = residual;

							iResSimple = i;
							jResSimple = j;
							kResSimple = k;
						}
					}

			timer.Tock(tidCorrsResiduals);

			if (showText)
			{
				std::cout << "##########################################\n";
				std::cout << "##########################################\n";

				std::cout << "resSimple: " << resSimple << " ures "
					<< ures << " vres " << vres << " wres " << wres << " pres " << pres << std::endl;
				std::cout << "maior mudanca em pres i " << iPres << " j " << jPres << " k " << kPres << " pij " << pij << std::endl;
				std::cout << "pres acumulado " << presAcum << std::endl;
				std::cout << "maior mudanca em res simple i " << iResSimple << " j " << jResSimple << " k " << kResSimple << std::endl;
				std::cout << "CFL medio " << CFLMean << " maior CFL " << CFLmax << std::endl;
			}

#ifdef _DEBUG
			field.SaveIKCutFieldToCSVFile(NY / 2, "OUT_DEBUG/poscorrlong");
			field.SaveJKCutFieldToCSVFile(NX - 2, "OUT_DEBUG/poscorrtrans");
#endif

			if ((resSimple < maxContinuityResidual && pres < maxpResidual) && ures < maxuResidual
				&& vres < maxvResidual && wres < maxwResidual)
			{
				std::cout << "Simple converged; time: " << tempoPassado << std::endl;
				std::cout << "Tempo para convergir simple: " << time(nullptr) - tSimple << " \n";

				std::cout << "##########################################\n";
				std::cout << "##########################################\n";

				if (((itTime % nItToSave) == 0) || (tempoPassado > 12.0 && (itTime % (nItToSave/6)) == 0))
				{
					field.SaveIKCutFieldToCSVFile(NY / 2, "turbFlatPlate/corteLong" + std::to_string(tempoPassado));
					field.SaveJKCutFieldToCSVFile(NX / 2, "turbFlatPlate/corteTrans" + std::to_string(tempoPassado));

					field.SaveField("output/field" + std::to_string(tempoPassado));
				}

				break;
			}
			if (showText)
			{
				std::cout << "Last simple iteration spent " << time(nullptr) - t1 << " seconds\n";
			}

			timer.Tock(tidSimple);

			if (showText)
			{
				timer.WriteToCoutAllTickTock();
			}

		}

		u0Field = field.cUVelField;
		v0Field = field.cVVelField;
		w0Field = field.cWVelField;

		tempoPassado += dt;
		if (showText)
		{
			std::cout << "Tempo aprox restante para conclusao : " << (NITMAX - itTime) * (time(nullptr) - tSimple) << " segs  ou "
				<< ((NITMAX - itTime) * (time(nullptr) - tSimple)) / 3600 << "horas" << std::endl;
		}

	}

	return 0;
}


int REMIXCSIMPROVED2D()
{
	Timer timer;
	const int tidSimple = timer.SetTimer("SIMPLE iteration");

	const int tidCalcHatValues = timer.SetTimer("Calc. of hat values");

	const int tidCreatPress = timer.SetTimer("Creation of pressure LS");
	const int tidPress = timer.SetTimer("Pressure LS solver");

	const int tidCorrsResiduals = timer.SetTimer("Corrections and residuals calcs");

	constexpr int NITMAXSIMPLE = 1000;

	double dd = 0.002/4;

	double l = dd * 10 * 4;

	double L = 26.0 * l;
	double H = 15.0 * l;

	double h = l;

	double dt = 0.001/4;

	double maxTime = 15.0;

	double rho = 1.0;
	double mu = 3.0e-5;

	double ufreestream = 0.15;
	double wfreestream = 0.0;

	double maxuResidual = 5e-5;
	double maxwResidual = 5e-5;
	double maxpResidual = 5e-5;

	double maxContinuityResidual = 1e-5;

	const char answer = GetValue<char>("Recomecar? (y/n): ");
	std::string baseFileName = "";
	double initialTime = 0.0;

	if (answer == 'y' || answer == 'Y')
	{
		baseFileName = GetValue<std::string>("Nome base arqs .dat (output/field10.010000): ");
		initialTime = GetValue<double>("Tempo inicial representado pelo arquivo (10.01): ");

		l = GetValue<double>("Diametro cilindro (l): ");

		const double nL = GetValue<double>("Comp. relativo do campo (L/l)[19.0-21.0]: ");
		const double nH = GetValue<double>("Alt. relativa do campo (H/l)[15.0-18.0]: ");

		L = nL * l;
		H = nH * l;

		h = l;

		dd = GetValue<double>("dd (10 % * l): ");

		dt = GetValue<double>("dt (lookout for 0.05<CFL<0.30)");

		maxTime = GetValue<double>("Max time (maxTime): ");

		rho = GetValue<double>("Densidade (1.0): ");
		mu = GetValue<double>("Mu [(3/2)*1e-5]: ");

		ufreestream = GetValue<double>("U freestream (0.15): ");
		wfreestream = GetValue<double>("W freestream (0.0): ");

		const double CFL = 3.0 * ufreestream * dt / dd;
		const double ReD = rho * ufreestream * l / mu;

		std::cout << "CFL " << CFL << std::endl;
		std::cout << "ReD " << ReD << std::endl;
	}


	const int NITMAX = int(maxTime / dt) + 2;

	std::cout << "Max. iteracoes tempo: " << NITMAX << std::endl;


	// FOR DEBUGGING PURPOSES

	const int64_t NX = int64_t(L / dd) + 2;
	const int64_t NZ = int64_t(H / dd) + 2;

	const int64_t NN = NX * NZ;

	const double ReL = rho * ufreestream * L / mu;
	const double ReD = rho * ufreestream * l / mu;

	const double delta = 5.0 * L / sqrt(ReL);

	const double Peclet = rho * ufreestream / (mu / dd);

	const double CFL = 3.0 * ufreestream * dt / dd;

	const double f0 = 0.16 * ufreestream / l;
	const double f1 = 0.2 * ufreestream / l;

	const double T0 = 1.0 / f0;
	const double T1 = 1.0 / f1;

	const double saveEverySec = (T1 / 10.0);

	const int nItToSave = 366*4;// (saveEverySec / dt);

	std::cout << "Atencao: o valor de pontos eh de " << NN << " sera necessario aprox. " << 0.7 * double(NN) / 1000000 << "GB de memoria" << std::endl;

	if (NN > 500000)
	{
		std::cout << "Continuar? (y/n): ";
		char b;
		std::cin >> b;
		if (b != 'y')
			return 0;
		else
			std::cout << "Continuando....\n";
	}

	// FOR DEBUGGING PURPOSES

	FieldRectCylinder2D field(L, H, l, dd, rho, mu, ufreestream, wfreestream);

	field.SetBCFlatPlate();

	if (answer == 'y' || answer == 'Y')
	{
		field.ReadField(baseFileName);

		std::cout << "Leu " << baseFileName << std::endl;
	}


	double tempoPassado = initialTime + dt;

#ifdef _DEBUG
	//field.SaveIKCutFieldToCSVFile(NY / 2, "OUT_DEBUG/initlong");
	//field.SaveIJCutFieldToCSVFile(NZ / 2, "OUT_DEBUG/initlat");
#endif

	//FieldRectCylinder fieldnm = field;

	doubleField2D ulastIField = field.cUVelField;
	doubleField2D wlastIField = field.cWVelField;
	doubleField2D pLastIField = field.cPressureField;

	doubleField2D u0Field = field.cUVelField;
	doubleField2D w0Field = field.cWVelField;

	std::vector<int>   ptr, col;
	std::vector<double> val, rhs;

	LSSolver linearSystemSolver;

	bool showText = true;

	int nItPressure = 100;
	const int maxNItPressureToReCondition = 40;

	for (int itTime = (tempoPassado / dt); itTime < NITMAX; itTime++)
	{
		const time_t tSimple = time(nullptr);

		// SET INLET PERTUBATION
		{
			for (int k = 1; k < NZ - 1; k++)
			{
				field.cUVelField[0 + (k * NX)] = ufreestream +
					(sin(tempoPassado + (k * dd)) * ufreestream * 0.01);

				field.cUVelField[1 + (k * NX)] = ufreestream +
					(sin(tempoPassado + (k * dd)) * ufreestream * 0.01);
			}
		}

		for (int itSimple = 0; itSimple < NITMAXSIMPLE; itSimple++)
		{
			timer.Tick(tidSimple);

			if (GetAsyncKeyState('H') & 0x8000)
			{
				showText = false;
				std::cout << "@@@@@ CHANGED ALL TEXT TO HIDE (HOLD H TO HIDE - S TO SHOW) @@@@@\n";
			}
			else if (GetAsyncKeyState('S') & 0x8000)
			{
				showText = true;
				std::cout << "@@@@@ CHANGED ALL TEXT TO SHOW (HOLD H TO HIDE - S TO SHOW) @@@@@\n";
			}

#ifdef _DDEBUG
			double resSimplePosMom = 0.0;
			for (int k = 2; k < field.cNZ - 1; k++)
				for (int i = 2; i < field.cNX - 1; i++)
				{
					resSimplePosMom = fmax(
						fabs(
						(
							field.cUVelField[(j * NZ * NX) + (k * NX) + i] - field.cUVelField[(j * NZ * NX) + (k * NX) + i + 1] +
							field.cWVelField[(j * NZ * NX) + (k * NX) + i] - field.cWVelField[(j * NZ * NX) + ((k + 1) * NX) + i]
							)
						), resSimplePosMom);
				}
				

			std::cout << "Residual simple pre momentum eqts. " << resSimplePosMom << std::endl;
#endif

			//ulastIField = field.cUVelField;
			//wlastIField = field.cWVelField;
			//pLastIField = field.cPressureField;

			memcpy(ulastIField.data(), field.cUVelField.data(), field.cUVelField.size() * sizeof(double));
			memcpy(wlastIField.data(), field.cWVelField.data(), field.cWVelField.size() * sizeof(double));
			memcpy(pLastIField.data(), field.cPressureField.data(), field.cPressureField.size() * sizeof(double));

			time_t t1 = time(nullptr);

			// Calculate field hat values
			timer.Tick(tidCalcHatValues);

			field.CalcUHatValuesCN(u0Field, dt);
			field.CalcWHatValuesCN(w0Field, dt);

			timer.Tock(tidCalcHatValues);

			// Solve pressure equation
			timer.Tick(tidPress);

			timer.Tick(tidCreatPress);
			int nn = field.CreatePressureLSCSRFI(ptr, col, val, rhs);
			timer.Tock(tidCreatPress);

			double errPressLS = 0.0;
			int nitPressLS = 0;
			// solve LS Pressure eqt

			if (nItPressure > maxNItPressureToReCondition)
			{
				linearSystemSolver.PrecondtionCRS(nn, ptr, col, val, maxContinuityResidual * 1e-3);
			}

			linearSystemSolver.SolvePreconditionedCRS(std::move(field.cPressureField), nn, ptr, col,
				val, rhs, errPressLS, nitPressLS);

			//linearSystemSolver.SolveSparseCRS(std::move(field.cPressureField), nn, ptr, col, val, rhs, errPressLS, nitPressLS);

			nItPressure = nitPressLS;


			if (isnan(errPressLS))
			{
				std::cout << "Erro ls Pressure - nao foi possivel resolver ls\n";
				int x = 0;
				std::cin >> x;
			}

			timer.Tock(tidPress);

			if (showText)
			{
				std::cout << "Solved Pressure LS it " << nitPressLS << " residual " << errPressLS << std::endl;
			}

#ifdef _DEBUG
			field.SaveIKCutFieldToCSVFile(NY / 2, "OUT_DEBUG/posmomlong");
			field.SaveJKCutFieldToCSVFile(NX - 2, "OUT_DEBUG/posmomtrans");
#endif   



			timer.Tick(tidCorrsResiduals);

			// CALCULATE TRUE VELOCITIES
			// u
			for (int K = 1; K < NZ - 1; K++)
			{
				for (int i = 2; i < NX - 1; i++)
				{
					if (i >= field.ciCyInit && i < field.ciCyEnd + 1 &&
						K >= field.ckCyInit && K < field.ckCyEnd)
					{
						// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO
					}
					else
					{
						field.cUVelField[(K * NX) + i] = field.cUHatVelField[(K * NX) + i]
							+ ((field.cPressureField[(K * NX) + i - 1] - field.cPressureField[(K * NX) + i])
								*(dd / field.caijkU[(K * NX) + i]));

						/*(dd * dd / field.caijkU[(K * NX) + i]) *
						(field.cPressureCorrField[(K * NX) + i - 1] - field.cPressureCorrField[(K * NX) + i]);
						*/
					}
				}
			}

			//w
			for (int K = 2; K < NZ - 1; K++)
			{
				for (int i = 1; i < NX - 1; i++)
				{
					if (i >= field.ciCyInit && i < field.ciCyEnd &&
						K >= field.ckCyInit && K < field.ckCyEnd + 1)
					{
						// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO
					}
					else
					{
						//field.cWVelField[(J * NZ * NX) + (K * NX) + i] += (dd * dd / field.caijkW[(J * NZ * NX) + (K * NX) + i]) *
						//	(field.cPressureCorrField[(J * NZ * NX) + ((K - 1) * NX) + i] - field.cPressureCorrField[(J * NZ * NX) + (K * NX) + i]);

						field.cWVelField[(K * NX) + i] = field.cWHatVelField[(K * NX) + i]
							+ ((field.cPressureField[((K - 1) * NX) + i] - field.cPressureField[(K * NX) + i])
								* (dd / field.caijkW[(K * NX) + i]));
					}
				}
			}

			double ures = 0.0;
			double vres = 0.0;
			double wres = 0.0;
			double pres = 0.0;

			double sumCFL = 0.0;

			double CFLmax = 0.0;

			// calc residuals velocities
			for (int k = 0; k < NZ; k++)
				for (int i = 0; i < NX; i++)
				{
					const double uResidual = fabs(field.cUVelField[(k * NX) + i] - ulastIField[(k * NX) + i]);

					if (uResidual > ures)
					{
						ures = uResidual;
					}
					
					const double wResidual = fabs(field.cWVelField[(k * NX) + i] - wlastIField[(k * NX) + i]);

					if (wResidual > wres)
					{
						wres = wResidual;
					}

					const double cfl = (fabs(field.cUVelField[(k * NX) + i])
						+ fabs(field.cWVelField[(k * NX) + i])) * dt / dd;

					sumCFL += cfl;

					if (cfl > CFLmax)
						CFLmax = cfl;
				}

			const double CFLMean = sumCFL / (field.cNX * field.cNZ);

			int iPres = 0;
			int jPres = 0;
			int kPres = 0;
			double pij = 0.0;
			double presAcum = 0.0;

			for (int k = 2; k < NZ - 2; k++)
				for (int i = 2; i < NX - 2; i++)
					{
						const double residual = fabs((field.cPressureField[(k * NX) + i] -
							pLastIField[(k * NX) + i]));

						if (residual > pres)
						{
							pres = residual;
							iPres = i;
							kPres = k;
							pij = field.cPressureField[(k * NX) + i];
						}
					}

			for (int k = 2; k < NZ - 2; k++)
				for (int i = 2; i < NX - 2; i++)
					{
						const double residual = fabs((field.cPressureField[(k * NX) + i] -
							pLastIField[(k * NX) + i]));

						presAcum += residual;
					}

			double resSimple = 0.0;
			int iResSimple = 0;
			int jResSimple = 0;
			int kResSimple = 0;

			// SET/EXTRAPOLATE/ENFORCES OUTLET B.C.
			{
				for (int64_t k = 1; k < NZ - 1; k++)
				{
					const int64_t i = NX - 2;

					field.cUVelField[(k * NX) + i + 1] = field.cUVelField[(k * NX) + i] +
						field.cWVelField[(k * NX) + i] - field.cWVelField[((k + 1) * NX) + i];
				}

				field.ExtrapolateBackwardWJk(NX - 1);
			}

			for (int k = 2; k < field.cNZ - 1; k++)
				for (int i = 2; i < field.cNX - 1; i++)
				{
					const double residual = fabs(
						field.cUVelField[(k * NX) + i] - field.cUVelField[(k * NX) + i + 1] +
						field.cWVelField[(k * NX) + i] - field.cWVelField[((k + 1) * NX) + i]
					) * field.cRHO * field.cdA;

					if (residual > resSimple)
					{
						resSimple = residual;

						iResSimple = i;
						kResSimple = k;
					}
				}

			timer.Tock(tidCorrsResiduals);

			if (showText)
			{
				std::cout << "##########################################\n";
				std::cout << "##########################################\n";

				std::cout << "resSimple: " << resSimple << " ures "
					<< ures << " vres " << vres << " wres " << wres << " pres " << pres << std::endl;
				std::cout << "maior mudanca em pres i " << iPres << " j " << jPres << " k " << kPres << " pij " << pij << std::endl;
				std::cout << "pres acumulado " << presAcum << std::endl;
				std::cout << "maior mudanca em res simple i " << iResSimple << " j " << jResSimple << " k " << kResSimple << std::endl;
				std::cout << "CFL medio " << CFLMean << " maior CFL " << CFLmax << std::endl;
			}

#ifdef _DEBUG
			field.SaveIKCutFieldToCSVFile(NY / 2, "OUT_DEBUG/poscorrlong");
			field.SaveJKCutFieldToCSVFile(NX - 2, "OUT_DEBUG/poscorrtrans");
#endif

			if ((resSimple < maxContinuityResidual && pres < maxpResidual) && ures < maxuResidual
				&& wres < maxwResidual)
			{
				std::cout << "Simple converged; time: " << tempoPassado << std::endl;
				std::cout << "Tempo para convergir simple: " << time(nullptr) - tSimple << " \n";

				std::cout << "##########################################\n";
				std::cout << "##########################################\n";

				if (((itTime % nItToSave) == 0) || (tempoPassado > 12.0 && (itTime % (nItToSave / 6)) == 0))
				{
					field.SaveField("output/field" + std::to_string(tempoPassado));
					field.SaveIKCutFieldToCSVFile("turbFlatPlate/corteLong" + std::to_string(tempoPassado));
				}

				break;
			}
			if (showText)
			{
				std::cout << "Last simple iteration spent " << time(nullptr) - t1 << " seconds\n";
			}

			timer.Tock(tidSimple);

			if (showText)
			{
				timer.WriteToCoutAllTickTock();
			}

		}

		u0Field = field.cUVelField;
		w0Field = field.cWVelField;

		tempoPassado += dt;
		if (showText)
		{
			std::cout << "Tempo aprox restante para conclusao : " << (NITMAX - itTime) * (time(nullptr) - tSimple) << " segs  ou "
				<< ((NITMAX - itTime) * (time(nullptr) - tSimple)) / 3600 << "horas" << std::endl;
		}

	}

	return 0;
}

template<class T>
T GetValue(const char* message)
{
	T v;

	std::cout << message;
	std::cin >> v;
	std::cout << "\n";

	return v;
}