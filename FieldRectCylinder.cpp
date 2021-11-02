#include "FieldRectCylinder.h"
#include <fstream>
#include <iostream>
#include <ctype.h>

FieldRectCylinder::FieldRectCylinder(double L, double H, double W, double l, double deltaSize, double rho, double mu, double Ufarfield,
	double Vfarfield, double Wfarfield)
	:
	ctidPushback(cTimer.SetTimer("Push back time (U mom)")),
	ctidCalcCoefs(cTimer.SetTimer("Calc coefs time (U mom)")),
	cdd(deltaSize),
	cdA(cdd* cdd),
	cdV(cdd* cdd* cdd),
	cNY(int64_t(W / cdd) + 2),
	cNX(int64_t(L / cdd) + 2),
	cNZ(int64_t(H / cdd) + 2),
	cNXCylinder(int64_t(0.5 + (l / cdd))),
	cNYCylinder(cNY),
	cNZCylinder(cNXCylinder),
	ciCyInit((cNXCylinder * 8) - cNXCylinder),
	cjCyInit(0),
	ckCyInit((cNZ/2) - (cNZCylinder/2)),
	ciCyEnd(cNXCylinder * 8),
	cjCyEnd(cNY),
	ckCyEnd((cNZ/2) + (cNZCylinder/2)),
	cRHO(rho),
	cMU(mu),
	cD(cdA * mu / cdd)
{
	cPressureField.resize(cNY * cNZ * cNX, 0.0);

	cPressureCorrField.resize(cNY * cNZ * cNX, 0.0);

	cUVelField.resize(cNY * cNZ * cNX, Ufarfield);
	cVVelField.resize(cNY * cNZ * cNX, Vfarfield);	
	cWVelField.resize(cNY * cNZ * cNX, Wfarfield);
	
	cSpUField.resize(cNY * cNZ * cNX, 0.0);
	cSpVField.resize(cNY * cNZ * cNX, 0.0);
	cSpWField.resize(cNY * cNZ * cNX, 0.0);
	
	cUHatVelField = cPressureCorrField;
	cVHatVelField = cPressureCorrField;
	cWHatVelField = cPressureCorrField;

	caijkU = cPressureCorrField;
	caijkV = cPressureCorrField;
	caijkW = cPressureCorrField;
}

FieldRectCylinder::~FieldRectCylinder()
{}

FieldRectCylinder::Intervalo::Intervalo(int64_t v1, int64_t v2)
	:
	valor1(v1),
	valor2(v2)
{}

int64_t FieldRectCylinder::Intervalo::Tamanho() const
{
	return valor2 - valor1;
}

void FieldRectCylinder::SetBCFlatPlate()
{
	// SET INITIAL VALUES
	SetUValueInterval({ ciCyInit, ciCyEnd+1 }, { cjCyInit, cjCyEnd }, { ckCyInit, ckCyEnd }, 0.0);
	SetVValueInterval({ ciCyInit, ciCyEnd }, { cjCyInit, cjCyEnd }, { ckCyInit, ckCyEnd }, 0.0);
	SetWValueInterval({ ciCyInit, ciCyEnd }, { cjCyInit, cjCyEnd }, { ckCyInit, ckCyEnd + 1 }, 0.0);
	SetPValueInterval({ ciCyInit, ciCyEnd }, { cjCyInit, cjCyEnd }, { ckCyInit, ckCyEnd }, 0.0);

	// INLET; OUTLET TUDO OK
	// VALORES PADROES

	// WALL (1)
	SetSpVValueInterval({ ciCyInit-1, ciCyInit }, { cjCyInit, cjCyEnd }, { ckCyInit, ckCyEnd }, -cMU * cdA / (cdd / 2.0));
	SetSpWValueInterval({ ciCyInit-1, ciCyInit }, { cjCyInit, cjCyEnd }, { ckCyInit, ckCyEnd + 1 }, -cMU * cdA / (cdd / 2.0));

	// WALL (2)
	SetSpVValueInterval({ ciCyEnd, ciCyEnd + 1 }, { cjCyInit, cjCyEnd }, { ckCyInit, ckCyEnd }, -cMU * cdA / (cdd / 2.0));
	SetSpWValueInterval({ ciCyEnd, ciCyEnd + 1 }, { cjCyInit, cjCyEnd }, { ckCyInit, ckCyEnd + 1 }, -cMU * cdA / (cdd / 2.0));

	// WALL (4)
	SetSpUValueInterval({ ciCyInit, ciCyEnd + 1 }, { cjCyInit, cjCyEnd }, { ckCyInit - 1, ckCyInit }, -cMU * cdA / (cdd / 2.0));
	SetSpVValueInterval({ ciCyInit, ciCyEnd }, { cjCyInit, cjCyEnd }, { ckCyInit - 1, ckCyInit }, -cMU * cdA / (cdd / 2.0));

	// WALL (3)
	SetSpUValueInterval({ ciCyInit, ciCyEnd + 1 }, { cjCyInit, cjCyEnd }, { ckCyEnd, ckCyEnd + 1 }, -cMU * cdA / (cdd / 2.0));
	SetSpVValueInterval({ ciCyInit, ciCyEnd }, { cjCyInit, cjCyEnd }, { ckCyEnd, ckCyEnd + 1 }, -cMU * cdA / (cdd / 2.0));
}

int FieldRectCylinder::CreateUMomentumLSCSR(const doubleField3D& uField0, double dt,
	std::vector<int>& ptr, std::vector<int>& col, std::vector<double>& val, std::vector<double>& rhs)
{
	double totalTimeCalcs = 0.0;
	double totalTimePushs = 0.0;

	const int n2 = static_cast<int64_t>(cNX * cNY * cNZ);        // Number of points in the grid.
		
	ptr.clear();
	if (ptr.capacity() != n2 + 1)
		ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	if (col.capacity() != n2 * 7)
		col.reserve(n2 * 7);

	val.clear();
	if (val.capacity() != n2 * 7)
		val.reserve(n2 * 7);

	rhs.clear();
	rhs.resize(n2);

	for (int j = 0; j < cNY; j++) 
	{
		for (int k = 0; k < cNZ; k++)
		{
			for (int i = 0; i < cNX; i++) 
			{
				const int64_t index = static_cast<int64_t>(i + (k * cNX) + (j * cNX * cNZ));

				if (i == 1 || i == 0 || i == cNX - 1 || j == 0 || j == cNY - 1 || k == 0 || k == cNZ - 1)
				{
					// Boundary point. Use Dirichlet condition. (seta o valor da variavel diretamente)

					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = cUVelField[(j * cNZ * cNX) + (k * cNX) + i];
				}
				else if (i >= ciCyInit && i < ciCyEnd + 1 && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA U = 0 )
					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else 
				{
					//if (isnan(aijkU.aijk) || isnan(aijkU.source))
					//{
					//	int x = 0;
					//}

					// Interior point. Use 5-point finite difference stencil.

					// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]

					cTimer.Tick(ctidPushback);

					// jki
					col.push_back(index + static_cast<int64_t>((0 * cNX) + (0 * cNX * cNZ)));
					val.push_back(NAN);

					// jkmi
					col.push_back(index + static_cast<int64_t>((-1 * cNX) + (0 * cNX * cNZ)));
					val.push_back(-NAN);

					// jkim
					col.push_back(index - static_cast<int64_t>(1 + (0 * cNX) + (0 * cNX * cNZ)));
					val.push_back(-NAN);					

					// jkip
					col.push_back(index + 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-NAN);

					// jkpi
					col.push_back(index + static_cast<int64_t>((1 * cNX) + (0 * cNX * cNZ)));
					val.push_back(-NAN);

					// jmki
					col.push_back(index + static_cast<int64_t>((0 * cNX) + (-1 * cNX * cNZ)));
					val.push_back(-NAN);

					// jpki
					col.push_back(index + static_cast<int64_t>((0 * cNX) + (1 * cNX * cNZ)));
					val.push_back(-NAN);

					rhs[index] = NAN;

					cTimer.Tock(ctidPushback);

					totalTimePushs += cTimer.GetLastTickTock(ctidPushback);
				}

				ptr.push_back(static_cast<int64_t>(col.size()));
			}
		}		
	}

	cTimer.Tick(ctidCalcCoefs);

#pragma omp parallel for
	for (int j = 0; j < cNY; j++)
	{
		for (int k = 0; k < cNZ; k++)
		{
			for (int i = 0; i < cNX; i++) 
			{
				const int64_t index = static_cast<int64_t>(i + (k * cNX) + (j * cNX * cNZ));

				if (i == 1 || i == 0 || i == cNX - 1 || j == 0 || j == cNY - 1 || k == 0 || k == cNZ - 1)
				{
					// do nothing, already did it before;
				}
				else if (i >= ciCyInit && i < ciCyEnd + 1 && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					// do nothing, already did it before;
				}
				else
				{
					CoefsData aijkU;

					const double Fw = GetFwInternalUMomentum(i, j, k);
					const double Fe = GetFeInternalUMomentum(i, j, k);
					const double Fs = GetFsInternalUMomentum(i, j, k);
					const double Fn = GetFnInternalUMomentum(i, j, k);
					const double Fb = GetFbInternalUMomentum(i, j, k);
					const double Ft = GetFtInternalUMomentum(i, j, k);

					const double deltaf = cdA * (Fe - Fw + Fn - Fs + Ft - Fb);

					aijkU.aimjk = cD + (cdA*fmax(Fw, 0.0));
					aijkU.aipjk = cD + (cdA*fmax(-Fe, 0.0));

					aijkU.aijmk = cD + (cdA*fmax(Fs, 0.0));
					aijkU.aijpk = cD + (cdA*fmax(-Fn, 0.0));

					aijkU.aijkm = cD + (cdA*fmax(Fb, 0.0));
					aijkU.aijkp = cD + (cdA*fmax(-Ft, 0.0));

					const double alfaw = (Fw > 0.0) ? 1.0 : 0.0;
					const double alfae = (Fe > 0.0) ? 1.0 : 0.0;
					const double alfas = (Fs > 0.0) ? 1.0 : 0.0;
					const double alfan = (Fn > 0.0) ? 1.0 : 0.0;
					const double alfat = (Ft > 0.0) ? 1.0 : 0.0;
					const double alfab = (Fb > 0.0) ? 1.0 : 0.0;

					const double reneg = ren(i, j, k, cUVelField);
					const double rwneg = rwn(i, j, k, cUVelField);
					const double rnneg = rnn(i, j, k, cUVelField);
					const double rsneg = rsn(i, j, k, cUVelField);
					const double rtneg = rtn(i, j, k, cUVelField);
					const double rbneg = rbn(i, j, k, cUVelField);

					const double repos = rep(i, j, k, cUVelField);
					const double rwpos = rwp(i, j, k, true, cUVelField);
					const double rnpos = rnp(i, j, k, cUVelField);
					const double rspos = rsp(i, j, k, false, cUVelField);
					const double rtpos = rtp(i, j, k, cUVelField);
					const double rbpos = rbp(i, j, k, false, cUVelField);

					const double phiE = cUVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phiP = cUVelField[index];
					const double phiW = cUVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phiN = cUVelField[((j + 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiS = cUVelField[((j - 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiT = cUVelField[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
					const double phiB = cUVelField[(j * cNZ * cNX) + ((k - 1) * cNX) + i];

					const double SuDc = cdA * 0.5*(
						(Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
						+ (Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
						+ (Fn * (((1.0 - alfan) * Psir(rnneg)) - (alfan * Psir(rnpos))) * (phiN - phiP))
						+ (Fs * (-((1.0 - alfas) * Psir(rsneg)) + (alfas * Psir(rspos))) * (phiP - phiS))
						+ (Ft * (((1.0 - alfat) * Psir(rtneg)) - (alfat * Psir(rtpos))) * (phiT - phiP))
						+ (Fb * (-((1.0 - alfab) * Psir(rbneg)) + (alfab * Psir(rbpos))) * (phiP - phiB))
						);


					// SERIA UM ERRO SETAR ESSES COEFS A ZERO (Os valores da fronteira sao utilizados)
					//if (i == 2)
					//	aijkU[J][K][i].aimjk = 0.0;
					//if (i == NXx - 3)
					//	aijkU[J][K][i].aipjk = 0.0;
					// ^^ SERIA UM ERRO SETAR ESSES COEFS A ZERO (Os valores da fronteira sao utilizados)
					// FIM
					// UNICA PAREDE ( k = 1 ) 
					//if (k == 1)
					//	aijkU.aijkm = 0.0;
					//else if (k == cNZ - 2)
					//	aijkU.aijkp = 0.0;
					//if (j == 1)
					//	aijkU.aijmk = 0.0;
					//else if (j == cNY - 2)
					//	aijkU.aijpk = 0.0;

					// PAREDES
					// WALL (1)
					// DESNECESSARIO (só resolver normalmente)

					// WALL (2)
					// DESNECESSARIO (só resolver normalmente)

					// WALL (4)
					if (k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd + 1)
					{
						aijkU.aijkp = 0.0;
					}
					// WALL (3)
					else if (k == ckCyEnd && i >= ciCyInit && i < ciCyEnd + 1)
					{
						aijkU.aijkm = 0.0;
					}					

					aijkU.aijk = aijkU.aimjk + aijkU.aipjk + aijkU.aijmk
						+ aijkU.aijpk + aijkU.aijkm + aijkU.aijkp + deltaf - cSpUField[index]
						+ (cRHO * cdV / dt);

					caijkU[(j * cNZ * cNX) + (k * cNX) + i] = aijkU.aijk;
										

					aijkU.source = SuDc +
						((cPressureField[(j * cNZ * cNX) + (k * cNX) + i - 1] - cPressureField[(j * cNZ * cNX) + (k * cNX) + i]) * cdA) +
						/*cSuUField[(j * cNZ * cNX) + (k * cNX) + i] +*/
						((cRHO * cdd * cdd * cdd / dt) * uField0[(j * cNZ * cNX) + (k * cNX) + i]);
					

					// jki
					val[ptr[index]] = aijkU.aijk;

					// jkmi
					val[ptr[index]+1] = -aijkU.aijkm;

					// jkim
					val[ptr[index] + 2] = -aijkU.aimjk;

					// jkip
					val[ptr[index] + 3]= -aijkU.aipjk;

					// jkpi
					val[ptr[index] + 4] = -aijkU.aijkp;

					// jmki
					val[ptr[index] + 5] = -aijkU.aijmk;

					// jpki
					val[ptr[index] + 6] = -aijkU.aijpk;

					rhs[index] = aijkU.source;									   			
				}
			}
		}
	}

	cTimer.Tock(ctidCalcCoefs);

	//cTimer.WriteToCoutTickTock(ctidCalcCoefs);
	//   
	//std::cout << "# Timer total time calcs " << totalTimeCalcs << " total time pushs " << totalTimePushs << std::endl;

	return n2;
}


int FieldRectCylinder::CreateVMomentumLSCSRParallel(const doubleField3D& vField0, double dt, std::vector<int>& ptr,
	std::vector<int>& col, std::vector<double>& val, std::vector<double>& rhs)
{
	const int n2 = cNX * cNY * cNZ; // Number of points in the grid.

	ptr.clear();
	if (ptr.capacity() != n2 + 1)
		ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	if (col.capacity() != n2 * 7)
		col.reserve(n2 * 7);

	val.clear();
	if (val.capacity() != n2 * 7)
		val.reserve(n2 * 7);

	rhs.clear();
	rhs.resize(n2);

	for (int j = 0, index = 0; j < cNY; ++j)
	{
		for (int k = 0; k < cNZ; k++)
		{
			for (int i = 0; i < cNX; i++, index++) 
			{

				if (i == 0 || i == cNX - 1 || j == 0 || j == 1 || j == cNY - 1 || k == 0 || k == cNZ - 1)
				{
					// Boundary point. Use Dirichlet condition. (seta o valor da variavel diretamente)

					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = cVVelField[(j * cNZ * cNX) + (k * cNX) + i];
				}
				else if (i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd + 1 && k >= ckCyInit && k < ckCyEnd)
				{
					// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA V = 0 )
					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else
				{					
					// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]

					// jki
					col.push_back(index + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(NAN);

					// jkmi
					col.push_back(index + (-1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-NAN);

					// jkim
					col.push_back(index - 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-NAN);					

					// jkip
					col.push_back(index + 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-NAN);

					// jkpi
					col.push_back(index + (1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-NAN);

					// jmki
					col.push_back(index + (0 * cNX) + (-1 * cNX * cNZ));
					val.push_back(-NAN);

					// jpki
					col.push_back(index + (0 * cNX) + (1 * cNX * cNZ));
					val.push_back(-NAN);

					//rhs[index] = NAN;

				}

				ptr.push_back(col.size());
			}
		}
	}

#pragma omp parallel for
	for (int j = 0; j < cNY; ++j)
	{
		for (int k = 0; k < cNZ; k++)
		{
			for (int i = 0; i < cNX; i++) 
			{
				const int index = i + (k * cNX) + (j * cNX * cNZ);

				if (i == 0 || i == cNX - 1 || j == 0 || j == 1 || j == cNY - 1 || k == 0 || k == cNZ - 1)
				{
					// do nothing, already did it
				}
				else if (i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd + 1 && k >= ckCyInit && k < ckCyEnd)
				{
					// do nothing, already did it before;
				}
				else
				{
					CoefsData aijkV;

					const double Fw = GetFwInternalVMomentum(i, j, k) * cdA;
					const double Fe = GetFeInternalVMomentum(i, j, k) * cdA;
					const double Fs = GetFsInternalVMomentum(i, j, k) * cdA;
					const double Fn = GetFnInternalVMomentum(i, j, k) * cdA;
					const double Fb = GetFbInternalVMomentum(i, j, k) * cdA;
					const double Ft = GetFtInternalVMomentum(i, j, k) * cdA;

					const double deltaf = Fe - Fw + Fn - Fs + Ft - Fb;

					aijkV.aimjk = cD + fmax(Fw, 0.0);
					aijkV.aipjk = cD + fmax(-Fe, 0.0);

					aijkV.aijmk = cD + fmax(Fs, 0.0);
					aijkV.aijpk = cD + fmax(-Fn, 0.0);

					aijkV.aijkm = cD + fmax(Fb, 0.0);
					aijkV.aijkp = cD + fmax(-Ft, 0.0);

					const double alfaw = (Fw > 0) ? 1.0 : 0.0;
					const double alfae = (Fe > 0) ? 1.0 : 0.0;
					const double alfas = (Fs > 0) ? 1.0 : 0.0;
					const double alfan = (Fn > 0) ? 1.0 : 0.0;
					const double alfab = (Fb > 0) ? 1.0 : 0.0;
					const double alfat = (Ft > 0) ? 1.0 : 0.0;

					const double reneg = ren(i, j, k, cVVelField);
					const double rwneg = rwn(i, j, k, cVVelField);
					const double rnneg = rnn(i, j, k, cVVelField);
					const double rsneg = rsn(i, j, k, cVVelField);
					const double rtneg = rtn(i, j, k, cVVelField);
					const double rbneg = rbn(i, j, k, cVVelField);

					const double repos = rep(i, j, k, cVVelField);
					const double rwpos = rwp(i, j, k, false, cVVelField);
					const double rnpos = rnp(i, j, k, cVVelField);
					const double rspos = rsp(i, j, k, true, cVVelField);
					const double rtpos = rtp(i, j, k, cVVelField);
					const double rbpos = rbp(i, j, k, false, cVVelField);

					const double phiE = cVVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phiP = cVVelField[(j * cNZ * cNX) + (k * cNX) + i];
					const double phiW = cVVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phiN = cVVelField[((j + 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiS = cVVelField[((j - 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiT = cVVelField[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
					const double phiB = cVVelField[(j * cNZ * cNX) + ((k - 1) * cNX) + i];


					const double SuDc = (0.5 * Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
						+ (0.5 * Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
						+ (0.5 * Fn * (((1.0 - alfan) * Psir(rnneg)) - (alfan * Psir(rnpos))) * (phiN - phiP))
						+ (0.5 * Fs * (-((1.0 - alfas) * Psir(rsneg)) + (alfas * Psir(rspos))) * (phiP - phiS))
						+ (0.5 * Ft * (((1.0 - alfat) * Psir(rtneg)) - (alfat * Psir(rtpos))) * (phiT - phiP))
						+ (0.5 * Fb * (-((1.0 - alfab) * Psir(rbneg)) + (alfab * Psir(rbpos))) * (phiP - phiB));
					
					// PAREDES
					// WALL (1)
					if (i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd)
					{
						aijkV.aipjk = 0.0;
					}
					// WALL (2)
					else if (i == ciCyEnd && k >= ckCyInit && k < ckCyEnd)
					{
						aijkV.aimjk = 0.0;
					}

					// WALL (4)
					if (k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd)
					{
						aijkV.aijkp = 0.0;
					}
					// WALL (3)
					else if (k == ckCyEnd && i >= ciCyInit && i < ciCyEnd)
					{
						aijkV.aijkm = 0.0;
					}

					aijkV.aijk = aijkV.aimjk + aijkV.aipjk + aijkV.aijmk
						+ aijkV.aijpk + aijkV.aijkm + aijkV.aijkp + deltaf - cSpVField[(j * cNZ * cNX) + (k * cNX) + i]
						+ (cRHO * cdV / dt);

					caijkV[(j * cNZ * cNX) + (k * cNX) + i] = aijkV.aijk;

					aijkV.source = SuDc +
						((cPressureField[((j - 1) * cNZ *cNX) + (k*cNX) + i] - cPressureField[(j * cNZ * cNX) + (k * cNX) + i])*cdA) +
						/*cSuVField[(j * cNZ * cNX) + (k * cNX) + i] + */
						((cRHO * cdd * cdd * cdd / dt) * vField0[(j * cNZ * cNX) + (k * cNX) + i]);

					//if (isnan(aijkV.aijk) || isnan(aijkV.source))
					//{
					//	int x = 0;
					//}

					// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]


					// jki
					val[ptr[index] + 0] = aijkV.aijk;

					// jkmi
					val[ptr[index] + 1] = -aijkV.aijkm;

					// jkim
					val[ptr[index] + 2] = -aijkV.aimjk;

					// jkip
					val[ptr[index] + 3] = -aijkV.aipjk;

					// jkpi
					val[ptr[index] + 4] = -aijkV.aijkp;

					// jmki
					val[ptr[index] + 5] = -aijkV.aijmk;

					// jpki
					val[ptr[index] + 6] = -aijkV.aijpk;

					rhs[index] = aijkV.source;
				}
			}
		}
	}

	return n2;
}

int FieldRectCylinder::CreateWMomentumLSCSRParallel(const doubleField3D& wField0, double dt, std::vector<int>& ptr,
	std::vector<int>& col, std::vector<double>& val, std::vector<double>& rhs)
{
	const int n2 = cNX * cNY * cNZ; // Number of points in the grid.

	ptr.clear();
	if (ptr.capacity() != n2 + 1)
		ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	if (col.capacity() != n2 * 7)
		col.reserve(n2 * 7);

	val.clear();
	if (val.capacity() != n2 * 7)
		val.reserve(n2 * 7);

	rhs.clear();
	rhs.resize(n2);

	for (int j = 0, index = 0; j < cNY; ++j)
	{
		for (int k = 0; k < cNZ; k++)
		{
			for (int i = 0; i < cNX; i++, index++) {

				if (i == 0 || i == cNX - 1 || j == 0 || j == cNY - 1 || k == 0 || k == 1 || k == cNZ - 1)
				{
					// Boundary point. Use Dirichlet condition. (seta o valor da variavel diretamente)

					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = cWVelField[(j * cNZ * cNX) + (k * cNX) + i];
				}
				else if (i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd + 1)
				{
					// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA V = 0 )
					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else
				{					
					// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]

					// jki
					col.push_back(index + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(NAN);

					// jkmi
					col.push_back(index + (-1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-NAN);

					// jkim
					col.push_back(index - 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-NAN);					

					// jkip
					col.push_back(index + 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-NAN);

					// jkpi
					col.push_back(index + (1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-NAN);

					// jmki
					col.push_back(index + (0 * cNX) + (-1 * cNX * cNZ));
					val.push_back(-NAN);

					// jpki
					col.push_back(index + (0 * cNX) + (1 * cNX * cNZ));
					val.push_back(-NAN);

					//rhs[index] = aijkW.source;
				}

				ptr.push_back(col.size());
			}
		}
	}

#pragma omp parallel for
	for (int j = 0; j < cNY; ++j)
	{
		for (int k = 0; k < cNZ; k++)
		{
			for (int i = 0; i < cNX; i++) 
			{
				const int index = i + (k * cNX) + (j * cNX * cNZ);

				if (i == 0 || i == cNX - 1 || j == 0 || j == cNY - 1 || k == 0 || k == 1 || k == cNZ - 1)
				{
					// do nothing, already did it
				}
				else if (i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd + 1)
				{
					// do nothing, already did it
				}
				else
				{
					CoefsData aijkW;

					const double Fw = GetFwInternalWMomentum(i, j, k) * cdA;
					const double Fe = GetFeInternalWMomentum(i, j, k) * cdA;
					const double Fs = GetFsInternalWMomentum(i, j, k) * cdA;
					const double Fn = GetFnInternalWMomentum(i, j, k) * cdA;
					const double Fb = GetFbInternalWMomentum(i, j, k) * cdA;
					const double Ft = GetFtInternalWMomentum(i, j, k) * cdA;

					const double deltaf = Fe - Fw + Fn - Fs + Ft - Fb;

					aijkW.aimjk = cD + fmax(Fw, 0.0);
					aijkW.aipjk = cD + fmax(-Fe, 0.0);

					aijkW.aijmk = cD + fmax(Fs, 0.0);
					aijkW.aijpk = cD + fmax(-Fn, 0.0);

					aijkW.aijkm = cD + fmax(Fb, 0.0);
					aijkW.aijkp = cD + fmax(-Ft, 0.0);

					const double alfaw = (Fw > 0) ? 1.0 : 0.0;
					const double alfae = (Fe > 0) ? 1.0 : 0.0;
					const double alfas = (Fs > 0) ? 1.0 : 0.0;
					const double alfan = (Fn > 0) ? 1.0 : 0.0;
					const double alfab = (Fb > 0) ? 1.0 : 0.0;
					const double alfat = (Ft > 0) ? 1.0 : 0.0;

					const double reneg = ren(i, j, k, cWVelField);
					const double rwneg = rwn(i, j, k, cWVelField);
					const double rnneg = rnn(i, j, k, cWVelField);
					const double rsneg = rsn(i, j, k, cWVelField);
					const double rtneg = rtn(i, j, k, cWVelField);
					const double rbneg = rbn(i, j, k, cWVelField);

					const double repos = rep(i, j, k, cWVelField);
					const double rwpos = rwp(i, j, k, false, cWVelField);
					const double rnpos = rnp(i, j, k, cWVelField);
					const double rspos = rsp(i, j, k, false, cWVelField);
					const double rtpos = rtp(i, j, k, cWVelField);
					const double rbpos = rbp(i, j, k, true, cWVelField);

					const double phiE = cWVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phiP = cWVelField[(j * cNZ * cNX) + (k * cNX) + i];
					const double phiW = cWVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phiN = cWVelField[((j + 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiS = cWVelField[((j - 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiT = cWVelField[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
					const double phiB = cWVelField[(j * cNZ * cNX) + ((k - 1) * cNX) + i];

					const double SuDc = (0.5 * Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
						+ (0.5 * Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
						+ (0.5 * Fn * (((1.0 - alfan) * Psir(rnneg)) - (alfan * Psir(rnpos))) * (phiN - phiP))
						+ (0.5 * Fs * (-((1.0 - alfas) * Psir(rsneg)) + (alfas * Psir(rspos))) * (phiP - phiS))
						+ (0.5 * Ft * (((1.0 - alfat) * Psir(rtneg)) - (alfat * Psir(rtpos))) * (phiT - phiP))
						+ (0.5 * Fb * (-((1.0 - alfab) * Psir(rbneg)) + (alfab * Psir(rbpos))) * (phiP - phiB));

					
					// PAREDES
					// WALL (1)
					if (i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd + 1)
					{
						aijkW.aipjk = 0.0;
					}
					// WALL (2)
					else if (i == ciCyEnd && k >= ckCyInit && k < ckCyEnd + 1)
					{
						aijkW.aimjk = 0.0;
					}

					// WALL (4)
					// RESOLVER NORMALMENTE

					// WALL (3)
					// RESOLVER NORMALMENTE

					aijkW.aijk = aijkW.aimjk + aijkW.aipjk + aijkW.aijmk
						+ aijkW.aijpk + aijkW.aijkm + aijkW.aijkp + deltaf - cSpWField[(j * cNZ * cNX) + (k * cNX) + i]
						+ (cRHO * cdV / dt);

					caijkW[(j * cNZ * cNX) + (k * cNX) + i] = aijkW.aijk;

					aijkW.source = SuDc + ((cPressureField[(j * cNZ * cNX) + ((k - 1) * cNX) + i] - cPressureField[(j * cNZ * cNX) + (k * cNX) + i]) * cdA) +
						/*cSuWField[(j * cNZ * cNX) + (k * cNX) + i] +*/ ((cRHO * cdd * cdd * cdd / dt) * wField0[(j * cNZ * cNX) + (k * cNX) + i]);

					// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]

					// jki
					val[ptr[index] + 0] = aijkW.aijk;

					// jkmi
					val[ptr[index] + 1] = -aijkW.aijkm;

					// jkim
					val[ptr[index] + 2] = -aijkW.aimjk;
					
					// jkip
					val[ptr[index] + 3] = -aijkW.aipjk;

					// jkpi
					val[ptr[index] + 4] = -aijkW.aijkp;

					// jmki
					val[ptr[index] + 5] = -aijkW.aijmk;

					// jpki
					val[ptr[index] + 6] = -aijkW.aijpk;

					rhs[index] = aijkW.source;
				}
			}
		}
	}
	return n2;
}
int FieldRectCylinder::CreatePressureCorrectionLSCSR(std::vector<int>& ptr, std::vector<int>& col, std::vector<double>& val,
	std::vector<double>& rhs) const
{
	const int n2 = cNX * cNY * cNZ; // Number of points in the grid.
	/*
	if (ptr.size() != n2 + 1)
	{
		ptr.clear();
		ptr.reserve(n2 + 1);
		ptr.push_back(0);
	}

	if (col.size() != n2 * 7)
	{
		col.clear();
		col.reserve(n2 * 7); // sete coeficientes
	}

	if (val.size() != n2 * 7)
	{
		val.clear();
		val.reserve(n2 * 7);  // sete coeficientes
	}

	if (rhs.size() != n2)
	{
		rhs.resize(n2);
	}
	*/

	ptr.clear();
	ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	col.reserve(n2 * 7);

	val.clear();
	val.reserve(n2 * 7);

	rhs.clear();
	rhs.resize(n2);
	for (int j = 0, index = 0; j < cNY; ++j)
	{
		for (int k = 0; k < cNZ; k++)
		{
			for (int i = 0; i < cNX; i++, index++) {

				if (i == 0 || i == cNX - 1 || j == 0 || j == cNY - 1 || k == 0 || k == cNZ - 1)
				{
					// Boundary point. Use Dirichlet condition. (seta o valor da variavel diretamente)

					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else if (i == cNX - 2)
				{
					// NO CORRECTION -> already known
					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else if (i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA PCORR = 0 )
					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else
				{
					CoefsData aijkPCorr;

					aijkPCorr.aimjk = cRHO * cdA * cdA / caijkU[(j * cNZ * cNX) + (k * cNX) + i];
					aijkPCorr.aipjk = cRHO * cdA * cdA / caijkU[(j * cNZ * cNX) + (k * cNX) + i + 1];

					aijkPCorr.aijmk = cRHO * cdA * cdA / caijkV[(j * cNZ * cNX) + (k * cNX) + i];
					aijkPCorr.aijpk = cRHO * cdA * cdA / caijkV[((j+1) * cNZ * cNX) + (k * cNX) + i];

					aijkPCorr.aijkm = cRHO * cdA * cdA / caijkW[(j * cNZ * cNX) + (k * cNX) + i];
					aijkPCorr.aijkp = cRHO * cdA * cdA / caijkW[(j * cNZ * cNX) + ((k+1) * cNX) + i];

					if (i == 1)
						aijkPCorr.aimjk = 0.0;
					else if (i == cNX - 2)
						aijkPCorr.aipjk = 0.0;
					if (j == 1)
						aijkPCorr.aijmk = 0.0;
					else if (j == cNY - 2)
						aijkPCorr.aijpk = 0.0;
					if (k == 1)
						aijkPCorr.aijkm = 0.0;
					else if (k == cNZ - 2)
						aijkPCorr.aijkp = 0.0;

					// PAREDES
					// WALL (1)
					if (i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd)
					{
						aijkPCorr.aipjk = 0.0;
					}
					// WALL (2)
					else if (i == ciCyEnd && k >= ckCyInit && k < ckCyEnd)
					{
						aijkPCorr.aimjk = 0.0;
					}

					// WALL (4)
					if (k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd)
					{
						aijkPCorr.aijkp = 0.0;
					}
					// WALL (3)
					else if (k == ckCyEnd && i >= ciCyInit && i < ciCyEnd)
					{
						aijkPCorr.aijkm = 0.0;
					}


					aijkPCorr.aijk = (aijkPCorr.aimjk + aijkPCorr.aipjk +
						aijkPCorr.aijmk + aijkPCorr.aijpk + aijkPCorr.aijkm +
						aijkPCorr.aijkp);

					aijkPCorr.source = cRHO * cdA * (cUVelField[(j * cNZ * cNX) + (k * cNX) + i] - cUVelField[(j * cNZ * cNX) + (k * cNX) + i + 1] +
						cVVelField[(j * cNZ * cNX) + (k * cNX) + i] - cVVelField[((j+1) * cNZ * cNX) + (k * cNX) + i] + 
						cWVelField[(j * cNZ * cNX) + (k * cNX) + i] - cWVelField[(j * cNZ * cNX) + ((k+1) * cNX) + i]);

					// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]
					
					// jki
					col.push_back(index + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(aijkPCorr.aijk);

					// jkmi
					col.push_back(index + (-1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkPCorr.aijkm);

					// jkim
					col.push_back(index - 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkPCorr.aimjk);
					
					// jkip
					col.push_back(index + 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkPCorr.aipjk);

					// jkpi
					col.push_back(index + (1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkPCorr.aijkp);

					// jmki
					col.push_back(index + (0 * cNX) + (-1 * cNX * cNZ));
					val.push_back(-aijkPCorr.aijmk);

					// jpki
					col.push_back(index + (0 * cNX) + (1 * cNX * cNZ));
					val.push_back(-aijkPCorr.aijpk);

					rhs[index] = aijkPCorr.source;
				}

				ptr.push_back(col.size());
			}
		}
	}

	return n2;
}

int FieldRectCylinder::CreatePressureLSCSRFI(std::vector<int>& ptr, std::vector<int>& col, std::vector<double>& val,
	std::vector<double>& rhs) const
{
	int n2 = cNX * cNY * cNZ; // total number of points in the grid.
	
	ptr.clear();
	ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	col.reserve(n2 * 7);

	val.clear();
	val.reserve(n2 * 7);

	rhs.clear();
	rhs.resize(n2);

	for (int j = 0; j < cNY; ++j)
	{
		for (int k = 0; k < cNZ; k++)
		{
			for (int i = 0; i < cNX; i++) {

				const int64_t index = i + (k * cNX) + (j * cNX * cNZ);
								
				if (i == 0 || i == cNX - 1 || j == 0 || j == cNY - 1 || k == 0 || k == cNZ - 1)
				{
					// Far flow boundary

					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else if (i == cNX - 2)
				{
					// Outlet => Pressure set to zero
					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else if (i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA P = 0 )
					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else
				{
					CoefsData aijkP;

					aijkP.aimjk = cRHO * cdA * cdA / caijkU[(j * cNZ * cNX) + (k * cNX) + i];
					aijkP.aipjk = cRHO * cdA * cdA / caijkU[(j * cNZ * cNX) + (k * cNX) + i + 1];

					aijkP.aijmk = cRHO * cdA * cdA / caijkV[(j * cNZ * cNX) + (k * cNX) + i];
					aijkP.aijpk = cRHO * cdA * cdA / caijkV[((j + 1) * cNZ * cNX) + (k * cNX) + i];

					aijkP.aijkm = cRHO * cdA * cdA / caijkW[(j * cNZ * cNX) + (k * cNX) + i];
					aijkP.aijkp = cRHO * cdA * cdA / caijkW[(j * cNZ * cNX) + ((k + 1) * cNX) + i];

					if (i == 1)
						aijkP.aimjk = 0.0;
					else if (i == cNX - 2)
						aijkP.aipjk = 0.0;
					if (j == 1)
						aijkP.aijmk = 0.0;
					else if (j == cNY - 2)
						aijkP.aijpk = 0.0;
					if (k == 1)
						aijkP.aijkm = 0.0;
					else if (k == cNZ - 2)
						aijkP.aijkp = 0.0;

					// PAREDES
					// WALL (1)
					if (i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd)
					{
						aijkP.aipjk = 0.0;
					}
					// WALL (2)
					else if (i == ciCyEnd && k >= ckCyInit && k < ckCyEnd)
					{
						aijkP.aimjk = 0.0;
					}

					// WALL (4)
					if (k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd)
					{
						aijkP.aijkp = 0.0;
					}
					// WALL (3)
					else if (k == ckCyEnd && i >= ciCyInit && i < ciCyEnd)
					{
						aijkP.aijkm = 0.0;
					}



					aijkP.aijk = (aijkP.aimjk + aijkP.aipjk +
						aijkP.aijmk + aijkP.aijpk + aijkP.aijkm +
						aijkP.aijkp);

					aijkP.source = cRHO * cdA * (cUHatVelField[(j * cNZ * cNX) + (k * cNX) + i] - cUHatVelField[(j * cNZ * cNX) + (k * cNX) + i + 1] +
						cVHatVelField[(j * cNZ * cNX) + (k * cNX) + i] - cVHatVelField[((j + 1) * cNZ * cNX) + (k * cNX) + i] +
						cWHatVelField[(j * cNZ * cNX) + (k * cNX) + i] - cWHatVelField[(j * cNZ * cNX) + ((k + 1) * cNX) + i]);

					// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]

					// jki
					col.push_back(index + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(aijkP.aijk);

					// jkmi
					col.push_back(index + (-1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkP.aijkm);

					// jkim
					col.push_back(index - 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkP.aimjk);

					// jkip
					col.push_back(index + 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkP.aipjk);

					// jkpi
					col.push_back(index + (1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkP.aijkp);

					// jmki
					col.push_back(index + (0 * cNX) + (-1 * cNX * cNZ));
					val.push_back(-aijkP.aijmk);

					// jpki
					col.push_back(index + (0 * cNX) + (1 * cNX * cNZ));
					val.push_back(-aijkP.aijpk);

					rhs[index] = aijkP.source;
				}

				ptr.push_back(col.size());
			}
		}
	}

	return n2;
}

void FieldRectCylinder::CalcUHatValuesFI(const doubleField3D& uField0, double dt)
{	
	// COPIAR TODOS OS VALORES DA CONDIÇÃO DE CONTORNO ( OBSERVE QUE DEPENDE DA GEOMETRIA SE QUISER SER MAIS EFICIENTE )
	//cUHatVelField = cUVelField;

	memcpy(cUHatVelField.data(), cUVelField.data(), cUVelField.size() * sizeof(double));

#pragma omp parallel for
	for (int j = 1; j < cNY - 1; j++)
	{
		for (int k = 1; k < cNZ - 1; k++)
		{
			for (int i = 2; i < cNX - 1; i++)
			{
				const int index = i + (k * cNX) + (j * cNX * cNZ);

				if (i >= ciCyInit && i < ciCyEnd + 1 && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA U = 0 )
					cUHatVelField[index] = 0.0;					
				}
				else
				{
					CoefsData aijkU;

					const double Fw = GetFwInternalUMomentum(i, j, k);
					const double Fe = GetFeInternalUMomentum(i, j, k);
					const double Fs = GetFsInternalUMomentum(i, j, k);
					const double Fn = GetFnInternalUMomentum(i, j, k);
					const double Fb = GetFbInternalUMomentum(i, j, k);
					const double Ft = GetFtInternalUMomentum(i, j, k);

					const double deltaf = cdA * (Fe - Fw + Fn - Fs + Ft - Fb);

					aijkU.aimjk = cD + (cdA * fmax(Fw, 0.0));
					aijkU.aipjk = cD + (cdA * fmax(-Fe, 0.0));

					aijkU.aijmk = cD + (cdA * fmax(Fs, 0.0));
					aijkU.aijpk = cD + (cdA * fmax(-Fn, 0.0));

					aijkU.aijkm = cD + (cdA * fmax(Fb, 0.0));
					aijkU.aijkp = cD + (cdA * fmax(-Ft, 0.0));

					const double alfaw = (Fw > 0.0) ? 1.0 : 0.0;
					const double alfae = (Fe > 0.0) ? 1.0 : 0.0;
					const double alfas = (Fs > 0.0) ? 1.0 : 0.0;
					const double alfan = (Fn > 0.0) ? 1.0 : 0.0;
					const double alfat = (Ft > 0.0) ? 1.0 : 0.0;
					const double alfab = (Fb > 0.0) ? 1.0 : 0.0;

					const double reneg = ren(i, j, k, cUVelField);
					const double rwneg = rwn(i, j, k, cUVelField);
					const double rnneg = rnn(i, j, k, cUVelField);
					const double rsneg = rsn(i, j, k, cUVelField);
					const double rtneg = rtn(i, j, k, cUVelField);
					const double rbneg = rbn(i, j, k, cUVelField);

					const double repos = rep(i, j, k, cUVelField);
					const double rwpos = rwp(i, j, k, true, cUVelField);
					const double rnpos = rnp(i, j, k, cUVelField);
					const double rspos = rsp(i, j, k, false, cUVelField);
					const double rtpos = rtp(i, j, k, cUVelField);
					const double rbpos = rbp(i, j, k, false, cUVelField);

					const double phiE = cUVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phiP = cUVelField[index];
					const double phiW = cUVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phiN = cUVelField[((j + 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiS = cUVelField[((j - 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiT = cUVelField[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
					const double phiB = cUVelField[(j * cNZ * cNX) + ((k - 1) * cNX) + i];

					const double SuDc = cdA * 0.5 * (
						(Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
						+ (Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
						+ (Fn * (((1.0 - alfan) * Psir(rnneg)) - (alfan * Psir(rnpos))) * (phiN - phiP))
						+ (Fs * (-((1.0 - alfas) * Psir(rsneg)) + (alfas * Psir(rspos))) * (phiP - phiS))
						+ (Ft * (((1.0 - alfat) * Psir(rtneg)) - (alfat * Psir(rtpos))) * (phiT - phiP))
						+ (Fb * (-((1.0 - alfab) * Psir(rbneg)) + (alfab * Psir(rbpos))) * (phiP - phiB))
						);
									   			

					// PAREDES
					// WALL (1)
					// DESNECESSARIO (só resolver normalmente)

					// WALL (2)
					// DESNECESSARIO (só resolver normalmente)

					// WALL (4)
					if (k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd + 1)
					{
						aijkU.aijkp = 0.0;
					}
					// WALL (3)
					else if (k == ckCyEnd && i >= ciCyInit && i < ciCyEnd + 1)
					{
						aijkU.aijkm = 0.0;
					}

					aijkU.aijk = aijkU.aimjk + aijkU.aipjk + aijkU.aijmk
						+ aijkU.aijpk + aijkU.aijkm + aijkU.aijkp + deltaf - cSpUField[index]
						+ (cRHO * cdV / dt);

					caijkU[index] = aijkU.aijk;

					// sem pressão (hat bruh)
					aijkU.source = SuDc +
						/*cSuUField[(j * cNZ * cNX) + (k * cNX) + i] +*/
						((cRHO * cdd * cdd * cdd / dt) * uField0[index]);

					cUHatVelField[index] = ((aijkU.aimjk * phiW) + (aijkU.aipjk * phiE)
						+ (aijkU.aijmk * phiS) + (aijkU.aijpk * phiN) + (aijkU.aijkm * phiB)
						+ (aijkU.aijkp * phiT) + aijkU.source) / aijkU.aijk;
				}				
			}
		}
	}

}
void FieldRectCylinder::CalcVHatValuesFI(const doubleField3D& vField0, double dt)
{
	// COPIAR TODOS OS VALORES DA CONDIÇÃO DE CONTORNO ( OBSERVE QUE DEPENDE DA GEOMETRIA SE QUISER SER MAIS EFICIENTE )
	//cVHatVelField = cVVelField;

	memcpy(cVHatVelField.data(), cVVelField.data(), cVVelField.size() * sizeof(double));

#pragma omp parallel for
	for (int j = 2; j < cNY - 1; j++)
	{
		for (int k = 1; k < cNZ - 1; k++)
		{
			for (int i = 1; i < cNX - 1; i++)
			{
				const int index = i + (k*cNX) + (j*cNX*cNZ);

				if (i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA V = 0 )
					cVHatVelField[index] = 0.0;
				}
				else
				{
					CoefsData aijkV;

					const double Fw = GetFwInternalVMomentum(i, j, k) * cdA;
					const double Fe = GetFeInternalVMomentum(i, j, k) * cdA;
					const double Fs = GetFsInternalVMomentum(i, j, k) * cdA;
					const double Fn = GetFnInternalVMomentum(i, j, k) * cdA;
					const double Fb = GetFbInternalVMomentum(i, j, k) * cdA;
					const double Ft = GetFtInternalVMomentum(i, j, k) * cdA;

					const double deltaf = Fe - Fw + Fn - Fs + Ft - Fb;

					aijkV.aimjk = cD + fmax(Fw, 0.0);
					aijkV.aipjk = cD + fmax(-Fe, 0.0);

					aijkV.aijmk = cD + fmax(Fs, 0.0);
					aijkV.aijpk = cD + fmax(-Fn, 0.0);

					aijkV.aijkm = cD + fmax(Fb, 0.0);
					aijkV.aijkp = cD + fmax(-Ft, 0.0);

					const double alfaw = (Fw > 0) ? 1.0 : 0.0;
					const double alfae = (Fe > 0) ? 1.0 : 0.0;
					const double alfas = (Fs > 0) ? 1.0 : 0.0;
					const double alfan = (Fn > 0) ? 1.0 : 0.0;
					const double alfab = (Fb > 0) ? 1.0 : 0.0;
					const double alfat = (Ft > 0) ? 1.0 : 0.0;

					const double reneg = ren(i, j, k, cVVelField);
					const double rwneg = rwn(i, j, k, cVVelField);
					const double rnneg = rnn(i, j, k, cVVelField);
					const double rsneg = rsn(i, j, k, cVVelField);
					const double rtneg = rtn(i, j, k, cVVelField);
					const double rbneg = rbn(i, j, k, cVVelField);

					const double repos = rep(i, j, k, cVVelField);
					const double rwpos = rwp(i, j, k, false, cVVelField);
					const double rnpos = rnp(i, j, k, cVVelField);
					const double rspos = rsp(i, j, k, true, cVVelField);
					const double rtpos = rtp(i, j, k, cVVelField);
					const double rbpos = rbp(i, j, k, false, cVVelField);

					const double phiE = cVVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phiP = cVVelField[(j * cNZ * cNX) + (k * cNX) + i];
					const double phiW = cVVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phiN = cVVelField[((j + 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiS = cVVelField[((j - 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiT = cVVelField[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
					const double phiB = cVVelField[(j * cNZ * cNX) + ((k - 1) * cNX) + i];


					const double SuDc = (0.5 * Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
						+ (0.5 * Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
						+ (0.5 * Fn * (((1.0 - alfan) * Psir(rnneg)) - (alfan * Psir(rnpos))) * (phiN - phiP))
						+ (0.5 * Fs * (-((1.0 - alfas) * Psir(rsneg)) + (alfas * Psir(rspos))) * (phiP - phiS))
						+ (0.5 * Ft * (((1.0 - alfat) * Psir(rtneg)) - (alfat * Psir(rtpos))) * (phiT - phiP))
						+ (0.5 * Fb * (-((1.0 - alfab) * Psir(rbneg)) + (alfab * Psir(rbpos))) * (phiP - phiB));

					// PAREDES
					// WALL (1)
					if (i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd)
					{
						aijkV.aipjk = 0.0;
					}
					// WALL (2)
					else if (i == ciCyEnd && k >= ckCyInit && k < ckCyEnd)
					{
						aijkV.aimjk = 0.0;
					}

					// WALL (4)
					if (k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd)
					{
						aijkV.aijkp = 0.0;
					}
					// WALL (3)
					else if (k == ckCyEnd && i >= ciCyInit && i < ciCyEnd)
					{
						aijkV.aijkm = 0.0;
					}
					

					aijkV.aijk = aijkV.aimjk + aijkV.aipjk + aijkV.aijmk
						+ aijkV.aijpk + aijkV.aijkm + aijkV.aijkp + deltaf - cSpVField[(j * cNZ * cNX) + (k * cNX) + i]
						+ (cRHO * cdV / dt);

					caijkV[(j * cNZ * cNX) + (k * cNX) + i] = aijkV.aijk;

					// no pressure comp. (hat bruh)
					aijkV.source = SuDc +
						/*cSuVField[(j * cNZ * cNX) + (k * cNX) + i] + */
						((cRHO * cdd * cdd * cdd / dt) * vField0[(j * cNZ * cNX) + (k * cNX) + i]);

					cVHatVelField[index] = ((aijkV.aimjk * phiW) + (aijkV.aipjk * phiE) + (aijkV.aijmk * phiS) +
						(aijkV.aijpk * phiN) + (aijkV.aijkm * phiB) + (aijkV.aijkp * phiT) + aijkV.source) / aijkV.aijk;
				}
			}
		}
	}
}
void FieldRectCylinder::CalcWHatValuesFI(const doubleField3D& wField0, double dt)
{
	// COPIAR TODOS OS VALORES DA CONDIÇÃO DE CONTORNO ( OBSERVE QUE DEPENDE DA GEOMETRIA SE QUISER SER MAIS EFICIENTE )
	//cWHatVelField = cWVelField;

	memcpy(cWHatVelField.data(), cWVelField.data(), cWVelField.size() * sizeof(double));

#pragma omp parallel for
	for (int j = 1; j < cNY - 1; j++)
	{
		for (int k = 2; k < cNZ - 1; k++)
		{
			for (int i = 1; i < cNX - 1; i++)
			{
				const int index = i + (k*cNX) + (j*cNX*cNZ);

				if (i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd + 1) 
				{
					// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA V = 0 )
					cWHatVelField[index] = 0.0;
				}
				else
				{
					CoefsData aijkW;

					const double Fw = GetFwInternalWMomentum(i, j, k) * cdA;
					const double Fe = GetFeInternalWMomentum(i, j, k) * cdA;
					const double Fs = GetFsInternalWMomentum(i, j, k) * cdA;
					const double Fn = GetFnInternalWMomentum(i, j, k) * cdA;
					const double Fb = GetFbInternalWMomentum(i, j, k) * cdA;
					const double Ft = GetFtInternalWMomentum(i, j, k) * cdA;

					const double deltaf = Fe - Fw + Fn - Fs + Ft - Fb;

					aijkW.aimjk = cD + fmax(Fw, 0.0);
					aijkW.aipjk = cD + fmax(-Fe, 0.0);

					aijkW.aijmk = cD + fmax(Fs, 0.0);
					aijkW.aijpk = cD + fmax(-Fn, 0.0);

					aijkW.aijkm = cD + fmax(Fb, 0.0);
					aijkW.aijkp = cD + fmax(-Ft, 0.0);

					const double alfaw = (Fw > 0) ? 1.0 : 0.0;
					const double alfae = (Fe > 0) ? 1.0 : 0.0;
					const double alfas = (Fs > 0) ? 1.0 : 0.0;
					const double alfan = (Fn > 0) ? 1.0 : 0.0;
					const double alfab = (Fb > 0) ? 1.0 : 0.0;
					const double alfat = (Ft > 0) ? 1.0 : 0.0;

					const double reneg = ren(i, j, k, cWVelField);
					const double rwneg = rwn(i, j, k, cWVelField);
					const double rnneg = rnn(i, j, k, cWVelField);
					const double rsneg = rsn(i, j, k, cWVelField);
					const double rtneg = rtn(i, j, k, cWVelField);
					const double rbneg = rbn(i, j, k, cWVelField);

					const double repos = rep(i, j, k, cWVelField);
					const double rwpos = rwp(i, j, k, false, cWVelField);
					const double rnpos = rnp(i, j, k, cWVelField);
					const double rspos = rsp(i, j, k, false, cWVelField);
					const double rtpos = rtp(i, j, k, cWVelField);
					const double rbpos = rbp(i, j, k, true, cWVelField);

					const double phiE = cWVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phiP = cWVelField[(j * cNZ * cNX) + (k * cNX) + i];
					const double phiW = cWVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phiN = cWVelField[((j + 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiS = cWVelField[((j - 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiT = cWVelField[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
					const double phiB = cWVelField[(j * cNZ * cNX) + ((k - 1) * cNX) + i];

					const double SuDc = (0.5 * Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
						+ (0.5 * Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
						+ (0.5 * Fn * (((1.0 - alfan) * Psir(rnneg)) - (alfan * Psir(rnpos))) * (phiN - phiP))
						+ (0.5 * Fs * (-((1.0 - alfas) * Psir(rsneg)) + (alfas * Psir(rspos))) * (phiP - phiS))
						+ (0.5 * Ft * (((1.0 - alfat) * Psir(rtneg)) - (alfat * Psir(rtpos))) * (phiT - phiP))
						+ (0.5 * Fb * (-((1.0 - alfab) * Psir(rbneg)) + (alfab * Psir(rbpos))) * (phiP - phiB));

					// PAREDES
					// WALL (1)
					if (i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd + 1)
					{
						aijkW.aipjk = 0.0;
					}
					// WALL (2)
					else if (i == ciCyEnd && k >= ckCyInit && k < ckCyEnd + 1)
					{
						aijkW.aimjk = 0.0;
					}

					// WALL (4)
					// RESOLVER NORMALMENTE

					// WALL (3)
					// RESOLVER NORMALMENTE


					aijkW.aijk = aijkW.aimjk + aijkW.aipjk + aijkW.aijmk
						+ aijkW.aijpk + aijkW.aijkm + aijkW.aijkp + deltaf - cSpWField[(j * cNZ * cNX) + (k * cNX) + i]
						+ (cRHO * cdV / dt);

					caijkW[(j * cNZ * cNX) + (k * cNX) + i] = aijkW.aijk;

					// no pressure comp. (hat bruh)
					aijkW.source = SuDc +
						((cRHO * cdd * cdd * cdd / dt) * wField0[(j * cNZ * cNX) + (k * cNX) + i]);

					cWHatVelField[index] = ((aijkW.aimjk * phiW) + (aijkW.aipjk * phiE) + (aijkW.aijmk * phiS) +
						(aijkW.aijpk * phiN) + (aijkW.aijkm * phiB) + (aijkW.aijkp * phiT) + aijkW.source) / aijkW.aijk;
				}
			}
		}
	}
}

void FieldRectCylinder::CalcUHatValuesCN(const doubleField3D& uField0, double dt)
{
	// COPIAR TODOS OS VALORES DA CONDIÇÃO DE CONTORNO ( OBSERVE QUE DEPENDE DA GEOMETRIA SE QUISER SER MAIS EFICIENTE )
	//cUHatVelField = cUVelField;

	memcpy(cUHatVelField.data(), cUVelField.data(), cUVelField.size() * sizeof(double));

#pragma omp parallel for
	for (int j = 1; j < cNY - 1; j++)
	{
		for (int k = 1; k < cNZ - 1; k++)
		{
			for (int i = 2; i < cNX - 1; i++)
			{
				const int index = i + (k * cNX) + (j * cNX * cNZ);

				if (i >= ciCyInit && i < ciCyEnd + 1 && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA U = 0 )
					cUHatVelField[index] = 0.0;
				}
				else
				{
					CoefsData aijkU;

					const double Fw = GetFwInternalUMomentum(i, j, k);
					const double Fe = GetFeInternalUMomentum(i, j, k);
					const double Fs = GetFsInternalUMomentum(i, j, k);
					const double Fn = GetFnInternalUMomentum(i, j, k);
					const double Fb = GetFbInternalUMomentum(i, j, k);
					const double Ft = GetFtInternalUMomentum(i, j, k);

					const double deltaf = cdA * (Fe - Fw + Fn - Fs + Ft - Fb);

					aijkU.aimjk = cD + (cdA * fmax(Fw, 0.0));
					aijkU.aipjk = cD + (cdA * fmax(-Fe, 0.0));

					aijkU.aijmk = cD + (cdA * fmax(Fs, 0.0));
					aijkU.aijpk = cD + (cdA * fmax(-Fn, 0.0));

					aijkU.aijkm = cD + (cdA * fmax(Fb, 0.0));
					aijkU.aijkp = cD + (cdA * fmax(-Ft, 0.0));

					const double alfaw = (Fw > 0.0) ? 1.0 : 0.0;
					const double alfae = (Fe > 0.0) ? 1.0 : 0.0;
					const double alfas = (Fs > 0.0) ? 1.0 : 0.0;
					const double alfan = (Fn > 0.0) ? 1.0 : 0.0;
					const double alfat = (Ft > 0.0) ? 1.0 : 0.0;
					const double alfab = (Fb > 0.0) ? 1.0 : 0.0;
					
					const double reneg = ren(i, j, k, cUVelField);
					const double rwneg = rwn(i, j, k, cUVelField);
					const double rnneg = rnn(i, j, k, cUVelField);
					const double rsneg = rsn(i, j, k, cUVelField);
					const double rtneg = rtn(i, j, k, cUVelField);
					const double rbneg = rbn(i, j, k, cUVelField);
					
					const double repos = rep(i, j, k, cUVelField);
					const double rwpos = rwp(i, j, k, true, cUVelField);
					const double rnpos = rnp(i, j, k, cUVelField);
					const double rspos = rsp(i, j, k, false, cUVelField);
					const double rtpos = rtp(i, j, k, cUVelField);
					const double rbpos = rbp(i, j, k, false, cUVelField);

					const double phiE = cUVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phiP = cUVelField[index];
					const double phiW = cUVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phiN = cUVelField[((j + 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiS = cUVelField[((j - 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiT = cUVelField[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
					const double phiB = cUVelField[(j * cNZ * cNX) + ((k - 1) * cNX) + i];

					const double phi0E = uField0[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phi0P = uField0[index];
					const double phi0W = uField0[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phi0N = uField0[((j + 1) * cNZ * cNX) + (k * cNX) + i];
					const double phi0S = uField0[((j - 1) * cNZ * cNX) + (k * cNX) + i];
					const double phi0T = uField0[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
					const double phi0B = uField0[(j * cNZ * cNX) + ((k - 1) * cNX) + i];

					const double SuDc = cdA * 0.5 * (
						(Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
						+ (Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
						+ (Fn * (((1.0 - alfan) * Psir(rnneg)) - (alfan * Psir(rnpos))) * (phiN - phiP))
						+ (Fs * (-((1.0 - alfas) * Psir(rsneg)) + (alfas * Psir(rspos))) * (phiP - phiS))
						+ (Ft * (((1.0 - alfat) * Psir(rtneg)) - (alfat * Psir(rtpos))) * (phiT - phiP))
						+ (Fb * (-((1.0 - alfab) * Psir(rbneg)) + (alfab * Psir(rbpos))) * (phiP - phiB))
						);


					// PAREDES
					// WALL (1)
					// DESNECESSARIO (só resolver normalmente)

					// WALL (2)
					// DESNECESSARIO (só resolver normalmente)

					// WALL (4)
					if (k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd + 1)
					{
						aijkU.aijkp = 0.0;
					}
					// WALL (3)
					else if (k == ckCyEnd && i >= ciCyInit && i < ciCyEnd + 1)
					{
						aijkU.aijkm = 0.0;
					}

					aijkU.aijk = 0.5*(aijkU.aimjk + aijkU.aipjk + aijkU.aijmk
						+ aijkU.aijpk + aijkU.aijkm + aijkU.aijkp) + deltaf - (0.5*cSpUField[index])
						+ (cRHO * cdV / dt);

					caijkU[index] = aijkU.aijk;

					// sem pressão (hat bruh)
					aijkU.source = SuDc + (0.5 * cSpUField[(j * cNZ * cNX) + (k * cNX) + i]*phi0P)+
						/*cSuUField[(j * cNZ * cNX) + (k * cNX) + i] +*/
						(( (cRHO * cdd * cdd * cdd / dt) 
							- (aijkU.aipjk/2.0) - (aijkU.aimjk / 2.0)
							- (aijkU.aijpk / 2.0) - (aijkU.aijmk / 2.0)
							- (aijkU.aijkp / 2.0) - (aijkU.aijkm / 2.0)) * phi0P);

					cUHatVelField[index] = ((aijkU.aimjk * 0.5*(phiW+phi0W)) + (aijkU.aipjk * 0.5*(phiE+phi0E))
						+ (aijkU.aijmk * 0.5*(phiS+phi0S)) + (aijkU.aijpk * 0.5*(phiN+phi0N)) + (aijkU.aijkm * 0.5*(phiB+phi0B))
						+ (aijkU.aijkp * 0.5*(phiT+phi0T)) + aijkU.source) / aijkU.aijk;
				}
			}
		}
	}

}
void FieldRectCylinder::CalcVHatValuesCN(const doubleField3D& vField0, double dt)
{
	// COPIAR TODOS OS VALORES DA CONDIÇÃO DE CONTORNO ( OBSERVE QUE DEPENDE DA GEOMETRIA SE QUISER SER MAIS EFICIENTE )
	//cVHatVelField = cVVelField;

	memcpy(cVHatVelField.data(), cVVelField.data(), cVVelField.size() * sizeof(double));

#pragma omp parallel for
	for (int j = 2; j < cNY - 1; j++)
	{
		for (int k = 1; k < cNZ - 1; k++)
		{
			for (int i = 1; i < cNX - 1; i++)
			{
				const int index = i + (k * cNX) + (j * cNX * cNZ);

				if (i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA V = 0 )
					cVHatVelField[index] = 0.0;
				}
				else
				{
					CoefsData aijkV;

					const double Fw = GetFwInternalVMomentum(i, j, k) * cdA;
					const double Fe = GetFeInternalVMomentum(i, j, k) * cdA;
					const double Fs = GetFsInternalVMomentum(i, j, k) * cdA;
					const double Fn = GetFnInternalVMomentum(i, j, k) * cdA;
					const double Fb = GetFbInternalVMomentum(i, j, k) * cdA;
					const double Ft = GetFtInternalVMomentum(i, j, k) * cdA;

					const double deltaf = Fe - Fw + Fn - Fs + Ft - Fb;

					aijkV.aimjk = cD + fmax(Fw, 0.0);
					aijkV.aipjk = cD + fmax(-Fe, 0.0);

					aijkV.aijmk = cD + fmax(Fs, 0.0);
					aijkV.aijpk = cD + fmax(-Fn, 0.0);

					aijkV.aijkm = cD + fmax(Fb, 0.0);
					aijkV.aijkp = cD + fmax(-Ft, 0.0);

					const double alfaw = (Fw > 0) ? 1.0 : 0.0;
					const double alfae = (Fe > 0) ? 1.0 : 0.0;
					const double alfas = (Fs > 0) ? 1.0 : 0.0;
					const double alfan = (Fn > 0) ? 1.0 : 0.0;
					const double alfab = (Fb > 0) ? 1.0 : 0.0;
					const double alfat = (Ft > 0) ? 1.0 : 0.0;
					
					const double reneg = ren(i, j, k, cVVelField);
					const double rwneg = rwn(i, j, k, cVVelField);
					const double rnneg = rnn(i, j, k, cVVelField);
					const double rsneg = rsn(i, j, k, cVVelField);
					const double rtneg = rtn(i, j, k, cVVelField);
					const double rbneg = rbn(i, j, k, cVVelField);
					
					const double repos = rep(i, j, k, cVVelField);
					const double rwpos = rwp(i, j, k, false, cVVelField);
					const double rnpos = rnp(i, j, k, cVVelField);
					const double rspos = rsp(i, j, k, true, cVVelField);
					const double rtpos = rtp(i, j, k, cVVelField);
					const double rbpos = rbp(i, j, k, false, cVVelField);

					const double phiE = cVVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phiP = cVVelField[(j * cNZ * cNX) + (k * cNX) + i];
					const double phiW = cVVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phiN = cVVelField[((j + 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiS = cVVelField[((j - 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiT = cVVelField[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
					const double phiB = cVVelField[(j * cNZ * cNX) + ((k - 1) * cNX) + i];

					const double phi0E = vField0[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phi0P = vField0[(j * cNZ * cNX) + (k * cNX) + i];
					const double phi0W = vField0[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phi0N = vField0[((j + 1) * cNZ * cNX) + (k * cNX) + i];
					const double phi0S = vField0[((j - 1) * cNZ * cNX) + (k * cNX) + i];
					const double phi0T = vField0[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
					const double phi0B = vField0[(j * cNZ * cNX) + ((k - 1) * cNX) + i];


					const double SuDc = (0.5 * Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
						+ (0.5 * Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
						+ (0.5 * Fn * (((1.0 - alfan) * Psir(rnneg)) - (alfan * Psir(rnpos))) * (phiN - phiP))
						+ (0.5 * Fs * (-((1.0 - alfas) * Psir(rsneg)) + (alfas * Psir(rspos))) * (phiP - phiS))
						+ (0.5 * Ft * (((1.0 - alfat) * Psir(rtneg)) - (alfat * Psir(rtpos))) * (phiT - phiP))
						+ (0.5 * Fb * (-((1.0 - alfab) * Psir(rbneg)) + (alfab * Psir(rbpos))) * (phiP - phiB));

					// PAREDES
					// WALL (1)
					if (i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd)
					{
						aijkV.aipjk = 0.0;
					}
					// WALL (2)
					else if (i == ciCyEnd && k >= ckCyInit && k < ckCyEnd)
					{
						aijkV.aimjk = 0.0;
					}

					// WALL (4)
					if (k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd)
					{
						aijkV.aijkp = 0.0;
					}
					// WALL (3)
					else if (k == ckCyEnd && i >= ciCyInit && i < ciCyEnd)
					{
						aijkV.aijkm = 0.0;
					}


					aijkV.aijk = 0.5*(aijkV.aimjk + aijkV.aipjk + aijkV.aijmk
						+ aijkV.aijpk + aijkV.aijkm + aijkV.aijkp) + deltaf - (0.5*cSpVField[(j * cNZ * cNX) + (k * cNX) + i])
						+ (cRHO * cdV / dt);

					caijkV[(j * cNZ * cNX) + (k * cNX) + i] = aijkV.aijk;

					// no pressure comp. (hat bruh)
					aijkV.source = SuDc +
						(0.5 * cSpVField[(j * cNZ * cNX) + (k * cNX) + i] * phi0P)  +
						(((cRHO * cdd * cdd * cdd / dt)
							- (aijkV.aipjk / 2.0) - (aijkV.aimjk / 2.0)
							- (aijkV.aijpk / 2.0) - (aijkV.aijmk / 2.0)
							- (aijkV.aijkp / 2.0) - (aijkV.aijkm / 2.0)) * phi0P);

					cVHatVelField[index] = ((aijkV.aimjk * 0.5*(phiW+phi0W)) + (aijkV.aipjk *0.5*(phiE+phi0E)) + 
						(aijkV.aijmk * 0.5*(phiS+phi0S)) + (aijkV.aijpk * 0.5*(phiN+phi0N)) + (aijkV.aijkm * 0.5*(phiB+phi0B)) + 
						(aijkV.aijkp * 0.5*(phiT+phi0T)) + aijkV.source) / aijkV.aijk;
				}
			}
		}
	}
}
void FieldRectCylinder::CalcWHatValuesCN(const doubleField3D& wField0, double dt)
{
	// COPIAR TODOS OS VALORES DA CONDIÇÃO DE CONTORNO ( OBSERVE QUE DEPENDE DA GEOMETRIA SE QUISER SER MAIS EFICIENTE )
	//cWHatVelField = cWVelField;

	memcpy(cWHatVelField.data(), cWVelField.data(), cWVelField.size() * sizeof(double));

#pragma omp parallel for
	for (int j = 1; j < cNY - 1; j++)
	{
		for (int k = 2; k < cNZ - 1; k++)
		{
			for (int i = 1; i < cNX - 1; i++)
			{
				const int index = i + (k * cNX) + (j * cNX * cNZ);

				if (i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd + 1)
				{
					// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA V = 0 )
					cWHatVelField[index] = 0.0;
				}
				else
				{
					CoefsData aijkW;

					const double Fw = GetFwInternalWMomentum(i, j, k) * cdA;
					const double Fe = GetFeInternalWMomentum(i, j, k) * cdA;
					const double Fs = GetFsInternalWMomentum(i, j, k) * cdA;
					const double Fn = GetFnInternalWMomentum(i, j, k) * cdA;
					const double Fb = GetFbInternalWMomentum(i, j, k) * cdA;
					const double Ft = GetFtInternalWMomentum(i, j, k) * cdA;

					const double deltaf = Fe - Fw + Fn - Fs + Ft - Fb;

					aijkW.aimjk = cD + fmax(Fw, 0.0);
					aijkW.aipjk = cD + fmax(-Fe, 0.0);

					aijkW.aijmk = cD + fmax(Fs, 0.0);
					aijkW.aijpk = cD + fmax(-Fn, 0.0);

					aijkW.aijkm = cD + fmax(Fb, 0.0);
					aijkW.aijkp = cD + fmax(-Ft, 0.0);

					const double alfaw = (Fw > 0) ? 1.0 : 0.0;
					const double alfae = (Fe > 0) ? 1.0 : 0.0;
					const double alfas = (Fs > 0) ? 1.0 : 0.0;
					const double alfan = (Fn > 0) ? 1.0 : 0.0;
					const double alfab = (Fb > 0) ? 1.0 : 0.0;
					const double alfat = (Ft > 0) ? 1.0 : 0.0;
					
					const double reneg = ren(i, j, k, cWVelField);
					const double rwneg = rwn(i, j, k, cWVelField);
					const double rnneg = rnn(i, j, k, cWVelField);
					const double rsneg = rsn(i, j, k, cWVelField);
					const double rtneg = rtn(i, j, k, cWVelField);
					const double rbneg = rbn(i, j, k, cWVelField);
					
					const double repos = rep(i, j, k, cWVelField);
					const double rwpos = rwp(i, j, k, false, cWVelField);
					const double rnpos = rnp(i, j, k, cWVelField);
					const double rspos = rsp(i, j, k, false, cWVelField);
					const double rtpos = rtp(i, j, k, cWVelField);
					const double rbpos = rbp(i, j, k, true, cWVelField);

					const double phiE = cWVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phiP = cWVelField[(j * cNZ * cNX) + (k * cNX) + i];
					const double phiW = cWVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phiN = cWVelField[((j + 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiS = cWVelField[((j - 1) * cNZ * cNX) + (k * cNX) + i];
					const double phiT = cWVelField[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
					const double phiB = cWVelField[(j * cNZ * cNX) + ((k - 1) * cNX) + i];

					const double phi0E = wField0[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phi0P = wField0[(j * cNZ * cNX) + (k * cNX) + i];
					const double phi0W = wField0[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phi0N = wField0[((j + 1) * cNZ * cNX) + (k * cNX) + i];
					const double phi0S = wField0[((j - 1) * cNZ * cNX) + (k * cNX) + i];
					const double phi0T = wField0[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
					const double phi0B = wField0[(j * cNZ * cNX) + ((k - 1) * cNX) + i];

					const double SuDc =  (0.5 * Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
						+ (0.5 * Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
						+ (0.5 * Fn * (((1.0 - alfan) * Psir(rnneg)) - (alfan * Psir(rnpos))) * (phiN - phiP))
						+ (0.5 * Fs * (-((1.0 - alfas) * Psir(rsneg)) + (alfas * Psir(rspos))) * (phiP - phiS))
						+ (0.5 * Ft * (((1.0 - alfat) * Psir(rtneg)) - (alfat * Psir(rtpos))) * (phiT - phiP))
						+ (0.5 * Fb * (-((1.0 - alfab) * Psir(rbneg)) + (alfab * Psir(rbpos))) * (phiP - phiB));

					// PAREDES
					// WALL (1)
					if (i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd + 1)
					{
						aijkW.aipjk = 0.0;
					}
					// WALL (2)
					else if (i == ciCyEnd && k >= ckCyInit && k < ckCyEnd + 1)
					{
						aijkW.aimjk = 0.0;
					}

					// WALL (4)
					// RESOLVER NORMALMENTE

					// WALL (3)
					// RESOLVER NORMALMENTE


					aijkW.aijk = 0.5*(aijkW.aimjk + aijkW.aipjk + aijkW.aijmk
						+ aijkW.aijpk + aijkW.aijkm + aijkW.aijkp) + deltaf - (0.5*cSpWField[(j * cNZ * cNX) + (k * cNX) + i])
						+ (cRHO * cdV / dt);

					caijkW[(j * cNZ * cNX) + (k * cNX) + i] = aijkW.aijk;

					// no pressure comp. (hat bruh)
					aijkW.source = SuDc +
						(0.5 * cSpWField[(j * cNZ * cNX) + (k * cNX) + i] * phi0P) +
						(((cRHO * cdd * cdd * cdd / dt)
							- (aijkW.aipjk / 2.0) - (aijkW.aimjk / 2.0)
							- (aijkW.aijpk / 2.0) - (aijkW.aijmk / 2.0)
							- (aijkW.aijkp / 2.0) - (aijkW.aijkm / 2.0)) * phi0P);

					cWHatVelField[index] = ((aijkW.aimjk * 0.5*(phiW+phi0W)) + (aijkW.aipjk * 0.5*(phiE+phi0E)) + 
						(aijkW.aijmk * 0.5*(phiS+phi0S)) + (aijkW.aijpk * 0.5*(phiN+phi0N)) + (aijkW.aijkm * 0.5*(phiB+phi0B)) + 
						(aijkW.aijkp * 0.5*(phiT+phi0T)) + aijkW.source) / aijkW.aijk;
				}
			}
		}
	}
}

// Set P value for all 'I' and 'J' in K defined;
void FieldRectCylinder::SetPValueIJ(int64_t K, double pValue)
{
	for (int64_t J = 0; J < cNY; J++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			cPressureField[(J * cNZ * cNX) + (K * cNX) + I] = pValue;
		}
	}
}

// Set P value for all 'J' and 'K' in I defined;
void FieldRectCylinder::SetPValueJK(int64_t I, double pValue)
{
	for (int64_t J = 0; J < cNY; J++)
	{
		for (int64_t K = 0; K < cNZ; K++)
		{
			cPressureField[(J * cNZ * cNX) + (K * cNX) + I] = pValue;
		}
	}
}

// Set P value for all 'K' and 'I' in J defined; 
void FieldRectCylinder::SetPValueIK(int64_t J, double pValue)
{
	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			cPressureField[(J * cNZ * cNX) + (K * cNX) + I] = pValue;
		}
	}
}

// Set U value for all 'i' and 'J' in K defined;
void FieldRectCylinder::SetUValueiJ(int64_t K, double uValue)
{
	for (int64_t J = 0; J < cNY; J++)
	{
		for (int64_t i = 0; i < cNX; i++)
		{
			cUVelField[(J * cNZ * cNX) + (K * cNX) + i] = uValue;
		}
	}
}

// Set U value for all 'J' and 'K' in i defined;
void FieldRectCylinder::SetUValueJK(int64_t i, double uValue)
{
	for (int64_t J = 0; J < cNY; J++)
	{
		for (int64_t K = 0; K < cNZ; K++)
		{
			cUVelField[(J * cNZ * cNX) + (K * cNX) + i] = uValue;
		}
	}
}

// Set U value for all 'K' and 'i' in J defined; 
void FieldRectCylinder::SetUValueiK(int64_t J, double uValue)
{
	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t i = 0; i < cNX; i++)
		{
			cUVelField[(J * cNZ * cNX) + (K * cNX) + i] = uValue;
		}
	}
}

// Set V value for all 'I' and 'K' in j defined;
void FieldRectCylinder::SetVValueIK(int64_t j, double uValue)
{
	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			cVVelField[(j * cNZ * cNX) + (K * cNX) + I] = uValue;
		}
	}
}

// Set V value for all 'j' and 'I' in K defined;
void FieldRectCylinder::SetVValueIj(int64_t K, double vValue)
{
	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			cVVelField[(j * cNZ * cNX) + (K * cNX) + I] = vValue;
		}
	}
}

// Set V value for all 'K' and 'j' in I defined; 
void FieldRectCylinder::SetVValuejK(int64_t I, double vValue)
{
	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t K = 0; K < cNZ; K++)
		{
			cVVelField[(j * cNZ * cNX) + (K * cNX) + I] = vValue;
		}
	}
}

// Set W value for all 'I' and 'J' in k defined;
void FieldRectCylinder::SetWValueIJ(int64_t k, double wValue)
{
	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = wValue;
		}
	}
}

// Set W value for all 'J' and 'k' in I defined;
void FieldRectCylinder::SetWValueJk(int64_t I, double wValue)
{
	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = wValue;
		}
	}
}

// Set W value for all 'k' and 'I' in J defined; 
void FieldRectCylinder::SetWValueIk(int64_t J, double wValue)
{
	for (int64_t k = 0; k < cNZ; k++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			cWVelField[(J * cNZ * cNX) + (k * cNX) + I] = wValue;
		}
	}
}

void FieldRectCylinder::SetPValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& KK, double pValue)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = KK.valor1; k < KK.valor2; k++)
		{
			for (int64_t i = II.valor1; i < II.valor2; i++)
			{
				cPressureField[(j * cNZ * cNX) + (k * cNX) + i] = pValue;
			}
		}
	}
}

void FieldRectCylinder::SetUValueInterval(const Intervalo& ii, const Intervalo& JJ, const Intervalo& KK, double uValue)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = KK.valor1; k < KK.valor2; k++)
		{
			for (int64_t i = ii.valor1; i < ii.valor2; i++)
			{
				cUVelField[(j * cNZ * cNX) + (k * cNX) + i] = uValue;
			}
		}
	}
}

void FieldRectCylinder::SetVValueInterval(const Intervalo& II, const Intervalo& jj, const Intervalo& KK, double vValue)
{
	for (int64_t j = jj.valor1; j < jj.valor2; j++)
	{
		for (int64_t k = KK.valor1; k < KK.valor2; k++)
		{
			for (int64_t i = II.valor1; i < II.valor2; i++)
			{
				cVVelField[(j * cNZ * cNX) + (k * cNX) + i] = vValue;
			}
		}
	}
}

void FieldRectCylinder::SetWValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double wValue)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (int64_t i = II.valor1; i < II.valor2; i++)
			{
				cWVelField[(j * cNZ * cNX) + (k * cNX) + i] = wValue;
			}
		}
	}
}

void FieldRectCylinder::ExtrapolateForwardPJK(int64_t I)
{
	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			cPressureField[(j * cNZ * cNX) + (k * cNX) + I] = cPressureField[(j * cNZ * cNX) + (k * cNX) + I + 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateForwardUJK(int64_t i)
{
	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			cUVelField[(j * cNZ * cNX) + (k * cNX) + i] = cUVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateForwardVjK(int64_t I)
{
	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			cVVelField[(j * cNZ * cNX) + (k * cNX) + I] = cVVelField[(j * cNZ * cNX) + (k * cNX) + I + 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateForwardWJk(int64_t I)
{
	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = cWVelField[(j * cNZ * cNX) + (k * cNX) + I + 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateBackwardPJK(int64_t I)
{
	for (int64_t j = 1; j < cNY - 1; j++)
	{
		for (int64_t k = 1; k < cNZ - 1; k++)
		{
			cPressureField[(j * cNZ * cNX) + (k * cNX) + I] = cPressureField[(j * cNZ * cNX) + (k * cNX) + I - 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateBackwardUJK(int64_t i)
{
	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			cUVelField[(j * cNZ * cNX) + (k * cNX) + i] = cUVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateBackwardVjK(int64_t I)
{
	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			cVVelField[(j * cNZ * cNX) + (k * cNX) + I] = cVVelField[(j * cNZ * cNX) + (k * cNX) + I - 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateBackwardWJk(int64_t I)
{
	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = cWVelField[(j * cNZ * cNX) + (k * cNX) + I - 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateBackwardWeightedUJK(int64_t i, double weight)
{
	for (int64_t j = 1; j < cNY - 1; j++)
	{
		for (int64_t k = 1; k < cNZ - 1; k++)
		{
			cUVelField[(j * cNZ * cNX) + (k * cNX) + i] = cUVelField[(j * cNZ * cNX) + (k * cNX) + i - 1] * weight;
		}
	}
}

void FieldRectCylinder::ExtrapolateBackwardWeightedVjK(int64_t I, double weight)
{
	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			cVVelField[(j * cNZ * cNX) + (k * cNX) + I] = cVVelField[(j * cNZ * cNX) + (k * cNX) + I - 1] * weight;
		}
	}
}

void FieldRectCylinder::ExtrapolateBackwardWeightedWJk(int64_t I, double weight)
{
	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = cWVelField[(j * cNZ * cNX) + (k * cNX) + I - 1] * weight;
		}
	}
}

void FieldRectCylinder::ExtrapolateForwardIntervalPJK(int64_t I, const Intervalo& JJ, const Intervalo& KK)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = KK.valor1; k < KK.valor2; k++)
		{
			cPressureField[(j * cNZ * cNX) + (k * cNX) + I] = cPressureField[(j * cNZ * cNX) + (k * cNX) + I + 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateForwardIntervalUJK(int64_t i, const Intervalo& JJ, const Intervalo& KK)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = KK.valor1; k < KK.valor2; k++)
		{
			cUVelField[(j * cNZ * cNX) + (k * cNX) + i] = cUVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateForwardIntervalVjK(int64_t I, const Intervalo& jj, const Intervalo& KK)
{
	for (int64_t j = jj.valor1; j < jj.valor2; j++)
	{
		for (int64_t k = KK.valor1; k < KK.valor2; k++)
		{
			cVVelField[(j * cNZ * cNX) + (k * cNX) + I] = cVVelField[(j * cNZ * cNX) + (k * cNX) + I + 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateForwardIntervalWJk(int64_t I, const Intervalo& JJ, const Intervalo& kk)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = kk.valor1; k < kk.valor2; k++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = cWVelField[(j * cNZ * cNX) + (k * cNX) + I + 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateBackwardIntervalPJK(int64_t I, const Intervalo& JJ, const Intervalo& KK)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = KK.valor1; k < KK.valor2; k++)
		{
			cPressureField[(j * cNZ * cNX) + (k * cNX) + I] = cPressureField[(j * cNZ * cNX) + (k * cNX) + I - 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateBackwardIntervalUJK(int64_t i, const Intervalo& JJ, const Intervalo& KK)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = KK.valor1; k < KK.valor2; k++)
		{
			cUVelField[(j * cNZ * cNX) + (k * cNX) + i] = cUVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateBackwardIntervalVjK(int64_t I, const Intervalo& jj, const Intervalo& KK)
{
	for (int64_t j = jj.valor1; j < jj.valor2; j++)
	{
		for (int64_t k = KK.valor1; k < KK.valor2; k++)
		{
			cVVelField[(j * cNZ * cNX) + (k * cNX) + I] = cVVelField[(j * cNZ * cNX) + (k * cNX) + I - 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateBackwardIntervalWJk(int64_t I, const Intervalo& JJ, const Intervalo& kk)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = kk.valor1; k < kk.valor2; k++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = cWVelField[(j * cNZ * cNX) + (k * cNX) + I - 1];
		}
	}
}

void FieldRectCylinder::ExtrapolateForwardIntervalPJI(int64_t K, const Intervalo& II, const Intervalo& JJ)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t i = II.valor1; i < II.valor2; i++)
		{
			cPressureField[(j * cNZ * cNX) + (K * cNX) + i] = cPressureField[(j * cNZ * cNX) + ((K+1) * cNX) + i];
		}
	}
}

void FieldRectCylinder::ExtrapolateForwardIntervalUJI(int64_t K, const Intervalo& ii, const Intervalo& JJ)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t i = ii.valor1; i < ii.valor2; i++)
		{
			cUVelField[(j * cNZ * cNX) + ((K ) * cNX) + i] = cUVelField[(j * cNZ * cNX) + ((K + 1) * cNX) + i];
		}
	}
}
void FieldRectCylinder::ExtrapolateForwardIntervalVJI(int64_t K, const Intervalo& II, const Intervalo& jj)
{
	for (int64_t j = jj.valor1; j < jj.valor2; j++)
	{
		for (int64_t i = II.valor1; i < II.valor2; i++)
		{
			cVVelField[(j * cNZ * cNX) + ((K) * cNX) + i] = cVVelField[(j * cNZ * cNX) + ((K + 1) * cNX) + i];
		}
	}
}
void FieldRectCylinder::ExtrapolateForwardIntervalWJI(int64_t k, const Intervalo& II, const Intervalo& JJ)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t i = II.valor1; i < II.valor2; i++)
		{
			cVVelField[(j * cNZ * cNX) + ((k) * cNX) + i] = cVVelField[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
		}
	}
}

void FieldRectCylinder::ExtrapolateBackwardIntervalPJI(int64_t K, const Intervalo& II, const Intervalo& JJ)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t i = II.valor1; i < II.valor2; i++)
		{
			cPressureField[(j * cNZ * cNX) + ((K) * cNX) + i] = cPressureField[(j * cNZ * cNX) + ((K - 1) * cNX) + i];
		}
	}
}
void FieldRectCylinder::ExtrapolateBackwardIntervalUJI(int64_t K, const Intervalo& ii, const Intervalo& JJ)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t i = ii.valor1; i < ii.valor2; i++)
		{
			cUVelField[(j * cNZ * cNX) + ((K + 1) * cNX) + i] = cUVelField[(j * cNZ * cNX) + ((K - 1) * cNX) + i];
		}
	}
}
void FieldRectCylinder::ExtrapolateBackwardIntervalVJI(int64_t K, const Intervalo& II, const Intervalo& jj)
{
	for (int64_t j = jj.valor1; j < jj.valor2; j++)
	{
		for (int64_t i = II.valor1; i < II.valor2; i++)
		{
			cVVelField[(j * cNZ * cNX) + ((K) * cNX) + i] = cVVelField[(j * cNZ * cNX) + ((K - 1) * cNX) + i];
		}
	}
}
void FieldRectCylinder::ExtrapolateBackwardIntervalWJI(int64_t k, const Intervalo& II, const Intervalo& JJ)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t i = II.valor1; i < II.valor2; i++)
		{
			cVVelField[(j * cNZ * cNX) + ((k) * cNX) + i] = cVVelField[(j * cNZ * cNX) + ((k - 1) * cNX) + i];
		}
	}
}

void FieldRectCylinder::SaveField(const std::string& baseFileName) const
{
	std::ofstream uFile(baseFileName + "U.dat", std::ios::binary);

	uFile.write(reinterpret_cast<const char*>(&cNX), sizeof(int64_t));
	uFile.write(reinterpret_cast<const char*>(&cNY), sizeof(int64_t));
	uFile.write(reinterpret_cast<const char*>(&cNZ), sizeof(int64_t));

	uFile.write(reinterpret_cast<const char*>(cUVelField.data()), cUVelField.size() * sizeof(double));

	uFile.close();

	std::ofstream vFile(baseFileName + "V.dat", std::ios::binary);

	vFile.write(reinterpret_cast<const char*>(&cNX), sizeof(int64_t));
	vFile.write(reinterpret_cast<const char*>(&cNY), sizeof(int64_t));
	vFile.write(reinterpret_cast<const char*>(&cNZ), sizeof(int64_t));

	vFile.write(reinterpret_cast<const char*>(cVVelField.data()), cVVelField.size() * sizeof(double));

	vFile.close();

	std::ofstream wFile(baseFileName + "W.dat", std::ios::binary);

	wFile.write(reinterpret_cast<const char*>(&cNX), sizeof(int64_t));
	wFile.write(reinterpret_cast<const char*>(&cNY), sizeof(int64_t));
	wFile.write(reinterpret_cast<const char*>(&cNZ), sizeof(int64_t));

	wFile.write(reinterpret_cast<const char*>(cWVelField.data()), cWVelField.size() * sizeof(double));

	wFile.close();

	std::ofstream pFile(baseFileName + "P.dat", std::ios::binary);

	pFile.write(reinterpret_cast<const char*>(&cNX), sizeof(int64_t));
	pFile.write(reinterpret_cast<const char*>(&cNY), sizeof(int64_t));
	pFile.write(reinterpret_cast<const char*>(&cNZ), sizeof(int64_t));

	pFile.write(reinterpret_cast<const char*>(cPressureField.data()), cPressureField.size() * sizeof(double));

	pFile.close();
}

void FieldRectCylinder::ReadField(const std::string& basefilename)
{

	std::ifstream uFile(basefilename + "U.dat", std::ios::binary);

	uFile.read(reinterpret_cast<char*>(&cNX), sizeof(int64_t));
	uFile.read(reinterpret_cast<char*>(&cNY), sizeof(int64_t));
	uFile.read(reinterpret_cast<char*>(&cNZ), sizeof(int64_t));
	
	cUVelField.clear();
	cUVelField.resize(cNX * cNY * cNZ);

	uFile.read(reinterpret_cast<char*>(cUVelField.data()), cUVelField.size() * sizeof(double));

	uFile.close();

	std::ifstream vFile(basefilename + "V.dat", std::ios::binary);

	vFile.read(reinterpret_cast<char*>(&cNX), sizeof(int64_t));
	vFile.read(reinterpret_cast<char*>(&cNY), sizeof(int64_t));
	vFile.read(reinterpret_cast<char*>(&cNZ), sizeof(int64_t));

	cVVelField.clear();
	cVVelField.resize(cNX * cNY * cNZ);

	vFile.read(reinterpret_cast<char*>(cVVelField.data()), cVVelField.size() * sizeof(double));

	vFile.close();

	std::ifstream wFile(basefilename + "W.dat", std::ios::binary);

	wFile.read(reinterpret_cast<char*>(&cNX), sizeof(int64_t));
	wFile.read(reinterpret_cast<char*>(&cNY), sizeof(int64_t));
	wFile.read(reinterpret_cast<char*>(&cNZ), sizeof(int64_t));

	cWVelField.clear();
	cWVelField.resize(cNX * cNY * cNZ);

	wFile.read(reinterpret_cast<char*>(cWVelField.data()), cWVelField.size() * sizeof(double));

	wFile.close();

	std::ifstream pFile(basefilename + "P.dat", std::ios::binary);

	pFile.read(reinterpret_cast<char*>(&cNX), sizeof(int64_t));
	pFile.read(reinterpret_cast<char*>(&cNY), sizeof(int64_t));
	pFile.read(reinterpret_cast<char*>(&cNZ), sizeof(int64_t));

	cPressureField.clear();
	cPressureField.resize(cNX * cNY * cNZ);

	pFile.read(reinterpret_cast<char*>(cPressureField.data()), cPressureField.size() * sizeof(double));

	pFile.close();
}

void FieldRectCylinder::SaveIJCutFieldToCSVFile(int64_t K, const std::string& baseFileName) const
{
	std::ofstream pFile(baseFileName + "P.csv");

	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			pFile << cPressureField[(j * cNZ * cNX) + ((K) * cNX) + I] << "\t";
		}
		pFile << "\n";
	}

	pFile.close();

	std::ofstream uFile(baseFileName + "U.csv");

	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			uFile << cUVelField[(j * cNZ * cNX) + ((K)*cNX) + I] << "\t";
		}
		uFile << "\n";
	}

	uFile.close();

	std::ofstream vFile(baseFileName + "V.csv");

	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			vFile << cVVelField[(j * cNZ * cNX) + ((K)*cNX) + I] << "\t";
		}
		vFile << "\n";
	}

	vFile.close();

	std::ofstream wFile(baseFileName + "W.csv");

	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			wFile << cWVelField[(j * cNZ * cNX) + ((K)*cNX) + I] << "\t";
		}
		wFile << "\n";
	}

	wFile.close();

	std::ofstream uspFile(baseFileName + "USP.csv");

	for (int64_t J = 0; J < cNY; J++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			uspFile << cSpUField[(J * cNZ * cNX) + ((K)*cNX) + I] << "\t";
		}
		uspFile << "\n";
	}

	uspFile.close();
}

void FieldRectCylinder::SaveIKCutFieldToCSVFile(int64_t J, const std::string& baseFileName) const
{
	std::ofstream pFile(baseFileName + "P.csv");

	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			pFile << cPressureField[(J * cNZ * cNX) + ((K)*cNX) + I] << "\t";
		}
		pFile << "\n";
	}

	pFile.close();

	std::ofstream pcorrFile(baseFileName + "Pcorr.csv");
	
	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			pcorrFile << cPressureCorrField[(J * cNX * cNZ) + (K * cNX) + I] << "\t";
		}
		pcorrFile << "\n";
	}
	
	pcorrFile.close();

	std::ofstream uFile(baseFileName + "U.csv");

	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			uFile << cUVelField[(J* cNZ * cNX) + ((K)*cNX) + I] << "\t";
		}
		uFile << "\n";
	}

	uFile.close();

	std::ofstream vFile(baseFileName + "V.csv");

	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			vFile << cVVelField[(J * cNZ * cNX) + ((K)*cNX) + I] << "\t";
		}
		vFile << "\n";
	}

	vFile.close();

	std::ofstream wFile(baseFileName + "W.csv");

	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			wFile << cWVelField[(J * cNZ * cNX) + ((K)*cNX) + I] << "\t";
		}
		wFile << "\n";
	}

	wFile.close();

	//std::ofstream uspFile(baseFileName + "USP.csv");
	//
	//for (int64_t K = 0; K < cSpUField[0].size(); K++)
	//{
	//	for (int64_t I = 0; I < cSpUField[0][0].size(); I++)
	//	{
	//		uspFile << cSpUField[J][K][I] << "\t";
	//	}
	//	uspFile << "\n";
	//}
	//
	//uspFile.close();
}

void FieldRectCylinder::SaveJKCutFieldToCSVFile(int64_t I, const std::string& baseFileName) const
{
	std::ofstream pFile(baseFileName + "P.csv");

	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			pFile << cPressureField[(j * cNZ * cNX) + ((k)*cNX) + I] << "\t";
		}
		pFile << "\n";
	}
	pFile.close();

	std::ofstream pcorrFile(baseFileName + "Pcorr.csv");

	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			pFile << cPressureCorrField[(j * cNZ * cNX) + ((k)*cNX) + I] << "\t";
		}
		pcorrFile << "\n";
	}

	pcorrFile.close();

	std::ofstream uFile(baseFileName + "U.csv");

	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			uFile << cUVelField[(j * cNZ * cNX) + ((k)*cNX) + I] << "\t";
		}
		uFile << "\n";
	}
	uFile.close();

	std::ofstream vFile(baseFileName + "V.csv");

	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			vFile << cVVelField[(j * cNZ * cNX) + ((k) * cNX) + I] << "\t";
		}
		vFile << "\n";
	}
	vFile.close();

	std::ofstream wFile(baseFileName + "W.csv");

	for (int64_t j = 0; j < cNY; j++)
	{
		for (int64_t k = 0; k < cNZ; k++)
		{
			wFile << cWVelField[(j * cNZ * cNX) + ((k)*cNX) + I] << "\t";
		}
		wFile << "\n";
	}
	wFile.close();

	//std::ofstream uspFile(baseFileName + "USP.csv");
	//for (int64_t j = 0; j < cSpUField.size(); j++)
	//{
	//	for (int64_t K = 0; K < cSpUField[0].size(); K++)
	//	{
	//		uspFile << cSpUField[j][K][I] << "\t";
	//	}
	//
	//	uspFile << "\n";
	//}
	//
	//
	//uspFile.close();
}

inline double FieldRectCylinder::GetFwInternalUMomentum(int64_t i, int64_t J, int64_t K) const
{
	return (cRHO / 2.0) * (cUVelField[(J * cNZ * cNX) + ((K)*cNX) + i] + cUVelField[(J * cNZ * cNX) + ((K)*cNX) + i - 1]);
}
inline double FieldRectCylinder::GetFeInternalUMomentum(int64_t i, int64_t J, int64_t K) const
{
	return (cRHO / 2.0) * (cUVelField[(J * cNZ * cNX) + ((K)*cNX) + i + 1] + cUVelField[(J * cNZ * cNX) + ((K)*cNX) + i]);
}
inline double FieldRectCylinder::GetFsInternalUMomentum(int64_t i, int64_t J, int64_t K) const
{
	return (cRHO / 2.0) * (cVVelField[(J * cNZ * cNX) + ((K)*cNX) + i] + cVVelField[(J * cNZ * cNX) + ((K)*cNX) + i - 1]);
}
inline double FieldRectCylinder::GetFnInternalUMomentum(int64_t i, int64_t J, int64_t K) const
{
	return (cRHO / 2.0) * (cVVelField[((J+1) * cNZ * cNX) + ((K)*cNX) + i] + cVVelField[((J+1) * cNZ * cNX) + ((K)*cNX) + i - 1]);
}
inline double FieldRectCylinder::GetFbInternalUMomentum(int64_t i, int64_t J, int64_t K) const
{
	return (cRHO / 2.0) * (cWVelField[(J * cNZ * cNX) + ((K)*cNX) + i] + cWVelField[(J * cNZ * cNX) + ((K)*cNX) + i - 1]);
}
inline double FieldRectCylinder::GetFtInternalUMomentum(int64_t i, int64_t J, int64_t K) const
{
	return (cRHO / 2.0) * (cWVelField[(J * cNZ * cNX) + ((K+1)*cNX) + i] + cWVelField[(J * cNZ * cNX) + ((K + 1)*cNX) + i - 1]);
}

inline double FieldRectCylinder::GetFwInternalVMomentum(int64_t I, int64_t j, int64_t K) const
{
	return (cRHO / 2.0) * (cUVelField[(j * cNZ * cNX) + ((K)*cNX) + I] + cUVelField[((j-1) * cNZ * cNX) + ((K)*cNX) + I]);
}
inline double FieldRectCylinder::GetFeInternalVMomentum(int64_t I, int64_t j, int64_t K) const
{
	return (cRHO / 2.0) * (cUVelField[(j * cNZ * cNX) + ((K)*cNX) + I + 1] + cUVelField[((j-1) * cNZ * cNX) + ((K)*cNX) + I + 1]);
}
inline double FieldRectCylinder::GetFsInternalVMomentum(int64_t I, int64_t j, int64_t K) const
{
	return (cRHO / 2.0) * (cVVelField[(j * cNZ * cNX) + ((K)*cNX) + I] + cVVelField[((j-1) * cNZ * cNX) + ((K)*cNX) + I]);
}
inline double FieldRectCylinder::GetFnInternalVMomentum(int64_t I, int64_t j, int64_t K) const
{
	return (cRHO / 2.0) * (cVVelField[(j * cNZ * cNX) + ((K)*cNX) + I] + cVVelField[((j+1) * cNZ * cNX) + ((K)*cNX) + I]);
}
inline double FieldRectCylinder::GetFbInternalVMomentum(int64_t I, int64_t j, int64_t K) const
{
	return (cRHO / 2.0) * (cWVelField[(j * cNZ * cNX) + ((K)*cNX) + I] + cWVelField[((j-1) * cNZ * cNX) + ((K)*cNX) + I]);
}
inline double FieldRectCylinder::GetFtInternalVMomentum(int64_t I, int64_t j, int64_t K) const
{
	return (cRHO / 2.0) * (cWVelField[(j * cNZ * cNX) + ((K+1)*cNX) + I] + cWVelField[((j-1) * cNZ * cNX) + ((K+1)*cNX) + I]);
}

inline double FieldRectCylinder::GetFwInternalWMomentum(int64_t I, int64_t J, int64_t k) const
{
	return (cRHO / 2.0) * (cUVelField[(J * cNZ * cNX) + ((k)*cNX) + I] + cUVelField[(J * cNZ * cNX) + ((k-1)*cNX) + I]);
}
inline double FieldRectCylinder::GetFeInternalWMomentum(int64_t I, int64_t J, int64_t k) const
{
	return (cRHO / 2.0) * (cUVelField[(J * cNZ * cNX) + ((k)*cNX) + I + 1] + cUVelField[(J * cNZ * cNX) + ((k-1)*cNX) + I + 1]);
}
inline double FieldRectCylinder::GetFsInternalWMomentum(int64_t I, int64_t J, int64_t k) const
{
	return (cRHO / 2.0) * (cVVelField[(J * cNZ * cNX) + ((k)*cNX) + I] + cVVelField[(J * cNZ * cNX) + ((k-1)*cNX) + I]);
}
inline double FieldRectCylinder::GetFnInternalWMomentum(int64_t I, int64_t J, int64_t k) const
{
	return (cRHO / 2.0) * (cVVelField[((J+1) * cNZ * cNX) + ((k)*cNX) + I] + cVVelField[((J+1) * cNZ * cNX) + ((k-1)*cNX) + I]);
}
inline double FieldRectCylinder::GetFbInternalWMomentum(int64_t I, int64_t J, int64_t k) const
{
	return (cRHO / 2.0) * (cWVelField[(J * cNZ * cNX) + ((k)*cNX) + I] + cWVelField[(J * cNZ * cNX) + ((k - 1)*cNX) + I]);
}
inline double FieldRectCylinder::GetFtInternalWMomentum(int64_t I, int64_t J, int64_t k) const
{
	return (cRHO / 2.0) * (cWVelField[(J * cNZ * cNX) + ((k)*cNX) + I] + cWVelField[(J * cNZ * cNX) + ((k+1)*cNX) + I]);
}

void FieldRectCylinder::SetSpUValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SpValue)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (int64_t i = II.valor1; i < II.valor2; i++)
			{
				cSpUField[(j * cNZ * cNX) + ((k)*cNX) + i] = SpValue;
			}
		}
	}
}
void FieldRectCylinder::SetSuUValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SuValue)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (int64_t i = II.valor1; i < II.valor2; i++)
			{
				//cSuUField[(j * cNZ * cNX) + ((k)*cNX) + i] = SuValue;
			}
		}
	}
}

void FieldRectCylinder::SetSpVValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SpValue)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (int64_t i = II.valor1; i < II.valor2; i++)
			{
				cSpVField[(j * cNZ * cNX) + ((k)*cNX) + i] = SpValue;
			}
		}
	}
}
void FieldRectCylinder::SetSuVValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SuValue)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (int64_t i = II.valor1; i < II.valor2; i++)
			{
				//cSuVField[(j * cNZ * cNX) + ((k)*cNX) + i] = SuValue;
			}
		}
	}
}

void FieldRectCylinder::SetSpWValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SpValue)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (int64_t i = II.valor1; i < II.valor2; i++)
			{
				cSpWField[(j * cNZ * cNX) + ((k)*cNX) + i] = SpValue;
			}
		}
	}
}
void FieldRectCylinder::SetSuWValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SuValue)
{
	for (int64_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (int64_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (int64_t i = II.valor1; i < II.valor2; i++)
			{
				//cSuWField[(j * cNZ * cNX) + ((k)*cNX) + i] = SuValue;
			}
		}
	}
}

inline double FieldRectCylinder::Psir(double r) const
{
	return (r + (r * r)) / (1.0 + (r * r));
}

inline double FieldRectCylinder::PsirSUPERBEE(double r) const
{
	const double i1 = std::fmin(2.0*r, 1.0);
	const double i2 = std::fmin(r, 2.0);
	const double res = std::fmax(0.0, std::fmax(i1, i2));

	return res;
}

double FieldRectCylinder::rep(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const
{
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiW = field[(J * cNZ * cNX) + ((K)*cNX) + I - 1];
	const double phiE = field[(J * cNZ * cNX) + ((K)*cNX) + I + 1];

	if (std::fabs(phiE-phiP) < 1e-5)
		return 10e30;
	
	return (phiP - phiW) / (phiE - phiP);
}
double FieldRectCylinder::ren(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const
{
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiE = field[(J * cNZ * cNX) + ((K)*cNX) + I + 1];
	double phiEE;

	if (I == cNX - 2)
		phiEE = (2.0 * phiE) - phiP;
	else
		phiEE = field[(J * cNZ * cNX) + ((K)*cNX) + I + 2];

	if (std::fabs(phiE - phiP) < 1e-5)
		return 10e30;

	return (phiEE - phiE) / (phiE - phiP);
}
	 
double FieldRectCylinder::rwp(int64_t I, int64_t J, int64_t K, bool ufield, const doubleField3D& field) const
{
	const double phiW = field[(J * cNZ * cNX) + ((K)*cNX) + I - 1];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];

	const int iInicial = (ufield) ? 2 : 1;

	double phiWW;
	if (I == iInicial)
		phiWW = (2.0 * phiW) - phiP;
	else
		phiWW = field[(J * cNZ * cNX) + ((K)*cNX) + I - 2];

	if (std::fabs(phiW - phiP) < 1e-5)
		return 10e30;

	return (phiW - phiWW) / (phiP - phiW);
}
double FieldRectCylinder::rwn(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const
{
	const double phiE = field[(J * cNZ * cNX) + ((K)*cNX) + I + 1];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiW = field[(J * cNZ * cNX) + ((K)*cNX) + I - 1];

	if (std::fabs(phiW - phiP) < 1e-5)
		return 10e30;

	return (phiE - phiP) / (phiP - phiW);
}
	 
double FieldRectCylinder::rnp(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const
{
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiS = field[((J-1) * cNZ * cNX) + ((K)*cNX) + I];
	const double phiN = field[((J+1) * cNZ * cNX) + ((K)*cNX) + I];

	if (std::fabs(phiN - phiP) < 1e-5)
		return 10e30;

	return (phiP - phiS) / (phiN - phiP);
}
double FieldRectCylinder::rnn(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const
{
	const double phiN = field[((J+1) * cNZ * cNX) + ((K)*cNX) + I];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];

	double phiNN;
	if (J == cNY - 2)
		phiNN = (2.0 * phiN) - phiP;
	else
		phiNN = field[((J+2) * cNZ * cNX) + ((K)*cNX) + I];

	if (std::fabs(phiN - phiP) < 1e-5)
		return 10e30;

	return (phiNN - phiN) / (phiN - phiP);
}

double FieldRectCylinder::rsp(int64_t I, int64_t J, int64_t K, bool vfield, const doubleField3D& field) const
{
	const double phiS = field[((J-1) * cNZ * cNX) + ((K)*cNX) + I];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];

	const int jInicial = (vfield) ? 2 : 1;

	double phiSS;
	if (J == jInicial)
		phiSS = (phiS * 2.0) - phiP;
	else
		phiSS = field[((J-2) * cNZ * cNX) + ((K)*cNX) + I];

	if (std::fabs(phiS - phiP) < 1e-5)
		return 10e30;

	return (phiS - phiSS) / (phiP - phiS);
}
double FieldRectCylinder::rsn(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const
{
	const double phiN = field[((J+1) * cNZ * cNX) + ((K)*cNX) + I];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiS = field[((J-1) * cNZ * cNX) + ((K)*cNX) + I];

	if (std::fabs(phiS - phiP) < 1e-5)
		return 10e30;

	return (phiN - phiP) / (phiP - phiS);
}
	
double FieldRectCylinder::rbp(int64_t I, int64_t J, int64_t K, bool wfield, const doubleField3D& field) const
{
	const double phiB = field[(J * cNZ * cNX) + ((K-1)*cNX) + I];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];

	const int kInicial = wfield ?  2 : 1;

	double phiBB;
	if (K == kInicial)
		phiBB = (2.0 * phiB) - phiP;
	else
		phiBB = field[(J * cNZ * cNX) + ((K-2)*cNX) + I];

	if (std::fabs(phiB - phiP) < 1e-5)
		return 10e30;

	return (phiB - phiBB) / (phiP - phiB);
}
double FieldRectCylinder::rbn(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const
{
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiT = field[(J * cNZ * cNX) + ((K + 1)*cNX) + I];
	const double phiB = field[(J * cNZ * cNX) + ((K-1)*cNX) + I];

	if (std::fabs(phiP - phiB) < 1e-5)
		return 10e30;
	
	return (phiT - phiP) / (phiP - phiB);
}

double FieldRectCylinder::rtp(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const
{
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiT = field[(J * cNZ * cNX) + ((K+1)*cNX) + I];
	const double phiB = field[(J * cNZ * cNX) + ((K-1)*cNX) + I];

	if (std::fabs(phiP - phiT) < 1e-5)
		return 10e30;
	
	return (phiP - phiB) / (phiT - phiP);
}
double FieldRectCylinder::rtn(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const
{
	const double phiT = field[(J * cNZ * cNX) + ((K+1)*cNX) + I];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];

	double phiTT;
	if (K == cNZ - 2)
		phiTT = (2.0 * phiT) - phiP;
	else
		phiTT = field[(J * cNZ * cNX) + ((K+2)*cNX) + I];

	if (std::fabs(phiP - phiT) < 1e-5)
		return 10e30;

	return (phiTT - phiT) / (phiT - phiP);
}

double FieldRectCylinder::PsirCD(double r) const
{
	return 1.0;
}