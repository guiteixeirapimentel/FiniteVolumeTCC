#include "FieldRectCylinder2D.h"
#include <fstream>
#include <iostream>
#include <ctype.h>

FieldRectCylinder2D::FieldRectCylinder2D(double L, double H, double l, double deltaSize, double rho, double mu, double Ufarfield,
	double Wfarfield)
	:
	ctidPushback(cTimer.SetTimer("Push back time (U mom)")),
	ctidCalcCoefs(cTimer.SetTimer("Calc coefs time (U mom)")),
	cdd(deltaSize),
	cdA(cdd * 1),
	cdV(cdd* cdd* 1),
	cNX(int64_t(L / cdd) + 2),
	cNZ(int64_t(H / cdd) + 2),
	cNXCylinder(int64_t(0.5 + (l / cdd))),
	cNZCylinder(cNXCylinder),
	ciCyInit((cNXCylinder * 8) - cNXCylinder),
	ckCyInit((cNZ / 2) - (cNZCylinder / 2)),
	ciCyEnd(cNXCylinder * 8),
	ckCyEnd((cNZ / 2) + (cNZCylinder / 2)),
	cRHO(rho),
	cMU(mu),
	cD(cdA * mu / cdd)
{
	cPressureField.resize(cNZ * cNX, 0.0);

	cPressureCorrField.resize(cNZ * cNX, 0.0);

	cUVelField.resize(cNZ * cNX, Ufarfield);
	cWVelField.resize(cNZ * cNX, Wfarfield);

	cSpUField.resize(cNZ * cNX, 0.0);
	cSpWField.resize(cNZ * cNX, 0.0);

	cUHatVelField = cPressureCorrField;
	cWHatVelField = cPressureCorrField;

	caijkU = cPressureCorrField;
	caijkW = cPressureCorrField;
}

FieldRectCylinder2D::~FieldRectCylinder2D()
{}

FieldRectCylinder2D::Intervalo::Intervalo(int64_t v1, int64_t v2)
	:
	valor1(v1),
	valor2(v2)
{}

int64_t FieldRectCylinder2D::Intervalo::Tamanho() const
{
	return valor2 - valor1;
}

void FieldRectCylinder2D::SetBCFlatPlate()
{
	// SET INITIAL VALUES
	SetUValueInterval({ ciCyInit, ciCyEnd + 1 },{ ckCyInit, ckCyEnd }, 0.0);
	SetWValueInterval({ ciCyInit, ciCyEnd },{ ckCyInit, ckCyEnd + 1 }, 0.0);
	SetPValueInterval({ ciCyInit, ciCyEnd },{ ckCyInit, ckCyEnd }, 0.0);

	// INLET; OUTLET TUDO OK
	// VALORES PADROES

	// WALL (1)
	SetSpWValueInterval({ ciCyInit - 1, ciCyInit },{ ckCyInit, ckCyEnd + 1 }, -cMU * cdA / (cdd / 2.0));

	// WALL (2)
	SetSpWValueInterval({ ciCyEnd, ciCyEnd + 1 },{ ckCyInit, ckCyEnd + 1 }, -cMU * cdA / (cdd / 2.0));

	// WALL (4)
	SetSpUValueInterval({ ciCyInit, ciCyEnd + 1 },{ ckCyInit - 1, ckCyInit }, -cMU * cdA / (cdd / 2.0));
	
	// WALL (3)
	SetSpUValueInterval({ ciCyInit, ciCyEnd + 1 },{ ckCyEnd, ckCyEnd + 1 }, -cMU * cdA / (cdd / 2.0));
}

int FieldRectCylinder2D::CreateUMomentumLSCSR(const doubleField2D& uField0, double dt,
	std::vector<int>& ptr, std::vector<int>& col, std::vector<double>& val, std::vector<double>& rhs)
{
	double totalTimeCalcs = 0.0;
	double totalTimePushs = 0.0;

	const int n2 = static_cast<int64_t>(cNX * cNZ);        // Number of points in the grid.

	ptr.clear();
	if (ptr.capacity() != n2 + 1)
		ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	if (col.capacity() != n2 * 5)
		col.reserve(n2 * 5);

	val.clear();
	if (val.capacity() != n2 * 5)
		val.reserve(n2 * 5);

	rhs.clear();
	rhs.resize(n2);

	for (int k = 0; k < cNZ; k++)
	{
		for (int i = 0; i < cNX; i++)
		{
			const int64_t index = static_cast<int64_t>(i + (k * cNX));

			if (i == 1 || i == 0 || i == cNX - 1 || k == 0 || k == cNZ - 1)
			{
				// Boundary point. Use Dirichlet condition. (seta o valor da variavel diretamente)

				col.push_back(index);
				val.push_back(1.0);

				rhs[index] = cUVelField[ (k * cNX) + i];
			}
			else if (i >= ciCyInit && i < ciCyEnd + 1 && k >= ckCyInit && k < ckCyEnd)
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
				col.push_back(index + static_cast<int64_t>((0 * cNX)));
				val.push_back(NAN);

				// jkmi
				col.push_back(index + static_cast<int64_t>((-1 * cNX)));
				val.push_back(-NAN);

				// jkim
				col.push_back(index - static_cast<int64_t>(1 + (0 * cNX)));
				val.push_back(-NAN);

				// jkip
				col.push_back(index + 1 + (0 * cNX));
				val.push_back(-NAN);

				// jkpi
				col.push_back(index + static_cast<int64_t>((1 * cNX)));
				val.push_back(-NAN);

				// jmki
				col.push_back(index + static_cast<int64_t>((0 * cNX)));
				val.push_back(-NAN);

				// jpki
				col.push_back(index + static_cast<int64_t>((0 * cNX)));
				val.push_back(-NAN);

				rhs[index] = NAN;

				cTimer.Tock(ctidPushback);

				totalTimePushs += cTimer.GetLastTickTock(ctidPushback);
			}

			ptr.push_back(static_cast<int64_t>(col.size()));
		}
	}

	cTimer.Tick(ctidCalcCoefs);

#pragma omp parallel for
	for (int k = 0; k < cNZ; k++)
	{
		for (int i = 0; i < cNX; i++)
		{
			const int64_t index = static_cast<int64_t>(i + (k * cNX));

			if (i == 1 || i == 0 || i == cNX - 1 || k == 0 || k == cNZ - 1)
			{
				// do nothing, already did it before;
			}
			else if (i >= ciCyInit && i < ciCyEnd + 1 && k >= ckCyInit && k < ckCyEnd)
			{
				// do nothing, already did it before;
			}
			else
			{
				CoefsData2D aijkU;

				const double Fw = GetFwInternalUMomentum(i, k);
				const double Fe = GetFeInternalUMomentum(i, k);
				const double Fb = GetFbInternalUMomentum(i, k);
				const double Ft = GetFtInternalUMomentum(i, k);

				const double deltaf = cdA * (Fe - Fw + Ft - Fb);

				aijkU.aimjk = cD + (cdA*fmax(Fw, 0.0));
				aijkU.aipjk = cD + (cdA*fmax(-Fe, 0.0));
				
				aijkU.aijkm = cD + (cdA*fmax(Fb, 0.0));
				aijkU.aijkp = cD + (cdA*fmax(-Ft, 0.0));

				const double alfaw = (Fw > 0.0) ? 1.0 : 0.0;
				const double alfae = (Fe > 0.0) ? 1.0 : 0.0;
				const double alfat = (Ft > 0.0) ? 1.0 : 0.0;
				const double alfab = (Fb > 0.0) ? 1.0 : 0.0;

				const double reneg = ren(i, k, cUVelField);
				const double rwneg = rwn(i, k, cUVelField);
				const double rtneg = rtn(i, k, cUVelField);
				const double rbneg = rbn(i, k, cUVelField);

				const double repos = rep(i, k, cUVelField);
				const double rwpos = rwp(i, k, true, cUVelField);
				const double rtpos = rtp(i, k, cUVelField);
				const double rbpos = rbp(i, k, false, cUVelField);

				const double phiE = cUVelField[(k * cNX) + i + 1];
				const double phiP = cUVelField[index];
				const double phiW = cUVelField[(k * cNX) + i - 1];
				const double phiT = cUVelField[((k + 1) * cNX) + i];
				const double phiB = cUVelField[((k - 1) * cNX) + i];

				const double SuDc = cdA * 0.5*(
					(Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
					+ (Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
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

				aijkU.aijk = aijkU.aimjk + aijkU.aipjk +
					aijkU.aijkm + aijkU.aijkp + deltaf - cSpUField[index]
					+ (cRHO * cdV / dt);

				caijkU[ (k * cNX) + i] = aijkU.aijk;


				aijkU.source = SuDc +
					((cPressureField[ (k * cNX) + i - 1] - cPressureField[ (k * cNX) + i]) * cdA) +
					/*cSuUField[ (k * cNX) + i] +*/
					((cRHO * cdd * cdd / dt) * uField0[ (k * cNX) + i]);


				// jki
				val[ptr[index]] = aijkU.aijk;

				// jkmi
				val[ptr[index] + 1] = -aijkU.aijkm;

				// jkim
				val[ptr[index] + 2] = -aijkU.aimjk;

				// jkip
				val[ptr[index] + 3] = -aijkU.aipjk;

				// jkpi
				val[ptr[index] + 4] = -aijkU.aijkp;

				rhs[index] = aijkU.source;
			}
		}
	}

	cTimer.Tock(ctidCalcCoefs);

	//cTimer.WriteToCoutTickTock(ctidCalcCoefs);
	//   
	//std::cout << "# Timer total time calcs " << totalTimeCalcs << " total time pushs " << totalTimePushs << std::endl;

	return n2;
}


int FieldRectCylinder2D::CreateWMomentumLSCSRParallel(const doubleField2D& wField0, double dt, std::vector<int>& ptr,
	std::vector<int>& col, std::vector<double>& val, std::vector<double>& rhs)
{
	const int n2 = cNX * cNZ; // Number of points in the grid.

	ptr.clear();
	if (ptr.capacity() != n2 + 1)
		ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	if (col.capacity() != n2 * 5)
		col.reserve(n2 * 5);

	val.clear();
	if (val.capacity() != n2 * 5)
		val.reserve(n2 * 5);

	rhs.clear();
	rhs.resize(n2);

	for (int k = 0, index = 0; k < cNZ; k++)
	{
		for (int i = 0; i < cNX; i++, index++) {

			if (i == 0 || i == cNX - 1 || k == 0 || k == 1 || k == cNZ - 1)
			{
				// Boundary point. Use Dirichlet condition. (seta o valor da variavel diretamente)

				col.push_back(index);
				val.push_back(1.0);

				rhs[index] = cWVelField[(k * cNX) + i];
			}
			else if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd + 1)
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
				col.push_back(index + (0 * cNX));
				val.push_back(NAN);

				// jkmi
				col.push_back(index + (-1 * cNX));
				val.push_back(-NAN);

				// jkim
				col.push_back(index - 1 + (0 * cNX));
				val.push_back(-NAN);

				// jkip
				col.push_back(index + 1 + (0 * cNX));
				val.push_back(-NAN);

				// jkpi
				col.push_back(index + (1 * cNX));
				val.push_back(-NAN);

				//rhs[index] = aijkW.source;
			}

			ptr.push_back(col.size());
		}
	}


#pragma omp parallel for
	for (int k = 0; k < cNZ; k++)
	{
		for (int i = 0; i < cNX; i++)
		{
			const int index = i + (k * cNX) ;

			if (i == 0 || i == cNX - 1 || k == 0 || k == 1 || k == cNZ - 1)
			{
				// do nothing, already did it
			}
			else if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd + 1)
			{
				// do nothing, already did it
			}
			else
			{
				CoefsData2D aijkW;

				const double Fw = GetFwInternalWMomentum(i, k) * cdA;
				const double Fe = GetFeInternalWMomentum(i, k) * cdA;
				const double Fb = GetFbInternalWMomentum(i, k) * cdA;
				const double Ft = GetFtInternalWMomentum(i, k) * cdA;

				const double deltaf = Fe - Fw + Ft - Fb;

				aijkW.aimjk = cD + fmax(Fw, 0.0);
				aijkW.aipjk = cD + fmax(-Fe, 0.0);
				
				aijkW.aijkm = cD + fmax(Fb, 0.0);
				aijkW.aijkp = cD + fmax(-Ft, 0.0);

				const double alfaw = (Fw > 0) ? 1.0 : 0.0;
				const double alfae = (Fe > 0) ? 1.0 : 0.0;
				const double alfab = (Fb > 0) ? 1.0 : 0.0;
				const double alfat = (Ft > 0) ? 1.0 : 0.0;

				const double reneg = ren(i, k, cWVelField);
				const double rwneg = rwn(i, k, cWVelField);
				const double rtneg = rtn(i, k, cWVelField);
				const double rbneg = rbn(i, k, cWVelField);

				const double repos = rep(i, k, cWVelField);
				const double rwpos = rwp(i, k, false, cWVelField);
				const double rtpos = rtp(i, k, cWVelField);
				const double rbpos = rbp(i, k, true, cWVelField);

				const double phiE = cWVelField[(k * cNX) + i + 1];
				const double phiP = cWVelField[(k * cNX) + i];
				const double phiW = cWVelField[(k * cNX) + i - 1];
				const double phiT = cWVelField[((k + 1) * cNX) + i];
				const double phiB = cWVelField[((k - 1) * cNX) + i];

				const double SuDc = (0.5 * Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
					+ (0.5 * Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
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

				aijkW.aijk = aijkW.aimjk + aijkW.aipjk + aijkW.aijkm + 
					aijkW.aijkp + deltaf - cSpWField[(k * cNX) + i]
					+ (cRHO * cdV / dt);

				caijkW[(k * cNX) + i] = aijkW.aijk;

				aijkW.source = SuDc + ((cPressureField[((k - 1) * cNX) + i] - cPressureField[(k * cNX) + i]) * cdA) +
					/*cSuWField[ (k * cNX) + i] +*/ ((cRHO * cdd * cdd / dt) * wField0[(k * cNX) + i]);

				// [i + (k * cNX) ] = [index global]

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

				rhs[index] = aijkW.source;
			}
		}
	}

	return n2;
}
int FieldRectCylinder2D::CreatePressureCorrectionLSCSR(std::vector<int>& ptr, std::vector<int>& col, std::vector<double>& val,
	std::vector<double>& rhs) const
{
	const int n2 = cNX * cNZ; // Number of points in the grid.
	/*
	if (ptr.size() != n2 + 1)
	{
		ptr.clear();
		ptr.reserve(n2 + 1);
		ptr.push_back(0);
	}

	if (col.size() != n2 * 5)
	{
		col.clear();
		col.reserve(n2 * 5); // sete coeficientes
	}

	if (val.size() != n2 * 5)
	{
		val.clear();
		val.reserve(n2 * 5);  // sete coeficientes
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
	col.reserve(n2 * 5);

	val.clear();
	val.reserve(n2 * 5);

	rhs.clear();
	rhs.resize(n2);
	for (int k = 0, index = 0; k < cNZ; k++)
	{
		for (int i = 0; i < cNX; i++, index++) {

			if (i == 0 || i == cNX - 1 || k == 0 || k == cNZ - 1)
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
			else if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd)
			{
				// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA PCORR = 0 )
				col.push_back(index);
				val.push_back(1.0);

				rhs[index] = 0.0;
			}
			else
			{
				CoefsData2D aijkPCorr;

				aijkPCorr.aimjk = cRHO * cdA * cdA / caijkU[(k * cNX) + i];
				aijkPCorr.aipjk = cRHO * cdA * cdA / caijkU[(k * cNX) + i + 1];
				
				aijkPCorr.aijkm = cRHO * cdA * cdA / caijkW[(k * cNX) + i];
				aijkPCorr.aijkp = cRHO * cdA * cdA / caijkW[((k + 1) * cNX) + i];

				if (i == 1)
					aijkPCorr.aimjk = 0.0;
				else if (i == cNX - 2)
					aijkPCorr.aipjk = 0.0;
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
					aijkPCorr.aijkm + aijkPCorr.aijkp);

				aijkPCorr.source = cRHO * cdA * (cUVelField[(k * cNX) + i] - cUVelField[(k * cNX) + i + 1] +
					cWVelField[(k * cNX) + i] - cWVelField[((k + 1) * cNX) + i]);

				// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]

				// jki
				col.push_back(index + (0 * cNX));
				val.push_back(aijkPCorr.aijk);

				// jkmi
				col.push_back(index + (-1 * cNX));
				val.push_back(-aijkPCorr.aijkm);

				// jkim
				col.push_back(index - 1 + (0 * cNX));
				val.push_back(-aijkPCorr.aimjk);

				// jkip
				col.push_back(index + 1 + (0 * cNX));
				val.push_back(-aijkPCorr.aipjk);

				// jkpi
				col.push_back(index + (1 * cNX));
				val.push_back(-aijkPCorr.aijkp);
				
				rhs[index] = aijkPCorr.source;
			}

			ptr.push_back(col.size());
		}
	}

	return n2;
}

int FieldRectCylinder2D::CreatePressureLSCSRFI(std::vector<int>& ptr, std::vector<int>& col, std::vector<double>& val,
	std::vector<double>& rhs) const
{
	int n2 = cNX * cNZ; // total number of points in the grid.

	ptr.clear();
	ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	col.reserve(n2 * 5);

	val.clear();
	val.reserve(n2 * 5);

	rhs.clear();
	rhs.resize(n2);

	for (int k = 0; k < cNZ; k++)
	{
		for (int i = 0; i < cNX; i++) {

			const int64_t index = i + (k * cNX);

			if (i == 0 || i == cNX - 1 || k == 0 || k == cNZ - 1)
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
			else if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd)
			{
				// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA P = 0 )
				col.push_back(index);
				val.push_back(1.0);

				rhs[index] = 0.0;
			}
			else
			{
				CoefsData2D aijkP;

				aijkP.aimjk = cRHO * cdA * cdA / caijkU[(k * cNX) + i];
				aijkP.aipjk = cRHO * cdA * cdA / caijkU[(k * cNX) + i + 1];
				
				aijkP.aijkm = cRHO * cdA * cdA / caijkW[(k * cNX) + i];
				aijkP.aijkp = cRHO * cdA * cdA / caijkW[((k + 1) * cNX) + i];

				if (i == 1)
					aijkP.aimjk = 0.0;
				else if (i == cNX - 2)
					aijkP.aipjk = 0.0;
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
					aijkP.aijkm + aijkP.aijkp);

				aijkP.source = cRHO * cdA * (cUHatVelField[(k * cNX) + i] - cUHatVelField[(k * cNX) + i + 1] +
					cWHatVelField[(k * cNX) + i] - cWHatVelField[((k + 1) * cNX) + i]);

				// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]

				// jki
				col.push_back(index + (0 * cNX));
				val.push_back(aijkP.aijk);

				// jkmi
				col.push_back(index + (-1 * cNX));
				val.push_back(-aijkP.aijkm);

				// jkim
				col.push_back(index - 1 + (0 * cNX));
				val.push_back(-aijkP.aimjk);

				// jkip
				col.push_back(index + 1 + (0 * cNX));
				val.push_back(-aijkP.aipjk);

				// jkpi
				col.push_back(index + (1 * cNX));
				val.push_back(-aijkP.aijkp);
				
				rhs[index] = aijkP.source;
			}

			ptr.push_back(col.size());
		}
	}

	return n2;
}

void FieldRectCylinder2D::CalcUHatValuesFI(const doubleField2D& uField0, double dt)
{
	// COPIAR TODOS OS VALORES DA CONDIÇÃO DE CONTORNO ( OBSERVE QUE DEPENDE DA GEOMETRIA SE QUISER SER MAIS EFICIENTE )
	//cUHatVelField = cUVelField;

	memcpy(cUHatVelField.data(), cUVelField.data(), cUVelField.size() * sizeof(double));

#pragma omp parallel for
	for (int k = 1; k < cNZ - 1; k++)
	{
		for (int i = 2; i < cNX - 1; i++)
		{
			const int index = i + (k * cNX);

			if (i >= ciCyInit && i < ciCyEnd + 1 && k >= ckCyInit && k < ckCyEnd)
			{
				// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA U = 0 )
				cUHatVelField[index] = 0.0;
			}
			else
			{
				CoefsData2D aijkU;

				const double Fw = GetFwInternalUMomentum(i, k);
				const double Fe = GetFeInternalUMomentum(i, k);
				const double Fb = GetFbInternalUMomentum(i, k);
				const double Ft = GetFtInternalUMomentum(i, k);

				const double deltaf = cdA * (Fe - Fw + Ft - Fb);

				aijkU.aimjk = cD + (cdA * fmax(Fw, 0.0));
				aijkU.aipjk = cD + (cdA * fmax(-Fe, 0.0));
				
				aijkU.aijkm = cD + (cdA * fmax(Fb, 0.0));
				aijkU.aijkp = cD + (cdA * fmax(-Ft, 0.0));

				const double alfaw = (Fw > 0.0) ? 1.0 : 0.0;
				const double alfae = (Fe > 0.0) ? 1.0 : 0.0;
				const double alfat = (Ft > 0.0) ? 1.0 : 0.0;
				const double alfab = (Fb > 0.0) ? 1.0 : 0.0;

				const double reneg = ren(i, k, cUVelField);
				const double rwneg = rwn(i, k, cUVelField);
				const double rtneg = rtn(i, k, cUVelField);
				const double rbneg = rbn(i, k, cUVelField);

				const double repos = rep(i, k, cUVelField);
				const double rwpos = rwp(i, k, true, cUVelField);
				const double rtpos = rtp(i, k, cUVelField);
				const double rbpos = rbp(i, k, false, cUVelField);

				const double phiE = cUVelField[(k * cNX) + i + 1];
				const double phiP = cUVelField[index];
				const double phiW = cUVelField[(k * cNX) + i - 1];
				const double phiT = cUVelField[((k + 1) * cNX) + i];
				const double phiB = cUVelField[((k - 1) * cNX) + i];

				const double SuDc = cdA * 0.5 * (
					(Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
					+ (Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
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

				aijkU.aijk = aijkU.aimjk + aijkU.aipjk + aijkU.aijkm + 
					aijkU.aijkp + deltaf - cSpUField[index]
					+ (cRHO * cdV / dt);

				caijkU[index] = aijkU.aijk;

				// sem pressão (hat bruh)
				aijkU.source = SuDc +
					/*cSuUField[ (k * cNX) + i] +*/
					((cRHO * cdd * cdd / dt) * uField0[index]);

				cUHatVelField[index] = ((aijkU.aimjk * phiW) + (aijkU.aipjk * phiE)
					+ (aijkU.aijkm * phiB)
					+ (aijkU.aijkp * phiT) + aijkU.source) / aijkU.aijk;
			}
		}
	}


}
void FieldRectCylinder2D::CalcWHatValuesFI(const doubleField2D& wField0, double dt)
{
	// COPIAR TODOS OS VALORES DA CONDIÇÃO DE CONTORNO ( OBSERVE QUE DEPENDE DA GEOMETRIA SE QUISER SER MAIS EFICIENTE )
	//cWHatVelField = cWVelField;

	memcpy(cWHatVelField.data(), cWVelField.data(), cWVelField.size() * sizeof(double));

#pragma omp parallel for
	for (int k = 2; k < cNZ - 1; k++)
	{
		for (int i = 1; i < cNX - 1; i++)
		{
			const int index = i + (k*cNX);

			if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd + 1)
			{
				// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA V = 0 )
				cWHatVelField[index] = 0.0;
			}
			else
			{
				CoefsData2D aijkW;

				const double Fw = GetFwInternalWMomentum(i, k) * cdA;
				const double Fe = GetFeInternalWMomentum(i, k) * cdA;
				const double Fb = GetFbInternalWMomentum(i, k) * cdA;
				const double Ft = GetFtInternalWMomentum(i, k) * cdA;

				const double deltaf = Fe - Fw + Ft - Fb;

				aijkW.aimjk = cD + fmax(Fw, 0.0);
				aijkW.aipjk = cD + fmax(-Fe, 0.0);
				
				aijkW.aijkm = cD + fmax(Fb, 0.0);
				aijkW.aijkp = cD + fmax(-Ft, 0.0);

				const double alfaw = (Fw > 0) ? 1.0 : 0.0;
				const double alfae = (Fe > 0) ? 1.0 : 0.0;
				const double alfab = (Fb > 0) ? 1.0 : 0.0;
				const double alfat = (Ft > 0) ? 1.0 : 0.0;

				const double reneg = ren(i, k, cWVelField);
				const double rwneg = rwn(i, k, cWVelField);
				const double rtneg = rtn(i, k, cWVelField);
				const double rbneg = rbn(i, k, cWVelField);

				const double repos = rep(i, k, cWVelField);
				const double rwpos = rwp(i, k, false, cWVelField);
				const double rtpos = rtp(i, k, cWVelField);
				const double rbpos = rbp(i, k, true, cWVelField);

				const double phiE = cWVelField[(k * cNX) + i + 1];
				const double phiP = cWVelField[(k * cNX) + i];
				const double phiW = cWVelField[(k * cNX) + i - 1];
				const double phiT = cWVelField[((k + 1) * cNX) + i];
				const double phiB = cWVelField[((k - 1) * cNX) + i];

				const double SuDc = (0.5 * Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
					+ (0.5 * Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
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


				aijkW.aijk = aijkW.aimjk + aijkW.aipjk + 
					aijkW.aijkm + aijkW.aijkp + deltaf - cSpWField[(k * cNX) + i]
					+ (cRHO * cdV / dt);

				caijkW[(k * cNX) + i] = aijkW.aijk;

				// no pressure comp. (hat bruh)
				aijkW.source = SuDc +
					((cRHO * cdd * cdd / dt) * wField0[(k * cNX) + i]);

				cWHatVelField[index] = ((aijkW.aimjk * phiW) + (aijkW.aipjk * phiE) +
					(aijkW.aijkm * phiB) + (aijkW.aijkp * phiT) + aijkW.source) / aijkW.aijk;
			}
		}
	}

}

void FieldRectCylinder2D::CalcUHatValuesCN(const doubleField2D& uField0, double dt)
{
	// COPIAR TODOS OS VALORES DA CONDIÇÃO DE CONTORNO ( OBSERVE QUE DEPENDE DA GEOMETRIA SE QUISER SER MAIS EFICIENTE )
	//cUHatVelField = cUVelField;

	memcpy(cUHatVelField.data(), cUVelField.data(), cUVelField.size() * sizeof(double));

#pragma omp parallel for
	for (int k = 1; k < cNZ - 1; k++)
	{
		for (int i = 2; i < cNX - 1; i++)
		{
			const int index = i + (k * cNX);

			if (i >= ciCyInit && i < ciCyEnd + 1 && k >= ckCyInit && k < ckCyEnd)
			{
				// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA U = 0 )
				cUHatVelField[index] = 0.0;
			}
			else
			{
				CoefsData2D aijkU;

				const double Fw = GetFwInternalUMomentum(i, k);
				const double Fe = GetFeInternalUMomentum(i, k);
				const double Fb = GetFbInternalUMomentum(i, k);
				const double Ft = GetFtInternalUMomentum(i, k);

				const double deltaf = cdA * (Fe - Fw + Ft - Fb);

				aijkU.aimjk = cD + (cdA * fmax(Fw, 0.0));
				aijkU.aipjk = cD + (cdA * fmax(-Fe, 0.0));
				
				aijkU.aijkm = cD + (cdA * fmax(Fb, 0.0));
				aijkU.aijkp = cD + (cdA * fmax(-Ft, 0.0));

				const double alfaw = (Fw > 0.0) ? 1.0 : 0.0;
				const double alfae = (Fe > 0.0) ? 1.0 : 0.0;
				const double alfat = (Ft > 0.0) ? 1.0 : 0.0;
				const double alfab = (Fb > 0.0) ? 1.0 : 0.0;

				const double reneg = ren(i, k, cUVelField);
				const double rwneg = rwn(i, k, cUVelField);
				const double rtneg = rtn(i, k, cUVelField);
				const double rbneg = rbn(i, k, cUVelField);

				const double repos = rep(i, k, cUVelField);
				const double rwpos = rwp(i, k, true, cUVelField);
				const double rtpos = rtp(i, k, cUVelField);
				const double rbpos = rbp(i, k, false, cUVelField);

				const double phiE = cUVelField[(k * cNX) + i + 1];
				const double phiP = cUVelField[index];
				const double phiW = cUVelField[(k * cNX) + i - 1];
				const double phiT = cUVelField[((k + 1) * cNX) + i];
				const double phiB = cUVelField[((k - 1) * cNX) + i];

				const double phi0E = uField0[(k * cNX) + i + 1];
				const double phi0P = uField0[index];
				const double phi0W = uField0[(k * cNX) + i - 1];
				const double phi0T = uField0[((k + 1) * cNX) + i];
				const double phi0B = uField0[((k - 1) * cNX) + i];

				const double SuDc = cdA * 0.5 * (
					(Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
					+ (Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
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

				aijkU.aijk = 0.5*(aijkU.aimjk + aijkU.aipjk + 
					aijkU.aijkm + aijkU.aijkp) + deltaf - (0.5*cSpUField[index])
					+ (cRHO * cdV / dt);

				caijkU[index] = aijkU.aijk;

				// sem pressão (hat bruh)
				aijkU.source = SuDc + (0.5 * cSpUField[(k * cNX) + i] * phi0P) +
					/*cSuUField[ (k * cNX) + i] +*/
					(((cRHO * cdd * cdd / dt)
						- (aijkU.aipjk / 2.0) - (aijkU.aimjk / 2.0)
						- (aijkU.aijkp / 2.0) - (aijkU.aijkm / 2.0)) * phi0P);

				cUHatVelField[index] = ((aijkU.aimjk * 0.5*(phiW + phi0W)) + (aijkU.aipjk * 0.5*(phiE + phi0E))
					+ (aijkU.aijkm * 0.5*(phiB + phi0B))
					+ (aijkU.aijkp * 0.5*(phiT + phi0T)) + aijkU.source) / aijkU.aijk;
			}
		}
	}
}
void FieldRectCylinder2D::CalcWHatValuesCN(const doubleField2D& wField0, double dt)
{
	// COPIAR TODOS OS VALORES DA CONDIÇÃO DE CONTORNO ( OBSERVE QUE DEPENDE DA GEOMETRIA SE QUISER SER MAIS EFICIENTE )
	//cWHatVelField = cWVelField;

	memcpy(cWHatVelField.data(), cWVelField.data(), cWVelField.size() * sizeof(double));

#pragma omp parallel for
	for (int k = 2; k < cNZ - 1; k++)
	{
		for (int i = 1; i < cNX - 1; i++)
		{
			const int index = i + (k * cNX);

			if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd + 1)
			{
				// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA V = 0 )
				cWHatVelField[index] = 0.0;
			}
			else
			{
				CoefsData2D aijkW;

				const double Fw = GetFwInternalWMomentum(i, k) * cdA;
				const double Fe = GetFeInternalWMomentum(i, k) * cdA;
				const double Fb = GetFbInternalWMomentum(i, k) * cdA;
				const double Ft = GetFtInternalWMomentum(i, k) * cdA;

				const double deltaf = Fe - Fw + Ft - Fb;

				aijkW.aimjk = cD + fmax(Fw, 0.0);
				aijkW.aipjk = cD + fmax(-Fe, 0.0);
				
				aijkW.aijkm = cD + fmax(Fb, 0.0);
				aijkW.aijkp = cD + fmax(-Ft, 0.0);

				const double alfaw = (Fw > 0) ? 1.0 : 0.0;
				const double alfae = (Fe > 0) ? 1.0 : 0.0;
				const double alfab = (Fb > 0) ? 1.0 : 0.0;
				const double alfat = (Ft > 0) ? 1.0 : 0.0;

				const double reneg = ren(i, k, cWVelField);
				const double rwneg = rwn(i, k, cWVelField);
				const double rtneg = rtn(i, k, cWVelField);
				const double rbneg = rbn(i, k, cWVelField);

				const double repos = rep(i, k, cWVelField);
				const double rwpos = rwp(i, k, false, cWVelField);
				const double rtpos = rtp(i, k, cWVelField);
				const double rbpos = rbp(i, k, true, cWVelField);

				const double phiE = cWVelField[(k * cNX) + i + 1];
				const double phiP = cWVelField[(k * cNX) + i];
				const double phiW = cWVelField[(k * cNX) + i - 1];
				const double phiT = cWVelField[((k + 1) * cNX) + i];
				const double phiB = cWVelField[((k - 1) * cNX) + i];

				const double phi0E = wField0[(k * cNX) + i + 1];
				const double phi0P = wField0[(k * cNX) + i];
				const double phi0W = wField0[(k * cNX) + i - 1];
				const double phi0T = wField0[((k + 1) * cNX) + i];
				const double phi0B = wField0[((k - 1) * cNX) + i];

				const double SuDc = (0.5 * Fe * (((1.0 - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
					+ (0.5 * Fw * (-((1.0 - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
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


				aijkW.aijk = 0.5*(aijkW.aimjk + aijkW.aipjk + 
					aijkW.aijkm + aijkW.aijkp) + deltaf - (0.5*cSpWField[(k * cNX) + i])
					+ (cRHO * cdV / dt);

				caijkW[(k * cNX) + i] = aijkW.aijk;

				// no pressure comp. (hat bruh)
				aijkW.source = SuDc +
					(0.5 * cSpWField[(k * cNX) + i] * phi0P) +
					(((cRHO * cdd * cdd / dt)
						- (aijkW.aipjk / 2.0) - (aijkW.aimjk / 2.0)
						- (aijkW.aijkp / 2.0) - (aijkW.aijkm / 2.0)) * phi0P);

				cWHatVelField[index] = ((aijkW.aimjk * 0.5*(phiW + phi0W)) + (aijkW.aipjk * 0.5*(phiE + phi0E)) +
					(aijkW.aijkm * 0.5*(phiB + phi0B)) +
					(aijkW.aijkp * 0.5*(phiT + phi0T)) + aijkW.source) / aijkW.aijk;
			}
		}
	}
}

// Set P value for all 'K' and 'I' in J defined; 
void FieldRectCylinder2D::SetPValueIK(double pValue)
{
	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			cPressureField[ (K * cNX) + I] = pValue;
		}
	}
}

// Set U value for all 'K' and 'i' in J defined; 
void FieldRectCylinder2D::SetUValueiK(double uValue)
{
	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t i = 0; i < cNX; i++)
		{
			cUVelField[ (K * cNX) + i] = uValue;
		}
	}
}

// Set W value for all 'k' and 'I' in J defined; 
void FieldRectCylinder2D::SetWValueIk(double wValue)
{
	for (int64_t k = 0; k < cNZ; k++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			cWVelField[ (k * cNX) + I] = wValue;
		}
	}
}

void FieldRectCylinder2D::SetPValueInterval(const Intervalo& II,const Intervalo& KK, double pValue)
{
	for (int64_t k = KK.valor1; k < KK.valor2; k++)
	{
		for (int64_t i = II.valor1; i < II.valor2; i++)
		{
			cPressureField[(k * cNX) + i] = pValue;
		}
	}

}

void FieldRectCylinder2D::SetUValueInterval(const Intervalo& ii,const Intervalo& KK, double uValue)
{
	for (int64_t k = KK.valor1; k < KK.valor2; k++)
	{
		for (int64_t i = ii.valor1; i < ii.valor2; i++)
		{
			cUVelField[(k * cNX) + i] = uValue;
		}
	}
}

void FieldRectCylinder2D::SetWValueInterval(const Intervalo& II,const Intervalo& kk, double wValue)
{
	for (int64_t k = kk.valor1; k < kk.valor2; k++)
	{
		for (int64_t i = II.valor1; i < II.valor2; i++)
		{
			cWVelField[(k * cNX) + i] = wValue;
		}
	}
}

void FieldRectCylinder2D::ExtrapolateForwardPJK(int64_t I)
{
	for (int64_t k = 0; k < cNZ; k++)
	{
		cPressureField[(k * cNX) + I] = cPressureField[(k * cNX) + I + 1];
	}
}

void FieldRectCylinder2D::ExtrapolateForwardUJK(int64_t i)
{
	for (int64_t k = 0; k < cNZ; k++)
	{
		cUVelField[(k * cNX) + i] = cUVelField[(k * cNX) + i + 1];
	}
}

void FieldRectCylinder2D::ExtrapolateForwardWJk(int64_t I)
{
	for (int64_t k = 0; k < cNZ; k++)
	{
		cWVelField[(k * cNX) + I] = cWVelField[(k * cNX) + I + 1];
	}
}

void FieldRectCylinder2D::ExtrapolateBackwardPJK(int64_t I)
{
	for (int64_t k = 1; k < cNZ - 1; k++)
	{
		cPressureField[(k * cNX) + I] = cPressureField[(k * cNX) + I - 1];
	}
}

void FieldRectCylinder2D::ExtrapolateBackwardUJK(int64_t i)
{
	for (int64_t k = 0; k < cNZ; k++)
	{
		cUVelField[(k * cNX) + i] = cUVelField[(k * cNX) + i - 1];
	}
}

void FieldRectCylinder2D::ExtrapolateBackwardWJk(int64_t I)
{
	for (int64_t k = 0; k < cNZ; k++)
	{
		cWVelField[(k * cNX) + I] = cWVelField[(k * cNX) + I - 1];
	}
}

void FieldRectCylinder2D::ExtrapolateBackwardWeightedUJK(int64_t i, double weight)
{
	for (int64_t k = 1; k < cNZ - 1; k++)
	{
		cUVelField[(k * cNX) + i] = cUVelField[(k * cNX) + i - 1] * weight;
	}
}

void FieldRectCylinder2D::ExtrapolateBackwardWeightedWJk(int64_t I, double weight)
{
	for (int64_t k = 0; k < cNZ; k++)
	{
		cWVelField[(k * cNX) + I] = cWVelField[(k * cNX) + I - 1] * weight;
	}
}

void FieldRectCylinder2D::ExtrapolateForwardIntervalPJK(int64_t I,const Intervalo& KK)
{
	for (int64_t k = KK.valor1; k < KK.valor2; k++)
	{
		cPressureField[(k * cNX) + I] = cPressureField[(k * cNX) + I + 1];
	}
}

void FieldRectCylinder2D::ExtrapolateForwardIntervalUJK(int64_t i,const Intervalo& KK)
{
	for (int64_t k = KK.valor1; k < KK.valor2; k++)
	{
		cUVelField[(k * cNX) + i] = cUVelField[(k * cNX) + i + 1];
	}
}

void FieldRectCylinder2D::ExtrapolateForwardIntervalWJk(int64_t I,const Intervalo& kk)
{
	for (int64_t k = kk.valor1; k < kk.valor2; k++)
	{
		cWVelField[(k * cNX) + I] = cWVelField[(k * cNX) + I + 1];
	}
}

void FieldRectCylinder2D::ExtrapolateBackwardIntervalPJK(int64_t I,const Intervalo& KK)
{
	for (int64_t k = KK.valor1; k < KK.valor2; k++)
	{
		cPressureField[(k * cNX) + I] = cPressureField[(k * cNX) + I - 1];
	}
}

void FieldRectCylinder2D::ExtrapolateBackwardIntervalUJK(int64_t i,const Intervalo& KK)
{
	for (int64_t k = KK.valor1; k < KK.valor2; k++)
	{
		cUVelField[(k * cNX) + i] = cUVelField[(k * cNX) + i - 1];
	}
}

void FieldRectCylinder2D::ExtrapolateBackwardIntervalWJk(int64_t I,const Intervalo& kk)
{
	for (int64_t k = kk.valor1; k < kk.valor2; k++)
	{
		cWVelField[(k * cNX) + I] = cWVelField[(k * cNX) + I - 1];
	}
}

void FieldRectCylinder2D::SaveField(const std::string& baseFileName) const
{
	std::ofstream uFile(baseFileName + "U.dat", std::ios::binary);

	uFile.write(reinterpret_cast<const char*>(&cNX), sizeof(int64_t));
	uFile.write(reinterpret_cast<const char*>(&cNZ), sizeof(int64_t));

	uFile.write(reinterpret_cast<const char*>(cUVelField.data()), cUVelField.size() * sizeof(double));

	uFile.close();
	
	std::ofstream wFile(baseFileName + "W.dat", std::ios::binary);

	wFile.write(reinterpret_cast<const char*>(&cNX), sizeof(int64_t));
	wFile.write(reinterpret_cast<const char*>(&cNZ), sizeof(int64_t));

	wFile.write(reinterpret_cast<const char*>(cWVelField.data()), cWVelField.size() * sizeof(double));

	wFile.close();

	std::ofstream pFile(baseFileName + "P.dat", std::ios::binary);

	pFile.write(reinterpret_cast<const char*>(&cNX), sizeof(int64_t));
	pFile.write(reinterpret_cast<const char*>(&cNZ), sizeof(int64_t));

	pFile.write(reinterpret_cast<const char*>(cPressureField.data()), cPressureField.size() * sizeof(double));

	pFile.close();
}

void FieldRectCylinder2D::ReadField(const std::string& basefilename)
{

	std::ifstream uFile(basefilename + "U.dat", std::ios::binary);

	uFile.read(reinterpret_cast<char*>(&cNX), sizeof(int64_t));
	uFile.read(reinterpret_cast<char*>(&cNZ), sizeof(int64_t));

	cUVelField.clear();
	cUVelField.resize(cNX * cNZ);

	uFile.read(reinterpret_cast<char*>(cUVelField.data()), cUVelField.size() * sizeof(double));

	uFile.close();
	
	std::ifstream wFile(basefilename + "W.dat", std::ios::binary);

	wFile.read(reinterpret_cast<char*>(&cNX), sizeof(int64_t));
	wFile.read(reinterpret_cast<char*>(&cNZ), sizeof(int64_t));

	cWVelField.clear();
	cWVelField.resize(cNX * cNZ);

	wFile.read(reinterpret_cast<char*>(cWVelField.data()), cWVelField.size() * sizeof(double));

	wFile.close();

	std::ifstream pFile(basefilename + "P.dat", std::ios::binary);

	pFile.read(reinterpret_cast<char*>(&cNX), sizeof(int64_t));
	pFile.read(reinterpret_cast<char*>(&cNZ), sizeof(int64_t));

	cPressureField.clear();
	cPressureField.resize(cNX * cNZ);

	pFile.read(reinterpret_cast<char*>(cPressureField.data()), cPressureField.size() * sizeof(double));

	pFile.close();
}

void FieldRectCylinder2D::SaveIKCutFieldToCSVFile(const std::string& baseFileName) const
{
	std::ofstream pFile(baseFileName + "P.csv");

	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			pFile << cPressureField[((K)*cNX) + I] << "\t";
		}
		pFile << "\n";
	}

	pFile.close();

	std::ofstream pcorrFile(baseFileName + "Pcorr.csv");

	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			pcorrFile << cPressureCorrField[(K * cNX) + I] << "\t";
		}
		pcorrFile << "\n";
	}

	pcorrFile.close();

	std::ofstream uFile(baseFileName + "U.csv");

	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			uFile << cUVelField[((K)*cNX) + I] << "\t";
		}
		uFile << "\n";
	}

	uFile.close();
	
	std::ofstream wFile(baseFileName + "W.csv");

	for (int64_t K = 0; K < cNZ; K++)
	{
		for (int64_t I = 0; I < cNX; I++)
		{
			wFile << cWVelField[ ((K)*cNX) + I] << "\t";
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

inline double FieldRectCylinder2D::GetFwInternalUMomentum(int64_t i, int64_t K) const
{
	return (cRHO / 2.0) * (cUVelField[ ((K)*cNX) + i] + cUVelField[ ((K)*cNX) + i - 1]);
}
inline double FieldRectCylinder2D::GetFeInternalUMomentum(int64_t i, int64_t K) const
{
	return (cRHO / 2.0) * (cUVelField[ ((K)*cNX) + i + 1] + cUVelField[ ((K)*cNX) + i]);
}
inline double FieldRectCylinder2D::GetFbInternalUMomentum(int64_t i, int64_t K) const
{
	return (cRHO / 2.0) * (cWVelField[ ((K)*cNX) + i] + cWVelField[ ((K)*cNX) + i - 1]);
}
inline double FieldRectCylinder2D::GetFtInternalUMomentum(int64_t i, int64_t K) const
{
	return (cRHO / 2.0) * (cWVelField[ ((K + 1)*cNX) + i] + cWVelField[ ((K + 1)*cNX) + i - 1]);
}

inline double FieldRectCylinder2D::GetFwInternalWMomentum(int64_t I, int64_t k) const
{
	return (cRHO / 2.0) * (cUVelField[ ((k)*cNX) + I] + cUVelField[ ((k - 1)*cNX) + I]);
}
inline double FieldRectCylinder2D::GetFeInternalWMomentum(int64_t I, int64_t k) const
{
	return (cRHO / 2.0) * (cUVelField[ ((k)*cNX) + I + 1] + cUVelField[ ((k - 1)*cNX) + I + 1]);
}
inline double FieldRectCylinder2D::GetFbInternalWMomentum(int64_t I, int64_t k) const
{
	return (cRHO / 2.0) * (cWVelField[ ((k)*cNX) + I] + cWVelField[ ((k - 1)*cNX) + I]);
}
inline double FieldRectCylinder2D::GetFtInternalWMomentum(int64_t I, int64_t k) const
{
	return (cRHO / 2.0) * (cWVelField[ ((k)*cNX) + I] + cWVelField[ ((k + 1)*cNX) + I]);
}

void FieldRectCylinder2D::SetSpUValueInterval(const Intervalo& II,const Intervalo& kk, double SpValue)
{
	for (int64_t k = kk.valor1; k < kk.valor2; k++)
	{
		for (int64_t i = II.valor1; i < II.valor2; i++)
		{
			cSpUField[((k)*cNX) + i] = SpValue;
		}
	}
}
void FieldRectCylinder2D::SetSuUValueInterval(const Intervalo& II,const Intervalo& kk, double SuValue)
{
	for (int64_t k = kk.valor1; k < kk.valor2; k++)
	{
		for (int64_t i = II.valor1; i < II.valor2; i++)
		{
			//cSuUField[ ((k)*cNX) + i] = SuValue;
		}
	}
}

void FieldRectCylinder2D::SetSpWValueInterval(const Intervalo& II,const Intervalo& kk, double SpValue)
{
	for (int64_t k = kk.valor1; k < kk.valor2; k++)
	{
		for (int64_t i = II.valor1; i < II.valor2; i++)
		{
			cSpWField[((k)*cNX) + i] = SpValue;
		}
	}
}
void FieldRectCylinder2D::SetSuWValueInterval(const Intervalo& II,const Intervalo& kk, double SuValue)
{
	for (int64_t k = kk.valor1; k < kk.valor2; k++)
	{
		for (int64_t i = II.valor1; i < II.valor2; i++)
		{
			//cSuWField[ ((k)*cNX) + i] = SuValue;
		}
	}
}

inline double FieldRectCylinder2D::Psir(double r) const
{
	return (r + (r * r)) / (1.0 + (r * r));
}

inline double FieldRectCylinder2D::PsirSUPERBEE(double r) const
{
	const double i1 = std::fmin(2.0*r, 1.0);
	const double i2 = std::fmin(r, 2.0);
	const double res = std::fmax(0.0, std::fmax(i1, i2));

	return res;
}

double FieldRectCylinder2D::rep(int64_t I, int64_t K, const doubleField2D& field) const
{
	const double phiP = field[ ((K)*cNX) + I];
	const double phiW = field[ ((K)*cNX) + I - 1];
	const double phiE = field[ ((K)*cNX) + I + 1];

	if (std::fabs(phiE - phiP) < 1e-5)
		return 10e30;

	return (phiP - phiW) / (phiE - phiP);
}
double FieldRectCylinder2D::ren(int64_t I, int64_t K, const doubleField2D& field) const
{
	const double phiP = field[ ((K)*cNX) + I];
	const double phiE = field[ ((K)*cNX) + I + 1];
	double phiEE;

	if (I == cNX - 2)
		phiEE = (2.0 * phiE) - phiP;
	else
		phiEE = field[ ((K)*cNX) + I + 2];

	if (std::fabs(phiE - phiP) < 1e-5)
		return 10e30;

	return (phiEE - phiE) / (phiE - phiP);
}

double FieldRectCylinder2D::rwp(int64_t I, int64_t K, bool ufield, const doubleField2D& field) const
{
	const double phiW = field[ ((K)*cNX) + I - 1];
	const double phiP = field[ ((K)*cNX) + I];

	const int iInicial = (ufield) ? 2 : 1;

	double phiWW;
	if (I == iInicial)
		phiWW = (2.0 * phiW) - phiP;
	else
		phiWW = field[ ((K)*cNX) + I - 2];

	if (std::fabs(phiW - phiP) < 1e-5)
		return 10e30;

	return (phiW - phiWW) / (phiP - phiW);
}
double FieldRectCylinder2D::rwn(int64_t I, int64_t K, const doubleField2D& field) const
{
	const double phiE = field[ ((K)*cNX) + I + 1];
	const double phiP = field[ ((K)*cNX) + I];
	const double phiW = field[ ((K)*cNX) + I - 1];

	if (std::fabs(phiW - phiP) < 1e-5)
		return 10e30;

	return (phiE - phiP) / (phiP - phiW);
}

double FieldRectCylinder2D::rbp(int64_t I, int64_t K, bool wfield, const doubleField2D& field) const
{
	const double phiB = field[ ((K - 1)*cNX) + I];
	const double phiP = field[ ((K)*cNX) + I];

	const int kInicial = wfield ? 2 : 1;

	double phiBB;
	if (K == kInicial)
		phiBB = (2.0 * phiB) - phiP;
	else
		phiBB = field[ ((K - 2)*cNX) + I];

	if (std::fabs(phiB - phiP) < 1e-5)
		return 10e30;

	return (phiB - phiBB) / (phiP - phiB);
}
double FieldRectCylinder2D::rbn(int64_t I, int64_t K, const doubleField2D& field) const
{
	const double phiP = field[((K)*cNX) + I];
	const double phiT = field[((K + 1)*cNX) + I];
	const double phiB = field[((K - 1)*cNX) + I];

	if (std::fabs(phiP - phiB) < 1e-5)
		return 10e30;

	return (phiT - phiP) / (phiP - phiB);
}

double FieldRectCylinder2D::rtp(int64_t I, int64_t K, const doubleField2D& field) const
{
	const double phiP = field[ ((K)*cNX) + I];
	const double phiT = field[ ((K + 1)*cNX) + I];
	const double phiB = field[ ((K - 1)*cNX) + I];

	if (std::fabs(phiP - phiT) < 1e-5)
		return 10e30;

	return (phiP - phiB) / (phiT - phiP);
}
double FieldRectCylinder2D::rtn(int64_t I, int64_t K, const doubleField2D& field) const
{
	const double phiT = field[ ((K + 1)*cNX) + I];
	const double phiP = field[ ((K)*cNX) + I];

	double phiTT;
	if (K == cNZ - 2)
		phiTT = (2.0 * phiT) - phiP;
	else
		phiTT = field[ ((K + 2)*cNX) + I];

	if (std::fabs(phiP - phiT) < 1e-5)
		return 10e30;

	return (phiTT - phiT) / (phiT - phiP);
}

double FieldRectCylinder2D::PsirCD(double r) const
{
	return 1.0;
}