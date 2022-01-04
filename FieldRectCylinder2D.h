#pragma once
#include <vector>
#include <string>

#include "Timer.h"

typedef std::vector<double> doubleField2D;

struct CoefsData2D
{
public:
	double aijk;
	double aipjk;
	double aimjk;
	double aijkp;
	double aijkm;

	double source;
};

class FieldRectCylinder2D
{
public:
	class Intervalo
	{
	public:
		Intervalo(int64_t v1, int64_t v2);
		int64_t valor1;
		int64_t valor2;

		int64_t Tamanho() const;
	};
public:
	// W = w = H
	// h = l
	FieldRectCylinder2D(double L, double H, double l, double deltaSize, double rho, double mu,
		double Ufarfield, double Wfarfield);
	~FieldRectCylinder2D();

	inline int64_t GetNX() const { return cNX; }
	inline int64_t GetNZ() const { return cNZ; }

	// RECTANGLE CYLINDER 
	void SetBCFlatPlate();

	int CreateUMomentumLSCSR(
		const doubleField2D& uvel0,
		double dt,
		std::vector<int>& ptr,
		std::vector<int>& col,
		std::vector<double>& val,
		std::vector<double>& rhs
	);
	int CreateWMomentumLSCSRParallel(
		const doubleField2D& Wvel0,
		double dt,
		std::vector<int>& ptr,
		std::vector<int>& col,
		std::vector<double>& val,
		std::vector<double>& rhs
	);
	int CreatePressureCorrectionLSCSR(
		std::vector<int>& ptr,
		std::vector<int>& col,
		std::vector<double>& val,
		std::vector<double>& rhs
	) const;

	int CreatePressureLSCSRFI(
		std::vector<int>& ptr,
		std::vector<int>& col,
		std::vector<double>& val,
		std::vector<double>& rhs
	) const;


	void CalcUHatValuesFI(const doubleField2D& uField0, double dt);
	void CalcWHatValuesFI(const doubleField2D& wField0, double dt);

	void CalcUHatValuesCN(const doubleField2D& uField0, double dt);
	void CalcWHatValuesCN(const doubleField2D& wField0, double dt);

	// Set P value for all 'K' and 'I' in J defined; 
	void SetPValueIK(double pValue);

	// Set U value for all 'K' and 'i' in J defined; 
	void SetUValueiK(double uValue);

	// Set W value for all 'k' and 'I' in J defined; 
	void SetWValueIk(double wValue);

	void SetPValueInterval(const Intervalo& II, const Intervalo& KK, double pValue);
	void SetUValueInterval(const Intervalo& ii, const Intervalo& KK, double uValue);
	void SetWValueInterval(const Intervalo& II, const Intervalo& kk, double wValue);

	void SetSpUValueInterval(const Intervalo& II, const Intervalo& kk, double SpValue);
	void SetSuUValueInterval(const Intervalo& II, const Intervalo& kk, double SuValue);
	
	void SetSpWValueInterval(const Intervalo& II, const Intervalo& kk, double SpValue);
	void SetSuWValueInterval(const Intervalo& II, const Intervalo& kk, double SuValue);

	void ExtrapolateForwardPJK(int64_t I);
	void ExtrapolateForwardUJK(int64_t i);
	void ExtrapolateForwardWJk(int64_t I);

	void ExtrapolateBackwardPJK(int64_t I);
	void ExtrapolateBackwardUJK(int64_t i);
	void ExtrapolateBackwardWJk(int64_t I);

	void ExtrapolateBackwardWeightedUJK(int64_t i, double weight);
	void ExtrapolateBackwardWeightedWJk(int64_t I, double weight);

	void ExtrapolateForwardIntervalPJK(int64_t I, const Intervalo& KK);
	void ExtrapolateForwardIntervalUJK(int64_t i, const Intervalo& KK);
	void ExtrapolateForwardIntervalWJk(int64_t I, const Intervalo& kk);

	void ExtrapolateBackwardIntervalPJK(int64_t I, const Intervalo& KK);
	void ExtrapolateBackwardIntervalUJK(int64_t i, const Intervalo& KK);
	void ExtrapolateBackwardIntervalWJk(int64_t I, const Intervalo& kk);
	
	void SaveIKCutFieldToCSVFile(const std::string& baseFileName) const;

	void SaveField(const std::string& baseFileName) const;

	void ReadField(const std::string& basefilename);

	double GetFwInternalUMomentum(int64_t i, int64_t K) const;
	double GetFeInternalUMomentum(int64_t i, int64_t K) const;
	double GetFbInternalUMomentum(int64_t i, int64_t K) const;
	double GetFtInternalUMomentum(int64_t i, int64_t K) const;
	
	double GetFwInternalWMomentum(int64_t I, int64_t k) const;
	double GetFeInternalWMomentum(int64_t I, int64_t k) const;
	double GetFbInternalWMomentum(int64_t I, int64_t k) const;
	double GetFtInternalWMomentum(int64_t I, int64_t k) const;

private:
	double Psir(double r) const;
	double PsirSUPERBEE(double r) const;

	double PsirCD(double r) const;

	double rep(int64_t I, int64_t K, const doubleField2D& field) const;
	double ren(int64_t I, int64_t K, const doubleField2D& field) const;

	double rwp(int64_t I, int64_t K, bool ufield, const doubleField2D& field) const;
	double rwn(int64_t I, int64_t K, const doubleField2D& field) const;

	double rbp(int64_t I, int64_t K, bool wfield, const doubleField2D& field) const;
	double rbn(int64_t I, int64_t K, const doubleField2D& field) const;

	double rtp(int64_t I, int64_t K, const doubleField2D& field) const;
	double rtn(int64_t I, int64_t K, const doubleField2D& field) const;

public:
	Timer cTimer;

	const int ctidPushback;
	const int ctidCalcCoefs;
	const double cdd;
	const double cdA;
	const double cdV;

	int64_t cNX;
	int64_t cNZ;

	const int64_t cNXCylinder;
	const int64_t cNZCylinder;

	const int64_t ciCyInit;
	const int64_t ckCyInit;

	const int64_t ciCyEnd;
	const int64_t ckCyEnd;

	const double cRHO;
	const double cMU;
	const double cD;

	doubleField2D cPressureField;
	doubleField2D cPressureCorrField;

	doubleField2D cUVelField;
	doubleField2D cWVelField;

	doubleField2D cSpUField;
	doubleField2D cSpWField;

	doubleField2D caijkU;
	doubleField2D caijkW;

	doubleField2D cUHatVelField;
	doubleField2D cWHatVelField;
};