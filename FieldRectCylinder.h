#pragma once
#include <vector>
#include <string>

#include "Timer.h"

typedef std::vector<double> doubleField3D;

struct CoefsData
{
public:
	double aijk;
	double aipjk;
	double aimjk;
	double aijpk;
	double aijmk;
	double aijkp;
	double aijkm;

	double source;
};

class FieldRectCylinder
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
	FieldRectCylinder(double L, double H, double W, double l, double deltaSize, double rho, double mu,
		double Ufarfield, double Vfarfield, double Wfarfield);
	~FieldRectCylinder();

	inline int64_t GetNX() const { return cNX; }
	inline int64_t GetNY() const { return cNY; }
	inline int64_t GetNZ() const { return cNZ; }
	
	// RECTANGLE CYLINDER 
	void SetBCFlatPlate();
	
	int CreateUMomentumLSCSR(
		const doubleField3D& uvel0, 
		double dt, 
		std::vector<int>& ptr, 
		std::vector<int>& col, 
		std::vector<double>& val, 
		std::vector<double>& rhs
	);
	int CreateVMomentumLSCSRParallel(
		const doubleField3D& vvel0,
		double dt,
		std::vector<int>& ptr,
		std::vector<int>& col,
		std::vector<double>& val,
		std::vector<double>& rhs
	);
	int CreateWMomentumLSCSRParallel(
		const doubleField3D& Wvel0,
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


	void CalcUHatValuesFI(const doubleField3D& uField0, double dt);
	void CalcVHatValuesFI(const doubleField3D& vField0, double dt);
	void CalcWHatValuesFI(const doubleField3D& wField0, double dt);

	void CalcUHatValuesCN(const doubleField3D& uField0, double dt);
	void CalcVHatValuesCN(const doubleField3D& vField0, double dt);
	void CalcWHatValuesCN(const doubleField3D& wField0, double dt);
	
	// Set P value for all 'I' and 'J' in K defined;
	void SetPValueIJ(int64_t K, double pValue);

	// Set P value for all 'J' and 'K' in I defined;
	void SetPValueJK(int64_t I, double pValue);

	// Set P value for all 'K' and 'I' in J defined; 
	void SetPValueIK(int64_t J, double pValue);

	// Set U value for all 'i' and 'J' in K defined;
	void SetUValueiJ(int64_t K, double uValue);

	// Set U value for all 'J' and 'K' in i defined;
	void SetUValueJK(int64_t i, double uValue);

	// Set U value for all 'K' and 'i' in J defined; 
	void SetUValueiK(int64_t J, double uValue);

	// Set V value for all 'I' and 'K' in j defined;
	void SetVValueIK(int64_t j, double vValue);

	// Set V value for all 'j' and 'I' in K defined;
	void SetVValueIj(int64_t K, double vValue);

	// Set V value for all 'K' and 'j' in I defined; 
	void SetVValuejK(int64_t I, double vValue);

	// Set W value for all 'I' and 'J' in k defined;
	void SetWValueIJ(int64_t k, double wValue);

	// Set W value for all 'J' and 'k' in I defined;
	void SetWValueJk(int64_t I, double wValue);

	// Set W value for all 'k' and 'I' in J defined; 
	void SetWValueIk(int64_t J, double wValue);

	void SetPValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& KK, double pValue);
	void SetUValueInterval(const Intervalo& ii, const Intervalo& JJ, const Intervalo& KK, double uValue);
	void SetVValueInterval(const Intervalo& II, const Intervalo& jj, const Intervalo& KK, double vValue);
	void SetWValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double wValue);

	void SetSpUValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SpValue);
	void SetSuUValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SuValue);

	void SetSpVValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SpValue);
	void SetSuVValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SuValue);

	void SetSpWValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SpValue);
	void SetSuWValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SuValue);

	void ExtrapolateForwardPJK(int64_t I);
	void ExtrapolateForwardUJK(int64_t i);
	void ExtrapolateForwardVjK(int64_t I);
	void ExtrapolateForwardWJk(int64_t I);

	void ExtrapolateBackwardPJK(int64_t I);
	void ExtrapolateBackwardUJK(int64_t i);
	void ExtrapolateBackwardVjK(int64_t I);
	void ExtrapolateBackwardWJk(int64_t I);

	void ExtrapolateBackwardWeightedUJK(int64_t i, double weight);
	void ExtrapolateBackwardWeightedVjK(int64_t I, double weight);
	void ExtrapolateBackwardWeightedWJk(int64_t I, double weight);

	void ExtrapolateForwardIntervalPJK(int64_t I, const Intervalo& JJ, const Intervalo& KK);
	void ExtrapolateForwardIntervalUJK(int64_t i, const Intervalo& JJ, const Intervalo& KK);
	void ExtrapolateForwardIntervalVjK(int64_t I, const Intervalo& jj, const Intervalo& KK);
	void ExtrapolateForwardIntervalWJk(int64_t I, const Intervalo& JJ, const Intervalo& kk);

	void ExtrapolateBackwardIntervalPJK(int64_t I, const Intervalo& JJ, const Intervalo& KK);
	void ExtrapolateBackwardIntervalUJK(int64_t i, const Intervalo& JJ, const Intervalo& KK);
	void ExtrapolateBackwardIntervalVjK(int64_t I, const Intervalo& jj, const Intervalo& KK);
	void ExtrapolateBackwardIntervalWJk(int64_t I, const Intervalo& JJ, const Intervalo& kk);

	void ExtrapolateForwardIntervalPJI(int64_t K, const Intervalo& II, const Intervalo& JJ);
	void ExtrapolateForwardIntervalUJI(int64_t K, const Intervalo& ii, const Intervalo& JJ);
	void ExtrapolateForwardIntervalVJI(int64_t K, const Intervalo& II, const Intervalo& jj);
	void ExtrapolateForwardIntervalWJI(int64_t k, const Intervalo& II, const Intervalo& JJ);

	void ExtrapolateBackwardIntervalPJI(int64_t K, const Intervalo& II, const Intervalo& JJ);
	void ExtrapolateBackwardIntervalUJI(int64_t K, const Intervalo& ii, const Intervalo& JJ);
	void ExtrapolateBackwardIntervalVJI(int64_t K, const Intervalo& II, const Intervalo& jj);
	void ExtrapolateBackwardIntervalWJI(int64_t k, const Intervalo& II, const Intervalo& JJ);

	void SaveIJCutFieldToCSVFile(int64_t K, const std::string& baseFileName) const;
	void SaveIKCutFieldToCSVFile(int64_t J, const std::string& baseFileName) const;
	void SaveJKCutFieldToCSVFile(int64_t I, const std::string& baseFileName) const;

	void SaveField(const std::string& baseFileName) const;

	void ReadField(const std::string& basefilename);

	double GetFwInternalUMomentum(int64_t i, int64_t J, int64_t K) const;
	double GetFeInternalUMomentum(int64_t i, int64_t J, int64_t K) const;
	double GetFsInternalUMomentum(int64_t i, int64_t J, int64_t K) const;
	double GetFnInternalUMomentum(int64_t i, int64_t J, int64_t K) const;
	double GetFbInternalUMomentum(int64_t i, int64_t J, int64_t K) const;
	double GetFtInternalUMomentum(int64_t i, int64_t J, int64_t K) const;

	double GetFwInternalVMomentum(int64_t I, int64_t j, int64_t K) const;
	double GetFeInternalVMomentum(int64_t I, int64_t j, int64_t K) const;
	double GetFsInternalVMomentum(int64_t I, int64_t j, int64_t K) const;
	double GetFnInternalVMomentum(int64_t I, int64_t j, int64_t K) const;
	double GetFbInternalVMomentum(int64_t I, int64_t j, int64_t K) const;
	double GetFtInternalVMomentum(int64_t I, int64_t j, int64_t K) const;

	double GetFwInternalWMomentum(int64_t I, int64_t J, int64_t k) const;
	double GetFeInternalWMomentum(int64_t I, int64_t J, int64_t k) const;
	double GetFsInternalWMomentum(int64_t I, int64_t J, int64_t k) const;
	double GetFnInternalWMomentum(int64_t I, int64_t J, int64_t k) const;
	double GetFbInternalWMomentum(int64_t I, int64_t J, int64_t k) const;
	double GetFtInternalWMomentum(int64_t I, int64_t J, int64_t k) const;

private:
	double Psir(double r) const;
	double PsirSUPERBEE(double r) const;

	double PsirCD(double r) const;

	double rep(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const;
	double ren(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const;
					 					  	
	double rwp(int64_t I, int64_t J, int64_t K, bool ufield, const doubleField3D& field) const;
	double rwn(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const;
					 					  		  
	double rnp(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const;
	double rnn(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const;
					 					  		 
	double rsp(int64_t I, int64_t J, int64_t K, bool vfield, const doubleField3D& field) const;
	double rsn(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const;
					 					  		 
	double rbp(int64_t I, int64_t J, int64_t K, bool wfield, const doubleField3D& field) const;
	double rbn(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const;
					 					 		
	double rtp(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const;
	double rtn(int64_t I, int64_t J, int64_t K, const doubleField3D& field) const;

public:
	Timer cTimer;

	const int ctidPushback;
	const int ctidCalcCoefs;
	const double cdd;
	const double cdA;
	const double cdV;

	int64_t cNX;
	int64_t cNY;
	int64_t cNZ;
	 
	const int64_t cNXCylinder;
	const int64_t cNYCylinder;
	const int64_t cNZCylinder;
	 
	const int64_t ciCyInit;
	const int64_t cjCyInit;
	const int64_t ckCyInit;
	 
	const int64_t ciCyEnd;
	const int64_t cjCyEnd;
	const int64_t ckCyEnd;
	 
	const double cRHO;
	const double cMU;
	const double cD;

	doubleField3D cPressureField;
	doubleField3D cPressureCorrField;

	doubleField3D cUVelField;
	doubleField3D cVVelField;
	doubleField3D cWVelField;

	doubleField3D cSpUField;	
	doubleField3D cSpVField;	
	doubleField3D cSpWField;

	doubleField3D caijkU;
	doubleField3D caijkV;
	doubleField3D caijkW;

	doubleField3D cUHatVelField;
	doubleField3D cVHatVelField;
	doubleField3D cWHatVelField;
};