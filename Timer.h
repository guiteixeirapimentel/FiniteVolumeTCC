#pragma once
#include <Windows.h>
#include <string>
#include <vector>
#include <cmath>

#define _DEBUGTIMER

class Timer
{
private:
	struct TickTockPair
	{
		LARGE_INTEGER tick;
		LARGE_INTEGER tock;
	};
public:
	Timer();
	~Timer();

	// returns the id of the timer
	int SetTimer(const std::string& label);

	void Tick(int id);
	void Tock(int id);

	// if id does not exist -> returns nan
	double GetLastTickTock(int id) const;
	// if id does not exist -> returns  "NA"
	std::string GetLabel(int id) const;

	void WriteToCoutTickTock(int id) const;

	void WriteToCoutTickTockMicrosec(int id) const;

	void WriteToCoutAllTickTock() const;

private:
	LARGE_INTEGER cCountsPerSec;

	std::vector<TickTockPair> cTimers;
	std::vector<std::string> cLabels;

	int cLastIdGiven;
};