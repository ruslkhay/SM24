#pragma once
#include "solution.hpp"
#include <array>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
// #include <mpi.h>

void debugReceivePrint(const std::vector<double> &storage, const int currRank,
                       const int prevRank, std::pair<int, int> buffSize);

void debugSendPrint(const std::vector<std::vector<double>> &grid,
                    const int currRank, const int nextRank, int x0, int xM,
                    int y0, int yN, const std::vector<double> &boarderVal);

std::array<int, 4> GetLimitsTwoProc(int rank, int M, int N);

std::array<int, 4> GetLimitsFourProc(int rank, int M, int N);

std::array<int, 4> GetSectors(int procNum, int rank, int M, int N);

std::vector<std::vector<double>>
prepareSubGrid(const std::vector<std::vector<double>> &grid, int x0, int xM,
               int y0, int yN);