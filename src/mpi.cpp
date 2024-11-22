#include "mpi.hpp"

void debugReceivePrint(const std::vector<double> &storage, const int currRank,
                       const int prevRank, std::pair<int, int> buffSize) {
  auto message = std::string("<< Process ");
  message = message + std::to_string(currRank) + " received from " +
            std::to_string(prevRank) + "\n";
  for (int i = 0; i < buffSize.first * buffSize.second; i++) {
    if (i == 0) {
      message = message + " ";
    } else if (i % buffSize.first == 0) {
      message = message + "\n ";
    }
    message = message + std::to_string(storage[i]) + "| ";
  }
  message = message + "\n==========\n";
  std::cout << message;
}

void debugSendPrint(const std::vector<std::vector<double>> &grid,
                    const int currRank, const int nextRank, int x0, int xM,
                    int y0, int yN, const std::vector<double> &boarder) {
  auto message =
      std::string(">> Process " + std::to_string(currRank) + " doing:\n");
  switch (currRank) {
  case 0:
    message += "  Top Left Sub-grid:\n";
    break;
  case 1:
    message += "  Bottom Left Sub-grid:\n";
    break;
  case 2:
    message += "  Bottom Right Sub-grid:\n";
    break;
  case 3:
    message += "  Top Right Sub-grid:\n";
    break;
  }
  for (int i = 0; i < static_cast<int>(grid[0].size()); ++i) {
    for (auto col : grid) {
      message += std::to_string(col[i]);
      message += "|";
    }
    message += '\n';
  }
  message += "  And sending [";
  for (auto elem : boarder) {
    message += std::to_string(elem) + " ";
  }
  message += "] to process " + std::to_string(nextRank) + "\n";
  message += "----------\n";
  std::cout << message;
}

std::array<int, 4> GetLimitsTwoProc(int rank, int M, int N) {
  // + 1 appears because I've chosen specific splitting for 2 proc case
  const int xMiddle = (M + 1) / 2 + 1;
  int x0 = 0, xM = 0;
  switch (rank) {
  case 0:
    x0 = 0, xM = xMiddle;
    break;
  case 1:
    x0 = M - xMiddle, xM = M;
    break;
  }
  std::array<int, 4> limits{x0, xM, 0, N};
  return limits;
}

std::array<int, 4> GetLimitsFourProc(int rank, int M, int N) {
  const int xMiddle = (M + 1) / 2, yMiddle = (N + 1) / 2;
  // TODO: Check if it works correct (old version below)
  // const int xMiddle = M / 2 + 1, yMiddle = N / 2 + 1;
  int x0 = 0, xM = 0, y0 = 0, yN = 0;
  switch (rank) {
  case 0:
    x0 = 0, xM = xMiddle, y0 = 0, yN = N - yMiddle;
    break;
  case 1:
    x0 = 0, xM = xMiddle, y0 = N - yMiddle, yN = N;
    break;
  case 2:
    x0 = xMiddle, xM = M, y0 = N - yMiddle, yN = N;
    break;
  case 3:
    x0 = xMiddle, xM = M, y0 = 0, yN = N - yMiddle;
    break;
  }
  std::array<int, 4> limits{x0, xM, y0, yN};
  return limits;
}

std::array<int, 4> GetSectors(int procNum, int rank, int M, int N) {
  switch (procNum) {
  case 2:
    return GetLimitsTwoProc(rank, M, N);
  case 4:
    return GetLimitsFourProc(rank, M, N);
  default:
    throw std::invalid_argument("Only 1, 2 or 4 processes could be launched\n");
  }
}
std::vector<std::vector<double>>
prepareSubGrid(const std::vector<std::vector<double>> &grid, int x0, int xM,
               int y0, int yN) {
  std::vector<std::vector<double>> tmpBuf(xM - x0,
                                          std::vector<double>(yN - y0, 0));
  // printf("(%d; %d), (%d, %d)\n", x0, xM, y0, yN);
  // printf("tmpgrid size: (%lu; %lu)\n", tmpBuf.size(), tmpBuf[0].size());
  for (int i = x0; i < xM; ++i) {
    for (int j = y0; j < yN; ++j) {
      tmpBuf[i - x0][j - y0] = 1.0 * grid[i][j];
    }
  }
  return tmpBuf;
}
