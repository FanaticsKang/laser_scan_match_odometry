#pragma once

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>

namespace Alpha {

class TicToc {
 public:
  TicToc() { Tic(); }

  void Tic() { start = std::chrono::system_clock::now(); }

  // return ms
  double TocMicroseconds() const {
    std::chrono::time_point<std::chrono::system_clock> end =
        std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    return elapsed_seconds.count() * 1000;
  }

  double TocSecond() const {
    std::chrono::time_point<std::chrono::system_clock> end =
        std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    return elapsed_seconds.count();
  }

  friend std::ostream& operator<<(std::ostream& os, const TicToc& tic) {
    return os << "Time Costs: " << tic.TocMicroseconds() << " ms";
  }

 private:
  std::chrono::time_point<std::chrono::system_clock> start;
};
}