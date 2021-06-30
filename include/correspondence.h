#pragma once

struct Correspondence {
  /** 1 if this correspondence is valid  */
  int valid;
  /** Closest point in the other scan.  */
  int j1;
  /** Second closest point in the other scan.  */
  int j2;
  /** Type of correspondence (point to point, or point to line) */
  enum { corr_pp = 0, corr_pl = 1 } type;
  /** Squared distance from p(i) to point j1 */
  double dist2_j1;
};