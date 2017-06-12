/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "PhredHelper.h"

#include <math.h>

phredConverter phredConv;

phredConverter::phredConverter()
{
    // Create a quick lookup table to speed up conversion of
    // base quality values stored as log10 (error rates) into
    // fractional error rates
  for (int i = 0; i <= 255; i++) {
    phred2Err[i] = (i > 1) ? pow(0.1, i * 0.1) : 0.75;
    phred2Prob[i] = pow(0.1, i * 0.1);
    phred2Mat[i] = 1.-phred2Err[i];
    phred2Mat3[i] = 1.-phred2Err[i]/3.;
    phred2LogMat[i] = log10(phred2Mat[i]);
    phred2LogMat3[i] = log10(phred2Mat3[i]);
    phred2HalfLogMat3[i] = log10(0.5-phred2Err[i]/3);
  }
  log3 = log10(3.0);
    // doubleLookup[255] = 0.0;
}

uint8_t phredConverter::err2Phred(double err) {
  double d = 0-10*log10(err);
  return ( d > 255 ? 255 : ( d < 0 ? 0 : (uint8_t)floor(d+.5) ) );
}

uint8_t phredConverter::mat2Phred(double mat) {
  double d = 0-10*log10(1.-mat);
  return ( d > 255 ? 255 : ( d < 0 ? 0 : (uint8_t)floor(d+.5) ) );
}

uint8_t phredConverter::mat32Phred(double mat3) {
  double d = 0-10*log10(1.-mat3);
  return ( d > 255 ? 255 : ( d < 0 ? 0 : (uint8_t)floor(d+.5) ) );
}

