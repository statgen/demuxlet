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

#ifndef __PHRED_HELPER_H__
#define __PHRED_HELPER_H__

#include <stdint.h>

class phredConverter
{
public:
  phredConverter();
  double phred2Err[256];
  double phred2Prob[256];
  double phred2Mat[256];
  double phred2Mat3[256];
  double phred2LogMat[256];
  double phred2LogMat3[256];
  double phred2HalfLogMat3[256];
  double log3;

  uint8_t err2Phred(double err);
  uint8_t mat2Phred(double mat);
  uint8_t mat32Phred(double mat3);
  
  inline double toProb(uint32_t phred) { return phred > 255 ? phred2Prob[255] : ( phred < 0 ? 1.0 : phred2Prob[phred]); }
  inline double toErr(uint32_t phred) { return phred > 255 ? phred2Err[255] : (phred < 0 ? 1.0 : phred2Err[phred]); }
  inline double toMat(uint32_t phred) { return phred > 255 ? phred2Mat[255] : (phred < 0 ? 0.0 : phred2Mat[phred]); }
  inline double toMat3(uint32_t phred) { return phred > 255 ? phred2Mat3[255] : (phred < 0 ? 0.0 : phred2Mat3[phred]); }
  inline double toLogMat(uint32_t phred) { return phred > 255 ? phred2LogMat[255] : (phred < 0 ? phred2LogMat[0] : phred2LogMat[phred]); }
  inline double toLogMat3(uint32_t phred) { return phred > 255 ? phred2LogMat3[255] : (phred < 0 ? phred2LogMat3[0] : phred2LogMat3[phred]); }
  inline double toHalfLogMat3(uint32_t phred) { return phred > 255 ? phred2HalfLogMat3[255] : (phred < 0 ? phred2HalfLogMat3[0] : phred2HalfLogMat3[phred]); }
};

extern phredConverter phredConv;

#endif


