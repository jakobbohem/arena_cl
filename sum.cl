/* Copyright (C) 2011 X. Andrade
**
** FortranCL is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** FortranCL is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
** $Id$
**/

//__kernel void coll_op(const int size, const __global float * peps, const __global float * pp, const __global float * plam,
//    const __global float * pi1, const __global float * pi2,
//    __global float * ppcdot, __global float * pacp, __global float * plcdot, __global float * paclam,
//    const __global float * temperature, const __global float * Z_eff, const __global float * coulomb_logarithm,
//    ){
__kernel void sum(const int size, const __global float * vec1, __global float * vec2, float multiplier, const float lightspeed) {
  int ii = get_global_id(0);

  if(ii < size) {
    //vec2[ii] += vec1[ii];
    float output;
    
    // ("return the square root of 1+2 to test library inclusion");
    if (lightspeed > 3e8)
        output = sqrtf(vec1[ii]+vec2[ii]) * multiplier;
    else
    {
        output = sqrtf(vec1[ii]+vec2[ii]) * lightspeed;
    }
    
    
    vec2[ii] = output;
    } // end if

} // end of 'sum' method
