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
** ADDITIONAL MODIFICATIONS, JAKOB RYDEN FOR ARENA ON OPEN-CL:
** The collision operator below was translated to C by Jakob Ryden from the original Fortran
** Program ARENA by Lars-Goran Eriksson, the original source can be obtained with gforge credentials
** which are attainable from the ITM task force www.efda-itm.eu
**
** the first three inputs are the ARENA phase space variables;
** eps = r/R: inverse aspect ratio
** p: normalised momentum, will grow to large when a runaway
** lam: the banana orbit indicator variable
** 
** the next two inputs are the pre-calculated integral values, which are common for all MC particles
**
** finally, the last four variables are the monte carlo outputs, p_dot and lam_dot and their respective variances:
** ppcdot: momentum time derivative
** pacp: momentum variance
** plcdot: lambda time derivative
** paclam: lambda value variance
** 
** 
** In the end are required non-global variables that should be sent or accessed separately in a perfect world,
** and also because it is nice to have a similar function call as in the current ARENA code.
** 
** NOTE: the peps input is not used for more than a temperature profile - here the temperature is assumed constant
** across the (small) simulation domain.
** 
** NOTE2: the input type here defines the variable, if input type is a pointer (*)
** it will point to the beginning of an array
** 
** 
** $Id$
**/
// in the input, the 'const' variables are INPUT ONLY

__kernel void coll_op(const int size, const __global float * peps, const __global float * pp, const __global float * plam,
    __global const float * I1, __global const float * I2,
    __global float * ppcdot, __global float * pacp, __global float * plcdot, __global float * paclam,
    const float temperature, const float Z_eff, const float coulomb_logarithm){
    
    // hard code all globals for now, this will have to be fixed!!


    // the OpenCL gloabal ID:
    // - - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -
    int ii = get_global_id(0);
  
    if(ii < size) {
    // collision operator here:
    // section of copied physical parameters to not have to send them across functions at this point:
    // - - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -

    //!  float pi =      3.14159265358979323846264338327950288419716939937510
    float pi =      3.1416;
    float pi_sqrt = 1.7724538509055160;
    float speed_of_light = 3.0e8;
    float c0 = 2.99792458e8;
    float k_boltz = 1.3806504e-23;
    float r_e = 2.8179e-15;
    float electron_mass_in_ev = 511.17;
    float electron_mass = 9.10938188e-31;
    float electron_charge = 1.60217646e-19;
    float phi_c1 =  0.254829592; // THIS IS ANTS OPERATOR INTERPOLATION STUFF <- what is this??
    float phi_c2 = -0.284496736;
    float phi_c3 =  1.421413741;
    float phi_c4 = -1.453152027;
    float phi_c5 =  1.061405429;
    float eps0 = 8.85418782e-12;

    
    // SPLINE OUTPUTS:
    float ppcdot1;
    float plcdot1;
    float pacp1;
    float paclam1;
    
    float ppcdot2;
    float plcdot2;
    float pacp2;
    float paclam2;
    
// - - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -

    // Main copied collision operator part: This could be updated to look nicer!
    // - - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -


    // this is really a VELOCITY(or normalised momentum perhaps)
    float pth = sqrt(temperature/electron_mass_in_ev);
    float pth2 = temperature / electron_mass_in_ev; // the square of pth
    float pth_real = pth * sqrt (2.0);// ! zzpth: more modern definition of therm. vel., used in the lowE coll.op
    
    float pi1 = I1[ii];
    float pi2 = I2[ii];
    
    float p = pp[ii] + 1.0E-6;
    float p2 = p * p;
    float p3 = p2 * p;
    float p4 = p3 * p;
    
    float pc = 0.4 * pth;
    float pc4 = pc * pc * pc * pc;
    float pcpa = 2.0 * pth;


    float pb;
    float pb2;
    float pb3;
    float pb4;
    
    if (p < pcpa){
        pb = pcpa;
        pb2 = pcpa * pcpa;
        pb3 = pb2 * pcpa;
        pb4 = pb3 * pcpa;
    }
    else {
        pb2 = p * p;
        pb3 = p2 * p;
        pb4 = p3 * p;
    }
    
    // ztheta = 1.9563E-3 * temperature // actually exactly the same as pth, not the 'more modern' one
    // REALLY NO NEED TO COMPUTE AGAIN
    float beta = 1.9563E-3 * temperature; // ztheta = T_e / (m_e c^2), = pth!! (verifyable)
    
    // temporary variables that return in many equations...
    float zt1 = (1.0+p2);
    float zt2 = (1.0+logf(p+1.0)/coulomb_logarithm);
    float zt3 = pc4 + p4;
    float zt4 = sqrtf(zt1);
    
    float zz = beta * (1.0+2.0*pb2) / (pb2*(1.0+pb2));
    float limit_no = 1.0;
    zz = fmin(zz, limit_no); // min() is minmag in C
    
    // old zbeta
    float nu_D = 2.0 * (1.0+ Z_eff - zz) * sqrt (1.0+pb2) / pb3 * (1.0+log(p+1.0)/coulomb_logarithm);

    float nu = zt1 * zt2 / p3;
    nu = nu * p4 / zt3;
    float dnu = (1.0+3.0*p2-4.0*(1.0+p2)*p4/zt3) / zt3;
    // High ENERGY COLL OP.
    // - - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -
    ppcdot1 = -p * nu + beta*(((2.0*nu/p+dnu)*zt4+nu*p/zt4)*zt2+nu*zt4/((p+1.0)*coulomb_logarithm));
    pacp1 = sqrt(2.0* nu * beta*zt4);

    
    //! MJ  LOW ENERGY COLL OP.
    // - - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -

    float x = p / pth_real; // the normalised momentum!
    float x2 = x * x;
    
    //THE CHANDRASEKHAR FUNCTION HAS TO BE IMPLEMENTED INLINE HERE, TRY TO SEPARATE!
    // - - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -
    // chandra = gstix (zx), old: zg

    float chandra;
    float zt;
    float phi;
    if (x < 0.01) {
        chandra = 1.1284 * (0.33333 * x + 0.2 * pow(x,3) );
    }   
    else
    {
      if( x2 < 20.0)
      {
        // THE WHOLE NESTED CHANDRASEKHAR FUNCTION, this was 'gstix()'
        zt = 1/(1+0.3275911*x); // input parameter, WHAT IS THE CONSTANTS??
        phi = ((((phi_c5*zt + phi_c4)*zt + phi_c3)*zt + phi_c2)*zt + phi_c1)*zt * exp(x2);
        chandra = phi; // chandra = psi(phi(... ))
      }
      else
        chandra = 1.0/2.0/x2;
        
    }
    // MORE, the thermal temperature?? the scaling...
    // - - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -

    float pth_real2 = pow(pth_real,2);

    float ve = speed_of_light * pth_real / sqrt(1.0 + pth_real2);
    float nu_ve = speed_of_light * speed_of_light * speed_of_light / (ve*ve*ve);
    float exp_x2 = exp(-x2);
    
    // OUTPUT:
    pacp2 = sqrt(2.0*nu_ve * pth_real2 * chandra / x);
    ppcdot2 = nu_ve * pth_real * (-chandra *(2.0+1.0/x2) + 1.1248 * exp_x2 / x);
    
    // MJ  THE DERRIVATIVE PART OF L
    float alpha = p / sqrt(2*pth2*(1+p2));  // note, the INPUT PARAMETER p! - also "not real" pth 
    // this is used as arg for chandrasekhar function
    // so is it some kind of momentum?
    
    // reuse the same variables, zt and phi:
    zt = 1/(1 + 0.3275911*alpha); // input parameter, WHAT IS THE CONSTANTS STILL??
    chandra = ((((phi_c5*zt + phi_c4)*zt + phi_c3)*zt + phi_c2)*zt + phi_c1)*zt * exp(pow(alpha,2));
    
    float b1 = pth2 * zt1 * zt4 / p3; // note: uses "OLD" p_thermal^2
    float b2 = pth2 * zt1 * zt4 * chandra / p;
    
    // LOGISTICS function evaluation to speed things up, this was a separate function
    // why the hell would that speed things up??
    // - - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -

    // logisv = logistic (zp, temperature) 
    // logder = logisv * (1-logisv) // !derrivative of logistics function
    float xarg = (p - pth) / (0.1*pth);
    float logisv = 1/(1 + exp(-xarg) );
    float logder = logisv * (1-logisv);
    
    // finally do the (output data) assignment:
    // - - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -
      
    plcdot[ii] = -nu_D * (0.5 * plam[ii] - pi2 / pi1);
    paclam[ii] = sqrt(2.0*nu_D * plam[ii] * pi2/pi1);

    pacp[ii] = pacp1 * logisv + pacp2 * (1-logisv);
    ppcdot[ii] = ppcdot1 * logisv + ppcdot2 * (1-logisv) + b1 * logder - b2 * logder;
  } // end of if i<size

} // end of coll_op method