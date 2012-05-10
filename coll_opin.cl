__kernel void coll_opin(const int size, const __global float * peps, const __global float * pp, const __global float * plam,
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

    
    // finally do the (output data) assignment:
    // - - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -- - - - - - - -
      
    plcdot[ii] = 1;
    paclam[ii] = 2;

    pacp[ii] = 0.3;
    ppcdot[ii] = 0.4;
  } // end of if i<size

} // end of coll_op method
