## Process this file with automake to produce Makefile.in

## Copyright (C) 2011 X. Andrade <xavier@tddft.org>
## 
## FortranCL is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## FortranCL is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
## $Id$

noinst_PROGRAMS = sum devices multiply

sum_SOURCES = sum.f90

sum_LDADD = $(top_builddir)/src/libfortrancl.la @CL_LIBS@
dist_noinst_DATA = sum.cl multiply.cl

multiply_SOURCES = numeric_constants.f90 modules.f90 collision_operator.f90 r1mach.f90 rf.f rd.f ellipf.f90 cl_collop.f90 multiply.f90
multiply_LDADD = $(top_builddir)/src/libfortrancl.la @CL_LIBS@

devices_SOURCES = devices.f90
devices_LDADD = $(top_builddir)/src/libfortrancl.la @CL_LIBS@

AM_FCFLAGS = @F90_MODULE_FLAG@$(top_builddir)/src -g -O0
# AM_FFLAGS = @F77_MODULE_FLAG@$(top_builddir)/src

CLEANFILES = *~ *.bak *.mod *.MOD *.il *.d *.pc* ifc* $(noinst_PROGRAMS)

