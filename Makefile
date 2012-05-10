# Makefile.in generated by automake 1.10 from Makefile.am.
# examples/Makefile.  Generated from Makefile.in by configure.

# Copyright (C) 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
# 2003, 2004, 2005, 2006  Free Software Foundation, Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.





pkgdatadir = $(datadir)/fortrancl
pkglibdir = $(libdir)/fortrancl
pkgincludedir = $(includedir)/fortrancl
am__cd = CDPATH="$${ZSH_VERSION+.}$(PATH_SEPARATOR)" && cd
install_sh_DATA = $(install_sh) -c -m 644
install_sh_PROGRAM = $(install_sh) -c
install_sh_SCRIPT = $(install_sh) -c
INSTALL_HEADER = $(INSTALL_DATA)
transform = $(program_transform_name)
NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
build_triplet = i386-apple-darwin11.3.0
host_triplet = i386-apple-darwin11.3.0
noinst_PROGRAMS = sum$(EXEEXT) devices$(EXEEXT) multiply$(EXEEXT)
subdir = examples
DIST_COMMON = $(dist_noinst_DATA) $(srcdir)/Makefile.am \
	$(srcdir)/Makefile.in
ACLOCAL_M4 = $(top_srcdir)/aclocal.m4
am__aclocal_m4_deps = $(top_srcdir)/m4/acx_pthread.m4 \
	$(top_srcdir)/m4/ax_check_cl.m4 \
	$(top_srcdir)/m4/ax_lang_compiler_ms.m4 \
	$(top_srcdir)/m4/f90_module_extension.m4 \
	$(top_srcdir)/m4/f90_module_flag.m4 \
	$(top_srcdir)/m4/fcflags.m4 $(top_srcdir)/m4/libtool.m4 \
	$(top_srcdir)/m4/ltoptions.m4 $(top_srcdir)/m4/ltsugar.m4 \
	$(top_srcdir)/m4/ltversion.m4 $(top_srcdir)/m4/lt~obsolete.m4 \
	$(top_srcdir)/configure.ac
am__configure_deps = $(am__aclocal_m4_deps) $(CONFIGURE_DEPENDENCIES) \
	$(ACLOCAL_M4)
mkinstalldirs = $(install_sh) -d
CONFIG_HEADER = $(top_builddir)/config.h
CONFIG_CLEAN_FILES =
PROGRAMS = $(noinst_PROGRAMS)
am_devices_OBJECTS = devices.$(OBJEXT)
devices_OBJECTS = $(am_devices_OBJECTS)
devices_DEPENDENCIES = $(top_builddir)/src/libfortrancl.la
am_multiply_OBJECTS = numeric_constants.$(OBJEXT) modules.$(OBJEXT) \
	collision_operator.$(OBJEXT) r1mach.$(OBJEXT) rf.$(OBJEXT) \
	rd.$(OBJEXT) ellipf.$(OBJEXT) cl_collop.$(OBJEXT) \
	multiply.$(OBJEXT)
multiply_OBJECTS = $(am_multiply_OBJECTS)
multiply_DEPENDENCIES = $(top_builddir)/src/libfortrancl.la
am_sum_OBJECTS = sum.$(OBJEXT)
sum_OBJECTS = $(am_sum_OBJECTS)
sum_DEPENDENCIES = $(top_builddir)/src/libfortrancl.la
DEFAULT_INCLUDES = -I. -I$(top_builddir)
F77COMPILE = $(F77) $(AM_FFLAGS) $(FFLAGS)
LTF77COMPILE = $(LIBTOOL) --tag=F77 $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=compile $(F77) $(AM_FFLAGS) $(FFLAGS)
F77LD = $(F77)
F77LINK = $(LIBTOOL) --tag=F77 $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=link $(F77LD) $(AM_FFLAGS) $(FFLAGS) $(AM_LDFLAGS) \
	$(LDFLAGS) -o $@
FCCOMPILE = $(FC) $(AM_FCFLAGS) $(FCFLAGS)
LTFCCOMPILE = $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=compile $(FC) $(AM_FCFLAGS) $(FCFLAGS)
FCLD = $(FC)
FCLINK = $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=link \
	$(FCLD) $(AM_FCFLAGS) $(FCFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o \
	$@
SOURCES = $(devices_SOURCES) $(multiply_SOURCES) $(sum_SOURCES)
DIST_SOURCES = $(devices_SOURCES) $(multiply_SOURCES) $(sum_SOURCES)
DATA = $(dist_noinst_DATA)
ETAGS = etags
CTAGS = ctags
DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
ACLOCAL = ${SHELL} /Users/jakob/OpenCL/fortrancl/missing --run aclocal-1.10
AMTAR = ${SHELL} /Users/jakob/OpenCL/fortrancl/missing --run tar
AM_CXXFLAGS =  -O3
AR = ar
AUTOCONF = ${SHELL} /Users/jakob/OpenCL/fortrancl/missing --run autoconf
AUTOHEADER = ${SHELL} /Users/jakob/OpenCL/fortrancl/missing --run autoheader
AUTOMAKE = ${SHELL} /Users/jakob/OpenCL/fortrancl/missing --run automake-1.10
AWK = awk
CC = gcc
CCDEPMODE = depmode=gcc3
CFLAGS = -I/usr/local/cuda/include/ -O3
CL_CFLAGS = -D_THREAD_SAFE 
CL_LIBS = -framework OpenCL -L/usr/X11/lib -lX11  -lm
CPP = gcc -E
CPPFLAGS =  -I/usr/include
CYGPATH_W = echo
DEFS = -DHAVE_CONFIG_H
DEPDIR = .deps
DSYMUTIL = dsymutil
DUMPBIN = 
ECHO_C = \c
ECHO_N = 
ECHO_T = 
EGREP = /usr/bin/grep -E
EXEEXT = 
F77 = gfortran
F90_MODULE_FLAG = -I 
FC = gfortran
FCFLAGS = -g -O2
FCFLAGS_f90 = 
FCLIBS =  -L/usr/X11/lib -L/usr/local/lib/gcc/x86_64-apple-darwin11.0.0/4.6.1 -L/usr/local/lib/gcc/x86_64-apple-darwin11.0.0/4.6.1/../../.. -lX11 -lgfortran -lquadmath -lm
FFLAGS = -g -O2
FGREP = /usr/bin/grep -F
FORTRANCL_SO_VERSION = 0:0:0
GREP = /usr/bin/grep
INSTALL = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_PROGRAM = ${INSTALL}
INSTALL_SCRIPT = ${INSTALL}
INSTALL_STRIP_PROGRAM = $(install_sh) -c -s
LD = /usr/llvm-gcc-4.2/libexec/gcc/i686-apple-darwin11/4.2.1/ld
LDFLAGS = 
LIBOBJS =  ${LIBOBJDIR}lstat$U.o
LIBS =   -L/usr/X11/lib -L/usr/local/lib/gcc/x86_64-apple-darwin11.0.0/4.6.1 -L/usr/local/lib/gcc/x86_64-apple-darwin11.0.0/4.6.1/../../.. -lX11 -lgfortran -lquadmath -lm
LIBTOOL = $(SHELL) $(top_builddir)/libtool
LIPO = lipo
LN_S = ln -s
LTLIBOBJS =  ${LIBOBJDIR}lstat$U.lo
MAKEINFO = ${SHELL} /Users/jakob/OpenCL/fortrancl/missing --run makeinfo
MKDIR_P = .././install-sh -c -d
NM = /usr/bin/nm
NMEDIT = nmedit
OBJDUMP = false
OBJEXT = o
OTOOL = otool
OTOOL64 = :
PACKAGE = fortrancl
PACKAGE_BUGREPORT = xavier@tddft.org
PACKAGE_NAME = FortranCL
PACKAGE_STRING = FortranCL 0.1alpha3
PACKAGE_TARNAME = fortrancl
PACKAGE_VERSION = 0.1alpha3
PATH_SEPARATOR = :
PTHREAD_CC = gcc
PTHREAD_CFLAGS = -D_THREAD_SAFE 
PTHREAD_LIBS = 
RANLIB = ranlib
SED = /usr/bin/sed
SET_MAKE = 
SHELL = /bin/sh
STRIP = strip
VERSION = 0.1alpha3
XMKMF = 
YACC = bison -y
YFLAGS = 
abs_builddir = /Users/jakob/OpenCL/fortrancl/examples
abs_srcdir = /Users/jakob/OpenCL/fortrancl/examples
abs_top_builddir = /Users/jakob/OpenCL/fortrancl
abs_top_srcdir = /Users/jakob/OpenCL/fortrancl
ac_ct_CC = gcc
ac_ct_DUMPBIN = 
ac_ct_F77 = gfortran
ac_ct_FC = gfortran
acx_pthread_config = 
am__include = include
am__leading_dot = .
am__quote = 
am__tar = ${AMTAR} chof - "$$tardir"
am__untar = ${AMTAR} xf -
ax_cv_f90_modext = mod
bindir = ${exec_prefix}/bin
build = i386-apple-darwin11.3.0
build_alias = 
build_cpu = i386
build_os = darwin11.3.0
build_vendor = apple
builddir = .
datadir = ${datarootdir}
datarootdir = ${prefix}/share
docdir = ${datarootdir}/doc/${PACKAGE_TARNAME}
dvidir = ${docdir}
exec_prefix = ${prefix}
host = i386-apple-darwin11.3.0
host_alias = 
host_cpu = i386
host_os = darwin11.3.0
host_vendor = apple
htmldir = ${docdir}
includedir = ${prefix}/include
infodir = ${datarootdir}/info
install_sh = $(SHELL) /Users/jakob/OpenCL/fortrancl/install-sh
libdir = ${exec_prefix}/lib
libexecdir = ${exec_prefix}/libexec
localedir = ${datarootdir}/locale
localstatedir = ${prefix}/var
mandir = ${datarootdir}/man
mkdir_p = $(top_builddir)/./install-sh -c -d
oldincludedir = /usr/include
pdfdir = ${docdir}
prefix = /usr/local
program_transform_name = s,x,x,
psdir = ${docdir}
sbindir = ${exec_prefix}/sbin
sharedstatedir = ${prefix}/com
srcdir = .
sysconfdir = ${prefix}/etc
target_alias = 
top_builddir = ..
top_srcdir = ..
sum_SOURCES = sum.f90
sum_LDADD = $(top_builddir)/src/libfortrancl.la -framework OpenCL -L/usr/X11/lib -lX11  -lm
dist_noinst_DATA = sum.cl multiply.cl
multiply_SOURCES = numeric_constants.f90 modules.f90 collision_operator.f90 r1mach.f90 rf.f rd.f ellipf.f90 cl_collop.f90 multiply.f90
multiply_LDADD = $(top_builddir)/src/libfortrancl.la -framework OpenCL -L/usr/X11/lib -lX11  -lm
devices_SOURCES = devices.f90
devices_LDADD = $(top_builddir)/src/libfortrancl.la -framework OpenCL -L/usr/X11/lib -lX11  -lm
AM_FCFLAGS = -I $(top_builddir)/src -g -O0
# AM_FFLAGS = @F77_MODULE_FLAG@$(top_builddir)/src
CLEANFILES = *~ *.bak *.mod *.MOD *.il *.d *.pc* ifc* $(noinst_PROGRAMS)
all: all-am

.SUFFIXES:
.SUFFIXES: .f .f90 .lo .o .obj
$(srcdir)/Makefile.in:  $(srcdir)/Makefile.am  $(am__configure_deps)
	@for dep in $?; do \
	  case '$(am__configure_deps)' in \
	    *$$dep*) \
	      cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh \
		&& exit 0; \
	      exit 1;; \
	  esac; \
	done; \
	echo ' cd $(top_srcdir) && $(AUTOMAKE) --gnu  examples/Makefile'; \
	cd $(top_srcdir) && \
	  $(AUTOMAKE) --gnu  examples/Makefile
.PRECIOUS: Makefile
Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@case '$?' in \
	  *config.status*) \
	    cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh;; \
	  *) \
	    echo ' cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe)'; \
	    cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe);; \
	esac;

$(top_builddir)/config.status: $(top_srcdir)/configure $(CONFIG_STATUS_DEPENDENCIES)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh

$(top_srcdir)/configure:  $(am__configure_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh
$(ACLOCAL_M4):  $(am__aclocal_m4_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh

clean-noinstPROGRAMS:
	@list='$(noinst_PROGRAMS)'; for p in $$list; do \
	  f=`echo $$p|sed 's/$(EXEEXT)$$//'`; \
	  echo " rm -f $$p $$f"; \
	  rm -f $$p $$f ; \
	done
devices$(EXEEXT): $(devices_OBJECTS) $(devices_DEPENDENCIES) 
	@rm -f devices$(EXEEXT)
	$(FCLINK) $(devices_OBJECTS) $(devices_LDADD) $(LIBS)
multiply$(EXEEXT): $(multiply_OBJECTS) $(multiply_DEPENDENCIES) 
	@rm -f multiply$(EXEEXT)
	$(F77LINK) $(multiply_OBJECTS) $(multiply_LDADD) $(LIBS)
sum$(EXEEXT): $(sum_OBJECTS) $(sum_DEPENDENCIES) 
	@rm -f sum$(EXEEXT)
	$(FCLINK) $(sum_OBJECTS) $(sum_LDADD) $(LIBS)

mostlyclean-compile:
	-rm -f *.$(OBJEXT)

distclean-compile:
	-rm -f *.tab.c

.f.o:
	$(F77COMPILE) -c -o $@ $<

.f.obj:
	$(F77COMPILE) -c -o $@ `$(CYGPATH_W) '$<'`

.f.lo:
	$(LTF77COMPILE) -c -o $@ $<

.f90.o:
	$(FCCOMPILE) -c -o $@ $<

.f90.obj:
	$(FCCOMPILE) -c -o $@ `$(CYGPATH_W) '$<'`

.f90.lo:
	$(LTFCCOMPILE) -c -o $@ $<

mostlyclean-libtool:
	-rm -f *.lo

clean-libtool:
	-rm -rf .libs _libs

ID: $(HEADERS) $(SOURCES) $(LISP) $(TAGS_FILES)
	list='$(SOURCES) $(HEADERS) $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	mkid -fID $$unique
tags: TAGS

TAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) \
		$(TAGS_FILES) $(LISP)
	tags=; \
	here=`pwd`; \
	list='$(SOURCES) $(HEADERS)  $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	if test -z "$(ETAGS_ARGS)$$tags$$unique"; then :; else \
	  test -n "$$unique" || unique=$$empty_fix; \
	  $(ETAGS) $(ETAGSFLAGS) $(AM_ETAGSFLAGS) $(ETAGS_ARGS) \
	    $$tags $$unique; \
	fi
ctags: CTAGS
CTAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) \
		$(TAGS_FILES) $(LISP)
	tags=; \
	here=`pwd`; \
	list='$(SOURCES) $(HEADERS)  $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	test -z "$(CTAGS_ARGS)$$tags$$unique" \
	  || $(CTAGS) $(CTAGSFLAGS) $(AM_CTAGSFLAGS) $(CTAGS_ARGS) \
	     $$tags $$unique

GTAGS:
	here=`$(am__cd) $(top_builddir) && pwd` \
	  && cd $(top_srcdir) \
	  && gtags -i $(GTAGS_ARGS) $$here

distclean-tags:
	-rm -f TAGS ID GTAGS GRTAGS GSYMS GPATH tags

distdir: $(DISTFILES)
	@srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	topsrcdirstrip=`echo "$(top_srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	list='$(DISTFILES)'; \
	  dist_files=`for file in $$list; do echo $$file; done | \
	  sed -e "s|^$$srcdirstrip/||;t" \
	      -e "s|^$$topsrcdirstrip/|$(top_builddir)/|;t"`; \
	case $$dist_files in \
	  */*) $(MKDIR_P) `echo "$$dist_files" | \
			   sed '/\//!d;s|^|$(distdir)/|;s,/[^/]*$$,,' | \
			   sort -u` ;; \
	esac; \
	for file in $$dist_files; do \
	  if test -f $$file || test -d $$file; then d=.; else d=$(srcdir); fi; \
	  if test -d $$d/$$file; then \
	    dir=`echo "/$$file" | sed -e 's,/[^/]*$$,,'`; \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -pR $(srcdir)/$$file $(distdir)$$dir || exit 1; \
	    fi; \
	    cp -pR $$d/$$file $(distdir)$$dir || exit 1; \
	  else \
	    test -f $(distdir)/$$file \
	    || cp -p $$d/$$file $(distdir)/$$file \
	    || exit 1; \
	  fi; \
	done
check-am: all-am
check: check-am
all-am: Makefile $(PROGRAMS) $(DATA)
installdirs:
install: install-am
install-exec: install-exec-am
install-data: install-data-am
uninstall: uninstall-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am

installcheck: installcheck-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) INSTALL_PROGRAM="$(INSTALL_STRIP_PROGRAM)" \
	  install_sh_PROGRAM="$(INSTALL_STRIP_PROGRAM)" INSTALL_STRIP_FLAG=-s \
	  `test -z '$(STRIP)' || \
	    echo "INSTALL_PROGRAM_ENV=STRIPPROG='$(STRIP)'"` install
mostlyclean-generic:

clean-generic:
	-test -z "$(CLEANFILES)" || rm -f $(CLEANFILES)

distclean-generic:
	-test -z "$(CONFIG_CLEAN_FILES)" || rm -f $(CONFIG_CLEAN_FILES)

maintainer-clean-generic:
	@echo "This command is intended for maintainers to use"
	@echo "it deletes files that may require special tools to rebuild."
clean: clean-am

clean-am: clean-generic clean-libtool clean-noinstPROGRAMS \
	mostlyclean-am

distclean: distclean-am
	-rm -f Makefile
distclean-am: clean-am distclean-compile distclean-generic \
	distclean-tags

dvi: dvi-am

dvi-am:

html: html-am

info: info-am

info-am:

install-data-am:

install-dvi: install-dvi-am

install-exec-am:

install-html: install-html-am

install-info: install-info-am

install-man:

install-pdf: install-pdf-am

install-ps: install-ps-am

installcheck-am:

maintainer-clean: maintainer-clean-am
	-rm -f Makefile
maintainer-clean-am: distclean-am maintainer-clean-generic

mostlyclean: mostlyclean-am

mostlyclean-am: mostlyclean-compile mostlyclean-generic \
	mostlyclean-libtool

pdf: pdf-am

pdf-am:

ps: ps-am

ps-am:

uninstall-am:

.MAKE: install-am install-strip

.PHONY: CTAGS GTAGS all all-am check check-am clean clean-generic \
	clean-libtool clean-noinstPROGRAMS ctags distclean \
	distclean-compile distclean-generic distclean-libtool \
	distclean-tags distdir dvi dvi-am html html-am info info-am \
	install install-am install-data install-data-am install-dvi \
	install-dvi-am install-exec install-exec-am install-html \
	install-html-am install-info install-info-am install-man \
	install-pdf install-pdf-am install-ps install-ps-am \
	install-strip installcheck installcheck-am installdirs \
	maintainer-clean maintainer-clean-generic mostlyclean \
	mostlyclean-compile mostlyclean-generic mostlyclean-libtool \
	pdf pdf-am ps ps-am tags uninstall uninstall-am

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
