Place the public header files in this directory. They will be
available to your code (and other modules) with

     #include <IMP/centrioles/myheader.h>

All headers should include `IMP/centrioles/centrioles_config.h` as their
first include and surround all code with `IMPCENTRIOLES_BEGIN_NAMESPACE`
and `IMPCENTRIOLES_END_NAMESPACE` to put it in the IMP::centrioles namespace
and manage compiler warnings.

Headers should also be exposed to SWIG in the `pyext/swig.i-in` file.
