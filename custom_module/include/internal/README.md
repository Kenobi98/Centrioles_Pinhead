Place the private header files in this directory. They will be
available to your code with

     #include <IMP/centrioles/internal/myheader.h>

All headers should include `IMP/centrioles/centrioles_config.h` as their
first include and surround all code with `IMPCENTRIOLES_BEGIN_INTERNAL_NAMESPACE`
and `IMPCENTRIOLES_END_INTERNAL_NAMESPACE` to put it in the
IMP::centrioles::internal namespace and manage compiler warnings.
