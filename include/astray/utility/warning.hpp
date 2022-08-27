#if defined(_MSC_VER)
#define DISABLE_WARNING_PUSH           __pragma(warning( push ))
#define DISABLE_WARNING_POP            __pragma(warning( pop )) 
#define DISABLE_WARNING(warningNumber) __pragma(warning( disable : warningNumber ))

#define DISABLE_WARNING_PRAGMAS
#define DISABLE_WARNING_NO_DEPRECATED             DISABLE_WARNING(4996)
#define DISABLE_WARNING_NO_DEPRECATED_GPU_TARGETS 

#elif defined(__GNUC__) || defined(__clang__)
#define DO_PRAGMA(X) _Pragma(#X)
#define DISABLE_WARNING_PUSH           DO_PRAGMA(GCC diagnostic push)
#define DISABLE_WARNING_POP            DO_PRAGMA(GCC diagnostic pop) 
#define DISABLE_WARNING(warningName)   DO_PRAGMA(GCC diagnostic ignored #warningName)

#define DISABLE_WARNING_PRAGMAS                   DISABLE_WARNING(-Wpragmas)
#define DISABLE_WARNING_NO_DEPRECATED             DISABLE_WARNING(-Wno-deprecated)
#define DISABLE_WARNING_NO_DEPRECATED_GPU_TARGETS DISABLE_WARNING(-Wno-deprecated-gpu-targets)

#else
#define DISABLE_WARNING_PUSH
#define DISABLE_WARNING_POP

#define DISABLE_WARNING_PRAGMAS
#define DISABLE_WARNING_NO_DEPRECATED
#define DISABLE_WARNING_NO_DEPRECATED_GPU_TARGETS

#endif
