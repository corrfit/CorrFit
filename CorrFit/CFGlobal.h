#ifndef _CF_GLOBAL_H_
#define _CF_GLOBAL_H_

// Define global types

using namespace std;

struct _FourVector 
{
  double x, y, z, t;
};

typedef struct _FourVector FourVector;

struct _VectorPair 
{
  FourVector part1, part2;
};

typedef struct _VectorPair VectorPair;

// Define compilation specific variables

#define PRINT_MESSAGE(_mes) cout << _mes << endl;
#define D_(_mes) cout << _mes << endl;

#ifdef _DEBUG_
#define PRINT_DEBUG(_mes) cerr << _mes << endl;
#else
#define PRINT_DEBUG(_mes) 
#endif

#if _DEBUG_LEVEL_==0
#define PRINT_DEBUG_3(_mes) 
#define PRINT_DEBUG_2(_mes) 
#define PRINT_DEBUG_1(_mes) 
#elif _DEBUG_LEVEL_==1
#define PRINT_DEBUG_3(_mes) 
#define PRINT_DEBUG_2(_mes) 
#define PRINT_DEBUG_1(_mes) cerr << _mes << endl;
#elif _DEBUG_LEVEL_==2
#define PRINT_DEBUG_3(_mes) 
#define PRINT_DEBUG_2(_mes) cerr << _mes << endl;
#define PRINT_DEBUG_1(_mes) cerr << _mes << endl;
#elif _DEBUG_LEVEL_==3
#define PRINT_DEBUG_3(_mes) cerr << _mes << endl;
#define PRINT_DEBUG_2(_mes) cerr << _mes << endl;
#define PRINT_DEBUG_1(_mes) cerr << _mes << endl;
#endif

#ifdef _GCC2_
#define STDIOS ios
#endif

#ifdef _GCC3_
#define STDIOS ios_base
#endif

#endif
