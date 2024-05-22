/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output, and Bison version.  */
#define YYBISON 30802

/* Bison version string.  */
#define YYBISON_VERSION "3.8.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         lefyyparse
#define yylex           lefyylex
#define yyerror         lefyyerror
#define yydebug         lefyydebug
#define yynerrs         lefyynerrs
#define yylval          lefyylval
#define yychar          lefyychar

/* First part of user prologue.  */
#line 52 "lef/lef.y"

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "lex.h"
#include "lefiDefs.hpp"
#include "lefiUser.hpp"
#include "lefiUtil.hpp"

#include "lefrData.hpp"
#include "lefrCallBacks.hpp"
#include "lefrSettings.hpp"

BEGIN_LEFDEF_PARSER_NAMESPACE

#define LYPROP_ECAP "EDGE_CAPACITANCE"

#define YYINITDEPTH 10000  // pcr 640902 - initialize the yystacksize to 300 
                           // this may need to increase in a design gets 
                           // larger and a polygon has around 300 sizes 
                           // 11/21/2003 - incrreased to 500, design from 
                           // Artisan is greater than 300, need to find a 
                           // way to dynamically increase the size 
                           // 2/10/2004 - increased to 1000 for pcr 686073 
                           // 3/22/2004 - increased to 2000 for pcr 695879 
                           // 9/29/2004 - double the size for pcr 746865 
                           // tried to overwrite the yyoverflow definition 
                           // it is impossible due to the union structure 
                           // 10/03/2006 - increased to 10000 for pcr 913695 

#define YYMAXDEPTH 300000  // 1/24/2008 - increased from 150000 
                           // This value has to be greater than YYINITDEPTH 


// Macro to describe how we handle a callback.
// If the function was set then call it.
// If the function returns non zero then there was an error
// so call the error routine and exit.
#define CALLBACK(func, typ, data) \
    if (!lefData->lef_errors) { \
      if (func) { \
        if ((lefData->lefRetVal = (*func)(typ, data, lefSettings->UserData)) == 0) { \
        } else { \
          return lefData->lefRetVal; \
        } \
      } \
    }

#define CHKERR() \
    if (lefData->lef_errors > 20) { \
      lefError(1020, "Too many syntax errors."); \
      lefData->lef_errors = 0; \
      return 1; \
    }

// **********************************************************************
// **********************************************************************

#define C_EQ 0
#define C_NE 1
#define C_LT 2
#define C_LE 3
#define C_GT 4
#define C_GE 5


int comp_str(char *s1, int op, char *s2)
{
    int k = strcmp(s1, s2);
    switch (op) {
        case C_EQ: return k == 0;
        case C_NE: return k != 0;
        case C_GT: return k >  0;
        case C_GE: return k >= 0;
        case C_LT: return k <  0;
        case C_LE: return k <= 0;
        }
    return 0;
}
int comp_num(double s1, int op, double s2)
{
    double k = s1 - s2;
    switch (op) {
        case C_EQ: return k == 0;
        case C_NE: return k != 0;
        case C_GT: return k >  0;
        case C_GE: return k >= 0;
        case C_LT: return k <  0;
        case C_LE: return k <= 0;
        }
    return 0;
}

int validNum(int values) {
    switch (values) {
        case 100:
        case 200:
        case 1000:
        case 2000:
             return 1;
        case 400:
        case 800:
        case 4000:
        case 8000:
        case 10000:
        case 20000:
             if (lefData->versionNum < 5.6) {
                if (lefCallbacks->UnitsCbk) {
                  if (lefData->unitsWarnings++ < lefSettings->UnitsWarnings) {
                    lefData->outMsg = (char*)lefMalloc(10000);
                    sprintf (lefData->outMsg,
                       "Error found when processing LEF file '%s'\nUnit %d is a version 5.6 or later syntax\nYour lef file is defined with version %g.",
                    lefData->lefrFileName, values, lefData->versionNum);
                    lefError(1501, lefData->outMsg);
                    lefFree(lefData->outMsg);
                  }
                }
                return 0;
             } else {
                return 1;
             }        
    }
    if (lefData->unitsWarnings++ < lefSettings->UnitsWarnings) {
       lefData->outMsg = (char*)lefMalloc(10000);
       sprintf (lefData->outMsg,
          "The value %d defined for LEF UNITS DATABASE MICRONS is invalid\n. Correct value is 100, 200, 400, 800, 1000, 2000, 4000, 8000, 10000, or 20000", values);
       lefError(1502, lefData->outMsg);
       lefFree(lefData->outMsg);
    }
    CHKERR();
    return 0;
}

int zeroOrGt(double values) {
    if (values < 0)
      return 0;
    return 1;
}


#line 220 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

#include "lef.tab.h"
/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_K_HISTORY = 3,                  /* K_HISTORY  */
  YYSYMBOL_K_ABUT = 4,                     /* K_ABUT  */
  YYSYMBOL_K_ABUTMENT = 5,                 /* K_ABUTMENT  */
  YYSYMBOL_K_ACTIVE = 6,                   /* K_ACTIVE  */
  YYSYMBOL_K_ANALOG = 7,                   /* K_ANALOG  */
  YYSYMBOL_K_ARRAY = 8,                    /* K_ARRAY  */
  YYSYMBOL_K_AREA = 9,                     /* K_AREA  */
  YYSYMBOL_K_BLOCK = 10,                   /* K_BLOCK  */
  YYSYMBOL_K_BOTTOMLEFT = 11,              /* K_BOTTOMLEFT  */
  YYSYMBOL_K_BOTTOMRIGHT = 12,             /* K_BOTTOMRIGHT  */
  YYSYMBOL_K_BY = 13,                      /* K_BY  */
  YYSYMBOL_K_CAPACITANCE = 14,             /* K_CAPACITANCE  */
  YYSYMBOL_K_CAPMULTIPLIER = 15,           /* K_CAPMULTIPLIER  */
  YYSYMBOL_K_CLASS = 16,                   /* K_CLASS  */
  YYSYMBOL_K_CLOCK = 17,                   /* K_CLOCK  */
  YYSYMBOL_K_CLOCKTYPE = 18,               /* K_CLOCKTYPE  */
  YYSYMBOL_K_COLUMNMAJOR = 19,             /* K_COLUMNMAJOR  */
  YYSYMBOL_K_DESIGNRULEWIDTH = 20,         /* K_DESIGNRULEWIDTH  */
  YYSYMBOL_K_INFLUENCE = 21,               /* K_INFLUENCE  */
  YYSYMBOL_K_CORE = 22,                    /* K_CORE  */
  YYSYMBOL_K_CORNER = 23,                  /* K_CORNER  */
  YYSYMBOL_K_COVER = 24,                   /* K_COVER  */
  YYSYMBOL_K_CPERSQDIST = 25,              /* K_CPERSQDIST  */
  YYSYMBOL_K_CURRENT = 26,                 /* K_CURRENT  */
  YYSYMBOL_K_CURRENTSOURCE = 27,           /* K_CURRENTSOURCE  */
  YYSYMBOL_K_CUT = 28,                     /* K_CUT  */
  YYSYMBOL_K_DEFAULT = 29,                 /* K_DEFAULT  */
  YYSYMBOL_K_DATABASE = 30,                /* K_DATABASE  */
  YYSYMBOL_K_DATA = 31,                    /* K_DATA  */
  YYSYMBOL_K_DIELECTRIC = 32,              /* K_DIELECTRIC  */
  YYSYMBOL_K_DIRECTION = 33,               /* K_DIRECTION  */
  YYSYMBOL_K_DO = 34,                      /* K_DO  */
  YYSYMBOL_K_EDGECAPACITANCE = 35,         /* K_EDGECAPACITANCE  */
  YYSYMBOL_K_EEQ = 36,                     /* K_EEQ  */
  YYSYMBOL_K_END = 37,                     /* K_END  */
  YYSYMBOL_K_ENDCAP = 38,                  /* K_ENDCAP  */
  YYSYMBOL_K_FALL = 39,                    /* K_FALL  */
  YYSYMBOL_K_FALLCS = 40,                  /* K_FALLCS  */
  YYSYMBOL_K_FALLT0 = 41,                  /* K_FALLT0  */
  YYSYMBOL_K_FALLSATT1 = 42,               /* K_FALLSATT1  */
  YYSYMBOL_K_FALLRS = 43,                  /* K_FALLRS  */
  YYSYMBOL_K_FALLSATCUR = 44,              /* K_FALLSATCUR  */
  YYSYMBOL_K_FALLTHRESH = 45,              /* K_FALLTHRESH  */
  YYSYMBOL_K_FEEDTHRU = 46,                /* K_FEEDTHRU  */
  YYSYMBOL_K_FIXED = 47,                   /* K_FIXED  */
  YYSYMBOL_K_FOREIGN = 48,                 /* K_FOREIGN  */
  YYSYMBOL_K_FROMPIN = 49,                 /* K_FROMPIN  */
  YYSYMBOL_K_GENERATE = 50,                /* K_GENERATE  */
  YYSYMBOL_K_GENERATOR = 51,               /* K_GENERATOR  */
  YYSYMBOL_K_GROUND = 52,                  /* K_GROUND  */
  YYSYMBOL_K_HEIGHT = 53,                  /* K_HEIGHT  */
  YYSYMBOL_K_HORIZONTAL = 54,              /* K_HORIZONTAL  */
  YYSYMBOL_K_INOUT = 55,                   /* K_INOUT  */
  YYSYMBOL_K_INPUT = 56,                   /* K_INPUT  */
  YYSYMBOL_K_INPUTNOISEMARGIN = 57,        /* K_INPUTNOISEMARGIN  */
  YYSYMBOL_K_COMPONENTPIN = 58,            /* K_COMPONENTPIN  */
  YYSYMBOL_K_INTRINSIC = 59,               /* K_INTRINSIC  */
  YYSYMBOL_K_INVERT = 60,                  /* K_INVERT  */
  YYSYMBOL_K_IRDROP = 61,                  /* K_IRDROP  */
  YYSYMBOL_K_ITERATE = 62,                 /* K_ITERATE  */
  YYSYMBOL_K_IV_TABLES = 63,               /* K_IV_TABLES  */
  YYSYMBOL_K_LAYER = 64,                   /* K_LAYER  */
  YYSYMBOL_K_LEAKAGE = 65,                 /* K_LEAKAGE  */
  YYSYMBOL_K_LEQ = 66,                     /* K_LEQ  */
  YYSYMBOL_K_LIBRARY = 67,                 /* K_LIBRARY  */
  YYSYMBOL_K_MACRO = 68,                   /* K_MACRO  */
  YYSYMBOL_K_MATCH = 69,                   /* K_MATCH  */
  YYSYMBOL_K_MAXDELAY = 70,                /* K_MAXDELAY  */
  YYSYMBOL_K_MAXLOAD = 71,                 /* K_MAXLOAD  */
  YYSYMBOL_K_METALOVERHANG = 72,           /* K_METALOVERHANG  */
  YYSYMBOL_K_MILLIAMPS = 73,               /* K_MILLIAMPS  */
  YYSYMBOL_K_MILLIWATTS = 74,              /* K_MILLIWATTS  */
  YYSYMBOL_K_MINFEATURE = 75,              /* K_MINFEATURE  */
  YYSYMBOL_K_MUSTJOIN = 76,                /* K_MUSTJOIN  */
  YYSYMBOL_K_NAMESCASESENSITIVE = 77,      /* K_NAMESCASESENSITIVE  */
  YYSYMBOL_K_NANOSECONDS = 78,             /* K_NANOSECONDS  */
  YYSYMBOL_K_NETS = 79,                    /* K_NETS  */
  YYSYMBOL_K_NEW = 80,                     /* K_NEW  */
  YYSYMBOL_K_NONDEFAULTRULE = 81,          /* K_NONDEFAULTRULE  */
  YYSYMBOL_K_NONINVERT = 82,               /* K_NONINVERT  */
  YYSYMBOL_K_NONUNATE = 83,                /* K_NONUNATE  */
  YYSYMBOL_K_OBS = 84,                     /* K_OBS  */
  YYSYMBOL_K_OHMS = 85,                    /* K_OHMS  */
  YYSYMBOL_K_OFFSET = 86,                  /* K_OFFSET  */
  YYSYMBOL_K_ORIENTATION = 87,             /* K_ORIENTATION  */
  YYSYMBOL_K_ORIGIN = 88,                  /* K_ORIGIN  */
  YYSYMBOL_K_OUTPUT = 89,                  /* K_OUTPUT  */
  YYSYMBOL_K_OUTPUTNOISEMARGIN = 90,       /* K_OUTPUTNOISEMARGIN  */
  YYSYMBOL_K_OVERHANG = 91,                /* K_OVERHANG  */
  YYSYMBOL_K_OVERLAP = 92,                 /* K_OVERLAP  */
  YYSYMBOL_K_OFF = 93,                     /* K_OFF  */
  YYSYMBOL_K_ON = 94,                      /* K_ON  */
  YYSYMBOL_K_OVERLAPS = 95,                /* K_OVERLAPS  */
  YYSYMBOL_K_PAD = 96,                     /* K_PAD  */
  YYSYMBOL_K_PATH = 97,                    /* K_PATH  */
  YYSYMBOL_K_PATTERN = 98,                 /* K_PATTERN  */
  YYSYMBOL_K_PICOFARADS = 99,              /* K_PICOFARADS  */
  YYSYMBOL_K_PIN = 100,                    /* K_PIN  */
  YYSYMBOL_K_PITCH = 101,                  /* K_PITCH  */
  YYSYMBOL_K_PLACED = 102,                 /* K_PLACED  */
  YYSYMBOL_K_POLYGON = 103,                /* K_POLYGON  */
  YYSYMBOL_K_PORT = 104,                   /* K_PORT  */
  YYSYMBOL_K_POST = 105,                   /* K_POST  */
  YYSYMBOL_K_POWER = 106,                  /* K_POWER  */
  YYSYMBOL_K_PRE = 107,                    /* K_PRE  */
  YYSYMBOL_K_PULLDOWNRES = 108,            /* K_PULLDOWNRES  */
  YYSYMBOL_K_RECT = 109,                   /* K_RECT  */
  YYSYMBOL_K_RESISTANCE = 110,             /* K_RESISTANCE  */
  YYSYMBOL_K_RESISTIVE = 111,              /* K_RESISTIVE  */
  YYSYMBOL_K_RING = 112,                   /* K_RING  */
  YYSYMBOL_K_RISE = 113,                   /* K_RISE  */
  YYSYMBOL_K_RISECS = 114,                 /* K_RISECS  */
  YYSYMBOL_K_RISERS = 115,                 /* K_RISERS  */
  YYSYMBOL_K_RISESATCUR = 116,             /* K_RISESATCUR  */
  YYSYMBOL_K_RISETHRESH = 117,             /* K_RISETHRESH  */
  YYSYMBOL_K_RISESATT1 = 118,              /* K_RISESATT1  */
  YYSYMBOL_K_RISET0 = 119,                 /* K_RISET0  */
  YYSYMBOL_K_RISEVOLTAGETHRESHOLD = 120,   /* K_RISEVOLTAGETHRESHOLD  */
  YYSYMBOL_K_FALLVOLTAGETHRESHOLD = 121,   /* K_FALLVOLTAGETHRESHOLD  */
  YYSYMBOL_K_ROUTING = 122,                /* K_ROUTING  */
  YYSYMBOL_K_ROWMAJOR = 123,               /* K_ROWMAJOR  */
  YYSYMBOL_K_RPERSQ = 124,                 /* K_RPERSQ  */
  YYSYMBOL_K_SAMENET = 125,                /* K_SAMENET  */
  YYSYMBOL_K_SCANUSE = 126,                /* K_SCANUSE  */
  YYSYMBOL_K_SHAPE = 127,                  /* K_SHAPE  */
  YYSYMBOL_K_SHRINKAGE = 128,              /* K_SHRINKAGE  */
  YYSYMBOL_K_SIGNAL = 129,                 /* K_SIGNAL  */
  YYSYMBOL_K_SITE = 130,                   /* K_SITE  */
  YYSYMBOL_K_SIZE = 131,                   /* K_SIZE  */
  YYSYMBOL_K_SOURCE = 132,                 /* K_SOURCE  */
  YYSYMBOL_K_SPACER = 133,                 /* K_SPACER  */
  YYSYMBOL_K_SPACING = 134,                /* K_SPACING  */
  YYSYMBOL_K_SPECIALNETS = 135,            /* K_SPECIALNETS  */
  YYSYMBOL_K_STACK = 136,                  /* K_STACK  */
  YYSYMBOL_K_START = 137,                  /* K_START  */
  YYSYMBOL_K_STEP = 138,                   /* K_STEP  */
  YYSYMBOL_K_STOP = 139,                   /* K_STOP  */
  YYSYMBOL_K_STRUCTURE = 140,              /* K_STRUCTURE  */
  YYSYMBOL_K_SYMMETRY = 141,               /* K_SYMMETRY  */
  YYSYMBOL_K_TABLE = 142,                  /* K_TABLE  */
  YYSYMBOL_K_THICKNESS = 143,              /* K_THICKNESS  */
  YYSYMBOL_K_TIEHIGH = 144,                /* K_TIEHIGH  */
  YYSYMBOL_K_TIELOW = 145,                 /* K_TIELOW  */
  YYSYMBOL_K_TIEOFFR = 146,                /* K_TIEOFFR  */
  YYSYMBOL_K_TIME = 147,                   /* K_TIME  */
  YYSYMBOL_K_TIMING = 148,                 /* K_TIMING  */
  YYSYMBOL_K_TO = 149,                     /* K_TO  */
  YYSYMBOL_K_TOPIN = 150,                  /* K_TOPIN  */
  YYSYMBOL_K_TOPLEFT = 151,                /* K_TOPLEFT  */
  YYSYMBOL_K_TOPRIGHT = 152,               /* K_TOPRIGHT  */
  YYSYMBOL_K_TOPOFSTACKONLY = 153,         /* K_TOPOFSTACKONLY  */
  YYSYMBOL_K_TRISTATE = 154,               /* K_TRISTATE  */
  YYSYMBOL_K_TYPE = 155,                   /* K_TYPE  */
  YYSYMBOL_K_UNATENESS = 156,              /* K_UNATENESS  */
  YYSYMBOL_K_UNITS = 157,                  /* K_UNITS  */
  YYSYMBOL_K_USE = 158,                    /* K_USE  */
  YYSYMBOL_K_VARIABLE = 159,               /* K_VARIABLE  */
  YYSYMBOL_K_VERTICAL = 160,               /* K_VERTICAL  */
  YYSYMBOL_K_VHI = 161,                    /* K_VHI  */
  YYSYMBOL_K_VIA = 162,                    /* K_VIA  */
  YYSYMBOL_K_VIARULE = 163,                /* K_VIARULE  */
  YYSYMBOL_K_VLO = 164,                    /* K_VLO  */
  YYSYMBOL_K_VOLTAGE = 165,                /* K_VOLTAGE  */
  YYSYMBOL_K_VOLTS = 166,                  /* K_VOLTS  */
  YYSYMBOL_K_WIDTH = 167,                  /* K_WIDTH  */
  YYSYMBOL_K_X = 168,                      /* K_X  */
  YYSYMBOL_K_Y = 169,                      /* K_Y  */
  YYSYMBOL_T_STRING = 170,                 /* T_STRING  */
  YYSYMBOL_QSTRING = 171,                  /* QSTRING  */
  YYSYMBOL_NUMBER = 172,                   /* NUMBER  */
  YYSYMBOL_K_N = 173,                      /* K_N  */
  YYSYMBOL_K_S = 174,                      /* K_S  */
  YYSYMBOL_K_E = 175,                      /* K_E  */
  YYSYMBOL_K_W = 176,                      /* K_W  */
  YYSYMBOL_K_FN = 177,                     /* K_FN  */
  YYSYMBOL_K_FS = 178,                     /* K_FS  */
  YYSYMBOL_K_FE = 179,                     /* K_FE  */
  YYSYMBOL_K_FW = 180,                     /* K_FW  */
  YYSYMBOL_K_R0 = 181,                     /* K_R0  */
  YYSYMBOL_K_R90 = 182,                    /* K_R90  */
  YYSYMBOL_K_R180 = 183,                   /* K_R180  */
  YYSYMBOL_K_R270 = 184,                   /* K_R270  */
  YYSYMBOL_K_MX = 185,                     /* K_MX  */
  YYSYMBOL_K_MY = 186,                     /* K_MY  */
  YYSYMBOL_K_MXR90 = 187,                  /* K_MXR90  */
  YYSYMBOL_K_MYR90 = 188,                  /* K_MYR90  */
  YYSYMBOL_K_USER = 189,                   /* K_USER  */
  YYSYMBOL_K_MASTERSLICE = 190,            /* K_MASTERSLICE  */
  YYSYMBOL_K_ENDMACRO = 191,               /* K_ENDMACRO  */
  YYSYMBOL_K_ENDMACROPIN = 192,            /* K_ENDMACROPIN  */
  YYSYMBOL_K_ENDVIARULE = 193,             /* K_ENDVIARULE  */
  YYSYMBOL_K_ENDVIA = 194,                 /* K_ENDVIA  */
  YYSYMBOL_K_ENDLAYER = 195,               /* K_ENDLAYER  */
  YYSYMBOL_K_ENDSITE = 196,                /* K_ENDSITE  */
  YYSYMBOL_K_CANPLACE = 197,               /* K_CANPLACE  */
  YYSYMBOL_K_CANNOTOCCUPY = 198,           /* K_CANNOTOCCUPY  */
  YYSYMBOL_K_TRACKS = 199,                 /* K_TRACKS  */
  YYSYMBOL_K_FLOORPLAN = 200,              /* K_FLOORPLAN  */
  YYSYMBOL_K_GCELLGRID = 201,              /* K_GCELLGRID  */
  YYSYMBOL_K_DEFAULTCAP = 202,             /* K_DEFAULTCAP  */
  YYSYMBOL_K_MINPINS = 203,                /* K_MINPINS  */
  YYSYMBOL_K_WIRECAP = 204,                /* K_WIRECAP  */
  YYSYMBOL_K_STABLE = 205,                 /* K_STABLE  */
  YYSYMBOL_K_SETUP = 206,                  /* K_SETUP  */
  YYSYMBOL_K_HOLD = 207,                   /* K_HOLD  */
  YYSYMBOL_K_DEFINE = 208,                 /* K_DEFINE  */
  YYSYMBOL_K_DEFINES = 209,                /* K_DEFINES  */
  YYSYMBOL_K_DEFINEB = 210,                /* K_DEFINEB  */
  YYSYMBOL_K_IF = 211,                     /* K_IF  */
  YYSYMBOL_K_THEN = 212,                   /* K_THEN  */
  YYSYMBOL_K_ELSE = 213,                   /* K_ELSE  */
  YYSYMBOL_K_FALSE = 214,                  /* K_FALSE  */
  YYSYMBOL_K_TRUE = 215,                   /* K_TRUE  */
  YYSYMBOL_K_EQ = 216,                     /* K_EQ  */
  YYSYMBOL_K_NE = 217,                     /* K_NE  */
  YYSYMBOL_K_LE = 218,                     /* K_LE  */
  YYSYMBOL_K_LT = 219,                     /* K_LT  */
  YYSYMBOL_K_GE = 220,                     /* K_GE  */
  YYSYMBOL_K_GT = 221,                     /* K_GT  */
  YYSYMBOL_K_OR = 222,                     /* K_OR  */
  YYSYMBOL_K_AND = 223,                    /* K_AND  */
  YYSYMBOL_K_NOT = 224,                    /* K_NOT  */
  YYSYMBOL_K_DELAY = 225,                  /* K_DELAY  */
  YYSYMBOL_K_TABLEDIMENSION = 226,         /* K_TABLEDIMENSION  */
  YYSYMBOL_K_TABLEAXIS = 227,              /* K_TABLEAXIS  */
  YYSYMBOL_K_TABLEENTRIES = 228,           /* K_TABLEENTRIES  */
  YYSYMBOL_K_TRANSITIONTIME = 229,         /* K_TRANSITIONTIME  */
  YYSYMBOL_K_EXTENSION = 230,              /* K_EXTENSION  */
  YYSYMBOL_K_PROPDEF = 231,                /* K_PROPDEF  */
  YYSYMBOL_K_STRING = 232,                 /* K_STRING  */
  YYSYMBOL_K_INTEGER = 233,                /* K_INTEGER  */
  YYSYMBOL_K_REAL = 234,                   /* K_REAL  */
  YYSYMBOL_K_RANGE = 235,                  /* K_RANGE  */
  YYSYMBOL_K_PROPERTY = 236,               /* K_PROPERTY  */
  YYSYMBOL_K_VIRTUAL = 237,                /* K_VIRTUAL  */
  YYSYMBOL_K_BUSBITCHARS = 238,            /* K_BUSBITCHARS  */
  YYSYMBOL_K_VERSION = 239,                /* K_VERSION  */
  YYSYMBOL_K_BEGINEXT = 240,               /* K_BEGINEXT  */
  YYSYMBOL_K_ENDEXT = 241,                 /* K_ENDEXT  */
  YYSYMBOL_K_UNIVERSALNOISEMARGIN = 242,   /* K_UNIVERSALNOISEMARGIN  */
  YYSYMBOL_K_EDGERATETHRESHOLD1 = 243,     /* K_EDGERATETHRESHOLD1  */
  YYSYMBOL_K_CORRECTIONTABLE = 244,        /* K_CORRECTIONTABLE  */
  YYSYMBOL_K_EDGERATESCALEFACTOR = 245,    /* K_EDGERATESCALEFACTOR  */
  YYSYMBOL_K_EDGERATETHRESHOLD2 = 246,     /* K_EDGERATETHRESHOLD2  */
  YYSYMBOL_K_VICTIMNOISE = 247,            /* K_VICTIMNOISE  */
  YYSYMBOL_K_NOISETABLE = 248,             /* K_NOISETABLE  */
  YYSYMBOL_K_EDGERATE = 249,               /* K_EDGERATE  */
  YYSYMBOL_K_OUTPUTRESISTANCE = 250,       /* K_OUTPUTRESISTANCE  */
  YYSYMBOL_K_VICTIMLENGTH = 251,           /* K_VICTIMLENGTH  */
  YYSYMBOL_K_CORRECTIONFACTOR = 252,       /* K_CORRECTIONFACTOR  */
  YYSYMBOL_K_OUTPUTPINANTENNASIZE = 253,   /* K_OUTPUTPINANTENNASIZE  */
  YYSYMBOL_K_INPUTPINANTENNASIZE = 254,    /* K_INPUTPINANTENNASIZE  */
  YYSYMBOL_K_INOUTPINANTENNASIZE = 255,    /* K_INOUTPINANTENNASIZE  */
  YYSYMBOL_K_CURRENTDEN = 256,             /* K_CURRENTDEN  */
  YYSYMBOL_K_PWL = 257,                    /* K_PWL  */
  YYSYMBOL_K_ANTENNALENGTHFACTOR = 258,    /* K_ANTENNALENGTHFACTOR  */
  YYSYMBOL_K_TAPERRULE = 259,              /* K_TAPERRULE  */
  YYSYMBOL_K_DIVIDERCHAR = 260,            /* K_DIVIDERCHAR  */
  YYSYMBOL_K_ANTENNASIZE = 261,            /* K_ANTENNASIZE  */
  YYSYMBOL_K_ANTENNAMETALLENGTH = 262,     /* K_ANTENNAMETALLENGTH  */
  YYSYMBOL_K_ANTENNAMETALAREA = 263,       /* K_ANTENNAMETALAREA  */
  YYSYMBOL_K_RISESLEWLIMIT = 264,          /* K_RISESLEWLIMIT  */
  YYSYMBOL_K_FALLSLEWLIMIT = 265,          /* K_FALLSLEWLIMIT  */
  YYSYMBOL_K_FUNCTION = 266,               /* K_FUNCTION  */
  YYSYMBOL_K_BUFFER = 267,                 /* K_BUFFER  */
  YYSYMBOL_K_INVERTER = 268,               /* K_INVERTER  */
  YYSYMBOL_K_NAMEMAPSTRING = 269,          /* K_NAMEMAPSTRING  */
  YYSYMBOL_K_NOWIREEXTENSIONATPIN = 270,   /* K_NOWIREEXTENSIONATPIN  */
  YYSYMBOL_K_WIREEXTENSION = 271,          /* K_WIREEXTENSION  */
  YYSYMBOL_K_MESSAGE = 272,                /* K_MESSAGE  */
  YYSYMBOL_K_CREATEFILE = 273,             /* K_CREATEFILE  */
  YYSYMBOL_K_OPENFILE = 274,               /* K_OPENFILE  */
  YYSYMBOL_K_CLOSEFILE = 275,              /* K_CLOSEFILE  */
  YYSYMBOL_K_WARNING = 276,                /* K_WARNING  */
  YYSYMBOL_K_ERROR = 277,                  /* K_ERROR  */
  YYSYMBOL_K_FATALERROR = 278,             /* K_FATALERROR  */
  YYSYMBOL_K_RECOVERY = 279,               /* K_RECOVERY  */
  YYSYMBOL_K_SKEW = 280,                   /* K_SKEW  */
  YYSYMBOL_K_ANYEDGE = 281,                /* K_ANYEDGE  */
  YYSYMBOL_K_POSEDGE = 282,                /* K_POSEDGE  */
  YYSYMBOL_K_NEGEDGE = 283,                /* K_NEGEDGE  */
  YYSYMBOL_K_SDFCONDSTART = 284,           /* K_SDFCONDSTART  */
  YYSYMBOL_K_SDFCONDEND = 285,             /* K_SDFCONDEND  */
  YYSYMBOL_K_SDFCOND = 286,                /* K_SDFCOND  */
  YYSYMBOL_K_MPWH = 287,                   /* K_MPWH  */
  YYSYMBOL_K_MPWL = 288,                   /* K_MPWL  */
  YYSYMBOL_K_PERIOD = 289,                 /* K_PERIOD  */
  YYSYMBOL_K_ACCURRENTDENSITY = 290,       /* K_ACCURRENTDENSITY  */
  YYSYMBOL_K_DCCURRENTDENSITY = 291,       /* K_DCCURRENTDENSITY  */
  YYSYMBOL_K_AVERAGE = 292,                /* K_AVERAGE  */
  YYSYMBOL_K_PEAK = 293,                   /* K_PEAK  */
  YYSYMBOL_K_RMS = 294,                    /* K_RMS  */
  YYSYMBOL_K_FREQUENCY = 295,              /* K_FREQUENCY  */
  YYSYMBOL_K_CUTAREA = 296,                /* K_CUTAREA  */
  YYSYMBOL_K_MEGAHERTZ = 297,              /* K_MEGAHERTZ  */
  YYSYMBOL_K_USELENGTHTHRESHOLD = 298,     /* K_USELENGTHTHRESHOLD  */
  YYSYMBOL_K_LENGTHTHRESHOLD = 299,        /* K_LENGTHTHRESHOLD  */
  YYSYMBOL_K_ANTENNAINPUTGATEAREA = 300,   /* K_ANTENNAINPUTGATEAREA  */
  YYSYMBOL_K_ANTENNAINOUTDIFFAREA = 301,   /* K_ANTENNAINOUTDIFFAREA  */
  YYSYMBOL_K_ANTENNAOUTPUTDIFFAREA = 302,  /* K_ANTENNAOUTPUTDIFFAREA  */
  YYSYMBOL_K_ANTENNAAREARATIO = 303,       /* K_ANTENNAAREARATIO  */
  YYSYMBOL_K_ANTENNADIFFAREARATIO = 304,   /* K_ANTENNADIFFAREARATIO  */
  YYSYMBOL_K_ANTENNACUMAREARATIO = 305,    /* K_ANTENNACUMAREARATIO  */
  YYSYMBOL_K_ANTENNACUMDIFFAREARATIO = 306, /* K_ANTENNACUMDIFFAREARATIO  */
  YYSYMBOL_K_ANTENNAAREAFACTOR = 307,      /* K_ANTENNAAREAFACTOR  */
  YYSYMBOL_K_ANTENNASIDEAREARATIO = 308,   /* K_ANTENNASIDEAREARATIO  */
  YYSYMBOL_K_ANTENNADIFFSIDEAREARATIO = 309, /* K_ANTENNADIFFSIDEAREARATIO  */
  YYSYMBOL_K_ANTENNACUMSIDEAREARATIO = 310, /* K_ANTENNACUMSIDEAREARATIO  */
  YYSYMBOL_K_ANTENNACUMDIFFSIDEAREARATIO = 311, /* K_ANTENNACUMDIFFSIDEAREARATIO  */
  YYSYMBOL_K_ANTENNASIDEAREAFACTOR = 312,  /* K_ANTENNASIDEAREAFACTOR  */
  YYSYMBOL_K_DIFFUSEONLY = 313,            /* K_DIFFUSEONLY  */
  YYSYMBOL_K_MANUFACTURINGGRID = 314,      /* K_MANUFACTURINGGRID  */
  YYSYMBOL_K_FIXEDMASK = 315,              /* K_FIXEDMASK  */
  YYSYMBOL_K_ANTENNACELL = 316,            /* K_ANTENNACELL  */
  YYSYMBOL_K_CLEARANCEMEASURE = 317,       /* K_CLEARANCEMEASURE  */
  YYSYMBOL_K_EUCLIDEAN = 318,              /* K_EUCLIDEAN  */
  YYSYMBOL_K_MAXXY = 319,                  /* K_MAXXY  */
  YYSYMBOL_K_USEMINSPACING = 320,          /* K_USEMINSPACING  */
  YYSYMBOL_K_ROWMINSPACING = 321,          /* K_ROWMINSPACING  */
  YYSYMBOL_K_ROWABUTSPACING = 322,         /* K_ROWABUTSPACING  */
  YYSYMBOL_K_FLIP = 323,                   /* K_FLIP  */
  YYSYMBOL_K_NONE = 324,                   /* K_NONE  */
  YYSYMBOL_K_ANTENNAPARTIALMETALAREA = 325, /* K_ANTENNAPARTIALMETALAREA  */
  YYSYMBOL_K_ANTENNAPARTIALMETALSIDEAREA = 326, /* K_ANTENNAPARTIALMETALSIDEAREA  */
  YYSYMBOL_K_ANTENNAGATEAREA = 327,        /* K_ANTENNAGATEAREA  */
  YYSYMBOL_K_ANTENNADIFFAREA = 328,        /* K_ANTENNADIFFAREA  */
  YYSYMBOL_K_ANTENNAMAXAREACAR = 329,      /* K_ANTENNAMAXAREACAR  */
  YYSYMBOL_K_ANTENNAMAXSIDEAREACAR = 330,  /* K_ANTENNAMAXSIDEAREACAR  */
  YYSYMBOL_K_ANTENNAPARTIALCUTAREA = 331,  /* K_ANTENNAPARTIALCUTAREA  */
  YYSYMBOL_K_ANTENNAMAXCUTCAR = 332,       /* K_ANTENNAMAXCUTCAR  */
  YYSYMBOL_K_SLOTWIREWIDTH = 333,          /* K_SLOTWIREWIDTH  */
  YYSYMBOL_K_SLOTWIRELENGTH = 334,         /* K_SLOTWIRELENGTH  */
  YYSYMBOL_K_SLOTWIDTH = 335,              /* K_SLOTWIDTH  */
  YYSYMBOL_K_SLOTLENGTH = 336,             /* K_SLOTLENGTH  */
  YYSYMBOL_K_MAXADJACENTSLOTSPACING = 337, /* K_MAXADJACENTSLOTSPACING  */
  YYSYMBOL_K_MAXCOAXIALSLOTSPACING = 338,  /* K_MAXCOAXIALSLOTSPACING  */
  YYSYMBOL_K_MAXEDGESLOTSPACING = 339,     /* K_MAXEDGESLOTSPACING  */
  YYSYMBOL_K_SPLITWIREWIDTH = 340,         /* K_SPLITWIREWIDTH  */
  YYSYMBOL_K_MINIMUMDENSITY = 341,         /* K_MINIMUMDENSITY  */
  YYSYMBOL_K_MAXIMUMDENSITY = 342,         /* K_MAXIMUMDENSITY  */
  YYSYMBOL_K_DENSITYCHECKWINDOW = 343,     /* K_DENSITYCHECKWINDOW  */
  YYSYMBOL_K_DENSITYCHECKSTEP = 344,       /* K_DENSITYCHECKSTEP  */
  YYSYMBOL_K_FILLACTIVESPACING = 345,      /* K_FILLACTIVESPACING  */
  YYSYMBOL_K_MINIMUMCUT = 346,             /* K_MINIMUMCUT  */
  YYSYMBOL_K_ADJACENTCUTS = 347,           /* K_ADJACENTCUTS  */
  YYSYMBOL_K_ANTENNAMODEL = 348,           /* K_ANTENNAMODEL  */
  YYSYMBOL_K_BUMP = 349,                   /* K_BUMP  */
  YYSYMBOL_K_ENCLOSURE = 350,              /* K_ENCLOSURE  */
  YYSYMBOL_K_FROMABOVE = 351,              /* K_FROMABOVE  */
  YYSYMBOL_K_FROMBELOW = 352,              /* K_FROMBELOW  */
  YYSYMBOL_K_IMPLANT = 353,                /* K_IMPLANT  */
  YYSYMBOL_K_LENGTH = 354,                 /* K_LENGTH  */
  YYSYMBOL_K_MAXVIASTACK = 355,            /* K_MAXVIASTACK  */
  YYSYMBOL_K_AREAIO = 356,                 /* K_AREAIO  */
  YYSYMBOL_K_BLACKBOX = 357,               /* K_BLACKBOX  */
  YYSYMBOL_K_MAXWIDTH = 358,               /* K_MAXWIDTH  */
  YYSYMBOL_K_MINENCLOSEDAREA = 359,        /* K_MINENCLOSEDAREA  */
  YYSYMBOL_K_MINSTEP = 360,                /* K_MINSTEP  */
  YYSYMBOL_K_ORIENT = 361,                 /* K_ORIENT  */
  YYSYMBOL_K_OXIDE1 = 362,                 /* K_OXIDE1  */
  YYSYMBOL_K_OXIDE2 = 363,                 /* K_OXIDE2  */
  YYSYMBOL_K_OXIDE3 = 364,                 /* K_OXIDE3  */
  YYSYMBOL_K_OXIDE4 = 365,                 /* K_OXIDE4  */
  YYSYMBOL_K_PARALLELRUNLENGTH = 366,      /* K_PARALLELRUNLENGTH  */
  YYSYMBOL_K_MINWIDTH = 367,               /* K_MINWIDTH  */
  YYSYMBOL_K_PROTRUSIONWIDTH = 368,        /* K_PROTRUSIONWIDTH  */
  YYSYMBOL_K_SPACINGTABLE = 369,           /* K_SPACINGTABLE  */
  YYSYMBOL_K_WITHIN = 370,                 /* K_WITHIN  */
  YYSYMBOL_K_ABOVE = 371,                  /* K_ABOVE  */
  YYSYMBOL_K_BELOW = 372,                  /* K_BELOW  */
  YYSYMBOL_K_CENTERTOCENTER = 373,         /* K_CENTERTOCENTER  */
  YYSYMBOL_K_CUTSIZE = 374,                /* K_CUTSIZE  */
  YYSYMBOL_K_CUTSPACING = 375,             /* K_CUTSPACING  */
  YYSYMBOL_K_DENSITY = 376,                /* K_DENSITY  */
  YYSYMBOL_K_DIAG45 = 377,                 /* K_DIAG45  */
  YYSYMBOL_K_DIAG135 = 378,                /* K_DIAG135  */
  YYSYMBOL_K_MASK = 379,                   /* K_MASK  */
  YYSYMBOL_K_DIAGMINEDGELENGTH = 380,      /* K_DIAGMINEDGELENGTH  */
  YYSYMBOL_K_DIAGSPACING = 381,            /* K_DIAGSPACING  */
  YYSYMBOL_K_DIAGPITCH = 382,              /* K_DIAGPITCH  */
  YYSYMBOL_K_DIAGWIDTH = 383,              /* K_DIAGWIDTH  */
  YYSYMBOL_K_GENERATED = 384,              /* K_GENERATED  */
  YYSYMBOL_K_GROUNDSENSITIVITY = 385,      /* K_GROUNDSENSITIVITY  */
  YYSYMBOL_K_HARDSPACING = 386,            /* K_HARDSPACING  */
  YYSYMBOL_K_INSIDECORNER = 387,           /* K_INSIDECORNER  */
  YYSYMBOL_K_LAYERS = 388,                 /* K_LAYERS  */
  YYSYMBOL_K_LENGTHSUM = 389,              /* K_LENGTHSUM  */
  YYSYMBOL_K_MICRONS = 390,                /* K_MICRONS  */
  YYSYMBOL_K_MINCUTS = 391,                /* K_MINCUTS  */
  YYSYMBOL_K_MINSIZE = 392,                /* K_MINSIZE  */
  YYSYMBOL_K_NETEXPR = 393,                /* K_NETEXPR  */
  YYSYMBOL_K_OUTSIDECORNER = 394,          /* K_OUTSIDECORNER  */
  YYSYMBOL_K_PREFERENCLOSURE = 395,        /* K_PREFERENCLOSURE  */
  YYSYMBOL_K_ROWCOL = 396,                 /* K_ROWCOL  */
  YYSYMBOL_K_ROWPATTERN = 397,             /* K_ROWPATTERN  */
  YYSYMBOL_K_SOFT = 398,                   /* K_SOFT  */
  YYSYMBOL_K_SUPPLYSENSITIVITY = 399,      /* K_SUPPLYSENSITIVITY  */
  YYSYMBOL_K_USEVIA = 400,                 /* K_USEVIA  */
  YYSYMBOL_K_USEVIARULE = 401,             /* K_USEVIARULE  */
  YYSYMBOL_K_WELLTAP = 402,                /* K_WELLTAP  */
  YYSYMBOL_K_ARRAYCUTS = 403,              /* K_ARRAYCUTS  */
  YYSYMBOL_K_ARRAYSPACING = 404,           /* K_ARRAYSPACING  */
  YYSYMBOL_K_ANTENNAAREADIFFREDUCEPWL = 405, /* K_ANTENNAAREADIFFREDUCEPWL  */
  YYSYMBOL_K_ANTENNAAREAMINUSDIFF = 406,   /* K_ANTENNAAREAMINUSDIFF  */
  YYSYMBOL_K_ANTENNACUMROUTINGPLUSCUT = 407, /* K_ANTENNACUMROUTINGPLUSCUT  */
  YYSYMBOL_K_ANTENNAGATEPLUSDIFF = 408,    /* K_ANTENNAGATEPLUSDIFF  */
  YYSYMBOL_K_ENDOFLINE = 409,              /* K_ENDOFLINE  */
  YYSYMBOL_K_ENDOFNOTCHWIDTH = 410,        /* K_ENDOFNOTCHWIDTH  */
  YYSYMBOL_K_EXCEPTEXTRACUT = 411,         /* K_EXCEPTEXTRACUT  */
  YYSYMBOL_K_EXCEPTSAMEPGNET = 412,        /* K_EXCEPTSAMEPGNET  */
  YYSYMBOL_K_EXCEPTPGNET = 413,            /* K_EXCEPTPGNET  */
  YYSYMBOL_K_LONGARRAY = 414,              /* K_LONGARRAY  */
  YYSYMBOL_K_MAXEDGES = 415,               /* K_MAXEDGES  */
  YYSYMBOL_K_NOTCHLENGTH = 416,            /* K_NOTCHLENGTH  */
  YYSYMBOL_K_NOTCHSPACING = 417,           /* K_NOTCHSPACING  */
  YYSYMBOL_K_ORTHOGONAL = 418,             /* K_ORTHOGONAL  */
  YYSYMBOL_K_PARALLELEDGE = 419,           /* K_PARALLELEDGE  */
  YYSYMBOL_K_PARALLELOVERLAP = 420,        /* K_PARALLELOVERLAP  */
  YYSYMBOL_K_PGONLY = 421,                 /* K_PGONLY  */
  YYSYMBOL_K_PRL = 422,                    /* K_PRL  */
  YYSYMBOL_K_TWOEDGES = 423,               /* K_TWOEDGES  */
  YYSYMBOL_K_TWOWIDTHS = 424,              /* K_TWOWIDTHS  */
  YYSYMBOL_IF = 425,                       /* IF  */
  YYSYMBOL_LNOT = 426,                     /* LNOT  */
  YYSYMBOL_427_ = 427,                     /* '-'  */
  YYSYMBOL_428_ = 428,                     /* '+'  */
  YYSYMBOL_429_ = 429,                     /* '*'  */
  YYSYMBOL_430_ = 430,                     /* '/'  */
  YYSYMBOL_UMINUS = 431,                   /* UMINUS  */
  YYSYMBOL_432_ = 432,                     /* ';'  */
  YYSYMBOL_433_ = 433,                     /* '('  */
  YYSYMBOL_434_ = 434,                     /* ')'  */
  YYSYMBOL_435_ = 435,                     /* '='  */
  YYSYMBOL_436_n_ = 436,                   /* '\n'  */
  YYSYMBOL_437_ = 437,                     /* '<'  */
  YYSYMBOL_438_ = 438,                     /* '>'  */
  YYSYMBOL_YYACCEPT = 439,                 /* $accept  */
  YYSYMBOL_lef_file = 440,                 /* lef_file  */
  YYSYMBOL_version = 441,                  /* version  */
  YYSYMBOL_442_1 = 442,                    /* $@1  */
  YYSYMBOL_int_number = 443,               /* int_number  */
  YYSYMBOL_dividerchar = 444,              /* dividerchar  */
  YYSYMBOL_busbitchars = 445,              /* busbitchars  */
  YYSYMBOL_rules = 446,                    /* rules  */
  YYSYMBOL_end_library = 447,              /* end_library  */
  YYSYMBOL_rule = 448,                     /* rule  */
  YYSYMBOL_case_sensitivity = 449,         /* case_sensitivity  */
  YYSYMBOL_wireextension = 450,            /* wireextension  */
  YYSYMBOL_fixedmask = 451,                /* fixedmask  */
  YYSYMBOL_manufacturing = 452,            /* manufacturing  */
  YYSYMBOL_useminspacing = 453,            /* useminspacing  */
  YYSYMBOL_clearancemeasure = 454,         /* clearancemeasure  */
  YYSYMBOL_clearance_type = 455,           /* clearance_type  */
  YYSYMBOL_spacing_type = 456,             /* spacing_type  */
  YYSYMBOL_spacing_value = 457,            /* spacing_value  */
  YYSYMBOL_units_section = 458,            /* units_section  */
  YYSYMBOL_start_units = 459,              /* start_units  */
  YYSYMBOL_units_rules = 460,              /* units_rules  */
  YYSYMBOL_units_rule = 461,               /* units_rule  */
  YYSYMBOL_layer_rule = 462,               /* layer_rule  */
  YYSYMBOL_start_layer = 463,              /* start_layer  */
  YYSYMBOL_464_2 = 464,                    /* $@2  */
  YYSYMBOL_end_layer = 465,                /* end_layer  */
  YYSYMBOL_466_3 = 466,                    /* $@3  */
  YYSYMBOL_layer_options = 467,            /* layer_options  */
  YYSYMBOL_layer_option = 468,             /* layer_option  */
  YYSYMBOL_469_4 = 469,                    /* $@4  */
  YYSYMBOL_470_5 = 470,                    /* $@5  */
  YYSYMBOL_471_6 = 471,                    /* $@6  */
  YYSYMBOL_472_7 = 472,                    /* $@7  */
  YYSYMBOL_473_8 = 473,                    /* $@8  */
  YYSYMBOL_474_9 = 474,                    /* $@9  */
  YYSYMBOL_475_10 = 475,                   /* $@10  */
  YYSYMBOL_476_11 = 476,                   /* $@11  */
  YYSYMBOL_477_12 = 477,                   /* $@12  */
  YYSYMBOL_478_13 = 478,                   /* $@13  */
  YYSYMBOL_479_14 = 479,                   /* $@14  */
  YYSYMBOL_480_15 = 480,                   /* $@15  */
  YYSYMBOL_481_16 = 481,                   /* $@16  */
  YYSYMBOL_482_17 = 482,                   /* $@17  */
  YYSYMBOL_483_18 = 483,                   /* $@18  */
  YYSYMBOL_484_19 = 484,                   /* $@19  */
  YYSYMBOL_485_20 = 485,                   /* $@20  */
  YYSYMBOL_486_21 = 486,                   /* $@21  */
  YYSYMBOL_487_22 = 487,                   /* $@22  */
  YYSYMBOL_488_23 = 488,                   /* $@23  */
  YYSYMBOL_489_24 = 489,                   /* $@24  */
  YYSYMBOL_490_25 = 490,                   /* $@25  */
  YYSYMBOL_491_26 = 491,                   /* $@26  */
  YYSYMBOL_492_27 = 492,                   /* $@27  */
  YYSYMBOL_493_28 = 493,                   /* $@28  */
  YYSYMBOL_494_29 = 494,                   /* $@29  */
  YYSYMBOL_layer_arraySpacing_long = 495,  /* layer_arraySpacing_long  */
  YYSYMBOL_layer_arraySpacing_width = 496, /* layer_arraySpacing_width  */
  YYSYMBOL_layer_arraySpacing_arraycuts = 497, /* layer_arraySpacing_arraycuts  */
  YYSYMBOL_layer_arraySpacing_arraycut = 498, /* layer_arraySpacing_arraycut  */
  YYSYMBOL_sp_options = 499,               /* sp_options  */
  YYSYMBOL_500_30 = 500,                   /* $@30  */
  YYSYMBOL_501_31 = 501,                   /* $@31  */
  YYSYMBOL_502_32 = 502,                   /* $@32  */
  YYSYMBOL_503_33 = 503,                   /* $@33  */
  YYSYMBOL_504_34 = 504,                   /* $@34  */
  YYSYMBOL_505_35 = 505,                   /* $@35  */
  YYSYMBOL_506_36 = 506,                   /* $@36  */
  YYSYMBOL_layer_spacingtable_opts = 507,  /* layer_spacingtable_opts  */
  YYSYMBOL_layer_spacingtable_opt = 508,   /* layer_spacingtable_opt  */
  YYSYMBOL_layer_enclosure_type_opt = 509, /* layer_enclosure_type_opt  */
  YYSYMBOL_layer_enclosure_width_opt = 510, /* layer_enclosure_width_opt  */
  YYSYMBOL_511_37 = 511,                   /* $@37  */
  YYSYMBOL_layer_enclosure_width_except_opt = 512, /* layer_enclosure_width_except_opt  */
  YYSYMBOL_layer_preferenclosure_width_opt = 513, /* layer_preferenclosure_width_opt  */
  YYSYMBOL_layer_minimumcut_within = 514,  /* layer_minimumcut_within  */
  YYSYMBOL_layer_minimumcut_from = 515,    /* layer_minimumcut_from  */
  YYSYMBOL_layer_minimumcut_length = 516,  /* layer_minimumcut_length  */
  YYSYMBOL_layer_minstep_options = 517,    /* layer_minstep_options  */
  YYSYMBOL_layer_minstep_option = 518,     /* layer_minstep_option  */
  YYSYMBOL_layer_minstep_type = 519,       /* layer_minstep_type  */
  YYSYMBOL_layer_antenna_pwl = 520,        /* layer_antenna_pwl  */
  YYSYMBOL_521_38 = 521,                   /* $@38  */
  YYSYMBOL_layer_diffusion_ratios = 522,   /* layer_diffusion_ratios  */
  YYSYMBOL_layer_diffusion_ratio = 523,    /* layer_diffusion_ratio  */
  YYSYMBOL_layer_antenna_duo = 524,        /* layer_antenna_duo  */
  YYSYMBOL_layer_table_type = 525,         /* layer_table_type  */
  YYSYMBOL_layer_frequency = 526,          /* layer_frequency  */
  YYSYMBOL_527_39 = 527,                   /* $@39  */
  YYSYMBOL_528_40 = 528,                   /* $@40  */
  YYSYMBOL_529_41 = 529,                   /* $@41  */
  YYSYMBOL_ac_layer_table_opt = 530,       /* ac_layer_table_opt  */
  YYSYMBOL_531_42 = 531,                   /* $@42  */
  YYSYMBOL_532_43 = 532,                   /* $@43  */
  YYSYMBOL_dc_layer_table = 533,           /* dc_layer_table  */
  YYSYMBOL_534_44 = 534,                   /* $@44  */
  YYSYMBOL_int_number_list = 535,          /* int_number_list  */
  YYSYMBOL_number_list = 536,              /* number_list  */
  YYSYMBOL_layer_prop_list = 537,          /* layer_prop_list  */
  YYSYMBOL_layer_prop = 538,               /* layer_prop  */
  YYSYMBOL_current_density_pwl_list = 539, /* current_density_pwl_list  */
  YYSYMBOL_current_density_pwl = 540,      /* current_density_pwl  */
  YYSYMBOL_cap_points = 541,               /* cap_points  */
  YYSYMBOL_cap_point = 542,                /* cap_point  */
  YYSYMBOL_res_points = 543,               /* res_points  */
  YYSYMBOL_res_point = 544,                /* res_point  */
  YYSYMBOL_layer_type = 545,               /* layer_type  */
  YYSYMBOL_layer_direction = 546,          /* layer_direction  */
  YYSYMBOL_layer_minen_width = 547,        /* layer_minen_width  */
  YYSYMBOL_layer_oxide = 548,              /* layer_oxide  */
  YYSYMBOL_layer_sp_parallel_widths = 549, /* layer_sp_parallel_widths  */
  YYSYMBOL_layer_sp_parallel_width = 550,  /* layer_sp_parallel_width  */
  YYSYMBOL_551_45 = 551,                   /* $@45  */
  YYSYMBOL_layer_sp_TwoWidths = 552,       /* layer_sp_TwoWidths  */
  YYSYMBOL_layer_sp_TwoWidth = 553,        /* layer_sp_TwoWidth  */
  YYSYMBOL_554_46 = 554,                   /* $@46  */
  YYSYMBOL_layer_sp_TwoWidthsPRL = 555,    /* layer_sp_TwoWidthsPRL  */
  YYSYMBOL_layer_sp_influence_widths = 556, /* layer_sp_influence_widths  */
  YYSYMBOL_layer_sp_influence_width = 557, /* layer_sp_influence_width  */
  YYSYMBOL_maxstack_via = 558,             /* maxstack_via  */
  YYSYMBOL_559_47 = 559,                   /* $@47  */
  YYSYMBOL_via = 560,                      /* via  */
  YYSYMBOL_561_48 = 561,                   /* $@48  */
  YYSYMBOL_via_keyword = 562,              /* via_keyword  */
  YYSYMBOL_start_via = 563,                /* start_via  */
  YYSYMBOL_via_viarule = 564,              /* via_viarule  */
  YYSYMBOL_565_49 = 565,                   /* $@49  */
  YYSYMBOL_566_50 = 566,                   /* $@50  */
  YYSYMBOL_567_51 = 567,                   /* $@51  */
  YYSYMBOL_via_viarule_options = 568,      /* via_viarule_options  */
  YYSYMBOL_via_viarule_option = 569,       /* via_viarule_option  */
  YYSYMBOL_570_52 = 570,                   /* $@52  */
  YYSYMBOL_via_option = 571,               /* via_option  */
  YYSYMBOL_via_other_options = 572,        /* via_other_options  */
  YYSYMBOL_via_more_options = 573,         /* via_more_options  */
  YYSYMBOL_via_other_option = 574,         /* via_other_option  */
  YYSYMBOL_575_53 = 575,                   /* $@53  */
  YYSYMBOL_via_prop_list = 576,            /* via_prop_list  */
  YYSYMBOL_via_name_value_pair = 577,      /* via_name_value_pair  */
  YYSYMBOL_via_foreign = 578,              /* via_foreign  */
  YYSYMBOL_start_foreign = 579,            /* start_foreign  */
  YYSYMBOL_580_54 = 580,                   /* $@54  */
  YYSYMBOL_orientation = 581,              /* orientation  */
  YYSYMBOL_via_layer_rule = 582,           /* via_layer_rule  */
  YYSYMBOL_via_layer = 583,                /* via_layer  */
  YYSYMBOL_584_55 = 584,                   /* $@55  */
  YYSYMBOL_via_geometries = 585,           /* via_geometries  */
  YYSYMBOL_via_geometry = 586,             /* via_geometry  */
  YYSYMBOL_587_56 = 587,                   /* $@56  */
  YYSYMBOL_end_via = 588,                  /* end_via  */
  YYSYMBOL_589_57 = 589,                   /* $@57  */
  YYSYMBOL_viarule_keyword = 590,          /* viarule_keyword  */
  YYSYMBOL_591_58 = 591,                   /* $@58  */
  YYSYMBOL_viarule = 592,                  /* viarule  */
  YYSYMBOL_viarule_generate = 593,         /* viarule_generate  */
  YYSYMBOL_594_59 = 594,                   /* $@59  */
  YYSYMBOL_viarule_generate_default = 595, /* viarule_generate_default  */
  YYSYMBOL_viarule_layer_list = 596,       /* viarule_layer_list  */
  YYSYMBOL_opt_viarule_props = 597,        /* opt_viarule_props  */
  YYSYMBOL_viarule_props = 598,            /* viarule_props  */
  YYSYMBOL_viarule_prop = 599,             /* viarule_prop  */
  YYSYMBOL_600_60 = 600,                   /* $@60  */
  YYSYMBOL_viarule_prop_list = 601,        /* viarule_prop_list  */
  YYSYMBOL_viarule_layer = 602,            /* viarule_layer  */
  YYSYMBOL_via_names = 603,                /* via_names  */
  YYSYMBOL_via_name = 604,                 /* via_name  */
  YYSYMBOL_viarule_layer_name = 605,       /* viarule_layer_name  */
  YYSYMBOL_606_61 = 606,                   /* $@61  */
  YYSYMBOL_viarule_layer_options = 607,    /* viarule_layer_options  */
  YYSYMBOL_viarule_layer_option = 608,     /* viarule_layer_option  */
  YYSYMBOL_end_viarule = 609,              /* end_viarule  */
  YYSYMBOL_610_62 = 610,                   /* $@62  */
  YYSYMBOL_spacing_rule = 611,             /* spacing_rule  */
  YYSYMBOL_start_spacing = 612,            /* start_spacing  */
  YYSYMBOL_end_spacing = 613,              /* end_spacing  */
  YYSYMBOL_spacings = 614,                 /* spacings  */
  YYSYMBOL_spacing = 615,                  /* spacing  */
  YYSYMBOL_samenet_keyword = 616,          /* samenet_keyword  */
  YYSYMBOL_maskColor = 617,                /* maskColor  */
  YYSYMBOL_irdrop = 618,                   /* irdrop  */
  YYSYMBOL_start_irdrop = 619,             /* start_irdrop  */
  YYSYMBOL_end_irdrop = 620,               /* end_irdrop  */
  YYSYMBOL_ir_tables = 621,                /* ir_tables  */
  YYSYMBOL_ir_table = 622,                 /* ir_table  */
  YYSYMBOL_ir_table_values = 623,          /* ir_table_values  */
  YYSYMBOL_ir_table_value = 624,           /* ir_table_value  */
  YYSYMBOL_ir_tablename = 625,             /* ir_tablename  */
  YYSYMBOL_minfeature = 626,               /* minfeature  */
  YYSYMBOL_dielectric = 627,               /* dielectric  */
  YYSYMBOL_nondefault_rule = 628,          /* nondefault_rule  */
  YYSYMBOL_629_63 = 629,                   /* $@63  */
  YYSYMBOL_630_64 = 630,                   /* $@64  */
  YYSYMBOL_631_65 = 631,                   /* $@65  */
  YYSYMBOL_end_nd_rule = 632,              /* end_nd_rule  */
  YYSYMBOL_nd_hardspacing = 633,           /* nd_hardspacing  */
  YYSYMBOL_nd_rules = 634,                 /* nd_rules  */
  YYSYMBOL_nd_rule = 635,                  /* nd_rule  */
  YYSYMBOL_usevia = 636,                   /* usevia  */
  YYSYMBOL_useviarule = 637,               /* useviarule  */
  YYSYMBOL_mincuts = 638,                  /* mincuts  */
  YYSYMBOL_nd_prop = 639,                  /* nd_prop  */
  YYSYMBOL_640_66 = 640,                   /* $@66  */
  YYSYMBOL_nd_prop_list = 641,             /* nd_prop_list  */
  YYSYMBOL_nd_layer = 642,                 /* nd_layer  */
  YYSYMBOL_643_67 = 643,                   /* $@67  */
  YYSYMBOL_644_68 = 644,                   /* $@68  */
  YYSYMBOL_645_69 = 645,                   /* $@69  */
  YYSYMBOL_646_70 = 646,                   /* $@70  */
  YYSYMBOL_nd_layer_stmts = 647,           /* nd_layer_stmts  */
  YYSYMBOL_nd_layer_stmt = 648,            /* nd_layer_stmt  */
  YYSYMBOL_site = 649,                     /* site  */
  YYSYMBOL_start_site = 650,               /* start_site  */
  YYSYMBOL_651_71 = 651,                   /* $@71  */
  YYSYMBOL_end_site = 652,                 /* end_site  */
  YYSYMBOL_653_72 = 653,                   /* $@72  */
  YYSYMBOL_site_options = 654,             /* site_options  */
  YYSYMBOL_site_option = 655,              /* site_option  */
  YYSYMBOL_site_class = 656,               /* site_class  */
  YYSYMBOL_site_symmetry_statement = 657,  /* site_symmetry_statement  */
  YYSYMBOL_site_symmetries = 658,          /* site_symmetries  */
  YYSYMBOL_site_symmetry = 659,            /* site_symmetry  */
  YYSYMBOL_site_rowpattern_statement = 660, /* site_rowpattern_statement  */
  YYSYMBOL_661_73 = 661,                   /* $@73  */
  YYSYMBOL_site_rowpatterns = 662,         /* site_rowpatterns  */
  YYSYMBOL_site_rowpattern = 663,          /* site_rowpattern  */
  YYSYMBOL_664_74 = 664,                   /* $@74  */
  YYSYMBOL_pt = 665,                       /* pt  */
  YYSYMBOL_macro = 666,                    /* macro  */
  YYSYMBOL_667_75 = 667,                   /* $@75  */
  YYSYMBOL_start_macro = 668,              /* start_macro  */
  YYSYMBOL_669_76 = 669,                   /* $@76  */
  YYSYMBOL_end_macro = 670,                /* end_macro  */
  YYSYMBOL_671_77 = 671,                   /* $@77  */
  YYSYMBOL_macro_options = 672,            /* macro_options  */
  YYSYMBOL_macro_option = 673,             /* macro_option  */
  YYSYMBOL_674_78 = 674,                   /* $@78  */
  YYSYMBOL_macro_prop_list = 675,          /* macro_prop_list  */
  YYSYMBOL_macro_symmetry_statement = 676, /* macro_symmetry_statement  */
  YYSYMBOL_macro_symmetries = 677,         /* macro_symmetries  */
  YYSYMBOL_macro_symmetry = 678,           /* macro_symmetry  */
  YYSYMBOL_macro_name_value_pair = 679,    /* macro_name_value_pair  */
  YYSYMBOL_macro_class = 680,              /* macro_class  */
  YYSYMBOL_class_type = 681,               /* class_type  */
  YYSYMBOL_pad_type = 682,                 /* pad_type  */
  YYSYMBOL_core_type = 683,                /* core_type  */
  YYSYMBOL_endcap_type = 684,              /* endcap_type  */
  YYSYMBOL_macro_generator = 685,          /* macro_generator  */
  YYSYMBOL_macro_generate = 686,           /* macro_generate  */
  YYSYMBOL_macro_source = 687,             /* macro_source  */
  YYSYMBOL_macro_power = 688,              /* macro_power  */
  YYSYMBOL_macro_origin = 689,             /* macro_origin  */
  YYSYMBOL_macro_foreign = 690,            /* macro_foreign  */
  YYSYMBOL_macro_fixedMask = 691,          /* macro_fixedMask  */
  YYSYMBOL_macro_eeq = 692,                /* macro_eeq  */
  YYSYMBOL_693_79 = 693,                   /* $@79  */
  YYSYMBOL_macro_leq = 694,                /* macro_leq  */
  YYSYMBOL_695_80 = 695,                   /* $@80  */
  YYSYMBOL_macro_site = 696,               /* macro_site  */
  YYSYMBOL_macro_site_word = 697,          /* macro_site_word  */
  YYSYMBOL_site_word = 698,                /* site_word  */
  YYSYMBOL_macro_size = 699,               /* macro_size  */
  YYSYMBOL_macro_pin = 700,                /* macro_pin  */
  YYSYMBOL_start_macro_pin = 701,          /* start_macro_pin  */
  YYSYMBOL_702_81 = 702,                   /* $@81  */
  YYSYMBOL_end_macro_pin = 703,            /* end_macro_pin  */
  YYSYMBOL_704_82 = 704,                   /* $@82  */
  YYSYMBOL_macro_pin_options = 705,        /* macro_pin_options  */
  YYSYMBOL_macro_pin_option = 706,         /* macro_pin_option  */
  YYSYMBOL_707_83 = 707,                   /* $@83  */
  YYSYMBOL_708_84 = 708,                   /* $@84  */
  YYSYMBOL_709_85 = 709,                   /* $@85  */
  YYSYMBOL_710_86 = 710,                   /* $@86  */
  YYSYMBOL_711_87 = 711,                   /* $@87  */
  YYSYMBOL_712_88 = 712,                   /* $@88  */
  YYSYMBOL_713_89 = 713,                   /* $@89  */
  YYSYMBOL_714_90 = 714,                   /* $@90  */
  YYSYMBOL_715_91 = 715,                   /* $@91  */
  YYSYMBOL_716_92 = 716,                   /* $@92  */
  YYSYMBOL_pin_layer_oxide = 717,          /* pin_layer_oxide  */
  YYSYMBOL_pin_prop_list = 718,            /* pin_prop_list  */
  YYSYMBOL_pin_name_value_pair = 719,      /* pin_name_value_pair  */
  YYSYMBOL_electrical_direction = 720,     /* electrical_direction  */
  YYSYMBOL_start_macro_port = 721,         /* start_macro_port  */
  YYSYMBOL_macro_port_class_option = 722,  /* macro_port_class_option  */
  YYSYMBOL_macro_pin_use = 723,            /* macro_pin_use  */
  YYSYMBOL_macro_scan_use = 724,           /* macro_scan_use  */
  YYSYMBOL_pin_shape = 725,                /* pin_shape  */
  YYSYMBOL_geometries = 726,               /* geometries  */
  YYSYMBOL_geometry = 727,                 /* geometry  */
  YYSYMBOL_728_93 = 728,                   /* $@93  */
  YYSYMBOL_729_94 = 729,                   /* $@94  */
  YYSYMBOL_geometry_options = 730,         /* geometry_options  */
  YYSYMBOL_layer_exceptpgnet = 731,        /* layer_exceptpgnet  */
  YYSYMBOL_layer_spacing = 732,            /* layer_spacing  */
  YYSYMBOL_firstPt = 733,                  /* firstPt  */
  YYSYMBOL_nextPt = 734,                   /* nextPt  */
  YYSYMBOL_otherPts = 735,                 /* otherPts  */
  YYSYMBOL_via_placement = 736,            /* via_placement  */
  YYSYMBOL_737_95 = 737,                   /* $@95  */
  YYSYMBOL_738_96 = 738,                   /* $@96  */
  YYSYMBOL_stepPattern = 739,              /* stepPattern  */
  YYSYMBOL_sitePattern = 740,              /* sitePattern  */
  YYSYMBOL_trackPattern = 741,             /* trackPattern  */
  YYSYMBOL_742_97 = 742,                   /* $@97  */
  YYSYMBOL_743_98 = 743,                   /* $@98  */
  YYSYMBOL_744_99 = 744,                   /* $@99  */
  YYSYMBOL_745_100 = 745,                  /* $@100  */
  YYSYMBOL_trackLayers = 746,              /* trackLayers  */
  YYSYMBOL_layer_name = 747,               /* layer_name  */
  YYSYMBOL_gcellPattern = 748,             /* gcellPattern  */
  YYSYMBOL_macro_obs = 749,                /* macro_obs  */
  YYSYMBOL_start_macro_obs = 750,          /* start_macro_obs  */
  YYSYMBOL_macro_density = 751,            /* macro_density  */
  YYSYMBOL_density_layers = 752,           /* density_layers  */
  YYSYMBOL_density_layer = 753,            /* density_layer  */
  YYSYMBOL_754_101 = 754,                  /* $@101  */
  YYSYMBOL_755_102 = 755,                  /* $@102  */
  YYSYMBOL_density_layer_rects = 756,      /* density_layer_rects  */
  YYSYMBOL_density_layer_rect = 757,       /* density_layer_rect  */
  YYSYMBOL_macro_clocktype = 758,          /* macro_clocktype  */
  YYSYMBOL_759_103 = 759,                  /* $@103  */
  YYSYMBOL_timing = 760,                   /* timing  */
  YYSYMBOL_start_timing = 761,             /* start_timing  */
  YYSYMBOL_end_timing = 762,               /* end_timing  */
  YYSYMBOL_timing_options = 763,           /* timing_options  */
  YYSYMBOL_timing_option = 764,            /* timing_option  */
  YYSYMBOL_765_104 = 765,                  /* $@104  */
  YYSYMBOL_766_105 = 766,                  /* $@105  */
  YYSYMBOL_767_106 = 767,                  /* $@106  */
  YYSYMBOL_one_pin_trigger = 768,          /* one_pin_trigger  */
  YYSYMBOL_two_pin_trigger = 769,          /* two_pin_trigger  */
  YYSYMBOL_from_pin_trigger = 770,         /* from_pin_trigger  */
  YYSYMBOL_to_pin_trigger = 771,           /* to_pin_trigger  */
  YYSYMBOL_delay_or_transition = 772,      /* delay_or_transition  */
  YYSYMBOL_list_of_table_entries = 773,    /* list_of_table_entries  */
  YYSYMBOL_table_entry = 774,              /* table_entry  */
  YYSYMBOL_list_of_table_axis_dnumbers = 775, /* list_of_table_axis_dnumbers  */
  YYSYMBOL_slew_spec = 776,                /* slew_spec  */
  YYSYMBOL_risefall = 777,                 /* risefall  */
  YYSYMBOL_unateness = 778,                /* unateness  */
  YYSYMBOL_list_of_from_strings = 779,     /* list_of_from_strings  */
  YYSYMBOL_list_of_to_strings = 780,       /* list_of_to_strings  */
  YYSYMBOL_array = 781,                    /* array  */
  YYSYMBOL_782_107 = 782,                  /* $@107  */
  YYSYMBOL_start_array = 783,              /* start_array  */
  YYSYMBOL_784_108 = 784,                  /* $@108  */
  YYSYMBOL_end_array = 785,                /* end_array  */
  YYSYMBOL_786_109 = 786,                  /* $@109  */
  YYSYMBOL_array_rules = 787,              /* array_rules  */
  YYSYMBOL_array_rule = 788,               /* array_rule  */
  YYSYMBOL_789_110 = 789,                  /* $@110  */
  YYSYMBOL_790_111 = 790,                  /* $@111  */
  YYSYMBOL_791_112 = 791,                  /* $@112  */
  YYSYMBOL_792_113 = 792,                  /* $@113  */
  YYSYMBOL_793_114 = 793,                  /* $@114  */
  YYSYMBOL_floorplan_start = 794,          /* floorplan_start  */
  YYSYMBOL_floorplan_list = 795,           /* floorplan_list  */
  YYSYMBOL_floorplan_element = 796,        /* floorplan_element  */
  YYSYMBOL_797_115 = 797,                  /* $@115  */
  YYSYMBOL_798_116 = 798,                  /* $@116  */
  YYSYMBOL_cap_list = 799,                 /* cap_list  */
  YYSYMBOL_one_cap = 800,                  /* one_cap  */
  YYSYMBOL_msg_statement = 801,            /* msg_statement  */
  YYSYMBOL_802_117 = 802,                  /* $@117  */
  YYSYMBOL_create_file_statement = 803,    /* create_file_statement  */
  YYSYMBOL_804_118 = 804,                  /* $@118  */
  YYSYMBOL_def_statement = 805,            /* def_statement  */
  YYSYMBOL_806_119 = 806,                  /* $@119  */
  YYSYMBOL_807_120 = 807,                  /* $@120  */
  YYSYMBOL_808_121 = 808,                  /* $@121  */
  YYSYMBOL_dtrm = 809,                     /* dtrm  */
  YYSYMBOL_then = 810,                     /* then  */
  YYSYMBOL_else = 811,                     /* else  */
  YYSYMBOL_expression = 812,               /* expression  */
  YYSYMBOL_b_expr = 813,                   /* b_expr  */
  YYSYMBOL_s_expr = 814,                   /* s_expr  */
  YYSYMBOL_relop = 815,                    /* relop  */
  YYSYMBOL_prop_def_section = 816,         /* prop_def_section  */
  YYSYMBOL_817_122 = 817,                  /* $@122  */
  YYSYMBOL_prop_stmts = 818,               /* prop_stmts  */
  YYSYMBOL_prop_stmt = 819,                /* prop_stmt  */
  YYSYMBOL_820_123 = 820,                  /* $@123  */
  YYSYMBOL_821_124 = 821,                  /* $@124  */
  YYSYMBOL_822_125 = 822,                  /* $@125  */
  YYSYMBOL_823_126 = 823,                  /* $@126  */
  YYSYMBOL_824_127 = 824,                  /* $@127  */
  YYSYMBOL_825_128 = 825,                  /* $@128  */
  YYSYMBOL_826_129 = 826,                  /* $@129  */
  YYSYMBOL_827_130 = 827,                  /* $@130  */
  YYSYMBOL_prop_define = 828,              /* prop_define  */
  YYSYMBOL_opt_range_second = 829,         /* opt_range_second  */
  YYSYMBOL_opt_endofline = 830,            /* opt_endofline  */
  YYSYMBOL_831_131 = 831,                  /* $@131  */
  YYSYMBOL_opt_endofline_twoedges = 832,   /* opt_endofline_twoedges  */
  YYSYMBOL_opt_samenetPGonly = 833,        /* opt_samenetPGonly  */
  YYSYMBOL_opt_def_range = 834,            /* opt_def_range  */
  YYSYMBOL_opt_def_value = 835,            /* opt_def_value  */
  YYSYMBOL_opt_def_dvalue = 836,           /* opt_def_dvalue  */
  YYSYMBOL_layer_spacing_opts = 837,       /* layer_spacing_opts  */
  YYSYMBOL_layer_spacing_opt = 838,        /* layer_spacing_opt  */
  YYSYMBOL_839_132 = 839,                  /* $@132  */
  YYSYMBOL_layer_spacing_cut_routing = 840, /* layer_spacing_cut_routing  */
  YYSYMBOL_841_133 = 841,                  /* $@133  */
  YYSYMBOL_842_134 = 842,                  /* $@134  */
  YYSYMBOL_843_135 = 843,                  /* $@135  */
  YYSYMBOL_844_136 = 844,                  /* $@136  */
  YYSYMBOL_845_137 = 845,                  /* $@137  */
  YYSYMBOL_spacing_cut_layer_opt = 846,    /* spacing_cut_layer_opt  */
  YYSYMBOL_opt_adjacentcuts_exceptsame = 847, /* opt_adjacentcuts_exceptsame  */
  YYSYMBOL_opt_layer_name = 848,           /* opt_layer_name  */
  YYSYMBOL_849_138 = 849,                  /* $@138  */
  YYSYMBOL_req_layer_name = 850,           /* req_layer_name  */
  YYSYMBOL_851_139 = 851,                  /* $@139  */
  YYSYMBOL_universalnoisemargin = 852,     /* universalnoisemargin  */
  YYSYMBOL_edgeratethreshold1 = 853,       /* edgeratethreshold1  */
  YYSYMBOL_edgeratethreshold2 = 854,       /* edgeratethreshold2  */
  YYSYMBOL_edgeratescalefactor = 855,      /* edgeratescalefactor  */
  YYSYMBOL_noisetable = 856,               /* noisetable  */
  YYSYMBOL_857_140 = 857,                  /* $@140  */
  YYSYMBOL_end_noisetable = 858,           /* end_noisetable  */
  YYSYMBOL_noise_table_list = 859,         /* noise_table_list  */
  YYSYMBOL_noise_table_entry = 860,        /* noise_table_entry  */
  YYSYMBOL_output_resistance_entry = 861,  /* output_resistance_entry  */
  YYSYMBOL_862_141 = 862,                  /* $@141  */
  YYSYMBOL_num_list = 863,                 /* num_list  */
  YYSYMBOL_victim_list = 864,              /* victim_list  */
  YYSYMBOL_victim = 865,                   /* victim  */
  YYSYMBOL_866_142 = 866,                  /* $@142  */
  YYSYMBOL_vnoiselist = 867,               /* vnoiselist  */
  YYSYMBOL_correctiontable = 868,          /* correctiontable  */
  YYSYMBOL_869_143 = 869,                  /* $@143  */
  YYSYMBOL_end_correctiontable = 870,      /* end_correctiontable  */
  YYSYMBOL_correction_table_list = 871,    /* correction_table_list  */
  YYSYMBOL_correction_table_item = 872,    /* correction_table_item  */
  YYSYMBOL_output_list = 873,              /* output_list  */
  YYSYMBOL_874_144 = 874,                  /* $@144  */
  YYSYMBOL_numo_list = 875,                /* numo_list  */
  YYSYMBOL_corr_victim_list = 876,         /* corr_victim_list  */
  YYSYMBOL_corr_victim = 877,              /* corr_victim  */
  YYSYMBOL_878_145 = 878,                  /* $@145  */
  YYSYMBOL_corr_list = 879,                /* corr_list  */
  YYSYMBOL_input_antenna = 880,            /* input_antenna  */
  YYSYMBOL_output_antenna = 881,           /* output_antenna  */
  YYSYMBOL_inout_antenna = 882,            /* inout_antenna  */
  YYSYMBOL_antenna_input = 883,            /* antenna_input  */
  YYSYMBOL_antenna_inout = 884,            /* antenna_inout  */
  YYSYMBOL_antenna_output = 885,           /* antenna_output  */
  YYSYMBOL_extension_opt = 886,            /* extension_opt  */
  YYSYMBOL_extension = 887                 /* extension  */
};
typedef enum yysymbol_kind_t yysymbol_kind_t;




#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

/* Work around bug in HP-UX 11.23, which defines these macros
   incorrectly for preprocessor constants.  This workaround can likely
   be removed in 2023, as HPE has promised support for HP-UX 11.23
   (aka HP-UX 11i v2) only through the end of 2022; see Table 2 of
   <https://h20195.www2.hpe.com/V2/getpdf.aspx/4AA4-7673ENW.pdf>.  */
#ifdef __hpux
# undef UINT_LEAST8_MAX
# undef UINT_LEAST16_MAX
# define UINT_LEAST8_MAX 255
# define UINT_LEAST16_MAX 65535
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))


/* Stored state numbers (used for stacks). */
typedef yytype_int16 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if !defined yyoverflow

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* !defined yyoverflow */

#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  4
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   2695

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  439
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  449
/* YYNRULES -- Number of rules.  */
#define YYNRULES  1020
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  2044

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   682


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK                     \
   ? YY_CAST (yysymbol_kind_t, yytranslate[YYX])        \
   : YYSYMBOL_YYUNDEF)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int16 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     436,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     433,   434,   429,   428,     2,   427,     2,   430,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,   432,
     437,   435,   438,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   171,   172,   173,   174,
     175,   176,   177,   178,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   188,   189,   190,   191,   192,   193,   194,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,   239,   240,   241,   242,   243,   244,
     245,   246,   247,   248,   249,   250,   251,   252,   253,   254,
     255,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,   423,   424,
     425,   426,   431
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   334,   334,   355,   355,   394,   411,   424,   437,   438,
     439,   443,   449,   459,   459,   459,   459,   460,   460,   460,
     460,   460,   461,   461,   462,   462,   462,   462,   462,   462,
     462,   463,   463,   463,   464,   464,   465,   465,   466,   466,
     466,   467,   467,   468,   468,   468,   468,   468,   469,   469,
     469,   470,   473,   487,   505,   515,   526,   539,   546,   561,
     565,   566,   569,   570,   573,   574,   576,   582,   605,   606,
     609,   611,   613,   615,   617,   619,   621,   628,   631,   638,
     638,   670,   670,   728,   729,   734,   740,   733,   757,   763,
     778,   783,   788,   792,   796,   800,   804,   808,   812,   817,
     834,   833,   870,   869,   887,   901,   913,   924,   936,   947,
     959,   971,   983,   995,  1007,  1020,  1047,  1066,  1084,  1103,
    1103,  1108,  1107,  1122,  1137,  1153,  1176,  1152,  1179,  1202,
    1178,  1206,  1245,  1244,  1283,  1322,  1321,  1361,  1360,  1375,
    1414,  1413,  1452,  1491,  1490,  1530,  1529,  1570,  1569,  1608,
    1629,  1650,  1672,  1693,  1671,  1709,  1733,  1757,  1781,  1805,
    1829,  1853,  1877,  1900,  1918,  1936,  1954,  1972,  1990,  2015,
    2041,  2040,  2058,  2057,  2073,  2072,  2079,  2096,  2095,  2122,
    2121,  2141,  2140,  2158,  2176,  2202,  2201,  2219,  2221,  2227,
    2229,  2235,  2237,  2240,  2257,  2278,  2283,  2289,  2256,  2303,
    2307,  2302,  2339,  2338,  2363,  2365,  2368,  2375,  2376,  2377,
    2379,  2381,  2380,  2387,  2403,  2404,  2420,  2421,  2428,  2429,
    2445,  2446,  2465,  2484,  2485,  2503,  2504,  2507,  2511,  2515,
    2529,  2530,  2531,  2534,  2538,  2537,  2557,  2558,  2562,  2567,
    2568,  2602,  2603,  2604,  2608,  2610,  2613,  2607,  2617,  2619,
    2618,  2633,  2632,  2649,  2648,  2653,  2654,  2657,  2658,  2662,
    2663,  2667,  2675,  2683,  2695,  2697,  2700,  2704,  2705,  2708,
    2712,  2713,  2716,  2720,  2721,  2722,  2723,  2724,  2725,  2728,
    2729,  2730,  2731,  2733,  2734,  2741,  2746,  2751,  2756,  2763,
    2764,  2768,  2767,  2777,  2778,  2782,  2781,  2792,  2796,  2803,
    2804,  2807,  2810,  2846,  2846,  2873,  2873,  2883,  2886,  2895,
    2903,  2912,  2914,  2917,  2912,  2938,  2939,  2942,  2946,  2950,
    2954,  2954,  2960,  2961,  2963,  2966,  2967,  2971,  2973,  2975,
    2977,  2977,  2980,  2991,  2992,  2996,  3006,  3014,  3024,  3033,
    3042,  3051,  3061,  3061,  3065,  3066,  3067,  3068,  3069,  3070,
    3071,  3072,  3073,  3074,  3075,  3076,  3077,  3078,  3079,  3080,
    3082,  3085,  3085,  3092,  3094,  3098,  3112,  3111,  3135,  3135,
    3169,  3169,  3179,  3198,  3197,  3225,  3226,  3244,  3245,  3248,
    3250,  3254,  3255,  3258,  3258,  3263,  3264,  3268,  3276,  3284,
    3295,  3326,  3328,  3331,  3334,  3334,  3340,  3342,  3346,  3365,
    3384,  3420,  3422,  3425,  3427,  3429,  3459,  3487,  3487,  3521,
    3524,  3546,  3559,  3561,  3564,  3578,  3593,  3599,  3600,  3603,
    3606,  3617,  3626,  3628,  3631,  3639,  3641,  3644,  3647,  3650,
    3664,  3675,  3676,  3686,  3675,  3716,  3721,  3742,  3744,  3763,
    3764,  3768,  3769,  3770,  3771,  3772,  3773,  3774,  3777,  3794,
    3813,  3832,  3832,  3837,  3838,  3842,  3850,  3858,  3869,  3870,
    3879,  3883,  3869,  3925,  3927,  3931,  3936,  3939,  3963,  3986,
    4009,  4028,  4034,  4034,  4044,  4044,  4075,  4077,  4081,  4087,
    4089,  4094,  4098,  4099,  4100,  4102,  4105,  4107,  4111,  4113,
    4115,  4118,  4118,  4122,  4124,  4127,  4127,  4131,  4133,  4137,
    4136,  4144,  4144,  4162,  4162,  4184,  4186,  4190,  4191,  4192,
    4193,  4194,  4195,  4197,  4199,  4201,  4203,  4204,  4205,  4207,
    4209,  4211,  4213,  4215,  4217,  4219,  4221,  4223,  4223,  4228,
    4229,  4232,  4243,  4245,  4249,  4251,  4253,  4257,  4267,  4275,
    4284,  4292,  4293,  4312,  4313,  4314,  4333,  4352,  4353,  4366,
    4367,  4368,  4391,  4392,  4398,  4401,  4406,  4407,  4408,  4409,
    4410,  4411,  4414,  4415,  4416,  4417,  4434,  4451,  4470,  4471,
    4472,  4473,  4474,  4475,  4477,  4480,  4484,  4493,  4502,  4512,
    4522,  4562,  4566,  4570,  4574,  4580,  4590,  4590,  4593,  4593,
    4604,  4610,  4620,  4624,  4627,  4643,  4650,  4650,  4656,  4656,
    4678,  4679,  4683,  4692,  4701,  4710,  4719,  4728,  4737,  4737,
    4746,  4755,  4757,  4759,  4761,  4770,  4779,  4788,  4797,  4806,
    4815,  4824,  4833,  4835,  4835,  4837,  4837,  4846,  4846,  4855,
    4855,  4864,  4873,  4875,  4877,  4886,  4895,  4904,  4913,  4922,
    4931,  4940,  4942,  4942,  4945,  4959,  4971,  4992,  5013,  5034,
    5036,  5038,  5068,  5098,  5128,  5158,  5188,  5218,  5248,  5279,
    5278,  5308,  5308,  5324,  5324,  5340,  5340,  5358,  5363,  5368,
    5373,  5380,  5381,  5385,  5395,  5403,  5413,  5414,  5415,  5416,
    5417,  5419,  5431,  5432,  5437,  5438,  5439,  5440,  5441,  5442,
    5445,  5446,  5447,  5448,  5451,  5452,  5453,  5454,  5456,  5459,
    5460,  5459,  5476,  5488,  5509,  5530,  5550,  5570,  5592,  5613,
    5616,  5617,  5619,  5620,  5635,  5636,  5650,  5665,  5669,  5673,
    5675,  5679,  5679,  5692,  5692,  5707,  5711,  5722,  5735,  5743,
    5734,  5746,  5754,  5745,  5756,  5765,  5775,  5777,  5780,  5783,
    5792,  5802,  5813,  5825,  5837,  5858,  5859,  5862,  5863,  5862,
    5869,  5870,  5873,  5879,  5879,  5882,  5885,  5888,  5902,  5904,
    5909,  5908,  5919,  5919,  5922,  5921,  5925,  5934,  5936,  5938,
    5940,  5942,  5944,  5946,  5948,  5950,  5952,  5954,  5956,  5958,
    5960,  5962,  5964,  5966,  5968,  5972,  5974,  5976,  5980,  5982,
    5984,  5986,  5990,  5992,  5994,  5998,  6000,  6002,  6006,  6008,
    6012,  6014,  6017,  6021,  6023,  6028,  6029,  6031,  6036,  6038,
    6042,  6044,  6046,  6050,  6052,  6056,  6058,  6062,  6061,  6071,
    6071,  6081,  6081,  6105,  6106,  6110,  6110,  6117,  6117,  6124,
    6124,  6131,  6131,  6137,  6140,  6140,  6146,  6152,  6155,  6160,
    6161,  6165,  6165,  6172,  6172,  6182,  6183,  6186,  6190,  6190,
    6194,  6194,  6198,  6198,  6207,  6207,  6216,  6216,  6227,  6228,
    6229,  6232,  6233,  6237,  6238,  6242,  6243,  6244,  6245,  6246,
    6247,  6248,  6250,  6253,  6254,  6255,  6256,  6257,  6258,  6259,
    6260,  6261,  6262,  6263,  6264,  6265,  6267,  6268,  6271,  6277,
    6279,  6288,  6292,  6293,  6294,  6295,  6296,  6297,  6298,  6299,
    6300,  6304,  6303,  6316,  6317,  6321,  6321,  6330,  6330,  6339,
    6339,  6349,  6349,  6358,  6358,  6367,  6367,  6376,  6376,  6385,
    6385,  6396,  6401,  6406,  6411,  6416,  6424,  6425,  6430,  6437,
    6444,  6452,  6454,  6453,  6462,  6463,  6471,  6472,  6480,  6481,
    6486,  6487,  6492,  6493,  6496,  6498,  6500,  6525,  6524,  6549,
    6573,  6575,  6576,  6575,  6592,  6591,  6616,  6640,  6639,  6645,
    6651,  6659,  6658,  6674,  6688,  6705,  6706,  6714,  6715,  6732,
    6733,  6733,  6738,  6738,  6742,  6756,  6769,  6782,  6796,  6795,
    6801,  6814,  6815,  6819,  6826,  6830,  6829,  6835,  6838,  6843,
    6844,  6848,  6847,  6854,  6857,  6862,  6861,  6868,  6881,  6882,
    6886,  6893,  6897,  6896,  6903,  6906,  6911,  6912,  6917,  6916,
    6923,  6926,  6932,  6956,  6980,  7004,  7038,  7072,  7106,  7107,
    7109
};
#endif

/** Accessing symbol of state STATE.  */
#define YY_ACCESSING_SYMBOL(State) YY_CAST (yysymbol_kind_t, yystos[State])

#if YYDEBUG || 0
/* The user-facing name of the symbol whose (internal) number is
   YYSYMBOL.  No bounds checking.  */
static const char *yysymbol_name (yysymbol_kind_t yysymbol) YY_ATTRIBUTE_UNUSED;

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "\"invalid token\"", "K_HISTORY", "K_ABUT",
  "K_ABUTMENT", "K_ACTIVE", "K_ANALOG", "K_ARRAY", "K_AREA", "K_BLOCK",
  "K_BOTTOMLEFT", "K_BOTTOMRIGHT", "K_BY", "K_CAPACITANCE",
  "K_CAPMULTIPLIER", "K_CLASS", "K_CLOCK", "K_CLOCKTYPE", "K_COLUMNMAJOR",
  "K_DESIGNRULEWIDTH", "K_INFLUENCE", "K_CORE", "K_CORNER", "K_COVER",
  "K_CPERSQDIST", "K_CURRENT", "K_CURRENTSOURCE", "K_CUT", "K_DEFAULT",
  "K_DATABASE", "K_DATA", "K_DIELECTRIC", "K_DIRECTION", "K_DO",
  "K_EDGECAPACITANCE", "K_EEQ", "K_END", "K_ENDCAP", "K_FALL", "K_FALLCS",
  "K_FALLT0", "K_FALLSATT1", "K_FALLRS", "K_FALLSATCUR", "K_FALLTHRESH",
  "K_FEEDTHRU", "K_FIXED", "K_FOREIGN", "K_FROMPIN", "K_GENERATE",
  "K_GENERATOR", "K_GROUND", "K_HEIGHT", "K_HORIZONTAL", "K_INOUT",
  "K_INPUT", "K_INPUTNOISEMARGIN", "K_COMPONENTPIN", "K_INTRINSIC",
  "K_INVERT", "K_IRDROP", "K_ITERATE", "K_IV_TABLES", "K_LAYER",
  "K_LEAKAGE", "K_LEQ", "K_LIBRARY", "K_MACRO", "K_MATCH", "K_MAXDELAY",
  "K_MAXLOAD", "K_METALOVERHANG", "K_MILLIAMPS", "K_MILLIWATTS",
  "K_MINFEATURE", "K_MUSTJOIN", "K_NAMESCASESENSITIVE", "K_NANOSECONDS",
  "K_NETS", "K_NEW", "K_NONDEFAULTRULE", "K_NONINVERT", "K_NONUNATE",
  "K_OBS", "K_OHMS", "K_OFFSET", "K_ORIENTATION", "K_ORIGIN", "K_OUTPUT",
  "K_OUTPUTNOISEMARGIN", "K_OVERHANG", "K_OVERLAP", "K_OFF", "K_ON",
  "K_OVERLAPS", "K_PAD", "K_PATH", "K_PATTERN", "K_PICOFARADS", "K_PIN",
  "K_PITCH", "K_PLACED", "K_POLYGON", "K_PORT", "K_POST", "K_POWER",
  "K_PRE", "K_PULLDOWNRES", "K_RECT", "K_RESISTANCE", "K_RESISTIVE",
  "K_RING", "K_RISE", "K_RISECS", "K_RISERS", "K_RISESATCUR",
  "K_RISETHRESH", "K_RISESATT1", "K_RISET0", "K_RISEVOLTAGETHRESHOLD",
  "K_FALLVOLTAGETHRESHOLD", "K_ROUTING", "K_ROWMAJOR", "K_RPERSQ",
  "K_SAMENET", "K_SCANUSE", "K_SHAPE", "K_SHRINKAGE", "K_SIGNAL", "K_SITE",
  "K_SIZE", "K_SOURCE", "K_SPACER", "K_SPACING", "K_SPECIALNETS",
  "K_STACK", "K_START", "K_STEP", "K_STOP", "K_STRUCTURE", "K_SYMMETRY",
  "K_TABLE", "K_THICKNESS", "K_TIEHIGH", "K_TIELOW", "K_TIEOFFR", "K_TIME",
  "K_TIMING", "K_TO", "K_TOPIN", "K_TOPLEFT", "K_TOPRIGHT",
  "K_TOPOFSTACKONLY", "K_TRISTATE", "K_TYPE", "K_UNATENESS", "K_UNITS",
  "K_USE", "K_VARIABLE", "K_VERTICAL", "K_VHI", "K_VIA", "K_VIARULE",
  "K_VLO", "K_VOLTAGE", "K_VOLTS", "K_WIDTH", "K_X", "K_Y", "T_STRING",
  "QSTRING", "NUMBER", "K_N", "K_S", "K_E", "K_W", "K_FN", "K_FS", "K_FE",
  "K_FW", "K_R0", "K_R90", "K_R180", "K_R270", "K_MX", "K_MY", "K_MXR90",
  "K_MYR90", "K_USER", "K_MASTERSLICE", "K_ENDMACRO", "K_ENDMACROPIN",
  "K_ENDVIARULE", "K_ENDVIA", "K_ENDLAYER", "K_ENDSITE", "K_CANPLACE",
  "K_CANNOTOCCUPY", "K_TRACKS", "K_FLOORPLAN", "K_GCELLGRID",
  "K_DEFAULTCAP", "K_MINPINS", "K_WIRECAP", "K_STABLE", "K_SETUP",
  "K_HOLD", "K_DEFINE", "K_DEFINES", "K_DEFINEB", "K_IF", "K_THEN",
  "K_ELSE", "K_FALSE", "K_TRUE", "K_EQ", "K_NE", "K_LE", "K_LT", "K_GE",
  "K_GT", "K_OR", "K_AND", "K_NOT", "K_DELAY", "K_TABLEDIMENSION",
  "K_TABLEAXIS", "K_TABLEENTRIES", "K_TRANSITIONTIME", "K_EXTENSION",
  "K_PROPDEF", "K_STRING", "K_INTEGER", "K_REAL", "K_RANGE", "K_PROPERTY",
  "K_VIRTUAL", "K_BUSBITCHARS", "K_VERSION", "K_BEGINEXT", "K_ENDEXT",
  "K_UNIVERSALNOISEMARGIN", "K_EDGERATETHRESHOLD1", "K_CORRECTIONTABLE",
  "K_EDGERATESCALEFACTOR", "K_EDGERATETHRESHOLD2", "K_VICTIMNOISE",
  "K_NOISETABLE", "K_EDGERATE", "K_OUTPUTRESISTANCE", "K_VICTIMLENGTH",
  "K_CORRECTIONFACTOR", "K_OUTPUTPINANTENNASIZE", "K_INPUTPINANTENNASIZE",
  "K_INOUTPINANTENNASIZE", "K_CURRENTDEN", "K_PWL",
  "K_ANTENNALENGTHFACTOR", "K_TAPERRULE", "K_DIVIDERCHAR", "K_ANTENNASIZE",
  "K_ANTENNAMETALLENGTH", "K_ANTENNAMETALAREA", "K_RISESLEWLIMIT",
  "K_FALLSLEWLIMIT", "K_FUNCTION", "K_BUFFER", "K_INVERTER",
  "K_NAMEMAPSTRING", "K_NOWIREEXTENSIONATPIN", "K_WIREEXTENSION",
  "K_MESSAGE", "K_CREATEFILE", "K_OPENFILE", "K_CLOSEFILE", "K_WARNING",
  "K_ERROR", "K_FATALERROR", "K_RECOVERY", "K_SKEW", "K_ANYEDGE",
  "K_POSEDGE", "K_NEGEDGE", "K_SDFCONDSTART", "K_SDFCONDEND", "K_SDFCOND",
  "K_MPWH", "K_MPWL", "K_PERIOD", "K_ACCURRENTDENSITY",
  "K_DCCURRENTDENSITY", "K_AVERAGE", "K_PEAK", "K_RMS", "K_FREQUENCY",
  "K_CUTAREA", "K_MEGAHERTZ", "K_USELENGTHTHRESHOLD", "K_LENGTHTHRESHOLD",
  "K_ANTENNAINPUTGATEAREA", "K_ANTENNAINOUTDIFFAREA",
  "K_ANTENNAOUTPUTDIFFAREA", "K_ANTENNAAREARATIO",
  "K_ANTENNADIFFAREARATIO", "K_ANTENNACUMAREARATIO",
  "K_ANTENNACUMDIFFAREARATIO", "K_ANTENNAAREAFACTOR",
  "K_ANTENNASIDEAREARATIO", "K_ANTENNADIFFSIDEAREARATIO",
  "K_ANTENNACUMSIDEAREARATIO", "K_ANTENNACUMDIFFSIDEAREARATIO",
  "K_ANTENNASIDEAREAFACTOR", "K_DIFFUSEONLY", "K_MANUFACTURINGGRID",
  "K_FIXEDMASK", "K_ANTENNACELL", "K_CLEARANCEMEASURE", "K_EUCLIDEAN",
  "K_MAXXY", "K_USEMINSPACING", "K_ROWMINSPACING", "K_ROWABUTSPACING",
  "K_FLIP", "K_NONE", "K_ANTENNAPARTIALMETALAREA",
  "K_ANTENNAPARTIALMETALSIDEAREA", "K_ANTENNAGATEAREA",
  "K_ANTENNADIFFAREA", "K_ANTENNAMAXAREACAR", "K_ANTENNAMAXSIDEAREACAR",
  "K_ANTENNAPARTIALCUTAREA", "K_ANTENNAMAXCUTCAR", "K_SLOTWIREWIDTH",
  "K_SLOTWIRELENGTH", "K_SLOTWIDTH", "K_SLOTLENGTH",
  "K_MAXADJACENTSLOTSPACING", "K_MAXCOAXIALSLOTSPACING",
  "K_MAXEDGESLOTSPACING", "K_SPLITWIREWIDTH", "K_MINIMUMDENSITY",
  "K_MAXIMUMDENSITY", "K_DENSITYCHECKWINDOW", "K_DENSITYCHECKSTEP",
  "K_FILLACTIVESPACING", "K_MINIMUMCUT", "K_ADJACENTCUTS",
  "K_ANTENNAMODEL", "K_BUMP", "K_ENCLOSURE", "K_FROMABOVE", "K_FROMBELOW",
  "K_IMPLANT", "K_LENGTH", "K_MAXVIASTACK", "K_AREAIO", "K_BLACKBOX",
  "K_MAXWIDTH", "K_MINENCLOSEDAREA", "K_MINSTEP", "K_ORIENT", "K_OXIDE1",
  "K_OXIDE2", "K_OXIDE3", "K_OXIDE4", "K_PARALLELRUNLENGTH", "K_MINWIDTH",
  "K_PROTRUSIONWIDTH", "K_SPACINGTABLE", "K_WITHIN", "K_ABOVE", "K_BELOW",
  "K_CENTERTOCENTER", "K_CUTSIZE", "K_CUTSPACING", "K_DENSITY", "K_DIAG45",
  "K_DIAG135", "K_MASK", "K_DIAGMINEDGELENGTH", "K_DIAGSPACING",
  "K_DIAGPITCH", "K_DIAGWIDTH", "K_GENERATED", "K_GROUNDSENSITIVITY",
  "K_HARDSPACING", "K_INSIDECORNER", "K_LAYERS", "K_LENGTHSUM",
  "K_MICRONS", "K_MINCUTS", "K_MINSIZE", "K_NETEXPR", "K_OUTSIDECORNER",
  "K_PREFERENCLOSURE", "K_ROWCOL", "K_ROWPATTERN", "K_SOFT",
  "K_SUPPLYSENSITIVITY", "K_USEVIA", "K_USEVIARULE", "K_WELLTAP",
  "K_ARRAYCUTS", "K_ARRAYSPACING", "K_ANTENNAAREADIFFREDUCEPWL",
  "K_ANTENNAAREAMINUSDIFF", "K_ANTENNACUMROUTINGPLUSCUT",
  "K_ANTENNAGATEPLUSDIFF", "K_ENDOFLINE", "K_ENDOFNOTCHWIDTH",
  "K_EXCEPTEXTRACUT", "K_EXCEPTSAMEPGNET", "K_EXCEPTPGNET", "K_LONGARRAY",
  "K_MAXEDGES", "K_NOTCHLENGTH", "K_NOTCHSPACING", "K_ORTHOGONAL",
  "K_PARALLELEDGE", "K_PARALLELOVERLAP", "K_PGONLY", "K_PRL", "K_TWOEDGES",
  "K_TWOWIDTHS", "IF", "LNOT", "'-'", "'+'", "'*'", "'/'", "UMINUS", "';'",
  "'('", "')'", "'='", "'\\n'", "'<'", "'>'", "$accept", "lef_file",
  "version", "$@1", "int_number", "dividerchar", "busbitchars", "rules",
  "end_library", "rule", "case_sensitivity", "wireextension", "fixedmask",
  "manufacturing", "useminspacing", "clearancemeasure", "clearance_type",
  "spacing_type", "spacing_value", "units_section", "start_units",
  "units_rules", "units_rule", "layer_rule", "start_layer", "$@2",
  "end_layer", "$@3", "layer_options", "layer_option", "$@4", "$@5", "$@6",
  "$@7", "$@8", "$@9", "$@10", "$@11", "$@12", "$@13", "$@14", "$@15",
  "$@16", "$@17", "$@18", "$@19", "$@20", "$@21", "$@22", "$@23", "$@24",
  "$@25", "$@26", "$@27", "$@28", "$@29", "layer_arraySpacing_long",
  "layer_arraySpacing_width", "layer_arraySpacing_arraycuts",
  "layer_arraySpacing_arraycut", "sp_options", "$@30", "$@31", "$@32",
  "$@33", "$@34", "$@35", "$@36", "layer_spacingtable_opts",
  "layer_spacingtable_opt", "layer_enclosure_type_opt",
  "layer_enclosure_width_opt", "$@37", "layer_enclosure_width_except_opt",
  "layer_preferenclosure_width_opt", "layer_minimumcut_within",
  "layer_minimumcut_from", "layer_minimumcut_length",
  "layer_minstep_options", "layer_minstep_option", "layer_minstep_type",
  "layer_antenna_pwl", "$@38", "layer_diffusion_ratios",
  "layer_diffusion_ratio", "layer_antenna_duo", "layer_table_type",
  "layer_frequency", "$@39", "$@40", "$@41", "ac_layer_table_opt", "$@42",
  "$@43", "dc_layer_table", "$@44", "int_number_list", "number_list",
  "layer_prop_list", "layer_prop", "current_density_pwl_list",
  "current_density_pwl", "cap_points", "cap_point", "res_points",
  "res_point", "layer_type", "layer_direction", "layer_minen_width",
  "layer_oxide", "layer_sp_parallel_widths", "layer_sp_parallel_width",
  "$@45", "layer_sp_TwoWidths", "layer_sp_TwoWidth", "$@46",
  "layer_sp_TwoWidthsPRL", "layer_sp_influence_widths",
  "layer_sp_influence_width", "maxstack_via", "$@47", "via", "$@48",
  "via_keyword", "start_via", "via_viarule", "$@49", "$@50", "$@51",
  "via_viarule_options", "via_viarule_option", "$@52", "via_option",
  "via_other_options", "via_more_options", "via_other_option", "$@53",
  "via_prop_list", "via_name_value_pair", "via_foreign", "start_foreign",
  "$@54", "orientation", "via_layer_rule", "via_layer", "$@55",
  "via_geometries", "via_geometry", "$@56", "end_via", "$@57",
  "viarule_keyword", "$@58", "viarule", "viarule_generate", "$@59",
  "viarule_generate_default", "viarule_layer_list", "opt_viarule_props",
  "viarule_props", "viarule_prop", "$@60", "viarule_prop_list",
  "viarule_layer", "via_names", "via_name", "viarule_layer_name", "$@61",
  "viarule_layer_options", "viarule_layer_option", "end_viarule", "$@62",
  "spacing_rule", "start_spacing", "end_spacing", "spacings", "spacing",
  "samenet_keyword", "maskColor", "irdrop", "start_irdrop", "end_irdrop",
  "ir_tables", "ir_table", "ir_table_values", "ir_table_value",
  "ir_tablename", "minfeature", "dielectric", "nondefault_rule", "$@63",
  "$@64", "$@65", "end_nd_rule", "nd_hardspacing", "nd_rules", "nd_rule",
  "usevia", "useviarule", "mincuts", "nd_prop", "$@66", "nd_prop_list",
  "nd_layer", "$@67", "$@68", "$@69", "$@70", "nd_layer_stmts",
  "nd_layer_stmt", "site", "start_site", "$@71", "end_site", "$@72",
  "site_options", "site_option", "site_class", "site_symmetry_statement",
  "site_symmetries", "site_symmetry", "site_rowpattern_statement", "$@73",
  "site_rowpatterns", "site_rowpattern", "$@74", "pt", "macro", "$@75",
  "start_macro", "$@76", "end_macro", "$@77", "macro_options",
  "macro_option", "$@78", "macro_prop_list", "macro_symmetry_statement",
  "macro_symmetries", "macro_symmetry", "macro_name_value_pair",
  "macro_class", "class_type", "pad_type", "core_type", "endcap_type",
  "macro_generator", "macro_generate", "macro_source", "macro_power",
  "macro_origin", "macro_foreign", "macro_fixedMask", "macro_eeq", "$@79",
  "macro_leq", "$@80", "macro_site", "macro_site_word", "site_word",
  "macro_size", "macro_pin", "start_macro_pin", "$@81", "end_macro_pin",
  "$@82", "macro_pin_options", "macro_pin_option", "$@83", "$@84", "$@85",
  "$@86", "$@87", "$@88", "$@89", "$@90", "$@91", "$@92",
  "pin_layer_oxide", "pin_prop_list", "pin_name_value_pair",
  "electrical_direction", "start_macro_port", "macro_port_class_option",
  "macro_pin_use", "macro_scan_use", "pin_shape", "geometries", "geometry",
  "$@93", "$@94", "geometry_options", "layer_exceptpgnet", "layer_spacing",
  "firstPt", "nextPt", "otherPts", "via_placement", "$@95", "$@96",
  "stepPattern", "sitePattern", "trackPattern", "$@97", "$@98", "$@99",
  "$@100", "trackLayers", "layer_name", "gcellPattern", "macro_obs",
  "start_macro_obs", "macro_density", "density_layers", "density_layer",
  "$@101", "$@102", "density_layer_rects", "density_layer_rect",
  "macro_clocktype", "$@103", "timing", "start_timing", "end_timing",
  "timing_options", "timing_option", "$@104", "$@105", "$@106",
  "one_pin_trigger", "two_pin_trigger", "from_pin_trigger",
  "to_pin_trigger", "delay_or_transition", "list_of_table_entries",
  "table_entry", "list_of_table_axis_dnumbers", "slew_spec", "risefall",
  "unateness", "list_of_from_strings", "list_of_to_strings", "array",
  "$@107", "start_array", "$@108", "end_array", "$@109", "array_rules",
  "array_rule", "$@110", "$@111", "$@112", "$@113", "$@114",
  "floorplan_start", "floorplan_list", "floorplan_element", "$@115",
  "$@116", "cap_list", "one_cap", "msg_statement", "$@117",
  "create_file_statement", "$@118", "def_statement", "$@119", "$@120",
  "$@121", "dtrm", "then", "else", "expression", "b_expr", "s_expr",
  "relop", "prop_def_section", "$@122", "prop_stmts", "prop_stmt", "$@123",
  "$@124", "$@125", "$@126", "$@127", "$@128", "$@129", "$@130",
  "prop_define", "opt_range_second", "opt_endofline", "$@131",
  "opt_endofline_twoedges", "opt_samenetPGonly", "opt_def_range",
  "opt_def_value", "opt_def_dvalue", "layer_spacing_opts",
  "layer_spacing_opt", "$@132", "layer_spacing_cut_routing", "$@133",
  "$@134", "$@135", "$@136", "$@137", "spacing_cut_layer_opt",
  "opt_adjacentcuts_exceptsame", "opt_layer_name", "$@138",
  "req_layer_name", "$@139", "universalnoisemargin", "edgeratethreshold1",
  "edgeratethreshold2", "edgeratescalefactor", "noisetable", "$@140",
  "end_noisetable", "noise_table_list", "noise_table_entry",
  "output_resistance_entry", "$@141", "num_list", "victim_list", "victim",
  "$@142", "vnoiselist", "correctiontable", "$@143", "end_correctiontable",
  "correction_table_list", "correction_table_item", "output_list", "$@144",
  "numo_list", "corr_victim_list", "corr_victim", "$@145", "corr_list",
  "input_antenna", "output_antenna", "inout_antenna", "antenna_input",
  "antenna_inout", "antenna_output", "extension_opt", "extension", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-1553)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-736)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    1383, -1553,   146,  2275, -1553, -1553,   -11, -1553, -1553, -1553,
     -11,   232, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553,    43, -1553, -1553,   -11,   -11,   -11,   -11,   -11,
     -11,   -11,   -11,   -11,    52,   417, -1553, -1553,     3,   103,
     120,   -11,  -354,   242,    79,   -11, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553,    60, -1553,    57, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553,   314, -1553,   167, -1553,
    -114,   193,   213,   -11,    16,    23,   218,   297,   305,   312,
     317,   332, -1553,    93,   375,   -11,   118,   157,   166,   202,
   -1553,   228,   234,   236,   238,   247,   256,   410,   455,   264,
     268,   279,   284, -1553, -1553, -1553,   287, -1553, -1553,   527,
    -163,   446,  1920,     2,   469,   633, -1553,   567, -1553, -1553,
     141,    47,     1,    50,   485,   607, -1553, -1553, -1553, -1553,
   -1553,   290, -1553, -1553, -1553, -1553, -1553,   295,   298,   342,
    1098, -1553,   308,   353, -1553, -1553, -1553, -1553,   412, -1553,
   -1553, -1553, -1553, -1553, -1553,   379,   411, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553,   424, -1553, -1553,   771,   799,   515,
     749,   833,   823,   831,   744,   614, -1553,   745,   891,   -11,
     -21,   -11, -1553,   -11,   -11,   -11,   -36,   -11,   -11,   -11,
      63,   -11, -1553,  -122,   -11,   -11,   315,   630,   -11, -1553,
     -11, -1553,   -11,   -11, -1553,   -11, -1553,   -11,   -11,   -11,
     -11,   -11,   -11,   -11,   -11,   -11,   -11,   -11,   -11,   -11,
     -11,   -11, -1553,   271,   -11,   751,   -11,   -11,   -11,   506,
     -11,   -11,   -11,   -11,   -11, -1553,   271, -1553,   492,   -11,
     496,   -11, -1553, -1553, -1553, -1553, -1553, -1553,   -11, -1553,
   -1553, -1553, -1553,   889, -1553, -1553, -1553,   578, -1553, -1553,
   -1553, -1553,   760, -1553,   350,   382,   797, -1553, -1553, -1553,
     770,   880,   772, -1553, -1553, -1553,    45, -1553,   -11, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553,   164, -1553, -1553,   774,
     775, -1553, -1553,   -78, -1553,   -11, -1553,   -11,    95, -1553,
   -1553, -1553,   385,   516,   888,   624,   916, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
     784, -1553, -1553, -1553, -1553,   494, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553,   785, -1553,   -11, -1553,   919, -1553,
   -1553, -1553, -1553, -1553,   571,   -49,  -116,   275,   727, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
     454,   519,  -116,  -116, -1553,   791,   -11,   -11,   -11, -1553,
     -11,   -11,   -11,   -11,   813,   531,   -61,   536, -1553, -1553,
   -1553, -1553,   556,   558,   826,   566,   -73,   -55,   -41,   572,
     574, -1553,   575, -1553, -1553, -1553, -1553, -1553, -1553,   576,
     580,   835,   598,   -11,   581,   608,   609, -1553, -1553, -1553,
     -11,    -1,   610,    22,   612,    22, -1553,   622,    22,   623,
      22, -1553,   626,   636,   637,   638,   639,   641,   643,   644,
     645,   648,   -11,   650,   651,   876,   374, -1553, -1553,   -11,
     661, -1553, -1553,   662,   724,   711,   -20,   663,   664,   665,
     -47,   667,   -78,   -11,   586,   -78,   668, -1553,   669,   932,
     933,   672,   953,   954, -1553, -1553,   404, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553,   -11,   -11,   693,   652,   493,   567,
     694,   611, -1553,   957,  1091,   143, -1553, -1553,   359,   -11,
     -11,   -78,   -11,   -11,   -11,   -11, -1553, -1553,   959, -1553,
   -1553,   -39,   698,   699,   700,   963,  1121,   -89, -1553,   148,
       7, -1553,   787,   573,    -4, -1553, -1553, -1553, -1553,   705,
     969,   971,   978,   717,   980,   723,   987,   726,  1146,   728,
     729,   731,   -79,   994,   735,   736, -1553, -1553, -1553, -1553,
     737,   796, -1553, -1553,     0,   738,  2296, -1553, -1553,   794,
     794,   794,   -30,   -11,  1137, -1553, -1553,  1557,  1006,  1006,
     654, -1553,   680, -1553,  1006, -1553, -1553,   243,   746, -1553,
     275,   -49,   -49, -1553,   299, -1553,   275,  -116,   127,   275,
   -1553, -1553,   275,   275,   843,   278,   869, -1553,  1026,  1027,
    1029,  1030,  1031,  1032,  1033,  1034,   -11, -1553,   169, -1553,
   -1553,   -11, -1553,   230, -1553, -1553,   127,   127,  1035,   776,
     777,   778,   779,   780,   790,   792,   793, -1553,   773,   800,
   -1553, -1553, -1553, -1553, -1553, -1553,   801, -1553,   802,   798,
     803, -1553, -1553,   -80, -1553, -1553, -1553,   646,   -87, -1553,
     805,   -11, -1553, -1553, -1553,   804,   882,   -11,  1045,   807,
   -1553,   808, -1553,   815, -1553,   817,   894, -1553,   820, -1553,
     824,   894, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553,   825, -1553, -1553,   -11, -1553, -1553, -1553, -1553,
     834,   -11, -1553,  1056, -1553, -1553,   -11,   -11,  1059,   -11,
    1073,   837, -1553, -1553, -1553, -1553,   844, -1553, -1553, -1553,
     -11, -1553,  1081,   -78, -1553, -1553, -1553,   845, -1553,   850,
     649,   -74, -1553,  1092, -1553,   -11, -1553, -1553, -1553,   851,
     794,   794, -1553,   300, -1553, -1553, -1553, -1553,   143,   852,
   -1553, -1553, -1553,   853,   854,   855,   857,   -78,   858,  1250,
    1126,   -11,   -11, -1553,   -11, -1553, -1553, -1553, -1553, -1553,
     -11, -1553, -1553, -1553, -1553, -1553,   -66, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553,   859,   860,   861, -1553,   862, -1553, -1553, -1553,
     -11, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,   671,
     -64, -1553, -1553, -1553,  1125,   100, -1553, -1553,   864,  1128,
   -1553,   -11, -1553,   -11,    56,   568, -1553,   -11,   -11, -1553,
    1129,   -11, -1553,   -11,   -11, -1553, -1553, -1553,   -11,   -11,
     -11,   -11,   -11,   -11,   -11,   462,   357,   -11,   763,   -11,
     -11, -1553, -1553,  1130,   -11,   -11,  1131,   -11,   -11,  1133,
    1136,  1138,  1139,  1140,  1141,  1143,  1144, -1553, -1553, -1553,
   -1553,   -91, -1553, -1553, -1553,   133,  1132,   -11,   -46,   -44,
     -25,   794,   -78,   877, -1553,   479,  1153, -1553,   -11,   -11,
     -11,   -11, -1553, -1553,   -11,   -11,   -11,   -11, -1553,   117,
    1111, -1553, -1553,   -11,   885,   890, -1553, -1553,  1149,  1150,
    1154, -1553, -1553, -1553, -1553, -1553,  1101,   594,   237,   -11,
     892,   896,   -11,   -11,   897,   -11,   -11,   898,    78,   899,
    1169,  1172, -1553, -1553, -1553, -1553,    58,   281, -1553,   199,
     -49,   -49,   -49,   -49, -1553, -1553, -1553,   281,   184,  -116,
   -1553,   281, -1553,   816,    81,   830, -1553, -1553, -1553, -1553,
   -1553, -1553,   -49,   -49, -1553, -1553, -1553,   -49,   275,   275,
     275,   275, -1553,  -116,  -116,  -116,   306,   306,   306,   306,
     306,   306,   306,   306,   911,   -11,  1100,    98, -1553,   913,
     -11,  1099,    98, -1553, -1553, -1553,   914, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553,   915, -1553, -1553, -1553,   917,
   -1553, -1553, -1553, -1553,    55,   -80, -1553, -1553, -1553, -1553,
   -1553,   -11,   420, -1553,   918, -1553,  1177, -1553, -1553, -1553,
   -1553,   -78, -1553, -1553, -1553,   921, -1553, -1553,   922, -1553,
   -1553, -1553, -1553,   -11,   923,   -65,  1184,  1223,   -11, -1553,
     -11, -1553, -1553,  -102, -1553,   -11,   984, -1553, -1553,   986,
   -1553, -1553, -1553, -1553, -1553, -1553,   927, -1553, -1553,   -78,
    1091, -1553,   -93, -1553,  1194, -1553, -1553, -1553, -1553,   935,
   -1553,   -11,   -11,   936,   -94, -1553,   937,  1007, -1553, -1553,
   -1553, -1553, -1553, -1553,   938, -1553, -1553, -1553, -1553, -1553,
     939, -1553, -1553, -1553, -1553,  1007,   943,   944,   945,   946,
     947,   948,  -100,  1195,   949,   950,   -11,  1215,   955,  1216,
     956,   958,  1219,   -11,   960,   961,   962,   964,   965,   966,
     967, -1553, -1553, -1553, -1553,   968, -1553, -1553, -1553,   970,
     973, -1553, -1553, -1553, -1553, -1553, -1553,   974,   975,   976,
    1225,   -11,   977,  1337,  1337,  1337,   981,   982,  1337,  1337,
    1337,  1337,  1339,  1339,  1337,  1339,   409,  1240,  1241,  1246,
     -85, -1553,   842,   164, -1553,   479, -1553, -1553,   -78, -1553,
     -78,   -78,   -78,   -78,   -78, -1553, -1553, -1553, -1553,   -11,
     -11,   -11,   -11,  1247,   -11,   -11,   -11,   -11,  1248, -1553,
   -1553, -1553,   989,   -11, -1553,     4,   -11,   461, -1553, -1553,
     990,   991,   992,   -11, -1553, -1553, -1553,   602,   -11, -1553,
   -1553,  1263, -1553, -1553,  1391,  1392, -1553,  1393,  1395, -1553,
    1228,   -11, -1553, -1553, -1553, -1553,  1006,  1006, -1553,   689,
   -1553,  1261,  1262,  1264, -1553, -1553,  1396, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553,  1224,   -49, -1553,   457,   457, -1553,
   -1553,  -116, -1553, -1553,   275, -1553,   363,   363,   363, -1553,
   -1553,   679,   378,  1009,  1009,  1009,  1267,  1200,  1200,  1269,
    1008,  1010,  1011,  1013,  1014,  1016,  1017,  1018, -1553, -1553,
      11, -1553, -1553, -1553, -1553,    18, -1553, -1553, -1553,   -11,
     464, -1553,   -11,   466, -1553,  1020,  1280, -1553,   -11,   -11,
     -11,   -11,   -11,   -11,  1021, -1553,   -11,  1022, -1553,  1023,
   -1553, -1553, -1553,   -78, -1553, -1553,  1086,   -40, -1553, -1553,
   -1553, -1553,   -11, -1553,   -11, -1553, -1553, -1553,   -11,   -11,
    1087, -1553,  1037, -1553, -1553, -1553,  1294, -1553,   -11, -1553,
     -11, -1553,   -78,   -78, -1553, -1553, -1553, -1553, -1553,  1036,
    1038, -1553,  1039, -1553, -1553, -1553, -1553, -1553,  1428, -1553,
   -1553, -1553, -1553, -1553, -1553,  1040, -1553, -1553, -1553, -1553,
     -11,  1041, -1553,  1042, -1553, -1553,  1043,   -11, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553,   719,   -62, -1553,   -11, -1553, -1553,  1044,  1046,  1047,
   -1553, -1553,  1048,  1049,  1050,  1051, -1553,  1053,  1054,  1055,
    1057, -1553, -1553, -1553, -1553,  1060,  1061,  1062,  1064, -1553,
     934, -1553,  1069,  1070,  1426,  1052, -1553,   -52,   -78,   -78,
     -78,  1071, -1553,  1296,  1072,  1074,  1075,  1076, -1553,   -60,
    1077,  1079,  1080,  1083, -1553,   -56, -1553,  1260, -1553, -1553,
     -11, -1553, -1553, -1553, -1553, -1553,   -11, -1553, -1553, -1553,
    1243,   -11,   117,   -11,   -11,   -11,   -11, -1553,  1273,  1084,
    1088,  1318, -1553, -1553, -1553,   326,   -11,  1089,  1093,  1321,
   -1553, -1553,  -121,   -29,   716,   331,   351, -1553,   -11,   -11,
    1333, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
    1268, -1553,  1275, -1553,   -11,  1095, -1553,   -11,  1096, -1553,
   -1553, -1553, -1553,  1348,   -11,  1289,  1159,  1160,  1114, -1553,
   -1553,  1102, -1553, -1553, -1553,    21,    25, -1553,   -11,   552,
     -11,   -11,  1103, -1553, -1553,  1105, -1553,   -11,   -11,   -11,
     -11,   -11,  1107, -1553,   -98,   -11,   -78,  1110, -1553, -1553,
   -1553, -1553,  1423,   -11, -1553,  1112, -1553, -1553, -1553,  1115,
   -1553, -1553, -1553, -1553, -1553,  1117,  1363, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553,  1364, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553,  1118, -1553, -1553, -1553, -1553,    39,
     -19, -1553,   -78, -1553,  1509, -1553,  1381,  1120, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
     -11,   -11,   -11,   -11, -1553,  1327,  1416,  1417,  1418,  1419,
     -11, -1553, -1553, -1553, -1553,   -88,  1134, -1553, -1553, -1553,
   -1553,  1345,   -49,  -116,   275,   -11, -1553, -1553, -1553, -1553,
     -11,  1268, -1553,   -11,  1275, -1553,  1127, -1553,  1135, -1553,
   -1553, -1553,   -11,   -11,   -11,   -11, -1553,    29, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553,  1206, -1553, -1553, -1553,
   -1553,  1193,  1433,  1403, -1553, -1553, -1553, -1553,  1173,  1145,
   -1553, -1553,  1148,   -78, -1553,   -78, -1553,  1565, -1553, -1553,
   -1553, -1553, -1553, -1553,   -11,   -11,  1151,   -11,  1156, -1553,
     -12,  1157,  1509, -1553,    91,  1147,  1158,   -11,   -11,   -11,
     -11,   -11,   -11,   -11,  1163,  1415, -1553, -1553, -1553, -1553,
     363,  1009,   490, -1553,  1170, -1553,  1171, -1553, -1553, -1553,
    1449,    35,   -11, -1553, -1553,  1185, -1553,  1358,  1358,   -59,
     -11,  1175,  1197,   -11,  1178,  1193,   -11,   -11, -1553,   -11,
    1179,  1173, -1553,  1217, -1553,   -78,  1423,   -11, -1553, -1553,
   -1553,  1591, -1553,   -19, -1553, -1553,  1180,  1181, -1553, -1553,
     -11,   -11,  1450,   -11,  1183,  1186, -1553, -1553, -1553,   -11,
   -1553, -1553, -1553, -1553,   -11,   -11, -1553, -1553, -1553,  1204,
    1198,   -11,   -13,   -11, -1553, -1553, -1553,  1249, -1553,   -11,
   -1553,  1486, -1553, -1553, -1553, -1553,   -11,  1490, -1553, -1553,
   -1553, -1553,    -3,   -11, -1553,  1492,   -11,  1201, -1553, -1553,
    1202,   -11,   -11,   -11,  1568,  1571,  1207,  1388,  1394,  1407,
     -11, -1553, -1553,   -11, -1553, -1553,   -11,  1476,  1421, -1553,
     -11, -1553,   -11, -1553, -1553,  1483,   -11,  1481, -1553,  1220,
     -11,  1519, -1553, -1553,   -11,   -11,  1226, -1553, -1553, -1553,
     -11,   -11,   -11, -1553,  1290, -1553, -1553,  1487, -1553, -1553,
   -1553,  1494,   -11,   -11, -1553,  1483, -1553,  1493, -1553,   -11,
     -11,   -11,  1230, -1553, -1553, -1553, -1553, -1553,    49, -1553,
      53,   -11,   -11, -1553, -1553, -1553,    74,   -11, -1553, -1553,
    1037, -1553,  1495, -1553,   -11,   -11, -1553,  1497,  1497,    34,
   -1553, -1553, -1553, -1553, -1553, -1553,    75,    89, -1553, -1553,
    1303,  1511,   -11,  1232, -1553,   -11, -1553, -1553,  1654,   -11,
   -1553,  1556,   -11,   -11,   -11, -1553,  1258, -1553, -1553,    90,
     -11,   -11, -1553, -1553,  1307, -1553,   -11,  1259,  1520,   -11,
    1270,  1272,  1274, -1553, -1553, -1553,  1558, -1553, -1553,   -11,
    1276, -1553, -1553,  1277, -1553, -1553, -1553,   -11, -1553,   -11,
     -11, -1553, -1553, -1553,   -11,  1278,  1343,   -11,   -11,   -11,
     -11,  1279, -1553, -1553,   -23,   -11,   -11, -1553,   -11, -1553,
     -11,   -11,  1524,   -11,   -11,  1282,  1283,  1284,   -11, -1553,
   -1553, -1553,  1285, -1553
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_int16 yydefact[] =
{
       0,    10,     0,  1018,     1,   819,     0,   420,    79,   501,
       0,     0,   431,   472,   410,    67,   307,   370,   852,   854,
     856,   901,     0,     3,  1020,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   848,   850,     0,     0,
       0,     0,     0,     0,     0,     0,    13,    21,    14,     9,
      15,    22,    47,    46,    48,    49,    16,    68,    17,    83,
      50,    18,     0,   305,     0,    19,    20,    24,   412,    27,
     422,    26,    25,    32,    28,   476,    29,   505,    30,   823,
      23,    51,    31,    33,    34,    35,    37,    36,    38,    39,
      40,    41,    42,    43,    44,    45,    11,  1019,     0,     5,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   903,     0,     0,     0,     0,     0,     0,     0,
     978,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    56,    61,    60,     0,    62,    63,     0,
       0,     0,     0,   308,     0,   375,   394,   391,   377,   396,
       0,     0,     0,   499,   817,     0,     2,   820,   430,    80,
     502,     0,    53,    52,   432,   473,   371,     0,     0,     0,
       0,     7,     0,     0,   975,   995,   977,   976,     0,  1013,
    1012,  1014,     6,    55,    54,     0,     0,  1015,  1016,  1017,
      57,    59,    65,    64,     0,   303,   302,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    69,     0,     0,     0,
       0,     0,    81,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   119,     0,     0,     0,     0,     0,     0,   132,
       0,   135,     0,     0,   140,     0,   143,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   147,   207,     0,     0,     0,     0,     0,   177,
       0,     0,     0,     0,     0,   185,   207,    85,     0,     0,
       0,     0,    78,    84,   309,   310,   342,   361,     0,   332,
     311,   330,   322,     0,   323,   325,   327,     0,   328,   363,
     376,   373,     0,   378,   379,   390,     0,   416,   409,   413,
       0,     0,     0,   419,   423,   425,     0,   474,     0,   486,
     491,   471,   477,   480,   479,   481,     0,   753,   586,     0,
       0,   588,   743,     0,   596,     0,   592,     0,     0,   532,
     756,   527,     0,     0,     0,     0,     0,   506,   511,   507,
     508,   509,   510,   514,   513,   515,   512,   516,   517,   519,
       0,   518,   520,   600,   523,     0,   524,   525,   526,   758,
     593,   827,   829,   831,     0,   834,     0,   825,     0,   824,
     839,   837,    12,   429,   437,     0,     0,     0,     0,   907,
     917,   905,   911,   919,   909,   913,   915,   904,     4,   974,
       0,     0,     0,     0,    58,     0,     0,     0,     0,    66,
       0,     0,     0,     0,     0,     0,     0,     0,   279,   280,
     281,   282,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   100,     0,   274,   275,   273,   276,   277,   278,     0,
       0,     0,     0,     0,     0,     0,     0,   242,   241,   243,
     121,     0,     0,     0,     0,     0,   137,     0,     0,     0,
       0,   145,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   208,   209,     0,
       0,   170,   174,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   187,     0,     0,   149,     0,     0,
       0,     0,     0,     0,   368,   306,   324,   344,   346,   347,
     345,   348,   350,   351,   349,   352,   353,   354,   355,   358,
     356,   359,   357,   338,     0,     0,     0,     0,   360,     0,
       0,     0,   383,     0,     0,   380,   381,   392,     0,     0,
       0,     0,     0,     0,     0,     0,   397,   411,     0,   421,
     428,     0,     0,     0,     0,     0,     0,     0,   493,   544,
     552,   553,   541,     0,   549,   543,   550,   547,   548,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   585,   747,   745,   581,
       0,     0,   503,   500,     0,     0,     0,   742,   699,   417,
     417,   417,   417,     0,     0,   710,   709,     0,     0,     0,
       0,   838,     0,   845,     0,   821,   818,     0,     0,   439,
       0,     0,     0,   872,   858,   891,     0,     0,   858,     0,
     887,   886,     0,     0,     0,   858,     0,   902,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1002,     0,   998,
    1001,     0,   985,     0,   981,   984,   858,   858,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    99,     0,     0,
     113,   104,   114,    82,   109,    94,     0,    90,     0,     0,
       0,   183,   112,   944,   111,    88,    98,     0,     0,   259,
       0,     0,   116,   115,   110,     0,     0,     0,     0,     0,
     131,     0,   233,     0,   134,     0,   239,   139,     0,   142,
       0,   239,   155,   156,   157,   158,   159,   160,   161,   162,
     163,   164,     0,   166,   167,     0,   285,   286,   287,   288,
       0,     0,   168,   283,   225,   169,     0,     0,     0,     0,
       0,     0,    89,   184,    97,    92,     0,    96,   717,   719,
       0,   188,   189,     0,   151,   150,   343,     0,   329,     0,
       0,     0,   333,     0,   326,     0,   497,   341,   339,     0,
     417,   417,   364,   379,   395,   387,   388,   389,     0,     0,
     407,   372,   382,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   424,     0,   426,   483,   482,   484,   475,
       0,   488,   489,   490,   485,   487,     0,   545,   546,   562,
     565,   563,   564,   566,   567,   554,   542,   572,   573,   569,
     568,   570,   571,   555,   558,   556,   557,   559,   560,   561,
     551,   540,     0,     0,     0,   574,     0,   580,   597,   579,
       0,   578,   577,   576,   534,   535,   536,   531,   533,     0,
       0,   529,   521,   522,     0,     0,   584,   582,     0,     0,
     590,     0,   591,     0,     0,     0,   598,     0,     0,   629,
       0,     0,   608,     0,     0,   623,   625,   681,     0,     0,
       0,     0,     0,     0,     0,     0,   694,     0,     0,     0,
       0,   642,   627,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   659,   665,   661,
     663,     0,   595,   601,   611,   682,     0,     0,     0,     0,
       0,   417,     0,     0,   741,   698,     0,   809,     0,     0,
       0,     0,   760,   808,     0,     0,     0,     0,   762,     0,
       0,   788,   789,     0,     0,     0,   790,   791,     0,     0,
       0,   785,   786,   787,   755,   759,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   841,   843,   840,   438,   433,     0,   869,     0,
       0,     0,     0,     0,   859,   860,   853,     0,     0,     0,
     855,     0,   883,     0,     0,     0,   896,   897,   892,   893,
     894,   895,     0,     0,   898,   899,   900,     0,     0,     0,
       0,     0,   857,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   858,   999,     0,
       0,     0,   858,   982,   849,   851,     0,    71,    74,    76,
      73,    72,    70,    75,    77,     0,   107,    95,    91,     0,
     105,   947,   946,   949,   950,   944,   261,   262,   263,   120,
     260,     0,     0,   264,     0,   123,     0,   122,   128,   125,
     124,     0,   133,   136,   240,     0,   141,   144,     0,   165,
     172,   148,   179,     0,     0,     0,     0,     0,     0,   194,
       0,   178,    93,     0,   181,     0,     0,   152,   362,     0,
     337,   336,   335,   331,   334,   369,     0,   340,   366,     0,
       0,   385,     0,   393,     0,   398,   399,   406,   405,     0,
     404,     0,     0,     0,     0,   427,     0,     0,   492,   494,
     754,   587,   575,   589,     0,   539,   538,   537,   528,   530,
       0,   744,   746,   583,   504,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   690,   691,   692,   693,     0,   695,   697,   696,     0,
       0,   685,   688,   689,   687,   686,   684,     0,     0,     0,
       0,     0,     0,   969,   969,   969,     0,     0,   969,   969,
     969,   969,     0,     0,   969,     0,     0,     0,     0,     0,
       0,   602,     0,     0,   645,     0,   700,   418,     0,   719,
       0,     0,     0,     0,     0,   721,   702,   711,   757,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   810,
     811,   812,     0,     0,   803,     0,     0,     0,   800,   784,
       0,     0,     0,     0,   792,   793,   794,     0,     0,   798,
     799,     0,   828,   830,     0,     0,   832,     0,     0,   835,
       0,     0,   846,   826,   822,   833,     0,     0,   458,     0,
     451,     0,     0,     0,   442,   443,     0,   440,   445,   446,
     447,   444,   441,   861,     0,     0,   870,   866,   865,   867,
     868,     0,   889,   888,     0,   884,   875,   874,   873,   879,
     880,   882,   881,   878,   877,   876,   923,   938,   938,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1000,  1004,
       0,   997,   996,   983,   987,     0,   980,   979,   304,     0,
       0,   267,     0,     0,   270,   936,     0,   951,     0,     0,
       0,     0,     0,     0,     0,   945,     0,     0,   265,     0,
     244,   255,   257,     0,   138,   146,   218,   210,   284,   171,
     232,   230,     0,   231,     0,   175,   226,   227,     0,     0,
       0,   255,   297,   186,   718,   720,   216,   190,     0,   236,
       0,   498,     0,     0,   374,   384,   386,   408,   402,     0,
       0,   400,     0,   414,   478,   495,   594,   748,   727,   631,
     636,   637,   680,   679,   676,     0,   677,   599,   618,   616,
       0,     0,   614,     0,   632,   633,     0,     0,   610,   635,
     634,   617,   615,   638,   639,   613,   622,   621,   612,   620,
     619,     0,     0,   671,     0,   641,   970,     0,     0,     0,
     649,   650,     0,     0,     0,     0,   972,     0,     0,     0,
       0,   667,   668,   669,   670,     0,     0,     0,     0,   605,
       0,   603,     0,     0,     0,   712,   719,     0,     0,     0,
       0,     0,   723,     0,     0,     0,     0,     0,   813,     0,
       0,     0,     0,     0,   815,     0,   777,     0,   767,   804,
       0,   768,   801,   781,   782,   783,     0,   795,   796,   797,
       0,     0,     0,     0,     0,     0,     0,   836,     0,     0,
       0,     0,   455,   456,   457,     0,     0,     0,     0,   435,
     434,   862,     0,     0,     0,     0,     0,   924,     0,   942,
     940,   925,   908,   918,   906,   912,   920,   910,   914,   916,
       0,  1005,     0,   988,     0,     0,   268,     0,     0,   271,
     937,   948,   956,     0,     0,   959,     0,     0,     0,   963,
     101,     0,   117,   118,   257,     0,     0,   234,     0,   220,
       0,     0,     0,   228,   229,     0,   102,     0,   195,     0,
       0,     0,     0,    86,     0,     0,     0,     0,   403,   401,
     415,   496,     0,     0,   678,     0,   640,   609,   624,     0,
     675,   674,   673,   643,   672,     0,     0,   646,   648,   647,
     651,   652,   655,   654,     0,   656,   657,   653,   658,   660,
     666,   662,   664,   606,     0,   604,   683,   644,   713,   714,
       0,   703,     0,   719,     0,   705,     0,     0,   772,   776,
     774,   770,   814,   761,   771,   769,   773,   775,   816,   763,
       0,     0,     0,     0,   764,     0,     0,     0,     0,     0,
       0,   842,   844,   459,   453,     0,     0,   448,   449,   436,
     863,     0,     0,     0,     0,     0,   943,   921,   941,   922,
       0,  1003,  1006,     0,   986,   989,     0,   108,     0,   106,
     952,   957,     0,     0,     0,     0,   266,     0,   129,   256,
     258,   126,   236,   219,   221,   222,   223,   211,   213,   180,
     176,   204,     0,     0,   298,   199,   217,   182,   191,     0,
     237,   238,     0,     0,   365,     0,   750,     0,   630,   626,
     628,   971,   973,   607,     0,     0,     0,     0,     0,   719,
       0,     0,     0,   722,     0,     0,     0,     0,   805,     0,
       0,     0,     0,     0,     0,     0,   452,   454,   450,   864,
     871,   890,   885,   939,     0,  1007,     0,   990,   269,   272,
     965,   926,     0,   954,   961,     0,   245,     0,     0,     0,
       0,     0,   214,     0,     0,   204,     0,     0,   255,     0,
       0,   191,   153,     0,   719,     0,   749,     0,   716,   715,
     701,     0,   704,     0,   707,   706,     0,     0,   802,   780,
       0,     0,     0,     0,   728,   731,   739,   740,   847,     0,
    1008,   991,   966,   953,     0,     0,   927,   958,   960,   967,
     931,     0,   248,     0,   130,   127,   235,     0,   173,     0,
     212,     0,   103,   205,   202,   196,   200,     0,    87,   192,
     154,   312,     0,     0,   751,     0,     0,     0,   724,   778,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   928,
       0,   968,   955,     0,   962,   964,     0,     0,     0,   253,
       0,   215,     0,   299,   255,   293,     0,     0,   367,     0,
       0,     0,   708,   779,     0,     0,     0,   729,   732,   460,
       0,     0,     0,   930,     0,   251,   249,     0,   255,   224,
     206,   203,   197,     0,   201,   293,   193,     0,   752,     0,
       0,   806,     0,   766,   736,   736,   463,  1010,     0,   993,
       0,     0,     0,   255,   257,   246,     0,     0,   300,   289,
     297,   294,     0,   726,     0,     0,   765,   730,   733,     0,
    1009,  1011,   992,   994,   929,   932,     0,     0,   257,   254,
       0,   198,     0,     0,   725,     0,   738,   737,     0,     0,
     461,     0,     0,     0,     0,   464,   934,   252,   250,     0,
       0,     0,   290,   295,     0,   807,     0,     0,     0,     0,
       0,     0,     0,   935,   933,   247,     0,   291,   255,     0,
       0,   469,   462,     0,   465,   466,   470,     0,   255,   296,
       0,   468,   467,   301,   292,     0,     0,     0,     0,     0,
       0,     0,   313,   315,   314,     0,     0,   320,     0,   316,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   318,
     321,   317,     0,   319
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -1553, -1553, -1553, -1553,    -6, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,   -96, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,   -86, -1553,
    1430, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553,   219, -1553,    10, -1553,  1000, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553,   -58, -1553, -1343, -1552, -1553,  1058,
   -1553,   653, -1553,   388, -1553,   386, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553,  -194, -1553, -1553,  -218, -1553, -1553, -1553,
   -1553,   757, -1553,  1431, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553,  1231, -1553, -1553,   979, -1553,
    -148, -1553,  -322, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553,  1205,   972, -1553,
    -519, -1553, -1553,  -144, -1553, -1553, -1553, -1553, -1553, -1553,
     632, -1553,   767, -1553, -1553, -1553, -1553, -1553,  -433, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1397, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553,  -211, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
     886, -1553,   534, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
     302, -1553, -1553, -1553, -1553, -1553, -1553,   535,   836, -1553,
   -1553, -1553, -1553, -1553,  -898, -1181, -1202, -1553, -1553, -1553,
   -1469,  -590, -1553, -1553, -1553, -1553, -1553,  -186, -1553, -1553,
   -1553, -1553, -1553, -1553,   903, -1553, -1553, -1553,   -43, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553,   512, -1553, -1553,     6,   250,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553,  1600, -1553, -1553, -1553,
    -589,  -400,  -834,  -334,  -576,  -332,  -624, -1553, -1553, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,  -265,
   -1553, -1553, -1553, -1553, -1553,   447, -1553, -1553,   712, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,  -540,
   -1553,  -894, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553,
   -1553,  1113, -1553, -1553, -1553, -1553,    76, -1553, -1553, -1553,
   -1553, -1553, -1553,  1123, -1553, -1553, -1553, -1553,    77, -1553,
   -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553, -1553
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
       0,     2,    46,   114,   515,    47,    48,     3,   156,    49,
      50,    51,    52,    53,    54,    55,   136,   139,   194,    56,
      57,   141,   206,    58,    59,   101,   272,   414,   142,   273,
     484,  1718,   673,  1711,   431,   686,  1352,  1778,  1351,  1777,
     443,   445,   696,   448,   450,   701,   466,  1379,  1850,   723,
    1356,   724,   476,  1357,  1376,   482,   742,  1086,  1790,  1791,
     731,  1371,  1713,  1884,  1939,  1788,  1885,  1883,  1784,  1785,
     469,  1572,  1782,  1840,  1582,  1569,  1706,  1781,  1075,  1366,
    1367,   693,  1702,  1584,  1720,  1065,   440,  1057,  1564,  1832,
    1958,  1878,  1934,  1933,  1834,  1908,  1565,  1566,   678,   679,
    1052,  1053,  1330,  1331,  1333,  1334,   429,   412,  1074,   720,
    1961,  1982,  2008,  1914,  1915,  1998,  1580,  1911,  1938,    60,
     395,    61,   144,    62,    63,   282,   492,  1887,  2023,  2024,
    2029,  2032,   283,   284,   496,   285,   493,   751,   752,   286,
     287,   489,   516,   288,   289,   490,   518,   762,  1382,   495,
     753,    64,   108,    65,    66,   519,   291,   147,   524,   525,
     526,   768,  1102,   148,   294,   527,   149,   292,   295,   536,
     771,  1104,    67,    68,   298,   150,   299,   300,   908,    69,
      70,   303,   151,   304,   541,   785,   305,    71,    72,    73,
     106,   374,  1276,  1520,   609,   966,  1277,  1278,  1279,  1280,
    1281,  1515,  1665,  1282,  1511,  1755,  1926,  1988,  1949,  1975,
      74,    75,   107,   311,   545,   152,   312,   313,   314,   547,
     795,   315,   548,   796,  1119,  1591,  1374,    76,   336,    77,
     102,   583,   849,   153,   337,   573,   840,   338,   572,   838,
     841,   339,   559,   820,   805,   813,   340,   341,   342,   343,
     344,   345,   346,   347,   561,   348,   564,   349,   350,   367,
     351,   352,   353,   566,   902,  1143,   586,   903,  1149,  1152,
    1153,  1181,  1146,  1180,  1196,  1198,  1199,  1197,  1455,  1432,
    1433,   904,   905,  1205,  1177,  1165,  1169,   594,   595,   906,
    1465,   915,  1629,  1736,   739,  1375,  1083,   596,  1473,  1636,
    1738,   585,   954,  1864,  1924,  1865,  1925,  1947,  1967,   957,
     354,   355,   356,   845,   578,   844,  1592,  1796,  1726,   357,
     560,   358,   359,   944,   597,   945,  1223,  1228,  1748,   946,
     947,  1247,  1500,  1251,  1237,  1238,  1235,  1812,   948,  1232,
    1479,  1485,    78,   368,    79,    98,   606,   960,   154,   369,
     604,   598,   599,   600,   602,   370,   607,   964,  1266,  1267,
     958,  1262,    80,   127,    81,   128,    82,   109,   110,   111,
     976,  1285,  1672,   624,   625,   626,   997,    83,   112,   170,
     387,   630,   628,   633,   631,   634,   635,   629,   632,  1310,
    1827,  1874,  1976,  1994,  1551,  1529,  1679,  1677,  1044,  1045,
    1335,  1344,  1553,  1770,  1829,  1771,  1830,  1823,  1872,  1437,
    1606,  1447,  1614,    84,    85,    86,    87,    88,   178,  1022,
     643,   644,   645,  1020,  1325,  1684,  1685,  1868,  1930,    89,
     390,  1017,   638,   639,   640,  1015,  1320,  1681,  1682,  1867,
    1928,    90,    91,    92,    93,    94,    95,    96,    97
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     100,   728,  1005,   293,   103,   335,   772,  1467,   950,   951,
    1209,  1211,  1697,   580,   959,  1737,  1208,   306,  1210,   115,
     116,   117,   118,   119,   120,   121,   122,   123,  1578,   980,
    1469,   274,   911,   408,   967,   132,  1002,  1212,   307,   140,
     977,   614,  1392,   981,   618,  1041,   982,   984,  1968,  1200,
      99,   814,   815,   799,  1405,   615,  1824,  1024,  1025,  1734,
     646,   647,  1137,  2025,  1336,  2026,   316,   542,   317,  1969,
      99,  1970,   195,  1360,    99,  2027,   517,   521,   133,   791,
     792,    99,  1269,   677,   301,   816,   318,    99,   418,   834,
     835,   423,  1670,   793,    99,   616,   750,   161,   276,    99,
     319,   320,   817,   836,  1117,   569,   839,   145,  1431,   173,
    1642,    99,   565,    99,  1648,  1260,   321,    99,  1664,  1337,
      99,   146,  1268,    99,   581,    99,    99,  1570,    99,   818,
     917,    99,   308,    99,   322,   432,    99,  1131,   323,   409,
     800,   543,   309,   522,  1971,   570,     4,    99,  1270,  1203,
     324,   801,   802,    99,  1876,   424,   325,   909,   910,   912,
      99,    99,   610,   137,   577,  1741,   687,  1138,  1972,    99,
    1204,    99,    99,  1735,   549,   129,    99,  1229,   296,   138,
     326,   327,   328,    99,  1670,   425,   550,   551,   552,   302,
      99,   329,    14,    99,    99,   759,   658,  1700,   330,  1230,
    1231,  1700,   553,   407,   923,   413,  1016,   415,   416,   417,
     419,   420,   421,   422,   113,   430,   669,   434,   435,   436,
      16,    99,   442,   124,   444,    99,   446,   447,  1269,   449,
     143,   451,   452,   453,   454,   455,   456,   457,   458,   459,
     460,   461,   462,   463,   464,   465,    99,    99,   470,  1101,
     472,   473,   474,   426,   477,   478,   479,   480,   481,   848,
     554,  1700,  1700,   486,  1630,   488,   297,  1021,  1757,   196,
    1825,   738,   491,  1806,   743,   130,   555,   968,   969,   691,
     961,  1261,   544,  1877,   571,   978,   331,  1632,  1633,   983,
    1338,   985,   131,  1042,  1270,   688,  1248,   998,   999,  1448,
     427,  1450,   546,  1000,  1001,  1973,   970,   971,   972,   973,
    1466,   433,  1468,   521,  1571,  1671,   332,   617,   158,   567,
     777,   568,  1361,   803,  1362,   104,   105,  1098,  1099,  1363,
    1373,   514,  1406,  1826,  1857,   514,  1719,   157,  1393,  1385,
    1043,  1201,   514,   794,  1756,  1049,   729,  1459,   514,   907,
    1364,   155,   819,   837,  1339,   514,   410,   411,  1093,   665,
     603,  1005,  1166,   159,   146,   333,  1118,  1365,  1128,   613,
    1603,   613,  1643,  2028,   514,  1836,  1649,   667,   611,   522,
    1631,   514,  1957,   160,   612,   735,   275,   514,   164,   514,
     649,   650,   651,   783,   652,   653,   654,   655,   310,   979,
     659,   556,  1340,  1167,   730,  1723,  1979,  1671,   514,   804,
     666,   668,   670,   773,   514,   528,   428,  1974,   636,   637,
    1804,   514,  1299,  1300,  1301,  1302,   334,   681,  1322,  1888,
     514,  1740,   850,  1327,   685,   689,  1488,   692,   901,   692,
     962,   963,   692,  1540,   692,  1846,   615,    99,   162,  1271,
    1542,  1739,   276,  1698,   529,   163,   712,  1701,  1272,  1273,
     197,  1776,  1249,   721,  1341,  1342,  1250,   165,   277,  1168,
     521,  1343,   198,   530,   736,   166,   199,   740,  1214,   641,
     642,  1950,   167,   200,  1586,  1952,   619,   168,   557,   620,
     621,   531,   532,  1283,   998,   999,  1269,   998,   999,   622,
    1000,  1001,   169,  1000,  1001,   797,  1959,  1977,   755,   756,
     125,   126,    16,   558,   278,  1295,   533,   276,  1161,   774,
     521,  1978,  1995,   775,   776,   171,   778,   779,   780,   781,
     974,   587,  1087,   277,   975,   784,   522,  1803,  1306,  1307,
    1308,  1912,  1794,   588,  1670,   172,   798,   998,   999,   534,
     174,  1162,   201,  1000,  1001,   979,   202,   279,   588,   974,
     134,   135,  1270,   975,  1670,  1936,  1109,   986,   987,   988,
     989,   990,   991,  1003,  1004,  1309,   589,  1291,   851,   278,
     185,  1294,   590,  1386,   807,   808,   522,   913,   591,   175,
    1956,   589,  1852,   203,   998,   999,   760,   590,   176,  1163,
    1000,  1164,   761,   591,   613,   613,   613,   437,   438,   439,
     613,   204,   979,   613,  1139,   360,   613,   613,  1292,   293,
     192,   193,   279,  1140,  1141,   186,   970,   971,   972,   973,
    1014,   146,   280,  1286,   177,  1019,  1287,  1288,  1289,  1290,
     281,   592,   467,   468,  1438,  1439,   593,  1293,  1442,  1443,
    1444,  1445,   574,   575,  1449,  2009,   592,  1142,  1296,  1297,
     179,   593,   290,  1298,   695,  2014,   180,   698,   181,   700,
     182,  1303,  1304,  1305,   372,  1054,  1509,  1510,   809,   183,
     810,  1058,   361,   362,   363,   364,   365,   366,   184,  1673,
    1202,  1674,  1673,    18,    19,    20,   187,   738,   738,  1213,
     188,  1215,   611,   636,   637,   281,   998,   999,   623,  1070,
     974,   189,  1000,  1001,   975,  1072,   190,  1284,  1525,   191,
    1076,  1077,   373,  1079,   811,   812,   970,   971,   972,   973,
     375,   974,   535,   376,  1084,   975,   716,   717,   718,   719,
     388,   205,  1311,  1312,  1313,  1314,  1315,  1316,  1317,  1096,
      99,   497,   498,   499,   500,   501,   502,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,  1671,   641,   642,
    1171,  1451,  1452,  1453,  1454,  1113,  1114,   377,  1115,   979,
    1172,   765,   766,   767,  1116,   389,   994,  1671,   995,   996,
     970,   971,   972,   973,  1173,  1395,    99,   497,   498,   499,
     500,   501,   502,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,  1398,   392,  1174,  1046,  1047,  1048,  1090,
    1091,  1092,   952,   953,  1124,   497,   498,   499,   500,   501,
     502,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,  1125,  1126,  1127,   391,  1135,   393,  1136,   955,   956,
    1353,  1144,  1145,  1051,  1347,  1148,   394,  1150,  1151,  1512,
    1513,  1514,  1154,  1155,  1156,  1157,  1158,  1159,  1160,  1175,
     396,  1170,   397,  1178,  1179,  1244,  1245,  1246,  1183,  1184,
    1462,  1186,  1187,  1497,  1498,  1499,   972,   973,  1383,  1600,
    1601,  1602,  1176,  1491,  1236,   998,   999,  1329,  1545,  1332,
    1548,  1207,  1005,  1704,  1705,   398,   399,   400,   401,   402,
     403,   404,  1219,  1220,  1221,  1222,   406,   405,  1224,  1225,
    1226,  1227,   441,   471,   475,   485,   494,  1234,   487,  1670,
     520,   537,   986,   987,   988,   989,   990,   991,   992,   993,
     538,   539,   540,   851,   562,   563,  1254,  1255,   576,  1257,
    1258,  1522,   577,   582,   584,   601,   605,   608,   627,  1523,
    1524,   648,  1526,   657,   613,   613,   613,   613,   660,   497,
     498,   499,   500,   501,   502,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   656,   613,   613,   661,  1460,
     662,   613,   613,   613,   613,   613,   663,   738,   664,   738,
     741,  1470,  1471,  1472,   671,   677,   672,   674,   675,  1319,
     513,   514,   676,   682,  1324,   497,   498,   499,   500,   501,
     502,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   680,   986,   987,   988,   989,   990,   991,   992,   993,
     683,   684,   690,   715,   694,  1346,   986,   987,   988,   989,
     990,   991,  1003,  1004,   697,   699,   579,   514,   702,   986,
     987,   988,   989,   990,   991,   992,   993,  1358,   703,   704,
     705,   706,  1370,   707,  1372,   708,   709,   710,   726,  1377,
     711,   727,   713,   714,   758,   986,   987,   988,   989,   990,
     991,  1003,  1004,   722,   725,   732,   733,   734,  1762,   737,
     744,   745,   746,   747,   748,  1389,  1390,   497,   498,   499,
     500,   501,   502,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   749,   750,   757,   764,   769,   770,   782,
     786,   787,   788,   789,   790,   378,   806,   821,  1624,   822,
    1410,   823,  1567,   970,   971,   972,   973,  1417,   824,   825,
     826,   994,  1671,   995,   996,   827,   379,   828,   829,   830,
     831,   832,   380,   833,   839,   381,   382,   842,   843,   846,
     852,   738,  1587,   907,   914,  1434,   949,  1056,   965,   383,
     497,   498,   499,   500,   501,   502,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,  1006,  1007,   384,  1008,
    1009,  1010,  1011,  1012,  1013,  1026,  1035,  1064,  1027,  1028,
    1029,  1030,  1031,  1474,  1475,  1476,  1477,  1059,  1480,  1481,
    1482,  1483,  1032,  1073,  1033,  1034,  1078,  1487,   847,  1489,
    1490,  1039,  1036,  1037,  1038,  1040,  1055,  1496,  1051,  1060,
    1080,  1061,  1501,   970,   971,   972,   973,  1062,  1085,  1063,
    1286,   994,  1066,   995,   996,  1508,  1067,  1069,   979,  1634,
     385,   386,  1095,  1111,  1292,   994,  1071,   995,   996,  1081,
     970,   971,   972,   973,  1461,  1112,  1082,  1088,   994,   613,
     995,   996,  1089,  1097,  1103,  1105,  1106,  1107,   613,  1108,
    1110,  1120,  1121,  1122,  1123,  1130,  1133,   979,  1134,  1147,
    1182,  1218,  1206,  1185,   994,  1188,   995,   996,  1189,  1216,
    1190,  1191,  1192,  1193,  1541,  1194,  1195,  1233,  1236,  1543,
    1240,  1241,  1239,  1544,  1252,  1242,  1547,  1243,  1253,  1256,
    1259,  1263,  1554,  1555,  1556,  1557,  1558,  1559,  1760,  1264,
    1561,  1761,  1265,  1318,  1321,  1323,  1328,  1326,  1329,  1350,
    1332,  1368,  1349,  1354,  1355,  1359,  1573,  1369,  1574,  1378,
    1380,  1381,  1575,  1576,  1387,  1407,  1623,  1388,  1391,  1394,
    1396,  1397,  1583,  1721,  1585,  1399,  1400,  1401,  1402,  1403,
    1404,  1408,  1409,    -8,     1,  1411,  1413,  1412,  1414,  1416,
    1415,    -8,  1418,  1419,  1420,  1431,  1421,  1422,  1423,  1424,
    1425,  1436,  1426,  1446,  1595,  1427,  1428,  1429,  1430,  1435,
    1456,  1599,  1457,  1440,  1441,    -8,  1458,  1478,  1484,  1502,
      -8,  1486,  1493,  1494,  1495,  1503,  1504,  1505,  1605,  1506,
    1507,  1516,  1517,  1519,  1518,  1528,  1521,   979,  1527,  1531,
    1532,  1550,  1533,  1534,    -8,  1535,  1536,    -8,  1537,  1538,
    1539,    -8,  1552,  1560,  1562,  1563,  1568,  1577,    -8,  1579,
      -8,  1581,  1593,  1627,    -8,  1628,  1637,  1650,  1588,  1653,
    1589,  1590,  1594,  1596,  1597,  1598,  1607,  1660,  1608,  1609,
    1610,  1611,  1612,  1613,  1651,  1615,  1616,  1617,  1663,  1618,
    1652,  1669,  1619,  1620,  1621,  1654,  1622,  1656,  1657,  1658,
    1659,  1625,  1626,  1635,  1638,  1678,  1639,  1640,  1641,  1644,
    1666,  1645,  1646,    -8,  1795,  1647,  1661,    -8,  1690,  1680,
    1662,  1667,  1675,  1676,  1692,  1668,  1683,  1687,  1689,  1693,
    1694,  1695,  1725,  1731,  1732,  1709,  1696,  1710,  1686,  1717,
      -8,  1688,  1724,  1737,  1728,    -8,    -8,  1729,  1691,  1730,
    1733,  1742,  1743,  1749,  1750,  1751,  1752,  1753,  1759,  1699,
    1780,  1768,  1703,  1783,  1707,  1708,  1758,  1786,  1721,  1769,
    1787,  1712,  1699,  1714,  1715,  1716,  1789,  1792,  1797,  1722,
    1793,  1808,  1819,  1800,  1853,  1822,  1833,  1727,  1802,  1805,
    1809,    -8,    -8,    -8,   916,  1818,   917,   918,   919,   920,
     921,  1831,  1820,  1821,  1856,  1851,   922,  1838,  1839,  1862,
    1842,  1848,  1858,  1859,    -8,  -734,  1871,  1873,  -735,  1880,
    1882,    -8,    -8,    -8,  1886,    -8,    -8,    -8,    -8,    -8,
    1890,    -8,  1897,  1892,  1893,  1898,    -8,    -8,    -8,  1899,
    1900,  1901,  1902,    -8,  1744,  1745,  1746,  1747,  1906,  1907,
    1913,  1917,  1918,    -8,  1754,    -8,    -8,  1920,  1923,  1935,
    1932,  1937,  1946,  1942,  1984,  1963,   613,  1966,   613,  1763,
     923,   924,   925,  1980,  1764,   926,   927,  1766,  1981,  1986,
    1989,  1993,  1999,    -8,    -8,    -8,  1772,  1773,  1774,  1775,
    2002,  2001,  2007,  2017,  2036,  1849,   483,    -8,    -8,  1843,
      -8,  1068,  2004,    -8,  2005,  1348,  2006,   928,  2011,  2012,
    2016,  2022,  1779,   929,  2039,  2040,  2041,  2043,  1546,  1549,
    1835,  1941,  1962,  1274,   763,   523,  1129,   754,  1798,  1799,
    1094,  1801,  1384,  1275,  1604,  1100,  1050,  1463,    -8,  1948,
    1464,  1810,  1811,  1813,  1814,  1815,  1816,  1817,  1132,  1492,
    1807,  1217,  1655,  1854,   371,  1530,  1023,  1345,  1765,     0,
    1767,  1018,   930,   931,   932,     0,  1828,     0,     0,     0,
       0,     0,     0,     0,  1837,     0,     0,  1841,     0,     0,
    1844,  1845,     0,  1847,   933,   934,     0,   935,     0,     0,
       0,  1855,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1860,  1861,     0,  1863,     0,     0,
       0,     0,     0,  1866,     0,     0,     0,     0,  1869,  1870,
       0,     0,     0,     0,     0,  1875,     0,  1879,     0,     0,
       0,     0,     0,  1881,     0,     0,   936,   937,     0,     0,
    1699,   938,   939,   940,   941,   942,   943,  1889,     0,     0,
    1891,     0,     0,     0,     0,  1894,  1895,  1896,     0,     0,
       0,     0,     0,     0,  1903,     0,     0,  1904,     0,     0,
    1905,     0,     0,     0,  1909,     0,  1910,     0,     0,     0,
    1916,     0,     0,     0,  1919,     0,     0,     0,  1921,  1922,
       0,     0,     0,     0,  1927,  1929,  1931,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1699,  1940,     0,     0,
       0,     0,     0,  1943,  1944,  1945,     0,     0,     0,     0,
       0,     0,  1951,     0,  1953,  1954,  1955,     0,     0,   207,
    1699,  1960,     0,     0,   208,   209,     0,     0,  1964,  1965,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1699,     0,     0,   210,     0,   211,  1983,   212,     0,  1985,
       0,     0,     0,  1987,     0,     0,  1990,  1991,  1992,     0,
       0,     0,     0,   213,  1996,  1997,     0,     0,     0,     0,
    2000,     0,     0,  2003,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2010,     0,     0,     0,     0,     0,     0,
       0,  2013,     0,  1699,  2015,     0,   214,     0,  1699,     0,
       0,  2018,  2019,  2020,  2021,     0,     0,     0,     0,  2030,
    2031,   215,  2033,     0,  2034,  2035,     0,  2037,  2038,     0,
     216,     0,  2042,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   217,     0,
       0,     0,     0,     0,   218,     0,     0,     0,     0,     0,
       0,     0,     0,   219,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   220,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   221,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   222,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   223,     0,   224,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   225,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     226,   227,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   238,   239,   240,   241,   242,   243,   244,
     245,   246,   247,   248,   249,   250,   251,     0,   252,     0,
     253,     0,     0,     0,     0,     0,     0,     0,   254,   255,
     256,     0,     0,     5,     0,     0,     0,   257,   258,   259,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   260,
     261,   262,   263,   264,     0,     0,     0,     6,     0,     0,
     853,     0,   265,     0,     0,   266,     0,     0,     0,     0,
       0,     0,     0,   854,   267,   268,   269,   270,   271,   855,
       0,     0,     0,   856,     0,     0,     7,     0,     0,     8,
     857,   858,     0,     9,   276,     0,     0,     0,     0,     0,
      10,     0,    11,   859,     0,     0,    12,     0,     0,   860,
       0,   861,   862,     0,     0,     0,   863,   864,     0,     0,
       0,     0,   865,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   866,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     867,     0,   868,     0,   869,    13,   870,     0,     0,    14,
       0,     0,   871,   872,     0,     0,   873,   874,     0,     0,
       0,     0,   875,   876,     0,     0,     0,     0,     0,     0,
       0,     0,    15,     0,     0,     0,     0,    16,    17,     0,
       0,     0,   877,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   878,     0,     0,   879,     0,     0,
     880,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    18,    19,    20,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    21,     0,     0,     0,
       0,     0,     0,    22,    23,    24,     0,    25,    26,    27,
      28,    29,     0,    30,     0,     0,     0,     0,    31,    32,
      33,     0,   881,     0,     0,    34,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    35,   882,    36,    37,     0,
       0,     0,     0,     0,     0,   883,     0,   884,   885,   886,
     887,   888,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    38,    39,    40,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    41,
      42,     0,    43,     0,     0,    44,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   889,   890,   891,   892,   893,   894,   895,   896,     0,
      45,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   897,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   898,     0,     0,     0,     0,     0,     0,     0,   899,
       0,     0,     0,     0,     0,   900
};

static const yytype_int16 yycheck[] =
{
       6,    21,   626,   147,    10,   153,   525,  1209,   598,   599,
     908,   909,  1564,   335,   604,    34,    62,    16,    62,    25,
      26,    27,    28,    29,    30,    31,    32,    33,  1371,   618,
    1211,    29,    62,    54,   610,    41,   625,    62,    37,    45,
     616,   375,   136,   619,   376,   125,   622,   623,    14,   140,
     172,    55,    56,    46,   154,   171,    21,   646,   647,    20,
     392,   393,     6,    86,     9,    88,    16,    22,    18,    35,
     172,    37,   235,   138,   172,    98,   287,   170,   432,   168,
     169,   172,   170,   170,    37,    89,    36,   172,   124,   168,
     169,    28,   213,   182,   172,   211,   170,   103,    48,   172,
      50,    51,   106,   182,   170,    10,   170,    50,   170,   115,
     170,   172,   323,   172,   170,    37,    66,   172,  1515,    64,
     172,    64,    64,   172,   335,   172,   172,   167,   172,   133,
      39,   172,   131,   172,    84,   257,   172,    37,    88,   160,
     133,    96,   141,   236,   110,    50,     0,   172,   236,    16,
     100,   144,   145,   172,   167,    92,   106,   590,   591,   592,
     172,   172,   211,    84,    64,  1634,   167,   111,   134,   172,
      37,   172,   172,   134,    10,   172,   172,    60,    37,   100,
     130,   131,   132,   172,   213,   122,    22,    23,    24,   142,
     172,   141,   134,   172,   172,   517,   257,   172,   148,    82,
      83,   172,    38,   209,   113,   211,    37,   213,   214,   215,
     216,   217,   218,   219,   171,   221,   257,   223,   224,   225,
     162,   172,   228,   171,   230,   172,   232,   233,   170,   235,
     170,   237,   238,   239,   240,   241,   242,   243,   244,   245,
     246,   247,   248,   249,   250,   251,   172,   172,   254,   768,
     256,   257,   258,   190,   260,   261,   262,   263,   264,   581,
      96,   172,   172,   269,  1466,   271,   125,    37,  1665,   432,
     235,   482,   278,  1742,   485,   172,   112,   611,   612,   257,
      37,   203,   237,   296,   189,   617,   236,  1468,  1469,   623,
     235,   623,   172,   373,   236,   296,    59,   216,   217,  1193,
     237,  1195,   308,   222,   223,   271,   427,   428,   429,   430,
    1208,   433,  1210,   170,   354,   436,   266,   433,   432,   325,
     531,   327,   387,   316,   389,    93,    94,   760,   761,   394,
     432,   433,   432,   298,  1803,   433,   434,   170,   432,   432,
     420,   432,   433,   432,   432,   432,   366,   432,   433,   379,
     415,    37,   356,   432,   299,   433,   377,   378,   432,   432,
     366,   985,     5,   170,    64,   315,   432,   432,   432,   375,
     432,   377,   432,   396,   433,   434,   432,   432,   427,   236,
     432,   433,  1934,   170,   433,   432,   384,   433,   170,   433,
     396,   397,   398,   432,   400,   401,   402,   403,   397,   428,
     406,   237,   347,    46,   424,  1586,  1958,   436,   433,   402,
     416,   417,   418,    54,   433,    33,   353,   383,   249,   250,
     432,   433,   998,   999,  1000,  1001,   376,   433,  1017,   432,
     433,  1633,   432,  1022,   440,   441,   432,   443,   586,   445,
     197,   198,   448,   432,   450,  1788,   171,   172,   432,   391,
     432,  1632,    48,   432,    72,   432,   462,   432,   400,   401,
      14,   432,   225,   469,   409,   410,   229,   170,    64,   112,
     170,   416,    26,    91,   480,   170,    30,   483,   911,   249,
     250,   432,   170,    37,  1382,   432,   211,   170,   324,   214,
     215,   109,   110,   212,   216,   217,   170,   216,   217,   224,
     222,   223,   170,   222,   223,   357,   432,   432,   514,   515,
      93,    94,   162,   349,   110,   434,   134,    48,    56,   160,
     170,   432,   432,   529,   530,   432,   532,   533,   534,   535,
     432,    37,   743,    64,   436,   541,   236,  1739,   232,   233,
     234,  1884,  1723,    64,   213,   170,   398,   216,   217,   167,
     432,    89,   106,   222,   223,   428,   110,   153,    64,   432,
     318,   319,   236,   436,   213,  1908,   777,   216,   217,   218,
     219,   220,   221,   222,   223,   269,    97,   977,   584,   110,
     170,   981,   103,  1102,    11,    12,   236,   593,   109,   432,
    1933,    97,  1794,   147,   216,   217,   103,   103,   432,   137,
     222,   139,   109,   109,   610,   611,   612,   292,   293,   294,
     616,   165,   428,   619,    46,   130,   622,   623,   434,   763,
      93,    94,   153,    55,    56,   170,   427,   428,   429,   430,
     636,    64,   163,   434,   432,   641,   970,   971,   972,   973,
     236,   162,   371,   372,  1184,  1185,   167,   979,  1188,  1189,
    1190,  1191,   267,   268,  1194,  1998,   162,    89,   992,   993,
     432,   167,    29,   997,   445,  2008,   432,   448,   432,   450,
     432,  1003,  1004,  1005,    67,   681,  1266,  1267,   105,   432,
     107,   687,   197,   198,   199,   200,   201,   202,   432,  1523,
     901,  1525,  1526,   208,   209,   210,   432,   908,   909,   910,
     432,   912,   427,   249,   250,   236,   216,   217,   433,   715,
     432,   432,   222,   223,   436,   721,   432,   436,  1294,   432,
     726,   727,   432,   729,   151,   152,   427,   428,   429,   430,
     435,   432,   350,   435,   740,   436,   362,   363,   364,   365,
     432,   295,  1007,  1008,  1009,  1010,  1011,  1012,  1013,   755,
     172,   173,   174,   175,   176,   177,   178,   179,   180,   181,
     182,   183,   184,   185,   186,   187,   188,   436,   249,   250,
       7,   362,   363,   364,   365,   781,   782,   435,   784,   428,
      17,   170,   171,   172,   790,   432,   435,   436,   437,   438,
     427,   428,   429,   430,    31,  1117,   172,   173,   174,   175,
     176,   177,   178,   179,   180,   181,   182,   183,   184,   185,
     186,   187,   188,  1135,   435,    52,   170,   171,   172,   170,
     171,   172,   168,   169,   830,   173,   174,   175,   176,   177,
     178,   179,   180,   181,   182,   183,   184,   185,   186,   187,
     188,   170,   171,   172,   432,   851,   435,   853,   168,   169,
    1061,   857,   858,   433,   434,   861,   432,   863,   864,   170,
     171,   172,   868,   869,   870,   871,   872,   873,   874,   106,
      99,   877,    73,   879,   880,   281,   282,   283,   884,   885,
    1202,   887,   888,   281,   282,   283,   429,   430,  1099,   170,
     171,   172,   129,   432,   433,   216,   217,   433,   434,   433,
     434,   907,  1526,   351,   352,   390,   157,    74,    85,    78,
     166,   297,   918,   919,   920,   921,    25,   172,   924,   925,
     926,   927,   292,   172,   418,   433,    37,   933,   432,   213,
     170,   134,   216,   217,   218,   219,   220,   221,   222,   223,
     170,    61,   170,   949,   170,   170,   952,   953,   432,   955,
     956,  1285,    64,    37,   170,   170,    37,   386,   231,  1291,
    1294,   170,  1294,   432,   970,   971,   972,   973,   432,   173,
     174,   175,   176,   177,   178,   179,   180,   181,   182,   183,
     184,   185,   186,   187,   188,   172,   992,   993,   432,  1200,
     432,   997,   998,   999,  1000,  1001,   170,  1208,   432,  1210,
     414,  1212,  1213,  1214,   432,   170,   432,   432,   432,  1015,
     432,   433,   432,   432,  1020,   173,   174,   175,   176,   177,
     178,   179,   180,   181,   182,   183,   184,   185,   186,   187,
     188,   433,   216,   217,   218,   219,   220,   221,   222,   223,
     432,   432,   432,   167,   432,  1051,   216,   217,   218,   219,
     220,   221,   222,   223,   432,   432,   432,   433,   432,   216,
     217,   218,   219,   220,   221,   222,   223,  1073,   432,   432,
     432,   432,  1078,   432,  1080,   432,   432,   432,   354,  1085,
     432,   370,   432,   432,   432,   216,   217,   218,   219,   220,
     221,   222,   223,   432,   432,   432,   432,   432,  1674,   432,
     432,   432,   170,   170,   432,  1111,  1112,   173,   174,   175,
     176,   177,   178,   179,   180,   181,   182,   183,   184,   185,
     186,   187,   188,   170,   170,   432,   432,   170,    37,   170,
     432,   432,   432,   170,    13,    37,   349,   432,  1460,   170,
    1146,   170,  1353,   427,   428,   429,   430,  1153,   170,   432,
     170,   435,   436,   437,   438,   432,    58,   170,   432,    13,
     432,   432,    64,   432,   170,    67,    68,   432,   432,   432,
     432,  1382,  1383,   379,    37,  1181,   170,   295,   432,    81,
     173,   174,   175,   176,   177,   178,   179,   180,   181,   182,
     183,   184,   185,   186,   187,   188,   170,   170,   100,   170,
     170,   170,   170,   170,   170,   170,   433,   313,   432,   432,
     432,   432,   432,  1219,  1220,  1221,  1222,   172,  1224,  1225,
    1226,  1227,   432,   167,   432,   432,   167,  1233,   432,  1235,
    1236,   433,   432,   432,   432,   432,   432,  1243,   433,   432,
     167,   433,  1248,   427,   428,   429,   430,   432,   167,   432,
     434,   435,   432,   437,   438,  1261,   432,   432,   428,  1470,
     162,   163,   170,    13,   434,   435,   432,   437,   438,   432,
     427,   428,   429,   430,   432,   149,   432,   432,   435,  1285,
     437,   438,   432,   432,   432,   432,   432,   432,  1294,   432,
     432,   432,   432,   432,   432,   170,   432,   428,   170,   170,
     170,   148,   170,   172,   435,   172,   437,   438,   172,   432,
     172,   172,   172,   172,  1320,   172,   172,   206,   433,  1325,
     171,   171,   432,  1329,   432,   171,  1332,   226,   432,   432,
     432,   432,  1338,  1339,  1340,  1341,  1342,  1343,  1672,   170,
    1346,  1673,   170,   432,   244,   432,   432,   248,   433,   172,
     433,   167,   434,   432,   432,   432,  1362,   134,  1364,   375,
     374,   434,  1368,  1369,   170,   170,   432,   432,   432,   432,
     432,   432,  1378,  1584,  1380,   432,   432,   432,   432,   432,
     432,   432,   432,     0,     1,   170,   170,   432,   432,   170,
     432,     8,   432,   432,   432,   170,   432,   432,   432,   432,
     432,    64,   432,    64,  1410,   432,   432,   432,   432,   432,
     170,  1417,   171,   432,   432,    32,   170,   170,   170,   156,
      37,   432,   432,   432,   432,    34,    34,    34,  1434,    34,
     202,   170,   170,    37,   170,   235,   212,   428,   171,   170,
     432,   421,   432,   432,    61,   432,   432,    64,   432,   432,
     432,    68,   172,   432,   432,   432,   370,   370,    75,   422,
      77,   167,    34,    37,    81,   413,   170,   207,   432,   226,
     432,   432,   432,   432,   432,   432,   432,   204,   432,   432,
     432,   432,   432,   432,  1490,   432,   432,   432,   170,   432,
    1496,   170,   432,   432,   432,  1501,   432,  1503,  1504,  1505,
    1506,   432,   432,   432,   432,   172,   432,   432,   432,   432,
    1516,   432,   432,   130,  1725,   432,   432,   134,   170,   251,
     432,   432,  1528,  1529,   235,   432,   251,   432,   432,   370,
     370,   417,   109,   170,   170,   432,   434,   432,  1544,   432,
     157,  1547,   432,    34,   432,   162,   163,   432,  1554,   432,
     432,   170,   432,   226,   138,   138,   138,   138,   213,  1565,
     354,   434,  1568,   370,  1570,  1571,   432,   134,  1779,   434,
     167,  1577,  1578,  1579,  1580,  1581,   403,   432,    13,  1585,
     432,   434,   167,   432,  1795,   136,   228,  1593,   432,   432,
     432,   208,   209,   210,    37,   432,    39,    40,    41,    42,
      43,   416,   432,   432,    13,   388,    49,   432,   411,   159,
     432,   432,   432,   432,   231,   432,   412,   419,   432,   370,
     134,   238,   239,   240,   134,   242,   243,   244,   245,   246,
     138,   248,    64,   432,   432,    64,   253,   254,   255,   432,
     252,   247,   235,   260,  1650,  1651,  1652,  1653,   172,   228,
     167,   170,   432,   270,  1660,   272,   273,   138,   432,   172,
     370,   167,   432,   170,   432,   170,  1672,   170,  1674,  1675,
     113,   114,   115,   370,  1680,   118,   119,  1683,   167,    25,
     124,   423,   375,   300,   301,   302,  1692,  1693,  1694,  1695,
     170,   432,   134,   350,   170,  1791,   266,   314,   315,  1785,
     317,   701,   432,   320,   432,  1052,   432,   150,   432,   432,
     432,   432,  1702,   156,   432,   432,   432,   432,  1330,  1333,
    1778,  1915,  1940,   966,   519,   294,   840,   496,  1734,  1735,
     751,  1737,  1100,   966,  1432,   763,   678,  1203,   355,  1925,
    1205,  1747,  1748,  1749,  1750,  1751,  1752,  1753,   845,  1237,
    1744,   915,  1502,  1796,   154,  1308,   643,  1045,  1681,    -1,
    1684,   638,   205,   206,   207,    -1,  1772,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1780,    -1,    -1,  1783,    -1,    -1,
    1786,  1787,    -1,  1789,   227,   228,    -1,   230,    -1,    -1,
      -1,  1797,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1810,  1811,    -1,  1813,    -1,    -1,
      -1,    -1,    -1,  1819,    -1,    -1,    -1,    -1,  1824,  1825,
      -1,    -1,    -1,    -1,    -1,  1831,    -1,  1833,    -1,    -1,
      -1,    -1,    -1,  1839,    -1,    -1,   279,   280,    -1,    -1,
    1846,   284,   285,   286,   287,   288,   289,  1853,    -1,    -1,
    1856,    -1,    -1,    -1,    -1,  1861,  1862,  1863,    -1,    -1,
      -1,    -1,    -1,    -1,  1870,    -1,    -1,  1873,    -1,    -1,
    1876,    -1,    -1,    -1,  1880,    -1,  1882,    -1,    -1,    -1,
    1886,    -1,    -1,    -1,  1890,    -1,    -1,    -1,  1894,  1895,
      -1,    -1,    -1,    -1,  1900,  1901,  1902,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1912,  1913,    -1,    -1,
      -1,    -1,    -1,  1919,  1920,  1921,    -1,    -1,    -1,    -1,
      -1,    -1,  1928,    -1,  1930,  1931,  1932,    -1,    -1,     9,
    1936,  1937,    -1,    -1,    14,    15,    -1,    -1,  1944,  1945,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1956,    -1,    -1,    33,    -1,    35,  1962,    37,    -1,  1965,
      -1,    -1,    -1,  1969,    -1,    -1,  1972,  1973,  1974,    -1,
      -1,    -1,    -1,    53,  1980,  1981,    -1,    -1,    -1,    -1,
    1986,    -1,    -1,  1989,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1999,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2007,    -1,  2009,  2010,    -1,    86,    -1,  2014,    -1,
      -1,  2017,  2018,  2019,  2020,    -1,    -1,    -1,    -1,  2025,
    2026,   101,  2028,    -1,  2030,  2031,    -1,  2033,  2034,    -1,
     110,    -1,  2038,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   128,    -1,
      -1,    -1,    -1,    -1,   134,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   143,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   155,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   167,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   236,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   256,    -1,   258,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   271,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     290,   291,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   303,   304,   305,   306,   307,   308,   309,
     310,   311,   312,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   333,   334,   335,   336,   337,   338,   339,
     340,   341,   342,   343,   344,   345,   346,    -1,   348,    -1,
     350,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   358,   359,
     360,    -1,    -1,     8,    -1,    -1,    -1,   367,   368,   369,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   379,
     380,   381,   382,   383,    -1,    -1,    -1,    32,    -1,    -1,
      14,    -1,   392,    -1,    -1,   395,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    27,   404,   405,   406,   407,   408,    33,
      -1,    -1,    -1,    37,    -1,    -1,    61,    -1,    -1,    64,
      44,    45,    -1,    68,    48,    -1,    -1,    -1,    -1,    -1,
      75,    -1,    77,    57,    -1,    -1,    81,    -1,    -1,    63,
      -1,    65,    66,    -1,    -1,    -1,    70,    71,    -1,    -1,
      -1,    -1,    76,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    90,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     104,    -1,   106,    -1,   108,   130,   110,    -1,    -1,   134,
      -1,    -1,   116,   117,    -1,    -1,   120,   121,    -1,    -1,
      -1,    -1,   126,   127,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   157,    -1,    -1,    -1,    -1,   162,   163,    -1,
      -1,    -1,   146,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   158,    -1,    -1,   161,    -1,    -1,
     164,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   208,   209,   210,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   231,    -1,    -1,    -1,
      -1,    -1,    -1,   238,   239,   240,    -1,   242,   243,   244,
     245,   246,    -1,   248,    -1,    -1,    -1,    -1,   253,   254,
     255,    -1,   236,    -1,    -1,   260,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   270,   250,   272,   273,    -1,
      -1,    -1,    -1,    -1,    -1,   259,    -1,   261,   262,   263,
     264,   265,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   300,   301,   302,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   314,
     315,    -1,   317,    -1,    -1,   320,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   325,   326,   327,   328,   329,   330,   331,   332,    -1,
     355,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   348,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   385,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   393,
      -1,    -1,    -1,    -1,    -1,   399
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_int16 yystos[] =
{
       0,     1,   440,   446,     0,     8,    32,    61,    64,    68,
      75,    77,    81,   130,   134,   157,   162,   163,   208,   209,
     210,   231,   238,   239,   240,   242,   243,   244,   245,   246,
     248,   253,   254,   255,   260,   270,   272,   273,   300,   301,
     302,   314,   315,   317,   320,   355,   441,   444,   445,   448,
     449,   450,   451,   452,   453,   454,   458,   459,   462,   463,
     558,   560,   562,   563,   590,   592,   593,   611,   612,   618,
     619,   626,   627,   628,   649,   650,   666,   668,   781,   783,
     801,   803,   805,   816,   852,   853,   854,   855,   856,   868,
     880,   881,   882,   883,   884,   885,   886,   887,   784,   172,
     443,   464,   669,   443,    93,    94,   629,   651,   591,   806,
     807,   808,   817,   171,   442,   443,   443,   443,   443,   443,
     443,   443,   443,   443,   171,    93,    94,   802,   804,   172,
     172,   172,   443,   432,   318,   319,   455,    84,   100,   456,
     443,   460,   467,   170,   561,    50,    64,   596,   602,   605,
     614,   621,   654,   672,   787,    37,   447,   170,   432,   170,
     170,   443,   432,   432,   170,   170,   170,   170,   170,   170,
     818,   432,   170,   443,   432,   432,   432,   432,   857,   432,
     432,   432,   432,   432,   432,   170,   170,   432,   432,   432,
     432,   432,    93,    94,   457,   235,   432,    14,    26,    30,
      37,   106,   110,   147,   165,   295,   461,     9,    14,    15,
      33,    35,    37,    53,    86,   101,   110,   128,   134,   143,
     155,   167,   236,   256,   258,   271,   290,   291,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   348,   350,   358,   359,   360,   367,   368,   369,
     379,   380,   381,   382,   383,   392,   395,   404,   405,   406,
     407,   408,   465,   468,    29,   384,    48,    64,   110,   153,
     163,   236,   564,   571,   572,   574,   578,   579,   582,   583,
      29,   595,   606,   602,   603,   607,    37,   125,   613,   615,
     616,    37,   142,   620,   622,   625,    16,    37,   131,   141,
     397,   652,   655,   656,   657,   660,    16,    18,    36,    50,
      51,    66,    84,    88,   100,   106,   130,   131,   132,   141,
     148,   236,   266,   315,   376,   579,   667,   673,   676,   680,
     685,   686,   687,   688,   689,   690,   691,   692,   694,   696,
     697,   699,   700,   701,   749,   750,   751,   758,   760,   761,
     130,   197,   198,   199,   200,   201,   202,   698,   782,   788,
     794,   805,    67,   432,   630,   435,   435,   435,    37,    58,
      64,    67,    68,    81,   100,   162,   163,   819,   432,   432,
     869,   432,   435,   435,   432,   559,    99,    73,   390,   157,
      74,    85,    78,   166,   297,   172,    25,   443,    54,   160,
     377,   378,   546,   443,   466,   443,   443,   443,   124,   443,
     443,   443,   443,    28,    92,   122,   190,   237,   353,   545,
     443,   473,   257,   433,   443,   443,   443,   292,   293,   294,
     525,   292,   443,   479,   443,   480,   443,   443,   482,   443,
     483,   443,   443,   443,   443,   443,   443,   443,   443,   443,
     443,   443,   443,   443,   443,   443,   485,   371,   372,   509,
     443,   172,   443,   443,   443,   418,   491,   443,   443,   443,
     443,   443,   494,   509,   469,   433,   443,   432,   443,   580,
     584,   443,   565,   575,    37,   588,   573,   173,   174,   175,
     176,   177,   178,   179,   180,   181,   182,   183,   184,   185,
     186,   187,   188,   432,   433,   443,   581,   665,   585,   594,
     170,   170,   236,   562,   597,   598,   599,   604,    33,    72,
      91,   109,   110,   134,   167,   350,   608,   134,   170,    61,
     170,   623,    22,    96,   237,   653,   443,   658,   661,    10,
      22,    23,    24,    38,    96,   112,   237,   324,   349,   681,
     759,   693,   170,   170,   695,   665,   702,   443,   443,    10,
      50,   189,   677,   674,   267,   268,   432,    64,   753,   432,
     581,   665,    37,   670,   170,   740,   705,    37,    64,    97,
     103,   109,   162,   167,   726,   727,   736,   763,   790,   791,
     792,   170,   793,   443,   789,    37,   785,   795,   386,   633,
     211,   427,   433,   443,   812,   171,   211,   433,   814,   211,
     214,   215,   224,   433,   812,   813,   814,   231,   821,   826,
     820,   823,   827,   822,   824,   825,   249,   250,   871,   872,
     873,   249,   250,   859,   860,   861,   814,   814,   170,   443,
     443,   443,   443,   443,   443,   443,   172,   432,   257,   443,
     432,   432,   432,   170,   432,   432,   443,   432,   443,   257,
     443,   432,   432,   471,   432,   432,   432,   170,   537,   538,
     433,   443,   432,   432,   432,   443,   474,   167,   296,   443,
     432,   257,   443,   520,   432,   520,   481,   432,   520,   432,
     520,   484,   432,   432,   432,   432,   432,   432,   432,   432,
     432,   432,   443,   432,   432,   167,   362,   363,   364,   365,
     548,   443,   432,   488,   490,   432,   354,   370,    21,   366,
     424,   499,   432,   432,   432,   432,   443,   432,   665,   733,
     443,   414,   495,   665,   432,   432,   170,   170,   432,   170,
     170,   576,   577,   589,   574,   443,   443,   432,   432,   581,
     103,   109,   586,   596,   432,   170,   171,   172,   600,   170,
      37,   609,   599,    54,   160,   443,   443,   665,   443,   443,
     443,   443,   170,   432,   443,   624,   432,   432,   432,   170,
      13,   168,   169,   182,   432,   659,   662,   357,   398,    46,
     133,   144,   145,   316,   402,   683,   349,    11,    12,   105,
     107,   151,   152,   684,    55,    56,    89,   106,   133,   356,
     682,   432,   170,   170,   170,   432,   170,   432,   170,   432,
      13,   432,   432,   432,   168,   169,   182,   432,   678,   170,
     675,   679,   432,   432,   754,   752,   432,   432,   581,   671,
     432,   443,   432,    14,    27,    33,    37,    44,    45,    57,
      63,    65,    66,    70,    71,    76,    90,   104,   106,   108,
     110,   116,   117,   120,   121,   126,   127,   146,   158,   161,
     164,   236,   250,   259,   261,   262,   263,   264,   265,   325,
     326,   327,   328,   329,   330,   331,   332,   348,   385,   393,
     399,   579,   703,   706,   720,   721,   728,   379,   617,   617,
     617,    62,   617,   443,    37,   730,    37,    39,    40,    41,
      42,    43,    49,   113,   114,   115,   118,   119,   150,   156,
     205,   206,   207,   227,   228,   230,   279,   280,   284,   285,
     286,   287,   288,   289,   762,   764,   768,   769,   777,   170,
     740,   740,   168,   169,   741,   168,   169,   748,   799,   740,
     786,    37,   197,   198,   796,   432,   634,   813,   812,   812,
     427,   428,   429,   430,   432,   436,   809,   813,   814,   428,
     809,   813,   813,   812,   813,   814,   216,   217,   218,   219,
     220,   221,   222,   223,   435,   437,   438,   815,   216,   217,
     222,   223,   809,   222,   223,   815,   170,   170,   170,   170,
     170,   170,   170,   170,   443,   874,    37,   870,   872,   443,
     862,    37,   858,   860,   809,   809,   170,   432,   432,   432,
     432,   432,   432,   432,   432,   433,   432,   432,   432,   433,
     432,   125,   373,   420,   837,   838,   170,   171,   172,   432,
     538,   433,   539,   540,   443,   432,   295,   526,   443,   172,
     432,   433,   432,   432,   313,   524,   432,   432,   524,   432,
     443,   432,   443,   167,   547,   517,   443,   443,   167,   443,
     167,   432,   432,   735,   443,   167,   496,   665,   432,   432,
     170,   171,   172,   432,   577,   170,   443,   432,   617,   617,
     597,   599,   601,   432,   610,   432,   432,   432,   432,   665,
     432,    13,   149,   443,   443,   443,   443,   170,   432,   663,
     432,   432,   432,   432,   443,   170,   171,   172,   432,   679,
     170,    37,   753,   432,   170,   443,   443,     6,   111,    46,
      55,    56,    89,   704,   443,   443,   711,   170,   443,   707,
     443,   443,   708,   709,   443,   443,   443,   443,   443,   443,
     443,    56,    89,   137,   139,   724,     5,    46,   112,   725,
     443,     7,    17,    31,    52,   106,   129,   723,   443,   443,
     712,   710,   170,   443,   443,   172,   443,   443,   172,   172,
     172,   172,   172,   172,   172,   172,   713,   716,   714,   715,
     140,   432,   665,    16,    37,   722,   170,   443,    62,   733,
      62,   733,    62,   665,   617,   665,   432,   727,   148,   443,
     443,   443,   443,   765,   443,   443,   443,   443,   766,    60,
      82,    83,   778,   206,   443,   775,   433,   773,   774,   432,
     171,   171,   171,   226,   281,   282,   283,   770,    59,   225,
     229,   772,   432,   432,   443,   443,   432,   443,   443,   432,
      37,   203,   800,   432,   170,   170,   797,   798,    64,   170,
     236,   391,   400,   401,   560,   611,   631,   635,   636,   637,
     638,   639,   642,   212,   436,   810,   434,   812,   812,   812,
     812,   810,   434,   814,   810,   434,   812,   812,   812,   813,
     813,   813,   813,   814,   814,   814,   232,   233,   234,   269,
     828,   828,   828,   828,   828,   828,   828,   828,   432,   443,
     875,   244,   809,   432,   443,   863,   248,   809,   432,   433,
     541,   542,   433,   543,   544,   839,     9,    64,   235,   299,
     347,   409,   410,   416,   840,   837,   443,   434,   540,   434,
     172,   477,   475,   665,   432,   432,   489,   492,   443,   432,
     138,   387,   389,   394,   415,   432,   518,   519,   167,   134,
     443,   500,   443,   432,   665,   734,   493,   443,   375,   486,
     374,   434,   587,   665,   609,   432,   599,   170,   432,   443,
     443,   432,   136,   432,   432,   581,   432,   432,   581,   432,
     432,   432,   432,   432,   432,   154,   432,   170,   432,   432,
     443,   170,   432,   170,   432,   432,   170,   443,   432,   432,
     432,   432,   432,   432,   432,   432,   432,   432,   432,   432,
     432,   170,   718,   719,   443,   432,    64,   848,   848,   848,
     432,   432,   848,   848,   848,   848,    64,   850,   850,   848,
     850,   362,   363,   364,   365,   717,   170,   171,   170,   432,
     665,   432,   581,   681,   726,   729,   733,   735,   733,   734,
     665,   665,   665,   737,   443,   443,   443,   443,   170,   779,
     443,   443,   443,   443,   170,   780,   432,   443,   432,   443,
     443,   432,   774,   432,   432,   432,   443,   281,   282,   283,
     771,   443,   156,    34,    34,    34,    34,   202,   443,   740,
     740,   643,   170,   171,   172,   640,   170,   170,   170,    37,
     632,   212,   812,   814,   812,   813,   814,   171,   235,   834,
     834,   170,   432,   432,   432,   432,   432,   432,   432,   432,
     432,   443,   432,   443,   443,   434,   542,   443,   434,   544,
     421,   833,   172,   841,   443,   443,   443,   443,   443,   443,
     432,   443,   432,   432,   527,   535,   536,   665,   370,   514,
     167,   354,   510,   443,   443,   443,   443,   370,   535,   422,
     555,   167,   513,   443,   522,   443,   733,   665,   432,   432,
     432,   664,   755,    34,   432,   443,   432,   432,   432,   443,
     170,   171,   172,   432,   719,   443,   849,   432,   432,   432,
     432,   432,   432,   432,   851,   432,   432,   432,   432,   432,
     432,   432,   432,   432,   581,   432,   432,    37,   413,   731,
     735,   432,   734,   734,   665,   432,   738,   170,   432,   432,
     432,   432,   170,   432,   432,   432,   432,   432,   170,   432,
     207,   443,   443,   226,   443,   778,   443,   443,   443,   443,
     204,   432,   432,   170,   639,   641,   443,   432,   432,   170,
     213,   436,   811,   811,   811,   443,   443,   836,   172,   835,
     251,   876,   877,   251,   864,   865,   443,   432,   443,   432,
     170,   443,   235,   370,   370,   417,   434,   536,   432,   443,
     172,   432,   521,   443,   351,   352,   515,   443,   443,   432,
     432,   472,   443,   501,   443,   443,   443,   432,   470,   434,
     523,   665,   443,   734,   432,   109,   757,   443,   432,   432,
     432,   170,   170,   432,    20,   134,   732,    34,   739,   734,
     735,   739,   170,   432,   443,   443,   443,   443,   767,   226,
     138,   138,   138,   138,   443,   644,   432,   639,   432,   213,
     812,   814,   813,   443,   443,   877,   443,   865,   434,   434,
     842,   844,   443,   443,   443,   443,   432,   478,   476,   522,
     354,   516,   511,   370,   507,   508,   134,   167,   504,   403,
     497,   498,   432,   432,   734,   665,   756,    13,   443,   443,
     432,   443,   432,   735,   432,   432,   739,   777,   434,   432,
     443,   443,   776,   443,   443,   443,   443,   443,   432,   167,
     432,   432,   136,   846,    21,   235,   298,   829,   443,   843,
     845,   416,   528,   228,   533,   533,   434,   443,   432,   411,
     512,   443,   432,   507,   443,   443,   535,   443,   432,   497,
     487,   388,   735,   665,   757,   443,    13,   739,   432,   432,
     443,   443,   159,   443,   742,   744,   443,   878,   866,   443,
     443,   412,   847,   419,   830,   443,   167,   296,   530,   443,
     370,   443,   134,   506,   502,   505,   134,   566,   432,   443,
     138,   443,   432,   432,   443,   443,   443,    64,    64,   432,
     252,   247,   235,   443,   443,   443,   172,   228,   534,   443,
     443,   556,   535,   167,   552,   553,   443,   170,   432,   443,
     138,   443,   443,   432,   743,   745,   645,   443,   879,   443,
     867,   443,   370,   532,   531,   172,   535,   167,   557,   503,
     443,   552,   170,   443,   443,   443,   432,   746,   746,   647,
     432,   443,   432,   443,   443,   443,   535,   536,   529,   432,
     443,   549,   555,   170,   443,   443,   170,   747,    14,    35,
      37,   110,   134,   271,   383,   648,   831,   432,   432,   536,
     370,   167,   550,   443,   432,   443,    25,   443,   646,   124,
     443,   443,   443,   423,   832,   432,   443,   443,   554,   375,
     443,   432,   170,   443,   432,   432,   432,   134,   551,   535,
     443,   432,   432,   443,   535,   443,   432,   350,   443,   443,
     443,   443,   432,   567,   568,    86,    88,    98,   396,   569,
     443,   443,   570,   443,   443,   443,   170,   443,   443,   432,
     432,   432,   443,   432
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int16 yyr1[] =
{
       0,   439,   440,   442,   441,   443,   444,   445,   446,   446,
     446,   447,   447,   448,   448,   448,   448,   448,   448,   448,
     448,   448,   448,   448,   448,   448,   448,   448,   448,   448,
     448,   448,   448,   448,   448,   448,   448,   448,   448,   448,
     448,   448,   448,   448,   448,   448,   448,   448,   448,   448,
     448,   448,   449,   449,   450,   450,   451,   452,   453,   454,
     455,   455,   456,   456,   457,   457,   458,   459,   460,   460,
     461,   461,   461,   461,   461,   461,   461,   461,   462,   464,
     463,   466,   465,   467,   467,   469,   470,   468,   468,   468,
     468,   468,   468,   468,   468,   468,   468,   468,   468,   468,
     471,   468,   472,   468,   468,   468,   468,   468,   468,   468,
     468,   468,   468,   468,   468,   468,   468,   468,   468,   473,
     468,   474,   468,   468,   468,   475,   476,   468,   477,   478,
     468,   468,   479,   468,   468,   480,   468,   481,   468,   468,
     482,   468,   468,   483,   468,   484,   468,   485,   468,   468,
     468,   468,   486,   487,   468,   468,   468,   468,   468,   468,
     468,   468,   468,   468,   468,   468,   468,   468,   468,   468,
     488,   468,   489,   468,   490,   468,   468,   491,   468,   492,
     468,   493,   468,   468,   468,   494,   468,   495,   495,   496,
     496,   497,   497,   498,   500,   501,   502,   503,   499,   504,
     505,   499,   506,   499,   507,   507,   508,   509,   509,   509,
     510,   511,   510,   510,   512,   512,   513,   513,   514,   514,
     515,   515,   515,   516,   516,   517,   517,   518,   518,   518,
     519,   519,   519,   520,   521,   520,   522,   522,   523,   524,
     524,   525,   525,   525,   527,   528,   529,   526,   530,   531,
     530,   532,   530,   534,   533,   535,   535,   536,   536,   537,
     537,   538,   538,   538,   539,   539,   540,   541,   541,   542,
     543,   543,   544,   545,   545,   545,   545,   545,   545,   546,
     546,   546,   546,   547,   547,   548,   548,   548,   548,   549,
     549,   551,   550,   552,   552,   554,   553,   555,   555,   556,
     556,   557,   558,   559,   558,   561,   560,   562,   563,   563,
     563,   565,   566,   567,   564,   568,   568,   569,   569,   569,
     570,   569,   571,   571,   572,   573,   573,   574,   574,   574,
     575,   574,   574,   576,   576,   577,   577,   577,   578,   578,
     578,   578,   580,   579,   581,   581,   581,   581,   581,   581,
     581,   581,   581,   581,   581,   581,   581,   581,   581,   581,
     582,   584,   583,   585,   585,   586,   587,   586,   589,   588,
     591,   590,   592,   594,   593,   595,   595,   596,   596,   597,
     597,   598,   598,   600,   599,   601,   601,   599,   599,   599,
     602,   603,   603,   604,   606,   605,   607,   607,   608,   608,
     608,   608,   608,   608,   608,   608,   608,   610,   609,   611,
     612,   613,   614,   614,   615,   615,   616,   617,   617,   618,
     619,   620,   621,   621,   622,   623,   623,   624,   625,   626,
     627,   629,   630,   631,   628,   632,   632,   633,   633,   634,
     634,   635,   635,   635,   635,   635,   635,   635,   636,   637,
     638,   640,   639,   641,   641,   639,   639,   639,   643,   644,
     645,   646,   642,   647,   647,   648,   648,   648,   648,   648,
     648,   649,   651,   650,   653,   652,   654,   654,   655,   655,
     655,   655,   656,   656,   656,   657,   658,   658,   659,   659,
     659,   661,   660,   662,   662,   664,   663,   665,   665,   667,
     666,   669,   668,   671,   670,   672,   672,   673,   673,   673,
     673,   673,   673,   673,   673,   673,   673,   673,   673,   673,
     673,   673,   673,   673,   673,   673,   673,   674,   673,   675,
     675,   676,   677,   677,   678,   678,   678,   679,   679,   679,
     680,   681,   681,   681,   681,   681,   681,   681,   681,   681,
     681,   681,   681,   681,   681,   681,   682,   682,   682,   682,
     682,   682,   683,   683,   683,   683,   683,   683,   684,   684,
     684,   684,   684,   684,   685,   686,   687,   687,   687,   688,
     689,   690,   690,   690,   690,   691,   693,   692,   695,   694,
     696,   696,   697,   698,   699,   700,   702,   701,   704,   703,
     705,   705,   706,   706,   706,   706,   706,   706,   707,   706,
     706,   706,   706,   706,   706,   706,   706,   706,   706,   706,
     706,   706,   706,   708,   706,   709,   706,   710,   706,   711,
     706,   706,   706,   706,   706,   706,   706,   706,   706,   706,
     706,   706,   712,   706,   706,   706,   706,   706,   706,   706,
     706,   706,   706,   706,   706,   706,   706,   706,   706,   713,
     706,   714,   706,   715,   706,   716,   706,   717,   717,   717,
     717,   718,   718,   719,   719,   719,   720,   720,   720,   720,
     720,   721,   722,   722,   723,   723,   723,   723,   723,   723,
     724,   724,   724,   724,   725,   725,   725,   725,   726,   728,
     729,   727,   727,   727,   727,   727,   727,   727,   727,   727,
     730,   730,   731,   731,   732,   732,   732,   733,   734,   735,
     735,   737,   736,   738,   736,   739,   740,   740,   742,   743,
     741,   744,   745,   741,   741,   741,   746,   746,   747,   748,
     748,   749,   749,   750,   751,   752,   752,   754,   755,   753,
     756,   756,   757,   759,   758,   760,   761,   762,   763,   763,
     765,   764,   766,   764,   767,   764,   764,   764,   764,   764,
     764,   764,   764,   764,   764,   764,   764,   764,   764,   764,
     764,   764,   764,   764,   764,   768,   768,   768,   769,   769,
     769,   769,   770,   770,   770,   771,   771,   771,   772,   772,
     773,   773,   774,   775,   775,   776,   776,   776,   777,   777,
     778,   778,   778,   779,   779,   780,   780,   782,   781,   784,
     783,   786,   785,   787,   787,   789,   788,   790,   788,   791,
     788,   792,   788,   788,   793,   788,   788,   788,   794,   795,
     795,   797,   796,   798,   796,   799,   799,   800,   802,   801,
     804,   803,   806,   805,   807,   805,   808,   805,   809,   809,
     809,   810,   810,   811,   811,   812,   812,   812,   812,   812,
     812,   812,   812,   813,   813,   813,   813,   813,   813,   813,
     813,   813,   813,   813,   813,   813,   813,   813,   814,   814,
     814,   814,   815,   815,   815,   815,   815,   815,   815,   815,
     815,   817,   816,   818,   818,   820,   819,   821,   819,   822,
     819,   823,   819,   824,   819,   825,   819,   826,   819,   827,
     819,   828,   828,   828,   828,   828,   829,   829,   829,   829,
     829,   830,   831,   830,   832,   832,   833,   833,   834,   834,
     835,   835,   836,   836,   837,   837,   838,   839,   838,   838,
     840,   841,   842,   840,   843,   840,   840,   844,   840,   840,
     840,   845,   840,   840,   840,   846,   846,   847,   847,   848,
     849,   848,   851,   850,   852,   853,   854,   855,   857,   856,
     858,   859,   859,   860,   860,   862,   861,   863,   863,   864,
     864,   866,   865,   867,   867,   869,   868,   870,   871,   871,
     872,   872,   874,   873,   875,   875,   876,   876,   878,   877,
     879,   879,   880,   881,   882,   883,   884,   885,   886,   886,
     887
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     3,     0,     4,     1,     3,     3,     0,     2,
       1,     0,     2,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     3,     3,     3,     3,     2,     3,     4,     3,
       1,     1,     1,     1,     1,     1,     4,     1,     0,     2,
       4,     4,     4,     4,     4,     4,     4,     4,     3,     0,
       3,     0,     3,     0,     2,     0,     0,     9,     3,     3,
       3,     4,     3,     4,     3,     4,     3,     3,     3,     3,
       0,     6,     0,     9,     3,     4,     7,     4,     7,     3,
       3,     3,     3,     3,     3,     3,     3,     6,     6,     0,
       4,     0,     4,     4,     4,     0,     0,     9,     0,     0,
       9,     3,     0,     4,     3,     0,     4,     0,     5,     3,
       0,     4,     3,     0,     4,     0,     5,     0,     4,     2,
       3,     3,     0,     0,     9,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     4,     3,     3,     3,     3,
       0,     5,     0,     9,     0,     5,     7,     0,     4,     0,
       7,     0,     7,     3,     3,     0,     5,     0,     1,     0,
       2,     0,     2,     4,     0,     0,     0,     0,    11,     0,
       0,     9,     0,     9,     0,     2,     4,     0,     1,     1,
       0,     0,     4,     2,     0,     2,     0,     2,     0,     2,
       0,     1,     1,     0,     4,     0,     2,     1,     2,     2,
       1,     1,     1,     1,     0,     7,     0,     2,     1,     0,
       1,     1,     1,     1,     0,     0,     0,    12,     0,     0,
       5,     0,     5,     0,     5,     0,     2,     0,     2,     1,
       2,     2,     2,     2,     1,     2,     4,     1,     2,     4,
       1,     2,     4,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     0,     2,     1,     1,     1,     1,     0,
       2,     0,     4,     0,     2,     0,     6,     0,     2,     0,
       2,     6,     3,     0,     7,     0,     4,     1,     2,     3,
       3,     0,     0,     0,    26,     0,     2,     4,     4,     6,
       0,     4,     1,     1,     2,     0,     2,     1,     1,     3,
       0,     4,     1,     1,     2,     2,     2,     2,     2,     3,
       4,     3,     0,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       2,     0,     4,     0,     2,     5,     0,     8,     0,     3,
       0,     3,     5,     0,     7,     0,     1,     1,     2,     0,
       1,     1,     2,     0,     4,     1,     2,     2,     2,     2,
       2,     0,     2,     3,     0,     4,     0,     2,     3,     3,
       4,     5,     4,     5,     3,     3,     3,     0,     3,     3,
       1,     2,     0,     2,     5,     6,     1,     0,     2,     3,
       1,     2,     0,     2,     3,     0,     2,     2,     2,     4,
       3,     0,     0,     0,     8,     1,     2,     0,     2,     0,
       2,     1,     1,     1,     1,     1,     1,     1,     3,     3,
       4,     0,     4,     1,     2,     2,     2,     2,     0,     0,
       0,     0,    12,     0,     2,     3,     3,     4,     4,     3,
       3,     3,     0,     3,     0,     3,     0,     2,     5,     1,
       1,     1,     3,     3,     3,     3,     0,     2,     1,     1,
       1,     0,     4,     0,     2,     0,     3,     2,     4,     0,
       4,     0,     3,     0,     3,     0,     2,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     3,     3,     1,     1,     1,     1,     0,     4,     1,
       2,     3,     0,     2,     1,     1,     1,     2,     2,     2,
       3,     1,     2,     1,     1,     2,     2,     1,     1,     1,
       1,     2,     1,     1,     2,     2,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     4,     3,     3,     3,     3,
       3,     2,     3,     4,     3,     2,     0,     4,     0,     4,
       3,     3,     1,     1,     5,     3,     0,     3,     0,     3,
       0,     2,     2,     3,     4,     3,     4,     5,     0,     4,
       3,     1,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     0,     4,     0,     5,     0,     5,     0,
       5,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       4,     3,     0,     4,     4,     2,     4,     4,     4,     3,
       3,     4,     4,     4,     4,     4,     4,     4,     4,     0,
       4,     0,     4,     0,     4,     0,     4,     1,     1,     1,
       1,     1,     2,     2,     2,     2,     3,     3,     4,     3,
       3,     1,     0,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     0,     1,     1,     1,     2,     0,
       0,     7,     3,     5,     7,     5,     7,     7,     9,     1,
       0,     2,     0,     1,     0,     2,     2,     1,     1,     0,
       2,     0,     6,     0,     8,     7,    11,     4,     0,     0,
      10,     0,     0,    10,     6,     6,     0,     2,     1,     6,
       6,     3,     2,     1,     4,     0,     2,     0,     0,     7,
       0,     2,     5,     0,     4,     3,     1,     2,     0,     2,
       0,     4,     0,     4,     0,    10,     9,     3,     3,     4,
       4,     4,     4,     4,     4,     4,     4,     3,     7,     8,
       6,     3,     3,     3,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     2,     5,     1,     2,     0,     4,     7,     1,     1,
       1,     1,     1,     1,     2,     1,     2,     0,     4,     0,
       3,     0,     3,     0,     2,     0,     4,     0,     4,     0,
       4,     0,     4,     4,     0,     4,     5,     1,     2,     0,
       2,     0,     4,     0,     4,     0,     2,     5,     0,     6,
       0,     6,     0,     6,     0,     6,     0,     6,     0,     1,
       1,     1,     2,     1,     2,     3,     3,     3,     3,     2,
       3,     6,     1,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     2,     3,     6,     1,     1,     3,     3,
       6,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     0,     5,     0,     2,     0,     5,     0,     5,     0,
       5,     0,     5,     0,     5,     0,     5,     0,     5,     0,
       5,     3,     3,     1,     2,     2,     0,     1,     2,     5,
       3,     0,     0,     6,     0,     1,     0,     1,     0,     3,
       0,     1,     0,     1,     0,     2,     1,     0,     3,     1,
       0,     0,     0,     5,     0,     6,     2,     0,     5,     2,
       5,     0,     6,     2,     6,     0,     1,     0,     1,     0,
       0,     3,     0,     3,     4,     3,     3,     3,     0,     7,
       2,     1,     2,     3,     1,     0,     5,     1,     2,     1,
       2,     0,     7,     1,     2,     0,     7,     2,     1,     2,
       3,     1,     0,     5,     1,     2,     1,     2,     0,     7,
       1,     2,     3,     3,     3,     3,     3,     3,     0,     1,
       1
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYNOMEM         goto yyexhaustedlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Backward compatibility with an undocumented macro.
   Use YYerror or YYUNDEF. */
#define YYERRCODE YYUNDEF


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)




# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Kind, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  if (!yyvaluep)
    return;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo,
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp,
                 int yyrule)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       YY_ACCESSING_SYMBOL (+yyssp[yyi + 1 - yynrhs]),
                       &yyvsp[(yyi + 1) - (yynrhs)]);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) ((void) 0)
# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif






/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg,
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep)
{
  YY_USE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/* Lookahead token kind.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;




/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    yy_state_fast_t yystate = 0;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus = 0;

    /* Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* Their size.  */
    YYPTRDIFF_T yystacksize = YYINITDEPTH;

    /* The state stack: array, bottom, top.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss = yyssa;
    yy_state_t *yyssp = yyss;

    /* The semantic value stack: array, bottom, top.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp = yyvs;

  int yyn;
  /* The return value of yyparse.  */
  int yyresult;
  /* Lookahead symbol kind.  */
  yysymbol_kind_t yytoken = YYSYMBOL_YYEMPTY;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yychar = YYEMPTY; /* Cause a token to be read.  */

  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END
  YY_STACK_PRINT (yyss, yyssp);

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    YYNOMEM;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        YYNOMEM;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          YYNOMEM;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */


  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token\n"));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = YYEOF;
      yytoken = YYSYMBOL_YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else if (yychar == YYerror)
    {
      /* The scanner already issued an error message, process directly
         to error recovery.  But do not keep the error token as
         lookahead, it is too special and may lead us to an endless
         loop in error recovery. */
      yychar = YYUNDEF;
      yytoken = YYSYMBOL_YYerror;
      goto yyerrlab1;
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 2: /* lef_file: rules extension_opt end_library  */
#line 335 "lef/lef.y"
      {
        // 11/16/2001 - Wanda da Rosa - pcr 408334
        // Return 1 if there are errors
        if (lefData->lef_errors)
           return 1;
        if (!lefData->hasVer) {
              char temp[300];
              sprintf(temp, "No VERSION statement found, using the default value %2g.", lefData->versionNum);
              lefWarning(2001, temp);            
        }        
        //only pre 5.6, 5.6 it is obsolete
        if (!lefData->hasNameCase && lefData->versionNum < 5.6)
           lefWarning(2002, "NAMESCASESENSITIVE is a required statement on LEF file with version 5.5 and earlier.\nWithout NAMESCASESENSITIVE defined, the LEF file is technically incorrect.\nRefer the LEF/DEF 5.5 or earlier Language Referece manual on how to define this statement.");
        if (!lefData->hasBusBit && lefData->versionNum < 5.6)
           lefWarning(2003, "BUSBITCHARS is a required statement on LEF file with version 5.5 and earlier.\nWithout BUSBITCHARS defined, the LEF file is technically incorrect.\nRefer the LEF/DEF 5.5 or earlier Language Referece manual on how to define this statement.");
        if (!lefData->hasDivChar && lefData->versionNum < 5.6)
           lefWarning(2004, "DIVIDERCHAR is a required statementon LEF file with version 5.5 and earlier.\nWithout DIVIDECHAR defined, the LEF file is technically incorrect.\nRefer the LEF/DEF 5.5 or earlier Language Referece manual on how to define this statement.");

      }
#line 3879 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 3: /* $@1: %empty  */
#line 355 "lef/lef.y"
                   { lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 3885 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 4: /* version: K_VERSION $@1 T_STRING ';'  */
#line 356 "lef/lef.y"
      { 
		 // More than 1 VERSION in lef file within the open file - It's wrong syntax, 
		 // but copy old behavior - initialize lef reading.
         if (lefData->hasVer)     
         {
			lefData->initRead();
		 }

         lefData->versionNum = convert_name2num((yyvsp[-1].string));
         if (lefData->versionNum > CURRENT_VERSION) {
            char temp[120];
            sprintf(temp,
               "Lef parser %.1f does not support lef file with version %s. Parser will stop processing.", CURRENT_VERSION, (yyvsp[-1].string));
            lefError(1503, temp);
            return 1;
         }

         if (lefCallbacks->VersionStrCbk) {
            CALLBACK(lefCallbacks->VersionStrCbk, lefrVersionStrCbkType, (yyvsp[-1].string));
         } else {
            if (lefCallbacks->VersionCbk)
               CALLBACK(lefCallbacks->VersionCbk, lefrVersionCbkType, lefData->versionNum);
         }
         if (lefData->versionNum > 5.3 && lefData->versionNum < 5.4) {
            lefData->ignoreVersion = 1;
         }
         lefData->use5_3 = lefData->use5_4 = 0;
         lefData->lef_errors = 0;
         lefData->hasVer = 1;
         if (lefData->versionNum < 5.6) {
            lefData->doneLib = 0;
            lefData->namesCaseSensitive = lefSettings->CaseSensitive;
         } else {
            lefData->doneLib = 1;
            lefData->namesCaseSensitive = 1;
         }
      }
#line 3927 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 5: /* int_number: NUMBER  */
#line 395 "lef/lef.y"
      {
         // int_number represent 'integer-like' type. It can have fraction and exponent part 
         // but the value shouldn't exceed the 64-bit integer limit. 
         if (!(( yylval.dval >= lefData->leflVal) && ( yylval.dval <= lefData->lefrVal))) { // YES, it isn't really a number 
            char *str = (char*) lefMalloc(strlen(lefData->current_token) + strlen(lefData->lefrFileName) + 350);
            sprintf(str, "ERROR (LEFPARS-203) Number has exceeded the limit for an integer. See file %s at line %d.\n",
                    lefData->lefrFileName, lefData->lef_nlines);
            fflush(stdout);
            lefiError(0, 203, str);
            free(str);
            lefData->lef_errors++;
        }

        (yyval.dval) = yylval.dval ;
      }
#line 3947 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 6: /* dividerchar: K_DIVIDERCHAR QSTRING ';'  */
#line 412 "lef/lef.y"
      {
        if (lefCallbacks->DividerCharCbk) {
          if (strcmp((yyvsp[-1].string), "") != 0) {
             CALLBACK(lefCallbacks->DividerCharCbk, lefrDividerCharCbkType, (yyvsp[-1].string));
          } else {
             CALLBACK(lefCallbacks->DividerCharCbk, lefrDividerCharCbkType, "/");
             lefWarning(2005, "DIVIDERCHAR has an invalid null value. Value is set to default /");
          }
        }
        lefData->hasDivChar = 1;
      }
#line 3963 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 7: /* busbitchars: K_BUSBITCHARS QSTRING ';'  */
#line 425 "lef/lef.y"
      {
        if (lefCallbacks->BusBitCharsCbk) {
          if (strcmp((yyvsp[-1].string), "") != 0) {
             CALLBACK(lefCallbacks->BusBitCharsCbk, lefrBusBitCharsCbkType, (yyvsp[-1].string)); 
          } else {
             CALLBACK(lefCallbacks->BusBitCharsCbk, lefrBusBitCharsCbkType, "[]"); 
             lefWarning(2006, "BUSBITCHAR has an invalid null value. Value is set to default []");
          }
        }
        lefData->hasBusBit = 1;
      }
#line 3979 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 10: /* rules: error  */
#line 440 "lef/lef.y"
            { }
#line 3985 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 11: /* end_library: %empty  */
#line 443 "lef/lef.y"
      {
        if (lefData->versionNum >= 5.6) {
           lefData->doneLib = 1;
           lefData->ge56done = 1;
        }
      }
#line 3996 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 12: /* end_library: K_END K_LIBRARY  */
#line 450 "lef/lef.y"
      {
        lefData->doneLib = 1;
        lefData->ge56done = 1;
        if (lefCallbacks->LibraryEndCbk)
          CALLBACK(lefCallbacks->LibraryEndCbk, lefrLibraryEndCbkType, 0);
        // 11/16/2001 - Wanda da Rosa - pcr 408334
        // Return 1 if there are errors
      }
#line 4009 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 52: /* case_sensitivity: K_NAMESCASESENSITIVE K_ON ';'  */
#line 474 "lef/lef.y"
          {
            if (lefData->versionNum < 5.6) {
              lefData->namesCaseSensitive = TRUE;
              if (lefCallbacks->CaseSensitiveCbk)
                CALLBACK(lefCallbacks->CaseSensitiveCbk, 
                         lefrCaseSensitiveCbkType,
                         lefData->namesCaseSensitive);
              lefData->hasNameCase = 1;
            } else
              if (lefCallbacks->CaseSensitiveCbk) // write warning only if cbk is set 
                 if (lefData->caseSensitiveWarnings++ < lefSettings->CaseSensitiveWarnings)
                   lefWarning(2007, "NAMESCASESENSITIVE statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
          }
#line 4027 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 53: /* case_sensitivity: K_NAMESCASESENSITIVE K_OFF ';'  */
#line 488 "lef/lef.y"
          {
            if (lefData->versionNum < 5.6) {
              lefData->namesCaseSensitive = FALSE;
              if (lefCallbacks->CaseSensitiveCbk)
                CALLBACK(lefCallbacks->CaseSensitiveCbk, lefrCaseSensitiveCbkType,
                               lefData->namesCaseSensitive);
              lefData->hasNameCase = 1;
            } else {
              if (lefCallbacks->CaseSensitiveCbk) { // write error only if cbk is set 
                if (lefData->caseSensitiveWarnings++ < lefSettings->CaseSensitiveWarnings) {
                  lefError(1504, "NAMESCASESENSITIVE statement is set with OFF.\nStarting version 5.6, NAMESCASENSITIVE is obsolete,\nif it is defined, it has to have the ON value.\nParser will stop processing.");
                  CHKERR();
                }
              }
            }
          }
#line 4048 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 54: /* wireextension: K_NOWIREEXTENSIONATPIN K_ON ';'  */
#line 506 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->NoWireExtensionCbk)
          CALLBACK(lefCallbacks->NoWireExtensionCbk, lefrNoWireExtensionCbkType, "ON");
      } else
        if (lefCallbacks->NoWireExtensionCbk) // write warning only if cbk is set 
           if (lefData->noWireExtensionWarnings++ < lefSettings->NoWireExtensionWarnings)
             lefWarning(2008, "NOWIREEXTENSIONATPIN statement is obsolete in version 5.6 or later.\nThe NOWIREEXTENSIONATPIN statement will be ignored.");
    }
#line 4062 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 55: /* wireextension: K_NOWIREEXTENSIONATPIN K_OFF ';'  */
#line 516 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->NoWireExtensionCbk)
          CALLBACK(lefCallbacks->NoWireExtensionCbk, lefrNoWireExtensionCbkType, "OFF");
      } else
        if (lefCallbacks->NoWireExtensionCbk) // write warning only if cbk is set 
           if (lefData->noWireExtensionWarnings++ < lefSettings->NoWireExtensionWarnings)
             lefWarning(2008, "NOWIREEXTENSIONATPIN statement is obsolete in version 5.6 or later.\nThe NOWIREEXTENSIONATPIN statement will be ignored.");
    }
#line 4076 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 56: /* fixedmask: K_FIXEDMASK ';'  */
#line 527 "lef/lef.y"
    { 
       if (lefData->versionNum >= 5.8) {
       
          if (lefCallbacks->FixedMaskCbk) {
            lefData->lefFixedMask = 1;
            CALLBACK(lefCallbacks->FixedMaskCbk, lefrFixedMaskCbkType, lefData->lefFixedMask);
          }
          
          lefData->hasFixedMask = 1;
       }
    }
#line 4092 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 57: /* manufacturing: K_MANUFACTURINGGRID int_number ';'  */
#line 540 "lef/lef.y"
    {
      if (lefCallbacks->ManufacturingCbk)
        CALLBACK(lefCallbacks->ManufacturingCbk, lefrManufacturingCbkType, (yyvsp[-1].dval));
      lefData->hasManufactur = 1;
    }
#line 4102 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 58: /* useminspacing: K_USEMINSPACING spacing_type spacing_value ';'  */
#line 547 "lef/lef.y"
  {
    if ((strcmp((yyvsp[-2].string), "PIN") == 0) && (lefData->versionNum >= 5.6)) {
      if (lefCallbacks->UseMinSpacingCbk) // write warning only if cbk is set 
         if (lefData->useMinSpacingWarnings++ < lefSettings->UseMinSpacingWarnings)
            lefWarning(2009, "USEMINSPACING PIN statement is obsolete in version 5.6 or later.\n The USEMINSPACING PIN statement will be ignored.");
    } else {
        if (lefCallbacks->UseMinSpacingCbk) {
          lefData->lefrUseMinSpacing.set((yyvsp[-2].string), (yyvsp[-1].integer));
          CALLBACK(lefCallbacks->UseMinSpacingCbk, lefrUseMinSpacingCbkType,
                   &lefData->lefrUseMinSpacing);
      }
    }
  }
#line 4120 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 59: /* clearancemeasure: K_CLEARANCEMEASURE clearance_type ';'  */
#line 562 "lef/lef.y"
    { CALLBACK(lefCallbacks->ClearanceMeasureCbk, lefrClearanceMeasureCbkType, (yyvsp[-1].string)); }
#line 4126 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 60: /* clearance_type: K_MAXXY  */
#line 565 "lef/lef.y"
            {(yyval.string) = (char*)"MAXXY";}
#line 4132 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 61: /* clearance_type: K_EUCLIDEAN  */
#line 566 "lef/lef.y"
                  {(yyval.string) = (char*)"EUCLIDEAN";}
#line 4138 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 62: /* spacing_type: K_OBS  */
#line 569 "lef/lef.y"
            {(yyval.string) = (char*)"OBS";}
#line 4144 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 63: /* spacing_type: K_PIN  */
#line 570 "lef/lef.y"
            {(yyval.string) = (char*)"PIN";}
#line 4150 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 64: /* spacing_value: K_ON  */
#line 573 "lef/lef.y"
            {(yyval.integer) = 1;}
#line 4156 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 65: /* spacing_value: K_OFF  */
#line 574 "lef/lef.y"
            {(yyval.integer) = 0;}
#line 4162 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 66: /* units_section: start_units units_rules K_END K_UNITS  */
#line 577 "lef/lef.y"
    { 
      if (lefCallbacks->UnitsCbk)
        CALLBACK(lefCallbacks->UnitsCbk, lefrUnitsCbkType, &lefData->lefrUnits);
    }
#line 4171 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 67: /* start_units: K_UNITS  */
#line 583 "lef/lef.y"
    {
      lefData->lefrUnits.clear();
      if (lefData->hasManufactur) {
        if (lefData->unitsWarnings++ < lefSettings->UnitsWarnings) {
          lefError(1505, "MANUFACTURINGGRID statement was defined before UNITS.\nRefer the LEF Language Reference manual for the order of LEF statements.");
          CHKERR();
        }
      }
      if (lefData->hasMinfeature) {
        if (lefData->unitsWarnings++ < lefSettings->UnitsWarnings) {
          lefError(1712, "MINFEATURE statement was defined before UNITS.\nRefer the LEF Language Reference manual for the order of LEF statements.");
          CHKERR();
        }
      }
      if (lefData->versionNum < 5.6) {
        if (lefData->hasSite) {//SITE is defined before UNIT and is illegal in pre 5.6
          lefError(1713, "SITE statement was defined before UNITS.\nRefer the LEF Language Reference manual for the order of LEF statements.");
          CHKERR();
        }
      }
    }
#line 4197 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 70: /* units_rule: K_TIME K_NANOSECONDS int_number ';'  */
#line 610 "lef/lef.y"
    { if (lefCallbacks->UnitsCbk) lefData->lefrUnits.setTime((yyvsp[-1].dval)); }
#line 4203 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 71: /* units_rule: K_CAPACITANCE K_PICOFARADS int_number ';'  */
#line 612 "lef/lef.y"
    { if (lefCallbacks->UnitsCbk) lefData->lefrUnits.setCapacitance((yyvsp[-1].dval)); }
#line 4209 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 72: /* units_rule: K_RESISTANCE K_OHMS int_number ';'  */
#line 614 "lef/lef.y"
    { if (lefCallbacks->UnitsCbk) lefData->lefrUnits.setResistance((yyvsp[-1].dval)); }
#line 4215 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 73: /* units_rule: K_POWER K_MILLIWATTS int_number ';'  */
#line 616 "lef/lef.y"
    { if (lefCallbacks->UnitsCbk) lefData->lefrUnits.setPower((yyvsp[-1].dval)); }
#line 4221 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 74: /* units_rule: K_CURRENT K_MILLIAMPS int_number ';'  */
#line 618 "lef/lef.y"
    { if (lefCallbacks->UnitsCbk) lefData->lefrUnits.setCurrent((yyvsp[-1].dval)); }
#line 4227 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 75: /* units_rule: K_VOLTAGE K_VOLTS int_number ';'  */
#line 620 "lef/lef.y"
    { if (lefCallbacks->UnitsCbk) lefData->lefrUnits.setVoltage((yyvsp[-1].dval)); }
#line 4233 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 76: /* units_rule: K_DATABASE K_MICRONS int_number ';'  */
#line 622 "lef/lef.y"
    { 
      if(validNum((int)(yyvsp[-1].dval))) {
         if (lefCallbacks->UnitsCbk)
            lefData->lefrUnits.setDatabase("MICRONS", (yyvsp[-1].dval));
      }
    }
#line 4244 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 77: /* units_rule: K_FREQUENCY K_MEGAHERTZ NUMBER ';'  */
#line 629 "lef/lef.y"
    { if (lefCallbacks->UnitsCbk) lefData->lefrUnits.setFrequency((yyvsp[-1].dval)); }
#line 4250 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 78: /* layer_rule: start_layer layer_options end_layer  */
#line 633 "lef/lef.y"
    { 
      if (lefCallbacks->LayerCbk)
        CALLBACK(lefCallbacks->LayerCbk, lefrLayerCbkType, &lefData->lefrLayer);
    }
#line 4259 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 79: /* $@2: %empty  */
#line 638 "lef/lef.y"
                     {lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 4265 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 80: /* start_layer: K_LAYER $@2 T_STRING  */
#line 639 "lef/lef.y"
    { 
      if (lefData->lefrHasMaxVS) {   // 5.5 
        if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
          if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
            lefError(1506, "A MAXVIASTACK statement is defined before the LAYER statement.\nRefer to the LEF Language Reference manual for the order of LEF statements.");
            CHKERR();
          }
        }
      }
      if (lefCallbacks->LayerCbk)
        lefData->lefrLayer.setName((yyvsp[0].string));
      lefData->useLenThr = 0;
      lefData->layerCut = 0;
      lefData->layerMastOver = 0;
      lefData->layerRout = 0;
      lefData->layerDir = 0;
      lefData->lefrHasLayer = 1;
      //strcpy(lefData->layerName, $3);
      lefData->layerName = strdup((yyvsp[0].string));
      lefData->hasType = 0;
      lefData->hasMask = 0;
      lefData->hasPitch = 0;
      lefData->hasWidth = 0;
      lefData->hasDirection = 0;
      lefData->hasParallel = 0;
      lefData->hasInfluence = 0;
      lefData->hasTwoWidths = 0;
      lefData->lefrHasSpacingTbl = 0;
      lefData->lefrHasSpacing = 0;
    }
#line 4300 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 81: /* $@3: %empty  */
#line 670 "lef/lef.y"
                 {lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 4306 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 82: /* end_layer: K_END $@3 T_STRING  */
#line 671 "lef/lef.y"
    { 
      if (strcmp(lefData->layerName, (yyvsp[0].string)) != 0) {
        if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
          if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
             lefData->outMsg = (char*)lefMalloc(10000);
             sprintf (lefData->outMsg,
                "END LAYER name %s is different from the LAYER name %s.\nCorrect the LEF file before rerunning it through the LEF parser.", (yyvsp[0].string), lefData->layerName);
             lefError(1507, lefData->outMsg);
             lefFree(lefData->outMsg);
             lefFree(lefData->layerName);
             CHKERR(); 
          } else
             lefFree(lefData->layerName);
        } else
          lefFree(lefData->layerName);
      } else
        lefFree(lefData->layerName);
      if (!lefSettings->RelaxMode) {
        if (lefData->hasType == 0) {
          if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefError(1508, "TYPE statement is a required statement in a LAYER and it is not defined.");
               CHKERR(); 
            }
          }
        }
        if ((lefData->layerRout == 1) && (lefData->hasPitch == 0)) {
          if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1509, "PITCH statement is a required statement in a LAYER with type ROUTING and it is not defined.");
              CHKERR(); 
            }
          }
        }
        if ((lefData->layerRout == 1) && (lefData->hasWidth == 0)) {
          if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1510, "WIDTH statement is a required statement in a LAYER with type ROUTING and it is not defined.");
              CHKERR(); 
            }
          }
        }
        if ((lefData->layerRout == 1) && (lefData->hasDirection == 0)) {
          if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg, "The DIRECTION statement which is required in a LAYER with TYPE ROUTING is not defined in LAYER %s.\nUpdate your lef file and add the DIRECTION statement for layer %s.", (yyvsp[0].string), (yyvsp[0].string));
              lefError(1511, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR(); 
            }
          }
        }
      }
    }
#line 4366 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 83: /* layer_options: %empty  */
#line 728 "lef/lef.y"
    { }
#line 4372 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 84: /* layer_options: layer_options layer_option  */
#line 730 "lef/lef.y"
    { }
#line 4378 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 85: /* $@4: %empty  */
#line 734 "lef/lef.y"
    {
       // let setArraySpacingCutSpacing to set the data 
    }
#line 4386 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 86: /* $@5: %empty  */
#line 740 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
         lefData->lefrLayer.setArraySpacingCut((yyvsp[0].dval));
         lefData->arrayCutsVal = 0;
      }
    }
#line 4397 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 87: /* layer_option: K_ARRAYSPACING $@4 layer_arraySpacing_long layer_arraySpacing_width K_CUTSPACING int_number $@5 layer_arraySpacing_arraycuts ';'  */
#line 747 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
         lefData->outMsg = (char*)lefMalloc(10000);
         sprintf(lefData->outMsg,
           "ARRAYSPACING is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
         lefError(1685, lefData->outMsg);
         lefFree(lefData->outMsg);
         CHKERR();
      }
    }
#line 4412 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 88: /* layer_option: K_TYPE layer_type ';'  */
#line 758 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
         lefData->lefrLayer.setType((yyvsp[-1].string));
      lefData->hasType = 1;
    }
#line 4422 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 89: /* layer_option: K_MASK int_number ';'  */
#line 764 "lef/lef.y"
    {
      if (lefData->versionNum < 5.8) {
          if (lefData->layerWarnings++ < lefSettings->ViaWarnings) {
              lefError(2081, "MASK information can only be defined with version 5.8");
              CHKERR(); 
          }           
      } else {
          if (lefCallbacks->LayerCbk) {
            lefData->lefrLayer.setMask((int)(yyvsp[-1].dval));
          }
          
          lefData->hasMask = 1;
      }
    }
#line 4441 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 90: /* layer_option: K_PITCH int_number ';'  */
#line 779 "lef/lef.y"
    { 
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setPitch((yyvsp[-1].dval));
      lefData->hasPitch = 1;  
    }
#line 4450 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 91: /* layer_option: K_PITCH int_number int_number ';'  */
#line 784 "lef/lef.y"
    { 
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setPitchXY((yyvsp[-2].dval), (yyvsp[-1].dval));
      lefData->hasPitch = 1;  
    }
#line 4459 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 92: /* layer_option: K_DIAGPITCH int_number ';'  */
#line 789 "lef/lef.y"
    { 
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setDiagPitch((yyvsp[-1].dval));
    }
#line 4467 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 93: /* layer_option: K_DIAGPITCH int_number int_number ';'  */
#line 793 "lef/lef.y"
    { 
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setDiagPitchXY((yyvsp[-2].dval), (yyvsp[-1].dval));
    }
#line 4475 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 94: /* layer_option: K_OFFSET int_number ';'  */
#line 797 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setOffset((yyvsp[-1].dval));
    }
#line 4483 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 95: /* layer_option: K_OFFSET int_number int_number ';'  */
#line 801 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setOffsetXY((yyvsp[-2].dval), (yyvsp[-1].dval));
    }
#line 4491 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 96: /* layer_option: K_DIAGWIDTH int_number ';'  */
#line 805 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setDiagWidth((yyvsp[-1].dval));
    }
#line 4499 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 97: /* layer_option: K_DIAGSPACING int_number ';'  */
#line 809 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setDiagSpacing((yyvsp[-1].dval));
    }
#line 4507 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 98: /* layer_option: K_WIDTH int_number ';'  */
#line 813 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setWidth((yyvsp[-1].dval));
      lefData->hasWidth = 1;  
    }
#line 4516 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 99: /* layer_option: K_AREA NUMBER ';'  */
#line 818 "lef/lef.y"
    {
      // Issue an error is this is defined in masterslice
      if (lefData->layerMastOver) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefError(1715, "It is incorrect to define an AREA statement in LAYER with TYPE MASTERSLICE or OVERLAP. Parser will stop processing.");
               CHKERR();
            }
         }
      }

      if (lefCallbacks->LayerCbk) {
         lefData->lefrLayer.setArea((yyvsp[-1].dval));
      }
    }
#line 4536 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 100: /* $@6: %empty  */
#line 834 "lef/lef.y"
    {
      lefData->hasSpCenter = 0;       // reset to 0, only once per spacing is allowed 
      lefData->hasSpSamenet = 0;
      lefData->hasSpParallel = 0;
      lefData->hasSpLayer = 0;
      lefData->layerCutSpacing = (yyvsp[0].dval);  // for error message purpose
      // 11/22/99 - Wanda da Rosa, PCR 283762
      //            Issue an error is this is defined in masterslice
      if (lefData->layerMastOver) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1512, "It is incorrect to define a SPACING statement in LAYER with TYPE MASTERSLICE or OVERLAP. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      // 5.5 either SPACING or SPACINGTABLE, not both for routing layer only
      if (lefData->layerRout) {
        if (lefData->lefrHasSpacingTbl && lefData->versionNum < 5.7) {
           if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
              if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                lefWarning(2010, "It is incorrect to have both SPACING rules & SPACINGTABLE rules within a ROUTING layer");
              }
           }
        }
        if (lefCallbacks->LayerCbk)
           lefData->lefrLayer.setSpacingMin((yyvsp[0].dval));
        lefData->lefrHasSpacing = 1;
      } else { 
        if (lefCallbacks->LayerCbk)
           lefData->lefrLayer.setSpacingMin((yyvsp[0].dval));
      }
    }
#line 4574 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 101: /* layer_option: K_SPACING int_number $@6 layer_spacing_opts layer_spacing_cut_routing ';'  */
#line 868 "lef/lef.y"
                                  {}
#line 4580 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 102: /* $@7: %empty  */
#line 870 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
         lefData->lefrLayer.setSpacingTableOrtho();
      if (lefCallbacks->LayerCbk) // due to converting to C, else, convertor produce 
         lefData->lefrLayer.addSpacingTableOrthoWithin((yyvsp[-2].dval), (yyvsp[0].dval));//bad code
    }
#line 4591 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 103: /* layer_option: K_SPACINGTABLE K_ORTHOGONAL K_WITHIN int_number K_SPACING int_number $@7 layer_spacingtable_opts ';'  */
#line 877 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
         lefData->outMsg = (char*)lefMalloc(10000);
         sprintf(lefData->outMsg,
           "SPACINGTABLE ORTHOGONAL is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
         lefError(1694, lefData->outMsg);
         lefFree(lefData->outMsg);
         CHKERR();
      }
    }
#line 4606 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 104: /* layer_option: K_DIRECTION layer_direction ';'  */
#line 888 "lef/lef.y"
    {
      lefData->layerDir = 1;
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1513, "DIRECTION statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setDirection((yyvsp[-1].string));
      lefData->hasDirection = 1;  
    }
#line 4624 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 105: /* layer_option: K_RESISTANCE K_RPERSQ int_number ';'  */
#line 902 "lef/lef.y"
    {
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1514, "RESISTANCE statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setResistance((yyvsp[-1].dval));
    }
#line 4640 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 106: /* layer_option: K_RESISTANCE K_RPERSQ K_PWL '(' res_points ')' ';'  */
#line 914 "lef/lef.y"
    {
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1515, "RESISTANCE statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      }
    }
#line 4655 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 107: /* layer_option: K_CAPACITANCE K_CPERSQDIST int_number ';'  */
#line 925 "lef/lef.y"
    {
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1516, "CAPACITANCE statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setCapacitance((yyvsp[-1].dval));
    }
#line 4671 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 108: /* layer_option: K_CAPACITANCE K_CPERSQDIST K_PWL '(' cap_points ')' ';'  */
#line 937 "lef/lef.y"
    {
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1517, "CAPACITANCE statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      }
    }
#line 4686 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 109: /* layer_option: K_HEIGHT int_number ';'  */
#line 948 "lef/lef.y"
    {
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1518, "HEIGHT statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setHeight((yyvsp[-1].dval));
    }
#line 4702 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 110: /* layer_option: K_WIREEXTENSION int_number ';'  */
#line 960 "lef/lef.y"
    {
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1519, "WIREEXTENSION statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setWireExtension((yyvsp[-1].dval));
    }
#line 4718 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 111: /* layer_option: K_THICKNESS int_number ';'  */
#line 972 "lef/lef.y"
    {
      if (!lefData->layerRout && (lefData->layerCut || lefData->layerMastOver)) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1520, "THICKNESS statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setThickness((yyvsp[-1].dval));
    }
#line 4734 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 112: /* layer_option: K_SHRINKAGE int_number ';'  */
#line 984 "lef/lef.y"
    {
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1521, "SHRINKAGE statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setShrinkage((yyvsp[-1].dval));
    }
#line 4750 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 113: /* layer_option: K_CAPMULTIPLIER int_number ';'  */
#line 996 "lef/lef.y"
    {
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1522, "CAPMULTIPLIER statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setCapMultiplier((yyvsp[-1].dval));
    }
#line 4766 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 114: /* layer_option: K_EDGECAPACITANCE int_number ';'  */
#line 1008 "lef/lef.y"
    {
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1523, "EDGECAPACITANCE statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setEdgeCap((yyvsp[-1].dval));
    }
#line 4782 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 115: /* layer_option: K_ANTENNALENGTHFACTOR int_number ';'  */
#line 1021 "lef/lef.y"
    { // 5.3 syntax 
      lefData->use5_3 = 1;
      if (!lefData->layerRout && (lefData->layerCut || lefData->layerMastOver)) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1525, "ANTENNALENGTHFACTOR statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      } else if (lefData->versionNum >= 5.4) {
         if (lefData->use5_4) {
            if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
               if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                  lefData->outMsg = (char*)lefMalloc(10000);
                  sprintf (lefData->outMsg,
                    "ANTENNALENGTHFACTOR statement is a version 5.3 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNALENGTHFACTOR syntax, which is incorrect.", lefData->versionNum);
                  lefError(1526, lefData->outMsg);
                  lefFree(lefData->outMsg);
                  CHKERR();
               }
            }
         }
      }

      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setAntennaLength((yyvsp[-1].dval));
    }
#line 4813 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 116: /* layer_option: K_CURRENTDEN int_number ';'  */
#line 1048 "lef/lef.y"
    {
      if (lefData->versionNum < 5.2) {
         if (!lefData->layerRout) {
            if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
               if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                 lefError(1702, "CURRENTDEN statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
                 CHKERR();
               }
            }
         }
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setCurrentDensity((yyvsp[-1].dval));
      } else {
         if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
            lefWarning(2079, "CURRENTDEN statement is obsolete in version 5.2 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.2 or later.");
            CHKERR();
         }
      }
    }
#line 4836 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 117: /* layer_option: K_CURRENTDEN K_PWL '(' current_density_pwl_list ')' ';'  */
#line 1067 "lef/lef.y"
    { 
      if (lefData->versionNum < 5.2) {
         if (!lefData->layerRout) {
            if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
               if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                 lefError(1702, "CURRENTDEN statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
                 CHKERR();
               }
            }
         }
      } else {
         if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
            lefWarning(2079, "CURRENTDEN statement is obsolete in version 5.2 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.2 or later.");
            CHKERR();
         }
      }
    }
#line 4858 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 118: /* layer_option: K_CURRENTDEN '(' int_number int_number ')' ';'  */
#line 1085 "lef/lef.y"
    {
      if (lefData->versionNum < 5.2) {
         if (!lefData->layerRout) {
            if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
               if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                 lefError(1702, "CURRENTDEN statement can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
                 CHKERR();
               }
            }
         }
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setCurrentPoint((yyvsp[-3].dval), (yyvsp[-2].dval));
      } else {
         if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
            lefWarning(2079, "CURRENTDEN statement is obsolete in version 5.2 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.2 or later.");
            CHKERR();
         }
      }
    }
#line 4881 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 119: /* $@8: %empty  */
#line 1103 "lef/lef.y"
               { lefData->lefDumbMode = 10000000;}
#line 4887 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 120: /* layer_option: K_PROPERTY $@8 layer_prop_list ';'  */
#line 1104 "lef/lef.y"
    {
      lefData->lefDumbMode = 0;
    }
#line 4895 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 121: /* $@9: %empty  */
#line 1108 "lef/lef.y"
    {
      if (lefData->layerMastOver) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1527, "ACCURRENTDENSITY statement can't be defined in LAYER with TYPE MASTERSLICE or OVERLAP. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.addAccurrentDensity((yyvsp[0].string));      
    }
#line 4911 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 122: /* layer_option: K_ACCURRENTDENSITY layer_table_type $@9 layer_frequency  */
#line 1119 "lef/lef.y"
                    {

    }
#line 4919 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 123: /* layer_option: K_ACCURRENTDENSITY layer_table_type int_number ';'  */
#line 1123 "lef/lef.y"
    {
      if (lefData->layerMastOver) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1527, "ACCURRENTDENSITY statement can't be defined in LAYER with TYPE MASTERSLICE or OVERLAP. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) {
           lefData->lefrLayer.addAccurrentDensity((yyvsp[-2].string));
           lefData->lefrLayer.setAcOneEntry((yyvsp[-1].dval));
      }
    }
#line 4938 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 124: /* layer_option: K_DCCURRENTDENSITY K_AVERAGE int_number ';'  */
#line 1138 "lef/lef.y"
    {
      if (lefData->layerMastOver) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1528, "DCCURRENTDENSITY statement can't be defined in LAYER with TYPE MASTERSLICE or OVERLAP. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) {
         lefData->lefrLayer.addDccurrentDensity("AVERAGE");
         lefData->lefrLayer.setDcOneEntry((yyvsp[-1].dval));
      }
    }
#line 4957 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 125: /* $@10: %empty  */
#line 1153 "lef/lef.y"
    {
      if (lefData->layerMastOver) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1528, "DCCURRENTDENSITY statement can't be defined in LAYER with TYPE MASTERSLICE or OVERLAP. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (!lefData->layerCut) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1529, "CUTAREA statement can only be defined in LAYER with type CUT. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) {
         lefData->lefrLayer.addDccurrentDensity("AVERAGE");
         lefData->lefrLayer.addNumber((yyvsp[0].dval));
      }
    }
#line 4984 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 126: /* $@11: %empty  */
#line 1176 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addDcCutarea(); }
#line 4990 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 127: /* layer_option: K_DCCURRENTDENSITY K_AVERAGE K_CUTAREA NUMBER $@10 number_list ';' $@11 dc_layer_table  */
#line 1177 "lef/lef.y"
                   {}
#line 4996 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 128: /* $@12: %empty  */
#line 1179 "lef/lef.y"
    {
      if (lefData->layerMastOver) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1528, "DCCURRENTDENSITY can't be defined in LAYER with TYPE MASTERSLICE or OVERLAP. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1530, "WIDTH statement can only be defined in LAYER with type ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) {
         lefData->lefrLayer.addDccurrentDensity("AVERAGE");
         lefData->lefrLayer.addNumber((yyvsp[0].dval));
      }
    }
#line 5023 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 129: /* $@13: %empty  */
#line 1202 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addDcWidth(); }
#line 5029 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 130: /* layer_option: K_DCCURRENTDENSITY K_AVERAGE K_WIDTH int_number $@12 int_number_list ';' $@13 dc_layer_table  */
#line 1203 "lef/lef.y"
                   {}
#line 5035 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 131: /* layer_option: K_ANTENNAAREARATIO int_number ';'  */
#line 1207 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "ANTENNAAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1531, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else if (lefData->use5_3) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "ANTENNADIFFAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNAAREARATIO syntax, which is incorrect.", lefData->versionNum);
               lefError(1704, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (!lefData->layerRout && !lefData->layerCut && lefData->layerMastOver) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1533, "ANTENNAAREARATIO statement can only be defined in LAYER with TYPE ROUTING or CUT. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setAntennaAreaRatio((yyvsp[-1].dval));
    }
#line 5077 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 132: /* $@14: %empty  */
#line 1245 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "ANTENNADIFFAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1532, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else if (lefData->use5_3) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "ANTENNADIFFAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNADIFFAREARATIO syntax, which is incorrect.", lefData->versionNum);
               lefError(1704, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (!lefData->layerRout && !lefData->layerCut && lefData->layerMastOver) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefError(1534, "ANTENNADIFFAREARATIO statement can only be defined in LAYER with TYPE ROUTING or CUT. Parser will stop processing.");
               CHKERR();
            }
         }
      }
      lefData->antennaType = lefiAntennaDAR; 
    }
#line 5119 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 133: /* layer_option: K_ANTENNADIFFAREARATIO $@14 layer_antenna_pwl ';'  */
#line 1282 "lef/lef.y"
                          {}
#line 5125 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 134: /* layer_option: K_ANTENNACUMAREARATIO int_number ';'  */
#line 1284 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "ANTENNACUMAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1535, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else if (lefData->use5_3) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "ANTENNACUMAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNACUMAREARATIO syntax, which is incorrect.", lefData->versionNum);
               lefError(1536, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (!lefData->layerRout && !lefData->layerCut && lefData->layerMastOver) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefError(1537, "ANTENNACUMAREARATIO statement can only be defined in LAYER with TYPE ROUTING or CUT. Parser will stop processing.");
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setAntennaCumAreaRatio((yyvsp[-1].dval));
    }
#line 5167 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 135: /* $@15: %empty  */
#line 1322 "lef/lef.y"
    {  // 5.4 syntax 
      lefData->use5_4 = 1;
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "ANTENNACUMDIFFAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1538, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else if (lefData->use5_3) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "ANTENNACUMDIFFAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNACUMDIFFAREARATIO syntax, which is incorrect.", lefData->versionNum);
               lefError(1539, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (!lefData->layerRout && !lefData->layerCut && lefData->layerMastOver) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1540, "ANTENNACUMDIFFAREARATIO statement can only be defined in LAYER with TYPE ROUTING or CUT. Parser will stop processing.");
              CHKERR();
            }
         }
      }
      lefData->antennaType = lefiAntennaCDAR;
    }
#line 5209 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 136: /* layer_option: K_ANTENNACUMDIFFAREARATIO $@15 layer_antenna_pwl ';'  */
#line 1359 "lef/lef.y"
                          {}
#line 5215 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 137: /* $@16: %empty  */
#line 1361 "lef/lef.y"
    { // both 5.3  & 5.4 syntax 
      if (!lefData->layerRout && !lefData->layerCut && lefData->layerMastOver) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefError(1541, "ANTENNAAREAFACTOR can only be defined in LAYER with TYPE ROUTING or CUT. Parser will stop processing.");
               CHKERR();
            }
         }
      }
      // this does not need to check, since syntax is in both 5.3 & 5.4 
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setAntennaAreaFactor((yyvsp[0].dval));
      lefData->antennaType = lefiAntennaAF;
    }
#line 5233 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 138: /* layer_option: K_ANTENNAAREAFACTOR int_number $@16 layer_antenna_duo ';'  */
#line 1374 "lef/lef.y"
                          {}
#line 5239 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 139: /* layer_option: K_ANTENNASIDEAREARATIO int_number ';'  */
#line 1376 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (!lefData->layerRout && (lefData->layerCut || lefData->layerMastOver)) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefError(1542, "ANTENNASIDEAREARATIO can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
               CHKERR();
            }
         }
      }
      if (lefData->ignoreVersion) {
        // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "ANTENNASIDEAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1543, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else if (lefData->use5_3) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "ANTENNASIDEAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNASIDEAREARATIO syntax, which is incorrect.", lefData->versionNum);
               lefError(1544, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setAntennaSideAreaRatio((yyvsp[-1].dval));
    }
#line 5281 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 140: /* $@17: %empty  */
#line 1414 "lef/lef.y"
    {  // 5.4 syntax 
      lefData->use5_4 = 1;
      if (!lefData->layerRout && (lefData->layerCut || lefData->layerMastOver)) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefError(1545, "ANTENNADIFFSIDEAREARATIO can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
               CHKERR();
            }
         }
      }
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "ANTENNADIFFSIDEAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1546, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else if (lefData->use5_3) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "ANTENNADIFFSIDEAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNADIFFSIDEAREARATIO syntax, which is incorrect.", lefData->versionNum);
               lefError(1547, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      lefData->antennaType = lefiAntennaDSAR;
    }
#line 5323 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 141: /* layer_option: K_ANTENNADIFFSIDEAREARATIO $@17 layer_antenna_pwl ';'  */
#line 1451 "lef/lef.y"
                          {}
#line 5329 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 142: /* layer_option: K_ANTENNACUMSIDEAREARATIO int_number ';'  */
#line 1453 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (!lefData->layerRout && (lefData->layerCut || lefData->layerMastOver)) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefError(1548, "ANTENNACUMSIDEAREARATIO can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
               CHKERR();
            }
         }
      }
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "ANTENNACUMSIDEAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1549, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else if (lefData->use5_3) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "ANTENNACUMSIDEAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNACUMSIDEAREARATIO syntax, which is incorrect.", lefData->versionNum);
               lefError(1550, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setAntennaCumSideAreaRatio((yyvsp[-1].dval));
    }
#line 5371 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 143: /* $@18: %empty  */
#line 1491 "lef/lef.y"
    {  // 5.4 syntax 
      lefData->use5_4 = 1;
      if (!lefData->layerRout && (lefData->layerCut || lefData->layerMastOver)) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefError(1551, "ANTENNACUMDIFFSIDEAREARATIO can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
               CHKERR();
            }
         }
      }
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "ANTENNACUMDIFFSIDEAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1552, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else if (lefData->use5_3) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "ANTENNACUMDIFFSIDEAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNACUMDIFFSIDEAREARATIO syntax, which is incorrect.", lefData->versionNum);
               lefError(1553, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      lefData->antennaType = lefiAntennaCDSAR;
    }
#line 5413 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 144: /* layer_option: K_ANTENNACUMDIFFSIDEAREARATIO $@18 layer_antenna_pwl ';'  */
#line 1528 "lef/lef.y"
                          {}
#line 5419 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 145: /* $@19: %empty  */
#line 1530 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (!lefData->layerRout && (lefData->layerCut || lefData->layerMastOver)) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefError(1554, "ANTENNASIDEAREAFACTOR can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
               CHKERR();
            }
         }
      }
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "ANTENNASIDEAREAFACTOR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1555, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else if (lefData->use5_3) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "ANTENNASIDEAREAFACTOR statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNASIDEAREAFACTOR syntax, which is incorrect.", lefData->versionNum);
               lefError(1556, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setAntennaSideAreaFactor((yyvsp[0].dval));
      lefData->antennaType = lefiAntennaSAF;
    }
#line 5462 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 146: /* layer_option: K_ANTENNASIDEAREAFACTOR int_number $@19 layer_antenna_duo ';'  */
#line 1568 "lef/lef.y"
                          {}
#line 5468 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 147: /* $@20: %empty  */
#line 1570 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (!lefData->layerRout && !lefData->layerCut && lefData->layerMastOver) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefError(1557, "ANTENNAMODEL can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
               CHKERR();
            }
         }
      }
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.5) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "ANTENNAMODEL statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1558, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else if (lefData->use5_3) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "ANTENNAMODEL statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNAMODEL syntax, which is incorrect.", lefData->versionNum);
               lefError(1559, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      lefData->antennaType = lefiAntennaO;
    }
#line 5510 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 148: /* layer_option: K_ANTENNAMODEL $@20 layer_oxide ';'  */
#line 1607 "lef/lef.y"
                    {}
#line 5516 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 149: /* layer_option: K_ANTENNACUMROUTINGPLUSCUT ';'  */
#line 1609 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
         lefData->outMsg = (char*)lefMalloc(10000);
         sprintf(lefData->outMsg,
           "ANTENNACUMROUTINGPLUSCUT is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
         lefError(1686, lefData->outMsg);
         lefFree(lefData->outMsg);
         CHKERR();
      } else {
         if (!lefData->layerRout && !lefData->layerCut) {
            if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
               if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                  lefError(1560, "ANTENNACUMROUTINGPLUSCUT can only be defined in LAYER with type ROUTING or CUT. Parser will stop processing.");
                  CHKERR();
               }
            }
         }
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setAntennaCumRoutingPlusCut();
      }
    }
#line 5541 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 150: /* layer_option: K_ANTENNAGATEPLUSDIFF int_number ';'  */
#line 1630 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
         lefData->outMsg = (char*)lefMalloc(10000);
         sprintf(lefData->outMsg,
           "ANTENNAGATEPLUSDIFF is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
         lefError(1687, lefData->outMsg);
         lefFree(lefData->outMsg);
         CHKERR();
      } else {
         if (!lefData->layerRout && !lefData->layerCut) {
            if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
               if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                  lefError(1561, "ANTENNAGATEPLUSDIFF can only be defined in LAYER with type ROUTING or CUT. Parser will stop processing.");
                  CHKERR();
               }
            }
         }
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setAntennaGatePlusDiff((yyvsp[-1].dval));
      }
    }
#line 5566 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 151: /* layer_option: K_ANTENNAAREAMINUSDIFF int_number ';'  */
#line 1651 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
         lefData->outMsg = (char*)lefMalloc(10000);
         sprintf(lefData->outMsg,
           "ANTENNAAREAMINUSDIFF is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
         lefError(1688, lefData->outMsg);
         lefFree(lefData->outMsg);
         CHKERR();
      } else {
         if (!lefData->layerRout && !lefData->layerCut) {
            if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
               if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                  lefError(1562, "ANTENNAAREAMINUSDIFF can only be defined in LAYER with type ROUTING or CUT. Parser will stop processing.");
                  CHKERR();
               }
            }
         }
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setAntennaAreaMinusDiff((yyvsp[-1].dval));
      }
    }
#line 5591 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 152: /* $@21: %empty  */
#line 1672 "lef/lef.y"
    {
      if (!lefData->layerRout && !lefData->layerCut) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefError(1563, "ANTENNAAREADIFFREDUCEPWL can only be defined in LAYER with type ROUTING or CUT. Parser will stop processing.");
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) { // require min 2 points, set the 1st 2 
         if (lefData->lefrAntennaPWLPtr) {
            lefData->lefrAntennaPWLPtr->Destroy();
            lefFree(lefData->lefrAntennaPWLPtr);
         }

         lefData->lefrAntennaPWLPtr = lefiAntennaPWL::create();
         lefData->lefrAntennaPWLPtr->addAntennaPWL((yyvsp[-1].pt).x, (yyvsp[-1].pt).y);
         lefData->lefrAntennaPWLPtr->addAntennaPWL((yyvsp[0].pt).x, (yyvsp[0].pt).y);
      }
    }
#line 5616 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 153: /* $@22: %empty  */
#line 1693 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
        lefData->lefrLayer.setAntennaPWL(lefiAntennaADR, lefData->lefrAntennaPWLPtr);
        lefData->lefrAntennaPWLPtr = NULL;
      }
    }
#line 5627 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 154: /* layer_option: K_ANTENNAAREADIFFREDUCEPWL '(' pt pt $@21 layer_diffusion_ratios ')' ';' $@22  */
#line 1699 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
        lefData->outMsg = (char*)lefMalloc(10000);
        sprintf(lefData->outMsg,
          "ANTENNAAREADIFFREDUCEPWL is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
        lefError(1689, lefData->outMsg);
        lefFree(lefData->outMsg);
        CHKERR();
      }
    }
#line 5642 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 155: /* layer_option: K_SLOTWIREWIDTH int_number ';'  */
#line 1710 "lef/lef.y"
    { // 5.4 syntax 
      if (lefData->ignoreVersion) {
         // do nothing 
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setSlotWireWidth((yyvsp[-1].dval));
      } else if (lefData->versionNum >= 5.7) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings)
               lefWarning(2011, "SLOTWIREWIDTH statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "SLOTWIREWIDTH statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1564, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setSlotWireWidth((yyvsp[-1].dval));
    }
#line 5670 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 156: /* layer_option: K_SLOTWIRELENGTH int_number ';'  */
#line 1734 "lef/lef.y"
    { // 5.4 syntax 
      if (lefData->ignoreVersion) {
         // do nothing 
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setSlotWireLength((yyvsp[-1].dval));
      } else if (lefData->versionNum >= 5.7) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings)
               lefWarning(2012, "SLOTWIRELENGTH statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "SLOTWIRELENGTH statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1565, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setSlotWireLength((yyvsp[-1].dval));
    }
#line 5698 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 157: /* layer_option: K_SLOTWIDTH int_number ';'  */
#line 1758 "lef/lef.y"
    { // 5.4 syntax 
      if (lefData->ignoreVersion) {
         // do nothing 
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setSlotWidth((yyvsp[-1].dval));
      } else if (lefData->versionNum >= 5.7) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings)
               lefWarning(2013, "SLOTWIDTH statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "SLOTWIDTH statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1566, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setSlotWidth((yyvsp[-1].dval));
    }
#line 5726 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 158: /* layer_option: K_SLOTLENGTH int_number ';'  */
#line 1782 "lef/lef.y"
    { // 5.4 syntax 
      if (lefData->ignoreVersion) {
         // do nothing 
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setSlotLength((yyvsp[-1].dval));
      } else if (lefData->versionNum >= 5.7) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings)
               lefWarning(2014, "SLOTLENGTH statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "SLOTLENGTH statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1567, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setSlotLength((yyvsp[-1].dval));
    }
#line 5754 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 159: /* layer_option: K_MAXADJACENTSLOTSPACING int_number ';'  */
#line 1806 "lef/lef.y"
    { // 5.4 syntax 
      if (lefData->ignoreVersion) {
         // do nothing 
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setMaxAdjacentSlotSpacing((yyvsp[-1].dval));
      } else if (lefData->versionNum >= 5.7) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings)
               lefWarning(2015, "MAXADJACENTSLOTSPACING statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "MAXADJACENTSLOTSPACING statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1568, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setMaxAdjacentSlotSpacing((yyvsp[-1].dval));
    }
#line 5782 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 160: /* layer_option: K_MAXCOAXIALSLOTSPACING int_number ';'  */
#line 1830 "lef/lef.y"
    { // 5.4 syntax 
      if (lefData->ignoreVersion) {
         // do nothing 
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setMaxCoaxialSlotSpacing((yyvsp[-1].dval));
      } else if (lefData->versionNum >= 5.7) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings)
                lefWarning(2016, "MAXCOAXIALSLOTSPACING statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "MAXCOAXIALSLOTSPACING statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1569, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setMaxCoaxialSlotSpacing((yyvsp[-1].dval));
    }
#line 5810 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 161: /* layer_option: K_MAXEDGESLOTSPACING int_number ';'  */
#line 1854 "lef/lef.y"
    { // 5.4 syntax 
      if (lefData->ignoreVersion) {
         // do nothing 
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setMaxEdgeSlotSpacing((yyvsp[-1].dval));
      } else if (lefData->versionNum >= 5.7) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings)
               lefWarning(2017, "MAXEDGESLOTSPACING statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "MAXEDGESLOTSPACING statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1570, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else
         if (lefCallbacks->LayerCbk) lefData->lefrLayer.setMaxEdgeSlotSpacing((yyvsp[-1].dval));
    }
#line 5838 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 162: /* layer_option: K_SPLITWIREWIDTH int_number ';'  */
#line 1878 "lef/lef.y"
    { // 5.4 syntax 
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum >= 5.7) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings)
               lefWarning(2018, "SPLITWIREWIDTH statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "SPLITWIREWIDTH statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1571, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setSplitWireWidth((yyvsp[-1].dval));
    }
#line 5865 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 163: /* layer_option: K_MINIMUMDENSITY int_number ';'  */
#line 1901 "lef/lef.y"
    { // 5.4 syntax, pcr 394389 
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "MINIMUMDENSITY statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1572, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setMinimumDensity((yyvsp[-1].dval));
    }
#line 5887 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 164: /* layer_option: K_MAXIMUMDENSITY int_number ';'  */
#line 1919 "lef/lef.y"
    { // 5.4 syntax, pcr 394389 
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "MAXIMUMDENSITY statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1573, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setMaximumDensity((yyvsp[-1].dval));
    }
#line 5909 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 165: /* layer_option: K_DENSITYCHECKWINDOW int_number int_number ';'  */
#line 1937 "lef/lef.y"
    { // 5.4 syntax, pcr 394389 
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "DENSITYCHECKWINDOW statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1574, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setDensityCheckWindow((yyvsp[-2].dval), (yyvsp[-1].dval));
    }
#line 5931 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 166: /* layer_option: K_DENSITYCHECKSTEP int_number ';'  */
#line 1955 "lef/lef.y"
    { // 5.4 syntax, pcr 394389 
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "DENSITYCHECKSTEP statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1575, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setDensityCheckStep((yyvsp[-1].dval));
    }
#line 5953 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 167: /* layer_option: K_FILLACTIVESPACING int_number ';'  */
#line 1973 "lef/lef.y"
    { // 5.4 syntax, pcr 394389 
      if (lefData->ignoreVersion) {
         // do nothing 
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "FILLACTIVESPACING statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1576, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setFillActiveSpacing((yyvsp[-1].dval));
    }
#line 5975 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 168: /* layer_option: K_MAXWIDTH int_number ';'  */
#line 1991 "lef/lef.y"
    {
      // 5.5 MAXWIDTH, is for routing layer only
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefError(1577, "MAXWIDTH statement can only be defined in LAYER with TYPE ROUTING.  Parser will stop processing.");
               CHKERR();
            }
         }
      }
      if (lefData->versionNum < 5.5) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "MAXWIDTH statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1578, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setMaxwidth((yyvsp[-1].dval));
    }
#line 6004 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 169: /* layer_option: K_MINWIDTH int_number ';'  */
#line 2016 "lef/lef.y"
    {
      // 5.5 MINWIDTH, is for routing layer only
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1579, "MINWIDTH statement can only be defined in LAYER with TYPE ROUTING.  Parser will stop processing.");
              CHKERR();
            }
         }
      }
      if (lefData->versionNum < 5.5) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "MINWIDTH statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1580, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setMinwidth((yyvsp[-1].dval));
    }
#line 6033 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 170: /* $@23: %empty  */
#line 2041 "lef/lef.y"
    {
      if (lefData->versionNum < 5.5) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "MINENCLOSEDAREA statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1581, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.addMinenclosedarea((yyvsp[0].dval));
    }
#line 6053 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 171: /* layer_option: K_MINENCLOSEDAREA NUMBER $@23 layer_minen_width ';'  */
#line 2056 "lef/lef.y"
                          {}
#line 6059 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 172: /* $@24: %empty  */
#line 2058 "lef/lef.y"
    { // pcr 409334 
      if (lefCallbacks->LayerCbk)
         lefData->lefrLayer.addMinimumcut((int)(yyvsp[-2].dval), (yyvsp[0].dval)); 
      lefData->hasLayerMincut = 0;
    }
#line 6069 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 173: /* layer_option: K_MINIMUMCUT int_number K_WIDTH int_number $@24 layer_minimumcut_within layer_minimumcut_from layer_minimumcut_length ';'  */
#line 2066 "lef/lef.y"
    {
      if (!lefData->hasLayerMincut) {   // FROMABOVE nor FROMBELOW is set 
         if (lefCallbacks->LayerCbk)
             lefData->lefrLayer.addMinimumcutConnect((char*)"");
      }
    }
#line 6080 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 174: /* $@25: %empty  */
#line 2073 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.addMinstep((yyvsp[0].dval));
    }
#line 6088 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 175: /* layer_option: K_MINSTEP int_number $@25 layer_minstep_options ';'  */
#line 2077 "lef/lef.y"
    {
    }
#line 6095 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 176: /* layer_option: K_PROTRUSIONWIDTH int_number K_LENGTH int_number K_WIDTH int_number ';'  */
#line 2080 "lef/lef.y"
    {
      if (lefData->versionNum < 5.5) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "PROTRUSION RULE statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1582, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.setProtrusion((yyvsp[-5].dval), (yyvsp[-3].dval), (yyvsp[-1].dval));
    }
#line 6115 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 177: /* $@26: %empty  */
#line 2096 "lef/lef.y"
    {
      if (lefData->versionNum < 5.5) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "SPACINGTABLE statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1583, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      // 5.5 either SPACING or SPACINGTABLE in a layer, not both
      if (lefData->lefrHasSpacing && lefData->layerRout && lefData->versionNum < 5.7) {
         if (lefCallbacks->LayerCbk)  // write warning only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefWarning(2010, "It is incorrect to have both SPACING rules & SPACINGTABLE rules within a ROUTING layer");
            }
      } 
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.addSpacingTable();
      lefData->lefrHasSpacingTbl = 1;
    }
#line 6143 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 178: /* layer_option: K_SPACINGTABLE $@26 sp_options ';'  */
#line 2119 "lef/lef.y"
                   {}
#line 6149 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 179: /* $@27: %empty  */
#line 2122 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "ENCLOSURE statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1584, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk)
         lefData->lefrLayer.addEnclosure((yyvsp[-2].string), (yyvsp[-1].dval), (yyvsp[0].dval));
    }
#line 6170 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 180: /* layer_option: K_ENCLOSURE layer_enclosure_type_opt int_number int_number $@27 layer_enclosure_width_opt ';'  */
#line 2138 "lef/lef.y"
                                  {}
#line 6176 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 181: /* $@28: %empty  */
#line 2141 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "PREFERENCLOSURE statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1585, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk)
         lefData->lefrLayer.addPreferEnclosure((yyvsp[-2].string), (yyvsp[-1].dval), (yyvsp[0].dval));
    }
#line 6197 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 182: /* layer_option: K_PREFERENCLOSURE layer_enclosure_type_opt int_number int_number $@28 layer_preferenclosure_width_opt ';'  */
#line 2157 "lef/lef.y"
                                        {}
#line 6203 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 183: /* layer_option: K_RESISTANCE int_number ';'  */
#line 2159 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "RESISTANCE statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1586, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else {
         if (lefCallbacks->LayerCbk)
            lefData->lefrLayer.setResPerCut((yyvsp[-1].dval));
      }
    }
#line 6225 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 184: /* layer_option: K_DIAGMINEDGELENGTH int_number ';'  */
#line 2177 "lef/lef.y"
    {
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1587, "DIAGMINEDGELENGTH can only be defined in LAYER with TYPE ROUTING. Parser will stop processing.");
              CHKERR();
            }
         }
      } else if (lefData->versionNum < 5.6) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "DIAGMINEDGELENGTH statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1588, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else {
         if (lefCallbacks->LayerCbk)
            lefData->lefrLayer.setDiagMinEdgeLength((yyvsp[-1].dval));
      }
    }
#line 6254 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 185: /* $@29: %empty  */
#line 2202 "lef/lef.y"
    {
      // Use the polygon code to retrieve the points for MINSIZE
      lefData->lefrGeometriesPtr = (lefiGeometries*)lefMalloc(sizeof(lefiGeometries));
      lefData->lefrGeometriesPtr->Init();
      lefData->lefrDoGeometries = 1;
    }
#line 6265 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 186: /* layer_option: K_MINSIZE $@29 firstPt otherPts ';'  */
#line 2209 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
         lefData->lefrGeometriesPtr->addPolygon();
         lefData->lefrLayer.setMinSize(lefData->lefrGeometriesPtr);
      }
     lefData->lefrDoGeometries = 0;
      lefData->lefrGeometriesPtr->Destroy();
      lefFree(lefData->lefrGeometriesPtr);
    }
#line 6279 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 188: /* layer_arraySpacing_long: K_LONGARRAY  */
#line 2222 "lef/lef.y"
    {
        if (lefCallbacks->LayerCbk)
           lefData->lefrLayer.setArraySpacingLongArray();
    }
#line 6288 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 190: /* layer_arraySpacing_width: K_WIDTH int_number  */
#line 2230 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
         lefData->lefrLayer.setArraySpacingWidth((yyvsp[0].dval));
    }
#line 6297 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 193: /* layer_arraySpacing_arraycut: K_ARRAYCUTS int_number K_SPACING int_number  */
#line 2241 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
         lefData->lefrLayer.addArraySpacingArray((int)(yyvsp[-2].dval), (yyvsp[0].dval));
         if (lefData->arrayCutsVal > (int)(yyvsp[-2].dval)) {
            // Mulitiple ARRAYCUTS value needs to me in ascending order 
            if (!lefData->arrayCutsWar) {
               if (lefData->layerWarnings++ < lefSettings->LayerWarnings)
                  lefWarning(2080, "The number of cut values in multiple ARRAYSPACING ARRAYCUTS are not in increasing order.\nTo be consistent with the documentation, update the cut values to increasing order.");
               lefData->arrayCutsWar = 1;
            }
         }
         lefData->arrayCutsVal = (int)(yyvsp[-2].dval);
    }
#line 6315 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 194: /* $@30: %empty  */
#line 2257 "lef/lef.y"
    { 
      if (lefData->hasInfluence) {  // 5.5 - INFLUENCE table must follow a PARALLEL
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1589, "An INFLUENCE table statement was defined before the PARALLELRUNLENGTH table statement.\nINFLUENCE table statement should be defined following the PARALLELRUNLENGTH.\nChange the LEF file and rerun the parser.");
              CHKERR();
            }
         }
      }
      if (lefData->hasParallel) { // 5.5 - Only one PARALLEL table is allowed per layer
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1590, "There is multiple PARALLELRUNLENGTH table statements are defined within a layer.\nAccording to the LEF Reference Manual, only one PARALLELRUNLENGTH table statement is allowed per layer.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.addNumber((yyvsp[0].dval));
      lefData->hasParallel = 1;
    }
#line 6340 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 195: /* $@31: %empty  */
#line 2278 "lef/lef.y"
    {
      lefData->spParallelLength = lefData->lefrLayer.getNumber();
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.addSpParallelLength();
    }
#line 6349 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 196: /* $@32: %empty  */
#line 2283 "lef/lef.y"
    { 
      if (lefCallbacks->LayerCbk) {
         lefData->lefrLayer.addSpParallelWidth((yyvsp[0].dval));
      }
    }
#line 6359 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 197: /* $@33: %empty  */
#line 2289 "lef/lef.y"
    { 
      if (lefData->lefrLayer.getNumber() != lefData->spParallelLength) {
         if (lefCallbacks->LayerCbk) {
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1591, "The number of length in the PARALLELRUNLENGTH statement is not equal to\nthe total number of spacings defined in the WIDTH statement in the SPACINGTABLE.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.addSpParallelWidthSpacing();
    }
#line 6375 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 199: /* $@34: %empty  */
#line 2303 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.addNumber((yyvsp[0].dval));
    }
#line 6383 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 200: /* $@35: %empty  */
#line 2307 "lef/lef.y"
    {
      if (lefData->hasParallel) { // 5.7 - Either PARALLEL OR TWOWIDTHS per layer
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1592, "A PARALLELRUNLENGTH statement was already defined in the layer.\nIt is PARALLELRUNLENGTH or TWOWIDTHS is allowed per layer.");
              CHKERR();
            }
         }
      }
      if (lefData->hasTwoWidths) { // 5.7 - only 1 TWOWIDTHS per layer
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1593, "A TWOWIDTHS table statement was already defined in the layer.\nOnly one TWOWIDTHS statement is allowed per layer.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.addSpTwoWidths((yyvsp[-4].dval), (yyvsp[-3].dval));
      lefData->hasTwoWidths = 1;
    }
#line 6408 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 201: /* sp_options: K_TWOWIDTHS K_WIDTH int_number layer_sp_TwoWidthsPRL int_number $@34 int_number_list $@35 layer_sp_TwoWidths  */
#line 2328 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
        lefData->outMsg = (char*)lefMalloc(10000);
        sprintf(lefData->outMsg,
          "TWOWIDTHS is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
        lefError(1697, lefData->outMsg);
        lefFree(lefData->outMsg);
        CHKERR();
      } 
    }
#line 6423 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 202: /* $@36: %empty  */
#line 2339 "lef/lef.y"
    {
      if (lefData->hasInfluence) {  // 5.5 - INFLUENCE table must follow a PARALLEL
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1594, "A INFLUENCE table statement was already defined in the layer.\nOnly one INFLUENCE statement is allowed per layer.");
              CHKERR();
            }
         }
      }
      if (!lefData->hasParallel) {  // 5.5 - INFLUENCE must follow a PARALLEL
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1595, "An INFLUENCE table statement was already defined before the layer.\nINFLUENCE statement has to be defined after the PARALLELRUNLENGTH table statement in the layer.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) {
         lefData->lefrLayer.setInfluence();
         lefData->lefrLayer.addSpInfluence((yyvsp[-4].dval), (yyvsp[-2].dval), (yyvsp[0].dval));
      }
    }
#line 6450 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 206: /* layer_spacingtable_opt: K_WITHIN int_number K_SPACING int_number  */
#line 2369 "lef/lef.y"
  {
    if (lefCallbacks->LayerCbk)
       lefData->lefrLayer.addSpacingTableOrthoWithin((yyvsp[-2].dval), (yyvsp[0].dval));
  }
#line 6459 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 207: /* layer_enclosure_type_opt: %empty  */
#line 2375 "lef/lef.y"
    {(yyval.string) = (char*)"NULL";}
#line 6465 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 208: /* layer_enclosure_type_opt: K_ABOVE  */
#line 2376 "lef/lef.y"
             {(yyval.string) = (char*)"ABOVE";}
#line 6471 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 209: /* layer_enclosure_type_opt: K_BELOW  */
#line 2377 "lef/lef.y"
             {(yyval.string) = (char*)"BELOW";}
#line 6477 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 211: /* $@37: %empty  */
#line 2381 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
         lefData->lefrLayer.addEnclosureWidth((yyvsp[0].dval));
      }
    }
#line 6487 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 213: /* layer_enclosure_width_opt: K_LENGTH int_number  */
#line 2388 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
         lefData->outMsg = (char*)lefMalloc(10000);
         sprintf(lefData->outMsg,
           "LENGTH is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
         lefError(1691, lefData->outMsg);
         lefFree(lefData->outMsg);
         CHKERR();
      } else {
         if (lefCallbacks->LayerCbk) {
            lefData->lefrLayer.addEnclosureLength((yyvsp[0].dval));
         }
      }
    }
#line 6506 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 215: /* layer_enclosure_width_except_opt: K_EXCEPTEXTRACUT int_number  */
#line 2405 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
         lefData->outMsg = (char*)lefMalloc(10000);
         sprintf(lefData->outMsg,
           "EXCEPTEXTRACUT is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
         lefError(1690, lefData->outMsg);
         lefFree(lefData->outMsg);
         CHKERR();
      } else {
         if (lefCallbacks->LayerCbk) {
            lefData->lefrLayer.addEnclosureExceptEC((yyvsp[0].dval));
         }
      }
    }
#line 6525 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 217: /* layer_preferenclosure_width_opt: K_WIDTH int_number  */
#line 2422 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
         lefData->lefrLayer.addPreferEnclosureWidth((yyvsp[0].dval));
      }
    }
#line 6535 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 219: /* layer_minimumcut_within: K_WITHIN int_number  */
#line 2430 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
        lefData->outMsg = (char*)lefMalloc(10000);
        sprintf(lefData->outMsg,
          "MINIMUMCUT WITHIN is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
        lefError(1700, lefData->outMsg);
        lefFree(lefData->outMsg);
        CHKERR();
      } else {
         if (lefCallbacks->LayerCbk) {
            lefData->lefrLayer.addMinimumcutWithin((yyvsp[0].dval));
         }
      }
    }
#line 6554 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 221: /* layer_minimumcut_from: K_FROMABOVE  */
#line 2447 "lef/lef.y"
    {
      if (lefData->versionNum < 5.5) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "FROMABOVE statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1596, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      lefData->hasLayerMincut = 1;
      if (lefCallbacks->LayerCbk)
         lefData->lefrLayer.addMinimumcutConnect((char*)"FROMABOVE");

    }
#line 6577 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 222: /* layer_minimumcut_from: K_FROMBELOW  */
#line 2466 "lef/lef.y"
    {
      if (lefData->versionNum < 5.5) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "FROMBELOW statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1597, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      }
      lefData->hasLayerMincut = 1;
      if (lefCallbacks->LayerCbk)
         lefData->lefrLayer.addMinimumcutConnect((char*)"FROMBELOW");
    }
#line 6599 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 224: /* layer_minimumcut_length: K_LENGTH int_number K_WITHIN int_number  */
#line 2486 "lef/lef.y"
    {   
      if (lefData->versionNum < 5.5) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "LENGTH WITHIN statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1598, lefData->outMsg);
               lefFree(lefData->outMsg);
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk)
         lefData->lefrLayer.addMinimumcutLengDis((yyvsp[-2].dval), (yyvsp[0].dval));
    }
#line 6620 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 227: /* layer_minstep_option: layer_minstep_type  */
#line 2508 "lef/lef.y"
  {
    if (lefCallbacks->LayerCbk) lefData->lefrLayer.addMinstepType((yyvsp[0].string));
  }
#line 6628 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 228: /* layer_minstep_option: K_LENGTHSUM int_number  */
#line 2512 "lef/lef.y"
  {
    if (lefCallbacks->LayerCbk) lefData->lefrLayer.addMinstepLengthsum((yyvsp[0].dval));
  }
#line 6636 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 229: /* layer_minstep_option: K_MAXEDGES int_number  */
#line 2516 "lef/lef.y"
  {
    if (lefData->versionNum < 5.7) {
      lefData->outMsg = (char*)lefMalloc(10000);
      sprintf(lefData->outMsg,
        "MAXEDGES is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
      lefError(1710, lefData->outMsg);
      lefFree(lefData->outMsg);
      CHKERR();
    } else
       if (lefCallbacks->LayerCbk) lefData->lefrLayer.addMinstepMaxedges((int)(yyvsp[0].dval));
  }
#line 6652 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 230: /* layer_minstep_type: K_INSIDECORNER  */
#line 2529 "lef/lef.y"
                 {(yyval.string) = (char*)"INSIDECORNER";}
#line 6658 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 231: /* layer_minstep_type: K_OUTSIDECORNER  */
#line 2530 "lef/lef.y"
                    {(yyval.string) = (char*)"OUTSIDECORNER";}
#line 6664 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 232: /* layer_minstep_type: K_STEP  */
#line 2531 "lef/lef.y"
           {(yyval.string) = (char*)"STEP";}
#line 6670 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 233: /* layer_antenna_pwl: int_number  */
#line 2535 "lef/lef.y"
      { if (lefCallbacks->LayerCbk)
          lefData->lefrLayer.setAntennaValue(lefData->antennaType, (yyvsp[0].dval)); }
#line 6677 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 234: /* $@38: %empty  */
#line 2538 "lef/lef.y"
      { if (lefCallbacks->LayerCbk) { // require min 2 points, set the 1st 2 
          if (lefData->lefrAntennaPWLPtr) {
            lefData->lefrAntennaPWLPtr->Destroy();
            lefFree(lefData->lefrAntennaPWLPtr);
          }

          lefData->lefrAntennaPWLPtr = lefiAntennaPWL::create();
          lefData->lefrAntennaPWLPtr->addAntennaPWL((yyvsp[-1].pt).x, (yyvsp[-1].pt).y);
          lefData->lefrAntennaPWLPtr->addAntennaPWL((yyvsp[0].pt).x, (yyvsp[0].pt).y);
        }
      }
#line 6693 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 235: /* layer_antenna_pwl: K_PWL '(' pt pt $@38 layer_diffusion_ratios ')'  */
#line 2550 "lef/lef.y"
      { 
        if (lefCallbacks->LayerCbk) {
          lefData->lefrLayer.setAntennaPWL(lefData->antennaType, lefData->lefrAntennaPWLPtr);
          lefData->lefrAntennaPWLPtr = NULL;
        }
      }
#line 6704 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 238: /* layer_diffusion_ratio: pt  */
#line 2563 "lef/lef.y"
  { if (lefCallbacks->LayerCbk)
      lefData->lefrAntennaPWLPtr->addAntennaPWL((yyvsp[0].pt).x, (yyvsp[0].pt).y);
  }
#line 6712 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 240: /* layer_antenna_duo: K_DIFFUSEONLY  */
#line 2569 "lef/lef.y"
      { 
        lefData->use5_4 = 1;
        if (lefData->ignoreVersion) {
           // do nothing 
        }
        else if ((lefData->antennaType == lefiAntennaAF) && (lefData->versionNum <= 5.3)) {
           if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
              if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                 lefData->outMsg = (char*)lefMalloc(10000);
                 sprintf (lefData->outMsg,
                   "ANTENNAAREAFACTOR with DIFFUSEONLY statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
                 lefError(1599, lefData->outMsg);
                 lefFree(lefData->outMsg);
                 CHKERR();
              }
           }
        } else if (lefData->use5_3) {
           if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
              if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                 lefData->outMsg = (char*)lefMalloc(10000);
                 sprintf (lefData->outMsg,
                   "ANTENNAAREAFACTOR with DIFFUSEONLY statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
                 lefError(1599, lefData->outMsg);
                 lefFree(lefData->outMsg);
                 CHKERR();
              }
           }
        }
        if (lefCallbacks->LayerCbk)
          lefData->lefrLayer.setAntennaDUO(lefData->antennaType);
      }
#line 6748 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 241: /* layer_table_type: K_PEAK  */
#line 2602 "lef/lef.y"
               {(yyval.string) = (char*)"PEAK";}
#line 6754 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 242: /* layer_table_type: K_AVERAGE  */
#line 2603 "lef/lef.y"
               {(yyval.string) = (char*)"AVERAGE";}
#line 6760 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 243: /* layer_table_type: K_RMS  */
#line 2604 "lef/lef.y"
               {(yyval.string) = (char*)"RMS";}
#line 6766 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 244: /* $@39: %empty  */
#line 2608 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addNumber((yyvsp[0].dval)); }
#line 6772 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 245: /* $@40: %empty  */
#line 2610 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addAcFrequency(); }
#line 6778 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 246: /* $@41: %empty  */
#line 2613 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addNumber((yyvsp[0].dval)); }
#line 6784 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 247: /* layer_frequency: K_FREQUENCY NUMBER $@39 number_list ';' $@40 ac_layer_table_opt K_TABLEENTRIES NUMBER $@41 number_list ';'  */
#line 2615 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addAcTableEntry(); }
#line 6790 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 249: /* $@42: %empty  */
#line 2619 "lef/lef.y"
    {
      if (!lefData->layerCut) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1600, "CUTAREA statement can only be defined in LAYER with TYPE CUT.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.addNumber((yyvsp[0].dval));
    }
#line 6806 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 250: /* ac_layer_table_opt: K_CUTAREA NUMBER $@42 number_list ';'  */
#line 2631 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addAcCutarea(); }
#line 6812 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 251: /* $@43: %empty  */
#line 2633 "lef/lef.y"
    {
      if (!lefData->layerRout) {
         if (lefCallbacks->LayerCbk) { // write error only if cbk is set 
            if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1601, "WIDTH can only be defined in LAYER with TYPE ROUTING.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->LayerCbk) lefData->lefrLayer.addNumber((yyvsp[0].dval));
    }
#line 6828 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 252: /* ac_layer_table_opt: K_WIDTH int_number $@43 int_number_list ';'  */
#line 2645 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addAcWidth(); }
#line 6834 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 253: /* $@44: %empty  */
#line 2649 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addNumber((yyvsp[0].dval)); }
#line 6840 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 254: /* dc_layer_table: K_TABLEENTRIES int_number $@44 int_number_list ';'  */
#line 2651 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addDcTableEntry(); }
#line 6846 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 256: /* int_number_list: int_number_list int_number  */
#line 2655 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addNumber((yyvsp[0].dval)); }
#line 6852 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 258: /* number_list: number_list NUMBER  */
#line 2659 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addNumber((yyvsp[0].dval)); }
#line 6858 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 261: /* layer_prop: T_STRING T_STRING  */
#line 2668 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
        char propTp;
        propTp = lefSettings->lefProps.lefrLayerProp.propType((yyvsp[-1].string));
        lefData->lefrLayer.addProp((yyvsp[-1].string), (yyvsp[0].string), propTp);
      }
    }
#line 6870 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 262: /* layer_prop: T_STRING QSTRING  */
#line 2676 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
        char propTp;
        propTp = lefSettings->lefProps.lefrLayerProp.propType((yyvsp[-1].string));
        lefData->lefrLayer.addProp((yyvsp[-1].string), (yyvsp[0].string), propTp);
      }
    }
#line 6882 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 263: /* layer_prop: T_STRING NUMBER  */
#line 2684 "lef/lef.y"
    {
      char temp[32];
      sprintf(temp, "%.11g", (yyvsp[0].dval));
      if (lefCallbacks->LayerCbk) {
        char propTp;
        propTp = lefSettings->lefProps.lefrLayerProp.propType((yyvsp[-1].string));
        lefData->lefrLayer.addNumProp((yyvsp[-1].string), (yyvsp[0].dval), temp, propTp);
      }
    }
#line 6896 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 264: /* current_density_pwl_list: current_density_pwl  */
#line 2696 "lef/lef.y"
    { }
#line 6902 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 265: /* current_density_pwl_list: current_density_pwl_list current_density_pwl  */
#line 2698 "lef/lef.y"
    { }
#line 6908 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 266: /* current_density_pwl: '(' int_number int_number ')'  */
#line 2701 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.setCurrentPoint((yyvsp[-2].dval), (yyvsp[-1].dval)); }
#line 6914 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 269: /* cap_point: '(' int_number int_number ')'  */
#line 2709 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.setCapacitancePoint((yyvsp[-2].dval), (yyvsp[-1].dval)); }
#line 6920 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 271: /* res_points: res_points res_point  */
#line 2714 "lef/lef.y"
    { }
#line 6926 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 272: /* res_point: '(' int_number int_number ')'  */
#line 2717 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.setResistancePoint((yyvsp[-2].dval), (yyvsp[-1].dval)); }
#line 6932 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 273: /* layer_type: K_ROUTING  */
#line 2720 "lef/lef.y"
                  {(yyval.string) = (char*)"ROUTING"; lefData->layerRout = 1;}
#line 6938 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 274: /* layer_type: K_CUT  */
#line 2721 "lef/lef.y"
                  {(yyval.string) = (char*)"CUT"; lefData->layerCut = 1;}
#line 6944 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 275: /* layer_type: K_OVERLAP  */
#line 2722 "lef/lef.y"
                  {(yyval.string) = (char*)"OVERLAP"; lefData->layerMastOver = 1;}
#line 6950 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 276: /* layer_type: K_MASTERSLICE  */
#line 2723 "lef/lef.y"
                  {(yyval.string) = (char*)"MASTERSLICE"; lefData->layerMastOver = 1;}
#line 6956 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 277: /* layer_type: K_VIRTUAL  */
#line 2724 "lef/lef.y"
                  {(yyval.string) = (char*)"VIRTUAL";}
#line 6962 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 278: /* layer_type: K_IMPLANT  */
#line 2725 "lef/lef.y"
                  {(yyval.string) = (char*)"IMPLANT";}
#line 6968 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 279: /* layer_direction: K_HORIZONTAL  */
#line 2728 "lef/lef.y"
                    {(yyval.string) = (char*)"HORIZONTAL";}
#line 6974 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 280: /* layer_direction: K_VERTICAL  */
#line 2729 "lef/lef.y"
                    {(yyval.string) = (char*)"VERTICAL";}
#line 6980 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 281: /* layer_direction: K_DIAG45  */
#line 2730 "lef/lef.y"
                    {(yyval.string) = (char*)"DIAG45";}
#line 6986 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 282: /* layer_direction: K_DIAG135  */
#line 2731 "lef/lef.y"
                    {(yyval.string) = (char*)"DIAG135";}
#line 6992 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 284: /* layer_minen_width: K_WIDTH int_number  */
#line 2735 "lef/lef.y"
    {
    if (lefCallbacks->LayerCbk)
       lefData->lefrLayer.addMinenclosedareaWidth((yyvsp[0].dval));
    }
#line 7001 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 285: /* layer_oxide: K_OXIDE1  */
#line 2742 "lef/lef.y"
    {
    if (lefCallbacks->LayerCbk)
       lefData->lefrLayer.addAntennaModel(1);
    }
#line 7010 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 286: /* layer_oxide: K_OXIDE2  */
#line 2747 "lef/lef.y"
    {
    if (lefCallbacks->LayerCbk)
       lefData->lefrLayer.addAntennaModel(2);
    }
#line 7019 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 287: /* layer_oxide: K_OXIDE3  */
#line 2752 "lef/lef.y"
    {
    if (lefCallbacks->LayerCbk)
       lefData->lefrLayer.addAntennaModel(3);
    }
#line 7028 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 288: /* layer_oxide: K_OXIDE4  */
#line 2757 "lef/lef.y"
    {
    if (lefCallbacks->LayerCbk)
       lefData->lefrLayer.addAntennaModel(4);
    }
#line 7037 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 289: /* layer_sp_parallel_widths: %empty  */
#line 2763 "lef/lef.y"
    { }
#line 7043 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 290: /* layer_sp_parallel_widths: layer_sp_parallel_widths layer_sp_parallel_width  */
#line 2765 "lef/lef.y"
    { }
#line 7049 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 291: /* $@45: %empty  */
#line 2768 "lef/lef.y"
    { 
      if (lefCallbacks->LayerCbk) {
         lefData->lefrLayer.addSpParallelWidth((yyvsp[0].dval));
      }
    }
#line 7059 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 292: /* layer_sp_parallel_width: K_WIDTH int_number $@45 int_number_list  */
#line 2774 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addSpParallelWidthSpacing(); }
#line 7065 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 293: /* layer_sp_TwoWidths: %empty  */
#line 2777 "lef/lef.y"
    { }
#line 7071 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 294: /* layer_sp_TwoWidths: layer_sp_TwoWidth layer_sp_TwoWidths  */
#line 2779 "lef/lef.y"
    { }
#line 7077 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 295: /* $@46: %empty  */
#line 2782 "lef/lef.y"
    {
       if (lefCallbacks->LayerCbk) lefData->lefrLayer.addNumber((yyvsp[0].dval));
    }
#line 7085 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 296: /* layer_sp_TwoWidth: K_WIDTH int_number layer_sp_TwoWidthsPRL int_number $@46 int_number_list  */
#line 2786 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
         lefData->lefrLayer.addSpTwoWidths((yyvsp[-4].dval), (yyvsp[-3].dval));
    }
#line 7094 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 297: /* layer_sp_TwoWidthsPRL: %empty  */
#line 2792 "lef/lef.y"
    { 
        (yyval.dval) = -1; // cannot use 0, since PRL number can be 0 
        lefData->lefrLayer.setSpTwoWidthsHasPRL(0);
    }
#line 7103 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 298: /* layer_sp_TwoWidthsPRL: K_PRL int_number  */
#line 2797 "lef/lef.y"
    { 
        (yyval.dval) = (yyvsp[0].dval); 
        lefData->lefrLayer.setSpTwoWidthsHasPRL(1);
    }
#line 7112 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 299: /* layer_sp_influence_widths: %empty  */
#line 2803 "lef/lef.y"
    { }
#line 7118 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 300: /* layer_sp_influence_widths: layer_sp_influence_widths layer_sp_influence_width  */
#line 2805 "lef/lef.y"
    { }
#line 7124 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 301: /* layer_sp_influence_width: K_WIDTH int_number K_WITHIN int_number K_SPACING int_number  */
#line 2808 "lef/lef.y"
    { if (lefCallbacks->LayerCbk) lefData->lefrLayer.addSpInfluence((yyvsp[-4].dval), (yyvsp[-2].dval), (yyvsp[0].dval)); }
#line 7130 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 302: /* maxstack_via: K_MAXVIASTACK int_number ';'  */
#line 2811 "lef/lef.y"
    {
      if (!lefData->lefrHasLayer) {  // 5.5 
        if (lefCallbacks->MaxStackViaCbk) { // write error only if cbk is set 
           if (lefData->maxStackViaWarnings++ < lefSettings->MaxStackViaWarnings) {
             lefError(1602, "MAXVIASTACK statement has to be defined after the LAYER statement.");
             CHKERR();
           }
        }
      } else if (lefData->lefrHasMaxVS) {
        if (lefCallbacks->MaxStackViaCbk) { // write error only if cbk is set 
           if (lefData->maxStackViaWarnings++ < lefSettings->MaxStackViaWarnings) {
             lefError(1603, "A MAXVIASTACK was already defined.\nOnly one MAXVIASTACK is allowed per lef file.");
             CHKERR();
           }
        }
      } else {
        if (lefCallbacks->MaxStackViaCbk) {
           lefData->lefrMaxStackVia.setMaxStackVia((int)(yyvsp[-1].dval));
           CALLBACK(lefCallbacks->MaxStackViaCbk, lefrMaxStackViaCbkType, &lefData->lefrMaxStackVia);
        }
      }
      if (lefData->versionNum < 5.5) {
        if (lefCallbacks->MaxStackViaCbk) { // write error only if cbk is set 
           if (lefData->maxStackViaWarnings++ < lefSettings->MaxStackViaWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                "MAXVIASTACK statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1604, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
      lefData->lefrHasMaxVS = 1;
    }
#line 7170 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 303: /* $@47: %empty  */
#line 2846 "lef/lef.y"
                                     {lefData->lefDumbMode = 2; lefData->lefNoNum= 2;}
#line 7176 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 304: /* maxstack_via: K_MAXVIASTACK int_number K_RANGE $@47 T_STRING T_STRING ';'  */
#line 2848 "lef/lef.y"
    {
      if (!lefData->lefrHasLayer) {  // 5.5 
        if (lefCallbacks->MaxStackViaCbk) { // write error only if cbk is set 
           if (lefData->maxStackViaWarnings++ < lefSettings->MaxStackViaWarnings) {
              lefError(1602, "MAXVIASTACK statement has to be defined after the LAYER statement.");
              CHKERR();
           }
        }
      } else if (lefData->lefrHasMaxVS) {
        if (lefCallbacks->MaxStackViaCbk) { // write error only if cbk is set 
           if (lefData->maxStackViaWarnings++ < lefSettings->MaxStackViaWarnings) {
             lefError(1603, "A MAXVIASTACK was already defined.\nOnly one MAXVIASTACK is allowed per lef file.");
             CHKERR();
           }
        }
      } else {
        if (lefCallbacks->MaxStackViaCbk) {
           lefData->lefrMaxStackVia.setMaxStackVia((int)(yyvsp[-5].dval));
           lefData->lefrMaxStackVia.setMaxStackViaRange((yyvsp[-2].string), (yyvsp[-1].string));
           CALLBACK(lefCallbacks->MaxStackViaCbk, lefrMaxStackViaCbkType, &lefData->lefrMaxStackVia);
        }
      }
      lefData->lefrHasMaxVS = 1;
    }
#line 7205 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 305: /* $@48: %empty  */
#line 2873 "lef/lef.y"
                { lefData->hasViaRule_layer = 0; }
#line 7211 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 306: /* via: start_via $@48 via_option end_via  */
#line 2874 "lef/lef.y"
    { 
      if (lefCallbacks->ViaCbk) {
        if (lefData->ndRule) 
            lefData->nd->addViaRule(&lefData->lefrVia);
         else 
            CALLBACK(lefCallbacks->ViaCbk, lefrViaCbkType, &lefData->lefrVia);
       }
    }
#line 7224 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 307: /* via_keyword: K_VIA  */
#line 2884 "lef/lef.y"
     { lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 7230 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 308: /* start_via: via_keyword T_STRING  */
#line 2887 "lef/lef.y"
    {
      // 0 is nodefault 
      if (lefCallbacks->ViaCbk) lefData->lefrVia.setName((yyvsp[0].string), 0);
      lefData->viaLayer = 0;
      lefData->numVia++;
      //strcpy(lefData->viaName, $2);
      lefData->viaName = strdup((yyvsp[0].string));
    }
#line 7243 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 309: /* start_via: via_keyword T_STRING K_DEFAULT  */
#line 2896 "lef/lef.y"
    {
      // 1 is default 
      if (lefCallbacks->ViaCbk) lefData->lefrVia.setName((yyvsp[-1].string), 1);
      lefData->viaLayer = 0;
      //strcpy(lefData->viaName, $2);
      lefData->viaName = strdup((yyvsp[-1].string));
    }
#line 7255 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 310: /* start_via: via_keyword T_STRING K_GENERATED  */
#line 2904 "lef/lef.y"
    {
      // 2 is generated 
      if (lefCallbacks->ViaCbk) lefData->lefrVia.setName((yyvsp[-1].string), 2);
      lefData->viaLayer = 0;
      //strcpy(lefData->viaName, $2);
      lefData->viaName = strdup((yyvsp[-1].string));
    }
#line 7267 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 311: /* $@49: %empty  */
#line 2912 "lef/lef.y"
                       {lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 7273 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 312: /* $@50: %empty  */
#line 2914 "lef/lef.y"
           {lefData->lefDumbMode = 3; lefData->lefNoNum = 1; }
#line 7279 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 313: /* $@51: %empty  */
#line 2917 "lef/lef.y"
    {
       if (lefData->versionNum < 5.6) {
         if (lefCallbacks->ViaCbk) { // write error only if cbk is set 
            if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                "VIARULE statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1709, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
            }
         }
       }  else
          if (lefCallbacks->ViaCbk) lefData->lefrVia.setViaRule((yyvsp[-21].string), (yyvsp[-18].dval), (yyvsp[-17].dval), (yyvsp[-13].string), (yyvsp[-12].string), (yyvsp[-11].string),
                          (yyvsp[-8].dval), (yyvsp[-7].dval), (yyvsp[-4].dval), (yyvsp[-3].dval), (yyvsp[-2].dval), (yyvsp[-1].dval));
       lefData->viaLayer++;
       lefData->hasViaRule_layer = 1;
    }
#line 7302 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 317: /* via_viarule_option: K_ROWCOL int_number int_number ';'  */
#line 2943 "lef/lef.y"
    {
       if (lefCallbacks->ViaCbk) lefData->lefrVia.setRowCol((int)(yyvsp[-2].dval), (int)(yyvsp[-1].dval));
    }
#line 7310 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 318: /* via_viarule_option: K_ORIGIN int_number int_number ';'  */
#line 2947 "lef/lef.y"
    {
       if (lefCallbacks->ViaCbk) lefData->lefrVia.setOrigin((yyvsp[-2].dval), (yyvsp[-1].dval));
    }
#line 7318 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 319: /* via_viarule_option: K_OFFSET int_number int_number int_number int_number ';'  */
#line 2951 "lef/lef.y"
    {
       if (lefCallbacks->ViaCbk) lefData->lefrVia.setOffset((yyvsp[-4].dval), (yyvsp[-3].dval), (yyvsp[-2].dval), (yyvsp[-1].dval));
    }
#line 7326 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 320: /* $@52: %empty  */
#line 2954 "lef/lef.y"
              {lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 7332 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 321: /* via_viarule_option: K_PATTERN $@52 T_STRING ';'  */
#line 2955 "lef/lef.y"
    {
       if (lefCallbacks->ViaCbk) lefData->lefrVia.setPattern((yyvsp[-1].string));
    }
#line 7340 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 327: /* via_other_option: via_foreign  */
#line 2972 "lef/lef.y"
    { }
#line 7346 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 328: /* via_other_option: via_layer_rule  */
#line 2974 "lef/lef.y"
    { }
#line 7352 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 329: /* via_other_option: K_RESISTANCE int_number ';'  */
#line 2976 "lef/lef.y"
    { if (lefCallbacks->ViaCbk) lefData->lefrVia.setResistance((yyvsp[-1].dval)); }
#line 7358 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 330: /* $@53: %empty  */
#line 2977 "lef/lef.y"
               { lefData->lefDumbMode = 1000000; }
#line 7364 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 331: /* via_other_option: K_PROPERTY $@53 via_prop_list ';'  */
#line 2978 "lef/lef.y"
    { lefData->lefDumbMode = 0;
    }
#line 7371 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 332: /* via_other_option: K_TOPOFSTACKONLY  */
#line 2981 "lef/lef.y"
    { 
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->ViaCbk) lefData->lefrVia.setTopOfStack();
      } else
        if (lefCallbacks->ViaCbk)  // write warning only if cbk is set 
           if (lefData->viaWarnings++ < lefSettings->ViaWarnings)
              lefWarning(2019, "TOPOFSTACKONLY statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later");
    }
#line 7384 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 335: /* via_name_value_pair: T_STRING NUMBER  */
#line 2997 "lef/lef.y"
    { 
      char temp[32];
      sprintf(temp, "%.11g", (yyvsp[0].dval));
      if (lefCallbacks->ViaCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrViaProp.propType((yyvsp[-1].string));
         lefData->lefrVia.addNumProp((yyvsp[-1].string), (yyvsp[0].dval), temp, propTp);
      }
    }
#line 7398 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 336: /* via_name_value_pair: T_STRING QSTRING  */
#line 3007 "lef/lef.y"
    {
      if (lefCallbacks->ViaCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrViaProp.propType((yyvsp[-1].string));
         lefData->lefrVia.addProp((yyvsp[-1].string), (yyvsp[0].string), propTp);
      }
    }
#line 7410 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 337: /* via_name_value_pair: T_STRING T_STRING  */
#line 3015 "lef/lef.y"
    {
      if (lefCallbacks->ViaCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrViaProp.propType((yyvsp[-1].string));
         lefData->lefrVia.addProp((yyvsp[-1].string), (yyvsp[0].string), propTp);
      }
    }
#line 7422 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 338: /* via_foreign: start_foreign ';'  */
#line 3025 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->ViaCbk) lefData->lefrVia.setForeign((yyvsp[-1].string), 0, 0.0, 0.0, -1);
      } else
        if (lefCallbacks->ViaCbk)  // write warning only if cbk is set 
           if (lefData->viaWarnings++ < lefSettings->ViaWarnings)
             lefWarning(2020, "FOREIGN statement in VIA is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 7435 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 339: /* via_foreign: start_foreign pt ';'  */
#line 3034 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->ViaCbk) lefData->lefrVia.setForeign((yyvsp[-2].string), 1, (yyvsp[-1].pt).x, (yyvsp[-1].pt).y, -1);
      } else
        if (lefCallbacks->ViaCbk)  // write warning only if cbk is set 
           if (lefData->viaWarnings++ < lefSettings->ViaWarnings)
             lefWarning(2020, "FOREIGN statement in VIA is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 7448 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 340: /* via_foreign: start_foreign pt orientation ';'  */
#line 3043 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->ViaCbk) lefData->lefrVia.setForeign((yyvsp[-3].string), 1, (yyvsp[-2].pt).x, (yyvsp[-2].pt).y, (yyvsp[-1].integer));
      } else
        if (lefCallbacks->ViaCbk)  // write warning only if cbk is set 
           if (lefData->viaWarnings++ < lefSettings->ViaWarnings)
             lefWarning(2020, "FOREIGN statement in VIA is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 7461 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 341: /* via_foreign: start_foreign orientation ';'  */
#line 3052 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->ViaCbk) lefData->lefrVia.setForeign((yyvsp[-2].string), 0, 0.0, 0.0, (yyvsp[-1].integer));
      } else
        if (lefCallbacks->ViaCbk)  // write warning only if cbk is set 
           if (lefData->viaWarnings++ < lefSettings->ViaWarnings)
             lefWarning(2020, "FOREIGN statement in VIA is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 7474 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 342: /* $@54: %empty  */
#line 3061 "lef/lef.y"
                                {lefData->lefDumbMode = 1; lefData->lefNoNum= 1;}
#line 7480 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 343: /* start_foreign: K_FOREIGN $@54 T_STRING  */
#line 3062 "lef/lef.y"
    { (yyval.string) = (yyvsp[0].string); }
#line 7486 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 344: /* orientation: K_N  */
#line 3065 "lef/lef.y"
              {(yyval.integer) = 0;}
#line 7492 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 345: /* orientation: K_W  */
#line 3066 "lef/lef.y"
              {(yyval.integer) = 1;}
#line 7498 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 346: /* orientation: K_S  */
#line 3067 "lef/lef.y"
              {(yyval.integer) = 2;}
#line 7504 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 347: /* orientation: K_E  */
#line 3068 "lef/lef.y"
              {(yyval.integer) = 3;}
#line 7510 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 348: /* orientation: K_FN  */
#line 3069 "lef/lef.y"
              {(yyval.integer) = 4;}
#line 7516 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 349: /* orientation: K_FW  */
#line 3070 "lef/lef.y"
              {(yyval.integer) = 5;}
#line 7522 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 350: /* orientation: K_FS  */
#line 3071 "lef/lef.y"
              {(yyval.integer) = 6;}
#line 7528 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 351: /* orientation: K_FE  */
#line 3072 "lef/lef.y"
              {(yyval.integer) = 7;}
#line 7534 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 352: /* orientation: K_R0  */
#line 3073 "lef/lef.y"
              {(yyval.integer) = 0;}
#line 7540 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 353: /* orientation: K_R90  */
#line 3074 "lef/lef.y"
              {(yyval.integer) = 1;}
#line 7546 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 354: /* orientation: K_R180  */
#line 3075 "lef/lef.y"
              {(yyval.integer) = 2;}
#line 7552 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 355: /* orientation: K_R270  */
#line 3076 "lef/lef.y"
              {(yyval.integer) = 3;}
#line 7558 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 356: /* orientation: K_MY  */
#line 3077 "lef/lef.y"
              {(yyval.integer) = 4;}
#line 7564 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 357: /* orientation: K_MYR90  */
#line 3078 "lef/lef.y"
              {(yyval.integer) = 5;}
#line 7570 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 358: /* orientation: K_MX  */
#line 3079 "lef/lef.y"
              {(yyval.integer) = 6;}
#line 7576 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 359: /* orientation: K_MXR90  */
#line 3080 "lef/lef.y"
              {(yyval.integer) = 7;}
#line 7582 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 360: /* via_layer_rule: via_layer via_geometries  */
#line 3083 "lef/lef.y"
    { }
#line 7588 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 361: /* $@55: %empty  */
#line 3085 "lef/lef.y"
                   {lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 7594 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 362: /* via_layer: K_LAYER $@55 T_STRING ';'  */
#line 3086 "lef/lef.y"
    {
      if (lefCallbacks->ViaCbk) lefData->lefrVia.addLayer((yyvsp[-1].string));
      lefData->viaLayer++;
      lefData->hasViaRule_layer = 1;
    }
#line 7604 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 365: /* via_geometry: K_RECT maskColor pt pt ';'  */
#line 3099 "lef/lef.y"
    { 
      if (lefCallbacks->ViaCbk) {
        if (lefData->versionNum < 5.8 && (int)(yyvsp[-3].integer) > 0) {
          if (lefData->viaWarnings++ < lefSettings->ViaWarnings) {
              lefError(2081, "MASK information can only be defined with version 5.8");
              CHKERR(); 
            }           
        } else {
          lefData->lefrVia.addRectToLayer((int)(yyvsp[-3].integer), (yyvsp[-2].pt).x, (yyvsp[-2].pt).y, (yyvsp[-1].pt).x, (yyvsp[-1].pt).y);
        }
      }
    }
#line 7621 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 366: /* $@56: %empty  */
#line 3112 "lef/lef.y"
    {
      lefData->lefrGeometriesPtr = (lefiGeometries*)lefMalloc(sizeof(lefiGeometries));
      lefData->lefrGeometriesPtr->Init();
      lefData->lefrDoGeometries = 1;
    }
#line 7631 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 367: /* via_geometry: K_POLYGON maskColor $@56 firstPt nextPt nextPt otherPts ';'  */
#line 3118 "lef/lef.y"
    { 
      if (lefCallbacks->ViaCbk) {
        if (lefData->versionNum < 5.8 && (yyvsp[-6].integer) > 0) {
          if (lefData->viaWarnings++ < lefSettings->ViaWarnings) {
              lefError(2083, "Color mask information can only be defined with version 5.8.");
              CHKERR(); 
            }           
        } else {
            lefData->lefrGeometriesPtr->addPolygon((int)(yyvsp[-6].integer));
            lefData->lefrVia.addPolyToLayer((int)(yyvsp[-6].integer), lefData->lefrGeometriesPtr);   // 5.6
        }
      }
      lefData->lefrGeometriesPtr->clearPolyItems(); // free items fields
      lefFree((char*)(lefData->lefrGeometriesPtr)); // Don't need anymore, poly data has
      lefData->lefrDoGeometries = 0;                // copied
    }
#line 7652 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 368: /* $@57: %empty  */
#line 3135 "lef/lef.y"
               {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 7658 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 369: /* end_via: K_END $@57 T_STRING  */
#line 3136 "lef/lef.y"
    { 
      // 10/17/2001 - Wanda da Rosa, PCR 404149
      //              Error if no layer in via
      if (!lefData->viaLayer) {
         if (lefCallbacks->ViaCbk) {  // write error only if cbk is set 
            if (lefData->viaWarnings++ < lefSettings->ViaWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                "A LAYER statement is missing in the VIA %s.\nAt least one LAYERis required per VIA statement.", (yyvsp[0].string));
              lefError(1606, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
            }
         }
      }
      if (strcmp(lefData->viaName, (yyvsp[0].string)) != 0) {
         if (lefCallbacks->ViaCbk) { // write error only if cbk is set 
            if (lefData->viaWarnings++ < lefSettings->ViaWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                "END VIA name %s is different from the VIA name %s.\nCorrect the LEF file before rerunning it through the LEF parser.", (yyvsp[0].string), lefData->viaName);
              lefError(1607, lefData->outMsg);
              lefFree(lefData->outMsg);
              lefFree(lefData->viaName);
              CHKERR();
            } else
              lefFree(lefData->viaName);
         } else
            lefFree(lefData->viaName);
      } else
         lefFree(lefData->viaName);
    }
#line 7695 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 370: /* $@58: %empty  */
#line 3169 "lef/lef.y"
                            { lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 7701 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 371: /* viarule_keyword: K_VIARULE $@58 T_STRING  */
#line 3170 "lef/lef.y"
    { 
      if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.setName((yyvsp[0].string));
      lefData->viaRuleLayer = 0;
      //strcpy(lefData->viaRuleName, $3);
      lefData->viaRuleName = strdup((yyvsp[0].string));
      lefData->isGenerate = 0;
    }
#line 7713 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 372: /* viarule: viarule_keyword viarule_layer_list via_names opt_viarule_props end_viarule  */
#line 3180 "lef/lef.y"
    {
      if (lefData->viaRuleLayer == 0 || lefData->viaRuleLayer > 2) {
         if (lefCallbacks->ViaRuleCbk) {  // write error only if cbk is set 
            if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
              lefError(1608, "A VIARULE statement requires two layers.");
              CHKERR();
            }
         }
      }
      if (lefCallbacks->ViaRuleCbk)
        CALLBACK(lefCallbacks->ViaRuleCbk, lefrViaRuleCbkType, &lefData->lefrViaRule);
      // 2/19/2004 - reset the ENCLOSURE overhang values which may be
      // set by the old syntax OVERHANG -- Not necessary, but just incase
      if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.clearLayerOverhang();
    }
#line 7733 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 373: /* $@59: %empty  */
#line 3198 "lef/lef.y"
    {
      lefData->isGenerate = 1;
    }
#line 7741 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 374: /* viarule_generate: viarule_keyword K_GENERATE viarule_generate_default $@59 viarule_layer_list opt_viarule_props end_viarule  */
#line 3202 "lef/lef.y"
    {
      if (lefData->viaRuleLayer == 0) {
         if (lefCallbacks->ViaRuleCbk) {  // write error only if cbk is set 
            if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
              lefError(1708, "A VIARULE GENERATE requires three layers.");
              CHKERR();
            }
         }
      } else if ((lefData->viaRuleLayer < 3) && (lefData->versionNum >= 5.6)) {
         if (lefCallbacks->ViaRuleCbk)  // write warning only if cbk is set 
            if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings)
              lefWarning(2021, "turn-via is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
      } else {
         if (lefCallbacks->ViaRuleCbk) {
            lefData->lefrViaRule.setGenerate();
            CALLBACK(lefCallbacks->ViaRuleCbk, lefrViaRuleCbkType, &lefData->lefrViaRule);
         }
      }
      // 2/19/2004 - reset the ENCLOSURE overhang values which may be
      // set by the old syntax OVERHANG
      if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.clearLayerOverhang();
    }
#line 7768 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 376: /* viarule_generate_default: K_DEFAULT  */
#line 3227 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
         if (lefCallbacks->ViaRuleCbk) {  // write error only if cbk is set 
            if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                "DEFAULT statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1605, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
            }
         }
      } else
        if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.setDefault();
    }
#line 7788 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 383: /* $@60: %empty  */
#line 3258 "lef/lef.y"
                         { lefData->lefDumbMode = 10000000;}
#line 7794 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 384: /* viarule_prop: K_PROPERTY $@60 viarule_prop_list ';'  */
#line 3259 "lef/lef.y"
    { lefData->lefDumbMode = 0;
    }
#line 7801 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 387: /* viarule_prop: T_STRING T_STRING  */
#line 3269 "lef/lef.y"
    {
      if (lefCallbacks->ViaRuleCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrViaRuleProp.propType((yyvsp[-1].string));
         lefData->lefrViaRule.addProp((yyvsp[-1].string), (yyvsp[0].string), propTp);
      }
    }
#line 7813 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 388: /* viarule_prop: T_STRING QSTRING  */
#line 3277 "lef/lef.y"
    {
      if (lefCallbacks->ViaRuleCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrViaRuleProp.propType((yyvsp[-1].string));
         lefData->lefrViaRule.addProp((yyvsp[-1].string), (yyvsp[0].string), propTp);
      }
    }
#line 7825 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 389: /* viarule_prop: T_STRING NUMBER  */
#line 3285 "lef/lef.y"
    {
      char temp[32];
      sprintf(temp, "%.11g", (yyvsp[0].dval));
      if (lefCallbacks->ViaRuleCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrViaRuleProp.propType((yyvsp[-1].string));
         lefData->lefrViaRule.addNumProp((yyvsp[-1].string), (yyvsp[0].dval), temp, propTp);
      }
    }
#line 7839 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 390: /* viarule_layer: viarule_layer_name viarule_layer_options  */
#line 3296 "lef/lef.y"
    {
      // 10/18/2001 - Wanda da Rosa PCR 404181
      //              Make sure the 1st 2 layers in viarule has direction
      // 04/28/2004 - PCR 704072 - DIRECTION in viarule generate is
      //              obsolete in 5.6
      if (lefData->versionNum >= 5.6) {
         if (lefData->viaRuleLayer < 2 && !lefData->viaRuleHasDir && !lefData->viaRuleHasEnc &&
             !lefData->isGenerate) {
            if (lefCallbacks->ViaRuleCbk) {  // write error only if cbk is set 
               if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
                  lefError(1705, "VIARULE statement in a layer, requires a DIRECTION construct statement.");
                  CHKERR(); 
               }
            }
         }
      } else {
         if (lefData->viaRuleLayer < 2 && !lefData->viaRuleHasDir && !lefData->viaRuleHasEnc &&
             lefData->isGenerate) {
            if (lefCallbacks->ViaRuleCbk) {  // write error only if cbk is set 
               if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
                  lefError(1705, "VIARULE statement in a layer, requires a DIRECTION construct statement.");
                  CHKERR(); 
               }
            }
         }
      }
      lefData->viaRuleLayer++;
    }
#line 7872 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 393: /* via_name: via_keyword T_STRING ';'  */
#line 3332 "lef/lef.y"
    { if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.addViaName((yyvsp[-1].string)); }
#line 7878 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 394: /* $@61: %empty  */
#line 3334 "lef/lef.y"
                            {lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 7884 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 395: /* viarule_layer_name: K_LAYER $@61 T_STRING ';'  */
#line 3335 "lef/lef.y"
    { if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.setLayer((yyvsp[-1].string));
      lefData->viaRuleHasDir = 0;
      lefData->viaRuleHasEnc = 0;
    }
#line 7893 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 398: /* viarule_layer_option: K_DIRECTION K_HORIZONTAL ';'  */
#line 3347 "lef/lef.y"
    {
      if (lefData->viaRuleHasEnc) {
        if (lefCallbacks->ViaRuleCbk) {  // write error only if cbk is set 
           if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
              lefError(1706, "An ENCLOSRE statement was already defined in the layer.\nIt is DIRECTION or ENCLOSURE can be specified in a layer.");
              CHKERR();
           }
        }
      } else {
        if ((lefData->versionNum < 5.6) || (!lefData->isGenerate)) {
          if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.setHorizontal();
        } else
          if (lefCallbacks->ViaRuleCbk)  // write warning only if cbk is set 
             if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings)
               lefWarning(2022, "DIRECTION statement in VIARULE is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
      }
      lefData->viaRuleHasDir = 1;
    }
#line 7916 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 399: /* viarule_layer_option: K_DIRECTION K_VERTICAL ';'  */
#line 3366 "lef/lef.y"
    { 
      if (lefData->viaRuleHasEnc) {
        if (lefCallbacks->ViaRuleCbk) { // write error only if cbk is set 
           if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
              lefError(1706, "An ENCLOSRE statement was already defined in the layer.\nIt is DIRECTION or ENCLOSURE can be specified in a layer.");
              CHKERR();
           }
        }
      } else {
        if ((lefData->versionNum < 5.6) || (!lefData->isGenerate)) {
          if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.setVertical();
        } else
          if (lefCallbacks->ViaRuleCbk) // write warning only if cbk is set 
           if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings)
              lefWarning(2022, "DIRECTION statement in VIARULE is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
      }
      lefData->viaRuleHasDir = 1;
    }
#line 7939 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 400: /* viarule_layer_option: K_ENCLOSURE int_number int_number ';'  */
#line 3385 "lef/lef.y"
    {
      if (lefData->versionNum < 5.5) {
         if (lefCallbacks->ViaRuleCbk) { // write error only if cbk is set 
           if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                "ENCLOSURE statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1707, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
         }
      }
      // 2/19/2004 - Enforced the rule that ENCLOSURE can only be defined
      // in VIARULE GENERATE
      if (!lefData->isGenerate) {
         if (lefCallbacks->ViaRuleCbk) { // write error only if cbk is set 
           if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
              lefError(1614, "An ENCLOSURE statement is defined in a VIARULE statement only.\nOVERHANG statement can only be defined in VIARULE GENERATE.");
              CHKERR();
           }
         }
      }
      if (lefData->viaRuleHasDir) {
         if (lefCallbacks->ViaRuleCbk) { // write error only if cbk is set 
           if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
              lefError(1609, "A DIRECTION statement was already defined in the layer.\nIt is DIRECTION or ENCLOSURE can be specified in a layer.");
              CHKERR();
           }
         }
      } else {
         if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.setEnclosure((yyvsp[-2].dval), (yyvsp[-1].dval));
      }
      lefData->viaRuleHasEnc = 1;
    }
#line 7979 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 401: /* viarule_layer_option: K_WIDTH int_number K_TO int_number ';'  */
#line 3421 "lef/lef.y"
    { if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.setWidth((yyvsp[-3].dval),(yyvsp[-1].dval)); }
#line 7985 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 402: /* viarule_layer_option: K_RECT pt pt ';'  */
#line 3423 "lef/lef.y"
    { if (lefCallbacks->ViaRuleCbk)
        lefData->lefrViaRule.setRect((yyvsp[-2].pt).x, (yyvsp[-2].pt).y, (yyvsp[-1].pt).x, (yyvsp[-1].pt).y); }
#line 7992 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 403: /* viarule_layer_option: K_SPACING int_number K_BY int_number ';'  */
#line 3426 "lef/lef.y"
    { if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.setSpacing((yyvsp[-3].dval),(yyvsp[-1].dval)); }
#line 7998 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 404: /* viarule_layer_option: K_RESISTANCE int_number ';'  */
#line 3428 "lef/lef.y"
    { if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.setResistance((yyvsp[-1].dval)); }
#line 8004 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 405: /* viarule_layer_option: K_OVERHANG int_number ';'  */
#line 3430 "lef/lef.y"
    {
      if (!lefData->viaRuleHasDir) {
         if (lefCallbacks->ViaRuleCbk) {  // write error only if cbk is set 
            if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
               lefError(1610, "An OVERHANG statement is defined, but the required DIRECTION statement is not yet defined.\nUpdate the LEF file to define the DIRECTION statement before the OVERHANG.");
               CHKERR();
            }
         }
      }
      // 2/19/2004 - Enforced the rule that OVERHANG can only be defined
      // in VIARULE GENERATE after 5.3
      if ((lefData->versionNum > 5.3) && (!lefData->isGenerate)) {
         if (lefCallbacks->ViaRuleCbk) {  // write error only if cbk is set 
            if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
               lefError(1611, "An OVERHANG statement is defined in a VIARULE statement only.\nOVERHANG statement can only be defined in VIARULE GENERATE.");
               CHKERR();
            }
         }
      }
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.setOverhang((yyvsp[-1].dval));
      } else {
        if (lefCallbacks->ViaRuleCbk)  // write warning only if cbk is set 
           if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings)
              lefWarning(2023, "OVERHANG statement will be translated into similar ENCLOSURE rule");
        // In 5.6 & later, set it to either ENCLOSURE overhang1 or overhang2
        if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.setOverhangToEnclosure((yyvsp[-1].dval));
      }
    }
#line 8038 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 406: /* viarule_layer_option: K_METALOVERHANG int_number ';'  */
#line 3460 "lef/lef.y"
    {
      // 2/19/2004 - Enforced the rule that METALOVERHANG can only be defined
      // in VIARULE GENERATE
      if ((lefData->versionNum > 5.3) && (!lefData->isGenerate)) {
         if (lefCallbacks->ViaRuleCbk) {  // write error only if cbk is set 
            if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
               lefError(1612, "An METALOVERHANG statement is defined in a VIARULE statement only.\nOVERHANG statement can only be defined in VIARULE GENERATE.");
               CHKERR();
            }
         }
      }
      if (lefData->versionNum < 5.6) {
        if (!lefData->viaRuleHasDir) {
           if (lefCallbacks->ViaRuleCbk) {  // write error only if cbk is set 
             if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
                lefError(1613, "An METALOVERHANG statement is defined, but the required DIRECTION statement is not yet defined.\nUpdate the LEF file to define the DIRECTION statement before the OVERHANG.");
                CHKERR();
             } 
           }
        }
        if (lefCallbacks->ViaRuleCbk) lefData->lefrViaRule.setMetalOverhang((yyvsp[-1].dval));
      } else
        if (lefCallbacks->ViaRuleCbk)  // write warning only if cbk is set 
           if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings)
             lefWarning(2024, "METALOVERHANG statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 8069 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 407: /* $@62: %empty  */
#line 3487 "lef/lef.y"
                   {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 8075 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 408: /* end_viarule: K_END $@62 T_STRING  */
#line 3488 "lef/lef.y"
    {
      if ((lefData->isGenerate) && (lefCallbacks->ViaRuleCbk) && lefData->lefrViaRule.numLayers() >= 3) {         
        if (!lefData->lefrViaRule.layer(0)->hasRect() &&
            !lefData->lefrViaRule.layer(1)->hasRect() &&
            !lefData->lefrViaRule.layer(2)->hasRect()) {
            lefData->outMsg = (char*)lefMalloc(10000);
            sprintf (lefData->outMsg, 
                     "VIARULE GENERATE '%s' cut layer definition should have RECT statement.\nCorrect the LEF file before rerunning it through the LEF parser.", 
                      lefData->viaRuleName);
            lefWarning(1714, lefData->outMsg); 
            lefFree(lefData->outMsg);            
            CHKERR();                
        }
      }

      if (strcmp(lefData->viaRuleName, (yyvsp[0].string)) != 0) {
        if (lefCallbacks->ViaRuleCbk) {  // write error only if cbk is set 
           if (lefData->viaRuleWarnings++ < lefSettings->ViaRuleWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "END VIARULE name %s is different from the VIARULE name %s.\nCorrect the LEF file before rerunning it through the LEF parser.", (yyvsp[0].string), lefData->viaRuleName);
              lefError(1615, lefData->outMsg);
              lefFree(lefData->outMsg);
              lefFree(lefData->viaRuleName);
              CHKERR();
           } else
              lefFree(lefData->viaRuleName);
        } else
           lefFree(lefData->viaRuleName);
      } else
        lefFree(lefData->viaRuleName);
    }
#line 8112 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 409: /* spacing_rule: start_spacing spacings end_spacing  */
#line 3522 "lef/lef.y"
    { }
#line 8118 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 410: /* start_spacing: K_SPACING  */
#line 3525 "lef/lef.y"
    {
      lefData->hasSamenet = 0;
      if ((lefData->versionNum < 5.6) || (!lefData->ndRule)) {
        // if 5.6 and in nondefaultrule, it should not get in here, 
        // it should go to the else statement to write out a warning 
        // if 5.6, not in nondefaultrule, it will get in here 
        // if 5.5 and earlier in nondefaultrule is ok to get in here 
        if (lefData->versionNum >= 5.7) { // will get to this if statement if  
                           // lefData->versionNum is 5.6 and higher but lefData->ndRule = 0 
           if (lefData->spacingWarnings == 0) {  // only print once 
              lefWarning(2077, "A SPACING SAMENET section is defined but it is not legal in a LEF 5.7 version file.\nIt will be ignored which will probably cause real DRC violations to be ignored, and may\ncause false DRC violations to occur.\n\nTo avoid this warning, and correctly handle these DRC rules, you should modify your\nLEF to use the appropriate SAMENET keywords as described in the LEF/DEF 5.7\nmanual under the SPACING statements in the LAYER (Routing) and LAYER (Cut)\nsections listed in the LEF Table of Contents.");
              lefData->spacingWarnings++;
           }
        } else if (lefCallbacks->SpacingBeginCbk && !lefData->ndRule)
          CALLBACK(lefCallbacks->SpacingBeginCbk, lefrSpacingBeginCbkType, 0);
      } else
        if (lefCallbacks->SpacingBeginCbk && !lefData->ndRule)  // write warning only if cbk is set 
           if (lefData->spacingWarnings++ < lefSettings->SpacingWarnings)
             lefWarning(2025, "SAMENET statement in NONDEFAULTRULE is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 8143 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 411: /* end_spacing: K_END K_SPACING  */
#line 3547 "lef/lef.y"
    {
      if ((lefData->versionNum < 5.6) || (!lefData->ndRule)) {
        if ((lefData->versionNum <= 5.4) && (!lefData->hasSamenet)) {
           lefError(1616, "SAMENET statement is required inside SPACING for any lef file with version 5.4 and earlier, but is not defined in the parsed lef file.");
           CHKERR();
        } else if (lefData->versionNum < 5.7) { // obsolete in 5.7 and later 
           if (lefCallbacks->SpacingEndCbk && !lefData->ndRule)
             CALLBACK(lefCallbacks->SpacingEndCbk, lefrSpacingEndCbkType, 0);
        }
      }
    }
#line 8159 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 414: /* spacing: samenet_keyword T_STRING T_STRING int_number ';'  */
#line 3565 "lef/lef.y"
    {
      if ((lefData->versionNum < 5.6) || (!lefData->ndRule)) {
        if (lefData->versionNum < 5.7) {
          if (lefCallbacks->SpacingCbk) {
            lefData->lefrSpacing.set((yyvsp[-3].string), (yyvsp[-2].string), (yyvsp[-1].dval), 0);
            if (lefData->ndRule)
                lefData->nd->addSpacingRule(&lefData->lefrSpacing);
            else 
                CALLBACK(lefCallbacks->SpacingCbk, lefrSpacingCbkType, &lefData->lefrSpacing);            
          }
        }
      }
    }
#line 8177 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 415: /* spacing: samenet_keyword T_STRING T_STRING int_number K_STACK ';'  */
#line 3579 "lef/lef.y"
    {
      if ((lefData->versionNum < 5.6) || (!lefData->ndRule)) {
        if (lefData->versionNum < 5.7) {
          if (lefCallbacks->SpacingCbk) {
            lefData->lefrSpacing.set((yyvsp[-4].string), (yyvsp[-3].string), (yyvsp[-2].dval), 1);
            if (lefData->ndRule)
                lefData->nd->addSpacingRule(&lefData->lefrSpacing);
            else 
                CALLBACK(lefCallbacks->SpacingCbk, lefrSpacingCbkType, &lefData->lefrSpacing);    
          }
        }
      }
    }
#line 8195 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 416: /* samenet_keyword: K_SAMENET  */
#line 3595 "lef/lef.y"
    { lefData->lefDumbMode = 2; lefData->lefNoNum = 2; lefData->hasSamenet = 1; }
#line 8201 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 417: /* maskColor: %empty  */
#line 3599 "lef/lef.y"
    { (yyval.integer) = 0; }
#line 8207 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 418: /* maskColor: K_MASK int_number  */
#line 3601 "lef/lef.y"
    { (yyval.integer) = (int)(yyvsp[0].dval); }
#line 8213 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 419: /* irdrop: start_irdrop ir_tables end_irdrop  */
#line 3604 "lef/lef.y"
    { }
#line 8219 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 420: /* start_irdrop: K_IRDROP  */
#line 3607 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->IRDropBeginCbk) 
          CALLBACK(lefCallbacks->IRDropBeginCbk, lefrIRDropBeginCbkType, 0);
      } else
        if (lefCallbacks->IRDropBeginCbk) // write warning only if cbk is set 
          if (lefData->iRDropWarnings++ < lefSettings->IRDropWarnings)
            lefWarning(2026, "IRDROP statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 8233 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 421: /* end_irdrop: K_END K_IRDROP  */
#line 3618 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->IRDropEndCbk)
          CALLBACK(lefCallbacks->IRDropEndCbk, lefrIRDropEndCbkType, 0);
      }
    }
#line 8244 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 424: /* ir_table: ir_tablename ir_table_values ';'  */
#line 3632 "lef/lef.y"
    { 
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->IRDropCbk)
          CALLBACK(lefCallbacks->IRDropCbk, lefrIRDropCbkType, &lefData->lefrIRDrop);
      }
    }
#line 8255 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 427: /* ir_table_value: int_number int_number  */
#line 3645 "lef/lef.y"
  { if (lefCallbacks->IRDropCbk) lefData->lefrIRDrop.setValues((yyvsp[-1].dval), (yyvsp[0].dval)); }
#line 8261 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 428: /* ir_tablename: K_TABLE T_STRING  */
#line 3648 "lef/lef.y"
  { if (lefCallbacks->IRDropCbk) lefData->lefrIRDrop.setTableName((yyvsp[0].string)); }
#line 8267 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 429: /* minfeature: K_MINFEATURE int_number int_number ';'  */
#line 3651 "lef/lef.y"
  {
    lefData->hasMinfeature = 1;
    if (lefData->versionNum < 5.4) {
       if (lefCallbacks->MinFeatureCbk) {
         lefData->lefrMinFeature.set((yyvsp[-2].dval), (yyvsp[-1].dval));
         CALLBACK(lefCallbacks->MinFeatureCbk, lefrMinFeatureCbkType, &lefData->lefrMinFeature);
       }
    } else
       if (lefCallbacks->MinFeatureCbk) // write warning only if cbk is set 
          if (lefData->minFeatureWarnings++ < lefSettings->MinFeatureWarnings)
            lefWarning(2027, "MINFEATURE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
  }
#line 8284 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 430: /* dielectric: K_DIELECTRIC int_number ';'  */
#line 3665 "lef/lef.y"
  {
    if (lefData->versionNum < 5.4) {
       if (lefCallbacks->DielectricCbk)
         CALLBACK(lefCallbacks->DielectricCbk, lefrDielectricCbkType, (yyvsp[-1].dval));
    } else
       if (lefCallbacks->DielectricCbk) // write warning only if cbk is set 
         if (lefData->dielectricWarnings++ < lefSettings->DielectricWarnings)
           lefWarning(2028, "DIELECTRIC statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
  }
#line 8298 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 431: /* $@63: %empty  */
#line 3675 "lef/lef.y"
                                  {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 8304 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 432: /* $@64: %empty  */
#line 3676 "lef/lef.y"
  {
    (void)lefSetNonDefault((yyvsp[0].string));
    if (lefCallbacks->NonDefaultCbk) lefData->lefrNonDefault.setName((yyvsp[0].string));
    lefData->ndLayer = 0;
    lefData->ndRule = 1;
    lefData->numVia = 0;
    //strcpy(lefData->nonDefaultRuleName, $3);
    lefData->nonDefaultRuleName = strdup((yyvsp[0].string));
  }
#line 8318 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 433: /* $@65: %empty  */
#line 3686 "lef/lef.y"
           {lefData->lefNdRule = 1;}
#line 8324 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 434: /* nondefault_rule: K_NONDEFAULTRULE $@63 T_STRING $@64 nd_hardspacing nd_rules $@65 end_nd_rule  */
#line 3687 "lef/lef.y"
  {
    // 10/18/2001 - Wanda da Rosa, PCR 404189
    //              At least 1 layer is required
    if ((!lefData->ndLayer) && (!lefSettings->RelaxMode)) {
       if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
         if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings) {
            lefError(1617, "NONDEFAULTRULE statement requires at least one LAYER statement.");
            CHKERR();
         }
       }
    }
    if ((!lefData->numVia) && (!lefSettings->RelaxMode) && (lefData->versionNum < 5.6)) {
       // VIA is no longer a required statement in 5.6
       if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
         if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings) {
            lefError(1618, "NONDEFAULTRULE statement requires at least one VIA statement.");
            CHKERR();
         }
       }
    }
    if (lefCallbacks->NonDefaultCbk) {
      lefData->lefrNonDefault.end();
      CALLBACK(lefCallbacks->NonDefaultCbk, lefrNonDefaultCbkType, &lefData->lefrNonDefault);
    }
    lefData->ndRule = 0;
    lefData->lefDumbMode = 0;
    (void)lefUnsetNonDefault();
  }
#line 8357 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 435: /* end_nd_rule: K_END  */
#line 3717 "lef/lef.y"
    {
      if ((lefData->nonDefaultRuleName) && (*lefData->nonDefaultRuleName != '\0'))
        lefFree(lefData->nonDefaultRuleName);
    }
#line 8366 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 436: /* end_nd_rule: K_END T_STRING  */
#line 3722 "lef/lef.y"
    {
      if (strcmp(lefData->nonDefaultRuleName, (yyvsp[0].string)) != 0) {
        if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
          if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings) {
             lefData->outMsg = (char*)lefMalloc(10000);
             sprintf (lefData->outMsg,
                "END NONDEFAULTRULE name %s is different from the NONDEFAULTRULE name %s.\nCorrect the LEF file before rerunning it through the LEF parser.", (yyvsp[0].string), lefData->nonDefaultRuleName);
             lefError(1619, lefData->outMsg);
             lefFree(lefData->nonDefaultRuleName);
             lefFree(lefData->outMsg);
             CHKERR();
          } else
             lefFree(lefData->nonDefaultRuleName);
        } else
           lefFree(lefData->nonDefaultRuleName);
      } else
        lefFree(lefData->nonDefaultRuleName);
    }
#line 8389 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 438: /* nd_hardspacing: K_HARDSPACING ';'  */
#line 3745 "lef/lef.y"
    {
       if (lefData->versionNum < 5.6) {
          if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
            if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "HARDSPACING statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1620, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
          }
       } else 
          if (lefCallbacks->NonDefaultCbk)
             lefData->lefrNonDefault.setHardspacing();
    }
#line 8410 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 448: /* usevia: K_USEVIA T_STRING ';'  */
#line 3778 "lef/lef.y"
    {
       if (lefData->versionNum < 5.6) {
          if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
             lefData->outMsg = (char*)lefMalloc(10000);
             sprintf (lefData->outMsg,
               "USEVIA statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
             lefError(1621, lefData->outMsg);
             lefFree(lefData->outMsg);
             CHKERR();
          }
       } else {
          if (lefCallbacks->NonDefaultCbk)
             lefData->lefrNonDefault.addUseVia((yyvsp[-1].string));
       }
    }
#line 8430 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 449: /* useviarule: K_USEVIARULE T_STRING ';'  */
#line 3795 "lef/lef.y"
    {
       if (lefData->versionNum < 5.6) {
          if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
             if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings) {
                lefData->outMsg = (char*)lefMalloc(10000);
                sprintf (lefData->outMsg,
                  "USEVIARULE statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
                lefError(1622, lefData->outMsg);
                lefFree(lefData->outMsg);
                CHKERR();
             }
          }
       } else {
          if (lefCallbacks->NonDefaultCbk)
             lefData->lefrNonDefault.addUseViaRule((yyvsp[-1].string));
       }
    }
#line 8452 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 450: /* mincuts: K_MINCUTS T_STRING int_number ';'  */
#line 3814 "lef/lef.y"
    {
       if (lefData->versionNum < 5.6) {
          if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
             if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings) {
                lefData->outMsg = (char*)lefMalloc(10000);
                sprintf (lefData->outMsg,
                  "MINCUTS statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
                lefError(1623, lefData->outMsg);
                lefFree(lefData->outMsg);
                CHKERR();
             }
          }
       } else {
          if (lefCallbacks->NonDefaultCbk)
             lefData->lefrNonDefault.addMinCuts((yyvsp[-2].string), (int)(yyvsp[-1].dval));
       }
    }
#line 8474 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 451: /* $@66: %empty  */
#line 3832 "lef/lef.y"
                    { lefData->lefDumbMode = 10000000;}
#line 8480 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 452: /* nd_prop: K_PROPERTY $@66 nd_prop_list ';'  */
#line 3833 "lef/lef.y"
    { lefData->lefDumbMode = 0;
    }
#line 8487 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 455: /* nd_prop: T_STRING T_STRING  */
#line 3843 "lef/lef.y"
    {
      if (lefCallbacks->NonDefaultCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrNondefProp.propType((yyvsp[-1].string));
         lefData->lefrNonDefault.addProp((yyvsp[-1].string), (yyvsp[0].string), propTp);
      }
    }
#line 8499 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 456: /* nd_prop: T_STRING QSTRING  */
#line 3851 "lef/lef.y"
    {
      if (lefCallbacks->NonDefaultCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrNondefProp.propType((yyvsp[-1].string));
         lefData->lefrNonDefault.addProp((yyvsp[-1].string), (yyvsp[0].string), propTp);
      }
    }
#line 8511 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 457: /* nd_prop: T_STRING NUMBER  */
#line 3859 "lef/lef.y"
    {
      if (lefCallbacks->NonDefaultCbk) {
         char temp[32];
         char propTp;
         sprintf(temp, "%.11g", (yyvsp[0].dval));
         propTp = lefSettings->lefProps.lefrNondefProp.propType((yyvsp[-1].string));
         lefData->lefrNonDefault.addNumProp((yyvsp[-1].string), (yyvsp[0].dval), temp, propTp);
      }
    }
#line 8525 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 458: /* $@67: %empty  */
#line 3869 "lef/lef.y"
                  {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 8531 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 459: /* $@68: %empty  */
#line 3870 "lef/lef.y"
  {
    if (lefCallbacks->NonDefaultCbk) lefData->lefrNonDefault.addLayer((yyvsp[0].string));
    lefData->ndLayer++;
    //strcpy(lefData->layerName, $3);
    lefData->layerName = strdup((yyvsp[0].string));
    lefData->ndLayerWidth = 0;
    lefData->ndLayerSpace = 0;
  }
#line 8544 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 460: /* $@69: %empty  */
#line 3879 "lef/lef.y"
  { 
    lefData->ndLayerWidth = 1;
    if (lefCallbacks->NonDefaultCbk) lefData->lefrNonDefault.addWidth((yyvsp[-1].dval));
  }
#line 8553 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 461: /* $@70: %empty  */
#line 3883 "lef/lef.y"
                       {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 8559 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 462: /* nd_layer: K_LAYER $@67 T_STRING $@68 K_WIDTH int_number ';' $@69 nd_layer_stmts K_END $@70 T_STRING  */
#line 3884 "lef/lef.y"
  {
    if (strcmp(lefData->layerName, (yyvsp[0].string)) != 0) {
      if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
         if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings) {
            lefData->outMsg = (char*)lefMalloc(10000);
            sprintf (lefData->outMsg,
               "END LAYER name %s is different from the LAYER name %s.\nCorrect the LEF file before rerunning it through the LEF parser.", (yyvsp[-9].string), lefData->layerName);
            lefError(1624, lefData->outMsg);
            lefFree(lefData->outMsg);
            lefFree(lefData->layerName);
            CHKERR();
         } else
            lefFree(lefData->layerName);
      } else
         lefFree(lefData->layerName);
    } else
      lefFree(lefData->layerName);
    if (!lefData->ndLayerWidth) {
      if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
         if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings) {
            lefError(1625, "A WIDTH statement is required in the LAYER statement in NONDEFULTRULE.");
            CHKERR();
         }
      }
    }
    if (!lefData->ndLayerSpace && lefData->versionNum < 5.6) {   // 5.6, SPACING is optional
      if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
         if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings) {
            lefData->outMsg = (char*)lefMalloc(10000);
            sprintf (lefData->outMsg,
               "A SPACING statement is required in the LAYER statement in NONDEFAULTRULE for lef file with version 5.5 and earlier.\nYour lef file is defined with version %g. Update your lef to add a LAYER statement and try again.",
                lefData->versionNum);
            lefError(1626, lefData->outMsg);
            lefFree(lefData->outMsg);
            CHKERR();
         }
      }
    }
  }
#line 8603 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 465: /* nd_layer_stmt: K_SPACING int_number ';'  */
#line 3932 "lef/lef.y"
    {
      lefData->ndLayerSpace = 1;
      if (lefCallbacks->NonDefaultCbk) lefData->lefrNonDefault.addSpacing((yyvsp[-1].dval));
    }
#line 8612 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 466: /* nd_layer_stmt: K_WIREEXTENSION int_number ';'  */
#line 3937 "lef/lef.y"
    { if (lefCallbacks->NonDefaultCbk)
         lefData->lefrNonDefault.addWireExtension((yyvsp[-1].dval)); }
#line 8619 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 467: /* nd_layer_stmt: K_RESISTANCE K_RPERSQ int_number ';'  */
#line 3940 "lef/lef.y"
    {
      if (lefData->ignoreVersion) {
         if (lefCallbacks->NonDefaultCbk)
            lefData->lefrNonDefault.addResistance((yyvsp[-1].dval));
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
            if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "RESISTANCE RPERSQ statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1627, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else if (lefData->versionNum > 5.5) {  // obsolete in 5.6
         if (lefCallbacks->NonDefaultCbk) // write warning only if cbk is set 
            if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings)
              lefWarning(2029, "RESISTANCE RPERSQ statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
      } else if (lefCallbacks->NonDefaultCbk)
         lefData->lefrNonDefault.addResistance((yyvsp[-1].dval));
    }
#line 8646 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 468: /* nd_layer_stmt: K_CAPACITANCE K_CPERSQDIST int_number ';'  */
#line 3964 "lef/lef.y"
    {
      if (lefData->ignoreVersion) {
         if (lefCallbacks->NonDefaultCbk)
            lefData->lefrNonDefault.addCapacitance((yyvsp[-1].dval));
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
            if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "CAPACITANCE CPERSQDIST statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1628, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
            }
         }
      } else if (lefData->versionNum > 5.5) { // obsolete in 5.6
         if (lefCallbacks->NonDefaultCbk) // write warning only if cbk is set 
            if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings)
              lefWarning(2030, "CAPACITANCE CPERSQDIST statement is obsolete in version 5.6. and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
      } else if (lefCallbacks->NonDefaultCbk)
         lefData->lefrNonDefault.addCapacitance((yyvsp[-1].dval));
    }
#line 8673 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 469: /* nd_layer_stmt: K_EDGECAPACITANCE int_number ';'  */
#line 3987 "lef/lef.y"
    {
      if (lefData->ignoreVersion) {
         if (lefCallbacks->NonDefaultCbk)
            lefData->lefrNonDefault.addEdgeCap((yyvsp[-1].dval));
      } else if (lefData->versionNum < 5.4) {
         if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
            if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "EDGECAPACITANCE statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1629, lefData->outMsg);
               lefFree(lefData->outMsg);
              CHKERR();
            }
         }
      } else if (lefData->versionNum > 5.5) {  // obsolete in 5.6
         if (lefCallbacks->NonDefaultCbk) // write warning only if cbk is set 
            if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings)
              lefWarning(2031, "EDGECAPACITANCE statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
      } else if (lefCallbacks->NonDefaultCbk)
         lefData->lefrNonDefault.addEdgeCap((yyvsp[-1].dval));
    }
#line 8700 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 470: /* nd_layer_stmt: K_DIAGWIDTH int_number ';'  */
#line 4010 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {  // 5.6 syntax
         if (lefCallbacks->NonDefaultCbk) { // write error only if cbk is set 
            if (lefData->nonDefaultWarnings++ < lefSettings->NonDefaultWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                 "DIAGWIDTH statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
               lefError(1630, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR(); 
            }
         }
      } else {
         if (lefCallbacks->NonDefaultCbk)
            lefData->lefrNonDefault.addDiagWidth((yyvsp[-1].dval));
      }
    }
#line 8722 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 471: /* site: start_site site_options end_site  */
#line 4029 "lef/lef.y"
    { 
      if (lefCallbacks->SiteCbk)
        CALLBACK(lefCallbacks->SiteCbk, lefrSiteCbkType, &lefData->lefrSite);
    }
#line 8731 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 472: /* $@71: %empty  */
#line 4034 "lef/lef.y"
                   {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 8737 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 473: /* start_site: K_SITE $@71 T_STRING  */
#line 4035 "lef/lef.y"
    { 
      if (lefCallbacks->SiteCbk) lefData->lefrSite.setName((yyvsp[0].string));
      //strcpy(lefData->siteName, $3);
      lefData->siteName = strdup((yyvsp[0].string));
      lefData->hasSiteClass = 0;
      lefData->hasSiteSize = 0;
      lefData->hasSite = 1;
    }
#line 8750 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 474: /* $@72: %empty  */
#line 4044 "lef/lef.y"
                {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 8756 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 475: /* end_site: K_END $@72 T_STRING  */
#line 4045 "lef/lef.y"
    {
      if (strcmp(lefData->siteName, (yyvsp[0].string)) != 0) {
        if (lefCallbacks->SiteCbk) { // write error only if cbk is set 
           if (lefData->siteWarnings++ < lefSettings->SiteWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "END SITE name %s is different from the SITE name %s.\nCorrect the LEF file before rerunning it through the LEF parser.", (yyvsp[0].string), lefData->siteName);
              lefError(1631, lefData->outMsg);
              lefFree(lefData->outMsg);
              lefFree(lefData->siteName);
              CHKERR();
           } else
              lefFree(lefData->siteName);
        } else
           lefFree(lefData->siteName);
      } else {
        lefFree(lefData->siteName);
        if (lefCallbacks->SiteCbk) { // write error only if cbk is set 
          if (lefData->hasSiteClass == 0) {
             lefError(1632, "A CLASS statement is required in the SITE statement.");
             CHKERR();
          }
          if (lefData->hasSiteSize == 0) {
             lefError(1633, "A SIZE  statement is required in the SITE statement.");
             CHKERR();
          }
        }
      }
    }
#line 8790 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 478: /* site_option: K_SIZE int_number K_BY int_number ';'  */
#line 4082 "lef/lef.y"
    {

      if (lefCallbacks->SiteCbk) lefData->lefrSite.setSize((yyvsp[-3].dval),(yyvsp[-1].dval));
      lefData->hasSiteSize = 1;
    }
#line 8800 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 479: /* site_option: site_symmetry_statement  */
#line 4088 "lef/lef.y"
    { }
#line 8806 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 480: /* site_option: site_class  */
#line 4090 "lef/lef.y"
    { 
      if (lefCallbacks->SiteCbk) lefData->lefrSite.setClass((yyvsp[0].string));
      lefData->hasSiteClass = 1;
    }
#line 8815 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 481: /* site_option: site_rowpattern_statement  */
#line 4095 "lef/lef.y"
    { }
#line 8821 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 482: /* site_class: K_CLASS K_PAD ';'  */
#line 4098 "lef/lef.y"
                    {(yyval.string) = (char*)"PAD"; }
#line 8827 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 483: /* site_class: K_CLASS K_CORE ';'  */
#line 4099 "lef/lef.y"
                        {(yyval.string) = (char*)"CORE"; }
#line 8833 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 484: /* site_class: K_CLASS K_VIRTUAL ';'  */
#line 4100 "lef/lef.y"
                           {(yyval.string) = (char*)"VIRTUAL"; }
#line 8839 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 485: /* site_symmetry_statement: K_SYMMETRY site_symmetries ';'  */
#line 4103 "lef/lef.y"
    { }
#line 8845 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 488: /* site_symmetry: K_X  */
#line 4112 "lef/lef.y"
    { if (lefCallbacks->SiteCbk) lefData->lefrSite.setXSymmetry(); }
#line 8851 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 489: /* site_symmetry: K_Y  */
#line 4114 "lef/lef.y"
    { if (lefCallbacks->SiteCbk) lefData->lefrSite.setYSymmetry(); }
#line 8857 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 490: /* site_symmetry: K_R90  */
#line 4116 "lef/lef.y"
    { if (lefCallbacks->SiteCbk) lefData->lefrSite.set90Symmetry(); }
#line 8863 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 491: /* $@73: %empty  */
#line 4118 "lef/lef.y"
                                        {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 8869 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 492: /* site_rowpattern_statement: K_ROWPATTERN $@73 site_rowpatterns ';'  */
#line 4120 "lef/lef.y"
    { }
#line 8875 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 495: /* $@74: %empty  */
#line 4127 "lef/lef.y"
                                      {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 8881 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 496: /* site_rowpattern: T_STRING orientation $@74  */
#line 4128 "lef/lef.y"
    { if (lefCallbacks->SiteCbk) lefData->lefrSite.addRowPattern((yyvsp[-2].string), (yyvsp[-1].integer)); }
#line 8887 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 497: /* pt: int_number int_number  */
#line 4132 "lef/lef.y"
    { (yyval.pt).x = (yyvsp[-1].dval); (yyval.pt).y = (yyvsp[0].dval); }
#line 8893 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 498: /* pt: '(' int_number int_number ')'  */
#line 4134 "lef/lef.y"
    { (yyval.pt).x = (yyvsp[-2].dval); (yyval.pt).y = (yyvsp[-1].dval); }
#line 8899 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 499: /* $@75: %empty  */
#line 4137 "lef/lef.y"
    { 
      if (lefCallbacks->MacroCbk)
        CALLBACK(lefCallbacks->MacroCbk, lefrMacroCbkType, &lefData->lefrMacro);
      lefData->lefrDoSite = 0;
    }
#line 8909 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 501: /* $@76: %empty  */
#line 4144 "lef/lef.y"
                     {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 8915 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 502: /* start_macro: K_MACRO $@76 T_STRING  */
#line 4145 "lef/lef.y"
    {
      lefData->siteDef = 0;
      lefData->symDef = 0;
      lefData->sizeDef = 0; 
      lefData->pinDef = 0; 
      lefData->obsDef = 0; 
      lefData->origDef = 0;
      lefData->lefrMacro.clear();      
      if (lefCallbacks->MacroBeginCbk || lefCallbacks->MacroCbk) {
        // some reader may not have MacroBeginCB, but has MacroCB set
        lefData->lefrMacro.setName((yyvsp[0].string));
        CALLBACK(lefCallbacks->MacroBeginCbk, lefrMacroBeginCbkType, (yyvsp[0].string));
      }
      //strcpy(lefData->macroName, $3);
      lefData->macroName = strdup((yyvsp[0].string));
    }
#line 8936 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 503: /* $@77: %empty  */
#line 4162 "lef/lef.y"
                 {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 8942 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 504: /* end_macro: K_END $@77 T_STRING  */
#line 4163 "lef/lef.y"
    {
      if (strcmp(lefData->macroName, (yyvsp[0].string)) != 0) {
        if (lefCallbacks->MacroEndCbk) { // write error only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "END MACRO name %s is different from the MACRO name %s.\nCorrect the LEF file before rerunning it through the LEF parser.", (yyvsp[0].string), lefData->macroName);
              lefError(1634, lefData->outMsg);
              lefFree(lefData->outMsg);
              lefFree(lefData->macroName);
              CHKERR();
           } else
              lefFree(lefData->macroName);
        } else
           lefFree(lefData->macroName);
      } else
        lefFree(lefData->macroName);
      if (lefCallbacks->MacroEndCbk)
        CALLBACK(lefCallbacks->MacroEndCbk, lefrMacroEndCbkType, (yyvsp[0].string));
    }
#line 8967 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 512: /* macro_option: macro_fixedMask  */
#line 4196 "lef/lef.y"
      { }
#line 8973 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 513: /* macro_option: macro_origin  */
#line 4198 "lef/lef.y"
      { }
#line 8979 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 514: /* macro_option: macro_power  */
#line 4200 "lef/lef.y"
      { }
#line 8985 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 515: /* macro_option: macro_foreign  */
#line 4202 "lef/lef.y"
      { }
#line 8991 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 518: /* macro_option: macro_size  */
#line 4206 "lef/lef.y"
      { }
#line 8997 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 519: /* macro_option: macro_site  */
#line 4208 "lef/lef.y"
      { }
#line 9003 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 520: /* macro_option: macro_pin  */
#line 4210 "lef/lef.y"
      { }
#line 9009 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 521: /* macro_option: K_FUNCTION K_BUFFER ';'  */
#line 4212 "lef/lef.y"
      { if (lefCallbacks->MacroCbk) lefData->lefrMacro.setBuffer(); }
#line 9015 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 522: /* macro_option: K_FUNCTION K_INVERTER ';'  */
#line 4214 "lef/lef.y"
      { if (lefCallbacks->MacroCbk) lefData->lefrMacro.setInverter(); }
#line 9021 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 523: /* macro_option: macro_obs  */
#line 4216 "lef/lef.y"
      { }
#line 9027 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 524: /* macro_option: macro_density  */
#line 4218 "lef/lef.y"
      { }
#line 9033 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 525: /* macro_option: macro_clocktype  */
#line 4220 "lef/lef.y"
      { }
#line 9039 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 526: /* macro_option: timing  */
#line 4222 "lef/lef.y"
      { }
#line 9045 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 527: /* $@78: %empty  */
#line 4223 "lef/lef.y"
               {lefData->lefDumbMode = 1000000; }
#line 9051 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 528: /* macro_option: K_PROPERTY $@78 macro_prop_list ';'  */
#line 4224 "lef/lef.y"
      { lefData->lefDumbMode = 0;
      }
#line 9058 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 531: /* macro_symmetry_statement: K_SYMMETRY macro_symmetries ';'  */
#line 4233 "lef/lef.y"
    {
      if (lefData->siteDef) { // SITE is defined before SYMMETRY 
          // pcr 283846 suppress warning 
          if (lefCallbacks->MacroCbk) // write warning only if cbk is set 
             if (lefData->macroWarnings++ < lefSettings->MacroWarnings)
               lefWarning(2032, "A SITE statement is defined before SYMMETRY statement.\nTo avoid this warning in the future, define SITE after SYMMETRY");
      }
      lefData->symDef = 1;
    }
#line 9072 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 534: /* macro_symmetry: K_X  */
#line 4250 "lef/lef.y"
    { if (lefCallbacks->MacroCbk) lefData->lefrMacro.setXSymmetry(); }
#line 9078 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 535: /* macro_symmetry: K_Y  */
#line 4252 "lef/lef.y"
    { if (lefCallbacks->MacroCbk) lefData->lefrMacro.setYSymmetry(); }
#line 9084 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 536: /* macro_symmetry: K_R90  */
#line 4254 "lef/lef.y"
    { if (lefCallbacks->MacroCbk) lefData->lefrMacro.set90Symmetry(); }
#line 9090 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 537: /* macro_name_value_pair: T_STRING NUMBER  */
#line 4258 "lef/lef.y"
    {
      char temp[32];
      sprintf(temp, "%.11g", (yyvsp[0].dval));
      if (lefCallbacks->MacroCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrMacroProp.propType((yyvsp[-1].string));
         lefData->lefrMacro.setNumProperty((yyvsp[-1].string), (yyvsp[0].dval), temp,  propTp);
      }
    }
#line 9104 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 538: /* macro_name_value_pair: T_STRING QSTRING  */
#line 4268 "lef/lef.y"
    {
      if (lefCallbacks->MacroCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrMacroProp.propType((yyvsp[-1].string));
         lefData->lefrMacro.setProperty((yyvsp[-1].string), (yyvsp[0].string), propTp);
      }
    }
#line 9116 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 539: /* macro_name_value_pair: T_STRING T_STRING  */
#line 4276 "lef/lef.y"
    {
      if (lefCallbacks->MacroCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrMacroProp.propType((yyvsp[-1].string));
         lefData->lefrMacro.setProperty((yyvsp[-1].string), (yyvsp[0].string), propTp);
      }
    }
#line 9128 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 540: /* macro_class: K_CLASS class_type ';'  */
#line 4285 "lef/lef.y"
    {
       if (lefCallbacks->MacroCbk) lefData->lefrMacro.setClass((yyvsp[-1].string));
       if (lefCallbacks->MacroClassTypeCbk)
          CALLBACK(lefCallbacks->MacroClassTypeCbk, lefrMacroClassTypeCbkType, (yyvsp[-1].string));
    }
#line 9138 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 541: /* class_type: K_COVER  */
#line 4292 "lef/lef.y"
          {(yyval.string) = (char*)"COVER"; }
#line 9144 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 542: /* class_type: K_COVER K_BUMP  */
#line 4294 "lef/lef.y"
    { (yyval.string) = (char*)"COVER BUMP";
      if (lefData->versionNum < 5.5) {
        if (lefCallbacks->MacroCbk) { // write error only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              if (lefSettings->RelaxMode)
                 lefWarning(2033, "The statement COVER BUMP is a LEF verion 5.5 syntax.\nYour LEF file is version 5.4 or earlier which is incorrect but will be allowed\nbecause this application does not enforce strict version checking.\nOther tools that enforce strict checking will have a syntax error when reading this file.\nYou can change the VERSION statement in this LEF file to 5.5 or higher to stop this warning.");
              else {
                 lefData->outMsg = (char*)lefMalloc(10000);
                 sprintf (lefData->outMsg,
                    "COVER BUMP statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
                 lefError(1635, lefData->outMsg);
                 lefFree(lefData->outMsg);
                 CHKERR();
              }
           }
        }
      }
    }
#line 9167 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 543: /* class_type: K_RING  */
#line 4312 "lef/lef.y"
              {(yyval.string) = (char*)"RING"; }
#line 9173 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 544: /* class_type: K_BLOCK  */
#line 4313 "lef/lef.y"
              {(yyval.string) = (char*)"BLOCK"; }
#line 9179 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 545: /* class_type: K_BLOCK K_BLACKBOX  */
#line 4315 "lef/lef.y"
    { (yyval.string) = (char*)"BLOCK BLACKBOX";
      if (lefData->versionNum < 5.5) {
        if (lefCallbacks->MacroCbk) { // write error only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
             if (lefSettings->RelaxMode)
                lefWarning(2034, "The statement BLOCK BLACKBOX is a LEF verion 5.5 syntax.\nYour LEF file is version 5.4 or earlier which is incorrect but will be allowed\nbecause this application does not enforce strict version checking.\nOther tools that enforce strict checking will have a syntax error when reading this file.\nYou can change the VERSION statement in this LEF file to 5.5 or higher to stop this warning.");
              else {
                 lefData->outMsg = (char*)lefMalloc(10000);
                 sprintf (lefData->outMsg,
                    "BLOCK BLACKBOX statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
                 lefError(1636, lefData->outMsg);
                 lefFree(lefData->outMsg);
                 CHKERR();
              }
           }
        }
      }
    }
#line 9202 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 546: /* class_type: K_BLOCK K_SOFT  */
#line 4334 "lef/lef.y"
    {
      if (lefData->ignoreVersion) {
        (yyval.string) = (char*)"BLOCK SOFT";
      } else if (lefData->versionNum < 5.6) {
        if (lefCallbacks->MacroCbk) { // write error only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "BLOCK SOFT statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1637, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
      else
        (yyval.string) = (char*)"BLOCK SOFT";
    }
#line 9225 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 547: /* class_type: K_NONE  */
#line 4352 "lef/lef.y"
              {(yyval.string) = (char*)"NONE"; }
#line 9231 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 548: /* class_type: K_BUMP  */
#line 4354 "lef/lef.y"
      {
        if (lefData->versionNum < 5.7) {
          lefData->outMsg = (char*)lefMalloc(10000);
          sprintf(lefData->outMsg,
            "BUMP is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
          lefError(1698, lefData->outMsg);
          lefFree(lefData->outMsg);
          CHKERR();
        }
       
        (yyval.string) = (char*)"BUMP";
     }
#line 9248 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 549: /* class_type: K_PAD  */
#line 4366 "lef/lef.y"
              {(yyval.string) = (char*)"PAD"; }
#line 9254 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 550: /* class_type: K_VIRTUAL  */
#line 4367 "lef/lef.y"
              {(yyval.string) = (char*)"VIRTUAL"; }
#line 9260 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 551: /* class_type: K_PAD pad_type  */
#line 4369 "lef/lef.y"
      {  sprintf(lefData->temp_name, "PAD %s", (yyvsp[0].string));
        (yyval.string) = lefData->temp_name; 
        if (lefData->versionNum < 5.5) {
           if (strcmp("AREAIO", (yyvsp[0].string)) != 0) {
             sprintf(lefData->temp_name, "PAD %s", (yyvsp[0].string));
             (yyval.string) = lefData->temp_name; 
           } else if (lefCallbacks->MacroCbk) { 
             if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
               if (lefSettings->RelaxMode)
                  lefWarning(2035, "The statement PAD AREAIO is a LEF verion 5.5 syntax.\nYour LEF file is version 5.4 or earlier which is incorrect but will be allowed\nbecause this application does not enforce strict version checking.\nOther tools that enforce strict checking will have a syntax error when reading this file.\nYou can change the VERSION statement in this LEF file to 5.5 or higher to stop this warning.");
               else {
                  lefData->outMsg = (char*)lefMalloc(10000);
                  sprintf (lefData->outMsg,
                     "PAD AREAIO statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
                  lefError(1638, lefData->outMsg);
                  lefFree(lefData->outMsg);
                  CHKERR();
               }
            }
          }
        }
      }
#line 9287 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 552: /* class_type: K_CORE  */
#line 4391 "lef/lef.y"
              {(yyval.string) = (char*)"CORE"; }
#line 9293 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 553: /* class_type: K_CORNER  */
#line 4393 "lef/lef.y"
      {(yyval.string) = (char*)"CORNER";
      // This token is NOT in the spec but has shown up in 
      // some lef files.  This exception came from LEFOUT
      // in 'frameworks'
      }
#line 9303 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 554: /* class_type: K_CORE core_type  */
#line 4399 "lef/lef.y"
      {sprintf(lefData->temp_name, "CORE %s", (yyvsp[0].string));
      (yyval.string) = lefData->temp_name;}
#line 9310 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 555: /* class_type: K_ENDCAP endcap_type  */
#line 4402 "lef/lef.y"
      {sprintf(lefData->temp_name, "ENDCAP %s", (yyvsp[0].string));
      (yyval.string) = lefData->temp_name;}
#line 9317 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 556: /* pad_type: K_INPUT  */
#line 4406 "lef/lef.y"
                  {(yyval.string) = (char*)"INPUT";}
#line 9323 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 557: /* pad_type: K_OUTPUT  */
#line 4407 "lef/lef.y"
                    {(yyval.string) = (char*)"OUTPUT";}
#line 9329 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 558: /* pad_type: K_INOUT  */
#line 4408 "lef/lef.y"
                    {(yyval.string) = (char*)"INOUT";}
#line 9335 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 559: /* pad_type: K_POWER  */
#line 4409 "lef/lef.y"
                    {(yyval.string) = (char*)"POWER";}
#line 9341 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 560: /* pad_type: K_SPACER  */
#line 4410 "lef/lef.y"
                    {(yyval.string) = (char*)"SPACER";}
#line 9347 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 561: /* pad_type: K_AREAIO  */
#line 4411 "lef/lef.y"
                {(yyval.string) = (char*)"AREAIO";}
#line 9353 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 562: /* core_type: K_FEEDTHRU  */
#line 4414 "lef/lef.y"
                    {(yyval.string) = (char*)"FEEDTHRU";}
#line 9359 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 563: /* core_type: K_TIEHIGH  */
#line 4415 "lef/lef.y"
                    {(yyval.string) = (char*)"TIEHIGH";}
#line 9365 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 564: /* core_type: K_TIELOW  */
#line 4416 "lef/lef.y"
                    {(yyval.string) = (char*)"TIELOW";}
#line 9371 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 565: /* core_type: K_SPACER  */
#line 4418 "lef/lef.y"
    { 
      (yyval.string) = (char*)"SPACER";

      if (!lefData->ignoreVersion && lefData->versionNum < 5.4) {
        if (lefCallbacks->MacroCbk) { // write error only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "SPACER statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1639, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
    }
#line 9392 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 566: /* core_type: K_ANTENNACELL  */
#line 4435 "lef/lef.y"
    { 
      (yyval.string) = (char*)"ANTENNACELL";

      if (!lefData->ignoreVersion && lefData->versionNum < 5.4) {
        if (lefCallbacks->MacroCbk) { // write error only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNACELL statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1640, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
    }
#line 9413 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 567: /* core_type: K_WELLTAP  */
#line 4452 "lef/lef.y"
    { 
      (yyval.string) = (char*)"WELLTAP";

      if (!lefData->ignoreVersion && lefData->versionNum < 5.6) {
        if (lefCallbacks->MacroCbk) { // write error only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "WELLTAP statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1641, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
    }
#line 9434 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 568: /* endcap_type: K_PRE  */
#line 4470 "lef/lef.y"
                {(yyval.string) = (char*)"PRE";}
#line 9440 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 569: /* endcap_type: K_POST  */
#line 4471 "lef/lef.y"
                   {(yyval.string) = (char*)"POST";}
#line 9446 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 570: /* endcap_type: K_TOPLEFT  */
#line 4472 "lef/lef.y"
                      {(yyval.string) = (char*)"TOPLEFT";}
#line 9452 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 571: /* endcap_type: K_TOPRIGHT  */
#line 4473 "lef/lef.y"
                       {(yyval.string) = (char*)"TOPRIGHT";}
#line 9458 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 572: /* endcap_type: K_BOTTOMLEFT  */
#line 4474 "lef/lef.y"
                         {(yyval.string) = (char*)"BOTTOMLEFT";}
#line 9464 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 573: /* endcap_type: K_BOTTOMRIGHT  */
#line 4475 "lef/lef.y"
                         {(yyval.string) = (char*)"BOTTOMRIGHT";}
#line 9470 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 574: /* macro_generator: K_GENERATOR T_STRING ';'  */
#line 4478 "lef/lef.y"
    { if (lefCallbacks->MacroCbk) lefData->lefrMacro.setGenerator((yyvsp[-1].string)); }
#line 9476 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 575: /* macro_generate: K_GENERATE T_STRING T_STRING ';'  */
#line 4481 "lef/lef.y"
    { if (lefCallbacks->MacroCbk) lefData->lefrMacro.setGenerate((yyvsp[-2].string), (yyvsp[-1].string)); }
#line 9482 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 576: /* macro_source: K_SOURCE K_USER ';'  */
#line 4485 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->MacroCbk) lefData->lefrMacro.setSource("USER");
      } else
        if (lefCallbacks->MacroCbk) // write warning only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings)
             lefWarning(2036, "SOURCE statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 9495 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 577: /* macro_source: K_SOURCE K_GENERATE ';'  */
#line 4494 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->MacroCbk) lefData->lefrMacro.setSource("GENERATE");
      } else
        if (lefCallbacks->MacroCbk) // write warning only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings)
             lefWarning(2037, "SOURCE statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 9508 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 578: /* macro_source: K_SOURCE K_BLOCK ';'  */
#line 4503 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->MacroCbk) lefData->lefrMacro.setSource("BLOCK");
      } else
        if (lefCallbacks->MacroCbk) // write warning only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings)
             lefWarning(2037, "SOURCE statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 9521 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 579: /* macro_power: K_POWER int_number ';'  */
#line 4513 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->MacroCbk) lefData->lefrMacro.setPower((yyvsp[-1].dval));
      } else
        if (lefCallbacks->MacroCbk) // write warning only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings)
             lefWarning(2038, "MACRO POWER statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9534 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 580: /* macro_origin: K_ORIGIN pt ';'  */
#line 4523 "lef/lef.y"
    { 
       if (lefData->origDef) { // Has multiple ORIGIN defined in a macro, stop parsing
          if (lefCallbacks->MacroCbk) { // write error only if cbk is set 
             if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
                lefError(1642, "ORIGIN statement has defined more than once in a MACRO statement.\nOnly one ORIGIN statement can be defined in a Macro.\nParser will stop processing.");
               CHKERR();
             }
          }
       }
       lefData->origDef = 1;
       if (lefData->siteDef) { // SITE is defined before ORIGIN 
          // pcr 283846 suppress warning 
          if (lefCallbacks->MacroCbk) // write warning only if cbk is set 
             if (lefData->macroWarnings++ < lefSettings->MacroWarnings)
               lefWarning(2039, "A SITE statement is defined before ORIGIN statement.\nTo avoid this warning in the future, define SITE after ORIGIN");
       }
       if (lefData->pinDef) { // PIN is defined before ORIGIN 
          // pcr 283846 suppress warning 
          if (lefCallbacks->MacroCbk) // write warning only if cbk is set 
             if (lefData->macroWarnings++ < lefSettings->MacroWarnings)
               lefWarning(2040, "A PIN statement is defined before ORIGIN statement.\nTo avoid this warning in the future, define PIN after ORIGIN");
       }
       if (lefData->obsDef) { // OBS is defined before ORIGIN 
          // pcr 283846 suppress warning 
          if (lefCallbacks->MacroCbk) // write warning only if cbk is set 
             if (lefData->macroWarnings++ < lefSettings->MacroWarnings)
               lefWarning(2041, "A OBS statement is defined before ORIGIN statement.\nTo avoid this warning in the future, define OBS after ORIGIN");
       }
      
       // Workaround for pcr 640902 
       if (lefCallbacks->MacroCbk) lefData->lefrMacro.setOrigin((yyvsp[-1].pt).x, (yyvsp[-1].pt).y);
       if (lefCallbacks->MacroOriginCbk) {
          lefData->macroNum.x = (yyvsp[-1].pt).x; 
          lefData->macroNum.y = (yyvsp[-1].pt).y; 
          CALLBACK(lefCallbacks->MacroOriginCbk, lefrMacroOriginCbkType, lefData->macroNum);
       }
    }
#line 9576 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 581: /* macro_foreign: start_foreign ';'  */
#line 4563 "lef/lef.y"
    { if (lefCallbacks->MacroCbk)
      lefData->lefrMacro.addForeign((yyvsp[-1].string), 0, 0.0, 0.0, -1);
    }
#line 9584 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 582: /* macro_foreign: start_foreign pt ';'  */
#line 4567 "lef/lef.y"
    { if (lefCallbacks->MacroCbk)
      lefData->lefrMacro.addForeign((yyvsp[-2].string), 1, (yyvsp[-1].pt).x, (yyvsp[-1].pt).y, -1);
    }
#line 9592 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 583: /* macro_foreign: start_foreign pt orientation ';'  */
#line 4571 "lef/lef.y"
    { if (lefCallbacks->MacroCbk)
      lefData->lefrMacro.addForeign((yyvsp[-3].string), 1, (yyvsp[-2].pt).x, (yyvsp[-2].pt).y, (yyvsp[-1].integer));
    }
#line 9600 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 584: /* macro_foreign: start_foreign orientation ';'  */
#line 4575 "lef/lef.y"
    { if (lefCallbacks->MacroCbk)
      lefData->lefrMacro.addForeign((yyvsp[-2].string), 0, 0.0, 0.0, (yyvsp[-1].integer));
    }
#line 9608 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 585: /* macro_fixedMask: K_FIXEDMASK ';'  */
#line 4581 "lef/lef.y"
   {   
       if (lefCallbacks->MacroCbk && lefData->versionNum >= 5.8) {
          lefData->lefrMacro.setFixedMask(1);
       }
       if (lefCallbacks->MacroFixedMaskCbk) {
          CALLBACK(lefCallbacks->MacroFixedMaskCbk, lefrMacroFixedMaskCbkType, 1);
       }        
    }
#line 9621 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 586: /* $@79: %empty  */
#line 4590 "lef/lef.y"
                 { lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 9627 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 587: /* macro_eeq: K_EEQ $@79 T_STRING ';'  */
#line 4591 "lef/lef.y"
    { if (lefCallbacks->MacroCbk) lefData->lefrMacro.setEEQ((yyvsp[-1].string)); }
#line 9633 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 588: /* $@80: %empty  */
#line 4593 "lef/lef.y"
                 { lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 9639 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 589: /* macro_leq: K_LEQ $@80 T_STRING ';'  */
#line 4594 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->MacroCbk) lefData->lefrMacro.setLEQ((yyvsp[-1].string));
      } else
        if (lefCallbacks->MacroCbk) // write warning only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings)
             lefWarning(2042, "LEQ statement in MACRO is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 9652 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 590: /* macro_site: macro_site_word T_STRING ';'  */
#line 4605 "lef/lef.y"
    {
      if (lefCallbacks->MacroCbk) {
        lefData->lefrMacro.setSiteName((yyvsp[-1].string));
      }
    }
#line 9662 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 591: /* macro_site: macro_site_word sitePattern ';'  */
#line 4611 "lef/lef.y"
    {
      if (lefCallbacks->MacroCbk) {
        // also set site name in the variable siteName_ in lefiMacro 
        // this, if user wants to use method lefData->siteName will get the name also 
        lefData->lefrMacro.setSitePattern(lefData->lefrSitePatternPtr);
        lefData->lefrSitePatternPtr = 0;
      }
    }
#line 9675 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 592: /* macro_site_word: K_SITE  */
#line 4621 "lef/lef.y"
    { lefData->lefDumbMode = 1; lefData->lefNoNum = 1; lefData->siteDef = 1;
        if (lefCallbacks->MacroCbk) lefData->lefrDoSite = 1; }
#line 9682 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 593: /* site_word: K_SITE  */
#line 4625 "lef/lef.y"
    { lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 9688 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 594: /* macro_size: K_SIZE int_number K_BY int_number ';'  */
#line 4628 "lef/lef.y"
    { 
      if (lefData->siteDef) { // SITE is defined before SIZE 
      }
      lefData->sizeDef = 1;
      if (lefCallbacks->MacroCbk) lefData->lefrMacro.setSize((yyvsp[-3].dval), (yyvsp[-1].dval));
      if (lefCallbacks->MacroSizeCbk) {
         lefData->macroNum.x = (yyvsp[-3].dval); 
         lefData->macroNum.y = (yyvsp[-1].dval); 
         CALLBACK(lefCallbacks->MacroSizeCbk, lefrMacroSizeCbkType, lefData->macroNum);
      }
    }
#line 9704 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 595: /* macro_pin: start_macro_pin macro_pin_options end_macro_pin  */
#line 4644 "lef/lef.y"
    { 
      if (lefCallbacks->PinCbk)
        CALLBACK(lefCallbacks->PinCbk, lefrPinCbkType, &lefData->lefrPin);
      lefData->lefrPin.clear();
    }
#line 9714 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 596: /* $@81: %empty  */
#line 4650 "lef/lef.y"
                       {lefData->lefDumbMode = 1; lefData->lefNoNum = 1; lefData->pinDef = 1;}
#line 9720 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 597: /* start_macro_pin: K_PIN $@81 T_STRING  */
#line 4651 "lef/lef.y"
    { if (lefCallbacks->PinCbk) lefData->lefrPin.setName((yyvsp[0].string));
      //strcpy(lefData->pinName, $3);
      lefData->pinName = strdup((yyvsp[0].string));
    }
#line 9729 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 598: /* $@82: %empty  */
#line 4656 "lef/lef.y"
                     {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 9735 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 599: /* end_macro_pin: K_END $@82 T_STRING  */
#line 4657 "lef/lef.y"
    {
      if (strcmp(lefData->pinName, (yyvsp[0].string)) != 0) {
        if (lefCallbacks->MacroCbk) { // write error only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "END PIN name %s is different from the PIN name %s.\nCorrect the LEF file before rerunning it through the LEF parser.", (yyvsp[0].string), lefData->pinName);
              lefError(1643, lefData->outMsg);
              lefFree(lefData->outMsg);
              lefFree(lefData->pinName);
              CHKERR();
           } else
              lefFree(lefData->pinName);
        } else
           lefFree(lefData->pinName);
      } else
        lefFree(lefData->pinName);
    }
#line 9758 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 600: /* macro_pin_options: %empty  */
#line 4678 "lef/lef.y"
    { }
#line 9764 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 601: /* macro_pin_options: macro_pin_options macro_pin_option  */
#line 4680 "lef/lef.y"
    { }
#line 9770 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 602: /* macro_pin_option: start_foreign ';'  */
#line 4684 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.addForeign((yyvsp[-1].string), 0, 0.0, 0.0, -1);
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2043, "FOREIGN statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 9783 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 603: /* macro_pin_option: start_foreign pt ';'  */
#line 4693 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.addForeign((yyvsp[-2].string), 1, (yyvsp[-1].pt).x, (yyvsp[-1].pt).y, -1);
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2043, "FOREIGN statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 9796 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 604: /* macro_pin_option: start_foreign pt orientation ';'  */
#line 4702 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.addForeign((yyvsp[-3].string), 1, (yyvsp[-2].pt).x, (yyvsp[-2].pt).y, (yyvsp[-1].integer));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2043, "FOREIGN statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 9809 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 605: /* macro_pin_option: start_foreign K_STRUCTURE ';'  */
#line 4711 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.addForeign((yyvsp[-2].string), 0, 0.0, 0.0, -1);
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2043, "FOREIGN statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 9822 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 606: /* macro_pin_option: start_foreign K_STRUCTURE pt ';'  */
#line 4720 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.addForeign((yyvsp[-3].string), 1, (yyvsp[-1].pt).x, (yyvsp[-1].pt).y, -1);
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2043, "FOREIGN statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 9835 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 607: /* macro_pin_option: start_foreign K_STRUCTURE pt orientation ';'  */
#line 4729 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.addForeign((yyvsp[-4].string), 1, (yyvsp[-2].pt).x, (yyvsp[-2].pt).y, (yyvsp[-1].integer));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2043, "FOREIGN statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 9848 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 608: /* $@83: %empty  */
#line 4737 "lef/lef.y"
          { lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 9854 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 609: /* macro_pin_option: K_LEQ $@83 T_STRING ';'  */
#line 4738 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setLEQ((yyvsp[-1].string));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2044, "LEQ statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
   }
#line 9867 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 610: /* macro_pin_option: K_POWER int_number ';'  */
#line 4747 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setPower((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2045, "MACRO POWER statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9880 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 611: /* macro_pin_option: electrical_direction  */
#line 4756 "lef/lef.y"
    { if (lefCallbacks->PinCbk) lefData->lefrPin.setDirection((yyvsp[0].string)); }
#line 9886 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 612: /* macro_pin_option: K_USE macro_pin_use ';'  */
#line 4758 "lef/lef.y"
    { if (lefCallbacks->PinCbk) lefData->lefrPin.setUse((yyvsp[-1].string)); }
#line 9892 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 613: /* macro_pin_option: K_SCANUSE macro_scan_use ';'  */
#line 4760 "lef/lef.y"
    { }
#line 9898 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 614: /* macro_pin_option: K_LEAKAGE int_number ';'  */
#line 4762 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setLeakage((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2046, "MACRO LEAKAGE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, r emove this statement from the LEF file with version 5.4 or later.");
    }
#line 9911 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 615: /* macro_pin_option: K_RISETHRESH int_number ';'  */
#line 4771 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setRiseThresh((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2047, "MACRO RISETHRESH statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9924 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 616: /* macro_pin_option: K_FALLTHRESH int_number ';'  */
#line 4780 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setFallThresh((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2048, "MACRO FALLTHRESH statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9937 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 617: /* macro_pin_option: K_RISESATCUR int_number ';'  */
#line 4789 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setRiseSatcur((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2049, "MACRO RISESATCUR statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9950 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 618: /* macro_pin_option: K_FALLSATCUR int_number ';'  */
#line 4798 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setFallSatcur((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2050, "MACRO FALLSATCUR statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9963 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 619: /* macro_pin_option: K_VLO int_number ';'  */
#line 4807 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setVLO((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2051, "MACRO VLO statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9976 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 620: /* macro_pin_option: K_VHI int_number ';'  */
#line 4816 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setVHI((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2052, "MACRO VHI statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9989 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 621: /* macro_pin_option: K_TIEOFFR int_number ';'  */
#line 4825 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setTieoffr((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2053, "MACRO TIEOFFR statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 10002 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 622: /* macro_pin_option: K_SHAPE pin_shape ';'  */
#line 4834 "lef/lef.y"
    { if (lefCallbacks->PinCbk) lefData->lefrPin.setShape((yyvsp[-1].string)); }
#line 10008 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 623: /* $@84: %empty  */
#line 4835 "lef/lef.y"
               {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 10014 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 624: /* macro_pin_option: K_MUSTJOIN $@84 T_STRING ';'  */
#line 4836 "lef/lef.y"
    { if (lefCallbacks->PinCbk) lefData->lefrPin.setMustjoin((yyvsp[-1].string)); }
#line 10020 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 625: /* $@85: %empty  */
#line 4837 "lef/lef.y"
                        {lefData->lefDumbMode = 1;}
#line 10026 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 626: /* macro_pin_option: K_OUTPUTNOISEMARGIN $@85 int_number int_number ';'  */
#line 4838 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setOutMargin((yyvsp[-2].dval), (yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2054, "MACRO OUTPUTNOISEMARGIN statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 10039 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 627: /* $@86: %empty  */
#line 4846 "lef/lef.y"
                       {lefData->lefDumbMode = 1;}
#line 10045 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 628: /* macro_pin_option: K_OUTPUTRESISTANCE $@86 int_number int_number ';'  */
#line 4847 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setOutResistance((yyvsp[-2].dval), (yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2055, "MACRO OUTPUTRESISTANCE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 10058 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 629: /* $@87: %empty  */
#line 4855 "lef/lef.y"
                       {lefData->lefDumbMode = 1;}
#line 10064 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 630: /* macro_pin_option: K_INPUTNOISEMARGIN $@87 int_number int_number ';'  */
#line 4856 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setInMargin((yyvsp[-2].dval), (yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2056, "MACRO INPUTNOISEMARGIN statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 10077 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 631: /* macro_pin_option: K_CAPACITANCE int_number ';'  */
#line 4865 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setCapacitance((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2057, "MACRO CAPACITANCE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 10090 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 632: /* macro_pin_option: K_MAXDELAY int_number ';'  */
#line 4874 "lef/lef.y"
    { if (lefCallbacks->PinCbk) lefData->lefrPin.setMaxdelay((yyvsp[-1].dval)); }
#line 10096 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 633: /* macro_pin_option: K_MAXLOAD int_number ';'  */
#line 4876 "lef/lef.y"
    { if (lefCallbacks->PinCbk) lefData->lefrPin.setMaxload((yyvsp[-1].dval)); }
#line 10102 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 634: /* macro_pin_option: K_RESISTANCE int_number ';'  */
#line 4878 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setResistance((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2058, "MACRO RESISTANCE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 10115 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 635: /* macro_pin_option: K_PULLDOWNRES int_number ';'  */
#line 4887 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setPulldownres((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2059, "MACRO PULLDOWNRES statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 10128 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 636: /* macro_pin_option: K_CURRENTSOURCE K_ACTIVE ';'  */
#line 4896 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setCurrentSource("ACTIVE");
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2060, "MACRO CURRENTSOURCE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 10141 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 637: /* macro_pin_option: K_CURRENTSOURCE K_RESISTIVE ';'  */
#line 4905 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setCurrentSource("RESISTIVE");
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2061, "MACRO CURRENTSOURCE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 10154 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 638: /* macro_pin_option: K_RISEVOLTAGETHRESHOLD int_number ';'  */
#line 4914 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setRiseVoltage((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2062, "MACRO RISEVOLTAGETHRESHOLD statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 10167 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 639: /* macro_pin_option: K_FALLVOLTAGETHRESHOLD int_number ';'  */
#line 4923 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setFallVoltage((yyvsp[-1].dval));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2063, "MACRO FALLVOLTAGETHRESHOLD statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 10180 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 640: /* macro_pin_option: K_IV_TABLES T_STRING T_STRING ';'  */
#line 4932 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) lefData->lefrPin.setTables((yyvsp[-2].string), (yyvsp[-1].string));
      } else
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2064, "MACRO IV_TABLES statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 10193 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 641: /* macro_pin_option: K_TAPERRULE T_STRING ';'  */
#line 4941 "lef/lef.y"
    { if (lefCallbacks->PinCbk) lefData->lefrPin.setTaperRule((yyvsp[-1].string)); }
#line 10199 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 642: /* $@88: %empty  */
#line 4942 "lef/lef.y"
               {lefData->lefDumbMode = 1000000; }
#line 10205 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 643: /* macro_pin_option: K_PROPERTY $@88 pin_prop_list ';'  */
#line 4943 "lef/lef.y"
    { lefData->lefDumbMode = 0;
    }
#line 10212 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 644: /* macro_pin_option: start_macro_port macro_port_class_option geometries K_END  */
#line 4946 "lef/lef.y"
    {
      lefData->lefDumbMode = 0;
      lefData->hasGeoLayer = 0;
      if (lefCallbacks->PinCbk) {
        lefData->lefrPin.addPort(lefData->lefrGeometriesPtr);
        lefData->lefrGeometriesPtr = 0;
        lefData->lefrDoGeometries = 0;
      }
      if ((lefData->needGeometry) && (lefData->needGeometry != 2))  // if the lefData->last LAYER in PORT
        if (lefCallbacks->PinCbk) // write warning only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings)
             lefWarning(2065, "Either PATH, RECT or POLYGON statement is a required in MACRO/PIN/PORT.");
    }
#line 10230 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 645: /* macro_pin_option: start_macro_port K_END  */
#line 4960 "lef/lef.y"
    {
      // Since in start_macro_port it has call the Init method, here
      // we need to call the Destroy method.
      // Still add a null pointer to set the number of port
      if (lefCallbacks->PinCbk) {
        lefData->lefrPin.addPort(lefData->lefrGeometriesPtr);
        lefData->lefrGeometriesPtr = 0;
        lefData->lefrDoGeometries = 0;
      }
      lefData->hasGeoLayer = 0;
    }
#line 10246 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 646: /* macro_pin_option: K_ANTENNASIZE int_number opt_layer_name ';'  */
#line 4972 "lef/lef.y"
    {  // a pre 5.4 syntax 
      lefData->use5_3 = 1;
      if (lefData->ignoreVersion) {
        // do nothing 
      } else if (lefData->versionNum >= 5.4) {
        if (lefData->use5_4) {
           if (lefCallbacks->PinCbk) { // write error only if cbk is set 
             if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
                lefData->outMsg = (char*)lefMalloc(10000);
                sprintf (lefData->outMsg,
                   "ANTENNASIZE statement is a version 5.3 and earlier syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
                lefError(1644, lefData->outMsg);
                lefFree(lefData->outMsg);
                CHKERR();
             }
           }
        }
      }
      if (lefCallbacks->PinCbk) lefData->lefrPin.addAntennaSize((yyvsp[-2].dval), (yyvsp[-1].string));
    }
#line 10271 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 647: /* macro_pin_option: K_ANTENNAMETALAREA NUMBER opt_layer_name ';'  */
#line 4993 "lef/lef.y"
    {  // a pre 5.4 syntax 
      lefData->use5_3 = 1;
      if (lefData->ignoreVersion) {
        // do nothing 
      } else if (lefData->versionNum >= 5.4) {
        if (lefData->use5_4) {
           if (lefCallbacks->PinCbk) { // write error only if cbk is set 
              if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
                 lefData->outMsg = (char*)lefMalloc(10000);
                 sprintf (lefData->outMsg,
                    "ANTENNAMETALAREA statement is a version 5.3 and earlier syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
                 lefError(1645, lefData->outMsg);
                 lefFree(lefData->outMsg);
                 CHKERR();
              }
           }
        }
      }
      if (lefCallbacks->PinCbk) lefData->lefrPin.addAntennaMetalArea((yyvsp[-2].dval), (yyvsp[-1].string));
    }
#line 10296 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 648: /* macro_pin_option: K_ANTENNAMETALLENGTH int_number opt_layer_name ';'  */
#line 5014 "lef/lef.y"
    { // a pre 5.4 syntax  
      lefData->use5_3 = 1;
      if (lefData->ignoreVersion) {
        // do nothing 
      } else if (lefData->versionNum >= 5.4) {
        if (lefData->use5_4) {
           if (lefCallbacks->PinCbk) { // write error only if cbk is set 
              if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
                 lefData->outMsg = (char*)lefMalloc(10000);
                 sprintf (lefData->outMsg,
                    "ANTENNAMETALLENGTH statement is a version 5.3 and earlier syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
                 lefError(1646, lefData->outMsg);
                 lefFree(lefData->outMsg);
                 CHKERR();
              }
           }
        }
      }
      if (lefCallbacks->PinCbk) lefData->lefrPin.addAntennaMetalLength((yyvsp[-2].dval), (yyvsp[-1].string));
    }
#line 10321 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 649: /* macro_pin_option: K_RISESLEWLIMIT int_number ';'  */
#line 5035 "lef/lef.y"
    { if (lefCallbacks->PinCbk) lefData->lefrPin.setRiseSlewLimit((yyvsp[-1].dval)); }
#line 10327 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 650: /* macro_pin_option: K_FALLSLEWLIMIT int_number ';'  */
#line 5037 "lef/lef.y"
    { if (lefCallbacks->PinCbk) lefData->lefrPin.setFallSlewLimit((yyvsp[-1].dval)); }
#line 10333 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 651: /* macro_pin_option: K_ANTENNAPARTIALMETALAREA NUMBER opt_layer_name ';'  */
#line 5039 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (lefData->ignoreVersion) {
        // do nothing 
      } else if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAPARTIALMETALAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1647, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      } else if (lefData->use5_3) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAPARTIALMETALAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1647, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
      if (lefCallbacks->PinCbk) lefData->lefrPin.addAntennaPartialMetalArea((yyvsp[-2].dval), (yyvsp[-1].string));
    }
#line 10367 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 652: /* macro_pin_option: K_ANTENNAPARTIALMETALSIDEAREA NUMBER opt_layer_name ';'  */
#line 5069 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (lefData->ignoreVersion) {
        // do nothing 
      } else if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAPARTIALMETALSIDEAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1648, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      } else if (lefData->use5_3) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAPARTIALMETALSIDEAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1648, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
      if (lefCallbacks->PinCbk) lefData->lefrPin.addAntennaPartialMetalSideArea((yyvsp[-2].dval), (yyvsp[-1].string));
    }
#line 10401 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 653: /* macro_pin_option: K_ANTENNAPARTIALCUTAREA NUMBER opt_layer_name ';'  */
#line 5099 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (lefData->ignoreVersion) {
        // do nothing 
      } else if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAPARTIALCUTAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1649, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      } else if (lefData->use5_3) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAPARTIALCUTAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1649, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
      if (lefCallbacks->PinCbk) lefData->lefrPin.addAntennaPartialCutArea((yyvsp[-2].dval), (yyvsp[-1].string));
    }
#line 10435 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 654: /* macro_pin_option: K_ANTENNADIFFAREA NUMBER opt_layer_name ';'  */
#line 5129 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (lefData->ignoreVersion) {
        // do nothing 
      } else if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNADIFFAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1650, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      } else if (lefData->use5_3) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNADIFFAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1650, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
      if (lefCallbacks->PinCbk) lefData->lefrPin.addAntennaDiffArea((yyvsp[-2].dval), (yyvsp[-1].string));
    }
#line 10469 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 655: /* macro_pin_option: K_ANTENNAGATEAREA NUMBER opt_layer_name ';'  */
#line 5159 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (lefData->ignoreVersion) {
        // do nothing 
      } else if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAGATEAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1651, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      } else if (lefData->use5_3) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAGATEAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1651, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
      if (lefCallbacks->PinCbk) lefData->lefrPin.addAntennaGateArea((yyvsp[-2].dval), (yyvsp[-1].string));
    }
#line 10503 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 656: /* macro_pin_option: K_ANTENNAMAXAREACAR NUMBER req_layer_name ';'  */
#line 5189 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (lefData->ignoreVersion) {
        // do nothing 
      } else if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAMAXAREACAR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1652, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      } else if (lefData->use5_3) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAMAXAREACAR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1652, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
      if (lefCallbacks->PinCbk) lefData->lefrPin.addAntennaMaxAreaCar((yyvsp[-2].dval), (yyvsp[-1].string));
    }
#line 10537 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 657: /* macro_pin_option: K_ANTENNAMAXSIDEAREACAR NUMBER req_layer_name ';'  */
#line 5219 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (lefData->ignoreVersion) {
        // do nothing 
      } else if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAMAXSIDEAREACAR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1653, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      } else if (lefData->use5_3) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAMAXSIDEAREACAR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1653, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
      if (lefCallbacks->PinCbk) lefData->lefrPin.addAntennaMaxSideAreaCar((yyvsp[-2].dval), (yyvsp[-1].string));
    }
#line 10571 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 658: /* macro_pin_option: K_ANTENNAMAXCUTCAR NUMBER req_layer_name ';'  */
#line 5249 "lef/lef.y"
    { // 5.4 syntax 
      lefData->use5_4 = 1;
      if (lefData->ignoreVersion) {
        // do nothing 
      } else if (lefData->versionNum < 5.4) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAMAXCUTCAR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1654, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      } else if (lefData->use5_3) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAMAXCUTCAR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1654, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
      if (lefCallbacks->PinCbk) lefData->lefrPin.addAntennaMaxCutCar((yyvsp[-2].dval), (yyvsp[-1].string));
    }
#line 10605 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 659: /* $@89: %empty  */
#line 5279 "lef/lef.y"
    { // 5.5 syntax 
      lefData->use5_4 = 1;
      if (lefData->ignoreVersion) {
        // do nothing 
      } else if (lefData->versionNum < 5.5) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAMODEL statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1655, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      } else if (lefData->use5_3) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ANTENNAMODEL statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1655, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      }
    }
#line 10638 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 661: /* $@90: %empty  */
#line 5308 "lef/lef.y"
              {lefData->lefDumbMode = 2; lefData->lefNoNum = 2; }
#line 10644 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 662: /* macro_pin_option: K_NETEXPR $@90 QSTRING ';'  */
#line 5309 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "NETEXPR statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1656, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      } else
        if (lefCallbacks->PinCbk) lefData->lefrPin.setNetExpr((yyvsp[-1].string));
    }
#line 10664 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 663: /* $@91: %empty  */
#line 5324 "lef/lef.y"
                        {lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 10670 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 664: /* macro_pin_option: K_SUPPLYSENSITIVITY $@91 T_STRING ';'  */
#line 5325 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "SUPPLYSENSITIVITY statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1657, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      } else
        if (lefCallbacks->PinCbk) lefData->lefrPin.setSupplySensitivity((yyvsp[-1].string));
    }
#line 10690 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 665: /* $@92: %empty  */
#line 5340 "lef/lef.y"
                        {lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 10696 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 666: /* macro_pin_option: K_GROUNDSENSITIVITY $@92 T_STRING ';'  */
#line 5341 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->PinCbk) { // write error only if cbk is set 
           if (lefData->pinWarnings++ < lefSettings->PinWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "GROUNDSENSITIVITY statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1658, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      } else
        if (lefCallbacks->PinCbk) lefData->lefrPin.setGroundSensitivity((yyvsp[-1].string));
    }
#line 10716 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 667: /* pin_layer_oxide: K_OXIDE1  */
#line 5359 "lef/lef.y"
    {
    if (lefCallbacks->PinCbk)
       lefData->lefrPin.addAntennaModel(1);
    }
#line 10725 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 668: /* pin_layer_oxide: K_OXIDE2  */
#line 5364 "lef/lef.y"
    {
    if (lefCallbacks->PinCbk)
       lefData->lefrPin.addAntennaModel(2);
    }
#line 10734 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 669: /* pin_layer_oxide: K_OXIDE3  */
#line 5369 "lef/lef.y"
    {
    if (lefCallbacks->PinCbk)
       lefData->lefrPin.addAntennaModel(3);
    }
#line 10743 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 670: /* pin_layer_oxide: K_OXIDE4  */
#line 5374 "lef/lef.y"
    {
    if (lefCallbacks->PinCbk)
       lefData->lefrPin.addAntennaModel(4);
    }
#line 10752 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 673: /* pin_name_value_pair: T_STRING NUMBER  */
#line 5386 "lef/lef.y"
    { 
      char temp[32];
      sprintf(temp, "%.11g", (yyvsp[0].dval));
      if (lefCallbacks->PinCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrPinProp.propType((yyvsp[-1].string));
         lefData->lefrPin.setNumProperty((yyvsp[-1].string), (yyvsp[0].dval), temp, propTp);
      }
    }
#line 10766 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 674: /* pin_name_value_pair: T_STRING QSTRING  */
#line 5396 "lef/lef.y"
    {
      if (lefCallbacks->PinCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrPinProp.propType((yyvsp[-1].string));
         lefData->lefrPin.setProperty((yyvsp[-1].string), (yyvsp[0].string), propTp);
      }
    }
#line 10778 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 675: /* pin_name_value_pair: T_STRING T_STRING  */
#line 5404 "lef/lef.y"
    {
      if (lefCallbacks->PinCbk) {
         char propTp;
         propTp = lefSettings->lefProps.lefrPinProp.propType((yyvsp[-1].string));
         lefData->lefrPin.setProperty((yyvsp[-1].string), (yyvsp[0].string), propTp);
      }
    }
#line 10790 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 676: /* electrical_direction: K_DIRECTION K_INPUT ';'  */
#line 5413 "lef/lef.y"
                                       {(yyval.string) = (char*)"INPUT";}
#line 10796 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 677: /* electrical_direction: K_DIRECTION K_OUTPUT ';'  */
#line 5414 "lef/lef.y"
                                        {(yyval.string) = (char*)"OUTPUT";}
#line 10802 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 678: /* electrical_direction: K_DIRECTION K_OUTPUT K_TRISTATE ';'  */
#line 5415 "lef/lef.y"
                                        {(yyval.string) = (char*)"OUTPUT TRISTATE";}
#line 10808 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 679: /* electrical_direction: K_DIRECTION K_INOUT ';'  */
#line 5416 "lef/lef.y"
                                        {(yyval.string) = (char*)"INOUT";}
#line 10814 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 680: /* electrical_direction: K_DIRECTION K_FEEDTHRU ';'  */
#line 5417 "lef/lef.y"
                                        {(yyval.string) = (char*)"FEEDTHRU";}
#line 10820 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 681: /* start_macro_port: K_PORT  */
#line 5420 "lef/lef.y"
    {
      if (lefCallbacks->PinCbk) {
        lefData->lefrDoGeometries = 1;
        lefData->hasPRP = 0;
        lefData->lefrGeometriesPtr = (lefiGeometries*)lefMalloc( sizeof(lefiGeometries));
        lefData->lefrGeometriesPtr->Init();
      }
      lefData->needGeometry = 0;  // don't need rect/path/poly define yet
      lefData->hasGeoLayer = 0;   // make sure LAYER is set before geometry
    }
#line 10835 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 683: /* macro_port_class_option: K_CLASS class_type ';'  */
#line 5433 "lef/lef.y"
    { if (lefData->lefrDoGeometries)
        lefData->lefrGeometriesPtr->addClass((yyvsp[-1].string)); }
#line 10842 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 684: /* macro_pin_use: K_SIGNAL  */
#line 5437 "lef/lef.y"
                {(yyval.string) = (char*)"SIGNAL";}
#line 10848 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 685: /* macro_pin_use: K_ANALOG  */
#line 5438 "lef/lef.y"
                {(yyval.string) = (char*)"ANALOG";}
#line 10854 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 686: /* macro_pin_use: K_POWER  */
#line 5439 "lef/lef.y"
                {(yyval.string) = (char*)"POWER";}
#line 10860 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 687: /* macro_pin_use: K_GROUND  */
#line 5440 "lef/lef.y"
                {(yyval.string) = (char*)"GROUND";}
#line 10866 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 688: /* macro_pin_use: K_CLOCK  */
#line 5441 "lef/lef.y"
                {(yyval.string) = (char*)"CLOCK";}
#line 10872 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 689: /* macro_pin_use: K_DATA  */
#line 5442 "lef/lef.y"
                {(yyval.string) = (char*)"DATA";}
#line 10878 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 690: /* macro_scan_use: K_INPUT  */
#line 5445 "lef/lef.y"
          {(yyval.string) = (char*)"INPUT";}
#line 10884 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 691: /* macro_scan_use: K_OUTPUT  */
#line 5446 "lef/lef.y"
                {(yyval.string) = (char*)"OUTPUT";}
#line 10890 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 692: /* macro_scan_use: K_START  */
#line 5447 "lef/lef.y"
                {(yyval.string) = (char*)"START";}
#line 10896 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 693: /* macro_scan_use: K_STOP  */
#line 5448 "lef/lef.y"
                {(yyval.string) = (char*)"STOP";}
#line 10902 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 694: /* pin_shape: %empty  */
#line 5451 "lef/lef.y"
  {(yyval.string) = (char*)""; }
#line 10908 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 695: /* pin_shape: K_ABUTMENT  */
#line 5452 "lef/lef.y"
                {(yyval.string) = (char*)"ABUTMENT";}
#line 10914 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 696: /* pin_shape: K_RING  */
#line 5453 "lef/lef.y"
                {(yyval.string) = (char*)"RING";}
#line 10920 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 697: /* pin_shape: K_FEEDTHRU  */
#line 5454 "lef/lef.y"
                {(yyval.string) = (char*)"FEEDTHRU";}
#line 10926 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 699: /* $@93: %empty  */
#line 5459 "lef/lef.y"
          {lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 10932 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 700: /* $@94: %empty  */
#line 5460 "lef/lef.y"
    {
      if ((lefData->needGeometry) && (lefData->needGeometry != 2)) // 1 LAYER follow after another
        if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
          // geometries is called by MACRO/OBS & MACRO/PIN/PORT 
          if (lefData->obsDef)
             lefWarning(2076, "Either PATH, RECT or POLYGON statement is a required in MACRO/OBS.");
          else
             lefWarning(2065, "Either PATH, RECT or POLYGON statement is a required in MACRO/PIN/PORT.");
        }
      if (lefData->lefrDoGeometries)
        lefData->lefrGeometriesPtr->addLayer((yyvsp[0].string));
      lefData->needGeometry = 1;    // within LAYER it requires either path, rect, poly
      lefData->hasGeoLayer = 1;
    }
#line 10951 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 702: /* geometry: K_WIDTH int_number ';'  */
#line 5477 "lef/lef.y"
    { 
      if (lefData->lefrDoGeometries) {
        if (lefData->hasGeoLayer == 0) {   // LAYER statement is missing 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              CHKERR();
           }
        } else
           lefData->lefrGeometriesPtr->addWidth((yyvsp[-1].dval)); 
      } 
    }
#line 10967 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 703: /* geometry: K_PATH maskColor firstPt otherPts ';'  */
#line 5489 "lef/lef.y"
    { if (lefData->lefrDoGeometries) {
        if (lefData->hasGeoLayer == 0) {   // LAYER statement is missing 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              CHKERR();
           }
        } else {
           if (lefData->versionNum < 5.8 && (int)(yyvsp[-3].integer) > 0) {
              if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
                 lefError(2083, "Color mask information can only be defined with version 5.8.");
                 CHKERR(); 
              }           
           } else {
                lefData->lefrGeometriesPtr->addPath((int)(yyvsp[-3].integer));
           }
        }
      }
      lefData->hasPRP = 1;
      lefData->needGeometry = 2;
    }
#line 10992 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 704: /* geometry: K_PATH maskColor K_ITERATE firstPt otherPts stepPattern ';'  */
#line 5510 "lef/lef.y"
    { if (lefData->lefrDoGeometries) {
        if (lefData->hasGeoLayer == 0) {   // LAYER statement is missing 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              CHKERR();
           }
        } else {
           if (lefData->versionNum < 5.8 && (int)(yyvsp[-5].integer) > 0) {
              if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
                 lefError(2083, "Color mask information can only be defined with version 5.8.");
                 CHKERR(); 
              }           
           } else {
              lefData->lefrGeometriesPtr->addPathIter((int)(yyvsp[-5].integer));
            }
         }
      } 
      lefData->hasPRP = 1;
      lefData->needGeometry = 2;
    }
#line 11017 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 705: /* geometry: K_RECT maskColor pt pt ';'  */
#line 5531 "lef/lef.y"
    { if (lefData->lefrDoGeometries) {
        if (lefData->hasGeoLayer == 0) {   // LAYER statement is missing 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              CHKERR();
           }
        } else {
           if (lefData->versionNum < 5.8 && (int)(yyvsp[-3].integer) > 0) {
              if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
                 lefError(2083, "Color mask information can only be defined with version 5.8.");
                 CHKERR(); 
              }           
           } else {
              lefData->lefrGeometriesPtr->addRect((int)(yyvsp[-3].integer), (yyvsp[-2].pt).x, (yyvsp[-2].pt).y, (yyvsp[-1].pt).x, (yyvsp[-1].pt).y);
           }
        }
      }
      lefData->needGeometry = 2;
    }
#line 11041 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 706: /* geometry: K_RECT maskColor K_ITERATE pt pt stepPattern ';'  */
#line 5551 "lef/lef.y"
    { if (lefData->lefrDoGeometries) {
        if (lefData->hasGeoLayer == 0) {   // LAYER statement is missing 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              CHKERR();
           }
        } else {
           if (lefData->versionNum < 5.8 && (int)(yyvsp[-5].integer) > 0) {
              if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
                 lefError(2083, "Color mask information can only be defined with version 5.8.");
                 CHKERR(); 
              }           
           } else {
              lefData->lefrGeometriesPtr->addRectIter((int)(yyvsp[-5].integer), (yyvsp[-3].pt).x, (yyvsp[-3].pt).y, (yyvsp[-2].pt).x, (yyvsp[-2].pt).y);
           }
        }
      }
      lefData->needGeometry = 2;
    }
#line 11065 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 707: /* geometry: K_POLYGON maskColor firstPt nextPt nextPt otherPts ';'  */
#line 5571 "lef/lef.y"
    {
      if (lefData->lefrDoGeometries) {
        if (lefData->hasGeoLayer == 0) {   // LAYER statement is missing 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              CHKERR();
           }
        } else {
           if (lefData->versionNum < 5.8 && (int)(yyvsp[-5].integer) > 0) {
              if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
                 lefError(2083, "Color mask information can only be defined with version 5.8.");
                 CHKERR(); 
              }           
           } else {
              lefData->lefrGeometriesPtr->addPolygon((int)(yyvsp[-5].integer));
            }
           }
      }
      lefData->hasPRP = 1;
      lefData->needGeometry = 2;
    }
#line 11091 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 708: /* geometry: K_POLYGON maskColor K_ITERATE firstPt nextPt nextPt otherPts stepPattern ';'  */
#line 5593 "lef/lef.y"
    { if (lefData->lefrDoGeometries) {
        if (lefData->hasGeoLayer == 0) {   // LAYER statement is missing 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              CHKERR();
           }
        } else {
           if (lefData->versionNum < 5.8 && (int)(yyvsp[-7].integer) > 0) {
              if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
                 lefError(2083, "Color mask information can only be defined with version 5.8.");
                 CHKERR(); 
              }           
           } else {
              lefData->lefrGeometriesPtr->addPolygonIter((int)(yyvsp[-7].integer));
           }
         }
      }
      lefData->hasPRP = 1;
      lefData->needGeometry = 2;
    }
#line 11116 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 709: /* geometry: via_placement  */
#line 5614 "lef/lef.y"
    { }
#line 11122 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 713: /* layer_exceptpgnet: K_EXCEPTPGNET  */
#line 5621 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
        lefData->outMsg = (char*)lefMalloc(10000);
        sprintf(lefData->outMsg,
          "EXCEPTPGNET is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
        lefError(1699, lefData->outMsg);
        lefFree(lefData->outMsg);
        CHKERR();
      } else {
       if (lefData->lefrDoGeometries)
        lefData->lefrGeometriesPtr->addLayerExceptPgNet();
      }
    }
#line 11140 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 715: /* layer_spacing: K_SPACING int_number  */
#line 5637 "lef/lef.y"
    { if (lefData->lefrDoGeometries) {
        if (zeroOrGt((yyvsp[0].dval)))
           lefData->lefrGeometriesPtr->addLayerMinSpacing((yyvsp[0].dval));
        else {
           lefData->outMsg = (char*)lefMalloc(10000);
           sprintf (lefData->outMsg,
              "THE SPACING statement has the value %g in MACRO OBS.\nValue has to be 0 or greater.", (yyvsp[0].dval));
           lefError(1659, lefData->outMsg);
           lefFree(lefData->outMsg);
           CHKERR();
        }
      }
    }
#line 11158 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 716: /* layer_spacing: K_DESIGNRULEWIDTH int_number  */
#line 5651 "lef/lef.y"
    { if (lefData->lefrDoGeometries) {
        if (zeroOrGt((yyvsp[0].dval)))
           lefData->lefrGeometriesPtr->addLayerRuleWidth((yyvsp[0].dval));
        else {
           lefData->outMsg = (char*)lefMalloc(10000);
           sprintf (lefData->outMsg,
              "THE DESIGNRULEWIDTH statement has the value %g in MACRO OBS.\nValue has to be 0 or greater.", (yyvsp[0].dval));
           lefError(1660, lefData->outMsg);
           lefFree(lefData->outMsg);
           CHKERR();
        }
      }
    }
#line 11176 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 717: /* firstPt: pt  */
#line 5666 "lef/lef.y"
    { if (lefData->lefrDoGeometries)
        lefData->lefrGeometriesPtr->startList((yyvsp[0].pt).x, (yyvsp[0].pt).y); }
#line 11183 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 718: /* nextPt: pt  */
#line 5670 "lef/lef.y"
    { if (lefData->lefrDoGeometries)
        lefData->lefrGeometriesPtr->addToList((yyvsp[0].pt).x, (yyvsp[0].pt).y); }
#line 11190 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 721: /* $@95: %empty  */
#line 5679 "lef/lef.y"
                     {lefData->lefDumbMode = 1;}
#line 11196 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 722: /* via_placement: K_VIA maskColor pt $@95 T_STRING ';'  */
#line 5680 "lef/lef.y"
    { 
        if (lefData->lefrDoGeometries){
            if (lefData->versionNum < 5.8 && (int)(yyvsp[-4].integer) > 0) {
              if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
                 lefError(2083, "Color mask information can only be defined with version 5.8.");
                 CHKERR(); 
              }           
            } else {
                lefData->lefrGeometriesPtr->addVia((int)(yyvsp[-4].integer), (yyvsp[-3].pt).x, (yyvsp[-3].pt).y, (yyvsp[-1].string));
            }
        }
    }
#line 11213 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 723: /* $@96: %empty  */
#line 5692 "lef/lef.y"
                                 {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 11219 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 724: /* via_placement: K_VIA K_ITERATE maskColor pt $@96 T_STRING stepPattern ';'  */
#line 5694 "lef/lef.y"
    { 
        if (lefData->lefrDoGeometries) {
            if (lefData->versionNum < 5.8 && (int)(yyvsp[-5].integer) > 0) {
              if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
                 lefError(2083, "Color mask information can only be defined with version 5.8.");
                 CHKERR(); 
              }           
            } else {
              lefData->lefrGeometriesPtr->addViaIter((int)(yyvsp[-5].integer), (yyvsp[-4].pt).x, (yyvsp[-4].pt).y, (yyvsp[-2].string)); 
            }
        }
    }
#line 11236 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 725: /* stepPattern: K_DO int_number K_BY int_number K_STEP int_number int_number  */
#line 5708 "lef/lef.y"
     { if (lefData->lefrDoGeometries)
         lefData->lefrGeometriesPtr->addStepPattern((yyvsp[-5].dval), (yyvsp[-3].dval), (yyvsp[-1].dval), (yyvsp[0].dval)); }
#line 11243 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 726: /* sitePattern: T_STRING int_number int_number orientation K_DO int_number K_BY int_number K_STEP int_number int_number  */
#line 5713 "lef/lef.y"
    {
      if (lefData->lefrDoSite) {
        lefData->lefrSitePatternPtr = (lefiSitePattern*)lefMalloc(
                                   sizeof(lefiSitePattern));
        lefData->lefrSitePatternPtr->Init();
        lefData->lefrSitePatternPtr->set((yyvsp[-10].string), (yyvsp[-9].dval), (yyvsp[-8].dval), (yyvsp[-7].integer), (yyvsp[-5].dval), (yyvsp[-3].dval),
          (yyvsp[-1].dval), (yyvsp[0].dval));
        }
    }
#line 11257 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 727: /* sitePattern: T_STRING int_number int_number orientation  */
#line 5723 "lef/lef.y"
    {
      if (lefData->lefrDoSite) {
        lefData->lefrSitePatternPtr = (lefiSitePattern*)lefMalloc(
                                   sizeof(lefiSitePattern));
        lefData->lefrSitePatternPtr->Init();
        lefData->lefrSitePatternPtr->set((yyvsp[-3].string), (yyvsp[-2].dval), (yyvsp[-1].dval), (yyvsp[0].integer), -1, -1,
          -1, -1);
        }
    }
#line 11271 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 728: /* $@97: %empty  */
#line 5735 "lef/lef.y"
    { 
      if (lefData->lefrDoTrack) {
        lefData->lefrTrackPatternPtr = (lefiTrackPattern*)lefMalloc(
                                sizeof(lefiTrackPattern));
        lefData->lefrTrackPatternPtr->Init();
        lefData->lefrTrackPatternPtr->set("X", (yyvsp[-4].dval), (int)(yyvsp[-2].dval), (yyvsp[0].dval));
      }    
    }
#line 11284 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 729: /* $@98: %empty  */
#line 5743 "lef/lef.y"
            {lefData->lefDumbMode = 1000000000;}
#line 11290 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 730: /* trackPattern: K_X int_number K_DO int_number K_STEP int_number $@97 K_LAYER $@98 trackLayers  */
#line 5744 "lef/lef.y"
    { lefData->lefDumbMode = 0;}
#line 11296 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 731: /* $@99: %empty  */
#line 5746 "lef/lef.y"
    { 
      if (lefData->lefrDoTrack) {
        lefData->lefrTrackPatternPtr = (lefiTrackPattern*)lefMalloc(
                                    sizeof(lefiTrackPattern));
        lefData->lefrTrackPatternPtr->Init();
        lefData->lefrTrackPatternPtr->set("Y", (yyvsp[-4].dval), (int)(yyvsp[-2].dval), (yyvsp[0].dval));
      }    
    }
#line 11309 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 732: /* $@100: %empty  */
#line 5754 "lef/lef.y"
            {lefData->lefDumbMode = 1000000000;}
#line 11315 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 733: /* trackPattern: K_Y int_number K_DO int_number K_STEP int_number $@99 K_LAYER $@100 trackLayers  */
#line 5755 "lef/lef.y"
    { lefData->lefDumbMode = 0;}
#line 11321 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 734: /* trackPattern: K_X int_number K_DO int_number K_STEP int_number  */
#line 5757 "lef/lef.y"
    { 
      if (lefData->lefrDoTrack) {
        lefData->lefrTrackPatternPtr = (lefiTrackPattern*)lefMalloc(
                                    sizeof(lefiTrackPattern));
        lefData->lefrTrackPatternPtr->Init();
        lefData->lefrTrackPatternPtr->set("X", (yyvsp[-4].dval), (int)(yyvsp[-2].dval), (yyvsp[0].dval));
      }    
    }
#line 11334 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 735: /* trackPattern: K_Y int_number K_DO int_number K_STEP int_number  */
#line 5766 "lef/lef.y"
    { 
      if (lefData->lefrDoTrack) {
        lefData->lefrTrackPatternPtr = (lefiTrackPattern*)lefMalloc(
                                    sizeof(lefiTrackPattern));
        lefData->lefrTrackPatternPtr->Init();
        lefData->lefrTrackPatternPtr->set("Y", (yyvsp[-4].dval), (int)(yyvsp[-2].dval), (yyvsp[0].dval));
      }    
    }
#line 11347 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 738: /* layer_name: T_STRING  */
#line 5781 "lef/lef.y"
    { if (lefData->lefrDoTrack) lefData->lefrTrackPatternPtr->addLayer((yyvsp[0].string)); }
#line 11353 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 739: /* gcellPattern: K_X int_number K_DO int_number K_STEP int_number  */
#line 5784 "lef/lef.y"
    {
      if (lefData->lefrDoGcell) {
        lefData->lefrGcellPatternPtr = (lefiGcellPattern*)lefMalloc(
                                    sizeof(lefiGcellPattern));
        lefData->lefrGcellPatternPtr->Init();
        lefData->lefrGcellPatternPtr->set("X", (yyvsp[-4].dval), (int)(yyvsp[-2].dval), (yyvsp[0].dval));
      }    
    }
#line 11366 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 740: /* gcellPattern: K_Y int_number K_DO int_number K_STEP int_number  */
#line 5793 "lef/lef.y"
    {
      if (lefData->lefrDoGcell) {
        lefData->lefrGcellPatternPtr = (lefiGcellPattern*)lefMalloc(
                                    sizeof(lefiGcellPattern));
        lefData->lefrGcellPatternPtr->Init();
        lefData->lefrGcellPatternPtr->set("Y", (yyvsp[-4].dval), (int)(yyvsp[-2].dval), (yyvsp[0].dval));
      }    
    }
#line 11379 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 741: /* macro_obs: start_macro_obs geometries K_END  */
#line 5803 "lef/lef.y"
    { 
      if (lefCallbacks->ObstructionCbk) {
        lefData->lefrObstruction.setGeometries(lefData->lefrGeometriesPtr);
        lefData->lefrGeometriesPtr = 0;
        lefData->lefrDoGeometries = 0;
        CALLBACK(lefCallbacks->ObstructionCbk, lefrObstructionCbkType, &lefData->lefrObstruction);
      }
      lefData->lefDumbMode = 0;
      lefData->hasGeoLayer = 0;       // reset 
    }
#line 11394 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 742: /* macro_obs: start_macro_obs K_END  */
#line 5814 "lef/lef.y"
    {
       // The pointer has malloced in start, need to free manually 
       if (lefData->lefrGeometriesPtr) {
          lefData->lefrGeometriesPtr->Destroy();
          lefFree(lefData->lefrGeometriesPtr);
          lefData->lefrGeometriesPtr = 0;
          lefData->lefrDoGeometries = 0;
       }
       lefData->hasGeoLayer = 0;
    }
#line 11409 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 743: /* start_macro_obs: K_OBS  */
#line 5826 "lef/lef.y"
    {
      lefData->obsDef = 1;
      if (lefCallbacks->ObstructionCbk) {
        lefData->lefrDoGeometries = 1;
        lefData->lefrGeometriesPtr = (lefiGeometries*)lefMalloc(
            sizeof(lefiGeometries));
        lefData->lefrGeometriesPtr->Init();
        }
      lefData->hasGeoLayer = 0;
    }
#line 11424 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 744: /* macro_density: K_DENSITY density_layer density_layers K_END  */
#line 5838 "lef/lef.y"
    { 
      if (lefData->versionNum < 5.6) {
        if (lefCallbacks->DensityCbk) { // write error only if cbk is set 
           if (lefData->macroWarnings++ < lefSettings->MacroWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "DENSITY statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1661, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
      } 
      if (lefCallbacks->DensityCbk) {
        CALLBACK(lefCallbacks->DensityCbk, lefrDensityCbkType, &lefData->lefrDensity);
        lefData->lefrDensity.clear();
      }
      lefData->lefDumbMode = 0;
    }
#line 11448 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 747: /* $@101: %empty  */
#line 5862 "lef/lef.y"
                       { lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 11454 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 748: /* $@102: %empty  */
#line 5863 "lef/lef.y"
    {
      if (lefCallbacks->DensityCbk)
        lefData->lefrDensity.addLayer((yyvsp[-1].string));
    }
#line 11463 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 752: /* density_layer_rect: K_RECT pt pt int_number ';'  */
#line 5874 "lef/lef.y"
    {
      if (lefCallbacks->DensityCbk)
        lefData->lefrDensity.addRect((yyvsp[-3].pt).x, (yyvsp[-3].pt).y, (yyvsp[-2].pt).x, (yyvsp[-2].pt).y, (yyvsp[-1].dval)); 
    }
#line 11472 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 753: /* $@103: %empty  */
#line 5879 "lef/lef.y"
                             { lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 11478 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 754: /* macro_clocktype: K_CLOCKTYPE $@103 T_STRING ';'  */
#line 5880 "lef/lef.y"
    { if (lefCallbacks->MacroCbk) lefData->lefrMacro.setClockType((yyvsp[-1].string)); }
#line 11484 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 755: /* timing: start_timing timing_options end_timing  */
#line 5883 "lef/lef.y"
    { }
#line 11490 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 756: /* start_timing: K_TIMING  */
#line 5886 "lef/lef.y"
    { }
#line 11496 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 757: /* end_timing: K_END K_TIMING  */
#line 5889 "lef/lef.y"
  {
    if (lefData->versionNum < 5.4) {
      if (lefCallbacks->TimingCbk && lefData->lefrTiming.hasData())
        CALLBACK(lefCallbacks->TimingCbk, lefrTimingCbkType, &lefData->lefrTiming);
      lefData->lefrTiming.clear();
    } else {
      if (lefCallbacks->TimingCbk) // write warning only if cbk is set 
        if (lefData->timingWarnings++ < lefSettings->TimingWarnings)
          lefWarning(2066, "MACRO TIMING statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
      lefData->lefrTiming.clear();
    }
  }
#line 11513 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 760: /* $@104: %empty  */
#line 5909 "lef/lef.y"
    {
    if (lefData->versionNum < 5.4) {
      if (lefCallbacks->TimingCbk && lefData->lefrTiming.hasData())
        CALLBACK(lefCallbacks->TimingCbk, lefrTimingCbkType, &lefData->lefrTiming);
    }
    lefData->lefDumbMode = 1000000000;
    lefData->lefrTiming.clear();
    }
#line 11526 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 761: /* timing_option: K_FROMPIN $@104 list_of_from_strings ';'  */
#line 5918 "lef/lef.y"
    { lefData->lefDumbMode = 0;}
#line 11532 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 762: /* $@105: %empty  */
#line 5919 "lef/lef.y"
            {lefData->lefDumbMode = 1000000000;}
#line 11538 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 763: /* timing_option: K_TOPIN $@105 list_of_to_strings ';'  */
#line 5920 "lef/lef.y"
    { lefData->lefDumbMode = 0;}
#line 11544 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 764: /* $@106: %empty  */
#line 5922 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addRiseFall((yyvsp[-3].string),(yyvsp[-1].dval),(yyvsp[0].dval)); }
#line 11550 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 765: /* timing_option: risefall K_INTRINSIC int_number int_number $@106 slew_spec K_VARIABLE int_number int_number ';'  */
#line 5924 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addRiseFallVariable((yyvsp[-2].dval),(yyvsp[-1].dval)); }
#line 11556 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 766: /* timing_option: risefall delay_or_transition K_UNATENESS unateness K_TABLEDIMENSION int_number int_number int_number ';'  */
#line 5927 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) {
        if ((yyvsp[-7].string)[0] == 'D' || (yyvsp[-7].string)[0] == 'd') // delay 
          lefData->lefrTiming.addDelay((yyvsp[-8].string), (yyvsp[-5].string), (yyvsp[-3].dval), (yyvsp[-2].dval), (yyvsp[-1].dval));
        else
          lefData->lefrTiming.addTransition((yyvsp[-8].string), (yyvsp[-5].string), (yyvsp[-3].dval), (yyvsp[-2].dval), (yyvsp[-1].dval));
      }
    }
#line 11568 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 767: /* timing_option: K_TABLEAXIS list_of_table_axis_dnumbers ';'  */
#line 5935 "lef/lef.y"
    { }
#line 11574 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 768: /* timing_option: K_TABLEENTRIES list_of_table_entries ';'  */
#line 5937 "lef/lef.y"
    { }
#line 11580 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 769: /* timing_option: K_RISERS int_number int_number ';'  */
#line 5939 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.setRiseRS((yyvsp[-2].dval),(yyvsp[-1].dval)); }
#line 11586 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 770: /* timing_option: K_FALLRS int_number int_number ';'  */
#line 5941 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.setFallRS((yyvsp[-2].dval),(yyvsp[-1].dval)); }
#line 11592 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 771: /* timing_option: K_RISECS int_number int_number ';'  */
#line 5943 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.setRiseCS((yyvsp[-2].dval),(yyvsp[-1].dval)); }
#line 11598 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 772: /* timing_option: K_FALLCS int_number int_number ';'  */
#line 5945 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.setFallCS((yyvsp[-2].dval),(yyvsp[-1].dval)); }
#line 11604 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 773: /* timing_option: K_RISESATT1 int_number int_number ';'  */
#line 5947 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.setRiseAtt1((yyvsp[-2].dval),(yyvsp[-1].dval)); }
#line 11610 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 774: /* timing_option: K_FALLSATT1 int_number int_number ';'  */
#line 5949 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.setFallAtt1((yyvsp[-2].dval),(yyvsp[-1].dval)); }
#line 11616 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 775: /* timing_option: K_RISET0 int_number int_number ';'  */
#line 5951 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.setRiseTo((yyvsp[-2].dval),(yyvsp[-1].dval)); }
#line 11622 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 776: /* timing_option: K_FALLT0 int_number int_number ';'  */
#line 5953 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.setFallTo((yyvsp[-2].dval),(yyvsp[-1].dval)); }
#line 11628 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 777: /* timing_option: K_UNATENESS unateness ';'  */
#line 5955 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addUnateness((yyvsp[-1].string)); }
#line 11634 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 778: /* timing_option: K_STABLE K_SETUP int_number K_HOLD int_number risefall ';'  */
#line 5957 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.setStable((yyvsp[-4].dval),(yyvsp[-2].dval),(yyvsp[-1].string)); }
#line 11640 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 779: /* timing_option: two_pin_trigger from_pin_trigger to_pin_trigger K_TABLEDIMENSION int_number int_number int_number ';'  */
#line 5959 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addSDF2Pins((yyvsp[-7].string),(yyvsp[-6].string),(yyvsp[-5].string),(yyvsp[-3].dval),(yyvsp[-2].dval),(yyvsp[-1].dval)); }
#line 11646 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 780: /* timing_option: one_pin_trigger K_TABLEDIMENSION int_number int_number int_number ';'  */
#line 5961 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addSDF1Pin((yyvsp[-5].string),(yyvsp[-3].dval),(yyvsp[-2].dval),(yyvsp[-2].dval)); }
#line 11652 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 781: /* timing_option: K_SDFCONDSTART QSTRING ';'  */
#line 5963 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.setSDFcondStart((yyvsp[-1].string)); }
#line 11658 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 782: /* timing_option: K_SDFCONDEND QSTRING ';'  */
#line 5965 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.setSDFcondEnd((yyvsp[-1].string)); }
#line 11664 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 783: /* timing_option: K_SDFCOND QSTRING ';'  */
#line 5967 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.setSDFcond((yyvsp[-1].string)); }
#line 11670 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 784: /* timing_option: K_EXTENSION ';'  */
#line 5969 "lef/lef.y"
    { }
#line 11676 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 785: /* one_pin_trigger: K_MPWH  */
#line 5973 "lef/lef.y"
    { (yyval.string) = (char*)"MPWH";}
#line 11682 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 786: /* one_pin_trigger: K_MPWL  */
#line 5975 "lef/lef.y"
    { (yyval.string) = (char*)"MPWL";}
#line 11688 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 787: /* one_pin_trigger: K_PERIOD  */
#line 5977 "lef/lef.y"
    { (yyval.string) = (char*)"PERIOD";}
#line 11694 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 788: /* two_pin_trigger: K_SETUP  */
#line 5981 "lef/lef.y"
    { (yyval.string) = (char*)"SETUP";}
#line 11700 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 789: /* two_pin_trigger: K_HOLD  */
#line 5983 "lef/lef.y"
    { (yyval.string) = (char*)"HOLD";}
#line 11706 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 790: /* two_pin_trigger: K_RECOVERY  */
#line 5985 "lef/lef.y"
    { (yyval.string) = (char*)"RECOVERY";}
#line 11712 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 791: /* two_pin_trigger: K_SKEW  */
#line 5987 "lef/lef.y"
    { (yyval.string) = (char*)"SKEW";}
#line 11718 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 792: /* from_pin_trigger: K_ANYEDGE  */
#line 5991 "lef/lef.y"
    { (yyval.string) = (char*)"ANYEDGE";}
#line 11724 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 793: /* from_pin_trigger: K_POSEDGE  */
#line 5993 "lef/lef.y"
    { (yyval.string) = (char*)"POSEDGE";}
#line 11730 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 794: /* from_pin_trigger: K_NEGEDGE  */
#line 5995 "lef/lef.y"
    { (yyval.string) = (char*)"NEGEDGE";}
#line 11736 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 795: /* to_pin_trigger: K_ANYEDGE  */
#line 5999 "lef/lef.y"
    { (yyval.string) = (char*)"ANYEDGE";}
#line 11742 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 796: /* to_pin_trigger: K_POSEDGE  */
#line 6001 "lef/lef.y"
    { (yyval.string) = (char*)"POSEDGE";}
#line 11748 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 797: /* to_pin_trigger: K_NEGEDGE  */
#line 6003 "lef/lef.y"
    { (yyval.string) = (char*)"NEGEDGE";}
#line 11754 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 798: /* delay_or_transition: K_DELAY  */
#line 6007 "lef/lef.y"
    { (yyval.string) = (char*)"DELAY"; }
#line 11760 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 799: /* delay_or_transition: K_TRANSITIONTIME  */
#line 6009 "lef/lef.y"
    { (yyval.string) = (char*)"TRANSITION"; }
#line 11766 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 800: /* list_of_table_entries: table_entry  */
#line 6013 "lef/lef.y"
    { }
#line 11772 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 801: /* list_of_table_entries: list_of_table_entries table_entry  */
#line 6015 "lef/lef.y"
    { }
#line 11778 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 802: /* table_entry: '(' int_number int_number int_number ')'  */
#line 6018 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addTableEntry((yyvsp[-3].dval),(yyvsp[-2].dval),(yyvsp[-1].dval)); }
#line 11784 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 803: /* list_of_table_axis_dnumbers: int_number  */
#line 6022 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addTableAxisNumber((yyvsp[0].dval)); }
#line 11790 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 804: /* list_of_table_axis_dnumbers: list_of_table_axis_dnumbers int_number  */
#line 6024 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addTableAxisNumber((yyvsp[0].dval)); }
#line 11796 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 805: /* slew_spec: %empty  */
#line 6028 "lef/lef.y"
    { }
#line 11802 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 806: /* slew_spec: int_number int_number int_number int_number  */
#line 6030 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addRiseFallSlew((yyvsp[-3].dval),(yyvsp[-2].dval),(yyvsp[-1].dval),(yyvsp[0].dval)); }
#line 11808 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 807: /* slew_spec: int_number int_number int_number int_number int_number int_number int_number  */
#line 6032 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addRiseFallSlew((yyvsp[-6].dval),(yyvsp[-5].dval),(yyvsp[-4].dval),(yyvsp[-3].dval));
      if (lefCallbacks->TimingCbk) lefData->lefrTiming.addRiseFallSlew2((yyvsp[-2].dval),(yyvsp[-1].dval),(yyvsp[0].dval)); }
#line 11815 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 808: /* risefall: K_RISE  */
#line 6037 "lef/lef.y"
    { (yyval.string) = (char*)"RISE"; }
#line 11821 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 809: /* risefall: K_FALL  */
#line 6039 "lef/lef.y"
    { (yyval.string) = (char*)"FALL"; }
#line 11827 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 810: /* unateness: K_INVERT  */
#line 6043 "lef/lef.y"
    { (yyval.string) = (char*)"INVERT"; }
#line 11833 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 811: /* unateness: K_NONINVERT  */
#line 6045 "lef/lef.y"
    { (yyval.string) = (char*)"NONINVERT"; }
#line 11839 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 812: /* unateness: K_NONUNATE  */
#line 6047 "lef/lef.y"
    { (yyval.string) = (char*)"NONUNATE"; }
#line 11845 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 813: /* list_of_from_strings: T_STRING  */
#line 6051 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addFromPin((yyvsp[0].string)); }
#line 11851 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 814: /* list_of_from_strings: list_of_from_strings T_STRING  */
#line 6053 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addFromPin((yyvsp[0].string)); }
#line 11857 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 815: /* list_of_to_strings: T_STRING  */
#line 6057 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addToPin((yyvsp[0].string)); }
#line 11863 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 816: /* list_of_to_strings: list_of_to_strings T_STRING  */
#line 6059 "lef/lef.y"
    { if (lefCallbacks->TimingCbk) lefData->lefrTiming.addToPin((yyvsp[0].string)); }
#line 11869 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 817: /* $@107: %empty  */
#line 6062 "lef/lef.y"
    {
      if (lefCallbacks->ArrayCbk)
        CALLBACK(lefCallbacks->ArrayCbk, lefrArrayCbkType, &lefData->lefrArray);
      lefData->lefrArray.clear();
      lefData->lefrSitePatternPtr = 0;
      lefData->lefrDoSite = 0;
   }
#line 11881 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 819: /* $@108: %empty  */
#line 6071 "lef/lef.y"
                     {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 11887 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 820: /* start_array: K_ARRAY $@108 T_STRING  */
#line 6072 "lef/lef.y"
    {
      if (lefCallbacks->ArrayCbk) {
        lefData->lefrArray.setName((yyvsp[0].string));
        CALLBACK(lefCallbacks->ArrayBeginCbk, lefrArrayBeginCbkType, (yyvsp[0].string));
      }
      //strcpy(lefData->arrayName, $3);
      lefData->arrayName = strdup((yyvsp[0].string));
    }
#line 11900 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 821: /* $@109: %empty  */
#line 6081 "lef/lef.y"
                 {lefData->lefDumbMode = 1; lefData->lefNoNum = 1;}
#line 11906 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 822: /* end_array: K_END $@109 T_STRING  */
#line 6082 "lef/lef.y"
    {
      if (lefCallbacks->ArrayCbk && lefCallbacks->ArrayEndCbk)
        CALLBACK(lefCallbacks->ArrayEndCbk, lefrArrayEndCbkType, (yyvsp[0].string));
      if (strcmp(lefData->arrayName, (yyvsp[0].string)) != 0) {
        if (lefCallbacks->ArrayCbk) { // write error only if cbk is set 
           if (lefData->arrayWarnings++ < lefSettings->ArrayWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "END ARRAY name %s is different from the ARRAY name %s.\nCorrect the LEF file before rerunning it through the LEF parser.", (yyvsp[0].string), lefData->arrayName);
              lefError(1662, lefData->outMsg);
              lefFree(lefData->outMsg);
              lefFree(lefData->arrayName);
              CHKERR();
           } else
              lefFree(lefData->arrayName);
        } else
           lefFree(lefData->arrayName);
      } else
        lefFree(lefData->arrayName);
    }
#line 11931 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 823: /* array_rules: %empty  */
#line 6105 "lef/lef.y"
    { }
#line 11937 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 824: /* array_rules: array_rules array_rule  */
#line 6107 "lef/lef.y"
    { }
#line 11943 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 825: /* $@110: %empty  */
#line 6110 "lef/lef.y"
            { if (lefCallbacks->ArrayCbk) lefData->lefrDoSite = 1; lefData->lefDumbMode = 1; }
#line 11949 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 826: /* array_rule: site_word $@110 sitePattern ';'  */
#line 6112 "lef/lef.y"
    {
      if (lefCallbacks->ArrayCbk) {
        lefData->lefrArray.addSitePattern(lefData->lefrSitePatternPtr);
      }
    }
#line 11959 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 827: /* $@111: %empty  */
#line 6117 "lef/lef.y"
               {lefData->lefDumbMode = 1; if (lefCallbacks->ArrayCbk) lefData->lefrDoSite = 1; }
#line 11965 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 828: /* array_rule: K_CANPLACE $@111 sitePattern ';'  */
#line 6119 "lef/lef.y"
    {
      if (lefCallbacks->ArrayCbk) {
        lefData->lefrArray.addCanPlace(lefData->lefrSitePatternPtr);
      }
    }
#line 11975 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 829: /* $@112: %empty  */
#line 6124 "lef/lef.y"
                   {lefData->lefDumbMode = 1; if (lefCallbacks->ArrayCbk) lefData->lefrDoSite = 1; }
#line 11981 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 830: /* array_rule: K_CANNOTOCCUPY $@112 sitePattern ';'  */
#line 6126 "lef/lef.y"
    {
      if (lefCallbacks->ArrayCbk) {
        lefData->lefrArray.addCannotOccupy(lefData->lefrSitePatternPtr);
      }
    }
#line 11991 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 831: /* $@113: %empty  */
#line 6131 "lef/lef.y"
             { if (lefCallbacks->ArrayCbk) lefData->lefrDoTrack = 1; }
#line 11997 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 832: /* array_rule: K_TRACKS $@113 trackPattern ';'  */
#line 6132 "lef/lef.y"
    {
      if (lefCallbacks->ArrayCbk) {
        lefData->lefrArray.addTrack(lefData->lefrTrackPatternPtr);
      }
    }
#line 12007 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 833: /* array_rule: floorplan_start floorplan_list K_END T_STRING  */
#line 6138 "lef/lef.y"
    {
    }
#line 12014 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 834: /* $@114: %empty  */
#line 6140 "lef/lef.y"
                { if (lefCallbacks->ArrayCbk) lefData->lefrDoGcell = 1; }
#line 12020 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 835: /* array_rule: K_GCELLGRID $@114 gcellPattern ';'  */
#line 6141 "lef/lef.y"
    {
      if (lefCallbacks->ArrayCbk) {
        lefData->lefrArray.addGcell(lefData->lefrGcellPatternPtr);
      }
    }
#line 12030 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 836: /* array_rule: K_DEFAULTCAP int_number cap_list K_END K_DEFAULTCAP  */
#line 6147 "lef/lef.y"
    {
      if (lefCallbacks->ArrayCbk) {
        lefData->lefrArray.setTableSize((int)(yyvsp[-3].dval));
      }
    }
#line 12040 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 837: /* array_rule: def_statement  */
#line 6153 "lef/lef.y"
    { }
#line 12046 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 838: /* floorplan_start: K_FLOORPLAN T_STRING  */
#line 6156 "lef/lef.y"
    { if (lefCallbacks->ArrayCbk) lefData->lefrArray.addFloorPlan((yyvsp[0].string)); }
#line 12052 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 839: /* floorplan_list: %empty  */
#line 6160 "lef/lef.y"
    { }
#line 12058 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 840: /* floorplan_list: floorplan_list floorplan_element  */
#line 6162 "lef/lef.y"
    { }
#line 12064 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 841: /* $@115: %empty  */
#line 6165 "lef/lef.y"
             { lefData->lefDumbMode = 1; if (lefCallbacks->ArrayCbk) lefData->lefrDoSite = 1; }
#line 12070 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 842: /* floorplan_element: K_CANPLACE $@115 sitePattern ';'  */
#line 6167 "lef/lef.y"
    {
      if (lefCallbacks->ArrayCbk)
        lefData->lefrArray.addSiteToFloorPlan("CANPLACE",
        lefData->lefrSitePatternPtr);
    }
#line 12080 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 843: /* $@116: %empty  */
#line 6172 "lef/lef.y"
                   { if (lefCallbacks->ArrayCbk) lefData->lefrDoSite = 1; lefData->lefDumbMode = 1; }
#line 12086 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 844: /* floorplan_element: K_CANNOTOCCUPY $@116 sitePattern ';'  */
#line 6174 "lef/lef.y"
    {
      if (lefCallbacks->ArrayCbk)
        lefData->lefrArray.addSiteToFloorPlan("CANNOTOCCUPY",
        lefData->lefrSitePatternPtr);
     }
#line 12096 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 845: /* cap_list: %empty  */
#line 6182 "lef/lef.y"
    { }
#line 12102 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 846: /* cap_list: cap_list one_cap  */
#line 6184 "lef/lef.y"
    { }
#line 12108 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 847: /* one_cap: K_MINPINS int_number K_WIRECAP int_number ';'  */
#line 6187 "lef/lef.y"
    { if (lefCallbacks->ArrayCbk) lefData->lefrArray.addDefaultCap((int)(yyvsp[-3].dval), (yyvsp[-1].dval)); }
#line 12114 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 848: /* $@117: %empty  */
#line 6190 "lef/lef.y"
            {lefData->lefDumbMode=1;lefData->lefNlToken=TRUE;}
#line 12120 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 849: /* msg_statement: K_MESSAGE $@117 T_STRING '=' s_expr dtrm  */
#line 6191 "lef/lef.y"
    {  }
#line 12126 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 850: /* $@118: %empty  */
#line 6194 "lef/lef.y"
               {lefData->lefDumbMode=1;lefData->lefNlToken=TRUE;}
#line 12132 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 851: /* create_file_statement: K_CREATEFILE $@118 T_STRING '=' s_expr dtrm  */
#line 6195 "lef/lef.y"
    { }
#line 12138 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 852: /* $@119: %empty  */
#line 6198 "lef/lef.y"
           {lefData->lefDumbMode=1;lefData->lefNlToken=TRUE;}
#line 12144 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 853: /* def_statement: K_DEFINE $@119 T_STRING '=' expression dtrm  */
#line 6199 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6)
        lefAddNumDefine((yyvsp[-3].string), (yyvsp[-1].dval));
      else
        if (lefCallbacks->ArrayCbk) // write warning only if cbk is set 
           if (lefData->arrayWarnings++ < lefSettings->ArrayWarnings)
             lefWarning(2067, "DEFINE statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 12157 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 854: /* $@120: %empty  */
#line 6207 "lef/lef.y"
               {lefData->lefDumbMode=1;lefData->lefNlToken=TRUE;}
#line 12163 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 855: /* def_statement: K_DEFINES $@120 T_STRING '=' s_expr dtrm  */
#line 6208 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6)
        lefAddStringDefine((yyvsp[-3].string), (yyvsp[-1].string));
      else
        if (lefCallbacks->ArrayCbk) // write warning only if cbk is set 
           if (lefData->arrayWarnings++ < lefSettings->ArrayWarnings)
             lefWarning(2068, "DEFINES statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 12176 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 856: /* $@121: %empty  */
#line 6216 "lef/lef.y"
               {lefData->lefDumbMode=1;lefData->lefNlToken=TRUE;}
#line 12182 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 857: /* def_statement: K_DEFINEB $@121 T_STRING '=' b_expr dtrm  */
#line 6217 "lef/lef.y"
    {
      if (lefData->versionNum < 5.6)
        lefAddBooleanDefine((yyvsp[-3].string), (yyvsp[-1].integer));
      else
        if (lefCallbacks->ArrayCbk) // write warning only if cbk is set 
           if (lefData->arrayWarnings++ < lefSettings->ArrayWarnings)
             lefWarning(2069, "DEFINEB statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 12195 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 859: /* dtrm: ';'  */
#line 6228 "lef/lef.y"
         {lefData->lefNlToken = FALSE;}
#line 12201 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 860: /* dtrm: '\n'  */
#line 6229 "lef/lef.y"
                 {lefData->lefNlToken = FALSE;}
#line 12207 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 865: /* expression: expression '+' expression  */
#line 6242 "lef/lef.y"
                                {(yyval.dval) = (yyvsp[-2].dval) + (yyvsp[0].dval); }
#line 12213 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 866: /* expression: expression '-' expression  */
#line 6243 "lef/lef.y"
                                {(yyval.dval) = (yyvsp[-2].dval) - (yyvsp[0].dval); }
#line 12219 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 867: /* expression: expression '*' expression  */
#line 6244 "lef/lef.y"
                                {(yyval.dval) = (yyvsp[-2].dval) * (yyvsp[0].dval); }
#line 12225 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 868: /* expression: expression '/' expression  */
#line 6245 "lef/lef.y"
                                {(yyval.dval) = (yyvsp[-2].dval) / (yyvsp[0].dval); }
#line 12231 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 869: /* expression: '-' expression  */
#line 6246 "lef/lef.y"
                                {(yyval.dval) = -(yyvsp[0].dval);}
#line 12237 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 870: /* expression: '(' expression ')'  */
#line 6247 "lef/lef.y"
                                {(yyval.dval) = (yyvsp[-1].dval);}
#line 12243 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 871: /* expression: K_IF b_expr then expression else expression  */
#line 6249 "lef/lef.y"
                {(yyval.dval) = ((yyvsp[-4].integer) != 0) ? (yyvsp[-2].dval) : (yyvsp[0].dval);}
#line 12249 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 872: /* expression: int_number  */
#line 6250 "lef/lef.y"
                                     {(yyval.dval) = (yyvsp[0].dval);}
#line 12255 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 873: /* b_expr: expression relop expression  */
#line 6253 "lef/lef.y"
                              {(yyval.integer) = comp_num((yyvsp[-2].dval),(yyvsp[-1].integer),(yyvsp[0].dval));}
#line 12261 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 874: /* b_expr: expression K_AND expression  */
#line 6254 "lef/lef.y"
                                {(yyval.integer) = (yyvsp[-2].dval) != 0 && (yyvsp[0].dval) != 0;}
#line 12267 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 875: /* b_expr: expression K_OR expression  */
#line 6255 "lef/lef.y"
                                {(yyval.integer) = (yyvsp[-2].dval) != 0 || (yyvsp[0].dval) != 0;}
#line 12273 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 876: /* b_expr: s_expr relop s_expr  */
#line 6256 "lef/lef.y"
                              {(yyval.integer) = comp_str((yyvsp[-2].string),(yyvsp[-1].integer),(yyvsp[0].string));}
#line 12279 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 877: /* b_expr: s_expr K_AND s_expr  */
#line 6257 "lef/lef.y"
                              {(yyval.integer) = (yyvsp[-2].string)[0] != 0 && (yyvsp[0].string)[0] != 0;}
#line 12285 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 878: /* b_expr: s_expr K_OR s_expr  */
#line 6258 "lef/lef.y"
                              {(yyval.integer) = (yyvsp[-2].string)[0] != 0 || (yyvsp[0].string)[0] != 0;}
#line 12291 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 879: /* b_expr: b_expr K_EQ b_expr  */
#line 6259 "lef/lef.y"
                              {(yyval.integer) = (yyvsp[-2].integer) == (yyvsp[0].integer);}
#line 12297 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 880: /* b_expr: b_expr K_NE b_expr  */
#line 6260 "lef/lef.y"
                              {(yyval.integer) = (yyvsp[-2].integer) != (yyvsp[0].integer);}
#line 12303 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 881: /* b_expr: b_expr K_AND b_expr  */
#line 6261 "lef/lef.y"
                              {(yyval.integer) = (yyvsp[-2].integer) && (yyvsp[0].integer);}
#line 12309 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 882: /* b_expr: b_expr K_OR b_expr  */
#line 6262 "lef/lef.y"
                              {(yyval.integer) = (yyvsp[-2].integer) || (yyvsp[0].integer);}
#line 12315 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 883: /* b_expr: K_NOT b_expr  */
#line 6263 "lef/lef.y"
                                               {(yyval.integer) = !(yyval.integer);}
#line 12321 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 884: /* b_expr: '(' b_expr ')'  */
#line 6264 "lef/lef.y"
                              {(yyval.integer) = (yyvsp[-1].integer);}
#line 12327 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 885: /* b_expr: K_IF b_expr then b_expr else b_expr  */
#line 6266 "lef/lef.y"
        {(yyval.integer) = ((yyvsp[-4].integer) != 0) ? (yyvsp[-2].integer) : (yyvsp[0].integer);}
#line 12333 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 886: /* b_expr: K_TRUE  */
#line 6267 "lef/lef.y"
                              {(yyval.integer) = 1;}
#line 12339 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 887: /* b_expr: K_FALSE  */
#line 6268 "lef/lef.y"
                              {(yyval.integer) = 0;}
#line 12345 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 888: /* s_expr: s_expr '+' s_expr  */
#line 6272 "lef/lef.y"
    {
      (yyval.string) = (char*)lefMalloc(strlen((yyvsp[-2].string))+strlen((yyvsp[0].string))+1);
      strcpy((yyval.string),(yyvsp[-2].string));
      strcat((yyval.string),(yyvsp[0].string));
    }
#line 12355 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 889: /* s_expr: '(' s_expr ')'  */
#line 6278 "lef/lef.y"
    { (yyval.string) = (yyvsp[-1].string); }
#line 12361 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 890: /* s_expr: K_IF b_expr then s_expr else s_expr  */
#line 6280 "lef/lef.y"
    {
      lefData->lefDefIf = TRUE;
      if ((yyvsp[-4].integer) != 0) {
        (yyval.string) = (yyvsp[-2].string);        
      } else {
        (yyval.string) = (yyvsp[0].string);
      }
    }
#line 12374 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 891: /* s_expr: QSTRING  */
#line 6289 "lef/lef.y"
    { (yyval.string) = (yyvsp[0].string); }
#line 12380 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 892: /* relop: K_LE  */
#line 6292 "lef/lef.y"
       {(yyval.integer) = C_LE;}
#line 12386 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 893: /* relop: K_LT  */
#line 6293 "lef/lef.y"
         {(yyval.integer) = C_LT;}
#line 12392 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 894: /* relop: K_GE  */
#line 6294 "lef/lef.y"
         {(yyval.integer) = C_GE;}
#line 12398 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 895: /* relop: K_GT  */
#line 6295 "lef/lef.y"
         {(yyval.integer) = C_GT;}
#line 12404 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 896: /* relop: K_EQ  */
#line 6296 "lef/lef.y"
         {(yyval.integer) = C_EQ;}
#line 12410 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 897: /* relop: K_NE  */
#line 6297 "lef/lef.y"
         {(yyval.integer) = C_NE;}
#line 12416 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 898: /* relop: '='  */
#line 6298 "lef/lef.y"
         {(yyval.integer) = C_EQ;}
#line 12422 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 899: /* relop: '<'  */
#line 6299 "lef/lef.y"
         {(yyval.integer) = C_LT;}
#line 12428 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 900: /* relop: '>'  */
#line 6300 "lef/lef.y"
         {(yyval.integer) = C_GT;}
#line 12434 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 901: /* $@122: %empty  */
#line 6304 "lef/lef.y"
    { 
      if (lefCallbacks->PropBeginCbk)
        CALLBACK(lefCallbacks->PropBeginCbk, lefrPropBeginCbkType, 0);
    }
#line 12443 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 902: /* prop_def_section: K_PROPDEF $@122 prop_stmts K_END K_PROPDEF  */
#line 6309 "lef/lef.y"
    { 
      if (lefCallbacks->PropEndCbk)
        CALLBACK(lefCallbacks->PropEndCbk, lefrPropEndCbkType, 0);
    }
#line 12452 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 903: /* prop_stmts: %empty  */
#line 6316 "lef/lef.y"
    { }
#line 12458 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 904: /* prop_stmts: prop_stmts prop_stmt  */
#line 6318 "lef/lef.y"
    { }
#line 12464 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 905: /* $@123: %empty  */
#line 6321 "lef/lef.y"
            {lefData->lefDumbMode = 1; lefData->lefrProp.clear(); }
#line 12470 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 906: /* prop_stmt: K_LIBRARY $@123 T_STRING prop_define ';'  */
#line 6323 "lef/lef.y"
    { 
      if (lefCallbacks->PropCbk) {
        lefData->lefrProp.setPropType("library", (yyvsp[-2].string));
        CALLBACK(lefCallbacks->PropCbk, lefrPropCbkType, &lefData->lefrProp);
      }
      lefSettings->lefProps.lefrLibProp.setPropType((yyvsp[-2].string), lefData->lefPropDefType);
    }
#line 12482 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 907: /* $@124: %empty  */
#line 6330 "lef/lef.y"
                   {lefData->lefDumbMode = 1; lefData->lefrProp.clear(); }
#line 12488 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 908: /* prop_stmt: K_COMPONENTPIN $@124 T_STRING prop_define ';'  */
#line 6332 "lef/lef.y"
    { 
      if (lefCallbacks->PropCbk) {
        lefData->lefrProp.setPropType("componentpin", (yyvsp[-2].string));
        CALLBACK(lefCallbacks->PropCbk, lefrPropCbkType, &lefData->lefrProp);
      }
      lefSettings->lefProps.lefrCompProp.setPropType((yyvsp[-2].string), lefData->lefPropDefType);
    }
#line 12500 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 909: /* $@125: %empty  */
#line 6339 "lef/lef.y"
          {lefData->lefDumbMode = 1; lefData->lefrProp.clear(); }
#line 12506 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 910: /* prop_stmt: K_PIN $@125 T_STRING prop_define ';'  */
#line 6341 "lef/lef.y"
    { 
      if (lefCallbacks->PropCbk) {
        lefData->lefrProp.setPropType("pin", (yyvsp[-2].string));
        CALLBACK(lefCallbacks->PropCbk, lefrPropCbkType, &lefData->lefrProp);
      }
      lefSettings->lefProps.lefrPinProp.setPropType((yyvsp[-2].string), lefData->lefPropDefType);
      
    }
#line 12519 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 911: /* $@126: %empty  */
#line 6349 "lef/lef.y"
            {lefData->lefDumbMode = 1; lefData->lefrProp.clear(); }
#line 12525 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 912: /* prop_stmt: K_MACRO $@126 T_STRING prop_define ';'  */
#line 6351 "lef/lef.y"
    { 
      if (lefCallbacks->PropCbk) {
        lefData->lefrProp.setPropType("macro", (yyvsp[-2].string));
        CALLBACK(lefCallbacks->PropCbk, lefrPropCbkType, &lefData->lefrProp);
      }
      lefSettings->lefProps.lefrMacroProp.setPropType((yyvsp[-2].string), lefData->lefPropDefType);
    }
#line 12537 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 913: /* $@127: %empty  */
#line 6358 "lef/lef.y"
          {lefData->lefDumbMode = 1; lefData->lefrProp.clear(); }
#line 12543 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 914: /* prop_stmt: K_VIA $@127 T_STRING prop_define ';'  */
#line 6360 "lef/lef.y"
    { 
      if (lefCallbacks->PropCbk) {
        lefData->lefrProp.setPropType("via", (yyvsp[-2].string));
        CALLBACK(lefCallbacks->PropCbk, lefrPropCbkType, &lefData->lefrProp);
      }
      lefSettings->lefProps.lefrViaProp.setPropType((yyvsp[-2].string), lefData->lefPropDefType);
    }
#line 12555 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 915: /* $@128: %empty  */
#line 6367 "lef/lef.y"
              {lefData->lefDumbMode = 1; lefData->lefrProp.clear(); }
#line 12561 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 916: /* prop_stmt: K_VIARULE $@128 T_STRING prop_define ';'  */
#line 6369 "lef/lef.y"
    { 
      if (lefCallbacks->PropCbk) {
        lefData->lefrProp.setPropType("viarule", (yyvsp[-2].string));
        CALLBACK(lefCallbacks->PropCbk, lefrPropCbkType, &lefData->lefrProp);
      }
      lefSettings->lefProps.lefrViaRuleProp.setPropType((yyvsp[-2].string), lefData->lefPropDefType);
    }
#line 12573 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 917: /* $@129: %empty  */
#line 6376 "lef/lef.y"
            {lefData->lefDumbMode = 1; lefData->lefrProp.clear(); }
#line 12579 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 918: /* prop_stmt: K_LAYER $@129 T_STRING prop_define ';'  */
#line 6378 "lef/lef.y"
    { 
      if (lefCallbacks->PropCbk) {
        lefData->lefrProp.setPropType("layer", (yyvsp[-2].string));
        CALLBACK(lefCallbacks->PropCbk, lefrPropCbkType, &lefData->lefrProp);
      }
      lefSettings->lefProps.lefrLayerProp.setPropType((yyvsp[-2].string), lefData->lefPropDefType);
    }
#line 12591 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 919: /* $@130: %empty  */
#line 6385 "lef/lef.y"
                     {lefData->lefDumbMode = 1; lefData->lefrProp.clear(); }
#line 12597 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 920: /* prop_stmt: K_NONDEFAULTRULE $@130 T_STRING prop_define ';'  */
#line 6387 "lef/lef.y"
    { 
      if (lefCallbacks->PropCbk) {
        lefData->lefrProp.setPropType("nondefaultrule", (yyvsp[-2].string));
        CALLBACK(lefCallbacks->PropCbk, lefrPropCbkType, &lefData->lefrProp);
      }
      lefSettings->lefProps.lefrNondefProp.setPropType((yyvsp[-2].string), lefData->lefPropDefType);
    }
#line 12609 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 921: /* prop_define: K_INTEGER opt_def_range opt_def_dvalue  */
#line 6397 "lef/lef.y"
    { 
      if (lefCallbacks->PropCbk) lefData->lefrProp.setPropInteger();
      lefData->lefPropDefType = 'I';
    }
#line 12618 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 922: /* prop_define: K_REAL opt_def_range opt_def_value  */
#line 6402 "lef/lef.y"
    { 
      if (lefCallbacks->PropCbk) lefData->lefrProp.setPropReal();
      lefData->lefPropDefType = 'R';
    }
#line 12627 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 923: /* prop_define: K_STRING  */
#line 6407 "lef/lef.y"
    {
      if (lefCallbacks->PropCbk) lefData->lefrProp.setPropString();
      lefData->lefPropDefType = 'S';
    }
#line 12636 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 924: /* prop_define: K_STRING QSTRING  */
#line 6412 "lef/lef.y"
    {
      if (lefCallbacks->PropCbk) lefData->lefrProp.setPropQString((yyvsp[0].string));
      lefData->lefPropDefType = 'Q';
    }
#line 12645 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 925: /* prop_define: K_NAMEMAPSTRING T_STRING  */
#line 6417 "lef/lef.y"
    {
      if (lefCallbacks->PropCbk) lefData->lefrProp.setPropNameMapString((yyvsp[0].string));
      lefData->lefPropDefType = 'S';
    }
#line 12654 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 926: /* opt_range_second: %empty  */
#line 6424 "lef/lef.y"
    { }
#line 12660 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 927: /* opt_range_second: K_USELENGTHTHRESHOLD  */
#line 6426 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
        lefData->lefrLayer.setSpacingRangeUseLength();
    }
#line 12669 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 928: /* opt_range_second: K_INFLUENCE int_number  */
#line 6431 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
        lefData->lefrLayer.setSpacingRangeInfluence((yyvsp[0].dval));
        lefData->lefrLayer.setSpacingRangeInfluenceRange(-1, -1);
      }
    }
#line 12680 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 929: /* opt_range_second: K_INFLUENCE int_number K_RANGE int_number int_number  */
#line 6438 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
        lefData->lefrLayer.setSpacingRangeInfluence((yyvsp[-3].dval));
        lefData->lefrLayer.setSpacingRangeInfluenceRange((yyvsp[-1].dval), (yyvsp[0].dval));
      }
    }
#line 12691 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 930: /* opt_range_second: K_RANGE int_number int_number  */
#line 6445 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
        lefData->lefrLayer.setSpacingRangeRange((yyvsp[-1].dval), (yyvsp[0].dval));
    }
#line 12700 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 931: /* opt_endofline: %empty  */
#line 6452 "lef/lef.y"
    { }
#line 12706 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 932: /* $@131: %empty  */
#line 6454 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
        lefData->lefrLayer.setSpacingParSW((yyvsp[-2].dval), (yyvsp[0].dval));
    }
#line 12715 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 934: /* opt_endofline_twoedges: %empty  */
#line 6462 "lef/lef.y"
    { }
#line 12721 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 935: /* opt_endofline_twoedges: K_TWOEDGES  */
#line 6464 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
        lefData->lefrLayer.setSpacingParTwoEdges();
    }
#line 12730 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 936: /* opt_samenetPGonly: %empty  */
#line 6471 "lef/lef.y"
    { }
#line 12736 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 937: /* opt_samenetPGonly: K_PGONLY  */
#line 6473 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
        lefData->lefrLayer.setSpacingSamenetPGonly();
    }
#line 12745 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 938: /* opt_def_range: %empty  */
#line 6480 "lef/lef.y"
    { }
#line 12751 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 939: /* opt_def_range: K_RANGE int_number int_number  */
#line 6482 "lef/lef.y"
    {  if (lefCallbacks->PropCbk) lefData->lefrProp.setRange((yyvsp[-1].dval), (yyvsp[0].dval)); }
#line 12757 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 940: /* opt_def_value: %empty  */
#line 6486 "lef/lef.y"
    { }
#line 12763 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 941: /* opt_def_value: NUMBER  */
#line 6488 "lef/lef.y"
    { if (lefCallbacks->PropCbk) lefData->lefrProp.setNumber((yyvsp[0].dval)); }
#line 12769 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 942: /* opt_def_dvalue: %empty  */
#line 6492 "lef/lef.y"
    { }
#line 12775 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 943: /* opt_def_dvalue: int_number  */
#line 6494 "lef/lef.y"
    { if (lefCallbacks->PropCbk) lefData->lefrProp.setNumber((yyvsp[0].dval)); }
#line 12781 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 946: /* layer_spacing_opt: K_CENTERTOCENTER  */
#line 6501 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
         if (lefData->hasSpCenter) {
           if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1663, "A CENTERTOCENTER statement was already defined in SPACING\nCENTERTOCENTER can only be defined once per LAYER CUT SPACING.");
              CHKERR();
           }
        }
        lefData->hasSpCenter = 1;
        if (lefData->versionNum < 5.6) {
           if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "CENTERTOCENTER statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1664, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
        if (lefCallbacks->LayerCbk)
          lefData->lefrLayer.setSpacingCenterToCenter();
      }
    }
#line 12809 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 947: /* $@132: %empty  */
#line 6525 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
        if (lefData->hasSpSamenet) {
           if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefError(1665, "A SAMENET statement was already defined in SPACING\nSAMENET can only be defined once per LAYER CUT SPACING.");
              CHKERR();
           }
        }
        lefData->hasSpSamenet = 1;
        if (lefCallbacks->LayerCbk)
          lefData->lefrLayer.setSpacingSamenet();
       }
    }
#line 12827 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 948: /* layer_spacing_opt: K_SAMENET $@132 opt_samenetPGonly  */
#line 6539 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
        lefData->outMsg = (char*)lefMalloc(10000);
        sprintf(lefData->outMsg,
          "SAMENET is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
        lefError(1684, lefData->outMsg);
        lefFree(lefData->outMsg);
        CHKERR();
      }
    }
#line 12842 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 949: /* layer_spacing_opt: K_PARALLELOVERLAP  */
#line 6550 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
        lefData->outMsg = (char*)lefMalloc(10000);
        sprintf(lefData->outMsg,
          "PARALLELOVERLAP is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
        lefError(1680, lefData->outMsg);
        lefFree(lefData->outMsg);
        CHKERR(); 
      } else {
        if (lefCallbacks->LayerCbk) {
          if (lefData->hasSpParallel) {
             if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                lefError(1666, "A PARALLELOVERLAP statement was already defined in SPACING\nPARALLELOVERLAP can only be defined once per LAYER CUT SPACING.");
                CHKERR();
             }
          }
          lefData->hasSpParallel = 1;
          if (lefCallbacks->LayerCbk)
            lefData->lefrLayer.setSpacingParallelOverlap();
        }
      }
    }
#line 12869 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 951: /* $@133: %empty  */
#line 6575 "lef/lef.y"
            {lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 12875 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 952: /* $@134: %empty  */
#line 6576 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
{
        if (lefData->versionNum < 5.7) {
           if (lefData->hasSpSamenet) {    // 5.6 and earlier does not allow 
              if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                 lefError(1667, "A SAMENET statement was already defined in SPACING\nEither SAMENET or LAYER can be defined, but not both.");
                 CHKERR();
              }
           }
        }
        lefData->lefrLayer.setSpacingName((yyvsp[0].string));
      }
    }
#line 12894 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 954: /* $@135: %empty  */
#line 6592 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
        if (lefData->versionNum < 5.5) {
           if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
              lefData->outMsg = (char*)lefMalloc(10000);
              sprintf (lefData->outMsg,
                 "ADJACENTCUTS statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
              lefError(1668, lefData->outMsg);
              lefFree(lefData->outMsg);
              CHKERR();
           }
        }
        if (lefData->versionNum < 5.7) {
           if (lefData->hasSpSamenet) {    // 5.6 and earlier does not allow 
              if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                 lefError(1669, "A SAMENET statement was already defined in SPACING\nEither SAMENET or ADJACENTCUTS can be defined, but not both.");
                 CHKERR();
              }
           }
        }
        lefData->lefrLayer.setSpacingAdjacent((int)(yyvsp[-2].dval), (yyvsp[0].dval));
      }
    }
#line 12922 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 956: /* layer_spacing_cut_routing: K_AREA NUMBER  */
#line 6617 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
        lefData->outMsg = (char*)lefMalloc(10000);
        sprintf(lefData->outMsg,
          "AREA is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
        lefError(1693, lefData->outMsg);
        lefFree(lefData->outMsg);
        CHKERR();
      } else {
        if (lefCallbacks->LayerCbk) {
          if (lefData->versionNum < 5.7) {
             if (lefData->hasSpSamenet) {    // 5.6 and earlier does not allow 
                if (lefData->layerWarnings++ < lefSettings->LayerWarnings) {
                   lefError(1670, "A SAMENET statement was already defined in SPACING\nEither SAMENET or AREA can be defined, but not both.");
                   CHKERR();
                }
             }
          }
          lefData->lefrLayer.setSpacingArea((yyvsp[0].dval));
        }
      }
    }
#line 12949 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 957: /* $@136: %empty  */
#line 6640 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
        lefData->lefrLayer.setSpacingRange((yyvsp[-1].dval), (yyvsp[0].dval));
    }
#line 12958 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 959: /* layer_spacing_cut_routing: K_LENGTHTHRESHOLD int_number  */
#line 6646 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
        lefData->lefrLayer.setSpacingLength((yyvsp[0].dval));
      }
    }
#line 12968 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 960: /* layer_spacing_cut_routing: K_LENGTHTHRESHOLD int_number K_RANGE int_number int_number  */
#line 6652 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk) {
        lefData->lefrLayer.setSpacingLength((yyvsp[-3].dval));
        lefData->lefrLayer.setSpacingLengthRange((yyvsp[-1].dval), (yyvsp[0].dval));
      }
    }
#line 12979 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 961: /* $@137: %empty  */
#line 6659 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
        lefData->lefrLayer.setSpacingEol((yyvsp[-2].dval), (yyvsp[0].dval));
    }
#line 12988 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 962: /* layer_spacing_cut_routing: K_ENDOFLINE int_number K_WITHIN int_number $@137 opt_endofline  */
#line 6664 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
        lefData->outMsg = (char*)lefMalloc(10000);
        sprintf(lefData->outMsg,
          "ENDOFLINE is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
        lefError(1681, lefData->outMsg);
        lefFree(lefData->outMsg);
        CHKERR();
      }
    }
#line 13003 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 963: /* layer_spacing_cut_routing: K_NOTCHLENGTH int_number  */
#line 6675 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
        lefData->outMsg = (char*)lefMalloc(10000);
        sprintf(lefData->outMsg,
          "NOTCHLENGTH is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
        lefError(1682, lefData->outMsg);
        lefFree(lefData->outMsg);
        CHKERR();
      } else {
        if (lefCallbacks->LayerCbk)
          lefData->lefrLayer.setSpacingNotchLength((yyvsp[0].dval));
      }
    }
#line 13021 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 964: /* layer_spacing_cut_routing: K_ENDOFNOTCHWIDTH int_number K_NOTCHSPACING int_number K_NOTCHLENGTH int_number  */
#line 6689 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
        lefData->outMsg = (char*)lefMalloc(10000);
        sprintf(lefData->outMsg,
          "ENDOFNOTCHWIDTH is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
        lefError(1696, lefData->outMsg);
        lefFree(lefData->outMsg);
        CHKERR();
      } else {
        if (lefCallbacks->LayerCbk)
          lefData->lefrLayer.setSpacingEndOfNotchWidth((yyvsp[-4].dval), (yyvsp[-2].dval), (yyvsp[0].dval));
      }
    }
#line 13039 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 965: /* spacing_cut_layer_opt: %empty  */
#line 6705 "lef/lef.y"
    {}
#line 13045 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 966: /* spacing_cut_layer_opt: K_STACK  */
#line 6707 "lef/lef.y"
    {
      if (lefCallbacks->LayerCbk)
        lefData->lefrLayer.setSpacingLayerStack();
    }
#line 13054 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 967: /* opt_adjacentcuts_exceptsame: %empty  */
#line 6714 "lef/lef.y"
    {}
#line 13060 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 968: /* opt_adjacentcuts_exceptsame: K_EXCEPTSAMEPGNET  */
#line 6716 "lef/lef.y"
    {
      if (lefData->versionNum < 5.7) {
        lefData->outMsg = (char*)lefMalloc(10000);
        sprintf(lefData->outMsg,
          "EXCEPTSAMEPGNET is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", lefData->versionNum);
        lefError(1683, lefData->outMsg);
        lefFree(lefData->outMsg);
        CHKERR();
      } else {
        if (lefCallbacks->LayerCbk)
          lefData->lefrLayer.setSpacingAdjacentExcept();
      }
    }
#line 13078 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 969: /* opt_layer_name: %empty  */
#line 6732 "lef/lef.y"
    { (yyval.string) = 0; }
#line 13084 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 970: /* $@138: %empty  */
#line 6733 "lef/lef.y"
            {lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 13090 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 971: /* opt_layer_name: K_LAYER $@138 T_STRING  */
#line 6734 "lef/lef.y"
    { (yyval.string) = (yyvsp[0].string); }
#line 13096 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 972: /* $@139: %empty  */
#line 6738 "lef/lef.y"
           {lefData->lefDumbMode = 1; lefData->lefNoNum = 1; }
#line 13102 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 973: /* req_layer_name: K_LAYER $@139 T_STRING  */
#line 6739 "lef/lef.y"
    { (yyval.string) = (yyvsp[0].string); }
#line 13108 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 974: /* universalnoisemargin: K_UNIVERSALNOISEMARGIN int_number int_number ';'  */
#line 6743 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->NoiseMarginCbk) {
          lefData->lefrNoiseMargin.low = (yyvsp[-2].dval);
          lefData->lefrNoiseMargin.high = (yyvsp[-1].dval);
          CALLBACK(lefCallbacks->NoiseMarginCbk, lefrNoiseMarginCbkType, &lefData->lefrNoiseMargin);
        }
      } else
        if (lefCallbacks->NoiseMarginCbk) // write warning only if cbk is set 
          if (lefData->noiseMarginWarnings++ < lefSettings->NoiseMarginWarnings)
            lefWarning(2070, "UNIVERSALNOISEMARGIN statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 13125 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 975: /* edgeratethreshold1: K_EDGERATETHRESHOLD1 int_number ';'  */
#line 6757 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->EdgeRateThreshold1Cbk) {
          CALLBACK(lefCallbacks->EdgeRateThreshold1Cbk,
          lefrEdgeRateThreshold1CbkType, (yyvsp[-1].dval));
        }
      } else
        if (lefCallbacks->EdgeRateThreshold1Cbk) // write warning only if cbk is set 
          if (lefData->edgeRateThreshold1Warnings++ < lefSettings->EdgeRateThreshold1Warnings)
            lefWarning(2071, "EDGERATETHRESHOLD1 statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 13141 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 976: /* edgeratethreshold2: K_EDGERATETHRESHOLD2 int_number ';'  */
#line 6770 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->EdgeRateThreshold2Cbk) {
          CALLBACK(lefCallbacks->EdgeRateThreshold2Cbk,
          lefrEdgeRateThreshold2CbkType, (yyvsp[-1].dval));
        }
      } else
        if (lefCallbacks->EdgeRateThreshold2Cbk) // write warning only if cbk is set 
          if (lefData->edgeRateThreshold2Warnings++ < lefSettings->EdgeRateThreshold2Warnings)
            lefWarning(2072, "EDGERATETHRESHOLD2 statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 13157 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 977: /* edgeratescalefactor: K_EDGERATESCALEFACTOR int_number ';'  */
#line 6783 "lef/lef.y"
    {
      if (lefData->versionNum < 5.4) {
        if (lefCallbacks->EdgeRateScaleFactorCbk) {
          CALLBACK(lefCallbacks->EdgeRateScaleFactorCbk,
          lefrEdgeRateScaleFactorCbkType, (yyvsp[-1].dval));
        }
      } else
        if (lefCallbacks->EdgeRateScaleFactorCbk) // write warning only if cbk is set 
          if (lefData->edgeRateScaleFactorWarnings++ < lefSettings->EdgeRateScaleFactorWarnings)
            lefWarning(2073, "EDGERATESCALEFACTOR statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 13173 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 978: /* $@140: %empty  */
#line 6796 "lef/lef.y"
    { if (lefCallbacks->NoiseTableCbk) lefData->lefrNoiseTable.setup((int)(yyvsp[0].dval)); }
#line 13179 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 979: /* noisetable: K_NOISETABLE int_number $@140 ';' noise_table_list end_noisetable dtrm  */
#line 6798 "lef/lef.y"
    { }
#line 13185 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 980: /* end_noisetable: K_END K_NOISETABLE  */
#line 6802 "lef/lef.y"
  {
    if (lefData->versionNum < 5.4) {
      if (lefCallbacks->NoiseTableCbk)
        CALLBACK(lefCallbacks->NoiseTableCbk, lefrNoiseTableCbkType, &lefData->lefrNoiseTable);
    } else
      if (lefCallbacks->NoiseTableCbk) // write warning only if cbk is set 
        if (lefData->noiseTableWarnings++ < lefSettings->NoiseTableWarnings)
          lefWarning(2074, "NOISETABLE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
  }
#line 13199 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 983: /* noise_table_entry: K_EDGERATE int_number ';'  */
#line 6820 "lef/lef.y"
    { if (lefCallbacks->NoiseTableCbk)
         {
            lefData->lefrNoiseTable.newEdge();
            lefData->lefrNoiseTable.addEdge((yyvsp[-1].dval));
         }
    }
#line 13210 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 984: /* noise_table_entry: output_resistance_entry  */
#line 6827 "lef/lef.y"
    { }
#line 13216 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 985: /* $@141: %empty  */
#line 6830 "lef/lef.y"
    { if (lefCallbacks->NoiseTableCbk) lefData->lefrNoiseTable.addResistance(); }
#line 13222 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 987: /* num_list: int_number  */
#line 6836 "lef/lef.y"
    { if (lefCallbacks->NoiseTableCbk)
    lefData->lefrNoiseTable.addResistanceNumber((yyvsp[0].dval)); }
#line 13229 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 988: /* num_list: num_list int_number  */
#line 6839 "lef/lef.y"
    { if (lefCallbacks->NoiseTableCbk)
    lefData->lefrNoiseTable.addResistanceNumber((yyvsp[0].dval)); }
#line 13236 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 991: /* $@142: %empty  */
#line 6848 "lef/lef.y"
        { if (lefCallbacks->NoiseTableCbk)
        lefData->lefrNoiseTable.addVictimLength((yyvsp[-1].dval)); }
#line 13243 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 992: /* victim: K_VICTIMLENGTH int_number ';' $@142 K_VICTIMNOISE vnoiselist ';'  */
#line 6851 "lef/lef.y"
        { }
#line 13249 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 993: /* vnoiselist: int_number  */
#line 6855 "lef/lef.y"
    { if (lefCallbacks->NoiseTableCbk)
    lefData->lefrNoiseTable.addVictimNoise((yyvsp[0].dval)); }
#line 13256 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 994: /* vnoiselist: vnoiselist int_number  */
#line 6858 "lef/lef.y"
    { if (lefCallbacks->NoiseTableCbk)
    lefData->lefrNoiseTable.addVictimNoise((yyvsp[0].dval)); }
#line 13263 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 995: /* $@143: %empty  */
#line 6862 "lef/lef.y"
    { if (lefCallbacks->CorrectionTableCbk)
    lefData->lefrCorrectionTable.setup((int)(yyvsp[-1].dval)); }
#line 13270 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 996: /* correctiontable: K_CORRECTIONTABLE int_number ';' $@143 correction_table_list end_correctiontable dtrm  */
#line 6865 "lef/lef.y"
    { }
#line 13276 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 997: /* end_correctiontable: K_END K_CORRECTIONTABLE  */
#line 6869 "lef/lef.y"
  {
    if (lefData->versionNum < 5.4) {
      if (lefCallbacks->CorrectionTableCbk)
        CALLBACK(lefCallbacks->CorrectionTableCbk, lefrCorrectionTableCbkType,
               &lefData->lefrCorrectionTable);
    } else
      if (lefCallbacks->CorrectionTableCbk) // write warning only if cbk is set 
        if (lefData->correctionTableWarnings++ < lefSettings->CorrectionTableWarnings)
          lefWarning(2075, "CORRECTIONTABLE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
  }
#line 13291 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1000: /* correction_table_item: K_EDGERATE int_number ';'  */
#line 6887 "lef/lef.y"
    { if (lefCallbacks->CorrectionTableCbk)
         {
            lefData->lefrCorrectionTable.newEdge();
            lefData->lefrCorrectionTable.addEdge((yyvsp[-1].dval));
         }
    }
#line 13302 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1001: /* correction_table_item: output_list  */
#line 6894 "lef/lef.y"
    { }
#line 13308 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1002: /* $@144: %empty  */
#line 6897 "lef/lef.y"
  { if (lefCallbacks->CorrectionTableCbk)
  lefData->lefrCorrectionTable.addResistance(); }
#line 13315 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1003: /* output_list: K_OUTPUTRESISTANCE $@144 numo_list ';' corr_victim_list  */
#line 6900 "lef/lef.y"
  { }
#line 13321 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1004: /* numo_list: int_number  */
#line 6904 "lef/lef.y"
    { if (lefCallbacks->CorrectionTableCbk)
    lefData->lefrCorrectionTable.addResistanceNumber((yyvsp[0].dval)); }
#line 13328 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1005: /* numo_list: numo_list int_number  */
#line 6907 "lef/lef.y"
    { if (lefCallbacks->CorrectionTableCbk)
    lefData->lefrCorrectionTable.addResistanceNumber((yyvsp[0].dval)); }
#line 13335 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1008: /* $@145: %empty  */
#line 6917 "lef/lef.y"
     { if (lefCallbacks->CorrectionTableCbk)
     lefData->lefrCorrectionTable.addVictimLength((yyvsp[-1].dval)); }
#line 13342 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1009: /* corr_victim: K_VICTIMLENGTH int_number ';' $@145 K_CORRECTIONFACTOR corr_list ';'  */
#line 6920 "lef/lef.y"
     { }
#line 13348 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1010: /* corr_list: int_number  */
#line 6924 "lef/lef.y"
    { if (lefCallbacks->CorrectionTableCbk)
        lefData->lefrCorrectionTable.addVictimCorrection((yyvsp[0].dval)); }
#line 13355 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1011: /* corr_list: corr_list int_number  */
#line 6927 "lef/lef.y"
    { if (lefCallbacks->CorrectionTableCbk)
        lefData->lefrCorrectionTable.addVictimCorrection((yyvsp[0].dval)); }
#line 13362 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1012: /* input_antenna: K_INPUTPINANTENNASIZE int_number ';'  */
#line 6933 "lef/lef.y"
    { // 5.3 syntax 
        lefData->use5_3 = 1;
        if (lefData->ignoreVersion) {
           // do nothing 
        } else if (lefData->versionNum > 5.3) {
           // A 5.3 syntax in 5.4 
           if (lefData->use5_4) {
              if (lefCallbacks->InputAntennaCbk) { // write warning only if cbk is set 
                if (lefData->inputAntennaWarnings++ < lefSettings->InputAntennaWarnings) {
                   lefData->outMsg = (char*)lefMalloc(10000);
                   sprintf (lefData->outMsg,
                      "INPUTPINANTENNASIZE statement is a version 5.3 or earlier syntax.\nYour lef file with version %g, has both old and new INPUTPINANTENNASIZE syntax, which is incorrect.", lefData->versionNum);
                   lefError(1671, lefData->outMsg);
                   lefFree(lefData->outMsg);
                   CHKERR();
                }
              }
           }
        }
        if (lefCallbacks->InputAntennaCbk)
          CALLBACK(lefCallbacks->InputAntennaCbk, lefrInputAntennaCbkType, (yyvsp[-1].dval));
    }
#line 13389 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1013: /* output_antenna: K_OUTPUTPINANTENNASIZE int_number ';'  */
#line 6957 "lef/lef.y"
    { // 5.3 syntax 
        lefData->use5_3 = 1;
        if (lefData->ignoreVersion) {
           // do nothing 
        } else if (lefData->versionNum > 5.3) {
           // A 5.3 syntax in 5.4 
           if (lefData->use5_4) {
              if (lefCallbacks->OutputAntennaCbk) { // write warning only if cbk is set 
                if (lefData->outputAntennaWarnings++ < lefSettings->OutputAntennaWarnings) {
                   lefData->outMsg = (char*)lefMalloc(10000);
                   sprintf (lefData->outMsg,
                      "OUTPUTPINANTENNASIZE statement is a version 5.3 or earlier syntax.\nYour lef file with version %g, has both old and new OUTPUTPINANTENNASIZE syntax, which is incorrect.", lefData->versionNum);
                   lefError(1672, lefData->outMsg);
                   lefFree(lefData->outMsg);
                   CHKERR();
                }
              }
           }
        }
        if (lefCallbacks->OutputAntennaCbk)
          CALLBACK(lefCallbacks->OutputAntennaCbk, lefrOutputAntennaCbkType, (yyvsp[-1].dval));
    }
#line 13416 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1014: /* inout_antenna: K_INOUTPINANTENNASIZE int_number ';'  */
#line 6981 "lef/lef.y"
    { // 5.3 syntax 
        lefData->use5_3 = 1;
        if (lefData->ignoreVersion) {
           // do nothing 
        } else if (lefData->versionNum > 5.3) {
           // A 5.3 syntax in 5.4 
           if (lefData->use5_4) {
              if (lefCallbacks->InoutAntennaCbk) { // write warning only if cbk is set 
                if (lefData->inoutAntennaWarnings++ < lefSettings->InoutAntennaWarnings) {
                   lefData->outMsg = (char*)lefMalloc(10000);
                   sprintf (lefData->outMsg,
                      "INOUTPINANTENNASIZE statement is a version 5.3 or earlier syntax.\nYour lef file with version %g, has both old and new INOUTPINANTENNASIZE syntax, which is incorrect.", lefData->versionNum);
                   lefError(1673, lefData->outMsg);
                   lefFree(lefData->outMsg);
                   CHKERR();
                }
              }
           }
        }
        if (lefCallbacks->InoutAntennaCbk)
          CALLBACK(lefCallbacks->InoutAntennaCbk, lefrInoutAntennaCbkType, (yyvsp[-1].dval));
    }
#line 13443 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1015: /* antenna_input: K_ANTENNAINPUTGATEAREA NUMBER ';'  */
#line 7005 "lef/lef.y"
    { // 5.4 syntax 
        // 11/12/2002 - this is obsolete in 5.5, suppose should be ingored 
        // 12/16/2002 - talked to Dave Noice, leave them in here for debugging
        lefData->use5_4 = 1;
        if (lefData->ignoreVersion) {
           // do nothing 
        } else if (lefData->versionNum < 5.4) {
           if (lefCallbacks->AntennaInputCbk) { // write warning only if cbk is set 
             if (lefData->antennaInputWarnings++ < lefSettings->AntennaInputWarnings) {
               lefData->outMsg = (char*)lefMalloc(10000);
               sprintf (lefData->outMsg,
                  "ANTENNAINPUTGATEAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.\nEither update your VERSION number or use the 5.3 syntax.", lefData->versionNum);
               lefError(1674, lefData->outMsg);
               lefFree(lefData->outMsg);
               CHKERR();
             }
           }
        } else if (lefData->use5_3) {
           if (lefCallbacks->AntennaInputCbk) { // write warning only if cbk is set 
             if (lefData->antennaInputWarnings++ < lefSettings->AntennaInputWarnings) {
                lefData->outMsg = (char*)lefMalloc(10000);
                sprintf (lefData->outMsg,
                   "ANTENNAINPUTGATEAREA statement is a version 5.4 or later syntax.\nYour lef file with version %g, has both old and new ANTENNAINPUTGATEAREA syntax, which is incorrect.", lefData->versionNum);
                lefError(1675, lefData->outMsg);
                lefFree(lefData->outMsg);
               CHKERR();
             }
           }
        }
        if (lefCallbacks->AntennaInputCbk)
          CALLBACK(lefCallbacks->AntennaInputCbk, lefrAntennaInputCbkType, (yyvsp[-1].dval));
    }
#line 13480 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1016: /* antenna_inout: K_ANTENNAINOUTDIFFAREA NUMBER ';'  */
#line 7039 "lef/lef.y"
    { // 5.4 syntax 
        // 11/12/2002 - this is obsolete in 5.5, & will be ignored 
        // 12/16/2002 - talked to Dave Noice, leave them in here for debugging
        lefData->use5_4 = 1;
        if (lefData->ignoreVersion) {
           // do nothing 
        } else if (lefData->versionNum < 5.4) {
           if (lefCallbacks->AntennaInoutCbk) { // write warning only if cbk is set 
              if (lefData->antennaInoutWarnings++ < lefSettings->AntennaInoutWarnings) {
                 lefData->outMsg = (char*)lefMalloc(10000);
                 sprintf (lefData->outMsg,
                    "ANTENNAINOUTDIFFAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.\nEither update your VERSION number or use the 5.3 syntax.", lefData->versionNum);
                 lefError(1676, lefData->outMsg);
                 lefFree(lefData->outMsg);
                 CHKERR();
              }
           }
        } else if (lefData->use5_3) {
           if (lefCallbacks->AntennaInoutCbk) { // write warning only if cbk is set 
              if (lefData->antennaInoutWarnings++ < lefSettings->AntennaInoutWarnings) {
                 lefData->outMsg = (char*)lefMalloc(10000);
                 sprintf (lefData->outMsg,
                    "ANTENNAINOUTDIFFAREA statement is a version 5.4 or later syntax.\nYour lef file with version %g, has both old and new ANTENNAINOUTDIFFAREA syntax, which is incorrect.", lefData->versionNum);
                 lefError(1677, lefData->outMsg);
                 lefFree(lefData->outMsg);
                 CHKERR();
              }
           }
        }
        if (lefCallbacks->AntennaInoutCbk)
          CALLBACK(lefCallbacks->AntennaInoutCbk, lefrAntennaInoutCbkType, (yyvsp[-1].dval));
    }
#line 13517 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1017: /* antenna_output: K_ANTENNAOUTPUTDIFFAREA NUMBER ';'  */
#line 7073 "lef/lef.y"
    { // 5.4 syntax 
        // 11/12/2002 - this is obsolete in 5.5, & will be ignored 
        // 12/16/2002 - talked to Dave Noice, leave them in here for debugging
        lefData->use5_4 = 1;
        if (lefData->ignoreVersion) {
           // do nothing 
        } else if (lefData->versionNum < 5.4) {
           if (lefCallbacks->AntennaOutputCbk) { // write warning only if cbk is set 
              if (lefData->antennaOutputWarnings++ < lefSettings->AntennaOutputWarnings) {
                 lefData->outMsg = (char*)lefMalloc(10000);
                 sprintf (lefData->outMsg,
                    "ANTENNAOUTPUTDIFFAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.\nEither update your VERSION number or use the 5.3 syntax.", lefData->versionNum);
                 lefError(1678, lefData->outMsg);
                 lefFree(lefData->outMsg);
                 CHKERR();
              }
           }
        } else if (lefData->use5_3) {
           if (lefCallbacks->AntennaOutputCbk) { // write warning only if cbk is set 
              if (lefData->antennaOutputWarnings++ < lefSettings->AntennaOutputWarnings) {
                 lefData->outMsg = (char*)lefMalloc(10000);
                 sprintf (lefData->outMsg,
                    "ANTENNAOUTPUTDIFFAREA statement is a version 5.4 or later syntax.\nYour lef file with version %g, has both old and new ANTENNAOUTPUTDIFFAREA syntax, which is incorrect.", lefData->versionNum);
                 lefError(1679, lefData->outMsg);
                 lefFree(lefData->outMsg);
                 CHKERR();
              }
           }
        }
        if (lefCallbacks->AntennaOutputCbk)
          CALLBACK(lefCallbacks->AntennaOutputCbk, lefrAntennaOutputCbkType, (yyvsp[-1].dval));
    }
#line 13554 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;

  case 1020: /* extension: K_BEGINEXT  */
#line 7110 "lef/lef.y"
    { 
        if (lefCallbacks->ExtensionCbk)
          CALLBACK(lefCallbacks->ExtensionCbk, lefrExtensionCbkType, &lefData->Hist_text[0]);
        if (lefData->versionNum >= 5.6)
           lefData->ge56almostDone = 1;
    }
#line 13565 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"
    break;


#line 13569 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/thirdparty/lefdef/5.8/lef/lef.tab.cpp"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", YY_CAST (yysymbol_kind_t, yyr1[yyn]), &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYSYMBOL_YYEMPTY : YYTRANSLATE (yychar);
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
      yyerror (YY_("syntax error"));
    }

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;
  ++yynerrs;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  /* Pop stack until we find a state that shifts the error token.  */
  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYSYMBOL_YYerror;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYSYMBOL_YYerror)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  YY_ACCESSING_SYMBOL (yystate), yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", YY_ACCESSING_SYMBOL (yyn), yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturnlab;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturnlab;


/*-----------------------------------------------------------.
| yyexhaustedlab -- YYNOMEM (memory exhaustion) comes here.  |
`-----------------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturnlab;


/*----------------------------------------------------------.
| yyreturnlab -- parsing is finished, clean up and return.  |
`----------------------------------------------------------*/
yyreturnlab:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 7117 "lef/lef.y"


END_LEFDEF_PARSER_NAMESPACE
