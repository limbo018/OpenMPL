// A Bison parser, made by GNU Bison 3.8.2.

// Skeleton interface for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015, 2018-2021 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// As a special exception, you may create a larger work that contains
// part or all of the Bison parser skeleton and distribute that work
// under terms of your choice, so long as that work isn't itself a
// parser generator using the skeleton or a modified version thereof
// as a parser skeleton.  Alternatively, if you modify or redistribute
// the parser skeleton itself, you may (at your option) remove this
// special exception, which will cause the skeleton and the resulting
// Bison output files to be licensed under the GNU General Public
// License without this special exception.

// This special exception was added by the Free Software Foundation in
// version 2.2 of Bison.


/**
 ** \file /home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.h
 ** Define the LefParser::parser class.
 */

// C++ LALR(1) parser skeleton written by Akim Demaille.

// DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
// especially those whose name start with YY_ or yy_.  They are
// private implementation details that can be changed or removed.

#ifndef YY_LEFPARSER_HOME_FANTASTIC_NA_OPENMPL_OPENMPL_BUILD_THIRDPARTY_LIMBO_LIMBO_PARSERS_LEF_BISON_LEFPARSER_H_INCLUDED
# define YY_LEFPARSER_HOME_FANTASTIC_NA_OPENMPL_OPENMPL_BUILD_THIRDPARTY_LIMBO_LIMBO_PARSERS_LEF_BISON_LEFPARSER_H_INCLUDED


# include <cstdlib> // std::abort
# include <iostream>
# include <stdexcept>
# include <string>
# include <vector>

#if defined __cplusplus
# define YY_CPLUSPLUS __cplusplus
#else
# define YY_CPLUSPLUS 199711L
#endif

// Support move semantics when possible.
#if 201103L <= YY_CPLUSPLUS
# define YY_MOVE           std::move
# define YY_MOVE_OR_COPY   move
# define YY_MOVE_REF(Type) Type&&
# define YY_RVREF(Type)    Type&&
# define YY_COPY(Type)     Type
#else
# define YY_MOVE
# define YY_MOVE_OR_COPY   copy
# define YY_MOVE_REF(Type) Type&
# define YY_RVREF(Type)    const Type&
# define YY_COPY(Type)     const Type&
#endif

// Support noexcept when possible.
#if 201103L <= YY_CPLUSPLUS
# define YY_NOEXCEPT noexcept
# define YY_NOTHROW
#else
# define YY_NOEXCEPT
# define YY_NOTHROW throw ()
#endif

// Support constexpr when possible.
#if 201703 <= YY_CPLUSPLUS
# define YY_CONSTEXPR constexpr
#else
# define YY_CONSTEXPR
#endif
# include "location.hh"


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

/* Debug traces.  */
#ifndef LEFPARSERDEBUG
# if defined YYDEBUG
#if YYDEBUG
#   define LEFPARSERDEBUG 1
#  else
#   define LEFPARSERDEBUG 0
#  endif
# else /* ! defined YYDEBUG */
#  define LEFPARSERDEBUG 1
# endif /* ! defined YYDEBUG */
#endif  /* ! defined LEFPARSERDEBUG */

namespace LefParser {
#line 190 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.h"




  /// A Bison parser.
  class Parser
  {
  public:
#ifdef LEFPARSERSTYPE
# ifdef __GNUC__
#  pragma GCC message "bison: do not #define LEFPARSERSTYPE in C++, use %define api.value.type"
# endif
    typedef LEFPARSERSTYPE value_type;
#else
    /// Symbol semantic values.
    union value_type
    {
#line 62 "LefParser.yy"

    int  			integerVal;
    double 			doubleVal;
	/*char*           string;*/
    std::string*		stringVal;
	std::string*		qstringVal;
	std::string*		binaryVal;

	struct lefPoint* pt; 

#line 219 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.h"

    };
#endif
    /// Backward compatibility (Bison 3.8).
    typedef value_type semantic_type;

    /// Symbol locations.
    typedef location location_type;

    /// Syntax errors thrown from user actions.
    struct syntax_error : std::runtime_error
    {
      syntax_error (const location_type& l, const std::string& m)
        : std::runtime_error (m)
        , location (l)
      {}

      syntax_error (const syntax_error& s)
        : std::runtime_error (s.what ())
        , location (s.location)
      {}

      ~syntax_error () YY_NOEXCEPT YY_NOTHROW;

      location_type location;
    };

    /// Token kinds.
    struct token
    {
      enum token_kind_type
      {
        LEFPARSEREMPTY = -2,
    END = 0,                       // "end of file"
    LEFPARSERerror = 256,          // error
    LEFPARSERUNDEF = 257,          // "invalid token"
    K_DEFINE = 258,                // "DEFINE"
    K_DEFINEB = 259,               // "DEFINEB"
    K_DEFINES = 260,               // "DEFINES"
    K_MESSAGE = 261,               // "MESSAGE"
    K_CREATEFILE = 262,            // "CREATEFILE"
    K_OPENFILE = 263,              // "OPENFILE"
    K_CLOSEFILE = 264,             // "CLOSEFILE"
    K_WARNING = 265,               // "WARNING"
    K_ERROR = 266,                 // "ERROR"
    K_FATALERROR = 267,            // "FATALERROR"
    K_ABOVE = 268,                 // "ABOVE"
    K_ABUT = 269,                  // "ABUT"
    K_ABUTMENT = 270,              // "ABUTMENT"
    K_ACCURRENTDENSITY = 271,      // "ACCURRENTDENSITY"
    K_ACTIVE = 272,                // "ACTIVE"
    K_ADJACENTCUTS = 273,          // "ADJACENTCUTS"
    K_ANALOG = 274,                // "ANALOG"
    K_AND = 275,                   // "AND"
    K_ANTENNAAREAFACTOR = 276,     // "ANTENNAAREAFACTOR"
    K_ANTENNAAREADIFFREDUCEPWL = 277, // "ANTENNAAREADIFFREDUCEPWL"
    K_ANTENNAAREAMINUSDIFF = 278,  // "ANTENNAAREAMINUSDIFF"
    K_ANTENNAAREARATIO = 279,      // "ANTENNAAREARATIO"
    K_ANTENNACELL = 280,           // "ANTENNACELL"
    K_ANTENNACUMAREARATIO = 281,   // "ANTENNACUMAREARATIO"
    K_ANTENNACUMDIFFAREARATIO = 282, // "ANTENNACUMDIFFAREARATIO"
    K_ANTENNACUMDIFFSIDEAREARATIO = 283, // "ANTENNACUMDIFFSIDEAREARATIO"
    K_ANTENNACUMROUTINGPLUSCUT = 284, // "ANTENNACUMROUTINGPLUSCUT"
    K_ANTENNACUMSIDEAREARATIO = 285, // "ANTENNACUMSIDEAREARATIO"
    K_ANTENNADIFFAREA = 286,       // "ANTENNADIFFAREA"
    K_ANTENNADIFFAREARATIO = 287,  // "ANTENNADIFFAREARATIO"
    K_ANTENNADIFFSIDEAREARATIO = 288, // "ANTENNADIFFSIDEAREARATIO"
    K_ANTENNAGATEAREA = 289,       // "ANTENNAGATEAREA"
    K_ANTENNAGATEPLUSDIFF = 290,   // "ANTENNAGATEPLUSDIFF"
    K_ANTENNAINOUTDIFFAREA = 291,  // "ANTENNAINOUTDIFFAREA"
    K_ANTENNAINPUTGATEAREA = 292,  // "ANTENNAINPUTGATEAREA"
    K_ANTENNALENGTHFACTOR = 293,   // "ANTENNALENGTHFACTOR"
    K_ANTENNAMAXAREACAR = 294,     // "ANTENNAMAXAREACAR"
    K_ANTENNAMAXCUTCAR = 295,      // "ANTENNAMAXCUTCAR"
    K_ANTENNAMAXSIDEAREACAR = 296, // "ANTENNAMAXSIDEAREACAR"
    K_ANTENNAMETALAREA = 297,      // "ANTENNAMETALAREA"
    K_ANTENNAMETALLENGTH = 298,    // "ANTENNAMETALLENGTH"
    K_ANTENNAMODEL = 299,          // "ANTENNAMODEL"
    K_ANTENNAOUTPUTDIFFAREA = 300, // "ANTENNAOUTPUTDIFFAREA"
    K_ANTENNAPARTIALCUTAREA = 301, // "ANTENNAPARTIALCUTAREA"
    K_ANTENNAPARTIALMETALAREA = 302, // "ANTENNAPARTIALMETALAREA"
    K_ANTENNAPARTIALMETALSIDEAREA = 303, // "ANTENNAPARTIALMETALSIDEAREA"
    K_ANTENNASIDEAREARATIO = 304,  // "ANTENNASIDEAREARATIO"
    K_ANTENNASIZE = 305,           // "ANTENNASIZE"
    K_ANTENNASIDEAREAFACTOR = 306, // "ANTENNASIDEAREAFACTOR"
    K_ANYEDGE = 307,               // "ANYEDGE"
    K_AREA = 308,                  // "AREA"
    K_AREAIO = 309,                // "AREAIO"
    K_ARRAY = 310,                 // "ARRAY"
    K_ARRAYCUTS = 311,             // "ARRAYCUTS"
    K_ARRAYSPACING = 312,          // "ARRAYSPACING"
    K_AVERAGE = 313,               // "AVERAGE"
    K_BELOW = 314,                 // "BELOW"
    K_BEGINEXT = 315,              // "BEGINEXT"
    K_BLACKBOX = 316,              // "BLACKBOX"
    K_BLOCK = 317,                 // "BLOCK"
    K_BOTTOMLEFT = 318,            // "BOTTOMLEFT"
    K_BOTTOMRIGHT = 319,           // "BOTTOMRIGHT"
    K_BUMP = 320,                  // "BUMP"
    K_BUSBITCHARS = 321,           // "BUSBITCHARS"
    K_BUFFER = 322,                // "BUFFER"
    K_BY = 323,                    // "BY"
    K_CANNOTOCCUPY = 324,          // "CANNOTOCCUPY"
    K_CANPLACE = 325,              // "CANPLACE"
    K_CAPACITANCE = 326,           // "CAPACITANCE"
    K_CAPMULTIPLIER = 327,         // "CAPMULTIPLIER"
    K_CENTERTOCENTER = 328,        // "CENTERTOCENTER"
    K_CLASS = 329,                 // "CLASS"
    K_CLEARANCEMEASURE = 330,      // "CLEARANCEMEASURE"
    K_CLOCK = 331,                 // "CLOCK"
    K_CLOCKTYPE = 332,             // "CLOCKTYPE"
    K_COLUMNMAJOR = 333,           // "COLUMNMAJOR"
    K_CURRENTDEN = 334,            // "CURRENTDEN"
    K_COMPONENTPIN = 335,          // "COMPONENTPIN"
    K_CORE = 336,                  // "CORE"
    K_CORNER = 337,                // "CORNER"
    K_CORRECTIONFACTOR = 338,      // "CORRECTIONFACTOR"
    K_CORRECTIONTABLE = 339,       // "CORRECTIONTABLE"
    K_COVER = 340,                 // "COVER"
    K_CPERSQDIST = 341,            // "CPERSQDIST"
    K_CURRENT = 342,               // "CURRENT"
    K_CURRENTSOURCE = 343,         // "CURRENTSOURCE"
    K_CUT = 344,                   // "CUT"
    K_CUTAREA = 345,               // "CUTAREA"
    K_CUTSIZE = 346,               // "CUTSIZE"
    K_CUTSPACING = 347,            // "CUTSPACING"
    K_DATA = 348,                  // "DATA"
    K_DATABASE = 349,              // "DATABASE"
    K_DCCURRENTDENSITY = 350,      // "DCCURRENTDENSITY"
    K_DEFAULT = 351,               // "DEFAULT"
    K_DEFAULTCAP = 352,            // "DEFAULTCAP"
    K_DELAY = 353,                 // "DELAY"
    K_DENSITY = 354,               // "DENSITY"
    K_DENSITYCHECKSTEP = 355,      // "DENSITYCHECKSTEP"
    K_DENSITYCHECKWINDOW = 356,    // "DENSITYCHECKWINDOW"
    K_DESIGNRULEWIDTH = 357,       // "DESIGNRULEWIDTH"
    K_DIAG45 = 358,                // "DIAG45"
    K_DIAG135 = 359,               // "DIAG135"
    K_DIAGMINEDGELENGTH = 360,     // "DIAGMINEDGELENGTH"
    K_DIAGSPACING = 361,           // "DIAGSPACING"
    K_DIAGPITCH = 362,             // "DIAGPITCH"
    K_DIAGWIDTH = 363,             // "DIAGWIDTH"
    K_DIELECTRIC = 364,            // "DIELECTRIC"
    K_DIFFUSEONLY = 365,           // "DIFFUSEONLY"
    K_DIRECTION = 366,             // "DIRECTION"
    K_DIVIDERCHAR = 367,           // "DIVIDERCHAR"
    K_DO = 368,                    // "DO"
    K_E = 369,                     // "E"
    K_EDGECAPACITANCE = 370,       // "EDGECAPACITANCE"
    K_EDGERATE = 371,              // "EDGERATE"
    K_EDGERATESCALEFACTOR = 372,   // "EDGERATESCALEFACTOR"
    K_EDGERATETHRESHOLD1 = 373,    // "EDGERATETHRESHOLD1"
    K_EDGERATETHRESHOLD2 = 374,    // "EDGERATETHRESHOLD2"
    K_EEQ = 375,                   // "EEQ"
    K_ELSE = 376,                  // "ELSE"
    K_ENCLOSURE = 377,             // "ENCLOSURE"
    K_END = 378,                   // "END"
    K_ENDEXT = 379,                // "ENDEXT"
    K_ENDCAP = 380,                // "ENDCAP"
    K_ENDOFLINE = 381,             // "ENDOFLINE"
    K_ENDOFNOTCHWIDTH = 382,       // "ENDOFNOTCHWIDTH"
    K_EUCLIDEAN = 383,             // "EUCLIDEAN"
    K_EXCEPTEXTRACUT = 384,        // "EXCEPTEXTRACUT"
    K_EXCEPTSAMEPGNET = 385,       // "EXCEPTSAMEPGNET"
    K_EXCEPTPGNET = 386,           // "EXCEPTPGNET"
    K_EXTENSION = 387,             // "EXTENSION"
    K_FALL = 388,                  // "FALL"
    K_FALLCS = 389,                // "FALLCS"
    K_FALLRS = 390,                // "FALLRS"
    K_FALLSATCUR = 391,            // "FALLSATCUR"
    K_FALLSATT1 = 392,             // "FALLSATT1"
    K_FALLSLEWLIMIT = 393,         // "FALLSLEWLIMIT"
    K_FALLT0 = 394,                // "FALLT0"
    K_FALLTHRESH = 395,            // "FALLTHRESH"
    K_FALLVOLTAGETHRESHOLD = 396,  // "FALLVOLTAGETHRESHOLD"
    K_FALSE = 397,                 // "FALSE"
    K_FE = 398,                    // "FE"
    K_FEEDTHRU = 399,              // "FEEDTHRU"
    K_FILLACTIVESPACING = 400,     // "FILLACTIVESPACING"
    K_FIXED = 401,                 // "FIXED"
    K_FLIP = 402,                  // "FLIP"
    K_FLOORPLAN = 403,             // "FLOORPLAN"
    K_FN = 404,                    // "FN"
    K_FOREIGN = 405,               // "FOREIGN"
    K_FREQUENCY = 406,             // "FREQUENCY"
    K_FROMABOVE = 407,             // "FROMABOVE"
    K_FROMBELOW = 408,             // "FROMBELOW"
    K_FROMPIN = 409,               // "FROMPIN"
    K_FUNCTION = 410,              // "FUNCTION"
    K_FS = 411,                    // "FS"
    K_FW = 412,                    // "FW"
    K_GCELLGRID = 413,             // "GCELLGRID"
    K_GENERATE = 414,              // "GENERATE"
    K_GENERATED = 415,             // "GENERATED"
    K_GENERATOR = 416,             // "GENERATOR"
    K_GROUND = 417,                // "GROUND"
    K_GROUNDSENSITIVITY = 418,     // "GROUNDSENSITIVITY"
    K_HARDSPACING = 419,           // "HARDSPACING"
    K_HEIGHT = 420,                // "HEIGHT"
    K_HISTORY = 421,               // "HISTORY"
    K_HOLD = 422,                  // "HOLD"
    K_HORIZONTAL = 423,            // "HORIZONTAL"
    K_IF = 424,                    // "IF"
    K_IMPLANT = 425,               // "IMPLANT"
    K_INFLUENCE = 426,             // "INFLUENCE"
    K_INOUT = 427,                 // "INOUT"
    K_INOUTPINANTENNASIZE = 428,   // "INOUTPINANTENNASIZE"
    K_INPUT = 429,                 // "INPUT"
    K_INPUTPINANTENNASIZE = 430,   // "INPUTPINANTENNASIZE"
    K_INPUTNOISEMARGIN = 431,      // "INPUTNOISEMARGIN"
    K_INSIDECORNER = 432,          // "INSIDECORNER"
    K_INTEGER = 433,               // "INTEGER"
    K_INTRINSIC = 434,             // "INTRINSIC"
    K_INVERT = 435,                // "INVERT"
    K_INVERTER = 436,              // "INVERTER"
    K_IRDROP = 437,                // "IRDROP"
    K_ITERATE = 438,               // "ITERATE"
    K_IV_TABLES = 439,             // "IV_TABLES"
    K_LAYER = 440,                 // "LAYER"
    K_LAYERS = 441,                // "LAYERS"
    K_LEAKAGE = 442,               // "LEAKAGE"
    K_LENGTH = 443,                // "LENGTH"
    K_LENGTHSUM = 444,             // "LENGTHSUM"
    K_LENGTHTHRESHOLD = 445,       // "LENGTHTHRESHOLD"
    K_LEQ = 446,                   // "LEQ"
    K_LIBRARY = 447,               // "LIBRARY"
    K_LONGARRAY = 448,             // "LONGARRAY"
    K_MACRO = 449,                 // "MACRO"
    K_MANUFACTURINGGRID = 450,     // "MANUFACTURINGGRID"
    K_MASTERSLICE = 451,           // "MASTERSLICE"
    K_MATCH = 452,                 // "MATCH"
    K_MAXADJACENTSLOTSPACING = 453, // "MAXADJACENTSLOTSPACING"
    K_MAXCOAXIALSLOTSPACING = 454, // "MAXCOAXIALSLOTSPACING"
    K_MAXDELAY = 455,              // "MAXDELAY"
    K_MAXEDGES = 456,              // "MAXEDGES"
    K_MAXEDGESLOTSPACING = 457,    // "MAXEDGESLOTSPACING"
    K_MAXLOAD = 458,               // "MAXLOAD"
    K_MAXIMUMDENSITY = 459,        // "MAXIMUMDENSITY"
    K_MAXVIASTACK = 460,           // "MAXVIASTACK"
    K_MAXWIDTH = 461,              // "MAXWIDTH"
    K_MAXXY = 462,                 // "MAXXY"
    K_MEGAHERTZ = 463,             // "MEGAHERTZ"
    K_METALOVERHANG = 464,         // "METALOVERHANG"
    K_MICRONS = 465,               // "MICRONS"
    K_MILLIAMPS = 466,             // "MILLIAMPS"
    K_MILLIWATTS = 467,            // "MILLIWATTS"
    K_MINCUTS = 468,               // "MINCUTS"
    K_MINENCLOSEDAREA = 469,       // "MINENCLOSEDAREA"
    K_MINFEATURE = 470,            // "MINFEATURE"
    K_MINIMUMCUT = 471,            // "MINIMUMCUT"
    K_MINIMUMDENSITY = 472,        // "MINIMUMDENSITY"
    K_MINPINS = 473,               // "MINPINS"
    K_MINSIZE = 474,               // "MINSIZE"
    K_MINSTEP = 475,               // "MINSTEP"
    K_MINWIDTH = 476,              // "MINWIDTH"
    K_MPWH = 477,                  // "MPWH"
    K_MPWL = 478,                  // "MPWL"
    K_MUSTJOIN = 479,              // "MUSTJOIN"
    K_MX = 480,                    // "MX"
    K_MY = 481,                    // "MY"
    K_MXR90 = 482,                 // "MXR90"
    K_MYR90 = 483,                 // "MYR90"
    K_N = 484,                     // "N"
    K_NAMEMAPSTRING = 485,         // "NAMEMAPSTRING"
    K_NAMESCASESENSITIVE = 486,    // "NAMESCASESENSITIVE"
    K_NANOSECONDS = 487,           // "NANOSECONDS"
    K_NEGEDGE = 488,               // "NEGEDGE"
    K_NETEXPR = 489,               // "NETEXPR"
    K_NETS = 490,                  // "NETS"
    K_NEW = 491,                   // "NEW"
    K_NONDEFAULTRULE = 492,        // "NONDEFAULTRULE"
    K_NONE = 493,                  // "NONE"
    K_NONINVERT = 494,             // "NONINVERT"
    K_NONUNATE = 495,              // "NONUNATE"
    K_NOISETABLE = 496,            // "NOISETABLE"
    K_NOTCHLENGTH = 497,           // "NOTCHLENGTH"
    K_NOTCHSPACING = 498,          // "NOTCHSPACING"
    K_NOWIREEXTENSIONATPIN = 499,  // "NOWIREEXTENSIONATPIN"
    K_OBS = 500,                   // "OBS"
    K_OFF = 501,                   // "OFF"
    K_OFFSET = 502,                // "OFFSET"
    K_OHMS = 503,                  // "OHMS"
    K_ON = 504,                    // "ON"
    K_OR = 505,                    // "OR"
    K_ORIENT = 506,                // "ORIENT"
    K_ORIENTATION = 507,           // "ORIENTATION"
    K_ORIGIN = 508,                // "ORIGIN"
    K_ORTHOGONAL = 509,            // "ORTHOGONAL"
    K_OUTPUT = 510,                // "OUTPUT"
    K_OUTPUTPINANTENNASIZE = 511,  // "OUTPUTPINANTENNASIZE"
    K_OUTPUTNOISEMARGIN = 512,     // "OUTPUTNOISEMARGIN"
    K_OUTPUTRESISTANCE = 513,      // "OUTPUTRESISTANCE"
    K_OUTSIDECORNER = 514,         // "OUTSIDECORNER"
    K_OVERHANG = 515,              // "OVERHANG"
    K_OVERLAP = 516,               // "OVERLAP"
    K_OVERLAPS = 517,              // "OVERLAPS"
    K_OXIDE1 = 518,                // "OXIDE1"
    K_OXIDE2 = 519,                // "OXIDE2"
    K_OXIDE3 = 520,                // "OXIDE3"
    K_OXIDE4 = 521,                // "OXIDE4"
    K_PAD = 522,                   // "PAD"
    K_PARALLELEDGE = 523,          // "PARALLELEDGE"
    K_PARALLELOVERLAP = 524,       // "PARALLELOVERLAP"
    K_PARALLELRUNLENGTH = 525,     // "PARALLELRUNLENGTH"
    K_PATH = 526,                  // "PATH"
    K_PATTERN = 527,               // "PATTERN"
    K_PEAK = 528,                  // "PEAK"
    K_PERIOD = 529,                // "PERIOD"
    K_PGONLY = 530,                // "PGONLY"
    K_PICOFARADS = 531,            // "PICOFARADS"
    K_PIN = 532,                   // "PIN"
    K_PITCH = 533,                 // "PITCH"
    K_PLACED = 534,                // "PLACED"
    K_POLYGON = 535,               // "POLYGON"
    K_PORT = 536,                  // "PORT"
    K_POSEDGE = 537,               // "POSEDGE"
    K_POST = 538,                  // "POST"
    K_POWER = 539,                 // "POWER"
    K_PRE = 540,                   // "PRE"
    K_PREFERENCLOSURE = 541,       // "PREFERENCLOSURE"
    K_PRL = 542,                   // "PRL"
    K_PROPERTY = 543,              // "PROPERTY"
    K_PROPERTYDEFINITIONS = 544,   // "PROPERTYDEFINITIONS"
    K_PROTRUSIONWIDTH = 545,       // "PROTRUSIONWIDTH"
    K_PULLDOWNRES = 546,           // "PULLDOWNRES"
    K_PWL = 547,                   // "PWL"
    K_R0 = 548,                    // "R0"
    K_R90 = 549,                   // "R90"
    K_R180 = 550,                  // "R180"
    K_R270 = 551,                  // "R270"
    K_RANGE = 552,                 // "RANGE"
    K_REAL = 553,                  // "REAL"
    K_RECOVERY = 554,              // "RECOVERY"
    K_RECT = 555,                  // "RECT"
    K_RESISTANCE = 556,            // "RESISTANCE"
    K_RESISTIVE = 557,             // "RESISTIVE"
    K_RING = 558,                  // "RING"
    K_RISE = 559,                  // "RISE"
    K_RISECS = 560,                // "RISECS"
    K_RISERS = 561,                // "RISERS"
    K_RISESATCUR = 562,            // "RISESATCUR"
    K_RISESATT1 = 563,             // "RISESATT1"
    K_RISESLEWLIMIT = 564,         // "RISESLEWLIMIT"
    K_RISET0 = 565,                // "RISET0"
    K_RISETHRESH = 566,            // "RISETHRESH"
    K_RISEVOLTAGETHRESHOLD = 567,  // "RISEVOLTAGETHRESHOLD"
    K_RMS = 568,                   // "RMS"
    K_ROUTING = 569,               // "ROUTING"
    K_ROWABUTSPACING = 570,        // "ROWABUTSPACING"
    K_ROWCOL = 571,                // "ROWCOL"
    K_ROWMAJOR = 572,              // "ROWMAJOR"
    K_ROWMINSPACING = 573,         // "ROWMINSPACING"
    K_ROWPATTERN = 574,            // "ROWPATTERN"
    K_RPERSQ = 575,                // "RPERSQ"
    K_S = 576,                     // "S"
    K_SAMENET = 577,               // "SAMENET"
    K_SCANUSE = 578,               // "SCANUSE"
    K_SDFCOND = 579,               // "SDFCOND"
    K_SDFCONDEND = 580,            // "SDFCONDEND"
    K_SDFCONDSTART = 581,          // "SDFCONDSTART"
    K_SETUP = 582,                 // "SETUP"
    K_SHAPE = 583,                 // "SHAPE"
    K_SHRINKAGE = 584,             // "SHRINKAGE"
    K_SIGNAL = 585,                // "SIGNAL"
    K_SITE = 586,                  // "SITE"
    K_SIZE = 587,                  // "SIZE"
    K_SKEW = 588,                  // "SKEW"
    K_SLOTLENGTH = 589,            // "SLOTLENGTH"
    K_SLOTWIDTH = 590,             // "SLOTWIDTH"
    K_SLOTWIRELENGTH = 591,        // "SLOTWIRELENGTH"
    K_SLOTWIREWIDTH = 592,         // "SLOTWIREWIDTH"
    K_SPLITWIREWIDTH = 593,        // "SPLITWIREWIDTH"
    K_SOFT = 594,                  // "SOFT"
    K_SOURCE = 595,                // "SOURCE"
    K_SPACER = 596,                // "SPACER"
    K_SPACING = 597,               // "SPACING"
    K_SPACINGTABLE = 598,          // "SPACINGTABLE"
    K_SPECIALNETS = 599,           // "SPECIALNETS"
    K_STABLE = 600,                // "STABLE"
    K_STACK = 601,                 // "STACK"
    K_START = 602,                 // "START"
    K_STEP = 603,                  // "STEP"
    K_STOP = 604,                  // "STOP"
    K_STRING = 605,                // "STRING"
    K_STRUCTURE = 606,             // "STRUCTURE"
    K_SUPPLYSENSITIVITY = 607,     // "SUPPLYSENSITIVITY"
    K_SYMMETRY = 608,              // "SYMMETRY"
    K_TABLE = 609,                 // "TABLE"
    K_TABLEAXIS = 610,             // "TABLEAXIS"
    K_TABLEDIMENSION = 611,        // "TABLEDIMENSION"
    K_TABLEENTRIES = 612,          // "TABLEENTRIES"
    K_TAPERRULE = 613,             // "TAPERRULE"
    K_THEN = 614,                  // "THEN"
    K_THICKNESS = 615,             // "THICKNESS"
    K_TIEHIGH = 616,               // "TIEHIGH"
    K_TIELOW = 617,                // "TIELOW"
    K_TIEOFFR = 618,               // "TIEOFFR"
    K_TIME = 619,                  // "TIME"
    K_TIMING = 620,                // "TIMING"
    K_TO = 621,                    // "TO"
    K_TOPIN = 622,                 // "TOPIN"
    K_TOPLEFT = 623,               // "TOPLEFT"
    K_TOPOFSTACKONLY = 624,        // "TOPOFSTACKONLY"
    K_TOPRIGHT = 625,              // "TOPRIGHT"
    K_TRACKS = 626,                // "TRACKS"
    K_TRANSITIONTIME = 627,        // "TRANSITIONTIME"
    K_TRISTATE = 628,              // "TRISTATE"
    K_TRUE = 629,                  // "TRUE"
    K_TWOEDGES = 630,              // "TWOEDGES"
    K_TWOWIDTHS = 631,             // "TWOWIDTHS"
    K_TYPE = 632,                  // "TYPE"
    K_UNATENESS = 633,             // "UNATENESS"
    K_UNITS = 634,                 // "UNITS"
    K_UNIVERSALNOISEMARGIN = 635,  // "UNIVERSALNOISEMARGIN"
    K_USE = 636,                   // "USE"
    K_USELENGTHTHRESHOLD = 637,    // "USELENGTHTHRESHOLD"
    K_USEMINSPACING = 638,         // "USEMINSPACING"
    K_USER = 639,                  // "USER"
    K_USEVIA = 640,                // "USEVIA"
    K_USEVIARULE = 641,            // "USEVIARULE"
    K_VARIABLE = 642,              // "VARIABLE"
    K_VERSION = 643,               // "VERSION"
    K_VERTICAL = 644,              // "VERTICAL"
    K_VHI = 645,                   // "VHI"
    K_VIA = 646,                   // "VIA"
    K_VIARULE = 647,               // "VIARULE"
    K_VICTIMLENGTH = 648,          // "VICTIMLENGTH"
    K_VICTIMNOISE = 649,           // "VICTIMNOISE"
    K_VIRTUAL = 650,               // "VIRTUAL"
    K_VLO = 651,                   // "VLO"
    K_VOLTAGE = 652,               // "VOLTAGE"
    K_VOLTS = 653,                 // "VOLTS"
    K_W = 654,                     // "W"
    K_WELLTAP = 655,               // "WELLTAP"
    K_WIDTH = 656,                 // "WIDTH"
    K_WITHIN = 657,                // "WITHIN"
    K_WIRECAP = 658,               // "WIRECAP"
    K_WIREEXTENSION = 659,         // "WIREEXTENSION"
    K_X = 660,                     // "X"
    K_Y = 661,                     // "Y"
    K_EQ = 662,                    // K_EQ
    K_NE = 663,                    // K_NE
    K_LE = 664,                    // K_LE
    K_LT = 665,                    // K_LT
    K_GE = 666,                    // K_GE
    K_GT = 667,                    // K_GT
    K_NOT = 668,                   // K_NOT
    INTEGER = 669,                 // "integer"
    DOUBLE = 670,                  // "double"
    STRING = 671,                  // "string"
    QSTRING = 672,                 // "qstring"
    BINARY = 673,                  // "binary numbers"
    IF = 676,                      // IF
    LNOT = 677,                    // LNOT
    UMINUS = 678                   // UMINUS
      };
      /// Backward compatibility alias (Bison 3.6).
      typedef token_kind_type yytokentype;
    };

    /// Token kind, as returned by yylex.
    typedef token::token_kind_type token_kind_type;

    /// Backward compatibility alias (Bison 3.6).
    typedef token_kind_type token_type;

    /// Symbol kinds.
    struct symbol_kind
    {
      enum symbol_kind_type
      {
        YYNTOKENS = 435, ///< Number of tokens.
        S_YYEMPTY = -2,
        S_YYEOF = 0,                             // "end of file"
        S_YYerror = 1,                           // error
        S_YYUNDEF = 2,                           // "invalid token"
        S_K_DEFINE = 3,                          // "DEFINE"
        S_K_DEFINEB = 4,                         // "DEFINEB"
        S_K_DEFINES = 5,                         // "DEFINES"
        S_K_MESSAGE = 6,                         // "MESSAGE"
        S_K_CREATEFILE = 7,                      // "CREATEFILE"
        S_K_OPENFILE = 8,                        // "OPENFILE"
        S_K_CLOSEFILE = 9,                       // "CLOSEFILE"
        S_K_WARNING = 10,                        // "WARNING"
        S_K_ERROR = 11,                          // "ERROR"
        S_K_FATALERROR = 12,                     // "FATALERROR"
        S_K_ABOVE = 13,                          // "ABOVE"
        S_K_ABUT = 14,                           // "ABUT"
        S_K_ABUTMENT = 15,                       // "ABUTMENT"
        S_K_ACCURRENTDENSITY = 16,               // "ACCURRENTDENSITY"
        S_K_ACTIVE = 17,                         // "ACTIVE"
        S_K_ADJACENTCUTS = 18,                   // "ADJACENTCUTS"
        S_K_ANALOG = 19,                         // "ANALOG"
        S_K_AND = 20,                            // "AND"
        S_K_ANTENNAAREAFACTOR = 21,              // "ANTENNAAREAFACTOR"
        S_K_ANTENNAAREADIFFREDUCEPWL = 22,       // "ANTENNAAREADIFFREDUCEPWL"
        S_K_ANTENNAAREAMINUSDIFF = 23,           // "ANTENNAAREAMINUSDIFF"
        S_K_ANTENNAAREARATIO = 24,               // "ANTENNAAREARATIO"
        S_K_ANTENNACELL = 25,                    // "ANTENNACELL"
        S_K_ANTENNACUMAREARATIO = 26,            // "ANTENNACUMAREARATIO"
        S_K_ANTENNACUMDIFFAREARATIO = 27,        // "ANTENNACUMDIFFAREARATIO"
        S_K_ANTENNACUMDIFFSIDEAREARATIO = 28,    // "ANTENNACUMDIFFSIDEAREARATIO"
        S_K_ANTENNACUMROUTINGPLUSCUT = 29,       // "ANTENNACUMROUTINGPLUSCUT"
        S_K_ANTENNACUMSIDEAREARATIO = 30,        // "ANTENNACUMSIDEAREARATIO"
        S_K_ANTENNADIFFAREA = 31,                // "ANTENNADIFFAREA"
        S_K_ANTENNADIFFAREARATIO = 32,           // "ANTENNADIFFAREARATIO"
        S_K_ANTENNADIFFSIDEAREARATIO = 33,       // "ANTENNADIFFSIDEAREARATIO"
        S_K_ANTENNAGATEAREA = 34,                // "ANTENNAGATEAREA"
        S_K_ANTENNAGATEPLUSDIFF = 35,            // "ANTENNAGATEPLUSDIFF"
        S_K_ANTENNAINOUTDIFFAREA = 36,           // "ANTENNAINOUTDIFFAREA"
        S_K_ANTENNAINPUTGATEAREA = 37,           // "ANTENNAINPUTGATEAREA"
        S_K_ANTENNALENGTHFACTOR = 38,            // "ANTENNALENGTHFACTOR"
        S_K_ANTENNAMAXAREACAR = 39,              // "ANTENNAMAXAREACAR"
        S_K_ANTENNAMAXCUTCAR = 40,               // "ANTENNAMAXCUTCAR"
        S_K_ANTENNAMAXSIDEAREACAR = 41,          // "ANTENNAMAXSIDEAREACAR"
        S_K_ANTENNAMETALAREA = 42,               // "ANTENNAMETALAREA"
        S_K_ANTENNAMETALLENGTH = 43,             // "ANTENNAMETALLENGTH"
        S_K_ANTENNAMODEL = 44,                   // "ANTENNAMODEL"
        S_K_ANTENNAOUTPUTDIFFAREA = 45,          // "ANTENNAOUTPUTDIFFAREA"
        S_K_ANTENNAPARTIALCUTAREA = 46,          // "ANTENNAPARTIALCUTAREA"
        S_K_ANTENNAPARTIALMETALAREA = 47,        // "ANTENNAPARTIALMETALAREA"
        S_K_ANTENNAPARTIALMETALSIDEAREA = 48,    // "ANTENNAPARTIALMETALSIDEAREA"
        S_K_ANTENNASIDEAREARATIO = 49,           // "ANTENNASIDEAREARATIO"
        S_K_ANTENNASIZE = 50,                    // "ANTENNASIZE"
        S_K_ANTENNASIDEAREAFACTOR = 51,          // "ANTENNASIDEAREAFACTOR"
        S_K_ANYEDGE = 52,                        // "ANYEDGE"
        S_K_AREA = 53,                           // "AREA"
        S_K_AREAIO = 54,                         // "AREAIO"
        S_K_ARRAY = 55,                          // "ARRAY"
        S_K_ARRAYCUTS = 56,                      // "ARRAYCUTS"
        S_K_ARRAYSPACING = 57,                   // "ARRAYSPACING"
        S_K_AVERAGE = 58,                        // "AVERAGE"
        S_K_BELOW = 59,                          // "BELOW"
        S_K_BEGINEXT = 60,                       // "BEGINEXT"
        S_K_BLACKBOX = 61,                       // "BLACKBOX"
        S_K_BLOCK = 62,                          // "BLOCK"
        S_K_BOTTOMLEFT = 63,                     // "BOTTOMLEFT"
        S_K_BOTTOMRIGHT = 64,                    // "BOTTOMRIGHT"
        S_K_BUMP = 65,                           // "BUMP"
        S_K_BUSBITCHARS = 66,                    // "BUSBITCHARS"
        S_K_BUFFER = 67,                         // "BUFFER"
        S_K_BY = 68,                             // "BY"
        S_K_CANNOTOCCUPY = 69,                   // "CANNOTOCCUPY"
        S_K_CANPLACE = 70,                       // "CANPLACE"
        S_K_CAPACITANCE = 71,                    // "CAPACITANCE"
        S_K_CAPMULTIPLIER = 72,                  // "CAPMULTIPLIER"
        S_K_CENTERTOCENTER = 73,                 // "CENTERTOCENTER"
        S_K_CLASS = 74,                          // "CLASS"
        S_K_CLEARANCEMEASURE = 75,               // "CLEARANCEMEASURE"
        S_K_CLOCK = 76,                          // "CLOCK"
        S_K_CLOCKTYPE = 77,                      // "CLOCKTYPE"
        S_K_COLUMNMAJOR = 78,                    // "COLUMNMAJOR"
        S_K_CURRENTDEN = 79,                     // "CURRENTDEN"
        S_K_COMPONENTPIN = 80,                   // "COMPONENTPIN"
        S_K_CORE = 81,                           // "CORE"
        S_K_CORNER = 82,                         // "CORNER"
        S_K_CORRECTIONFACTOR = 83,               // "CORRECTIONFACTOR"
        S_K_CORRECTIONTABLE = 84,                // "CORRECTIONTABLE"
        S_K_COVER = 85,                          // "COVER"
        S_K_CPERSQDIST = 86,                     // "CPERSQDIST"
        S_K_CURRENT = 87,                        // "CURRENT"
        S_K_CURRENTSOURCE = 88,                  // "CURRENTSOURCE"
        S_K_CUT = 89,                            // "CUT"
        S_K_CUTAREA = 90,                        // "CUTAREA"
        S_K_CUTSIZE = 91,                        // "CUTSIZE"
        S_K_CUTSPACING = 92,                     // "CUTSPACING"
        S_K_DATA = 93,                           // "DATA"
        S_K_DATABASE = 94,                       // "DATABASE"
        S_K_DCCURRENTDENSITY = 95,               // "DCCURRENTDENSITY"
        S_K_DEFAULT = 96,                        // "DEFAULT"
        S_K_DEFAULTCAP = 97,                     // "DEFAULTCAP"
        S_K_DELAY = 98,                          // "DELAY"
        S_K_DENSITY = 99,                        // "DENSITY"
        S_K_DENSITYCHECKSTEP = 100,              // "DENSITYCHECKSTEP"
        S_K_DENSITYCHECKWINDOW = 101,            // "DENSITYCHECKWINDOW"
        S_K_DESIGNRULEWIDTH = 102,               // "DESIGNRULEWIDTH"
        S_K_DIAG45 = 103,                        // "DIAG45"
        S_K_DIAG135 = 104,                       // "DIAG135"
        S_K_DIAGMINEDGELENGTH = 105,             // "DIAGMINEDGELENGTH"
        S_K_DIAGSPACING = 106,                   // "DIAGSPACING"
        S_K_DIAGPITCH = 107,                     // "DIAGPITCH"
        S_K_DIAGWIDTH = 108,                     // "DIAGWIDTH"
        S_K_DIELECTRIC = 109,                    // "DIELECTRIC"
        S_K_DIFFUSEONLY = 110,                   // "DIFFUSEONLY"
        S_K_DIRECTION = 111,                     // "DIRECTION"
        S_K_DIVIDERCHAR = 112,                   // "DIVIDERCHAR"
        S_K_DO = 113,                            // "DO"
        S_K_E = 114,                             // "E"
        S_K_EDGECAPACITANCE = 115,               // "EDGECAPACITANCE"
        S_K_EDGERATE = 116,                      // "EDGERATE"
        S_K_EDGERATESCALEFACTOR = 117,           // "EDGERATESCALEFACTOR"
        S_K_EDGERATETHRESHOLD1 = 118,            // "EDGERATETHRESHOLD1"
        S_K_EDGERATETHRESHOLD2 = 119,            // "EDGERATETHRESHOLD2"
        S_K_EEQ = 120,                           // "EEQ"
        S_K_ELSE = 121,                          // "ELSE"
        S_K_ENCLOSURE = 122,                     // "ENCLOSURE"
        S_K_END = 123,                           // "END"
        S_K_ENDEXT = 124,                        // "ENDEXT"
        S_K_ENDCAP = 125,                        // "ENDCAP"
        S_K_ENDOFLINE = 126,                     // "ENDOFLINE"
        S_K_ENDOFNOTCHWIDTH = 127,               // "ENDOFNOTCHWIDTH"
        S_K_EUCLIDEAN = 128,                     // "EUCLIDEAN"
        S_K_EXCEPTEXTRACUT = 129,                // "EXCEPTEXTRACUT"
        S_K_EXCEPTSAMEPGNET = 130,               // "EXCEPTSAMEPGNET"
        S_K_EXCEPTPGNET = 131,                   // "EXCEPTPGNET"
        S_K_EXTENSION = 132,                     // "EXTENSION"
        S_K_FALL = 133,                          // "FALL"
        S_K_FALLCS = 134,                        // "FALLCS"
        S_K_FALLRS = 135,                        // "FALLRS"
        S_K_FALLSATCUR = 136,                    // "FALLSATCUR"
        S_K_FALLSATT1 = 137,                     // "FALLSATT1"
        S_K_FALLSLEWLIMIT = 138,                 // "FALLSLEWLIMIT"
        S_K_FALLT0 = 139,                        // "FALLT0"
        S_K_FALLTHRESH = 140,                    // "FALLTHRESH"
        S_K_FALLVOLTAGETHRESHOLD = 141,          // "FALLVOLTAGETHRESHOLD"
        S_K_FALSE = 142,                         // "FALSE"
        S_K_FE = 143,                            // "FE"
        S_K_FEEDTHRU = 144,                      // "FEEDTHRU"
        S_K_FILLACTIVESPACING = 145,             // "FILLACTIVESPACING"
        S_K_FIXED = 146,                         // "FIXED"
        S_K_FLIP = 147,                          // "FLIP"
        S_K_FLOORPLAN = 148,                     // "FLOORPLAN"
        S_K_FN = 149,                            // "FN"
        S_K_FOREIGN = 150,                       // "FOREIGN"
        S_K_FREQUENCY = 151,                     // "FREQUENCY"
        S_K_FROMABOVE = 152,                     // "FROMABOVE"
        S_K_FROMBELOW = 153,                     // "FROMBELOW"
        S_K_FROMPIN = 154,                       // "FROMPIN"
        S_K_FUNCTION = 155,                      // "FUNCTION"
        S_K_FS = 156,                            // "FS"
        S_K_FW = 157,                            // "FW"
        S_K_GCELLGRID = 158,                     // "GCELLGRID"
        S_K_GENERATE = 159,                      // "GENERATE"
        S_K_GENERATED = 160,                     // "GENERATED"
        S_K_GENERATOR = 161,                     // "GENERATOR"
        S_K_GROUND = 162,                        // "GROUND"
        S_K_GROUNDSENSITIVITY = 163,             // "GROUNDSENSITIVITY"
        S_K_HARDSPACING = 164,                   // "HARDSPACING"
        S_K_HEIGHT = 165,                        // "HEIGHT"
        S_K_HISTORY = 166,                       // "HISTORY"
        S_K_HOLD = 167,                          // "HOLD"
        S_K_HORIZONTAL = 168,                    // "HORIZONTAL"
        S_K_IF = 169,                            // "IF"
        S_K_IMPLANT = 170,                       // "IMPLANT"
        S_K_INFLUENCE = 171,                     // "INFLUENCE"
        S_K_INOUT = 172,                         // "INOUT"
        S_K_INOUTPINANTENNASIZE = 173,           // "INOUTPINANTENNASIZE"
        S_K_INPUT = 174,                         // "INPUT"
        S_K_INPUTPINANTENNASIZE = 175,           // "INPUTPINANTENNASIZE"
        S_K_INPUTNOISEMARGIN = 176,              // "INPUTNOISEMARGIN"
        S_K_INSIDECORNER = 177,                  // "INSIDECORNER"
        S_K_INTEGER = 178,                       // "INTEGER"
        S_K_INTRINSIC = 179,                     // "INTRINSIC"
        S_K_INVERT = 180,                        // "INVERT"
        S_K_INVERTER = 181,                      // "INVERTER"
        S_K_IRDROP = 182,                        // "IRDROP"
        S_K_ITERATE = 183,                       // "ITERATE"
        S_K_IV_TABLES = 184,                     // "IV_TABLES"
        S_K_LAYER = 185,                         // "LAYER"
        S_K_LAYERS = 186,                        // "LAYERS"
        S_K_LEAKAGE = 187,                       // "LEAKAGE"
        S_K_LENGTH = 188,                        // "LENGTH"
        S_K_LENGTHSUM = 189,                     // "LENGTHSUM"
        S_K_LENGTHTHRESHOLD = 190,               // "LENGTHTHRESHOLD"
        S_K_LEQ = 191,                           // "LEQ"
        S_K_LIBRARY = 192,                       // "LIBRARY"
        S_K_LONGARRAY = 193,                     // "LONGARRAY"
        S_K_MACRO = 194,                         // "MACRO"
        S_K_MANUFACTURINGGRID = 195,             // "MANUFACTURINGGRID"
        S_K_MASTERSLICE = 196,                   // "MASTERSLICE"
        S_K_MATCH = 197,                         // "MATCH"
        S_K_MAXADJACENTSLOTSPACING = 198,        // "MAXADJACENTSLOTSPACING"
        S_K_MAXCOAXIALSLOTSPACING = 199,         // "MAXCOAXIALSLOTSPACING"
        S_K_MAXDELAY = 200,                      // "MAXDELAY"
        S_K_MAXEDGES = 201,                      // "MAXEDGES"
        S_K_MAXEDGESLOTSPACING = 202,            // "MAXEDGESLOTSPACING"
        S_K_MAXLOAD = 203,                       // "MAXLOAD"
        S_K_MAXIMUMDENSITY = 204,                // "MAXIMUMDENSITY"
        S_K_MAXVIASTACK = 205,                   // "MAXVIASTACK"
        S_K_MAXWIDTH = 206,                      // "MAXWIDTH"
        S_K_MAXXY = 207,                         // "MAXXY"
        S_K_MEGAHERTZ = 208,                     // "MEGAHERTZ"
        S_K_METALOVERHANG = 209,                 // "METALOVERHANG"
        S_K_MICRONS = 210,                       // "MICRONS"
        S_K_MILLIAMPS = 211,                     // "MILLIAMPS"
        S_K_MILLIWATTS = 212,                    // "MILLIWATTS"
        S_K_MINCUTS = 213,                       // "MINCUTS"
        S_K_MINENCLOSEDAREA = 214,               // "MINENCLOSEDAREA"
        S_K_MINFEATURE = 215,                    // "MINFEATURE"
        S_K_MINIMUMCUT = 216,                    // "MINIMUMCUT"
        S_K_MINIMUMDENSITY = 217,                // "MINIMUMDENSITY"
        S_K_MINPINS = 218,                       // "MINPINS"
        S_K_MINSIZE = 219,                       // "MINSIZE"
        S_K_MINSTEP = 220,                       // "MINSTEP"
        S_K_MINWIDTH = 221,                      // "MINWIDTH"
        S_K_MPWH = 222,                          // "MPWH"
        S_K_MPWL = 223,                          // "MPWL"
        S_K_MUSTJOIN = 224,                      // "MUSTJOIN"
        S_K_MX = 225,                            // "MX"
        S_K_MY = 226,                            // "MY"
        S_K_MXR90 = 227,                         // "MXR90"
        S_K_MYR90 = 228,                         // "MYR90"
        S_K_N = 229,                             // "N"
        S_K_NAMEMAPSTRING = 230,                 // "NAMEMAPSTRING"
        S_K_NAMESCASESENSITIVE = 231,            // "NAMESCASESENSITIVE"
        S_K_NANOSECONDS = 232,                   // "NANOSECONDS"
        S_K_NEGEDGE = 233,                       // "NEGEDGE"
        S_K_NETEXPR = 234,                       // "NETEXPR"
        S_K_NETS = 235,                          // "NETS"
        S_K_NEW = 236,                           // "NEW"
        S_K_NONDEFAULTRULE = 237,                // "NONDEFAULTRULE"
        S_K_NONE = 238,                          // "NONE"
        S_K_NONINVERT = 239,                     // "NONINVERT"
        S_K_NONUNATE = 240,                      // "NONUNATE"
        S_K_NOISETABLE = 241,                    // "NOISETABLE"
        S_K_NOTCHLENGTH = 242,                   // "NOTCHLENGTH"
        S_K_NOTCHSPACING = 243,                  // "NOTCHSPACING"
        S_K_NOWIREEXTENSIONATPIN = 244,          // "NOWIREEXTENSIONATPIN"
        S_K_OBS = 245,                           // "OBS"
        S_K_OFF = 246,                           // "OFF"
        S_K_OFFSET = 247,                        // "OFFSET"
        S_K_OHMS = 248,                          // "OHMS"
        S_K_ON = 249,                            // "ON"
        S_K_OR = 250,                            // "OR"
        S_K_ORIENT = 251,                        // "ORIENT"
        S_K_ORIENTATION = 252,                   // "ORIENTATION"
        S_K_ORIGIN = 253,                        // "ORIGIN"
        S_K_ORTHOGONAL = 254,                    // "ORTHOGONAL"
        S_K_OUTPUT = 255,                        // "OUTPUT"
        S_K_OUTPUTPINANTENNASIZE = 256,          // "OUTPUTPINANTENNASIZE"
        S_K_OUTPUTNOISEMARGIN = 257,             // "OUTPUTNOISEMARGIN"
        S_K_OUTPUTRESISTANCE = 258,              // "OUTPUTRESISTANCE"
        S_K_OUTSIDECORNER = 259,                 // "OUTSIDECORNER"
        S_K_OVERHANG = 260,                      // "OVERHANG"
        S_K_OVERLAP = 261,                       // "OVERLAP"
        S_K_OVERLAPS = 262,                      // "OVERLAPS"
        S_K_OXIDE1 = 263,                        // "OXIDE1"
        S_K_OXIDE2 = 264,                        // "OXIDE2"
        S_K_OXIDE3 = 265,                        // "OXIDE3"
        S_K_OXIDE4 = 266,                        // "OXIDE4"
        S_K_PAD = 267,                           // "PAD"
        S_K_PARALLELEDGE = 268,                  // "PARALLELEDGE"
        S_K_PARALLELOVERLAP = 269,               // "PARALLELOVERLAP"
        S_K_PARALLELRUNLENGTH = 270,             // "PARALLELRUNLENGTH"
        S_K_PATH = 271,                          // "PATH"
        S_K_PATTERN = 272,                       // "PATTERN"
        S_K_PEAK = 273,                          // "PEAK"
        S_K_PERIOD = 274,                        // "PERIOD"
        S_K_PGONLY = 275,                        // "PGONLY"
        S_K_PICOFARADS = 276,                    // "PICOFARADS"
        S_K_PIN = 277,                           // "PIN"
        S_K_PITCH = 278,                         // "PITCH"
        S_K_PLACED = 279,                        // "PLACED"
        S_K_POLYGON = 280,                       // "POLYGON"
        S_K_PORT = 281,                          // "PORT"
        S_K_POSEDGE = 282,                       // "POSEDGE"
        S_K_POST = 283,                          // "POST"
        S_K_POWER = 284,                         // "POWER"
        S_K_PRE = 285,                           // "PRE"
        S_K_PREFERENCLOSURE = 286,               // "PREFERENCLOSURE"
        S_K_PRL = 287,                           // "PRL"
        S_K_PROPERTY = 288,                      // "PROPERTY"
        S_K_PROPERTYDEFINITIONS = 289,           // "PROPERTYDEFINITIONS"
        S_K_PROTRUSIONWIDTH = 290,               // "PROTRUSIONWIDTH"
        S_K_PULLDOWNRES = 291,                   // "PULLDOWNRES"
        S_K_PWL = 292,                           // "PWL"
        S_K_R0 = 293,                            // "R0"
        S_K_R90 = 294,                           // "R90"
        S_K_R180 = 295,                          // "R180"
        S_K_R270 = 296,                          // "R270"
        S_K_RANGE = 297,                         // "RANGE"
        S_K_REAL = 298,                          // "REAL"
        S_K_RECOVERY = 299,                      // "RECOVERY"
        S_K_RECT = 300,                          // "RECT"
        S_K_RESISTANCE = 301,                    // "RESISTANCE"
        S_K_RESISTIVE = 302,                     // "RESISTIVE"
        S_K_RING = 303,                          // "RING"
        S_K_RISE = 304,                          // "RISE"
        S_K_RISECS = 305,                        // "RISECS"
        S_K_RISERS = 306,                        // "RISERS"
        S_K_RISESATCUR = 307,                    // "RISESATCUR"
        S_K_RISESATT1 = 308,                     // "RISESATT1"
        S_K_RISESLEWLIMIT = 309,                 // "RISESLEWLIMIT"
        S_K_RISET0 = 310,                        // "RISET0"
        S_K_RISETHRESH = 311,                    // "RISETHRESH"
        S_K_RISEVOLTAGETHRESHOLD = 312,          // "RISEVOLTAGETHRESHOLD"
        S_K_RMS = 313,                           // "RMS"
        S_K_ROUTING = 314,                       // "ROUTING"
        S_K_ROWABUTSPACING = 315,                // "ROWABUTSPACING"
        S_K_ROWCOL = 316,                        // "ROWCOL"
        S_K_ROWMAJOR = 317,                      // "ROWMAJOR"
        S_K_ROWMINSPACING = 318,                 // "ROWMINSPACING"
        S_K_ROWPATTERN = 319,                    // "ROWPATTERN"
        S_K_RPERSQ = 320,                        // "RPERSQ"
        S_K_S = 321,                             // "S"
        S_K_SAMENET = 322,                       // "SAMENET"
        S_K_SCANUSE = 323,                       // "SCANUSE"
        S_K_SDFCOND = 324,                       // "SDFCOND"
        S_K_SDFCONDEND = 325,                    // "SDFCONDEND"
        S_K_SDFCONDSTART = 326,                  // "SDFCONDSTART"
        S_K_SETUP = 327,                         // "SETUP"
        S_K_SHAPE = 328,                         // "SHAPE"
        S_K_SHRINKAGE = 329,                     // "SHRINKAGE"
        S_K_SIGNAL = 330,                        // "SIGNAL"
        S_K_SITE = 331,                          // "SITE"
        S_K_SIZE = 332,                          // "SIZE"
        S_K_SKEW = 333,                          // "SKEW"
        S_K_SLOTLENGTH = 334,                    // "SLOTLENGTH"
        S_K_SLOTWIDTH = 335,                     // "SLOTWIDTH"
        S_K_SLOTWIRELENGTH = 336,                // "SLOTWIRELENGTH"
        S_K_SLOTWIREWIDTH = 337,                 // "SLOTWIREWIDTH"
        S_K_SPLITWIREWIDTH = 338,                // "SPLITWIREWIDTH"
        S_K_SOFT = 339,                          // "SOFT"
        S_K_SOURCE = 340,                        // "SOURCE"
        S_K_SPACER = 341,                        // "SPACER"
        S_K_SPACING = 342,                       // "SPACING"
        S_K_SPACINGTABLE = 343,                  // "SPACINGTABLE"
        S_K_SPECIALNETS = 344,                   // "SPECIALNETS"
        S_K_STABLE = 345,                        // "STABLE"
        S_K_STACK = 346,                         // "STACK"
        S_K_START = 347,                         // "START"
        S_K_STEP = 348,                          // "STEP"
        S_K_STOP = 349,                          // "STOP"
        S_K_STRING = 350,                        // "STRING"
        S_K_STRUCTURE = 351,                     // "STRUCTURE"
        S_K_SUPPLYSENSITIVITY = 352,             // "SUPPLYSENSITIVITY"
        S_K_SYMMETRY = 353,                      // "SYMMETRY"
        S_K_TABLE = 354,                         // "TABLE"
        S_K_TABLEAXIS = 355,                     // "TABLEAXIS"
        S_K_TABLEDIMENSION = 356,                // "TABLEDIMENSION"
        S_K_TABLEENTRIES = 357,                  // "TABLEENTRIES"
        S_K_TAPERRULE = 358,                     // "TAPERRULE"
        S_K_THEN = 359,                          // "THEN"
        S_K_THICKNESS = 360,                     // "THICKNESS"
        S_K_TIEHIGH = 361,                       // "TIEHIGH"
        S_K_TIELOW = 362,                        // "TIELOW"
        S_K_TIEOFFR = 363,                       // "TIEOFFR"
        S_K_TIME = 364,                          // "TIME"
        S_K_TIMING = 365,                        // "TIMING"
        S_K_TO = 366,                            // "TO"
        S_K_TOPIN = 367,                         // "TOPIN"
        S_K_TOPLEFT = 368,                       // "TOPLEFT"
        S_K_TOPOFSTACKONLY = 369,                // "TOPOFSTACKONLY"
        S_K_TOPRIGHT = 370,                      // "TOPRIGHT"
        S_K_TRACKS = 371,                        // "TRACKS"
        S_K_TRANSITIONTIME = 372,                // "TRANSITIONTIME"
        S_K_TRISTATE = 373,                      // "TRISTATE"
        S_K_TRUE = 374,                          // "TRUE"
        S_K_TWOEDGES = 375,                      // "TWOEDGES"
        S_K_TWOWIDTHS = 376,                     // "TWOWIDTHS"
        S_K_TYPE = 377,                          // "TYPE"
        S_K_UNATENESS = 378,                     // "UNATENESS"
        S_K_UNITS = 379,                         // "UNITS"
        S_K_UNIVERSALNOISEMARGIN = 380,          // "UNIVERSALNOISEMARGIN"
        S_K_USE = 381,                           // "USE"
        S_K_USELENGTHTHRESHOLD = 382,            // "USELENGTHTHRESHOLD"
        S_K_USEMINSPACING = 383,                 // "USEMINSPACING"
        S_K_USER = 384,                          // "USER"
        S_K_USEVIA = 385,                        // "USEVIA"
        S_K_USEVIARULE = 386,                    // "USEVIARULE"
        S_K_VARIABLE = 387,                      // "VARIABLE"
        S_K_VERSION = 388,                       // "VERSION"
        S_K_VERTICAL = 389,                      // "VERTICAL"
        S_K_VHI = 390,                           // "VHI"
        S_K_VIA = 391,                           // "VIA"
        S_K_VIARULE = 392,                       // "VIARULE"
        S_K_VICTIMLENGTH = 393,                  // "VICTIMLENGTH"
        S_K_VICTIMNOISE = 394,                   // "VICTIMNOISE"
        S_K_VIRTUAL = 395,                       // "VIRTUAL"
        S_K_VLO = 396,                           // "VLO"
        S_K_VOLTAGE = 397,                       // "VOLTAGE"
        S_K_VOLTS = 398,                         // "VOLTS"
        S_K_W = 399,                             // "W"
        S_K_WELLTAP = 400,                       // "WELLTAP"
        S_K_WIDTH = 401,                         // "WIDTH"
        S_K_WITHIN = 402,                        // "WITHIN"
        S_K_WIRECAP = 403,                       // "WIRECAP"
        S_K_WIREEXTENSION = 404,                 // "WIREEXTENSION"
        S_K_X = 405,                             // "X"
        S_K_Y = 406,                             // "Y"
        S_K_EQ = 407,                            // K_EQ
        S_K_NE = 408,                            // K_NE
        S_K_LE = 409,                            // K_LE
        S_K_LT = 410,                            // K_LT
        S_K_GE = 411,                            // K_GE
        S_K_GT = 412,                            // K_GT
        S_K_NOT = 413,                           // K_NOT
        S_INTEGER = 414,                         // "integer"
        S_DOUBLE = 415,                          // "double"
        S_STRING = 416,                          // "string"
        S_QSTRING = 417,                         // "qstring"
        S_BINARY = 418,                          // "binary numbers"
        S_419_number_ = 419,                     // "number"
        S_420_generalized_string_ = 420,         // "generalized string"
        S_IF = 421,                              // IF
        S_LNOT = 422,                            // LNOT
        S_423_ = 423,                            // '-'
        S_424_ = 424,                            // '+'
        S_425_ = 425,                            // '*'
        S_426_ = 426,                            // '/'
        S_UMINUS = 427,                          // UMINUS
        S_428_ = 428,                            // ';'
        S_429_ = 429,                            // '('
        S_430_ = 430,                            // ')'
        S_431_ = 431,                            // '='
        S_432_n_ = 432,                          // '\n'
        S_433_ = 433,                            // '<'
        S_434_ = 434,                            // '>'
        S_YYACCEPT = 435,                        // $accept
        S_NUMBER = 436,                          // NUMBER
        S_GSTRING = 437,                         // GSTRING
        S_lef_file = 438,                        // lef_file
        S_version = 439,                         // version
        S_440_1 = 440,                           // $@1
        S_441_2 = 441,                           // $@2
        S_dividerchar = 442,                     // dividerchar
        S_busbitchars = 443,                     // busbitchars
        S_rules = 444,                           // rules
        S_end_library = 445,                     // end_library
        S_rule = 446,                            // rule
        S_case_sensitivity = 447,                // case_sensitivity
        S_wireextension = 448,                   // wireextension
        S_manufacturing = 449,                   // manufacturing
        S_useminspacing = 450,                   // useminspacing
        S_clearancemeasure = 451,                // clearancemeasure
        S_clearance_type = 452,                  // clearance_type
        S_spacing_type = 453,                    // spacing_type
        S_spacing_value = 454,                   // spacing_value
        S_units_section = 455,                   // units_section
        S_start_units = 456,                     // start_units
        S_units_rules = 457,                     // units_rules
        S_units_rule = 458,                      // units_rule
        S_layer_rule = 459,                      // layer_rule
        S_start_layer = 460,                     // start_layer
        S_461_3 = 461,                           // $@3
        S_end_layer = 462,                       // end_layer
        S_463_4 = 463,                           // $@4
        S_layer_options = 464,                   // layer_options
        S_layer_option = 465,                    // layer_option
        S_466_5 = 466,                           // $@5
        S_467_6 = 467,                           // $@6
        S_468_7 = 468,                           // $@7
        S_469_8 = 469,                           // $@8
        S_470_9 = 470,                           // $@9
        S_471_10 = 471,                          // $@10
        S_472_11 = 472,                          // $@11
        S_473_12 = 473,                          // $@12
        S_474_13 = 474,                          // $@13
        S_475_14 = 475,                          // $@14
        S_476_15 = 476,                          // $@15
        S_477_16 = 477,                          // $@16
        S_478_17 = 478,                          // $@17
        S_479_18 = 479,                          // $@18
        S_480_19 = 480,                          // $@19
        S_481_20 = 481,                          // $@20
        S_482_21 = 482,                          // $@21
        S_483_22 = 483,                          // $@22
        S_484_23 = 484,                          // $@23
        S_485_24 = 485,                          // $@24
        S_486_25 = 486,                          // $@25
        S_487_26 = 487,                          // $@26
        S_488_27 = 488,                          // $@27
        S_489_28 = 489,                          // $@28
        S_490_29 = 490,                          // $@29
        S_491_30 = 491,                          // $@30
        S_layer_arraySpacing_long = 492,         // layer_arraySpacing_long
        S_layer_arraySpacing_width = 493,        // layer_arraySpacing_width
        S_layer_arraySpacing_arraycuts = 494,    // layer_arraySpacing_arraycuts
        S_layer_arraySpacing_arraycut = 495,     // layer_arraySpacing_arraycut
        S_sp_options = 496,                      // sp_options
        S_497_31 = 497,                          // $@31
        S_498_32 = 498,                          // $@32
        S_499_33 = 499,                          // $@33
        S_500_34 = 500,                          // $@34
        S_501_35 = 501,                          // $@35
        S_502_36 = 502,                          // $@36
        S_503_37 = 503,                          // $@37
        S_layer_spacingtable_opts = 504,         // layer_spacingtable_opts
        S_layer_spacingtable_opt = 505,          // layer_spacingtable_opt
        S_layer_enclosure_type_opt = 506,        // layer_enclosure_type_opt
        S_layer_enclosure_width_opt = 507,       // layer_enclosure_width_opt
        S_508_38 = 508,                          // $@38
        S_layer_enclosure_width_except_opt = 509, // layer_enclosure_width_except_opt
        S_layer_preferenclosure_width_opt = 510, // layer_preferenclosure_width_opt
        S_layer_minimumcut_within = 511,         // layer_minimumcut_within
        S_layer_minimumcut_from = 512,           // layer_minimumcut_from
        S_layer_minimumcut_length = 513,         // layer_minimumcut_length
        S_layer_minstep_options = 514,           // layer_minstep_options
        S_layer_minstep_option = 515,            // layer_minstep_option
        S_layer_minstep_type = 516,              // layer_minstep_type
        S_layer_antenna_pwl = 517,               // layer_antenna_pwl
        S_518_39 = 518,                          // $@39
        S_layer_diffusion_ratios = 519,          // layer_diffusion_ratios
        S_layer_diffusion_ratio = 520,           // layer_diffusion_ratio
        S_layer_antenna_duo = 521,               // layer_antenna_duo
        S_layer_table_type = 522,                // layer_table_type
        S_layer_frequency = 523,                 // layer_frequency
        S_524_40 = 524,                          // $@40
        S_525_41 = 525,                          // $@41
        S_526_42 = 526,                          // $@42
        S_ac_layer_table_opt = 527,              // ac_layer_table_opt
        S_528_43 = 528,                          // $@43
        S_529_44 = 529,                          // $@44
        S_dc_layer_table = 530,                  // dc_layer_table
        S_531_45 = 531,                          // $@45
        S_number_list = 532,                     // number_list
        S_layer_prop_list = 533,                 // layer_prop_list
        S_layer_prop = 534,                      // layer_prop
        S_current_density_pwl_list = 535,        // current_density_pwl_list
        S_current_density_pwl = 536,             // current_density_pwl
        S_cap_points = 537,                      // cap_points
        S_cap_point = 538,                       // cap_point
        S_res_points = 539,                      // res_points
        S_res_point = 540,                       // res_point
        S_layer_type = 541,                      // layer_type
        S_layer_direction = 542,                 // layer_direction
        S_layer_minen_width = 543,               // layer_minen_width
        S_layer_oxide = 544,                     // layer_oxide
        S_layer_sp_parallel_widths = 545,        // layer_sp_parallel_widths
        S_layer_sp_parallel_width = 546,         // layer_sp_parallel_width
        S_547_46 = 547,                          // $@46
        S_layer_sp_TwoWidths = 548,              // layer_sp_TwoWidths
        S_layer_sp_TwoWidth = 549,               // layer_sp_TwoWidth
        S_550_47 = 550,                          // $@47
        S_layer_sp_TwoWidthsPRL = 551,           // layer_sp_TwoWidthsPRL
        S_layer_sp_influence_widths = 552,       // layer_sp_influence_widths
        S_layer_sp_influence_width = 553,        // layer_sp_influence_width
        S_maxstack_via = 554,                    // maxstack_via
        S_555_48 = 555,                          // $@48
        S_via = 556,                             // via
        S_557_49 = 557,                          // $@49
        S_via_keyword = 558,                     // via_keyword
        S_start_via = 559,                       // start_via
        S_via_viarule = 560,                     // via_viarule
        S_561_50 = 561,                          // $@50
        S_562_51 = 562,                          // $@51
        S_563_52 = 563,                          // $@52
        S_via_viarule_options = 564,             // via_viarule_options
        S_via_viarule_option = 565,              // via_viarule_option
        S_566_53 = 566,                          // $@53
        S_via_option = 567,                      // via_option
        S_via_other_options = 568,               // via_other_options
        S_via_more_options = 569,                // via_more_options
        S_via_other_option = 570,                // via_other_option
        S_571_54 = 571,                          // $@54
        S_via_prop_list = 572,                   // via_prop_list
        S_via_name_value_pair = 573,             // via_name_value_pair
        S_via_foreign = 574,                     // via_foreign
        S_start_foreign = 575,                   // start_foreign
        S_576_55 = 576,                          // $@55
        S_orientation = 577,                     // orientation
        S_via_layer_rule = 578,                  // via_layer_rule
        S_via_layer = 579,                       // via_layer
        S_580_56 = 580,                          // $@56
        S_via_geometries = 581,                  // via_geometries
        S_via_geometry = 582,                    // via_geometry
        S_583_57 = 583,                          // $@57
        S_end_via = 584,                         // end_via
        S_585_58 = 585,                          // $@58
        S_viarule_keyword = 586,                 // viarule_keyword
        S_587_59 = 587,                          // $@59
        S_viarule = 588,                         // viarule
        S_viarule_generate = 589,                // viarule_generate
        S_590_60 = 590,                          // $@60
        S_viarule_generate_default = 591,        // viarule_generate_default
        S_viarule_layer_list = 592,              // viarule_layer_list
        S_opt_viarule_props = 593,               // opt_viarule_props
        S_viarule_props = 594,                   // viarule_props
        S_viarule_prop = 595,                    // viarule_prop
        S_596_61 = 596,                          // $@61
        S_viarule_prop_list = 597,               // viarule_prop_list
        S_viarule_layer = 598,                   // viarule_layer
        S_via_names = 599,                       // via_names
        S_via_name = 600,                        // via_name
        S_viarule_layer_name = 601,              // viarule_layer_name
        S_602_62 = 602,                          // $@62
        S_viarule_layer_options = 603,           // viarule_layer_options
        S_viarule_layer_option = 604,            // viarule_layer_option
        S_end_viarule = 605,                     // end_viarule
        S_606_63 = 606,                          // $@63
        S_spacing_rule = 607,                    // spacing_rule
        S_start_spacing = 608,                   // start_spacing
        S_end_spacing = 609,                     // end_spacing
        S_spacings = 610,                        // spacings
        S_spacing = 611,                         // spacing
        S_samenet_keyword = 612,                 // samenet_keyword
        S_irdrop = 613,                          // irdrop
        S_start_irdrop = 614,                    // start_irdrop
        S_end_irdrop = 615,                      // end_irdrop
        S_ir_tables = 616,                       // ir_tables
        S_ir_table = 617,                        // ir_table
        S_ir_table_values = 618,                 // ir_table_values
        S_ir_table_value = 619,                  // ir_table_value
        S_ir_tablename = 620,                    // ir_tablename
        S_minfeature = 621,                      // minfeature
        S_dielectric = 622,                      // dielectric
        S_nondefault_rule = 623,                 // nondefault_rule
        S_624_64 = 624,                          // $@64
        S_625_65 = 625,                          // $@65
        S_626_66 = 626,                          // $@66
        S_end_nd_rule = 627,                     // end_nd_rule
        S_nd_hardspacing = 628,                  // nd_hardspacing
        S_nd_rules = 629,                        // nd_rules
        S_nd_rule = 630,                         // nd_rule
        S_usevia = 631,                          // usevia
        S_useviarule = 632,                      // useviarule
        S_mincuts = 633,                         // mincuts
        S_nd_prop = 634,                         // nd_prop
        S_635_67 = 635,                          // $@67
        S_nd_prop_list = 636,                    // nd_prop_list
        S_nd_layer = 637,                        // nd_layer
        S_638_68 = 638,                          // $@68
        S_639_69 = 639,                          // $@69
        S_640_70 = 640,                          // $@70
        S_641_71 = 641,                          // $@71
        S_nd_layer_stmts = 642,                  // nd_layer_stmts
        S_nd_layer_stmt = 643,                   // nd_layer_stmt
        S_site = 644,                            // site
        S_start_site = 645,                      // start_site
        S_646_72 = 646,                          // $@72
        S_end_site = 647,                        // end_site
        S_648_73 = 648,                          // $@73
        S_site_options = 649,                    // site_options
        S_site_option = 650,                     // site_option
        S_site_class = 651,                      // site_class
        S_site_symmetry_statement = 652,         // site_symmetry_statement
        S_site_symmetries = 653,                 // site_symmetries
        S_site_symmetry = 654,                   // site_symmetry
        S_site_rowpattern_statement = 655,       // site_rowpattern_statement
        S_656_74 = 656,                          // $@74
        S_site_rowpatterns = 657,                // site_rowpatterns
        S_site_rowpattern = 658,                 // site_rowpattern
        S_659_75 = 659,                          // $@75
        S_pt = 660,                              // pt
        S_macro = 661,                           // macro
        S_662_76 = 662,                          // $@76
        S_start_macro = 663,                     // start_macro
        S_664_77 = 664,                          // $@77
        S_end_macro = 665,                       // end_macro
        S_666_78 = 666,                          // $@78
        S_macro_options = 667,                   // macro_options
        S_macro_option = 668,                    // macro_option
        S_669_79 = 669,                          // $@79
        S_macro_prop_list = 670,                 // macro_prop_list
        S_macro_symmetry_statement = 671,        // macro_symmetry_statement
        S_macro_symmetries = 672,                // macro_symmetries
        S_macro_symmetry = 673,                  // macro_symmetry
        S_macro_name_value_pair = 674,           // macro_name_value_pair
        S_macro_class = 675,                     // macro_class
        S_class_type = 676,                      // class_type
        S_pad_type = 677,                        // pad_type
        S_core_type = 678,                       // core_type
        S_endcap_type = 679,                     // endcap_type
        S_macro_generator = 680,                 // macro_generator
        S_macro_generate = 681,                  // macro_generate
        S_macro_source = 682,                    // macro_source
        S_macro_power = 683,                     // macro_power
        S_macro_origin = 684,                    // macro_origin
        S_macro_foreign = 685,                   // macro_foreign
        S_macro_eeq = 686,                       // macro_eeq
        S_687_80 = 687,                          // $@80
        S_macro_leq = 688,                       // macro_leq
        S_689_81 = 689,                          // $@81
        S_macro_site = 690,                      // macro_site
        S_macro_site_word = 691,                 // macro_site_word
        S_site_word = 692,                       // site_word
        S_macro_size = 693,                      // macro_size
        S_macro_pin = 694,                       // macro_pin
        S_start_macro_pin = 695,                 // start_macro_pin
        S_696_82 = 696,                          // $@82
        S_end_macro_pin = 697,                   // end_macro_pin
        S_698_83 = 698,                          // $@83
        S_macro_pin_options = 699,               // macro_pin_options
        S_macro_pin_option = 700,                // macro_pin_option
        S_701_84 = 701,                          // $@84
        S_702_85 = 702,                          // $@85
        S_703_86 = 703,                          // $@86
        S_704_87 = 704,                          // $@87
        S_705_88 = 705,                          // $@88
        S_706_89 = 706,                          // $@89
        S_707_90 = 707,                          // $@90
        S_708_91 = 708,                          // $@91
        S_709_92 = 709,                          // $@92
        S_710_93 = 710,                          // $@93
        S_pin_layer_oxide = 711,                 // pin_layer_oxide
        S_pin_prop_list = 712,                   // pin_prop_list
        S_pin_name_value_pair = 713,             // pin_name_value_pair
        S_electrical_direction = 714,            // electrical_direction
        S_start_macro_port = 715,                // start_macro_port
        S_macro_port_class_option = 716,         // macro_port_class_option
        S_macro_pin_use = 717,                   // macro_pin_use
        S_macro_scan_use = 718,                  // macro_scan_use
        S_pin_shape = 719,                       // pin_shape
        S_geometries = 720,                      // geometries
        S_geometry = 721,                        // geometry
        S_722_94 = 722,                          // $@94
        S_723_95 = 723,                          // $@95
        S_geometry_options = 724,                // geometry_options
        S_layer_exceptpgnet = 725,               // layer_exceptpgnet
        S_layer_spacing = 726,                   // layer_spacing
        S_firstPt = 727,                         // firstPt
        S_nextPt = 728,                          // nextPt
        S_otherPts = 729,                        // otherPts
        S_via_placement = 730,                   // via_placement
        S_731_96 = 731,                          // $@96
        S_732_97 = 732,                          // $@97
        S_stepPattern = 733,                     // stepPattern
        S_sitePattern = 734,                     // sitePattern
        S_trackPattern = 735,                    // trackPattern
        S_736_98 = 736,                          // $@98
        S_737_99 = 737,                          // $@99
        S_738_100 = 738,                         // $@100
        S_739_101 = 739,                         // $@101
        S_trackLayers = 740,                     // trackLayers
        S_layer_name = 741,                      // layer_name
        S_gcellPattern = 742,                    // gcellPattern
        S_macro_obs = 743,                       // macro_obs
        S_start_macro_obs = 744,                 // start_macro_obs
        S_macro_density = 745,                   // macro_density
        S_density_layers = 746,                  // density_layers
        S_density_layer = 747,                   // density_layer
        S_748_102 = 748,                         // $@102
        S_749_103 = 749,                         // $@103
        S_density_layer_rects = 750,             // density_layer_rects
        S_density_layer_rect = 751,              // density_layer_rect
        S_macro_clocktype = 752,                 // macro_clocktype
        S_753_104 = 753,                         // $@104
        S_timing = 754,                          // timing
        S_start_timing = 755,                    // start_timing
        S_end_timing = 756,                      // end_timing
        S_timing_options = 757,                  // timing_options
        S_timing_option = 758,                   // timing_option
        S_759_105 = 759,                         // $@105
        S_760_106 = 760,                         // $@106
        S_761_107 = 761,                         // $@107
        S_one_pin_trigger = 762,                 // one_pin_trigger
        S_two_pin_trigger = 763,                 // two_pin_trigger
        S_from_pin_trigger = 764,                // from_pin_trigger
        S_to_pin_trigger = 765,                  // to_pin_trigger
        S_delay_or_transition = 766,             // delay_or_transition
        S_list_of_table_entries = 767,           // list_of_table_entries
        S_table_entry = 768,                     // table_entry
        S_list_of_table_axis_numbers = 769,      // list_of_table_axis_numbers
        S_slew_spec = 770,                       // slew_spec
        S_risefall = 771,                        // risefall
        S_unateness = 772,                       // unateness
        S_list_of_from_strings = 773,            // list_of_from_strings
        S_list_of_to_strings = 774,              // list_of_to_strings
        S_array = 775,                           // array
        S_776_108 = 776,                         // $@108
        S_start_array = 777,                     // start_array
        S_778_109 = 778,                         // $@109
        S_end_array = 779,                       // end_array
        S_780_110 = 780,                         // $@110
        S_array_rules = 781,                     // array_rules
        S_array_rule = 782,                      // array_rule
        S_783_111 = 783,                         // $@111
        S_784_112 = 784,                         // $@112
        S_785_113 = 785,                         // $@113
        S_786_114 = 786,                         // $@114
        S_787_115 = 787,                         // $@115
        S_floorplan_start = 788,                 // floorplan_start
        S_floorplan_list = 789,                  // floorplan_list
        S_floorplan_element = 790,               // floorplan_element
        S_791_116 = 791,                         // $@116
        S_792_117 = 792,                         // $@117
        S_cap_list = 793,                        // cap_list
        S_one_cap = 794,                         // one_cap
        S_msg_statement = 795,                   // msg_statement
        S_796_118 = 796,                         // $@118
        S_create_file_statement = 797,           // create_file_statement
        S_798_119 = 798,                         // $@119
        S_def_statement = 799,                   // def_statement
        S_800_120 = 800,                         // $@120
        S_801_121 = 801,                         // $@121
        S_802_122 = 802,                         // $@122
        S_dtrm = 803,                            // dtrm
        S_then = 804,                            // then
        S_else = 805,                            // else
        S_expression = 806,                      // expression
        S_b_expr = 807,                          // b_expr
        S_s_expr = 808,                          // s_expr
        S_relop = 809,                           // relop
        S_prop_def_section = 810,                // prop_def_section
        S_811_123 = 811,                         // $@123
        S_prop_stmts = 812,                      // prop_stmts
        S_prop_stmt = 813,                       // prop_stmt
        S_814_124 = 814,                         // $@124
        S_815_125 = 815,                         // $@125
        S_816_126 = 816,                         // $@126
        S_817_127 = 817,                         // $@127
        S_818_128 = 818,                         // $@128
        S_819_129 = 819,                         // $@129
        S_820_130 = 820,                         // $@130
        S_821_131 = 821,                         // $@131
        S_prop_define = 822,                     // prop_define
        S_823_132 = 823,                         // $@132
        S_824_133 = 824,                         // $@133
        S_opt_range_second = 825,                // opt_range_second
        S_opt_endofline = 826,                   // opt_endofline
        S_827_134 = 827,                         // $@134
        S_opt_endofline_twoedges = 828,          // opt_endofline_twoedges
        S_opt_samenetPGonly = 829,               // opt_samenetPGonly
        S_opt_def_range = 830,                   // opt_def_range
        S_opt_def_value = 831,                   // opt_def_value
        S_layer_spacing_opts = 832,              // layer_spacing_opts
        S_layer_spacing_opt = 833,               // layer_spacing_opt
        S_834_135 = 834,                         // $@135
        S_layer_spacing_cut_routing = 835,       // layer_spacing_cut_routing
        S_836_136 = 836,                         // $@136
        S_837_137 = 837,                         // $@137
        S_838_138 = 838,                         // $@138
        S_839_139 = 839,                         // $@139
        S_840_140 = 840,                         // $@140
        S_spacing_cut_layer_opt = 841,           // spacing_cut_layer_opt
        S_opt_adjacentcuts_exceptsame = 842,     // opt_adjacentcuts_exceptsame
        S_opt_layer_name = 843,                  // opt_layer_name
        S_844_141 = 844,                         // $@141
        S_req_layer_name = 845,                  // req_layer_name
        S_846_142 = 846,                         // $@142
        S_universalnoisemargin = 847,            // universalnoisemargin
        S_edgeratethreshold1 = 848,              // edgeratethreshold1
        S_edgeratethreshold2 = 849,              // edgeratethreshold2
        S_edgeratescalefactor = 850,             // edgeratescalefactor
        S_noisetable = 851,                      // noisetable
        S_852_143 = 852,                         // $@143
        S_end_noisetable = 853,                  // end_noisetable
        S_noise_table_list = 854,                // noise_table_list
        S_noise_table_entry = 855,               // noise_table_entry
        S_output_resistance_entry = 856,         // output_resistance_entry
        S_857_144 = 857,                         // $@144
        S_num_list = 858,                        // num_list
        S_victim_list = 859,                     // victim_list
        S_victim = 860,                          // victim
        S_861_145 = 861,                         // $@145
        S_vnoiselist = 862,                      // vnoiselist
        S_correctiontable = 863,                 // correctiontable
        S_864_146 = 864,                         // $@146
        S_end_correctiontable = 865,             // end_correctiontable
        S_correction_table_list = 866,           // correction_table_list
        S_correction_table_item = 867,           // correction_table_item
        S_output_list = 868,                     // output_list
        S_869_147 = 869,                         // $@147
        S_numo_list = 870,                       // numo_list
        S_corr_victim_list = 871,                // corr_victim_list
        S_corr_victim = 872,                     // corr_victim
        S_873_148 = 873,                         // $@148
        S_corr_list = 874,                       // corr_list
        S_input_antenna = 875,                   // input_antenna
        S_output_antenna = 876,                  // output_antenna
        S_inout_antenna = 877,                   // inout_antenna
        S_antenna_input = 878,                   // antenna_input
        S_antenna_inout = 879,                   // antenna_inout
        S_antenna_output = 880,                  // antenna_output
        S_extension_opt = 881,                   // extension_opt
        S_extension = 882                        // extension
      };
    };

    /// (Internal) symbol kind.
    typedef symbol_kind::symbol_kind_type symbol_kind_type;

    /// The number of tokens.
    static const symbol_kind_type YYNTOKENS = symbol_kind::YYNTOKENS;

    /// A complete symbol.
    ///
    /// Expects its Base type to provide access to the symbol kind
    /// via kind ().
    ///
    /// Provide access to semantic value and location.
    template <typename Base>
    struct basic_symbol : Base
    {
      /// Alias to Base.
      typedef Base super_type;

      /// Default constructor.
      basic_symbol () YY_NOEXCEPT
        : value ()
        , location ()
      {}

#if 201103L <= YY_CPLUSPLUS
      /// Move constructor.
      basic_symbol (basic_symbol&& that)
        : Base (std::move (that))
        , value (std::move (that.value))
        , location (std::move (that.location))
      {}
#endif

      /// Copy constructor.
      basic_symbol (const basic_symbol& that);
      /// Constructor for valueless symbols.
      basic_symbol (typename Base::kind_type t,
                    YY_MOVE_REF (location_type) l);

      /// Constructor for symbols with semantic value.
      basic_symbol (typename Base::kind_type t,
                    YY_RVREF (value_type) v,
                    YY_RVREF (location_type) l);

      /// Destroy the symbol.
      ~basic_symbol ()
      {
        clear ();
      }



      /// Destroy contents, and record that is empty.
      void clear () YY_NOEXCEPT
      {
        Base::clear ();
      }

      /// The user-facing name of this symbol.
      std::string name () const YY_NOEXCEPT
      {
        return Parser::symbol_name (this->kind ());
      }

      /// Backward compatibility (Bison 3.6).
      symbol_kind_type type_get () const YY_NOEXCEPT;

      /// Whether empty.
      bool empty () const YY_NOEXCEPT;

      /// Destructive move, \a s is emptied into this.
      void move (basic_symbol& s);

      /// The semantic value.
      value_type value;

      /// The location.
      location_type location;

    private:
#if YY_CPLUSPLUS < 201103L
      /// Assignment operator.
      basic_symbol& operator= (const basic_symbol& that);
#endif
    };

    /// Type access provider for token (enum) based symbols.
    struct by_kind
    {
      /// The symbol kind as needed by the constructor.
      typedef token_kind_type kind_type;

      /// Default constructor.
      by_kind () YY_NOEXCEPT;

#if 201103L <= YY_CPLUSPLUS
      /// Move constructor.
      by_kind (by_kind&& that) YY_NOEXCEPT;
#endif

      /// Copy constructor.
      by_kind (const by_kind& that) YY_NOEXCEPT;

      /// Constructor from (external) token numbers.
      by_kind (kind_type t) YY_NOEXCEPT;



      /// Record that this symbol is empty.
      void clear () YY_NOEXCEPT;

      /// Steal the symbol kind from \a that.
      void move (by_kind& that);

      /// The (internal) type number (corresponding to \a type).
      /// \a empty when empty.
      symbol_kind_type kind () const YY_NOEXCEPT;

      /// Backward compatibility (Bison 3.6).
      symbol_kind_type type_get () const YY_NOEXCEPT;

      /// The symbol kind.
      /// \a S_YYEMPTY when empty.
      symbol_kind_type kind_;
    };

    /// Backward compatibility for a private implementation detail (Bison 3.6).
    typedef by_kind by_type;

    /// "External" symbols: returned by the scanner.
    struct symbol_type : basic_symbol<by_kind>
    {};

    /// Build a parser object.
    Parser (class Driver& driver_yyarg);
    virtual ~Parser ();

#if 201103L <= YY_CPLUSPLUS
    /// Non copyable.
    Parser (const Parser&) = delete;
    /// Non copyable.
    Parser& operator= (const Parser&) = delete;
#endif

    /// Parse.  An alias for parse ().
    /// \returns  0 iff parsing succeeded.
    int operator() ();

    /// Parse.
    /// \returns  0 iff parsing succeeded.
    virtual int parse ();

#if LEFPARSERDEBUG
    /// The current debugging stream.
    std::ostream& debug_stream () const YY_ATTRIBUTE_PURE;
    /// Set the current debugging stream.
    void set_debug_stream (std::ostream &);

    /// Type for debugging levels.
    typedef int debug_level_type;
    /// The current debugging level.
    debug_level_type debug_level () const YY_ATTRIBUTE_PURE;
    /// Set the current debugging level.
    void set_debug_level (debug_level_type l);
#endif

    /// Report a syntax error.
    /// \param loc    where the syntax error is found.
    /// \param msg    a description of the syntax error.
    virtual void error (const location_type& loc, const std::string& msg);

    /// Report a syntax error.
    void error (const syntax_error& err);

    /// The user-facing name of the symbol whose (internal) number is
    /// YYSYMBOL.  No bounds checking.
    static std::string symbol_name (symbol_kind_type yysymbol);



    class context
    {
    public:
      context (const Parser& yyparser, const symbol_type& yyla);
      const symbol_type& lookahead () const YY_NOEXCEPT { return yyla_; }
      symbol_kind_type token () const YY_NOEXCEPT { return yyla_.kind (); }
      const location_type& location () const YY_NOEXCEPT { return yyla_.location; }

      /// Put in YYARG at most YYARGN of the expected tokens, and return the
      /// number of tokens stored in YYARG.  If YYARG is null, return the
      /// number of expected tokens (guaranteed to be less than YYNTOKENS).
      int expected_tokens (symbol_kind_type yyarg[], int yyargn) const;

    private:
      const Parser& yyparser_;
      const symbol_type& yyla_;
    };

  private:
#if YY_CPLUSPLUS < 201103L
    /// Non copyable.
    Parser (const Parser&);
    /// Non copyable.
    Parser& operator= (const Parser&);
#endif


    /// Stored state numbers (used for stacks).
    typedef short state_type;

    /// The arguments of the error message.
    int yy_syntax_error_arguments_ (const context& yyctx,
                                    symbol_kind_type yyarg[], int yyargn) const;

    /// Generate an error message.
    /// \param yyctx     the context in which the error occurred.
    virtual std::string yysyntax_error_ (const context& yyctx) const;
    /// Compute post-reduction state.
    /// \param yystate   the current state
    /// \param yysym     the nonterminal to push on the stack
    static state_type yy_lr_goto_state_ (state_type yystate, int yysym);

    /// Whether the given \c yypact_ value indicates a defaulted state.
    /// \param yyvalue   the value to check
    static bool yy_pact_value_is_default_ (int yyvalue) YY_NOEXCEPT;

    /// Whether the given \c yytable_ value indicates a syntax error.
    /// \param yyvalue   the value to check
    static bool yy_table_value_is_error_ (int yyvalue) YY_NOEXCEPT;

    static const short yypact_ninf_;
    static const short yytable_ninf_;

    /// Convert a scanner token kind \a t to a symbol kind.
    /// In theory \a t should be a token_kind_type, but character literals
    /// are valid, yet not members of the token_kind_type enum.
    static symbol_kind_type yytranslate_ (int t) YY_NOEXCEPT;

    /// Convert the symbol name \a n to a form suitable for a diagnostic.
    static std::string yytnamerr_ (const char *yystr);

    /// For a symbol, its name in clear.
    static const char* const yytname_[];


    // Tables.
    // YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
    // STATE-NUM.
    static const short yypact_[];

    // YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
    // Performed when YYTABLE does not specify something else to do.  Zero
    // means the default is an error.
    static const short yydefact_[];

    // YYPGOTO[NTERM-NUM].
    static const short yypgoto_[];

    // YYDEFGOTO[NTERM-NUM].
    static const short yydefgoto_[];

    // YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
    // positive, shift that token.  If negative, reduce the rule whose
    // number is the opposite.  If YYTABLE_NINF, syntax error.
    static const short yytable_[];

    static const short yycheck_[];

    // YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
    // state STATE-NUM.
    static const short yystos_[];

    // YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.
    static const short yyr1_[];

    // YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.
    static const signed char yyr2_[];


#if LEFPARSERDEBUG
    // YYRLINE[YYN] -- Source line where rule number YYN was defined.
    static const short yyrline_[];
    /// Report on the debug stream that the rule \a r is going to be reduced.
    virtual void yy_reduce_print_ (int r) const;
    /// Print the state stack on the debug stream.
    virtual void yy_stack_print_ () const;

    /// Debugging level.
    int yydebug_;
    /// Debug stream.
    std::ostream* yycdebug_;

    /// \brief Display a symbol kind, value and location.
    /// \param yyo    The output stream.
    /// \param yysym  The symbol.
    template <typename Base>
    void yy_print_ (std::ostream& yyo, const basic_symbol<Base>& yysym) const;
#endif

    /// \brief Reclaim the memory associated to a symbol.
    /// \param yymsg     Why this token is reclaimed.
    ///                  If null, print nothing.
    /// \param yysym     The symbol.
    template <typename Base>
    void yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const;

  private:
    /// Type access provider for state based symbols.
    struct by_state
    {
      /// Default constructor.
      by_state () YY_NOEXCEPT;

      /// The symbol kind as needed by the constructor.
      typedef state_type kind_type;

      /// Constructor.
      by_state (kind_type s) YY_NOEXCEPT;

      /// Copy constructor.
      by_state (const by_state& that) YY_NOEXCEPT;

      /// Record that this symbol is empty.
      void clear () YY_NOEXCEPT;

      /// Steal the symbol kind from \a that.
      void move (by_state& that);

      /// The symbol kind (corresponding to \a state).
      /// \a symbol_kind::S_YYEMPTY when empty.
      symbol_kind_type kind () const YY_NOEXCEPT;

      /// The state number used to denote an empty symbol.
      /// We use the initial state, as it does not have a value.
      enum { empty_state = 0 };

      /// The state.
      /// \a empty when empty.
      state_type state;
    };

    /// "Internal" symbol: element of the stack.
    struct stack_symbol_type : basic_symbol<by_state>
    {
      /// Superclass.
      typedef basic_symbol<by_state> super_type;
      /// Construct an empty symbol.
      stack_symbol_type ();
      /// Move or copy construction.
      stack_symbol_type (YY_RVREF (stack_symbol_type) that);
      /// Steal the contents from \a sym to build this.
      stack_symbol_type (state_type s, YY_MOVE_REF (symbol_type) sym);
#if YY_CPLUSPLUS < 201103L
      /// Assignment, needed by push_back by some old implementations.
      /// Moves the contents of that.
      stack_symbol_type& operator= (stack_symbol_type& that);

      /// Assignment, needed by push_back by other implementations.
      /// Needed by some other old implementations.
      stack_symbol_type& operator= (const stack_symbol_type& that);
#endif
    };

    /// A stack with random access from its top.
    template <typename T, typename S = std::vector<T> >
    class stack
    {
    public:
      // Hide our reversed order.
      typedef typename S::iterator iterator;
      typedef typename S::const_iterator const_iterator;
      typedef typename S::size_type size_type;
      typedef typename std::ptrdiff_t index_type;

      stack (size_type n = 200) YY_NOEXCEPT
        : seq_ (n)
      {}

#if 201103L <= YY_CPLUSPLUS
      /// Non copyable.
      stack (const stack&) = delete;
      /// Non copyable.
      stack& operator= (const stack&) = delete;
#endif

      /// Random access.
      ///
      /// Index 0 returns the topmost element.
      const T&
      operator[] (index_type i) const
      {
        return seq_[size_type (size () - 1 - i)];
      }

      /// Random access.
      ///
      /// Index 0 returns the topmost element.
      T&
      operator[] (index_type i)
      {
        return seq_[size_type (size () - 1 - i)];
      }

      /// Steal the contents of \a t.
      ///
      /// Close to move-semantics.
      void
      push (YY_MOVE_REF (T) t)
      {
        seq_.push_back (T ());
        operator[] (0).move (t);
      }

      /// Pop elements from the stack.
      void
      pop (std::ptrdiff_t n = 1) YY_NOEXCEPT
      {
        for (; 0 < n; --n)
          seq_.pop_back ();
      }

      /// Pop all elements from the stack.
      void
      clear () YY_NOEXCEPT
      {
        seq_.clear ();
      }

      /// Number of elements on the stack.
      index_type
      size () const YY_NOEXCEPT
      {
        return index_type (seq_.size ());
      }

      /// Iterator on top of the stack (going downwards).
      const_iterator
      begin () const YY_NOEXCEPT
      {
        return seq_.begin ();
      }

      /// Bottom of the stack.
      const_iterator
      end () const YY_NOEXCEPT
      {
        return seq_.end ();
      }

      /// Present a slice of the top of a stack.
      class slice
      {
      public:
        slice (const stack& stack, index_type range) YY_NOEXCEPT
          : stack_ (stack)
          , range_ (range)
        {}

        const T&
        operator[] (index_type i) const
        {
          return stack_[range_ - i];
        }

      private:
        const stack& stack_;
        index_type range_;
      };

    private:
#if YY_CPLUSPLUS < 201103L
      /// Non copyable.
      stack (const stack&);
      /// Non copyable.
      stack& operator= (const stack&);
#endif
      /// The wrapped container.
      S seq_;
    };


    /// Stack type.
    typedef stack<stack_symbol_type> stack_type;

    /// The stack.
    stack_type yystack_;

    /// Push a new state on the stack.
    /// \param m    a debug message to display
    ///             if null, no trace is output.
    /// \param sym  the symbol
    /// \warning the contents of \a s.value is stolen.
    void yypush_ (const char* m, YY_MOVE_REF (stack_symbol_type) sym);

    /// Push a new look ahead token on the state on the stack.
    /// \param m    a debug message to display
    ///             if null, no trace is output.
    /// \param s    the state
    /// \param sym  the symbol (for its value and location).
    /// \warning the contents of \a sym.value is stolen.
    void yypush_ (const char* m, state_type s, YY_MOVE_REF (symbol_type) sym);

    /// Pop \a n symbols from the stack.
    void yypop_ (int n = 1) YY_NOEXCEPT;

    /// Constants.
    enum
    {
      yylast_ = 3683,     ///< Last index in yytable_.
      yynnts_ = 448,  ///< Number of nonterminal symbols.
      yyfinal_ = 4 ///< Termination state number.
    };


    // User arguments.
    class Driver& driver;

  };


} // LefParser
#line 2100 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.h"




#endif // !YY_LEFPARSER_HOME_FANTASTIC_NA_OPENMPL_OPENMPL_BUILD_THIRDPARTY_LIMBO_LIMBO_PARSERS_LEF_BISON_LEFPARSER_H_INCLUDED
