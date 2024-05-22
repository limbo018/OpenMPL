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
 ** \file /home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.h
 ** Define the DefParser::parser class.
 */

// C++ LALR(1) parser skeleton written by Akim Demaille.

// DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
// especially those whose name start with YY_ or yy_.  They are
// private implementation details that can be changed or removed.

#ifndef YY_DEFPARSER_HOME_FANTASTIC_NA_OPENMPL_OPENMPL_BUILD_THIRDPARTY_LIMBO_LIMBO_PARSERS_DEF_BISON_DEFPARSER_H_INCLUDED
# define YY_DEFPARSER_HOME_FANTASTIC_NA_OPENMPL_OPENMPL_BUILD_THIRDPARTY_LIMBO_LIMBO_PARSERS_DEF_BISON_DEFPARSER_H_INCLUDED


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
#ifndef DEFPARSERDEBUG
# if defined YYDEBUG
#if YYDEBUG
#   define DEFPARSERDEBUG 1
#  else
#   define DEFPARSERDEBUG 0
#  endif
# else /* ! defined YYDEBUG */
#  define DEFPARSERDEBUG 1
# endif /* ! defined YYDEBUG */
#endif  /* ! defined DEFPARSERDEBUG */

namespace DefParser {
#line 190 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.h"




  /// A Bison parser.
  class Parser
  {
  public:
#ifdef DEFPARSERSTYPE
# ifdef __GNUC__
#  pragma GCC message "bison: do not #define DEFPARSERSTYPE in C++, use %define api.value.type"
# endif
    typedef DEFPARSERSTYPE value_type;
#else
    /// Symbol semantic values.
    union value_type
    {
#line 56 "DefParser.yy"

    int  			integerVal;
    double 			doubleVal;
    std::string*		stringVal;
	std::string*		quoteVal;
	std::string*		binaryVal;

/*	class IntegerArray* integerArrayVal;*/
	class StringArray* stringArrayVal;

#line 219 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.h"

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
        DEFPARSEREMPTY = -2,
    END = 0,                       // "end of file"
    DEFPARSERerror = 256,          // error
    DEFPARSERUNDEF = 257,          // "invalid token"
    INTEGER = 258,                 // "integer"
    DOUBLE = 259,                  // "double"
    STRING = 260,                  // "string"
    QUOTE = 261,                   // "quoted chars"
    BINARY = 262,                  // "binary numbers"
    KWD_VERSION = 263,             // "VERSION"
    KWD_DIVIDERCHAR = 264,         // "DIVIDERCHAR"
    KWD_BUSBITCHARS = 265,         // "BUSBITCHARS"
    KWD_DESIGN = 266,              // "DESIGN"
    KWD_UNITS = 267,               // "UNITS"
    KWD_DISTANCE = 268,            // "DISTANCE"
    KWD_MICRONS = 269,             // "MICRONS"
    KWD_DIEAREA = 270,             // "DIEAREA"
    KWD_ROW = 271,                 // "ROW"
    KWD_DO = 272,                  // "DO"
    KWD_BY = 273,                  // "BY"
    KWD_COMPONENTS = 274,          // "COMPONENTS"
    KWD_PINS = 275,                // "PINS"
    KWD_NETS = 276,                // "NETS"
    KWD_NET = 277,                 // "NET"
    KWD_STEP = 278,                // "STEP"
    KWD_END = 279,                 // "END"
    KWD_DIRECTION = 280,           // "DIRECTION"
    KWD_LAYER = 281,               // "LAYER"
    KWD_PROPERTYDEFINITIONS = 282, // "PROPERTYDEFINITIONS"
    KWD_COMPONENTPIN = 283,        // "COMPONENTPIN"
    KWD_TRACKS = 284,              // "TRACKS"
    KWD_GCELLGRID = 285,           // "GCELLGRID"
    KWD_VIAS = 286,                // "VIAS"
    KWD_VIA = 287,                 // "VIA"
    KWD_VIARULE = 288,             // "VIARULE"
    KWD_CUTSIZE = 289,             // "CUTSIZE"
    KWD_LAYERS = 290,              // "LAYERS"
    KWD_ROWCOL = 291,              // "ROWCOL"
    KWD_ENCLOSURE = 292,           // "ENCLOSURE"
    KWD_CUTSPACING = 293,          // "CUTSPACING"
    KWD_USE = 294,                 // "USE"
    KWD_SPECIALNETS = 295,         // "SPECIALNETS"
    KWD_SHAPE = 296,               // "SHAPE"
    KWD_RECT = 297,                // "RECT"
    KWD_SOURCE = 298,              // "SOURCE"
    KWD_NONDEFAULTRULES = 299,     // "NONDEFAULTRULES"
    KWD_NONDEFAULTRULE = 300,      // "NONDEFAULTRULE"
    KWD_HARDSPACING = 301,         // "HARDSPACING"
    KWD_WIDTH = 302,               // "WIDTH"
    KWD_SPACING = 303,             // "SPACING"
    KWD_REGIONS = 304,             // "REGIONS"
    KWD_REGION = 305,              // "REGION"
    KWD_TYPE = 306,                // "TYPE"
    KWD_FENCE = 307,               // "FENCE"
    KWD_GROUPS = 308,              // "GROUPS"
    KWD_GROUP = 309,               // "GROUP"
    KWD_BLOCKAGES = 310,           // "BLOCKAGES"
    KWD_PLACEMENT = 311,           // "PLACEMENT"
    KWD_ROUTING = 312              // "ROUTING"
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
        YYNTOKENS = 64, ///< Number of tokens.
        S_YYEMPTY = -2,
        S_YYEOF = 0,                             // "end of file"
        S_YYerror = 1,                           // error
        S_YYUNDEF = 2,                           // "invalid token"
        S_INTEGER = 3,                           // "integer"
        S_DOUBLE = 4,                            // "double"
        S_STRING = 5,                            // "string"
        S_QUOTE = 6,                             // "quoted chars"
        S_BINARY = 7,                            // "binary numbers"
        S_KWD_VERSION = 8,                       // "VERSION"
        S_KWD_DIVIDERCHAR = 9,                   // "DIVIDERCHAR"
        S_KWD_BUSBITCHARS = 10,                  // "BUSBITCHARS"
        S_KWD_DESIGN = 11,                       // "DESIGN"
        S_KWD_UNITS = 12,                        // "UNITS"
        S_KWD_DISTANCE = 13,                     // "DISTANCE"
        S_KWD_MICRONS = 14,                      // "MICRONS"
        S_KWD_DIEAREA = 15,                      // "DIEAREA"
        S_KWD_ROW = 16,                          // "ROW"
        S_KWD_DO = 17,                           // "DO"
        S_KWD_BY = 18,                           // "BY"
        S_KWD_COMPONENTS = 19,                   // "COMPONENTS"
        S_KWD_PINS = 20,                         // "PINS"
        S_KWD_NETS = 21,                         // "NETS"
        S_KWD_NET = 22,                          // "NET"
        S_KWD_STEP = 23,                         // "STEP"
        S_KWD_END = 24,                          // "END"
        S_KWD_DIRECTION = 25,                    // "DIRECTION"
        S_KWD_LAYER = 26,                        // "LAYER"
        S_KWD_PROPERTYDEFINITIONS = 27,          // "PROPERTYDEFINITIONS"
        S_KWD_COMPONENTPIN = 28,                 // "COMPONENTPIN"
        S_KWD_TRACKS = 29,                       // "TRACKS"
        S_KWD_GCELLGRID = 30,                    // "GCELLGRID"
        S_KWD_VIAS = 31,                         // "VIAS"
        S_KWD_VIA = 32,                          // "VIA"
        S_KWD_VIARULE = 33,                      // "VIARULE"
        S_KWD_CUTSIZE = 34,                      // "CUTSIZE"
        S_KWD_LAYERS = 35,                       // "LAYERS"
        S_KWD_ROWCOL = 36,                       // "ROWCOL"
        S_KWD_ENCLOSURE = 37,                    // "ENCLOSURE"
        S_KWD_CUTSPACING = 38,                   // "CUTSPACING"
        S_KWD_USE = 39,                          // "USE"
        S_KWD_SPECIALNETS = 40,                  // "SPECIALNETS"
        S_KWD_SHAPE = 41,                        // "SHAPE"
        S_KWD_RECT = 42,                         // "RECT"
        S_KWD_SOURCE = 43,                       // "SOURCE"
        S_KWD_NONDEFAULTRULES = 44,              // "NONDEFAULTRULES"
        S_KWD_NONDEFAULTRULE = 45,               // "NONDEFAULTRULE"
        S_KWD_HARDSPACING = 46,                  // "HARDSPACING"
        S_KWD_WIDTH = 47,                        // "WIDTH"
        S_KWD_SPACING = 48,                      // "SPACING"
        S_KWD_REGIONS = 49,                      // "REGIONS"
        S_KWD_REGION = 50,                       // "REGION"
        S_KWD_TYPE = 51,                         // "TYPE"
        S_KWD_FENCE = 52,                        // "FENCE"
        S_KWD_GROUPS = 53,                       // "GROUPS"
        S_KWD_GROUP = 54,                        // "GROUP"
        S_KWD_BLOCKAGES = 55,                    // "BLOCKAGES"
        S_KWD_PLACEMENT = 56,                    // "PLACEMENT"
        S_KWD_ROUTING = 57,                      // "ROUTING"
        S_58_ = 58,                              // ';'
        S_59_ = 59,                              // '+'
        S_60_ = 60,                              // '('
        S_61_ = 61,                              // ')'
        S_62_ = 62,                              // '-'
        S_63_ = 63,                              // '*'
        S_YYACCEPT = 64,                         // $accept
        S_string_array = 65,                     // string_array
        S_block_other = 66,                      // block_other
        S_single_row = 67,                       // single_row
        S_block_rows = 68,                       // block_rows
        S_single_tracks = 69,                    // single_tracks
        S_block_tracks = 70,                     // block_tracks
        S_single_gcellgrid = 71,                 // single_gcellgrid
        S_block_gcellgrid = 72,                  // block_gcellgrid
        S_begin_vias = 73,                       // begin_vias
        S_end_vias = 74,                         // end_vias
        S_via_addon = 75,                        // via_addon
        S_single_via = 76,                       // single_via
        S_multiple_vias = 77,                    // multiple_vias
        S_block_vias = 78,                       // block_vias
        S_begin_nondefaultrules = 79,            // begin_nondefaultrules
        S_end_nondefaultrules = 80,              // end_nondefaultrules
        S_nondefaultrule_addon = 81,             // nondefaultrule_addon
        S_single_nondefaultrule = 82,            // single_nondefaultrule
        S_multiple_nondefaultrules = 83,         // multiple_nondefaultrules
        S_block_nondefaultrules = 84,            // block_nondefaultrules
        S_begin_regions = 85,                    // begin_regions
        S_end_regions = 86,                      // end_regions
        S_region_addon = 87,                     // region_addon
        S_region_points = 88,                    // region_points
        S_single_region = 89,                    // single_region
        S_multiple_regions = 90,                 // multiple_regions
        S_block_regions = 91,                    // block_regions
        S_begin_components = 92,                 // begin_components
        S_end_components = 93,                   // end_components
        S_component_addon = 94,                  // component_addon
        S_single_component = 95,                 // single_component
        S_multiple_components = 96,              // multiple_components
        S_block_components = 97,                 // block_components
        S_begin_pins = 98,                       // begin_pins
        S_end_pins = 99,                         // end_pins
        S_pin_addon = 100,                       // pin_addon
        S_single_pin = 101,                      // single_pin
        S_multiple_pins = 102,                   // multiple_pins
        S_block_pins = 103,                      // block_pins
        S_begin_blockages = 104,                 // begin_blockages
        S_end_blockages = 105,                   // end_blockages
        S_single_blockage = 106,                 // single_blockage
        S_multiple_blockages = 107,              // multiple_blockages
        S_block_blockages = 108,                 // block_blockages
        S_begin_specialnets = 109,               // begin_specialnets
        S_end_specialnets = 110,                 // end_specialnets
        S_specialnets_metal_layer = 111,         // specialnets_metal_layer
        S_specialnets_metal_shape = 112,         // specialnets_metal_shape
        S_specialnets_metal_array = 113,         // specialnets_metal_array
        S_specialnets_addon = 114,               // specialnets_addon
        S_single_specialnet = 115,               // single_specialnet
        S_multiple_specialnets = 116,            // multiple_specialnets
        S_block_specialnets = 117,               // block_specialnets
        S_begin_nets = 118,                      // begin_nets
        S_end_nets = 119,                        // end_nets
        S_node_pin_pair = 120,                   // node_pin_pair
        S_node_pin_pairs = 121,                  // node_pin_pairs
        S_net_addon = 122,                       // net_addon
        S_single_net = 123,                      // single_net
        S_multiple_nets = 124,                   // multiple_nets
        S_block_nets = 125,                      // block_nets
        S_single_propterty = 126,                // single_propterty
        S_multiple_property = 127,               // multiple_property
        S_block_propertydefinitions = 128,       // block_propertydefinitions
        S_begin_groups = 129,                    // begin_groups
        S_end_groups = 130,                      // end_groups
        S_group_addon = 131,                     // group_addon
        S_single_group = 132,                    // single_group
        S_multiple_groups = 133,                 // multiple_groups
        S_block_groups = 134,                    // block_groups
        S_begin_design = 135,                    // begin_design
        S_end_design = 136,                      // end_design
        S_block_unit = 137,                      // block_unit
        S_block_diearea = 138,                   // block_diearea
        S_block_option = 139,                    // block_option
        S_block_design = 140,                    // block_design
        S_start = 141                            // start
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

#if DEFPARSERDEBUG
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
    static const signed char yytable_ninf_;

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
    static const unsigned char yydefact_[];

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
    static const unsigned char yystos_[];

    // YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.
    static const unsigned char yyr1_[];

    // YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.
    static const signed char yyr2_[];


#if DEFPARSERDEBUG
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
      yylast_ = 528,     ///< Last index in yytable_.
      yynnts_ = 78,  ///< Number of nonterminal symbols.
      yyfinal_ = 2 ///< Termination state number.
    };


    // User arguments.
    class Driver& driver;

  };


} // DefParser
#line 995 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.h"




#endif // !YY_DEFPARSER_HOME_FANTASTIC_NA_OPENMPL_OPENMPL_BUILD_THIRDPARTY_LIMBO_LIMBO_PARSERS_DEF_BISON_DEFPARSER_H_INCLUDED
