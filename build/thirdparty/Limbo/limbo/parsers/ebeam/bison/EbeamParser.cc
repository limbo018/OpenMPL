// A Bison parser, made by GNU Bison 3.8.2.

// Skeleton implementation for Bison LALR(1) parsers in C++

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

// DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
// especially those whose name start with YY_ or yy_.  They are
// private implementation details that can be changed or removed.


// Take the name prefix into account.
#define yylex   EbeamParserlex

// First part of user prologue.
#line 4 "EbeamParser.yy"
 /*** C/C++ Declarations ***/

#include <stdio.h>
#include <string>
#include <vector>

#include "EbeamDataBase.h"


#line 53 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"


#include "EbeamParser.h"

// Second part of user prologue.
#line 108 "EbeamParser.yy"


#include "EbeamDriver.h"
#include "EbeamScanner.h"

/* this "connects" the bison parser in the driver to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the driver context. */
#undef yylex
#define yylex driver.lexer->lex


#line 72 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"



#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> // FIXME: INFRINGES ON USER NAME SPACE.
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif


// Whether we are compiled with exception support.
#ifndef YY_EXCEPTIONS
# if defined __GNUC__ && !defined __EXCEPTIONS
#  define YY_EXCEPTIONS 0
# else
#  define YY_EXCEPTIONS 1
# endif
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K].location)
/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

# ifndef YYLLOC_DEFAULT
#  define YYLLOC_DEFAULT(Current, Rhs, N)                               \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).begin  = YYRHSLOC (Rhs, 1).begin;                   \
          (Current).end    = YYRHSLOC (Rhs, N).end;                     \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).begin = (Current).end = YYRHSLOC (Rhs, 0).end;      \
        }                                                               \
    while (false)
# endif


// Enable debugging if requested.
#if EBEAMPARSERDEBUG

// A pseudo ostream that takes yydebug_ into account.
# define YYCDEBUG if (yydebug_) (*yycdebug_)

# define YY_SYMBOL_PRINT(Title, Symbol)         \
  do {                                          \
    if (yydebug_)                               \
    {                                           \
      *yycdebug_ << Title << ' ';               \
      yy_print_ (*yycdebug_, Symbol);           \
      *yycdebug_ << '\n';                       \
    }                                           \
  } while (false)

# define YY_REDUCE_PRINT(Rule)          \
  do {                                  \
    if (yydebug_)                       \
      yy_reduce_print_ (Rule);          \
  } while (false)

# define YY_STACK_PRINT()               \
  do {                                  \
    if (yydebug_)                       \
      yy_stack_print_ ();                \
  } while (false)

#else // !EBEAMPARSERDEBUG

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Symbol)  YY_USE (Symbol)
# define YY_REDUCE_PRINT(Rule)           static_cast<void> (0)
# define YY_STACK_PRINT()                static_cast<void> (0)

#endif // !EBEAMPARSERDEBUG

#define yyerrok         (yyerrstatus_ = 0)
#define yyclearin       (yyla.clear ())

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)

namespace EbeamParser {
#line 165 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"

  /// Build a parser object.
  Parser::Parser (class Driver& driver_yyarg)
#if EBEAMPARSERDEBUG
    : yydebug_ (false),
      yycdebug_ (&std::cerr),
#else
    :
#endif
      driver (driver_yyarg)
  {}

  Parser::~Parser ()
  {}

  Parser::syntax_error::~syntax_error () YY_NOEXCEPT YY_NOTHROW
  {}

  /*---------.
  | symbol.  |
  `---------*/

  // basic_symbol.
  template <typename Base>
  Parser::basic_symbol<Base>::basic_symbol (const basic_symbol& that)
    : Base (that)
    , value (that.value)
    , location (that.location)
  {}


  /// Constructor for valueless symbols.
  template <typename Base>
  Parser::basic_symbol<Base>::basic_symbol (typename Base::kind_type t, YY_MOVE_REF (location_type) l)
    : Base (t)
    , value ()
    , location (l)
  {}

  template <typename Base>
  Parser::basic_symbol<Base>::basic_symbol (typename Base::kind_type t, YY_RVREF (value_type) v, YY_RVREF (location_type) l)
    : Base (t)
    , value (YY_MOVE (v))
    , location (YY_MOVE (l))
  {}


  template <typename Base>
  Parser::symbol_kind_type
  Parser::basic_symbol<Base>::type_get () const YY_NOEXCEPT
  {
    return this->kind ();
  }


  template <typename Base>
  bool
  Parser::basic_symbol<Base>::empty () const YY_NOEXCEPT
  {
    return this->kind () == symbol_kind::S_YYEMPTY;
  }

  template <typename Base>
  void
  Parser::basic_symbol<Base>::move (basic_symbol& s)
  {
    super_type::move (s);
    value = YY_MOVE (s.value);
    location = YY_MOVE (s.location);
  }

  // by_kind.
  Parser::by_kind::by_kind () YY_NOEXCEPT
    : kind_ (symbol_kind::S_YYEMPTY)
  {}

#if 201103L <= YY_CPLUSPLUS
  Parser::by_kind::by_kind (by_kind&& that) YY_NOEXCEPT
    : kind_ (that.kind_)
  {
    that.clear ();
  }
#endif

  Parser::by_kind::by_kind (const by_kind& that) YY_NOEXCEPT
    : kind_ (that.kind_)
  {}

  Parser::by_kind::by_kind (token_kind_type t) YY_NOEXCEPT
    : kind_ (yytranslate_ (t))
  {}



  void
  Parser::by_kind::clear () YY_NOEXCEPT
  {
    kind_ = symbol_kind::S_YYEMPTY;
  }

  void
  Parser::by_kind::move (by_kind& that)
  {
    kind_ = that.kind_;
    that.clear ();
  }

  Parser::symbol_kind_type
  Parser::by_kind::kind () const YY_NOEXCEPT
  {
    return kind_;
  }


  Parser::symbol_kind_type
  Parser::by_kind::type_get () const YY_NOEXCEPT
  {
    return this->kind ();
  }



  // by_state.
  Parser::by_state::by_state () YY_NOEXCEPT
    : state (empty_state)
  {}

  Parser::by_state::by_state (const by_state& that) YY_NOEXCEPT
    : state (that.state)
  {}

  void
  Parser::by_state::clear () YY_NOEXCEPT
  {
    state = empty_state;
  }

  void
  Parser::by_state::move (by_state& that)
  {
    state = that.state;
    that.clear ();
  }

  Parser::by_state::by_state (state_type s) YY_NOEXCEPT
    : state (s)
  {}

  Parser::symbol_kind_type
  Parser::by_state::kind () const YY_NOEXCEPT
  {
    if (state == empty_state)
      return symbol_kind::S_YYEMPTY;
    else
      return YY_CAST (symbol_kind_type, yystos_[+state]);
  }

  Parser::stack_symbol_type::stack_symbol_type ()
  {}

  Parser::stack_symbol_type::stack_symbol_type (YY_RVREF (stack_symbol_type) that)
    : super_type (YY_MOVE (that.state), YY_MOVE (that.value), YY_MOVE (that.location))
  {
#if 201103L <= YY_CPLUSPLUS
    // that is emptied.
    that.state = empty_state;
#endif
  }

  Parser::stack_symbol_type::stack_symbol_type (state_type s, YY_MOVE_REF (symbol_type) that)
    : super_type (s, YY_MOVE (that.value), YY_MOVE (that.location))
  {
    // that is emptied.
    that.kind_ = symbol_kind::S_YYEMPTY;
  }

#if YY_CPLUSPLUS < 201103L
  Parser::stack_symbol_type&
  Parser::stack_symbol_type::operator= (const stack_symbol_type& that)
  {
    state = that.state;
    value = that.value;
    location = that.location;
    return *this;
  }

  Parser::stack_symbol_type&
  Parser::stack_symbol_type::operator= (stack_symbol_type& that)
  {
    state = that.state;
    value = that.value;
    location = that.location;
    // that is emptied.
    that.state = empty_state;
    return *this;
  }
#endif

  template <typename Base>
  void
  Parser::yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const
  {
    if (yymsg)
      YY_SYMBOL_PRINT (yymsg, yysym);

    // User destructor.
    switch (yysym.kind ())
    {
      case symbol_kind::S_STRING: // "string"
#line 99 "EbeamParser.yy"
                    { delete (yysym.value.stringVal); }
#line 377 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
        break;

      case symbol_kind::S_QUOTE: // "quoted chars"
#line 99 "EbeamParser.yy"
                    { delete (yysym.value.quoteVal); }
#line 383 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
        break;

      case symbol_kind::S_integer_array: // integer_array
#line 100 "EbeamParser.yy"
                    { delete (yysym.value.integerArrayVal); }
#line 389 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
        break;

      case symbol_kind::S_string_array: // string_array
#line 100 "EbeamParser.yy"
                    { delete (yysym.value.stringArrayVal); }
#line 395 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
        break;

      default:
        break;
    }
  }

#if EBEAMPARSERDEBUG
  template <typename Base>
  void
  Parser::yy_print_ (std::ostream& yyo, const basic_symbol<Base>& yysym) const
  {
    std::ostream& yyoutput = yyo;
    YY_USE (yyoutput);
    if (yysym.empty ())
      yyo << "empty symbol";
    else
      {
        symbol_kind_type yykind = yysym.kind ();
        yyo << (yykind < YYNTOKENS ? "token" : "nterm")
            << ' ' << yysym.name () << " ("
            << yysym.location << ": ";
        YY_USE (yykind);
        yyo << ')';
      }
  }
#endif

  void
  Parser::yypush_ (const char* m, YY_MOVE_REF (stack_symbol_type) sym)
  {
    if (m)
      YY_SYMBOL_PRINT (m, sym);
    yystack_.push (YY_MOVE (sym));
  }

  void
  Parser::yypush_ (const char* m, state_type s, YY_MOVE_REF (symbol_type) sym)
  {
#if 201103L <= YY_CPLUSPLUS
    yypush_ (m, stack_symbol_type (s, std::move (sym)));
#else
    stack_symbol_type ss (s, sym);
    yypush_ (m, ss);
#endif
  }

  void
  Parser::yypop_ (int n) YY_NOEXCEPT
  {
    yystack_.pop (n);
  }

#if EBEAMPARSERDEBUG
  std::ostream&
  Parser::debug_stream () const
  {
    return *yycdebug_;
  }

  void
  Parser::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


  Parser::debug_level_type
  Parser::debug_level () const
  {
    return yydebug_;
  }

  void
  Parser::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }
#endif // EBEAMPARSERDEBUG

  Parser::state_type
  Parser::yy_lr_goto_state_ (state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - YYNTOKENS] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - YYNTOKENS];
  }

  bool
  Parser::yy_pact_value_is_default_ (int yyvalue) YY_NOEXCEPT
  {
    return yyvalue == yypact_ninf_;
  }

  bool
  Parser::yy_table_value_is_error_ (int yyvalue) YY_NOEXCEPT
  {
    return yyvalue == yytable_ninf_;
  }

  int
  Parser::operator() ()
  {
    return parse ();
  }

  int
  Parser::parse ()
  {
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// The lookahead symbol.
    symbol_type yyla;

    /// The locations where the error started and ended.
    stack_symbol_type yyerror_range[3];

    /// The return value of parse ().
    int yyresult;

#if YY_EXCEPTIONS
    try
#endif // YY_EXCEPTIONS
      {
    YYCDEBUG << "Starting parse\n";


    // User initialization code.
#line 41 "EbeamParser.yy"
{
    // initialize the initial location object
    yyla.location.begin.filename = yyla.location.end.filename = &driver.streamname;
}

#line 538 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"


    /* Initialize the stack.  The initial state will be set in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystack_.clear ();
    yypush_ (YY_NULLPTR, 0, YY_MOVE (yyla));

  /*-----------------------------------------------.
  | yynewstate -- push a new symbol on the stack.  |
  `-----------------------------------------------*/
  yynewstate:
    YYCDEBUG << "Entering state " << int (yystack_[0].state) << '\n';
    YY_STACK_PRINT ();

    // Accept?
    if (yystack_[0].state == yyfinal_)
      YYACCEPT;

    goto yybackup;


  /*-----------.
  | yybackup.  |
  `-----------*/
  yybackup:
    // Try to take a decision without lookahead.
    yyn = yypact_[+yystack_[0].state];
    if (yy_pact_value_is_default_ (yyn))
      goto yydefault;

    // Read a lookahead token.
    if (yyla.empty ())
      {
        YYCDEBUG << "Reading a token\n";
#if YY_EXCEPTIONS
        try
#endif // YY_EXCEPTIONS
          {
            yyla.kind_ = yytranslate_ (yylex (&yyla.value, &yyla.location));
          }
#if YY_EXCEPTIONS
        catch (const syntax_error& yyexc)
          {
            YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
            error (yyexc);
            goto yyerrlab1;
          }
#endif // YY_EXCEPTIONS
      }
    YY_SYMBOL_PRINT ("Next token is", yyla);

    if (yyla.kind () == symbol_kind::S_YYerror)
    {
      // The scanner already issued an error message, process directly
      // to error recovery.  But do not keep the error token as
      // lookahead, it is too special and may lead us to an endless
      // loop in error recovery. */
      yyla.kind_ = symbol_kind::S_YYUNDEF;
      goto yyerrlab1;
    }

    /* If the proper action on seeing token YYLA.TYPE is to reduce or
       to detect an error, take that action.  */
    yyn += yyla.kind ();
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.kind ())
      {
        goto yydefault;
      }

    // Reduce or error.
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
        if (yy_table_value_is_error_ (yyn))
          goto yyerrlab;
        yyn = -yyn;
        goto yyreduce;
      }

    // Count tokens shifted since error; after three, turn off error status.
    if (yyerrstatus_)
      --yyerrstatus_;

    // Shift the lookahead token.
    yypush_ ("Shifting", state_type (yyn), YY_MOVE (yyla));
    goto yynewstate;


  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[+yystack_[0].state];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;


  /*-----------------------------.
  | yyreduce -- do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    {
      stack_symbol_type yylhs;
      yylhs.state = yy_lr_goto_state_ (yystack_[yylen].state, yyr1_[yyn]);
      /* If YYLEN is nonzero, implement the default value of the
         action: '$$ = $1'.  Otherwise, use the top of the stack.

         Otherwise, the following line sets YYLHS.VALUE to garbage.
         This behavior is undocumented and Bison users should not rely
         upon it.  */
      if (yylen)
        yylhs.value = yystack_[yylen - 1].value;
      else
        yylhs.value = yystack_[0].value;

      // Default location.
      {
        stack_type::slice range (yystack_, yylen);
        YYLLOC_DEFAULT (yylhs.location, range, yylen);
        yyerror_range[1].location = yylhs.location;
      }

      // Perform the reduction.
      YY_REDUCE_PRINT (yyn);
#if YY_EXCEPTIONS
      try
#endif // YY_EXCEPTIONS
        {
          switch (yyn)
            {
  case 2: // number: "integer"
#line 126 "EbeamParser.yy"
                 {
		(yylhs.value.doubleVal) = (yystack_[0].value.integerVal);
	   }
#line 678 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 3: // number: "double"
#line 129 "EbeamParser.yy"
                    {
		(yylhs.value.doubleVal) = (yystack_[0].value.doubleVal);
	   }
#line 686 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 4: // integer_array: "integer"
#line 132 "EbeamParser.yy"
                        {
				(yylhs.value.integerArrayVal) = new IntegerArray(1, (yystack_[0].value.integerVal));
			  }
#line 694 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 5: // integer_array: integer_array "integer"
#line 135 "EbeamParser.yy"
                                                  {
				(yystack_[1].value.integerArrayVal)->push_back((yystack_[0].value.integerVal));
				(yylhs.value.integerArrayVal) = (yystack_[1].value.integerArrayVal);
			  }
#line 703 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 6: // string_array: "string"
#line 139 "EbeamParser.yy"
                      {
				(yylhs.value.stringArrayVal) = new StringArray(1, *(yystack_[0].value.stringVal)); delete (yystack_[0].value.stringVal);
			  }
#line 711 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 7: // string_array: string_array "string"
#line 142 "EbeamParser.yy"
                                                {
				(yystack_[1].value.stringArrayVal)->push_back(*(yystack_[0].value.stringVal)); delete (yystack_[0].value.stringVal);
				(yylhs.value.stringArrayVal) = (yystack_[1].value.stringArrayVal);
			  }
#line 720 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 8: // block_units: "UNITS" "DATABASE" "MICRONS" "integer" ';' "END" "UNITS"
#line 150 "EbeamParser.yy"
                                     {
			driver.units_cbk((yystack_[3].value.integerVal));
		   }
#line 728 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 9: // entry_offset: "OFFSET" number ';'
#line 155 "EbeamParser.yy"
                                     {
			 driver.ebeam_cbk_offset((yystack_[1].value.doubleVal));
			 }
#line 736 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 10: // entry_width: "WIDTH" number ';'
#line 158 "EbeamParser.yy"
                                   {
			 driver.ebeam_cbk_width((yystack_[1].value.doubleVal));
			}
#line 744 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 11: // entry_step: "STEP" number ';'
#line 161 "EbeamParser.yy"
                                 {
			 driver.ebeam_cbk_step((yystack_[1].value.doubleVal));
		   }
#line 752 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 12: // entry_layer: "LAYERID" integer_array ';'
#line 164 "EbeamParser.yy"
                                            {
				driver.ebeam_cbk_layerid(*(yystack_[1].value.integerArrayVal));
                delete (yystack_[1].value.integerArrayVal);
			}
#line 761 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 13: // entry_layer: "LAYER" string_array ';'
#line 168 "EbeamParser.yy"
                                                     {
				driver.ebeam_cbk_layer(*(yystack_[1].value.stringArrayVal));
			}
#line 769 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 14: // block_ebeam: "EBEAMBOUNDARY" entry_offset entry_width entry_step entry_layer "END" "EBEAMBOUNDARY"
#line 176 "EbeamParser.yy"
                                                  {
				driver.ebeam_cbk();
			}
#line 777 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 15: // begin_macro: "MACRO" "string"
#line 181 "EbeamParser.yy"
                               {
				driver.macro_cbk_name(*(yystack_[0].value.stringVal));
                delete (yystack_[0].value.stringVal);
			}
#line 786 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 16: // end_macro: "END" "string"
#line 185 "EbeamParser.yy"
                           {
			driver.macro_cbk(*(yystack_[0].value.stringVal));
            delete (yystack_[0].value.stringVal);
		  }
#line 795 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 17: // single_confsite: "CONFLICTSITE" "string" "LAYERID" "integer" ';' "SITE" integer_array ';' "END" "string"
#line 192 "EbeamParser.yy"
                                          {
				driver.macro_cbk_confsite_name(*(yystack_[8].value.stringVal));
				driver.macro_cbk_confsite_layerid((yystack_[6].value.integerVal));
				driver.macro_cbk_confsite_site(*(yystack_[3].value.integerArrayVal));
				driver.macro_cbk_confsite(*(yystack_[0].value.stringVal));
                delete (yystack_[8].value.stringVal); delete (yystack_[3].value.integerArrayVal); delete (yystack_[0].value.stringVal);
			   }
#line 807 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;

  case 18: // single_confsite: "CONFLICTSITE" "string" "LAYER" "string" ';' "SITE" integer_array ';' "END" "string"
#line 202 "EbeamParser.yy"
                                          {
				driver.macro_cbk_confsite_name(*(yystack_[8].value.stringVal));
				driver.macro_cbk_confsite_layer(*(yystack_[6].value.stringVal));
				driver.macro_cbk_confsite_site(*(yystack_[3].value.integerArrayVal));
				driver.macro_cbk_confsite(*(yystack_[0].value.stringVal));
                delete (yystack_[8].value.stringVal); delete (yystack_[6].value.stringVal); delete (yystack_[3].value.integerArrayVal); delete (yystack_[0].value.stringVal);
			   }
#line 819 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"
    break;


#line 823 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"

            default:
              break;
            }
        }
#if YY_EXCEPTIONS
      catch (const syntax_error& yyexc)
        {
          YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
          error (yyexc);
          YYERROR;
        }
#endif // YY_EXCEPTIONS
      YY_SYMBOL_PRINT ("-> $$ =", yylhs);
      yypop_ (yylen);
      yylen = 0;

      // Shift the result of the reduction.
      yypush_ (YY_NULLPTR, YY_MOVE (yylhs));
    }
    goto yynewstate;


  /*--------------------------------------.
  | yyerrlab -- here on detecting error.  |
  `--------------------------------------*/
  yyerrlab:
    // If not already recovering from an error, report this error.
    if (!yyerrstatus_)
      {
        ++yynerrs_;
        context yyctx (*this, yyla);
        std::string msg = yysyntax_error_ (yyctx);
        error (yyla.location, YY_MOVE (msg));
      }


    yyerror_range[1].location = yyla.location;
    if (yyerrstatus_ == 3)
      {
        /* If just tried and failed to reuse lookahead token after an
           error, discard it.  */

        // Return failure if at end of input.
        if (yyla.kind () == symbol_kind::S_YYEOF)
          YYABORT;
        else if (!yyla.empty ())
          {
            yy_destroy_ ("Error: discarding", yyla);
            yyla.clear ();
          }
      }

    // Else will try to reuse lookahead token after shifting the error token.
    goto yyerrlab1;


  /*---------------------------------------------------.
  | yyerrorlab -- error raised explicitly by YYERROR.  |
  `---------------------------------------------------*/
  yyerrorlab:
    /* Pacify compilers when the user code never invokes YYERROR and
       the label yyerrorlab therefore never appears in user code.  */
    if (false)
      YYERROR;

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    YY_STACK_PRINT ();
    goto yyerrlab1;


  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;   // Each real token shifted decrements this.
    // Pop stack until we find a state that shifts the error token.
    for (;;)
      {
        yyn = yypact_[+yystack_[0].state];
        if (!yy_pact_value_is_default_ (yyn))
          {
            yyn += symbol_kind::S_YYerror;
            if (0 <= yyn && yyn <= yylast_
                && yycheck_[yyn] == symbol_kind::S_YYerror)
              {
                yyn = yytable_[yyn];
                if (0 < yyn)
                  break;
              }
          }

        // Pop the current state because it cannot handle the error token.
        if (yystack_.size () == 1)
          YYABORT;

        yyerror_range[1].location = yystack_[0].location;
        yy_destroy_ ("Error: popping", yystack_[0]);
        yypop_ ();
        YY_STACK_PRINT ();
      }
    {
      stack_symbol_type error_token;

      yyerror_range[2].location = yyla.location;
      YYLLOC_DEFAULT (error_token.location, yyerror_range, 2);

      // Shift the error token.
      error_token.state = state_type (yyn);
      yypush_ ("Shifting", YY_MOVE (error_token));
    }
    goto yynewstate;


  /*-------------------------------------.
  | yyacceptlab -- YYACCEPT comes here.  |
  `-------------------------------------*/
  yyacceptlab:
    yyresult = 0;
    goto yyreturn;


  /*-----------------------------------.
  | yyabortlab -- YYABORT comes here.  |
  `-----------------------------------*/
  yyabortlab:
    yyresult = 1;
    goto yyreturn;


  /*-----------------------------------------------------.
  | yyreturn -- parsing is finished, return the result.  |
  `-----------------------------------------------------*/
  yyreturn:
    if (!yyla.empty ())
      yy_destroy_ ("Cleanup: discarding lookahead", yyla);

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    YY_STACK_PRINT ();
    while (1 < yystack_.size ())
      {
        yy_destroy_ ("Cleanup: popping", yystack_[0]);
        yypop_ ();
      }

    return yyresult;
  }
#if YY_EXCEPTIONS
    catch (...)
      {
        YYCDEBUG << "Exception caught: cleaning lookahead and stack\n";
        // Do not try to display the values of the reclaimed symbols,
        // as their printers might throw an exception.
        if (!yyla.empty ())
          yy_destroy_ (YY_NULLPTR, yyla);

        while (1 < yystack_.size ())
          {
            yy_destroy_ (YY_NULLPTR, yystack_[0]);
            yypop_ ();
          }
        throw;
      }
#endif // YY_EXCEPTIONS
  }

  void
  Parser::error (const syntax_error& yyexc)
  {
    error (yyexc.location, yyexc.what ());
  }

  /* Return YYSTR after stripping away unnecessary quotes and
     backslashes, so that it's suitable for yyerror.  The heuristic is
     that double-quoting is unnecessary unless the string contains an
     apostrophe, a comma, or backslash (other than backslash-backslash).
     YYSTR is taken from yytname.  */
  std::string
  Parser::yytnamerr_ (const char *yystr)
  {
    if (*yystr == '"')
      {
        std::string yyr;
        char const *yyp = yystr;

        for (;;)
          switch (*++yyp)
            {
            case '\'':
            case ',':
              goto do_not_strip_quotes;

            case '\\':
              if (*++yyp != '\\')
                goto do_not_strip_quotes;
              else
                goto append;

            append:
            default:
              yyr += *yyp;
              break;

            case '"':
              return yyr;
            }
      do_not_strip_quotes: ;
      }

    return yystr;
  }

  std::string
  Parser::symbol_name (symbol_kind_type yysymbol)
  {
    return yytnamerr_ (yytname_[yysymbol]);
  }



  // Parser::context.
  Parser::context::context (const Parser& yyparser, const symbol_type& yyla)
    : yyparser_ (yyparser)
    , yyla_ (yyla)
  {}

  int
  Parser::context::expected_tokens (symbol_kind_type yyarg[], int yyargn) const
  {
    // Actual number of expected tokens
    int yycount = 0;

    const int yyn = yypact_[+yyparser_.yystack_[0].state];
    if (!yy_pact_value_is_default_ (yyn))
      {
        /* Start YYX at -YYN if negative to avoid negative indexes in
           YYCHECK.  In other words, skip the first -YYN actions for
           this state because they are default actions.  */
        const int yyxbegin = yyn < 0 ? -yyn : 0;
        // Stay within bounds of both yycheck and yytname.
        const int yychecklim = yylast_ - yyn + 1;
        const int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
        for (int yyx = yyxbegin; yyx < yyxend; ++yyx)
          if (yycheck_[yyx + yyn] == yyx && yyx != symbol_kind::S_YYerror
              && !yy_table_value_is_error_ (yytable_[yyx + yyn]))
            {
              if (!yyarg)
                ++yycount;
              else if (yycount == yyargn)
                return 0;
              else
                yyarg[yycount++] = YY_CAST (symbol_kind_type, yyx);
            }
      }

    if (yyarg && yycount == 0 && 0 < yyargn)
      yyarg[0] = symbol_kind::S_YYEMPTY;
    return yycount;
  }






  int
  Parser::yy_syntax_error_arguments_ (const context& yyctx,
                                                 symbol_kind_type yyarg[], int yyargn) const
  {
    /* There are many possibilities here to consider:
       - If this state is a consistent state with a default action, then
         the only way this function was invoked is if the default action
         is an error action.  In that case, don't check for expected
         tokens because there are none.
       - The only way there can be no lookahead present (in yyla) is
         if this state is a consistent state with a default action.
         Thus, detecting the absence of a lookahead is sufficient to
         determine that there is no unexpected or expected token to
         report.  In that case, just report a simple "syntax error".
       - Don't assume there isn't a lookahead just because this state is
         a consistent state with a default action.  There might have
         been a previous inconsistent state, consistent state with a
         non-default action, or user semantic action that manipulated
         yyla.  (However, yyla is currently not documented for users.)
       - Of course, the expected token list depends on states to have
         correct lookahead information, and it depends on the parser not
         to perform extra reductions after fetching a lookahead from the
         scanner and before detecting a syntax error.  Thus, state merging
         (from LALR or IELR) and default reductions corrupt the expected
         token list.  However, the list is correct for canonical LR with
         one exception: it will still contain any token that will not be
         accepted due to an error action in a later state.
    */

    if (!yyctx.lookahead ().empty ())
      {
        if (yyarg)
          yyarg[0] = yyctx.token ();
        int yyn = yyctx.expected_tokens (yyarg ? yyarg + 1 : yyarg, yyargn - 1);
        return yyn + 1;
      }
    return 0;
  }

  // Generate an error message.
  std::string
  Parser::yysyntax_error_ (const context& yyctx) const
  {
    // Its maximum.
    enum { YYARGS_MAX = 5 };
    // Arguments of yyformat.
    symbol_kind_type yyarg[YYARGS_MAX];
    int yycount = yy_syntax_error_arguments_ (yyctx, yyarg, YYARGS_MAX);

    char const* yyformat = YY_NULLPTR;
    switch (yycount)
      {
#define YYCASE_(N, S)                         \
        case N:                               \
          yyformat = S;                       \
        break
      default: // Avoid compiler warnings.
        YYCASE_ (0, YY_("syntax error"));
        YYCASE_ (1, YY_("syntax error, unexpected %s"));
        YYCASE_ (2, YY_("syntax error, unexpected %s, expecting %s"));
        YYCASE_ (3, YY_("syntax error, unexpected %s, expecting %s or %s"));
        YYCASE_ (4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
        YYCASE_ (5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
#undef YYCASE_
      }

    std::string yyres;
    // Argument number.
    std::ptrdiff_t yyi = 0;
    for (char const* yyp = yyformat; *yyp; ++yyp)
      if (yyp[0] == '%' && yyp[1] == 's' && yyi < yycount)
        {
          yyres += symbol_name (yyarg[yyi++]);
          ++yyp;
        }
      else
        yyres += *yyp;
    return yyres;
  }


  const signed char Parser::yypact_ninf_ = -32;

  const signed char Parser::yytable_ninf_ = -1;

  const signed char
  Parser::yypact_[] =
  {
     -32,     5,   -32,   -17,     7,     1,    19,   -32,   -32,    -8,
     -32,   -32,     3,   -32,     6,    20,    13,   -32,    -4,    29,
     -32,   -32,    14,     6,    21,    10,    32,   -32,   -32,    16,
     -32,    17,     6,    12,    31,    35,   -32,    34,   -32,    22,
      39,    38,    40,    23,    24,    41,   -32,   -32,    -3,   -32,
       2,    30,    33,    36,   -32,   -32,   -32,   -32,   -32,   -32,
      39,    39,    -2,    -1,    43,    44,    47,    49,   -32,   -32
  };

  const signed char
  Parser::yydefact_[] =
  {
      22,     0,     1,     0,     0,     0,     0,    23,    24,     0,
      25,    26,     0,    15,     0,     0,     0,    19,     0,     0,
       2,     3,     0,     0,     0,     0,     0,    21,    20,     0,
       9,     0,     0,     0,     0,     0,    16,     0,    10,     0,
       0,     0,     0,     0,     0,     0,    11,     4,     0,     6,
       0,     0,     0,     0,     8,     5,    12,     7,    13,    14,
       0,     0,     0,     0,     0,     0,     0,     0,    17,    18
  };

  const signed char
  Parser::yypgoto_[] =
  {
     -32,   -15,   -31,   -32,   -32,   -32,   -32,   -32,   -32,   -32,
     -32,   -32,    37,   -32,   -32,   -32
  };

  const signed char
  Parser::yydefgoto_[] =
  {
       0,    22,    48,    50,     7,    15,    24,    33,    42,     8,
       9,    27,    17,    18,    10,     1
  };

  const signed char
  Parser::yytable_[] =
  {
      55,    55,    55,    26,    11,     2,    13,    57,    31,    20,
      21,    16,     3,     4,    19,    16,    12,    39,    25,    56,
      64,    65,     5,     6,    58,    34,    35,    40,    41,    62,
      63,    14,    29,    23,    43,    32,    30,    36,    37,    38,
      44,    45,    47,    49,    46,    52,    53,    51,    59,    54,
      66,    67,    68,    60,    69,    28,    61
  };

  const signed char
  Parser::yycheck_[] =
  {
       3,     3,     3,     7,    21,     0,     5,     5,    23,     3,
       4,    19,     7,     8,    11,    19,     9,    32,     5,    22,
      22,    22,    17,    18,    22,    15,    16,    15,    16,    60,
      61,    12,     3,    13,     3,    14,    22,     5,    22,    22,
       5,     7,     3,     5,    22,    22,    22,     7,    18,     8,
       7,     7,     5,    20,     5,    18,    20
  };

  const signed char
  Parser::yystos_[] =
  {
       0,    38,     0,     7,     8,    17,    18,    27,    32,    33,
      37,    21,     9,     5,    12,    28,    19,    35,    36,    11,
       3,     4,    24,    13,    29,     5,     7,    34,    35,     3,
      22,    24,    14,    30,    15,    16,     5,    22,    22,    24,
      15,    16,    31,     3,     5,     7,    22,     3,    25,     5,
      26,     7,    22,    22,     8,     3,    22,     5,    22,    18,
      20,    20,    25,    25,    22,    22,     7,     7,     5,     5
  };

  const signed char
  Parser::yyr1_[] =
  {
       0,    23,    24,    24,    25,    25,    26,    26,    27,    28,
      29,    30,    31,    31,    32,    33,    34,    35,    35,    36,
      36,    37,    38,    38,    38,    38,    38
  };

  const signed char
  Parser::yyr2_[] =
  {
       0,     2,     1,     1,     1,     2,     1,     2,     7,     3,
       3,     3,     3,     3,     7,     2,     2,    10,    10,     1,
       2,     3,     0,     2,     2,     2,     3
  };


#if EBEAMPARSERDEBUG || 1
  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a YYNTOKENS, nonterminals.
  const char*
  const Parser::yytname_[] =
  {
  "\"end of file\"", "error", "\"invalid token\"", "\"integer\"",
  "\"double\"", "\"string\"", "\"quoted chars\"", "\"END\"", "\"UNITS\"",
  "\"DATABASE\"", "\"DISTANCE\"", "\"MICRONS\"", "\"OFFSET\"", "\"WIDTH\"",
  "\"STEP\"", "\"LAYERID\"", "\"LAYER\"", "\"MACRO\"", "\"EBEAMBOUNDARY\"",
  "\"CONFLICTSITE\"", "\"SITE\"", "\"LIBRARY\"", "';'", "$accept",
  "number", "integer_array", "string_array", "block_units", "entry_offset",
  "entry_width", "entry_step", "entry_layer", "block_ebeam", "begin_macro",
  "end_macro", "single_confsite", "multiple_confsites", "block_macro",
  "start", YY_NULLPTR
  };
#endif


#if EBEAMPARSERDEBUG
  const unsigned char
  Parser::yyrline_[] =
  {
       0,   126,   126,   129,   132,   135,   139,   142,   148,   155,
     158,   161,   164,   168,   171,   181,   185,   189,   199,   209,
     210,   212,   218,   219,   220,   221,   222
  };

  void
  Parser::yy_stack_print_ () const
  {
    *yycdebug_ << "Stack now";
    for (stack_type::const_iterator
           i = yystack_.begin (),
           i_end = yystack_.end ();
         i != i_end; ++i)
      *yycdebug_ << ' ' << int (i->state);
    *yycdebug_ << '\n';
  }

  void
  Parser::yy_reduce_print_ (int yyrule) const
  {
    int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
               << " (line " << yylno << "):\n";
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
                       yystack_[(yynrhs) - (yyi + 1)]);
  }
#endif // EBEAMPARSERDEBUG

  Parser::symbol_kind_type
  Parser::yytranslate_ (int t) YY_NOEXCEPT
  {
    // YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to
    // TOKEN-NUM as returned by yylex.
    static
    const signed char
    translate_table[] =
    {
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    22,
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
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21
    };
    // Last valid token kind.
    const int code_max = 276;

    if (t <= 0)
      return symbol_kind::S_YYEOF;
    else if (t <= code_max)
      return static_cast <symbol_kind_type> (translate_table[t]);
    else
      return symbol_kind::S_YYUNDEF;
  }

} // EbeamParser
#line 1373 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/ebeam/bison/EbeamParser.cc"

#line 227 "EbeamParser.yy"
 /*** Additional Code ***/

void EbeamParser::Parser::error(const Parser::location_type& l,
			    const std::string& m)
{
    driver.error(l, m);
}
