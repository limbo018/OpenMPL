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
#define yylex   DefParserlex

// First part of user prologue.
#line 4 "DefParser.yy"
 /*** C/C++ Declarations ***/

#include <stdio.h>
#include <string>
#include <vector>

/*#include "expression.h"*/


#line 53 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"


#include "DefParser.h"

// Second part of user prologue.
#line 141 "DefParser.yy"


#include "DefDriver.h"
#include "DefScanner.h"

/* this "connects" the bison parser in the driver to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the driver context. */
#undef yylex
#define yylex driver.lexer->lex


#line 72 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"



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
#if DEFPARSERDEBUG

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

#else // !DEFPARSERDEBUG

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Symbol)  YY_USE (Symbol)
# define YY_REDUCE_PRINT(Rule)           static_cast<void> (0)
# define YY_STACK_PRINT()                static_cast<void> (0)

#endif // !DEFPARSERDEBUG

#define yyerrok         (yyerrstatus_ = 0)
#define yyclearin       (yyla.clear ())

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)

namespace DefParser {
#line 165 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"

  /// Build a parser object.
  Parser::Parser (class Driver& driver_yyarg)
#if DEFPARSERDEBUG
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
#line 132 "DefParser.yy"
                    { delete (yysym.value.stringVal); }
#line 377 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
        break;

      case symbol_kind::S_QUOTE: // "quoted chars"
#line 132 "DefParser.yy"
                    { delete (yysym.value.quoteVal); }
#line 383 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
        break;

      case symbol_kind::S_BINARY: // "binary numbers"
#line 132 "DefParser.yy"
                    { delete (yysym.value.binaryVal); }
#line 389 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
        break;

      case symbol_kind::S_string_array: // string_array
#line 133 "DefParser.yy"
                    { delete (yysym.value.stringArrayVal); }
#line 395 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
        break;

      default:
        break;
    }
  }

#if DEFPARSERDEBUG
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

#if DEFPARSERDEBUG
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
#endif // DEFPARSERDEBUG

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
#line 41 "DefParser.yy"
{
    // initialize the initial location object
    yyla.location.begin.filename = yyla.location.end.filename = &driver.streamname;
}

#line 538 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"


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
  case 2: // string_array: "string"
#line 167 "DefParser.yy"
                      {
				(yylhs.value.stringArrayVal) = new StringArray(1, *(yystack_[0].value.stringVal));
                delete (yystack_[0].value.stringVal);
			  }
#line 679 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 3: // string_array: string_array "string"
#line 171 "DefParser.yy"
                                                {
				(yystack_[1].value.stringArrayVal)->push_back(*(yystack_[0].value.stringVal));
                delete (yystack_[0].value.stringVal);
				(yylhs.value.stringArrayVal) = (yystack_[1].value.stringArrayVal);
			  }
#line 689 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 4: // block_other: "VERSION" "double" ';'
#line 177 "DefParser.yy"
                                     {
				driver.version_cbk((yystack_[1].value.doubleVal));
			}
#line 697 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 5: // block_other: "DIVIDERCHAR" "quoted chars" ';'
#line 180 "DefParser.yy"
                                                    {
				driver.dividerchar_cbk(*(yystack_[1].value.quoteVal));
                delete (yystack_[1].value.quoteVal);
			}
#line 706 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 6: // block_other: "BUSBITCHARS" "quoted chars" ';'
#line 184 "DefParser.yy"
                                                    {
				driver.busbitchars_cbk(*(yystack_[1].value.quoteVal));
                delete (yystack_[1].value.quoteVal);
			}
#line 715 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 7: // single_row: "ROW" "string" "string" "integer" "integer" "string" "DO" "integer" "BY" "integer" "STEP" "integer" "integer" ';'
#line 190 "DefParser.yy"
                                                                                                                     {
				driver.row_cbk(*(yystack_[12].value.stringVal), *(yystack_[11].value.stringVal), (yystack_[10].value.integerVal), (yystack_[9].value.integerVal), *(yystack_[8].value.stringVal), (yystack_[6].value.integerVal), (yystack_[4].value.integerVal), (yystack_[2].value.integerVal), (yystack_[1].value.integerVal));
                delete (yystack_[12].value.stringVal);
                delete (yystack_[11].value.stringVal);
                delete (yystack_[8].value.stringVal);
			 }
#line 726 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 10: // single_tracks: "TRACKS" "string" "integer" "DO" "integer" "STEP" "integer" "LAYER" string_array ';'
#line 202 "DefParser.yy"
                                                                                                     {
				driver.track_cbk(*(yystack_[8].value.stringVal), (yystack_[7].value.integerVal), (yystack_[5].value.integerVal), (yystack_[3].value.integerVal), *(yystack_[1].value.stringArrayVal));
                delete (yystack_[8].value.stringVal);
                delete (yystack_[1].value.stringArrayVal);
			 }
#line 736 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 13: // single_gcellgrid: "GCELLGRID" "string" "integer" "DO" "integer" "STEP" "integer" ';'
#line 213 "DefParser.yy"
                                                                                    {
				driver.gcellgrid_cbk(*(yystack_[6].value.stringVal), (yystack_[5].value.integerVal), (yystack_[3].value.integerVal), (yystack_[1].value.integerVal));
                delete (yystack_[6].value.stringVal);
			 }
#line 745 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 19: // via_addon: via_addon '+' "VIARULE" "string"
#line 230 "DefParser.yy"
                                                     {delete (yystack_[0].value.stringVal);}
#line 751 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 21: // via_addon: via_addon '+' "LAYERS" string_array
#line 232 "DefParser.yy"
                                                          {delete (yystack_[0].value.stringArrayVal);}
#line 757 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 24: // via_addon: via_addon '+' "RECT" "string" '(' "integer" "integer" ')' '(' "integer" "integer" ')'
#line 235 "DefParser.yy"
                                                                                                  {delete (yystack_[8].value.stringVal);}
#line 763 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 26: // single_via: '-' "string" via_addon ';'
#line 239 "DefParser.yy"
                                      {delete (yystack_[2].value.stringVal);}
#line 769 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 35: // nondefaultrule_addon: nondefaultrule_addon '+' "LAYER" "string" "WIDTH" "integer" "SPACING" "integer"
#line 258 "DefParser.yy"
                                                                                             {delete (yystack_[4].value.stringVal);}
#line 775 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 36: // nondefaultrule_addon: nondefaultrule_addon '+' "VIA" "string"
#line 259 "DefParser.yy"
                                                     {delete (yystack_[0].value.stringVal);}
#line 781 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 37: // single_nondefaultrule: '-' "string" nondefaultrule_addon ';'
#line 262 "DefParser.yy"
                                                            {delete (yystack_[2].value.stringVal);}
#line 787 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 48: // single_region: '-' "string" region_points region_addon ';'
#line 287 "DefParser.yy"
                                                          {delete (yystack_[3].value.stringVal);}
#line 793 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 53: // begin_components: "COMPONENTS" "integer" ';'
#line 299 "DefParser.yy"
                                              {
					driver.component_cbk_size((yystack_[1].value.integerVal));
				 }
#line 801 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 56: // component_addon: component_addon '+' "string" '(' "integer" "integer" ')' "string"
#line 307 "DefParser.yy"
                                                                                            {
					driver.component_cbk_position(*(yystack_[5].value.stringVal), (yystack_[3].value.integerVal), (yystack_[2].value.integerVal), *(yystack_[0].value.stringVal));
                    delete (yystack_[5].value.stringVal);
                    delete (yystack_[0].value.stringVal);
				}
#line 811 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 57: // component_addon: component_addon '+' "string" '(' "double" "double" ')' "string"
#line 312 "DefParser.yy"
                                                                                          { /*it may be double in some benchmarks*/
					driver.component_cbk_position(*(yystack_[5].value.stringVal), (yystack_[3].value.doubleVal), (yystack_[2].value.doubleVal), *(yystack_[0].value.stringVal));
                    delete (yystack_[5].value.stringVal); 
                    delete (yystack_[0].value.stringVal);
				}
#line 821 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 58: // component_addon: component_addon '+' "SOURCE" "string"
#line 317 "DefParser.yy"
                                                                        {
					driver.component_cbk_source(*(yystack_[0].value.stringVal));
                    delete (yystack_[0].value.stringVal);
				}
#line 830 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 59: // component_addon: component_addon '+' "string"
#line 321 "DefParser.yy"
                                                             {
					driver.component_cbk_position(*(yystack_[0].value.stringVal));
                    delete (yystack_[0].value.stringVal);
				}
#line 839 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 60: // single_component: '-' "string" "string" component_addon ';'
#line 327 "DefParser.yy"
                                                         {
				driver.component_cbk(*(yystack_[3].value.stringVal), *(yystack_[2].value.stringVal));
                delete (yystack_[3].value.stringVal);
                delete (yystack_[2].value.stringVal);
			}
#line 849 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 65: // begin_pins: "PINS" "integer" ';'
#line 345 "DefParser.yy"
                                  {
					driver.pin_cbk_size((yystack_[1].value.integerVal));
				}
#line 857 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 68: // pin_addon: pin_addon '+' "NET" "string"
#line 353 "DefParser.yy"
                                                 {
			driver.pin_cbk_net(*(yystack_[0].value.stringVal));
            delete (yystack_[0].value.stringVal);
		  }
#line 866 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 69: // pin_addon: pin_addon '+' "DIRECTION" "string"
#line 357 "DefParser.yy"
                                                       {
			driver.pin_cbk_direction(*(yystack_[0].value.stringVal));
            delete (yystack_[0].value.stringVal);
		  }
#line 875 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 70: // pin_addon: pin_addon '+' "string" '(' "integer" "integer" ')' "string"
#line 361 "DefParser.yy"
                                                                        {
			driver.pin_cbk_position(*(yystack_[5].value.stringVal), (yystack_[3].value.integerVal), (yystack_[2].value.integerVal), *(yystack_[0].value.stringVal));
            delete (yystack_[5].value.stringVal);
            delete (yystack_[0].value.stringVal);
		  }
#line 885 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 71: // pin_addon: pin_addon '+' "LAYER" "string" '(' "integer" "integer" ')' '(' "integer" "integer" ')'
#line 366 "DefParser.yy"
                                                                                                   {
			driver.pin_cbk_bbox(*(yystack_[8].value.stringVal), (yystack_[6].value.integerVal), (yystack_[5].value.integerVal), (yystack_[2].value.integerVal), (yystack_[1].value.integerVal));
            delete (yystack_[8].value.stringVal);
		  }
#line 894 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 72: // pin_addon: pin_addon '+' "USE" "string"
#line 370 "DefParser.yy"
                                                 {
			driver.pin_cbk_use(*(yystack_[0].value.stringVal));
            delete (yystack_[0].value.stringVal);
		  }
#line 903 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 73: // single_pin: '-' "string" pin_addon ';'
#line 376 "DefParser.yy"
                                      {
			driver.pin_cbk(*(yystack_[2].value.stringVal));
            delete (yystack_[2].value.stringVal);
		   }
#line 912 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 78: // begin_blockages: "BLOCKAGES" "integer" ';'
#line 393 "DefParser.yy"
                                            {driver.blockage_cbk_size((yystack_[1].value.integerVal));}
#line 918 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 80: // single_blockage: '-' "PLACEMENT" "RECT" '(' "integer" "integer" ')' '(' "integer" "integer" ')' ';'
#line 397 "DefParser.yy"
                                                                                                 {driver.blockage_cbk_placement((yystack_[7].value.integerVal), (yystack_[6].value.integerVal), (yystack_[3].value.integerVal), (yystack_[2].value.integerVal));}
#line 924 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 81: // single_blockage: '-' "ROUTING" "RECT" '(' "integer" "integer" ')' '(' "integer" "integer" ')' ';'
#line 398 "DefParser.yy"
                                                                                               {driver.blockage_cbk_routing((yystack_[7].value.integerVal), (yystack_[6].value.integerVal), (yystack_[3].value.integerVal), (yystack_[2].value.integerVal));}
#line 930 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 82: // single_blockage: '-' "LAYER" "string" "RECT" '(' "integer" "integer" ')' '(' "integer" "integer" ')' ';'
#line 399 "DefParser.yy"
                                                                                                    {delete (yystack_[10].value.stringVal);}
#line 936 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 89: // specialnets_metal_layer: "string" "string" "integer"
#line 414 "DefParser.yy"
                                                {delete (yystack_[2].value.stringVal); delete (yystack_[1].value.stringVal);}
#line 942 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 90: // specialnets_metal_shape: '+' "SHAPE" "string" '(' "integer" "integer" ')' '(' "integer" '*' ')'
#line 416 "DefParser.yy"
                                                                                           {delete (yystack_[8].value.stringVal);}
#line 948 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 91: // specialnets_metal_shape: '+' "SHAPE" "string" '(' "integer" "integer" ')' '(' '*' "integer" ')'
#line 417 "DefParser.yy"
                                                                                                                   {delete (yystack_[8].value.stringVal);}
#line 954 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 92: // specialnets_metal_shape: '+' "SHAPE" "string" '(' "integer" "integer" ')' "string"
#line 418 "DefParser.yy"
                                                                                                      {delete (yystack_[5].value.stringVal); delete (yystack_[0].value.stringVal);}
#line 960 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 97: // specialnets_addon: specialnets_addon '+' "USE" "string"
#line 425 "DefParser.yy"
                                                                         {delete (yystack_[0].value.stringVal);}
#line 966 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 98: // specialnets_addon: specialnets_addon '+' "RECT" "string" '(' "integer" "integer" ')' '(' "integer" "integer" ')'
#line 426 "DefParser.yy"
                                                                                                          {delete (yystack_[8].value.stringVal);}
#line 972 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 99: // single_specialnet: '-' "binary numbers" specialnets_addon ';'
#line 428 "DefParser.yy"
                                                     {delete (yystack_[2].value.binaryVal);}
#line 978 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 100: // single_specialnet: '-' "string" specialnets_addon ';'
#line 429 "DefParser.yy"
                                                                     {delete (yystack_[2].value.stringVal);}
#line 984 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 105: // begin_nets: "NETS" "integer" ';'
#line 439 "DefParser.yy"
                                  {
					driver.net_cbk_size((yystack_[1].value.integerVal));
				}
#line 992 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 107: // node_pin_pair: '(' "string" "string" ')'
#line 446 "DefParser.yy"
                                      {
		   /** be careful, this callback will be invoked before net_cbk_name **/
				driver.net_cbk_pin(*(yystack_[2].value.stringVal), *(yystack_[1].value.stringVal)); 
                delete (yystack_[2].value.stringVal);
                delete (yystack_[1].value.stringVal);
			  }
#line 1003 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 111: // net_addon: net_addon '+' "USE" "string"
#line 458 "DefParser.yy"
                                         {delete (yystack_[0].value.stringVal);}
#line 1009 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 112: // net_addon: net_addon '+' "NONDEFAULTRULE" "string"
#line 459 "DefParser.yy"
                                                    {delete (yystack_[0].value.stringVal);}
#line 1015 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 113: // single_net: '-' "string" node_pin_pairs net_addon ';'
#line 462 "DefParser.yy"
                                                     {
				driver.net_cbk_name(*(yystack_[3].value.stringVal));
                delete (yystack_[3].value.stringVal);
		   }
#line 1024 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 114: // single_net: '-' "binary numbers" node_pin_pairs net_addon ';'
#line 466 "DefParser.yy"
                                                             {
				driver.net_cbk_name(*(yystack_[3].value.binaryVal));
                delete (yystack_[3].value.binaryVal);
		   }
#line 1033 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 119: // single_propterty: "COMPONENTPIN" "string" "string" ';'
#line 486 "DefParser.yy"
                                                      {delete (yystack_[2].value.stringVal); delete (yystack_[1].value.stringVal);}
#line 1039 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 120: // single_propterty: "DESIGN" "string" "string" "double" ';'
#line 487 "DefParser.yy"
                                                                       {delete (yystack_[3].value.stringVal); delete (yystack_[2].value.stringVal);}
#line 1045 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 121: // single_propterty: "NET" "string" "string" ';'
#line 488 "DefParser.yy"
                                             {delete (yystack_[2].value.stringVal); delete (yystack_[1].value.stringVal);}
#line 1051 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 129: // group_addon: group_addon '+' "REGION" "string"
#line 509 "DefParser.yy"
                                                {delete (yystack_[0].value.stringVal);}
#line 1057 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 130: // single_group: '-' "string" string_array group_addon ';'
#line 512 "DefParser.yy"
                                                       {delete (yystack_[3].value.stringVal); delete (yystack_[2].value.stringArrayVal);}
#line 1063 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 135: // begin_design: "DESIGN" "string" ';'
#line 525 "DefParser.yy"
                                     {
				driver.design_cbk(*(yystack_[1].value.stringVal));
                delete (yystack_[1].value.stringVal);
			 }
#line 1072 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 137: // block_unit: "UNITS" "DISTANCE" "MICRONS" "integer" ';'
#line 533 "DefParser.yy"
                                                            {
				driver.unit_cbk((yystack_[1].value.integerVal));
		   }
#line 1080 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;

  case 138: // block_diearea: "DIEAREA" '(' "integer" "integer" ')' '(' "integer" "integer" ')' ';'
#line 537 "DefParser.yy"
                                                                                {
				driver.diearea_cbk((yystack_[7].value.integerVal), (yystack_[6].value.integerVal), (yystack_[3].value.integerVal), (yystack_[2].value.integerVal));
			  }
#line 1088 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"
    break;


#line 1092 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"

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


  const short Parser::yypact_ninf_ = -257;

  const signed char Parser::yytable_ninf_ = -1;

  const short
  Parser::yypact_[] =
  {
    -257,    91,  -257,    10,    13,    37,    51,  -257,    33,     1,
      23,    27,    53,    56,    35,   104,   112,   120,   123,    65,
     124,   133,   136,   137,   138,   139,   140,   141,  -257,   129,
    -257,   117,  -257,   118,   -20,  -257,   -15,  -257,   -11,  -257,
      -7,  -257,     3,  -257,     4,  -257,     5,  -257,     6,  -257,
    -257,     8,  -257,  -257,  -257,  -257,    -9,  -257,  -257,  -257,
    -257,   135,   144,   145,    93,    94,    95,   149,   150,   130,
     151,  -257,    68,   155,   156,   102,   103,   105,   106,   107,
     108,  -257,  -257,  -257,   142,   157,  -257,  -257,   -20,   125,
     162,  -257,  -257,   -15,   119,   165,  -257,  -257,   -11,   152,
     167,  -257,  -257,    -7,   154,   170,  -257,  -257,     3,   121,
     -18,  -257,  -257,     4,   143,    18,  -257,  -257,     5,   158,
      19,  -257,  -257,     6,   127,   172,  -257,  -257,     8,   171,
    -257,  -257,   175,   178,   181,  -257,  -257,  -257,   180,   182,
    -257,   183,   159,  -257,   173,   174,  -257,  -257,  -257,  -257,
    -257,  -257,  -257,  -257,  -257,  -257,  -257,  -257,  -257,  -257,
    -257,   132,  -257,  -257,  -257,   184,  -257,  -257,  -257,  -257,
    -257,  -257,  -257,   188,   153,   160,  -257,  -257,  -257,  -257,
    -257,  -257,  -257,  -257,   134,   134,  -257,  -257,  -257,   191,
    -257,  -257,  -257,   146,   147,   194,   195,   148,   161,  -257,
     197,   198,    25,    39,   200,   163,  -257,    55,   168,   164,
     166,    59,    61,   193,  -257,   134,   134,  -257,   202,  -257,
     169,   204,   176,  -257,  -257,   189,   190,  -257,    70,  -257,
      48,   208,   211,    63,    66,  -257,    11,   177,   212,   213,
    -257,    36,  -257,   215,  -257,    69,    72,  -257,    74,   214,
     201,  -257,   218,   219,   220,   224,   191,   225,   227,   228,
     230,   231,   233,  -257,   179,   229,  -257,   192,  -257,    29,
     185,   234,   236,   237,   239,   243,   244,   245,   246,   247,
     248,   196,   249,   199,  -257,    71,  -257,  -257,   206,   254,
     255,   207,   203,  -257,   256,   202,   259,   260,   261,   205,
     221,  -257,  -257,   209,   217,   216,   262,   263,  -257,  -257,
     222,  -257,   268,   223,   226,   269,  -257,   232,   238,  -257,
     196,  -257,   270,   272,   273,   235,   187,   191,  -257,  -257,
    -257,   271,  -257,   277,   278,  -257,  -257,   131,  -257,   280,
     282,   240,   242,   250,  -257,   283,   284,  -257,  -257,  -257,
    -257,   241,   285,     0,   287,   288,   252,   290,   291,   251,
     294,   253,   295,   300,   301,   257,  -257,   286,  -257,  -257,
     258,   302,   264,   265,   289,   266,   303,   304,   305,   267,
     308,   311,   274,  -257,   310,   313,  -257,   275,   317,   276,
     279,   281,   318,   319,   320,  -257,  -257,   321,   292,   293,
     296,   326,   297,   298,   327,   328,   299,  -257,  -257,   329,
      -3,  -257,   306,   307,  -257,   309,  -257,    -2,  -257,  -257,
    -257,   186,   330,   312,   314,  -257,  -257
  };

  const unsigned char
  Parser::yydefact_[] =
  {
     156,     0,     1,     0,     0,     0,     0,   157,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     8,   142,
      11,   143,    14,   144,     0,   145,     0,   146,     0,   147,
       0,   148,     0,   149,     0,   150,     0,   151,     0,   152,
     140,     0,   153,   139,   141,   154,     0,     4,     5,     6,
     135,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   122,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     9,    12,    15,     0,     0,    30,    27,     0,     0,
       0,    41,    38,     0,     0,     0,    52,    49,     0,     0,
       0,    64,    61,     0,     0,     0,    77,    74,     0,     0,
       0,    86,    83,     0,     0,     0,   104,   101,     0,     0,
       0,   118,   115,     0,     0,     0,   134,   131,     0,     0,
     158,   155,     0,     0,     0,    53,    65,   105,     0,     0,
     125,     0,     0,   123,     0,     0,    16,    87,    31,    42,
     126,    78,    17,    18,    29,    28,    32,    33,    40,    39,
      43,     0,    51,    50,    54,     0,    63,    62,    66,    67,
      76,    75,    79,     0,     0,     0,    85,    84,    88,    95,
      95,   103,   102,   106,     0,     0,   117,   116,   127,     0,
     133,   132,   136,     0,     0,     0,     0,     0,     0,   124,
       0,     0,     0,     0,     0,    44,    55,     0,     0,     0,
       0,     0,     0,     0,   108,   110,   110,     2,   128,   137,
       0,     0,     0,   121,   119,     0,     0,    26,     0,    37,
       0,     0,     0,     0,     0,    73,     0,     0,     0,     0,
     100,     0,    99,     0,   109,     0,     0,     3,     0,     0,
       0,   120,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    34,     0,     0,    48,     0,    60,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    96,     0,   113,     0,   114,   130,     0,     0,
       0,     0,     0,    19,     0,    21,     0,     0,     0,     0,
       0,    36,    46,     0,     0,    59,     0,     0,    68,    69,
       0,    72,     0,     0,     0,     0,    97,     0,     0,    93,
       0,   107,     0,     0,     0,     0,     0,     0,    13,    20,
      25,     0,    22,     0,     0,    47,    45,     0,    58,     0,
       0,     0,     0,     0,    89,     0,     0,    94,   111,   112,
     129,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   138,     0,    10,    23,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    35,     0,     0,    70,     0,     0,     0,
       0,     0,     0,     0,     0,    56,    57,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    80,    81,     0,
       0,     7,     0,     0,    82,     0,    92,     0,    24,    71,
      98,     0,     0,     0,     0,    90,    91
  };

  const short
  Parser::yypgoto_[] =
  {
    -257,  -256,  -257,   315,  -257,   316,  -257,   322,  -257,  -257,
     323,  -257,   324,  -257,  -257,  -257,   325,  -257,   331,  -257,
    -257,  -257,   332,  -257,  -257,   333,  -257,  -257,  -257,   334,
    -257,   335,  -257,  -257,  -257,   336,  -257,   337,  -257,  -257,
    -257,   338,   339,  -257,  -257,  -257,   340,   -32,   -47,  -257,
     210,   341,  -257,  -257,  -257,   342,   -79,   343,   100,   344,
    -257,  -257,   345,  -257,  -257,  -257,   346,  -257,   347,  -257,
    -257,  -257,  -257,  -257,  -257,   348,  -257,  -257
  };

  const short
  Parser::yydefgoto_[] =
  {
       0,   218,     7,    28,    29,    30,    31,    32,    33,    34,
      86,   202,    87,    88,    35,    36,    91,   203,    92,    93,
      37,    38,    96,   233,   205,    97,    98,    39,    40,   101,
     234,   102,   103,    41,    42,   106,   207,   107,   108,    43,
      44,   111,   112,   113,    45,    46,   116,   281,   319,   282,
     211,   117,   118,    47,    48,   121,   214,   215,   245,   122,
     123,    49,    71,    72,    50,    51,   126,   248,   127,   128,
      52,     8,   130,    53,    54,    55,    56,     1
  };

  const short
  Parser::yytable_[] =
  {
     295,   421,   416,    13,    84,   247,    14,    15,   173,    89,
      16,    17,    18,    94,     9,   129,   270,    99,    19,    10,
      20,    21,    22,   179,   184,   180,   185,   104,   109,   114,
     119,    23,   124,   271,   305,    24,   272,   273,   174,   175,
      25,   278,    85,    11,    26,    13,    27,    90,    14,    15,
     274,    95,    16,    17,    18,   100,    12,   417,   368,    57,
      19,   422,    20,    21,    22,   105,   110,   115,   120,    61,
     125,   353,   306,    23,   261,   279,    67,    24,   280,    67,
     262,    58,    25,   227,   228,    59,    26,    68,    27,    69,
      68,     2,   142,    70,   263,    62,    70,   229,   230,     3,
       4,     5,     6,   254,   255,   256,   257,   258,   259,    63,
     322,    60,   260,   235,   236,    64,   323,   240,   241,   242,
     241,   266,   267,    65,   268,   269,    66,   284,   285,    73,
     286,   285,   287,   288,   357,   358,   244,   244,    74,    75,
      76,    77,    78,    79,    80,    15,    20,   133,    21,   132,
     134,   135,   136,   137,   138,   139,   141,   140,   144,   145,
     146,   147,   153,   148,   149,   150,   151,   157,   160,   156,
     161,   164,   165,   152,   168,   169,   172,   189,   193,   183,
     188,   194,   192,   178,   195,   196,   199,   197,   198,   206,
     200,   201,   204,   208,   213,   209,   217,   221,   243,   222,
     225,   226,   210,   231,   219,   352,   223,   247,   220,   250,
     237,   264,   252,   253,   265,   276,   277,   289,   290,   224,
     283,   291,   292,   232,   238,   293,   239,   294,   296,   249,
     297,   298,   303,   327,   251,   299,   300,   275,   301,   308,
     302,   309,   310,   304,   311,   307,   312,   313,   314,   423,
     320,   315,   316,   317,   278,   318,   324,   325,   326,   329,
     321,   328,   330,   331,   332,   333,   339,   338,   334,   336,
     335,   341,   344,   347,   354,   348,   337,   349,   350,   346,
     355,   356,   340,   359,   342,   360,   364,   343,   367,   365,
     369,   370,   345,   372,   386,   373,   351,   375,   377,   366,
     371,   361,   362,   378,   379,   383,   388,   389,   390,   381,
     363,   392,   374,   376,   393,   395,   246,   380,   396,   382,
     398,   402,   403,   404,   405,   384,   385,   387,   391,   409,
     412,   413,   415,   424,   394,   397,     0,   399,     0,     0,
     400,   401,     0,     0,    81,     0,     0,    82,     0,     0,
       0,   407,     0,   406,   408,    83,   411,   414,   410,     0,
       0,     0,     0,     0,     0,     0,     0,   418,   419,     0,
     420,     0,     0,   425,     0,   426,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     212,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   131,     0,     0,     0,     0,     0,
       0,   154,   155,     0,     0,     0,     0,   143,   158,     0,
       0,     0,     0,     0,   159,     0,     0,     0,     0,     0,
     162,   163,     0,     0,     0,     0,     0,   166,   167,     0,
       0,     0,     0,     0,   170,   171,     0,     0,     0,     0,
       0,   176,   177,     0,     0,     0,     0,     0,   181,   182,
       0,     0,     0,     0,     0,   186,     0,   187,     0,     0,
       0,     0,     0,     0,   190,   191,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   216
  };

  const short
  Parser::yycheck_[] =
  {
     256,     3,     5,    12,    24,     5,    15,    16,    26,    24,
      19,    20,    21,    24,     4,    24,     5,    24,    27,     6,
      29,    30,    31,     5,     5,     7,     7,    24,    24,    24,
      24,    40,    24,    22,     5,    44,    25,    26,    56,    57,
      49,     5,    62,     6,    53,    12,    55,    62,    15,    16,
      39,    62,    19,    20,    21,    62,     5,    60,    58,    58,
      27,    63,    29,    30,    31,    62,    62,    62,    62,    13,
      62,   327,    43,    40,    26,    39,    11,    44,    42,    11,
      32,    58,    49,    58,    59,    58,    53,    22,    55,    24,
      22,     0,    24,    28,    46,    60,    28,    58,    59,     8,
       9,    10,    11,    33,    34,    35,    36,    37,    38,     5,
      39,    58,    42,    58,    59,     3,    45,    58,    59,    58,
      59,    58,    59,     3,    58,    59,     3,    58,    59,     5,
      58,    59,    58,    59,     3,     4,   215,   216,     5,     3,
       3,     3,     3,     3,     3,    16,    29,     3,    30,    14,
       5,    58,    58,    58,     5,     5,     5,    27,     3,     3,
      58,    58,     5,    58,    58,    58,    58,     5,    49,    44,
       5,    19,     5,    31,    20,     5,    55,     5,     3,    21,
      53,     3,    11,    40,     3,     5,    27,     5,     5,     5,
      17,    17,    60,     5,    60,    42,     5,     3,     5,     4,
       3,     3,    42,     3,    58,    18,    58,     5,    61,     5,
      42,     3,    23,    23,     3,     3,     3,     3,    17,    58,
       5,     3,     3,    60,    60,     5,    60,     3,     3,    60,
       3,     3,     3,    26,    58,     5,     5,    60,     5,     5,
      61,     5,     5,    51,     5,    60,     3,     3,     3,    63,
     282,     5,     5,     5,     5,    59,    50,     3,     3,     3,
      61,    58,     3,     3,     3,    60,     3,     5,    47,    52,
      61,     3,     3,   320,     3,     5,    60,     5,     5,    41,
       3,     3,    60,     3,    61,     3,     3,    61,     3,     5,
       3,     3,    60,     3,     5,     4,    61,     3,     3,    58,
      48,    61,    60,     3,     3,     3,     3,     3,     3,    23,
      60,     3,    61,    60,     3,     5,   216,    60,     5,    61,
       3,     3,     3,     3,     3,    61,    61,    61,    61,     3,
       3,     3,     3,     3,    60,    60,    -1,    61,    -1,    -1,
      61,    60,    -1,    -1,    29,    -1,    -1,    31,    -1,    -1,
      -1,    58,    -1,    61,    58,    33,    58,    58,    61,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    61,    61,    -1,
      61,    -1,    -1,    61,    -1,    61,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     180,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    56,    -1,    -1,    -1,    -1,    -1,
      -1,    88,    88,    -1,    -1,    -1,    -1,    72,    93,    -1,
      -1,    -1,    -1,    -1,    93,    -1,    -1,    -1,    -1,    -1,
      98,    98,    -1,    -1,    -1,    -1,    -1,   103,   103,    -1,
      -1,    -1,    -1,    -1,   108,   108,    -1,    -1,    -1,    -1,
      -1,   113,   113,    -1,    -1,    -1,    -1,    -1,   118,   118,
      -1,    -1,    -1,    -1,    -1,   123,    -1,   123,    -1,    -1,
      -1,    -1,    -1,    -1,   128,   128,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   185
  };

  const unsigned char
  Parser::yystos_[] =
  {
       0,   141,     0,     8,     9,    10,    11,    66,   135,     4,
       6,     6,     5,    12,    15,    16,    19,    20,    21,    27,
      29,    30,    31,    40,    44,    49,    53,    55,    67,    68,
      69,    70,    71,    72,    73,    78,    79,    84,    85,    91,
      92,    97,    98,   103,   104,   108,   109,   117,   118,   125,
     128,   129,   134,   137,   138,   139,   140,    58,    58,    58,
      58,    13,    60,     5,     3,     3,     3,    11,    22,    24,
      28,   126,   127,     5,     5,     3,     3,     3,     3,     3,
       3,    67,    69,    71,    24,    62,    74,    76,    77,    24,
      62,    80,    82,    83,    24,    62,    86,    89,    90,    24,
      62,    93,    95,    96,    24,    62,    99,   101,   102,    24,
      62,   105,   106,   107,    24,    62,   110,   115,   116,    24,
      62,   119,   123,   124,    24,    62,   130,   132,   133,    24,
     136,   139,    14,     3,     5,    58,    58,    58,     5,     5,
      27,     5,    24,   126,     3,     3,    58,    58,    58,    58,
      58,    58,    31,     5,    74,    76,    44,     5,    80,    82,
      49,     5,    86,    89,    19,     5,    93,    95,    20,     5,
      99,   101,    55,    26,    56,    57,   105,   106,    40,     5,
       7,   110,   115,    21,     5,     7,   119,   123,    53,     5,
     130,   132,    11,     3,     3,     3,     5,     5,     5,    27,
      17,    17,    75,    81,    60,    88,     5,   100,     5,    42,
      42,   114,   114,    60,   120,   121,   121,     5,    65,    58,
      61,     3,     4,    58,    58,     3,     3,    58,    59,    58,
      59,     3,    60,    87,    94,    58,    59,    42,    60,    60,
      58,    59,    58,     5,   120,   122,   122,     5,   131,    60,
       5,    58,    23,    23,    33,    34,    35,    36,    37,    38,
      42,    26,    32,    46,     3,     3,    58,    59,    58,    59,
       5,    22,    25,    26,    39,    60,     3,     3,     5,    39,
      42,   111,   113,     5,    58,    59,    58,    58,    59,     3,
      17,     3,     3,     5,     3,    65,     3,     3,     3,     5,
       5,     5,    61,     3,    51,     5,    43,    60,     5,     5,
       5,     5,     3,     3,     3,     5,     5,     5,    59,   112,
     111,    61,    39,    45,    50,     3,     3,    26,    58,     3,
       3,     3,     3,    60,    47,    61,    52,    60,     5,     3,
      60,     3,    61,    61,     3,    60,    41,   112,     5,     5,
       5,    61,    18,    65,     3,     3,     3,     3,     4,     3,
       3,    61,    60,    60,     3,     5,    58,     3,    58,     3,
       3,    48,     3,     4,    61,     3,    60,     3,     3,     3,
      60,    23,    61,     3,    61,    61,     5,    61,     3,     3,
       3,    61,     3,     3,    60,     5,     5,    60,     3,    61,
      61,    60,     3,     3,     3,     3,    61,    58,    58,     3,
      61,    58,     3,     3,    58,     3,     5,    60,    61,    61,
      61,     3,    63,    63,     3,    61,    61
  };

  const unsigned char
  Parser::yyr1_[] =
  {
       0,    64,    65,    65,    66,    66,    66,    67,    68,    68,
      69,    70,    70,    71,    72,    72,    73,    74,    75,    75,
      75,    75,    75,    75,    75,    75,    76,    77,    77,    78,
      78,    79,    80,    81,    81,    81,    81,    82,    83,    83,
      84,    84,    85,    86,    87,    87,    88,    88,    89,    90,
      90,    91,    91,    92,    93,    94,    94,    94,    94,    94,
      95,    96,    96,    97,    97,    98,    99,   100,   100,   100,
     100,   100,   100,   101,   102,   102,   103,   103,   104,   105,
     106,   106,   106,   107,   107,   108,   108,   109,   110,   111,
     112,   112,   112,   113,   113,   114,   114,   114,   114,   115,
     115,   116,   116,   117,   117,   118,   119,   120,   121,   121,
     122,   122,   122,   123,   123,   124,   124,   125,   125,   126,
     126,   126,   127,   127,   128,   128,   129,   130,   131,   131,
     132,   133,   133,   134,   134,   135,   136,   137,   138,   139,
     139,   139,   139,   139,   139,   139,   139,   139,   139,   139,
     139,   139,   139,   139,   140,   140,   141,   141,   141
  };

  const signed char
  Parser::yyr2_[] =
  {
       0,     2,     1,     2,     3,     3,     3,    14,     1,     2,
      10,     1,     2,     8,     1,     2,     3,     2,     0,     4,
       5,     4,     5,     7,    12,     5,     4,     1,     2,     3,
       2,     3,     2,     0,     3,     8,     4,     4,     1,     2,
       3,     2,     3,     2,     0,     4,     4,     5,     5,     1,
       2,     3,     2,     3,     2,     0,     8,     8,     4,     3,
       5,     1,     2,     3,     2,     3,     2,     0,     4,     4,
       8,    12,     4,     4,     1,     2,     3,     2,     3,     2,
      12,    12,    13,     1,     2,     3,     2,     3,     2,     3,
      11,    11,     8,     2,     3,     0,     3,     4,    12,     4,
       4,     1,     2,     3,     2,     3,     2,     4,     1,     2,
       0,     4,     4,     5,     5,     1,     2,     3,     2,     4,
       5,     4,     1,     2,     4,     3,     3,     2,     0,     4,
       5,     1,     2,     3,     2,     3,     2,     5,    10,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     2,     0,     2,     4
  };


#if DEFPARSERDEBUG || 1
  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a YYNTOKENS, nonterminals.
  const char*
  const Parser::yytname_[] =
  {
  "\"end of file\"", "error", "\"invalid token\"", "\"integer\"",
  "\"double\"", "\"string\"", "\"quoted chars\"", "\"binary numbers\"",
  "\"VERSION\"", "\"DIVIDERCHAR\"", "\"BUSBITCHARS\"", "\"DESIGN\"",
  "\"UNITS\"", "\"DISTANCE\"", "\"MICRONS\"", "\"DIEAREA\"", "\"ROW\"",
  "\"DO\"", "\"BY\"", "\"COMPONENTS\"", "\"PINS\"", "\"NETS\"", "\"NET\"",
  "\"STEP\"", "\"END\"", "\"DIRECTION\"", "\"LAYER\"",
  "\"PROPERTYDEFINITIONS\"", "\"COMPONENTPIN\"", "\"TRACKS\"",
  "\"GCELLGRID\"", "\"VIAS\"", "\"VIA\"", "\"VIARULE\"", "\"CUTSIZE\"",
  "\"LAYERS\"", "\"ROWCOL\"", "\"ENCLOSURE\"", "\"CUTSPACING\"", "\"USE\"",
  "\"SPECIALNETS\"", "\"SHAPE\"", "\"RECT\"", "\"SOURCE\"",
  "\"NONDEFAULTRULES\"", "\"NONDEFAULTRULE\"", "\"HARDSPACING\"",
  "\"WIDTH\"", "\"SPACING\"", "\"REGIONS\"", "\"REGION\"", "\"TYPE\"",
  "\"FENCE\"", "\"GROUPS\"", "\"GROUP\"", "\"BLOCKAGES\"", "\"PLACEMENT\"",
  "\"ROUTING\"", "';'", "'+'", "'('", "')'", "'-'", "'*'", "$accept",
  "string_array", "block_other", "single_row", "block_rows",
  "single_tracks", "block_tracks", "single_gcellgrid", "block_gcellgrid",
  "begin_vias", "end_vias", "via_addon", "single_via", "multiple_vias",
  "block_vias", "begin_nondefaultrules", "end_nondefaultrules",
  "nondefaultrule_addon", "single_nondefaultrule",
  "multiple_nondefaultrules", "block_nondefaultrules", "begin_regions",
  "end_regions", "region_addon", "region_points", "single_region",
  "multiple_regions", "block_regions", "begin_components",
  "end_components", "component_addon", "single_component",
  "multiple_components", "block_components", "begin_pins", "end_pins",
  "pin_addon", "single_pin", "multiple_pins", "block_pins",
  "begin_blockages", "end_blockages", "single_blockage",
  "multiple_blockages", "block_blockages", "begin_specialnets",
  "end_specialnets", "specialnets_metal_layer", "specialnets_metal_shape",
  "specialnets_metal_array", "specialnets_addon", "single_specialnet",
  "multiple_specialnets", "block_specialnets", "begin_nets", "end_nets",
  "node_pin_pair", "node_pin_pairs", "net_addon", "single_net",
  "multiple_nets", "block_nets", "single_propterty", "multiple_property",
  "block_propertydefinitions", "begin_groups", "end_groups", "group_addon",
  "single_group", "multiple_groups", "block_groups", "begin_design",
  "end_design", "block_unit", "block_diearea", "block_option",
  "block_design", "start", YY_NULLPTR
  };
#endif


#if DEFPARSERDEBUG
  const short
  Parser::yyrline_[] =
  {
       0,   167,   167,   171,   177,   180,   184,   190,   197,   198,
     202,   208,   209,   213,   218,   219,   223,   226,   229,   230,
     231,   232,   233,   234,   235,   236,   239,   242,   243,   246,
     247,   251,   253,   256,   257,   258,   259,   262,   265,   266,
     269,   270,   274,   276,   279,   280,   283,   284,   287,   290,
     291,   294,   295,   299,   303,   306,   307,   312,   317,   321,
     327,   333,   334,   337,   340,   345,   349,   352,   353,   357,
     361,   366,   370,   376,   381,   382,   385,   388,   393,   395,
     397,   398,   399,   401,   402,   404,   405,   410,   412,   414,
     416,   417,   418,   420,   421,   423,   424,   425,   426,   428,
     429,   431,   432,   434,   435,   439,   443,   446,   453,   454,
     457,   458,   459,   462,   466,   472,   473,   476,   479,   486,
     487,   488,   491,   492,   495,   498,   503,   505,   508,   509,
     512,   515,   516,   519,   520,   525,   530,   533,   537,   541,
     542,   543,   544,   545,   546,   547,   548,   549,   550,   551,
     552,   553,   554,   555,   558,   559,   563,   564,   565
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
#endif // DEFPARSERDEBUG

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
      60,    61,    63,    59,     2,    62,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    58,
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
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57
    };
    // Last valid token kind.
    const int code_max = 312;

    if (t <= 0)
      return symbol_kind::S_YYEOF;
    else if (t <= code_max)
      return static_cast <symbol_kind_type> (translate_table[t]);
    else
      return symbol_kind::S_YYUNDEF;
  }

} // DefParser
#line 1926 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/def/bison/DefParser.cc"

#line 570 "DefParser.yy"
 /*** Additional Code ***/

void DefParser::Parser::error(const Parser::location_type& l,
			    const std::string& m)
{
    driver.error(l, m);
}
