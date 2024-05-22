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
#define yylex   LefParserlex

// First part of user prologue.
#line 4 "LefParser.yy"
 /*** C/C++ Declarations ***/

#include <stdio.h>
#include <string>
#include <vector>
#include <string.h>
#include "LefDataBase.h"


#line 53 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"


#include "LefParser.h"

// Second part of user prologue.
#line 526 "LefParser.yy"


#include "LefDriver.h"
#include "LefScanner.h"

/* this "connects" the bison parser in the driver to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the driver context. */
#undef yylex
#define yylex driver.lexer->lex


#line 72 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"



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
#if LEFPARSERDEBUG

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

#else // !LEFPARSERDEBUG

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Symbol)  YY_USE (Symbol)
# define YY_REDUCE_PRINT(Rule)           static_cast<void> (0)
# define YY_STACK_PRINT()                static_cast<void> (0)

#endif // !LEFPARSERDEBUG

#define yyerrok         (yyerrstatus_ = 0)
#define yyclearin       (yyla.clear ())

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)

namespace LefParser {
#line 165 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"

  /// Build a parser object.
  Parser::Parser (class Driver& driver_yyarg)
#if LEFPARSERDEBUG
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
#line 511 "LefParser.yy"
                    { delete (yysym.value.stringVal); }
#line 377 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
        break;

      case symbol_kind::S_QSTRING: // "qstring"
#line 511 "LefParser.yy"
                    { delete (yysym.value.qstringVal); }
#line 383 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
        break;

      case symbol_kind::S_BINARY: // "binary numbers"
#line 511 "LefParser.yy"
                    { delete (yysym.value.binaryVal); }
#line 389 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
        break;

      case symbol_kind::S_GSTRING: // GSTRING
#line 513 "LefParser.yy"
                    { delete (yysym.value.stringVal); }
#line 395 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
        break;

      case symbol_kind::S_pt: // pt
#line 512 "LefParser.yy"
                    { delete (yysym.value.pt); }
#line 401 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
        break;

      default:
        break;
    }
  }

#if LEFPARSERDEBUG
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

#if LEFPARSERDEBUG
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
#endif // LEFPARSERDEBUG

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
#line 47 "LefParser.yy"
{
    // initialize the initial location object
    yyla.location.begin.filename = yyla.location.end.filename = &driver.streamname;
}

#line 544 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"


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
  case 2: // NUMBER: "integer"
#line 543 "LefParser.yy"
                 {(yylhs.value.doubleVal)=(yystack_[0].value.integerVal);}
#line 682 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 3: // NUMBER: "double"
#line 544 "LefParser.yy"
                    {(yylhs.value.doubleVal)=(yystack_[0].value.doubleVal);}
#line 688 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 4: // GSTRING: "string"
#line 546 "LefParser.yy"
                 {(yylhs.value.stringVal) = (yystack_[0].value.stringVal);}
#line 694 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 5: // GSTRING: "DEFINE"
#line 547 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("DEFINE");}
#line 700 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 6: // GSTRING: "DEFINEB"
#line 548 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("DEFINEB");}
#line 706 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 7: // GSTRING: "DEFINES"
#line 549 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("DEFINES");}
#line 712 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 8: // GSTRING: "MESSAGE"
#line 550 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("MESSAGE");}
#line 718 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 9: // GSTRING: "CREATEFILE"
#line 551 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("CREATEFILE");}
#line 724 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 10: // GSTRING: "OPENFILE"
#line 552 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("OPENFILE");}
#line 730 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 11: // GSTRING: "CLOSEFILE"
#line 553 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("CLOSEFILE");}
#line 736 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 12: // GSTRING: "WARNING"
#line 554 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("WARNING");}
#line 742 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 13: // GSTRING: "ERROR"
#line 555 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("ERROR");}
#line 748 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 14: // GSTRING: "FATALERROR"
#line 556 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("FATALERROR");}
#line 754 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 15: // GSTRING: "ABOVE"
#line 557 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("ABOVE");}
#line 760 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 16: // GSTRING: "ABUT"
#line 558 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("ABUT");}
#line 766 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 17: // GSTRING: "ABUTMENT"
#line 559 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("ABUTMENT");}
#line 772 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 18: // GSTRING: "ACCURRENTDENSITY"
#line 560 "LefParser.yy"
                     {(yylhs.value.stringVal) = new std::string ("ACCURRENTDENSITY");}
#line 778 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 19: // GSTRING: "ACTIVE"
#line 561 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("ACTIVE");}
#line 784 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 20: // GSTRING: "ADJACENTCUTS"
#line 562 "LefParser.yy"
                 {(yylhs.value.stringVal) = new std::string ("ADJACENTCUTS");}
#line 790 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 21: // GSTRING: "ANALOG"
#line 563 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("ANALOG");}
#line 796 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 22: // GSTRING: "AND"
#line 564 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("AND");}
#line 802 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 23: // GSTRING: "ANTENNAAREAFACTOR"
#line 565 "LefParser.yy"
                      {(yylhs.value.stringVal) = new std::string ("ANTENNAAREAFACTOR");}
#line 808 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 24: // GSTRING: "ANTENNAAREADIFFREDUCEPWL"
#line 566 "LefParser.yy"
                             {(yylhs.value.stringVal) = new std::string ("ANTENNAAREADIFFREDUCEPWL");}
#line 814 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 25: // GSTRING: "ANTENNAAREAMINUSDIFF"
#line 567 "LefParser.yy"
                         {(yylhs.value.stringVal) = new std::string ("ANTENNAAREAMINUSDIFF");}
#line 820 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 26: // GSTRING: "ANTENNAAREARATIO"
#line 568 "LefParser.yy"
                     {(yylhs.value.stringVal) = new std::string ("ANTENNAAREARATIO");}
#line 826 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 27: // GSTRING: "ANTENNACELL"
#line 569 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("ANTENNACELL");}
#line 832 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 28: // GSTRING: "ANTENNACUMAREARATIO"
#line 570 "LefParser.yy"
                        {(yylhs.value.stringVal) = new std::string ("ANTENNACUMAREARATIO");}
#line 838 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 29: // GSTRING: "ANTENNACUMDIFFAREARATIO"
#line 571 "LefParser.yy"
                            {(yylhs.value.stringVal) = new std::string ("ANTENNACUMDIFFAREARATIO");}
#line 844 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 30: // GSTRING: "ANTENNACUMDIFFSIDEAREARATIO"
#line 572 "LefParser.yy"
                                {(yylhs.value.stringVal) = new std::string ("ANTENNACUMDIFFSIDEAREARATIO");}
#line 850 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 31: // GSTRING: "ANTENNACUMROUTINGPLUSCUT"
#line 573 "LefParser.yy"
                             {(yylhs.value.stringVal) = new std::string ("ANTENNACUMROUTINGPLUSCUT");}
#line 856 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 32: // GSTRING: "ANTENNACUMSIDEAREARATIO"
#line 574 "LefParser.yy"
                            {(yylhs.value.stringVal) = new std::string ("ANTENNACUMSIDEAREARATIO");}
#line 862 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 33: // GSTRING: "ANTENNADIFFAREA"
#line 575 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("ANTENNADIFFAREA");}
#line 868 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 34: // GSTRING: "ANTENNADIFFAREARATIO"
#line 576 "LefParser.yy"
                         {(yylhs.value.stringVal) = new std::string ("ANTENNADIFFAREARATIO");}
#line 874 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 35: // GSTRING: "ANTENNADIFFSIDEAREARATIO"
#line 577 "LefParser.yy"
                             {(yylhs.value.stringVal) = new std::string ("ANTENNADIFFSIDEAREARATIO");}
#line 880 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 36: // GSTRING: "ANTENNAGATEAREA"
#line 578 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("ANTENNAGATEAREA");}
#line 886 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 37: // GSTRING: "ANTENNAGATEPLUSDIFF"
#line 579 "LefParser.yy"
                        {(yylhs.value.stringVal) = new std::string ("ANTENNAGATEPLUSDIFF");}
#line 892 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 38: // GSTRING: "ANTENNAINOUTDIFFAREA"
#line 580 "LefParser.yy"
                         {(yylhs.value.stringVal) = new std::string ("ANTENNAINOUTDIFFAREA");}
#line 898 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 39: // GSTRING: "ANTENNAINPUTGATEAREA"
#line 581 "LefParser.yy"
                         {(yylhs.value.stringVal) = new std::string ("ANTENNAINPUTGATEAREA");}
#line 904 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 40: // GSTRING: "ANTENNALENGTHFACTOR"
#line 582 "LefParser.yy"
                        {(yylhs.value.stringVal) = new std::string ("ANTENNALENGTHFACTOR");}
#line 910 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 41: // GSTRING: "ANTENNAMAXAREACAR"
#line 583 "LefParser.yy"
                      {(yylhs.value.stringVal) = new std::string ("ANTENNAMAXAREACAR");}
#line 916 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 42: // GSTRING: "ANTENNAMAXCUTCAR"
#line 584 "LefParser.yy"
                     {(yylhs.value.stringVal) = new std::string ("ANTENNAMAXCUTCAR");}
#line 922 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 43: // GSTRING: "ANTENNAMAXSIDEAREACAR"
#line 585 "LefParser.yy"
                          {(yylhs.value.stringVal) = new std::string ("ANTENNAMAXSIDEAREACAR");}
#line 928 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 44: // GSTRING: "ANTENNAMETALAREA"
#line 586 "LefParser.yy"
                     {(yylhs.value.stringVal) = new std::string ("ANTENNAMETALAREA");}
#line 934 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 45: // GSTRING: "ANTENNAMETALLENGTH"
#line 587 "LefParser.yy"
                       {(yylhs.value.stringVal) = new std::string ("ANTENNAMETALLENGTH");}
#line 940 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 46: // GSTRING: "ANTENNAMODEL"
#line 588 "LefParser.yy"
                 {(yylhs.value.stringVal) = new std::string ("ANTENNAMODEL");}
#line 946 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 47: // GSTRING: "ANTENNAOUTPUTDIFFAREA"
#line 589 "LefParser.yy"
                          {(yylhs.value.stringVal) = new std::string ("ANTENNAOUTPUTDIFFAREA");}
#line 952 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 48: // GSTRING: "ANTENNAPARTIALCUTAREA"
#line 590 "LefParser.yy"
                          {(yylhs.value.stringVal) = new std::string ("ANTENNAPARTIALCUTAREA");}
#line 958 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 49: // GSTRING: "ANTENNAPARTIALMETALAREA"
#line 591 "LefParser.yy"
                            {(yylhs.value.stringVal) = new std::string ("ANTENNAPARTIALMETALAREA");}
#line 964 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 50: // GSTRING: "ANTENNAPARTIALMETALSIDEAREA"
#line 592 "LefParser.yy"
                                {(yylhs.value.stringVal) = new std::string ("ANTENNAPARTIALMETALSIDEAREA");}
#line 970 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 51: // GSTRING: "ANTENNASIDEAREARATIO"
#line 593 "LefParser.yy"
                         {(yylhs.value.stringVal) = new std::string ("ANTENNASIDEAREARATIO");}
#line 976 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 52: // GSTRING: "ANTENNASIZE"
#line 594 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("ANTENNASIZE");}
#line 982 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 53: // GSTRING: "ANTENNASIDEAREAFACTOR"
#line 595 "LefParser.yy"
                          {(yylhs.value.stringVal) = new std::string ("ANTENNASIDEAREAFACTOR");}
#line 988 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 54: // GSTRING: "ANYEDGE"
#line 596 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("ANYEDGE");}
#line 994 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 55: // GSTRING: "AREA"
#line 597 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("AREA");}
#line 1000 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 56: // GSTRING: "AREAIO"
#line 598 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("AREAIO");}
#line 1006 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 57: // GSTRING: "ARRAY"
#line 599 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("ARRAY");}
#line 1012 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 58: // GSTRING: "ARRAYCUTS"
#line 600 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("ARRAYCUTS");}
#line 1018 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 59: // GSTRING: "ARRAYSPACING"
#line 601 "LefParser.yy"
                 {(yylhs.value.stringVal) = new std::string ("ARRAYSPACING");}
#line 1024 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 60: // GSTRING: "AVERAGE"
#line 602 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("AVERAGE");}
#line 1030 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 61: // GSTRING: "BELOW"
#line 603 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("BELOW");}
#line 1036 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 62: // GSTRING: "BEGINEXT"
#line 604 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("BEGINEXT");}
#line 1042 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 63: // GSTRING: "BLACKBOX"
#line 605 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("BLACKBOX");}
#line 1048 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 64: // GSTRING: "BLOCK"
#line 606 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("BLOCK");}
#line 1054 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 65: // GSTRING: "BOTTOMLEFT"
#line 607 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("BOTTOMLEFT");}
#line 1060 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 66: // GSTRING: "BOTTOMRIGHT"
#line 608 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("BOTTOMRIGHT");}
#line 1066 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 67: // GSTRING: "BUMP"
#line 609 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("BUMP");}
#line 1072 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 68: // GSTRING: "BUSBITCHARS"
#line 610 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("BUSBITCHARS");}
#line 1078 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 69: // GSTRING: "BUFFER"
#line 611 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("BUFFER");}
#line 1084 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 70: // GSTRING: "BY"
#line 612 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("BY");}
#line 1090 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 71: // GSTRING: "CANNOTOCCUPY"
#line 613 "LefParser.yy"
                 {(yylhs.value.stringVal) = new std::string ("CANNOTOCCUPY");}
#line 1096 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 72: // GSTRING: "CANPLACE"
#line 614 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("CANPLACE");}
#line 1102 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 73: // GSTRING: "CAPACITANCE"
#line 615 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("CAPACITANCE");}
#line 1108 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 74: // GSTRING: "CAPMULTIPLIER"
#line 616 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("CAPMULTIPLIER");}
#line 1114 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 75: // GSTRING: "CENTERTOCENTER"
#line 617 "LefParser.yy"
                   {(yylhs.value.stringVal) = new std::string ("CENTERTOCENTER");}
#line 1120 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 76: // GSTRING: "CLASS"
#line 618 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("CLASS");}
#line 1126 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 77: // GSTRING: "CLEARANCEMEASURE"
#line 619 "LefParser.yy"
                     {(yylhs.value.stringVal) = new std::string ("CLEARANCEMEASURE");}
#line 1132 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 78: // GSTRING: "CLOCK"
#line 620 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("CLOCK");}
#line 1138 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 79: // GSTRING: "CLOCKTYPE"
#line 621 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("CLOCKTYPE");}
#line 1144 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 80: // GSTRING: "COLUMNMAJOR"
#line 622 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("COLUMNMAJOR");}
#line 1150 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 81: // GSTRING: "CURRENTDEN"
#line 623 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("CURRENTDEN");}
#line 1156 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 82: // GSTRING: "COMPONENTPIN"
#line 624 "LefParser.yy"
                 {(yylhs.value.stringVal) = new std::string ("COMPONENTPIN");}
#line 1162 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 83: // GSTRING: "CORE"
#line 625 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("CORE");}
#line 1168 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 84: // GSTRING: "CORNER"
#line 626 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("CORNER");}
#line 1174 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 85: // GSTRING: "CORRECTIONFACTOR"
#line 627 "LefParser.yy"
                     {(yylhs.value.stringVal) = new std::string ("CORRECTIONFACTOR");}
#line 1180 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 86: // GSTRING: "CORRECTIONTABLE"
#line 628 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("CORRECTIONTABLE");}
#line 1186 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 87: // GSTRING: "COVER"
#line 629 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("COVER");}
#line 1192 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 88: // GSTRING: "CPERSQDIST"
#line 630 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("CPERSQDIST");}
#line 1198 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 89: // GSTRING: "CURRENT"
#line 631 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("CURRENT");}
#line 1204 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 90: // GSTRING: "CURRENTSOURCE"
#line 632 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("CURRENTSOURCE");}
#line 1210 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 91: // GSTRING: "CUT"
#line 633 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("CUT");}
#line 1216 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 92: // GSTRING: "CUTAREA"
#line 634 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("CUTAREA");}
#line 1222 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 93: // GSTRING: "CUTSIZE"
#line 635 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("CUTSIZE");}
#line 1228 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 94: // GSTRING: "CUTSPACING"
#line 636 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("CUTSPACING");}
#line 1234 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 95: // GSTRING: "DATA"
#line 637 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("DATA");}
#line 1240 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 96: // GSTRING: "DATABASE"
#line 638 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("DATABASE");}
#line 1246 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 97: // GSTRING: "DCCURRENTDENSITY"
#line 639 "LefParser.yy"
                     {(yylhs.value.stringVal) = new std::string ("DCCURRENTDENSITY");}
#line 1252 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 98: // GSTRING: "DEFAULT"
#line 640 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("DEFAULT");}
#line 1258 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 99: // GSTRING: "DEFAULTCAP"
#line 641 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("DEFAULTCAP");}
#line 1264 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 100: // GSTRING: "DELAY"
#line 642 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("DELAY");}
#line 1270 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 101: // GSTRING: "DENSITY"
#line 643 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("DENSITY");}
#line 1276 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 102: // GSTRING: "DENSITYCHECKSTEP"
#line 644 "LefParser.yy"
                     {(yylhs.value.stringVal) = new std::string ("DENSITYCHECKSTEP");}
#line 1282 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 103: // GSTRING: "DENSITYCHECKWINDOW"
#line 645 "LefParser.yy"
                       {(yylhs.value.stringVal) = new std::string ("DENSITYCHECKWINDOW");}
#line 1288 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 104: // GSTRING: "DESIGNRULEWIDTH"
#line 646 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("DESIGNRULEWIDTH");}
#line 1294 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 105: // GSTRING: "DIAG45"
#line 647 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("DIAG45");}
#line 1300 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 106: // GSTRING: "DIAG135"
#line 648 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("DIAG135");}
#line 1306 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 107: // GSTRING: "DIAGMINEDGELENGTH"
#line 649 "LefParser.yy"
                      {(yylhs.value.stringVal) = new std::string ("DIAGMINEDGELENGTH");}
#line 1312 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 108: // GSTRING: "DIAGSPACING"
#line 650 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("DIAGSPACING");}
#line 1318 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 109: // GSTRING: "DIAGPITCH"
#line 651 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("DIAGPITCH");}
#line 1324 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 110: // GSTRING: "DIAGWIDTH"
#line 652 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("DIAGWIDTH");}
#line 1330 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 111: // GSTRING: "DIELECTRIC"
#line 653 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("DIELECTRIC");}
#line 1336 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 112: // GSTRING: "DIFFUSEONLY"
#line 654 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("DIFFUSEONLY");}
#line 1342 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 113: // GSTRING: "DIRECTION"
#line 655 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("DIRECTION");}
#line 1348 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 114: // GSTRING: "DIVIDERCHAR"
#line 656 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("DIVIDERCHAR");}
#line 1354 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 115: // GSTRING: "DO"
#line 657 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("DO");}
#line 1360 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 116: // GSTRING: "E"
#line 658 "LefParser.yy"
      {(yylhs.value.stringVal) = new std::string ("E");}
#line 1366 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 117: // GSTRING: "EDGECAPACITANCE"
#line 659 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("EDGECAPACITANCE");}
#line 1372 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 118: // GSTRING: "EDGERATE"
#line 660 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("EDGERATE");}
#line 1378 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 119: // GSTRING: "EDGERATESCALEFACTOR"
#line 661 "LefParser.yy"
                        {(yylhs.value.stringVal) = new std::string ("EDGERATESCALEFACTOR");}
#line 1384 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 120: // GSTRING: "EDGERATETHRESHOLD1"
#line 662 "LefParser.yy"
                       {(yylhs.value.stringVal) = new std::string ("EDGERATETHRESHOLD1");}
#line 1390 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 121: // GSTRING: "EDGERATETHRESHOLD2"
#line 663 "LefParser.yy"
                       {(yylhs.value.stringVal) = new std::string ("EDGERATETHRESHOLD2");}
#line 1396 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 122: // GSTRING: "EEQ"
#line 664 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("EEQ");}
#line 1402 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 123: // GSTRING: "ELSE"
#line 665 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("ELSE");}
#line 1408 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 124: // GSTRING: "ENCLOSURE"
#line 666 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("ENCLOSURE");}
#line 1414 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 125: // GSTRING: "END"
#line 667 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("END");}
#line 1420 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 126: // GSTRING: "ENDEXT"
#line 668 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("ENDEXT");}
#line 1426 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 127: // GSTRING: "ENDCAP"
#line 669 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("ENDCAP");}
#line 1432 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 128: // GSTRING: "ENDOFLINE"
#line 670 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("ENDOFLINE");}
#line 1438 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 129: // GSTRING: "ENDOFNOTCHWIDTH"
#line 671 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("ENDOFNOTCHWIDTH");}
#line 1444 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 130: // GSTRING: "EUCLIDEAN"
#line 672 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("EUCLIDEAN");}
#line 1450 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 131: // GSTRING: "EXCEPTEXTRACUT"
#line 673 "LefParser.yy"
                   {(yylhs.value.stringVal) = new std::string ("EXCEPTEXTRACUT");}
#line 1456 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 132: // GSTRING: "EXCEPTSAMEPGNET"
#line 674 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("EXCEPTSAMEPGNET");}
#line 1462 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 133: // GSTRING: "EXCEPTPGNET"
#line 675 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("EXCEPTPGNET");}
#line 1468 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 134: // GSTRING: "EXTENSION"
#line 676 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("EXTENSION");}
#line 1474 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 135: // GSTRING: "FALL"
#line 677 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("FALL");}
#line 1480 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 136: // GSTRING: "FALLCS"
#line 678 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("FALLCS");}
#line 1486 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 137: // GSTRING: "FALLRS"
#line 679 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("FALLRS");}
#line 1492 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 138: // GSTRING: "FALLSATCUR"
#line 680 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("FALLSATCUR");}
#line 1498 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 139: // GSTRING: "FALLSATT1"
#line 681 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("FALLSATT1");}
#line 1504 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 140: // GSTRING: "FALLSLEWLIMIT"
#line 682 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("FALLSLEWLIMIT");}
#line 1510 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 141: // GSTRING: "FALLT0"
#line 683 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("FALLT0");}
#line 1516 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 142: // GSTRING: "FALLTHRESH"
#line 684 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("FALLTHRESH");}
#line 1522 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 143: // GSTRING: "FALLVOLTAGETHRESHOLD"
#line 685 "LefParser.yy"
                         {(yylhs.value.stringVal) = new std::string ("FALLVOLTAGETHRESHOLD");}
#line 1528 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 144: // GSTRING: "FALSE"
#line 686 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("FALSE");}
#line 1534 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 145: // GSTRING: "FE"
#line 687 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("FE");}
#line 1540 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 146: // GSTRING: "FEEDTHRU"
#line 688 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("FEEDTHRU");}
#line 1546 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 147: // GSTRING: "FILLACTIVESPACING"
#line 689 "LefParser.yy"
                      {(yylhs.value.stringVal) = new std::string ("FILLACTIVESPACING");}
#line 1552 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 148: // GSTRING: "FIXED"
#line 690 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("FIXED");}
#line 1558 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 149: // GSTRING: "FLIP"
#line 691 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("FLIP");}
#line 1564 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 150: // GSTRING: "FLOORPLAN"
#line 692 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("FLOORPLAN");}
#line 1570 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 151: // GSTRING: "FN"
#line 693 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("FN");}
#line 1576 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 152: // GSTRING: "FOREIGN"
#line 694 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("FOREIGN");}
#line 1582 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 153: // GSTRING: "FREQUENCY"
#line 695 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("FREQUENCY");}
#line 1588 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 154: // GSTRING: "FROMABOVE"
#line 696 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("FROMABOVE");}
#line 1594 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 155: // GSTRING: "FROMBELOW"
#line 697 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("FROMBELOW");}
#line 1600 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 156: // GSTRING: "FROMPIN"
#line 698 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("FROMPIN");}
#line 1606 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 157: // GSTRING: "FUNCTION"
#line 699 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("FUNCTION");}
#line 1612 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 158: // GSTRING: "FS"
#line 700 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("FS");}
#line 1618 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 159: // GSTRING: "FW"
#line 701 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("FW");}
#line 1624 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 160: // GSTRING: "GCELLGRID"
#line 702 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("GCELLGRID");}
#line 1630 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 161: // GSTRING: "GENERATE"
#line 703 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("GENERATE");}
#line 1636 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 162: // GSTRING: "GENERATED"
#line 704 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("GENERATED");}
#line 1642 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 163: // GSTRING: "GENERATOR"
#line 705 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("GENERATOR");}
#line 1648 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 164: // GSTRING: "GROUND"
#line 706 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("GROUND");}
#line 1654 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 165: // GSTRING: "GROUNDSENSITIVITY"
#line 707 "LefParser.yy"
                      {(yylhs.value.stringVal) = new std::string ("GROUNDSENSITIVITY");}
#line 1660 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 166: // GSTRING: "HARDSPACING"
#line 708 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("HARDSPACING");}
#line 1666 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 167: // GSTRING: "HEIGHT"
#line 709 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("HEIGHT");}
#line 1672 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 168: // GSTRING: "HISTORY"
#line 710 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("HISTORY");}
#line 1678 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 169: // GSTRING: "HOLD"
#line 711 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("HOLD");}
#line 1684 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 170: // GSTRING: "HORIZONTAL"
#line 712 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("HORIZONTAL");}
#line 1690 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 171: // GSTRING: "IF"
#line 713 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("IF");}
#line 1696 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 172: // GSTRING: "IMPLANT"
#line 714 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("IMPLANT");}
#line 1702 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 173: // GSTRING: "INFLUENCE"
#line 715 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("INFLUENCE");}
#line 1708 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 174: // GSTRING: "INOUT"
#line 716 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("INOUT");}
#line 1714 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 175: // GSTRING: "INOUTPINANTENNASIZE"
#line 717 "LefParser.yy"
                        {(yylhs.value.stringVal) = new std::string ("INOUTPINANTENNASIZE");}
#line 1720 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 176: // GSTRING: "INPUT"
#line 718 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("INPUT");}
#line 1726 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 177: // GSTRING: "INPUTPINANTENNASIZE"
#line 719 "LefParser.yy"
                        {(yylhs.value.stringVal) = new std::string ("INPUTPINANTENNASIZE");}
#line 1732 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 178: // GSTRING: "INPUTNOISEMARGIN"
#line 720 "LefParser.yy"
                     {(yylhs.value.stringVal) = new std::string ("INPUTNOISEMARGIN");}
#line 1738 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 179: // GSTRING: "INSIDECORNER"
#line 721 "LefParser.yy"
                 {(yylhs.value.stringVal) = new std::string ("INSIDECORNER");}
#line 1744 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 180: // GSTRING: "INTEGER"
#line 722 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("INTEGER");}
#line 1750 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 181: // GSTRING: "INTRINSIC"
#line 723 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("INTRINSIC");}
#line 1756 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 182: // GSTRING: "INVERT"
#line 724 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("INVERT");}
#line 1762 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 183: // GSTRING: "INVERTER"
#line 725 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("INVERTER");}
#line 1768 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 184: // GSTRING: "IRDROP"
#line 726 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("IRDROP");}
#line 1774 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 185: // GSTRING: "ITERATE"
#line 727 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("ITERATE");}
#line 1780 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 186: // GSTRING: "IV_TABLES"
#line 728 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("IV_TABLES");}
#line 1786 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 187: // GSTRING: "LAYER"
#line 729 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("LAYER");}
#line 1792 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 188: // GSTRING: "LAYERS"
#line 730 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("LAYERS");}
#line 1798 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 189: // GSTRING: "LEAKAGE"
#line 731 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("LEAKAGE");}
#line 1804 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 190: // GSTRING: "LENGTH"
#line 732 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("LENGTH");}
#line 1810 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 191: // GSTRING: "LENGTHSUM"
#line 733 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("LENGTHSUM");}
#line 1816 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 192: // GSTRING: "LENGTHTHRESHOLD"
#line 734 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("LENGTHTHRESHOLD");}
#line 1822 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 193: // GSTRING: "LEQ"
#line 735 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("LEQ");}
#line 1828 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 194: // GSTRING: "LIBRARY"
#line 736 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("LIBRARY");}
#line 1834 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 195: // GSTRING: "LONGARRAY"
#line 737 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("LONGARRAY");}
#line 1840 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 196: // GSTRING: "MACRO"
#line 738 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("MACRO");}
#line 1846 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 197: // GSTRING: "MANUFACTURINGGRID"
#line 739 "LefParser.yy"
                      {(yylhs.value.stringVal) = new std::string ("MANUFACTURINGGRID");}
#line 1852 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 198: // GSTRING: "MASTERSLICE"
#line 740 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("MASTERSLICE");}
#line 1858 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 199: // GSTRING: "MATCH"
#line 741 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("MATCH");}
#line 1864 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 200: // GSTRING: "MAXADJACENTSLOTSPACING"
#line 742 "LefParser.yy"
                           {(yylhs.value.stringVal) = new std::string ("MAXADJACENTSLOTSPACING");}
#line 1870 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 201: // GSTRING: "MAXCOAXIALSLOTSPACING"
#line 743 "LefParser.yy"
                          {(yylhs.value.stringVal) = new std::string ("MAXCOAXIALSLOTSPACING");}
#line 1876 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 202: // GSTRING: "MAXDELAY"
#line 744 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("MAXDELAY");}
#line 1882 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 203: // GSTRING: "MAXEDGES"
#line 745 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("MAXEDGES");}
#line 1888 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 204: // GSTRING: "MAXEDGESLOTSPACING"
#line 746 "LefParser.yy"
                       {(yylhs.value.stringVal) = new std::string ("MAXEDGESLOTSPACING");}
#line 1894 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 205: // GSTRING: "MAXLOAD"
#line 747 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("MAXLOAD");}
#line 1900 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 206: // GSTRING: "MAXIMUMDENSITY"
#line 748 "LefParser.yy"
                   {(yylhs.value.stringVal) = new std::string ("MAXIMUMDENSITY");}
#line 1906 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 207: // GSTRING: "MAXVIASTACK"
#line 749 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("MAXVIASTACK");}
#line 1912 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 208: // GSTRING: "MAXWIDTH"
#line 750 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("MAXWIDTH");}
#line 1918 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 209: // GSTRING: "MAXXY"
#line 751 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("MAXXY");}
#line 1924 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 210: // GSTRING: "MEGAHERTZ"
#line 752 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("MEGAHERTZ");}
#line 1930 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 211: // GSTRING: "METALOVERHANG"
#line 753 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("METALOVERHANG");}
#line 1936 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 212: // GSTRING: "MICRONS"
#line 754 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("MICRONS");}
#line 1942 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 213: // GSTRING: "MILLIAMPS"
#line 755 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("MILLIAMPS");}
#line 1948 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 214: // GSTRING: "MILLIWATTS"
#line 756 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("MILLIWATTS");}
#line 1954 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 215: // GSTRING: "MINCUTS"
#line 757 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("MINCUTS");}
#line 1960 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 216: // GSTRING: "MINENCLOSEDAREA"
#line 758 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("MINENCLOSEDAREA");}
#line 1966 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 217: // GSTRING: "MINFEATURE"
#line 759 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("MINFEATURE");}
#line 1972 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 218: // GSTRING: "MINIMUMCUT"
#line 760 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("MINIMUMCUT");}
#line 1978 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 219: // GSTRING: "MINIMUMDENSITY"
#line 761 "LefParser.yy"
                   {(yylhs.value.stringVal) = new std::string ("MINIMUMDENSITY");}
#line 1984 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 220: // GSTRING: "MINPINS"
#line 762 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("MINPINS");}
#line 1990 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 221: // GSTRING: "MINSIZE"
#line 763 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("MINSIZE");}
#line 1996 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 222: // GSTRING: "MINSTEP"
#line 764 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("MINSTEP");}
#line 2002 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 223: // GSTRING: "MINWIDTH"
#line 765 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("MINWIDTH");}
#line 2008 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 224: // GSTRING: "MPWH"
#line 766 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("MPWH");}
#line 2014 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 225: // GSTRING: "MPWL"
#line 767 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("MPWL");}
#line 2020 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 226: // GSTRING: "MUSTJOIN"
#line 768 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("MUSTJOIN");}
#line 2026 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 227: // GSTRING: "MX"
#line 769 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("MX");}
#line 2032 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 228: // GSTRING: "MY"
#line 770 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("MY");}
#line 2038 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 229: // GSTRING: "MXR90"
#line 771 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("MXR90");}
#line 2044 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 230: // GSTRING: "MYR90"
#line 772 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("MYR90");}
#line 2050 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 231: // GSTRING: "N"
#line 773 "LefParser.yy"
      {(yylhs.value.stringVal) = new std::string ("N");}
#line 2056 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 232: // GSTRING: "NAMEMAPSTRING"
#line 774 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("NAMEMAPSTRING");}
#line 2062 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 233: // GSTRING: "NAMESCASESENSITIVE"
#line 775 "LefParser.yy"
                       {(yylhs.value.stringVal) = new std::string ("NAMESCASESENSITIVE");}
#line 2068 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 234: // GSTRING: "NANOSECONDS"
#line 776 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("NANOSECONDS");}
#line 2074 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 235: // GSTRING: "NEGEDGE"
#line 777 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("NEGEDGE");}
#line 2080 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 236: // GSTRING: "NETEXPR"
#line 778 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("NETEXPR");}
#line 2086 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 237: // GSTRING: "NETS"
#line 779 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("NETS");}
#line 2092 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 238: // GSTRING: "NEW"
#line 780 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("NEW");}
#line 2098 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 239: // GSTRING: "NONDEFAULTRULE"
#line 781 "LefParser.yy"
                   {(yylhs.value.stringVal) = new std::string ("NONDEFAULTRULE");}
#line 2104 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 240: // GSTRING: "NONE"
#line 782 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("NONE");}
#line 2110 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 241: // GSTRING: "NONINVERT"
#line 783 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("NONINVERT");}
#line 2116 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 242: // GSTRING: "NONUNATE"
#line 784 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("NONUNATE");}
#line 2122 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 243: // GSTRING: "NOISETABLE"
#line 785 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("NOISETABLE");}
#line 2128 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 244: // GSTRING: "NOTCHLENGTH"
#line 786 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("NOTCHLENGTH");}
#line 2134 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 245: // GSTRING: "NOTCHSPACING"
#line 787 "LefParser.yy"
                 {(yylhs.value.stringVal) = new std::string ("NOTCHSPACING");}
#line 2140 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 246: // GSTRING: "NOWIREEXTENSIONATPIN"
#line 788 "LefParser.yy"
                         {(yylhs.value.stringVal) = new std::string ("NOWIREEXTENSIONATPIN");}
#line 2146 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 247: // GSTRING: "OBS"
#line 789 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("OBS");}
#line 2152 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 248: // GSTRING: "OFF"
#line 790 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("OFF");}
#line 2158 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 249: // GSTRING: "OFFSET"
#line 791 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("OFFSET");}
#line 2164 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 250: // GSTRING: "OHMS"
#line 792 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("OHMS");}
#line 2170 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 251: // GSTRING: "ON"
#line 793 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("ON");}
#line 2176 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 252: // GSTRING: "OR"
#line 794 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("OR");}
#line 2182 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 253: // GSTRING: "ORIENT"
#line 795 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("ORIENT");}
#line 2188 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 254: // GSTRING: "ORIENTATION"
#line 796 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("ORIENTATION");}
#line 2194 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 255: // GSTRING: "ORIGIN"
#line 797 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("ORIGIN");}
#line 2200 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 256: // GSTRING: "ORTHOGONAL"
#line 798 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("ORTHOGONAL");}
#line 2206 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 257: // GSTRING: "OUTPUT"
#line 799 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("OUTPUT");}
#line 2212 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 258: // GSTRING: "OUTPUTPINANTENNASIZE"
#line 800 "LefParser.yy"
                         {(yylhs.value.stringVal) = new std::string ("OUTPUTPINANTENNASIZE");}
#line 2218 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 259: // GSTRING: "OUTPUTNOISEMARGIN"
#line 801 "LefParser.yy"
                      {(yylhs.value.stringVal) = new std::string ("OUTPUTNOISEMARGIN");}
#line 2224 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 260: // GSTRING: "OUTPUTRESISTANCE"
#line 802 "LefParser.yy"
                     {(yylhs.value.stringVal) = new std::string ("OUTPUTRESISTANCE");}
#line 2230 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 261: // GSTRING: "OUTSIDECORNER"
#line 803 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("OUTSIDECORNER");}
#line 2236 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 262: // GSTRING: "OVERHANG"
#line 804 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("OVERHANG");}
#line 2242 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 263: // GSTRING: "OVERLAP"
#line 805 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("OVERLAP");}
#line 2248 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 264: // GSTRING: "OVERLAPS"
#line 806 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("OVERLAPS");}
#line 2254 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 265: // GSTRING: "OXIDE1"
#line 807 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("OXIDE1");}
#line 2260 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 266: // GSTRING: "OXIDE2"
#line 808 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("OXIDE2");}
#line 2266 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 267: // GSTRING: "OXIDE3"
#line 809 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("OXIDE3");}
#line 2272 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 268: // GSTRING: "OXIDE4"
#line 810 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("OXIDE4");}
#line 2278 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 269: // GSTRING: "PAD"
#line 811 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("PAD");}
#line 2284 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 270: // GSTRING: "PARALLELEDGE"
#line 812 "LefParser.yy"
                 {(yylhs.value.stringVal) = new std::string ("PARALLELEDGE");}
#line 2290 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 271: // GSTRING: "PARALLELOVERLAP"
#line 813 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("PARALLELOVERLAP");}
#line 2296 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 272: // GSTRING: "PARALLELRUNLENGTH"
#line 814 "LefParser.yy"
                      {(yylhs.value.stringVal) = new std::string ("PARALLELRUNLENGTH");}
#line 2302 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 273: // GSTRING: "PATH"
#line 815 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("PATH");}
#line 2308 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 274: // GSTRING: "PATTERN"
#line 816 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("PATTERN");}
#line 2314 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 275: // GSTRING: "PEAK"
#line 817 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("PEAK");}
#line 2320 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 276: // GSTRING: "PERIOD"
#line 818 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("PERIOD");}
#line 2326 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 277: // GSTRING: "PGONLY"
#line 819 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("PGONLY");}
#line 2332 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 278: // GSTRING: "PICOFARADS"
#line 820 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("PICOFARADS");}
#line 2338 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 279: // GSTRING: "PIN"
#line 821 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("PIN");}
#line 2344 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 280: // GSTRING: "PITCH"
#line 822 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("PITCH");}
#line 2350 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 281: // GSTRING: "PLACED"
#line 823 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("PLACED");}
#line 2356 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 282: // GSTRING: "POLYGON"
#line 824 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("POLYGON");}
#line 2362 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 283: // GSTRING: "PORT"
#line 825 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("PORT");}
#line 2368 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 284: // GSTRING: "POSEDGE"
#line 826 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("POSEDGE");}
#line 2374 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 285: // GSTRING: "POST"
#line 827 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("POST");}
#line 2380 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 286: // GSTRING: "POWER"
#line 828 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("POWER");}
#line 2386 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 287: // GSTRING: "PRE"
#line 829 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("PRE");}
#line 2392 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 288: // GSTRING: "PREFERENCLOSURE"
#line 830 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("PREFERENCLOSURE");}
#line 2398 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 289: // GSTRING: "PRL"
#line 831 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("PRL");}
#line 2404 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 290: // GSTRING: "PROPERTY"
#line 832 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("PROPERTY");}
#line 2410 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 291: // GSTRING: "PROPERTYDEFINITIONS"
#line 833 "LefParser.yy"
                        {(yylhs.value.stringVal) = new std::string ("PROPERTYDEFINITIONS");}
#line 2416 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 292: // GSTRING: "PROTRUSIONWIDTH"
#line 834 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("PROTRUSIONWIDTH");}
#line 2422 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 293: // GSTRING: "PULLDOWNRES"
#line 835 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("PULLDOWNRES");}
#line 2428 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 294: // GSTRING: "PWL"
#line 836 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("PWL");}
#line 2434 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 295: // GSTRING: "R0"
#line 837 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("R0");}
#line 2440 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 296: // GSTRING: "R90"
#line 838 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("R90");}
#line 2446 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 297: // GSTRING: "R180"
#line 839 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("R180");}
#line 2452 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 298: // GSTRING: "R270"
#line 840 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("R270");}
#line 2458 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 299: // GSTRING: "RANGE"
#line 841 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("RANGE");}
#line 2464 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 300: // GSTRING: "REAL"
#line 842 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("REAL");}
#line 2470 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 301: // GSTRING: "RECOVERY"
#line 843 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("RECOVERY");}
#line 2476 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 302: // GSTRING: "RECT"
#line 844 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("RECT");}
#line 2482 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 303: // GSTRING: "RESISTANCE"
#line 845 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("RESISTANCE");}
#line 2488 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 304: // GSTRING: "RESISTIVE"
#line 846 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("RESISTIVE");}
#line 2494 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 305: // GSTRING: "RING"
#line 847 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("RING");}
#line 2500 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 306: // GSTRING: "RISE"
#line 848 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("RISE");}
#line 2506 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 307: // GSTRING: "RISECS"
#line 849 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("RISECS");}
#line 2512 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 308: // GSTRING: "RISERS"
#line 850 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("RISERS");}
#line 2518 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 309: // GSTRING: "RISESATCUR"
#line 851 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("RISESATCUR");}
#line 2524 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 310: // GSTRING: "RISESATT1"
#line 852 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("RISESATT1");}
#line 2530 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 311: // GSTRING: "RISESLEWLIMIT"
#line 853 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("RISESLEWLIMIT");}
#line 2536 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 312: // GSTRING: "RISET0"
#line 854 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("RISET0");}
#line 2542 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 313: // GSTRING: "RISETHRESH"
#line 855 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("RISETHRESH");}
#line 2548 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 314: // GSTRING: "RISEVOLTAGETHRESHOLD"
#line 856 "LefParser.yy"
                         {(yylhs.value.stringVal) = new std::string ("RISEVOLTAGETHRESHOLD");}
#line 2554 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 315: // GSTRING: "RMS"
#line 857 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("RMS");}
#line 2560 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 316: // GSTRING: "ROUTING"
#line 858 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("ROUTING");}
#line 2566 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 317: // GSTRING: "ROWABUTSPACING"
#line 859 "LefParser.yy"
                   {(yylhs.value.stringVal) = new std::string ("ROWABUTSPACING");}
#line 2572 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 318: // GSTRING: "ROWCOL"
#line 860 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("ROWCOL");}
#line 2578 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 319: // GSTRING: "ROWMAJOR"
#line 861 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("ROWMAJOR");}
#line 2584 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 320: // GSTRING: "ROWMINSPACING"
#line 862 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("ROWMINSPACING");}
#line 2590 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 321: // GSTRING: "ROWPATTERN"
#line 863 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("ROWPATTERN");}
#line 2596 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 322: // GSTRING: "RPERSQ"
#line 864 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("RPERSQ");}
#line 2602 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 323: // GSTRING: "S"
#line 865 "LefParser.yy"
      {(yylhs.value.stringVal) = new std::string ("S");}
#line 2608 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 324: // GSTRING: "SAMENET"
#line 866 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("SAMENET");}
#line 2614 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 325: // GSTRING: "SCANUSE"
#line 867 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("SCANUSE");}
#line 2620 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 326: // GSTRING: "SDFCOND"
#line 868 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("SDFCOND");}
#line 2626 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 327: // GSTRING: "SDFCONDEND"
#line 869 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("SDFCONDEND");}
#line 2632 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 328: // GSTRING: "SDFCONDSTART"
#line 870 "LefParser.yy"
                 {(yylhs.value.stringVal) = new std::string ("SDFCONDSTART");}
#line 2638 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 329: // GSTRING: "SETUP"
#line 871 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("SETUP");}
#line 2644 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 330: // GSTRING: "SHAPE"
#line 872 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("SHAPE");}
#line 2650 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 331: // GSTRING: "SHRINKAGE"
#line 873 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("SHRINKAGE");}
#line 2656 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 332: // GSTRING: "SIGNAL"
#line 874 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("SIGNAL");}
#line 2662 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 333: // GSTRING: "SITE"
#line 875 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("SITE");}
#line 2668 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 334: // GSTRING: "SIZE"
#line 876 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("SIZE");}
#line 2674 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 335: // GSTRING: "SKEW"
#line 877 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("SKEW");}
#line 2680 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 336: // GSTRING: "SLOTLENGTH"
#line 878 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("SLOTLENGTH");}
#line 2686 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 337: // GSTRING: "SLOTWIDTH"
#line 879 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("SLOTWIDTH");}
#line 2692 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 338: // GSTRING: "SLOTWIRELENGTH"
#line 880 "LefParser.yy"
                   {(yylhs.value.stringVal) = new std::string ("SLOTWIRELENGTH");}
#line 2698 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 339: // GSTRING: "SLOTWIREWIDTH"
#line 881 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("SLOTWIREWIDTH");}
#line 2704 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 340: // GSTRING: "SPLITWIREWIDTH"
#line 882 "LefParser.yy"
                   {(yylhs.value.stringVal) = new std::string ("SPLITWIREWIDTH");}
#line 2710 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 341: // GSTRING: "SOFT"
#line 883 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("SOFT");}
#line 2716 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 342: // GSTRING: "SOURCE"
#line 884 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("SOURCE");}
#line 2722 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 343: // GSTRING: "SPACER"
#line 885 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("SPACER");}
#line 2728 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 344: // GSTRING: "SPACING"
#line 886 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("SPACING");}
#line 2734 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 345: // GSTRING: "SPACINGTABLE"
#line 887 "LefParser.yy"
                 {(yylhs.value.stringVal) = new std::string ("SPACINGTABLE");}
#line 2740 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 346: // GSTRING: "SPECIALNETS"
#line 888 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("SPECIALNETS");}
#line 2746 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 347: // GSTRING: "STABLE"
#line 889 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("STABLE");}
#line 2752 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 348: // GSTRING: "STACK"
#line 890 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("STACK");}
#line 2758 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 349: // GSTRING: "START"
#line 891 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("START");}
#line 2764 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 350: // GSTRING: "STEP"
#line 892 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("STEP");}
#line 2770 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 351: // GSTRING: "STOP"
#line 893 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("STOP");}
#line 2776 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 352: // GSTRING: "STRING"
#line 894 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("STRING");}
#line 2782 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 353: // GSTRING: "STRUCTURE"
#line 895 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("STRUCTURE");}
#line 2788 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 354: // GSTRING: "SUPPLYSENSITIVITY"
#line 896 "LefParser.yy"
                      {(yylhs.value.stringVal) = new std::string ("SUPPLYSENSITIVITY");}
#line 2794 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 355: // GSTRING: "SYMMETRY"
#line 897 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("SYMMETRY");}
#line 2800 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 356: // GSTRING: "TABLE"
#line 898 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("TABLE");}
#line 2806 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 357: // GSTRING: "TABLEAXIS"
#line 899 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("TABLEAXIS");}
#line 2812 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 358: // GSTRING: "TABLEDIMENSION"
#line 900 "LefParser.yy"
                   {(yylhs.value.stringVal) = new std::string ("TABLEDIMENSION");}
#line 2818 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 359: // GSTRING: "TABLEENTRIES"
#line 901 "LefParser.yy"
                 {(yylhs.value.stringVal) = new std::string ("TABLEENTRIES");}
#line 2824 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 360: // GSTRING: "TAPERRULE"
#line 902 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("TAPERRULE");}
#line 2830 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 361: // GSTRING: "THEN"
#line 903 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("THEN");}
#line 2836 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 362: // GSTRING: "THICKNESS"
#line 904 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("THICKNESS");}
#line 2842 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 363: // GSTRING: "TIEHIGH"
#line 905 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("TIEHIGH");}
#line 2848 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 364: // GSTRING: "TIELOW"
#line 906 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("TIELOW");}
#line 2854 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 365: // GSTRING: "TIEOFFR"
#line 907 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("TIEOFFR");}
#line 2860 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 366: // GSTRING: "TIME"
#line 908 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("TIME");}
#line 2866 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 367: // GSTRING: "TIMING"
#line 909 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("TIMING");}
#line 2872 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 368: // GSTRING: "TO"
#line 910 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("TO");}
#line 2878 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 369: // GSTRING: "TOPIN"
#line 911 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("TOPIN");}
#line 2884 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 370: // GSTRING: "TOPLEFT"
#line 912 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("TOPLEFT");}
#line 2890 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 371: // GSTRING: "TOPOFSTACKONLY"
#line 913 "LefParser.yy"
                   {(yylhs.value.stringVal) = new std::string ("TOPOFSTACKONLY");}
#line 2896 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 372: // GSTRING: "TOPRIGHT"
#line 914 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("TOPRIGHT");}
#line 2902 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 373: // GSTRING: "TRACKS"
#line 915 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("TRACKS");}
#line 2908 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 374: // GSTRING: "TRANSITIONTIME"
#line 916 "LefParser.yy"
                   {(yylhs.value.stringVal) = new std::string ("TRANSITIONTIME");}
#line 2914 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 375: // GSTRING: "TRISTATE"
#line 917 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("TRISTATE");}
#line 2920 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 376: // GSTRING: "TRUE"
#line 918 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("TRUE");}
#line 2926 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 377: // GSTRING: "TWOEDGES"
#line 919 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("TWOEDGES");}
#line 2932 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 378: // GSTRING: "TWOWIDTHS"
#line 920 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("TWOWIDTHS");}
#line 2938 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 379: // GSTRING: "TYPE"
#line 921 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("TYPE");}
#line 2944 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 380: // GSTRING: "UNATENESS"
#line 922 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("UNATENESS");}
#line 2950 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 381: // GSTRING: "UNITS"
#line 923 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("UNITS");}
#line 2956 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 382: // GSTRING: "UNIVERSALNOISEMARGIN"
#line 924 "LefParser.yy"
                         {(yylhs.value.stringVal) = new std::string ("UNIVERSALNOISEMARGIN");}
#line 2962 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 383: // GSTRING: "USE"
#line 925 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("USE");}
#line 2968 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 384: // GSTRING: "USELENGTHTHRESHOLD"
#line 926 "LefParser.yy"
                       {(yylhs.value.stringVal) = new std::string ("USELENGTHTHRESHOLD");}
#line 2974 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 385: // GSTRING: "USEMINSPACING"
#line 927 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("USEMINSPACING");}
#line 2980 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 386: // GSTRING: "USER"
#line 928 "LefParser.yy"
         {(yylhs.value.stringVal) = new std::string ("USER");}
#line 2986 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 387: // GSTRING: "USEVIA"
#line 929 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("USEVIA");}
#line 2992 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 388: // GSTRING: "USEVIARULE"
#line 930 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("USEVIARULE");}
#line 2998 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 389: // GSTRING: "VARIABLE"
#line 931 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("VARIABLE");}
#line 3004 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 390: // GSTRING: "VERSION"
#line 932 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("VERSION");}
#line 3010 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 391: // GSTRING: "VERTICAL"
#line 933 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("VERTICAL");}
#line 3016 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 392: // GSTRING: "VHI"
#line 934 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("VHI");}
#line 3022 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 393: // GSTRING: "VIA"
#line 935 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("VIA");}
#line 3028 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 394: // GSTRING: "VIARULE"
#line 936 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("VIARULE");}
#line 3034 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 395: // GSTRING: "VICTIMLENGTH"
#line 937 "LefParser.yy"
                 {(yylhs.value.stringVal) = new std::string ("VICTIMLENGTH");}
#line 3040 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 396: // GSTRING: "VICTIMNOISE"
#line 938 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("VICTIMNOISE");}
#line 3046 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 397: // GSTRING: "VIRTUAL"
#line 939 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("VIRTUAL");}
#line 3052 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 398: // GSTRING: "VLO"
#line 940 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("VLO");}
#line 3058 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 399: // GSTRING: "VOLTAGE"
#line 941 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("VOLTAGE");}
#line 3064 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 400: // GSTRING: "VOLTS"
#line 942 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("VOLTS");}
#line 3070 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 401: // GSTRING: "W"
#line 943 "LefParser.yy"
      {(yylhs.value.stringVal) = new std::string ("W");}
#line 3076 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 402: // GSTRING: "WELLTAP"
#line 944 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("WELLTAP");}
#line 3082 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 403: // GSTRING: "WIDTH"
#line 945 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("WIDTH");}
#line 3088 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 404: // GSTRING: "WITHIN"
#line 946 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("WITHIN");}
#line 3094 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 405: // GSTRING: "WIRECAP"
#line 947 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("WIRECAP");}
#line 3100 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 406: // GSTRING: "WIREEXTENSION"
#line 948 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("WIREEXTENSION");}
#line 3106 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 407: // GSTRING: "X"
#line 949 "LefParser.yy"
      {(yylhs.value.stringVal) = new std::string ("X");}
#line 3112 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 408: // GSTRING: "Y"
#line 950 "LefParser.yy"
      {(yylhs.value.stringVal) = new std::string ("Y");}
#line 3118 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 409: // GSTRING: K_EQ
#line 951 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("EQ");}
#line 3124 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 410: // GSTRING: K_NE
#line 952 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("NE");}
#line 3130 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 411: // GSTRING: K_LE
#line 953 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("LE");}
#line 3136 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 412: // GSTRING: K_LT
#line 954 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("LT");}
#line 3142 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 413: // GSTRING: K_GE
#line 955 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("GE");}
#line 3148 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 414: // GSTRING: K_GT
#line 956 "LefParser.yy"
       {(yylhs.value.stringVal) = new std::string ("GT");}
#line 3154 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 415: // GSTRING: K_NOT
#line 957 "LefParser.yy"
        {(yylhs.value.stringVal) = new std::string ("NOT");}
#line 3160 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 416: // lef_file: rules extension_opt end_library
#line 960 "LefParser.yy"
      {
        // 11/16/2001 - Wanda da Rosa - pcr 408334
        // Return 1 if there are errors
        if (driver.lef_errors)
           return 1;
        if (!driver.hasVer && driver.versionNum < 5.6)
           driver.lefWarning(2001, "VERSION is a required statement on LEF file with version 5.5 and earlier.\nWithout VERSION defined, the LEF file is technically illegal.\nRefer the LEF/DEF 5.5 and earlier Language Reference manual on how to defined this statement.");
        //only pre 5.6, 5.6 it is obsolete
        if (!driver.hasNameCase && driver.versionNum < 5.6)
           driver.lefWarning(2002, "NAMESCASESENSITIVE is a required statement on LEF file with version 5.5 and earlier.\nWithout NAMESCASESENSITIVE defined, the LEF file is technically illegal.\nRefer the LEF/DEF 5.5 and earlier Language Referece manual on how to define this statement.");
        if (!driver.hasBusBit && driver.versionNum < 5.6)
           driver.lefWarning(2003, "BUSBITCHARS is a required statement on LEF file with version 5.5 and earlier.\nWithout BUSBITCHARS defined, the LEF file is technically illegal.\nRefer the LEF/DEF 5.5 and earlier Language Referece manual on how to define this statement.");
        if (!driver.hasDivChar && driver.versionNum < 5.6)
           driver.lefWarning(2004, "DIVIDERCHAR is a required statementon LEF file with version 5.5 and earlier.\nWithout DIVIDECHAR defined, the LEF file is technically illegal.\nRefer the LEF/DEF 5.5 and earlier Language Referece manual on how to define this statement.");

       driver.resetVars();
      }
#line 3182 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 417: // $@1: %empty
#line 978 "LefParser.yy"
                   { driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 3188 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 418: // version: "VERSION" $@1 "string" ';'
#line 979 "LefParser.yy"
      { 
         driver.versionNum = driver.convert_name2num((*(yystack_[1].value.stringVal)).c_str());
         if (driver.versionNum > 5.7) {
            char temp[120];
            sprintf(temp,
               "Lef parser 5.7 does not support lef file with version %s. Parser stops executions.", (*(yystack_[1].value.stringVal)).c_str());
            driver.lefError(1503, temp);
            //return 1;
         }
/*
         driver.versionNum = $3;         Save the version number for future use */
		 driver.lefrVersionStrCbk(*(yystack_[1].value.stringVal));
		 driver.lefrVersionCbk(driver.versionNum);
         if (driver.versionNum > 5.3 && driver.versionNum < 5.4) {
            driver.ignoreVersion = 1;
         }
         driver.use5_3 = driver.use5_4 = 0;
         driver.lef_errors = 0;
         if (driver.hasVer)     /* More than 1 lef file within the open file */
            driver.resetVars();
         driver.hasVer = 1;
         if (driver.versionNum < 5.6) {
            driver.doneLib = 0;
            driver.lefNamesCaseSensitive = driver.lefReaderCaseSensitive;
         } else {
            driver.doneLib = 1;
            driver.lefNamesCaseSensitive = 1;
         }
         delete (yystack_[1].value.stringVal);
      }
#line 3223 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 419: // $@2: %empty
#line 1009 "LefParser.yy"
                      { driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 3229 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 420: // version: "VERSION" $@2 "double" ';'
#line 1010 "LefParser.yy"
      { 
         driver.versionNum = (yystack_[1].value.doubleVal);
         if (driver.versionNum > 5.7) {
            char temp[120];
            sprintf(temp,
               "Lef parser 5.7 does not support lef file with version %f. Parser stops executions.", (yystack_[1].value.doubleVal));
            driver.lefError(1503, temp);
            //return 1;
         }
/*
         driver.versionNum = $3;         Save the version number for future use */
		 char tmp[20];
		 snprintf(tmp, sizeof(tmp), "%g", (yystack_[1].value.doubleVal));
		 driver.lefrVersionStrCbk(tmp);
		 driver.lefrVersionCbk(driver.versionNum);
         if (driver.versionNum > 5.3 && driver.versionNum < 5.4) {
            driver.ignoreVersion = 1;
         }
         driver.use5_3 = driver.use5_4 = 0;
         driver.lef_errors = 0;
         if (driver.hasVer)     /* More than 1 lef file within the open file */
            driver.resetVars();
         driver.hasVer = 1;
         if (driver.versionNum < 5.6) {
            driver.doneLib = 0;
            driver.lefNamesCaseSensitive = driver.lefReaderCaseSensitive;
         } else {
            driver.doneLib = 1;
            driver.lefNamesCaseSensitive = 1;
         }
      }
#line 3265 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 421: // dividerchar: "DIVIDERCHAR" "qstring" ';'
#line 1044 "LefParser.yy"
      {
		  if (!(*(yystack_[1].value.qstringVal)).empty()) 
		  {
			  driver.lefrDividerCharCbk(*(yystack_[1].value.qstringVal));
		  }
		  else 
		  {
			  driver.lefrDividerCharCbk("/");
			  driver.lefWarning(2005, "DIVIDERCHAR has an invalid null value. Value is set to default /");
		  }
		  driver.hasDivChar = 1;
          delete (yystack_[1].value.qstringVal);
      }
#line 3283 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 422: // busbitchars: "BUSBITCHARS" "qstring" ';'
#line 1059 "LefParser.yy"
   {
		if (!(*(yystack_[1].value.qstringVal)).empty()) 
	   {
		   driver.lefrBusBitCharsCbk(*(yystack_[1].value.qstringVal));
	   } 
	  else 
	  {
		  driver.lefrBusBitCharsCbk("[]");
		 driver.lefWarning(2006, "BUSBITCHAR has an invalid null value. Value is set to default []");
	  }
	driver.hasBusBit = 1;
    delete (yystack_[1].value.qstringVal);
  }
#line 3301 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 425: // rules: error
#line 1076 "LefParser.yy"
            { }
#line 3307 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 426: // end_library: %empty
#line 1079 "LefParser.yy"
      {
        if (driver.versionNum >= 5.6) {
           driver.doneLib = 1;
           driver.ge56done = 1;
        }
      }
#line 3318 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 427: // end_library: "END" "LIBRARY"
#line 1086 "LefParser.yy"
      {
        driver.doneLib = 1;
        driver.ge56done = 1;
		driver.lefrLibraryEndCbk();
        // 11/16/2001 - Wanda da Rosa - pcr 408334
        // Return 1 if there are errors
/*
        if (driver.lef_errors)
           return 1;
        if (!driver.hasVer)
           driver.lefWarning(2001, "VERSION is a required statement.");
        if (!driver.hasNameCase)
           driver.lefWarning(2002, "NAMESCASESENSITIVE is a required statement.");
        if (!driver.hasBusBit && driver.versionNum < 5.6)
           driver.lefWarning(2003, "BUSBITCHARS is a required statement.");
        if (!driver.hasDivChar && driver.versionNum < 5.6)
           driver.lefWarning(2004, "DIVIDERCHAR is a required statement.");
        driver.hasVer = 0;
        driver.hasNameCase = 0;
        driver.hasBusBit = 0;
        driver.hasDivChar = 0;
        driver.hasManufactur = 0;
        driver.hasMinfeature = 0;
        driver.hasSite = 0;
*/
      }
#line 3349 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 466: // case_sensitivity: "NAMESCASESENSITIVE" "ON" ';'
#line 1128 "LefParser.yy"
          {
            if (driver.versionNum < 5.6) {
              driver.lefNamesCaseSensitive = TRUE;
			  driver.lefrCaseSensitiveCbk(driver.lefNamesCaseSensitive);
              driver.hasNameCase = 1;
            } 
			else if (driver.caseSensitiveWarnings++ < driver.lefrCaseSensitiveWarnings)
				driver.lefWarning(2007, "NAMESCASESENSITIVE statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
	  }
#line 3363 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 467: // case_sensitivity: "NAMESCASESENSITIVE" "OFF" ';'
#line 1138 "LefParser.yy"
          {
            if (driver.versionNum < 5.6) 
			{
			driver.lefNamesCaseSensitive = FALSE;
			driver.lefrCaseSensitiveCbk(driver.lefNamesCaseSensitive);
			driver.hasNameCase = 1;
            } 
			else 
			{
                if (driver.caseSensitiveWarnings++ < driver.lefrCaseSensitiveWarnings) 
				{
                  driver.lefError(1504, "NAMESCASESENSITIVE statement is set with OFF.\nStarting version 5.6, NAMESCASENSITIVE is obsolete,\nif it is defined, it has to have the ON value.\nParser stops executions.");
              }
            }
	  }
#line 3383 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 468: // wireextension: "NOWIREEXTENSIONATPIN" "ON" ';'
#line 1155 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) 
	  {
	  driver.lefrNoWireExtensionCbk("ON");
      } 
	  else if (driver.noWireExtensionWarnings++ < driver.lefrNoWireExtensionWarnings)
	  {
	  driver.lefWarning(2008, "NOWIREEXTENSIONATPIN statement is obsolete in version 5.6 or later.\nThe NOWIREEXTENSIONATPIN statement will be ignored.");
	  }
    }
#line 3398 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 469: // wireextension: "NOWIREEXTENSIONATPIN" "OFF" ';'
#line 1166 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) 
	  {
	  driver.lefrNoWireExtensionCbk("OFF");
      } 
	  else if (driver.noWireExtensionWarnings++ < driver.lefrNoWireExtensionWarnings)
	  {
	  driver.lefWarning(2008, "NOWIREEXTENSIONATPIN statement is obsolete in version 5.6 or later.\nThe NOWIREEXTENSIONATPIN statement will be ignored.");
	  }
    }
#line 3413 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 470: // manufacturing: "MANUFACTURINGGRID" NUMBER ';'
#line 1178 "LefParser.yy"
    {
	driver.lefrManufacturingCbk((yystack_[1].value.doubleVal));
	  driver.hasManufactur = 1;
    }
#line 3422 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 471: // useminspacing: "USEMINSPACING" spacing_type spacing_value ';'
#line 1184 "LefParser.yy"
  {
    if ((strcmp((*(yystack_[2].value.stringVal)).c_str(), "PIN") == 0) && (driver.versionNum >= 5.6)) 
	{
         if (driver.useMinSpacingWarnings++ < driver.lefrUseMinSpacingWarnings)
            driver.lefWarning(2009, "USEMINSPACING PIN statement is obsolete in version 5.6 or later.\n The USEMINSPACING PIN statement will be ignored.");
    } 
	else 
	{
		driver.lefrUseMinSpacing.lefiUseMinSpacing::set((*(yystack_[2].value.stringVal)).c_str(), (yystack_[1].value.integerVal));
		driver.lefrUseMinSpacingCbk(driver.lefrUseMinSpacing);
    }
    delete (yystack_[2].value.stringVal);
  }
#line 3440 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 472: // clearancemeasure: "CLEARANCEMEASURE" clearance_type ';'
#line 1199 "LefParser.yy"
    {driver.lefrClearanceMeasureCbk(*(yystack_[1].value.stringVal)); delete (yystack_[1].value.stringVal);}
#line 3446 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 473: // clearance_type: "MAXXY"
#line 1202 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("MAXXY");}
#line 3452 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 474: // clearance_type: "EUCLIDEAN"
#line 1203 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("EUCLIDEAN");}
#line 3458 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 475: // spacing_type: "OBS"
#line 1206 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("OBS");}
#line 3464 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 476: // spacing_type: "PIN"
#line 1207 "LefParser.yy"
            {(yylhs.value.stringVal) = new std::string ("PIN");}
#line 3470 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 477: // spacing_value: "ON"
#line 1210 "LefParser.yy"
            {(yylhs.value.integerVal) = 1;}
#line 3476 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 478: // spacing_value: "OFF"
#line 1211 "LefParser.yy"
            {(yylhs.value.integerVal) = 0;}
#line 3482 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 479: // units_section: start_units units_rules "END" "UNITS"
#line 1214 "LefParser.yy"
    { driver.lefrUnitsCbk(driver.lefrUnits);}
#line 3488 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 480: // start_units: "UNITS"
#line 1217 "LefParser.yy"
    {
      driver.lefrUnits.lefiUnits::clear();
      if (driver.hasManufactur) 
	  {
        if (driver.unitsWarnings++ < driver.lefrUnitsWarnings) 
		{driver.lefError(1505, "MANUFACTURINGGRID statement was defined before UNITS.\nRefer the LEF Language Reference manual for the order of LEF statements.");}
      }
      if (driver.hasMinfeature) {
        if (driver.unitsWarnings++ < driver.lefrUnitsWarnings) 
		{
          driver.lefError(1712, "MINFEATURE statement was defined before UNITS.\nRefer the LEF Language Reference manual for the order of LEF statements.");
        }
      }
      if (driver.versionNum < 5.6) {
        if (driver.hasSite) {/*SITE is defined before UNIT and is illegal in pre 5.6*/
          driver.lefError(1713, "SITE statement was defined before UNITS.\nRefer the LEF Language Reference manual for the order of LEF statements.");
        }
      }
    }
#line 3512 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 483: // units_rule: "TIME" "NANOSECONDS" NUMBER ';'
#line 1242 "LefParser.yy"
    { driver.lefrUnits.lefiUnits::setTime((yystack_[1].value.doubleVal)); }
#line 3518 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 484: // units_rule: "CAPACITANCE" "PICOFARADS" NUMBER ';'
#line 1244 "LefParser.yy"
    { driver.lefrUnits.lefiUnits::setCapacitance((yystack_[1].value.doubleVal)); }
#line 3524 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 485: // units_rule: "RESISTANCE" "OHMS" NUMBER ';'
#line 1246 "LefParser.yy"
    { driver.lefrUnits.lefiUnits::setResistance((yystack_[1].value.doubleVal)); }
#line 3530 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 486: // units_rule: "POWER" "MILLIWATTS" NUMBER ';'
#line 1248 "LefParser.yy"
    { driver.lefrUnits.lefiUnits::setPower((yystack_[1].value.doubleVal)); }
#line 3536 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 487: // units_rule: "CURRENT" "MILLIAMPS" NUMBER ';'
#line 1250 "LefParser.yy"
    { driver.lefrUnits.lefiUnits::setCurrent((yystack_[1].value.doubleVal)); }
#line 3542 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 488: // units_rule: "VOLTAGE" "VOLTS" NUMBER ';'
#line 1252 "LefParser.yy"
    { driver.lefrUnits.lefiUnits::setVoltage((yystack_[1].value.doubleVal)); }
#line 3548 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 489: // units_rule: "DATABASE" "MICRONS" NUMBER ';'
#line 1254 "LefParser.yy"
    { 
      if(driver.validNum((int)(yystack_[1].value.doubleVal))) {
            driver.lefrUnits.lefiUnits::setDatabase("MICRONS", (yystack_[1].value.doubleVal));
      }
    }
#line 3558 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 490: // units_rule: "FREQUENCY" "MEGAHERTZ" NUMBER ';'
#line 1260 "LefParser.yy"
    { driver.lefrUnits.lefiUnits::setFrequency((yystack_[1].value.doubleVal)); }
#line 3564 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 491: // layer_rule: start_layer layer_options end_layer
#line 1264 "LefParser.yy"
    { driver.lefrLayerCbk(driver.lefrLayer);}
#line 3570 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 492: // $@3: %empty
#line 1266 "LefParser.yy"
                     {driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 3576 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 493: // start_layer: "LAYER" $@3 GSTRING
#line 1267 "LefParser.yy"
    { 
      if (driver.lefrHasMaxVS) {   /* 5.5 */
          if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
            driver.lefError(1506, "A MAXVIASTACK statement is defined before the LAYER statement.\nRefer to the LEF Language Reference manual for the order of LEF statements.");
          }
      }
	  driver.lefrLayer.lefiLayer::setName((*(yystack_[0].value.stringVal)).c_str());
      driver.useLenThr = 0;
      driver.layerCut = 0;
      driver.layerMastOver = 0;
      driver.layerRout = 0;
      driver.layerDir = 0;
      driver.lefrHasLayer = 1;
      driver.layerName = *(yystack_[0].value.stringVal);
      driver.hasType = 0;
      driver.hasPitch = 0;
      driver.hasWidth = 0;
      driver.hasDirection = 0;
      driver.hasParallel = 0;
      driver.hasInfluence = 0;
      driver.hasTwoWidths = 0;
      driver.lefrHasSpacingTbl = 0;
      driver.lefrHasSpacing = 0;
      delete (yystack_[0].value.stringVal);
    }
#line 3606 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 494: // $@4: %empty
#line 1293 "LefParser.yy"
                 {driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 3612 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 495: // end_layer: "END" $@4 GSTRING
#line 1294 "LefParser.yy"
    { 
      if (driver.layerName != *(yystack_[0].value.stringVal)) {
          if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
             driver.outMsg = (char*)lefMalloc(10000);
             sprintf (driver.outMsg,
                "END LAYER name %s is different from the LAYER name %s.\nCorrect the LEF file before rerun it through the LEF parser.", (*(yystack_[0].value.stringVal)).c_str(), driver.layerName.c_str());
             driver.lefError(1507, driver.outMsg);
             lefFree(driver.outMsg);
          } 
      } 
      if (!driver.lefrRelaxMode) {
        if (driver.hasType == 0) {
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.lefError(1508, "TYPE statement is a required statement in a LAYER and it is not defined.");
            }
        }
        if ((driver.layerRout == 1) && (driver.hasPitch == 0)) {
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1509, "PITCH statement is a required statement in a LAYER with type ROUTING and it is not defined.");
          }
        }
        if ((driver.layerRout == 1) && (driver.hasWidth == 0)) {
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1510, "WIDTH statement is a required statement in a LAYER with type ROUTING and it is not defined.");
          }
        }
        if ((driver.layerRout == 1) && (driver.hasDirection == 0)) {
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg, "The DIRECTION statement which is required in a LAYER with TYPE ROUTING is not defined in LAYER %s.\nUpdate your lef file and add the DIRECTION statement for layer %s.", (*(yystack_[0].value.stringVal)).c_str(), (*(yystack_[0].value.stringVal)).c_str());
              driver.lefError(1511, driver.outMsg);
              lefFree(driver.outMsg);
            }
        }
      }
      delete (yystack_[0].value.stringVal);
    }
#line 3654 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 496: // layer_options: %empty
#line 1334 "LefParser.yy"
    { }
#line 3660 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 497: // layer_options: layer_options layer_option
#line 1339 "LefParser.yy"
    { }
#line 3666 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 498: // $@5: %empty
#line 1343 "LefParser.yy"
    {
       /* let setArraySpacingCutSpacing to set the data */
    }
#line 3674 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 499: // $@6: %empty
#line 1349 "LefParser.yy"
    {
         driver.lefrLayer.lefiLayer::setArraySpacingCut((yystack_[4].value.doubleVal));
         driver.lefrLayer.lefiLayer::addArraySpacingArray((int)(yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
         driver.arrayCutsVal = (int)(yystack_[2].value.doubleVal);  /* set the value */
         driver.arrayCutsWar = 0;
    }
#line 3685 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 500: // layer_option: "ARRAYSPACING" $@5 layer_arraySpacing_long layer_arraySpacing_width "CUTSPACING" NUMBER "ARRAYCUTS" NUMBER "SPACING" NUMBER $@6 layer_arraySpacing_arraycuts ';'
#line 1356 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
         driver.outMsg = (char*)lefMalloc(10000);
         sprintf(driver.outMsg,
           "ARRAYSPACING is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
         driver.lefError(1685, driver.outMsg);
         lefFree(driver.outMsg);
      }
    }
#line 3699 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 501: // layer_option: "TYPE" layer_type ';'
#line 1366 "LefParser.yy"
    {
         driver.lefrLayer.lefiLayer::setType((*(yystack_[1].value.stringVal)).c_str());
      driver.hasType = 1;
      delete (yystack_[1].value.stringVal);
    }
#line 3709 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 502: // layer_option: "PITCH" NUMBER ';'
#line 1372 "LefParser.yy"
    { 
      driver.lefrLayer.lefiLayer::setPitch((yystack_[1].value.doubleVal));
      driver.hasPitch = 1;  
    }
#line 3718 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 503: // layer_option: "PITCH" NUMBER NUMBER ';'
#line 1377 "LefParser.yy"
    { 
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setPitchXY((yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
      driver.hasPitch = 1;  
    }
#line 3727 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 504: // layer_option: "DIAGPITCH" NUMBER ';'
#line 1382 "LefParser.yy"
    { 
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setDiagPitch((yystack_[1].value.doubleVal));
    }
#line 3735 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 505: // layer_option: "DIAGPITCH" NUMBER NUMBER ';'
#line 1386 "LefParser.yy"
    { 
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setDiagPitchXY((yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
    }
#line 3743 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 506: // layer_option: "OFFSET" NUMBER ';'
#line 1390 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setOffset((yystack_[1].value.doubleVal));
    }
#line 3751 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 507: // layer_option: "OFFSET" NUMBER NUMBER ';'
#line 1394 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setOffsetXY((yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
    }
#line 3759 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 508: // layer_option: "DIAGWIDTH" NUMBER ';'
#line 1398 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setDiagWidth((yystack_[1].value.doubleVal));
    }
#line 3767 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 509: // layer_option: "DIAGSPACING" NUMBER ';'
#line 1402 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setDiagSpacing((yystack_[1].value.doubleVal));
    }
#line 3775 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 510: // layer_option: "WIDTH" NUMBER ';'
#line 1406 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setWidth((yystack_[1].value.doubleVal));
      driver.hasWidth = 1;  
    }
#line 3784 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 511: // layer_option: "AREA" NUMBER ';'
#line 1411 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setArea((yystack_[1].value.doubleVal));
    }
#line 3792 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 512: // $@7: %empty
#line 1415 "LefParser.yy"
    {
      driver.hasSpCenter = 0;       /* reset to 0, only once per spacing is allowed */
      driver.hasSpSamenet = 0;
      driver.hasSpParallel = 0;
      driver.hasSpLayer = 0;
      driver.layerCutSpacing = (yystack_[0].value.doubleVal);  /* for error message purpose */
      // 11/22/99 - Wanda da Rosa, PCR 283762
      //            Issue an error is this is defined in masterslice
      if (driver.layerMastOver) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1512, "It is illegal to define a SPACING statement in LAYER with TYPE MASTERSLICE or OVERLAP. Parser stops executions.");
            }
         }
      }
      // 5.5 either SPACING or SPACINGTABLE, not both for routing layer only
      if (driver.layerRout) {
        if (driver.lefrHasSpacingTbl && driver.versionNum < 5.7) {
           if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
              if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                driver.lefWarning(2010, "It is illegal to have both SPACING rules & SPACINGTABLE rules within a ROUTING layer");
              }
           }
        }
        if (/*driver.lefrLayerCbk*/ 1)
           driver.lefrLayer.lefiLayer::setSpacingMin((yystack_[0].value.doubleVal));
        driver.lefrHasSpacing = 1;
      } else { 
        if (/*driver.lefrLayerCbk*/ 1)
           driver.lefrLayer.lefiLayer::setSpacingMin((yystack_[0].value.doubleVal));
      }
    }
#line 3829 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 514: // $@8: %empty
#line 1450 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
         driver.lefrLayer.lefiLayer::setSpacingTableOrtho();
      if (/*driver.lefrLayerCbk*/ 1) /* due to converting to C, else, convertor produce */
         driver.lefrLayer.lefiLayer::addSpacingTableOrthoWithin((yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));/*bad code*/
    }
#line 3840 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 515: // layer_option: "SPACINGTABLE" "ORTHOGONAL" "WITHIN" NUMBER "SPACING" NUMBER $@8 layer_spacingtable_opts ';'
#line 1457 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
         driver.outMsg = (char*)lefMalloc(10000);
         sprintf(driver.outMsg,
           "SPACINGTABLE ORTHOGONAL is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
         driver.lefError(1694, driver.outMsg);
         lefFree(driver.outMsg);
      }
    }
#line 3854 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 516: // layer_option: "DIRECTION" layer_direction ';'
#line 1467 "LefParser.yy"
    {
      driver.layerDir = 1;
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1513, "DIRECTION statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setDirection(*(yystack_[1].value.stringVal));
      driver.hasDirection = 1;  
      delete (yystack_[1].value.stringVal);
    }
#line 3873 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 517: // layer_option: "RESISTANCE" "RPERSQ" NUMBER ';'
#line 1482 "LefParser.yy"
    {
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1514, "RESISTANCE statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setResistance((yystack_[1].value.doubleVal));
    }
#line 3889 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 518: // layer_option: "RESISTANCE" "RPERSQ" "PWL" '(' res_points ')' ';'
#line 1494 "LefParser.yy"
    {
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1515, "RESISTANCE statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
    }
#line 3904 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 519: // layer_option: "CAPACITANCE" "CPERSQDIST" NUMBER ';'
#line 1505 "LefParser.yy"
    {
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1516, "CAPACITANCE statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setCapacitance((yystack_[1].value.doubleVal));
    }
#line 3920 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 520: // layer_option: "CAPACITANCE" "CPERSQDIST" "PWL" '(' cap_points ')' ';'
#line 1517 "LefParser.yy"
    {
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1517, "CAPACITANCE statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
    }
#line 3935 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 521: // layer_option: "HEIGHT" NUMBER ';'
#line 1528 "LefParser.yy"
    {
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1518, "HEIGHT statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setHeight((yystack_[1].value.doubleVal));
    }
#line 3951 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 522: // layer_option: "WIREEXTENSION" NUMBER ';'
#line 1540 "LefParser.yy"
    {
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1519, "WIREEXTENSION statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setWireExtension((yystack_[1].value.doubleVal));
    }
#line 3967 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 523: // layer_option: "THICKNESS" NUMBER ';'
#line 1552 "LefParser.yy"
    {
      if (!driver.layerRout && (driver.layerCut || driver.layerMastOver)) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1520, "THICKNESS statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setThickness((yystack_[1].value.doubleVal));
    }
#line 3983 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 524: // layer_option: "SHRINKAGE" NUMBER ';'
#line 1564 "LefParser.yy"
    {
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1521, "SHRINKAGE statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setShrinkage((yystack_[1].value.doubleVal));
    }
#line 3999 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 525: // layer_option: "CAPMULTIPLIER" NUMBER ';'
#line 1576 "LefParser.yy"
    {
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1522, "CAPMULTIPLIER statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setCapMultiplier((yystack_[1].value.doubleVal));
    }
#line 4015 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 526: // layer_option: "EDGECAPACITANCE" NUMBER ';'
#line 1588 "LefParser.yy"
    {
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1523, "EDGECAPACITANCE statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setEdgeCap((yystack_[1].value.doubleVal));
    }
#line 4031 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 527: // layer_option: "ANTENNALENGTHFACTOR" NUMBER ';'
#line 1613 "LefParser.yy"
    { /* 5.3 syntax */
      driver.use5_3 = 1;
      if (!driver.layerRout && (driver.layerCut || driver.layerMastOver)) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1525, "ANTENNALENGTHFACTOR statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      } else if (driver.versionNum >= 5.4) {
         if (driver.use5_4) {
            if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
               if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                  driver.outMsg = (char*)lefMalloc(10000);
                  sprintf (driver.outMsg,
                    "ANTENNALENGTHFACTOR statement is a version 5.3 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNALENGTHFACTOR syntax, which is illegal.", driver.versionNum);
                  driver.lefError(1526, driver.outMsg);
                  lefFree(driver.outMsg);
                  /*CHKERR();*/
               }
            }
         }
      }

      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setAntennaLength((yystack_[1].value.doubleVal));
    }
#line 4062 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 528: // layer_option: "CURRENTDEN" NUMBER ';'
#line 1640 "LefParser.yy"
    {
      if (driver.versionNum < 5.2) {
         if (!driver.layerRout) {
            if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
               if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                 driver.lefError(1526, "CURRENTDEN statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
                 /*CHKERR();*/
               }
            }
         }
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setCurrentDensity((yystack_[1].value.doubleVal));
      } else {
         if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
            driver.lefWarning(2079, "CURRENTDEN statement is obsolete in version 5.2 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.2 or later.");
            /*CHKERR();*/
         }
      }
    }
#line 4085 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 529: // layer_option: "CURRENTDEN" "PWL" '(' current_density_pwl_list ')' ';'
#line 1659 "LefParser.yy"
    { 
      if (driver.versionNum < 5.2) {
         if (!driver.layerRout) {
            if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
               if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                 driver.lefError(1526, "CURRENTDEN statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
                 /*CHKERR();*/
               }
            }
         }
      } else {
         if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
            driver.lefWarning(2079, "CURRENTDEN statement is obsolete in version 5.2 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.2 or later.");
            /*CHKERR();*/
         }
      }
    }
#line 4107 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 530: // layer_option: "CURRENTDEN" '(' NUMBER NUMBER ')' ';'
#line 1677 "LefParser.yy"
    {
      if (driver.versionNum < 5.2) {
         if (!driver.layerRout) {
            if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
               if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                 driver.lefError(1526, "CURRENTDEN statement can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
                 /*CHKERR();*/
               }
            }
         }
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setCurrentPoint((yystack_[3].value.doubleVal), (yystack_[2].value.doubleVal));
      } else {
         if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
            driver.lefWarning(2079, "CURRENTDEN statement is obsolete in version 5.2 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.2 or later.");
            /*CHKERR();*/
         }
      }
    }
#line 4130 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 531: // $@9: %empty
#line 1695 "LefParser.yy"
               { driver.lefDumbMode = 10000000; driver.lefRealNum = 1; driver.lefInProp = 1; }
#line 4136 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 532: // layer_option: "PROPERTY" $@9 layer_prop_list ';'
#line 1696 "LefParser.yy"
    {
      driver.lefDumbMode = 0;
      driver.lefRealNum = 0;
      driver.lefInProp = 0;
    }
#line 4146 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 533: // $@10: %empty
#line 1702 "LefParser.yy"
    {
      if (driver.layerMastOver) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1527, "ACCURRENTDENSITY statement can't be defined in LAYER with TYPE MASTERSLICE or OVERLAP. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addAccurrentDensity((*(yystack_[0].value.stringVal)).c_str());
      delete (yystack_[0].value.stringVal);
    }
#line 4163 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 535: // layer_option: "ACCURRENTDENSITY" layer_table_type NUMBER ';'
#line 1716 "LefParser.yy"
    {
      if (driver.layerMastOver) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1527, "ACCURRENTDENSITY statement can't be defined in LAYER with TYPE MASTERSLICE or OVERLAP. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) {
           driver.lefrLayer.lefiLayer::addAccurrentDensity((*(yystack_[2].value.stringVal)).c_str());
           driver.lefrLayer.lefiLayer::setAcOneEntry((yystack_[1].value.doubleVal));
           delete (yystack_[2].value.stringVal);
      }
    }
#line 4183 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 536: // layer_option: "DCCURRENTDENSITY" "AVERAGE" NUMBER ';'
#line 1732 "LefParser.yy"
    {
      if (driver.layerMastOver) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1528, "DCCURRENTDENSITY statement can't be defined in LAYER with TYPE MASTERSLICE or OVERLAP. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) {
         driver.lefrLayer.lefiLayer::addDccurrentDensity("AVERAGE");
         driver.lefrLayer.lefiLayer::setDcOneEntry((yystack_[1].value.doubleVal));
      }
    }
#line 4202 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 537: // $@11: %empty
#line 1747 "LefParser.yy"
    {
      if (driver.layerMastOver) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1528, "DCCURRENTDENSITY statement can't be defined in LAYER with TYPE MASTERSLICE or OVERLAP. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (!driver.layerCut) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1529, "CUTAREA statement can only be defined in LAYER with type CUT. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) {
         driver.lefrLayer.lefiLayer::addDccurrentDensity("AVERAGE");
         driver.lefrLayer.lefiLayer::addNumber((yystack_[0].value.doubleVal));
      }
    }
#line 4229 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 538: // $@12: %empty
#line 1770 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addDcCutarea(); }
#line 4235 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 540: // $@13: %empty
#line 1773 "LefParser.yy"
    {
      if (driver.layerMastOver) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1528, "DCCURRENTDENSITY can't be defined in LAYER with TYPE MASTERSLICE or OVERLAP. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1530, "WIDTH statement can only be defined in LAYER with type ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) {
         driver.lefrLayer.lefiLayer::addDccurrentDensity("AVERAGE");
         driver.lefrLayer.lefiLayer::addNumber((yystack_[0].value.doubleVal));
      }
    }
#line 4262 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 541: // $@14: %empty
#line 1796 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addDcWidth(); }
#line 4268 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 543: // layer_option: "ANTENNAAREARATIO" NUMBER ';'
#line 1801 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "ANTENNAAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1531, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else if (driver.use5_3) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "ANTENNAAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNAAREARATIO syntax, which is illegal.", driver.versionNum);
               driver.lefError(1532, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (!driver.layerRout && !driver.layerCut && driver.layerMastOver) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1533, "ANTENNAAREARATIO statement can only be defined in LAYER with TYPE ROUTING or CUT. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setAntennaAreaRatio((yystack_[1].value.doubleVal));
    }
#line 4310 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 544: // $@15: %empty
#line 1839 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "ANTENNADIFFAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1532, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else if (driver.use5_3) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "ANTENNADIFFAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNADIFFAREARATIO syntax, which is illegal.", driver.versionNum);
               driver.lefError(1533, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (!driver.layerRout && !driver.layerCut && driver.layerMastOver) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.lefError(1534, "ANTENNADIFFAREARATIO statement can only be defined in LAYER with TYPE ROUTING or CUT. Parser stops executions.");
               /*CHKERR();*/
            }
         }
      }
      driver.antennaType = lefiAntennaDAR; 
    }
#line 4352 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 546: // layer_option: "ANTENNACUMAREARATIO" NUMBER ';'
#line 1878 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "ANTENNACUMAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1535, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else if (driver.use5_3) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "ANTENNACUMAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNACUMAREARATIO syntax, which is illegal.", driver.versionNum);
               driver.lefError(1536, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (!driver.layerRout && !driver.layerCut && driver.layerMastOver) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.lefError(1537, "ANTENNACUMAREARATIO statement can only be defined in LAYER with TYPE ROUTING or CUT. Parser stops executions.");
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setAntennaCumAreaRatio((yystack_[1].value.doubleVal));
    }
#line 4394 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 547: // $@16: %empty
#line 1916 "LefParser.yy"
    {  /* 5.4 syntax */
      driver.use5_4 = 1;
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "ANTENNACUMDIFFAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1538, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else if (driver.use5_3) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "ANTENNACUMDIFFAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNACUMDIFFAREARATIO syntax, which is illegal.", driver.versionNum);
               driver.lefError(1539, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (!driver.layerRout && !driver.layerCut && driver.layerMastOver) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1540, "ANTENNACUMDIFFAREARATIO statement can only be defined in LAYER with TYPE ROUTING or CUT. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      driver.antennaType = lefiAntennaCDAR;
    }
#line 4436 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 549: // $@17: %empty
#line 1955 "LefParser.yy"
    { /* both 5.3  & 5.4 syntax */
      if (!driver.layerRout && !driver.layerCut && driver.layerMastOver) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.lefError(1541, "ANTENNAAREAFACTOR can only be defined in LAYER with TYPE ROUTING or CUT. Parser stops executions.");
               /*CHKERR();*/
            }
         }
      }
      /* this does not need to check, since syntax is in both 5.3 & 5.4 */
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setAntennaAreaFactor((yystack_[0].value.doubleVal));
      driver.antennaType = lefiAntennaAF;
    }
#line 4454 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 551: // layer_option: "ANTENNASIDEAREARATIO" NUMBER ';'
#line 1970 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (!driver.layerRout && (driver.layerCut || driver.layerMastOver)) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.lefError(1542, "ANTENNASIDEAREARATIO can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
               /*CHKERR();*/
            }
         }
      }
      if (driver.ignoreVersion) {
        /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "ANTENNASIDEAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1543, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else if (driver.use5_3) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "ANTENNASIDEAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNASIDEAREARATIO syntax, which is illegal.", driver.versionNum);
               driver.lefError(1544, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setAntennaSideAreaRatio((yystack_[1].value.doubleVal));
    }
#line 4496 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 552: // $@18: %empty
#line 2008 "LefParser.yy"
    {  /* 5.4 syntax */
      driver.use5_4 = 1;
      if (!driver.layerRout && (driver.layerCut || driver.layerMastOver)) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.lefError(1545, "ANTENNADIFFSIDEAREARATIO can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
               /*CHKERR();*/
            }
         }
      }
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "ANTENNADIFFSIDEAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1546, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else if (driver.use5_3) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "ANTENNADIFFSIDEAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNADIFFSIDEAREARATIO syntax, which is illegal.", driver.versionNum);
               driver.lefError(1547, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      driver.antennaType = lefiAntennaDSAR;
    }
#line 4538 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 554: // layer_option: "ANTENNACUMSIDEAREARATIO" NUMBER ';'
#line 2047 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (!driver.layerRout && (driver.layerCut || driver.layerMastOver)) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.lefError(1548, "ANTENNACUMSIDEAREARATIO can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
               /*CHKERR();*/
            }
         }
      }
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "ANTENNACUMSIDEAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1549, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else if (driver.use5_3) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "ANTENNACUMSIDEAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNACUMSIDEAREARATIO syntax, which is illegal.", driver.versionNum);
               driver.lefError(1550, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setAntennaCumSideAreaRatio((yystack_[1].value.doubleVal));
    }
#line 4580 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 555: // $@19: %empty
#line 2085 "LefParser.yy"
    {  /* 5.4 syntax */
      driver.use5_4 = 1;
      if (!driver.layerRout && (driver.layerCut || driver.layerMastOver)) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.lefError(1551, "ANTENNACUMDIFFSIDEAREARATIO can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
               /*CHKERR();*/
            }
         }
      }
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "ANTENNACUMDIFFSIDEAREARATIO statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1552, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else if (driver.use5_3) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "ANTENNACUMDIFFSIDEAREARATIO statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNACUMDIFFSIDEAREARATIO syntax, which is illegal.", driver.versionNum);
               driver.lefError(1553, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      driver.antennaType = lefiAntennaCDSAR;
    }
#line 4622 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 557: // $@20: %empty
#line 2124 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (!driver.layerRout && (driver.layerCut || driver.layerMastOver)) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.lefError(1554, "ANTENNASIDEAREAFACTOR can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
               /*CHKERR();*/
            }
         }
      }
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "ANTENNASIDEAREAFACTOR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1555, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else if (driver.use5_3) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "ANTENNASIDEAREAFACTOR statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNASIDEAREAFACTOR syntax, which is illegal.", driver.versionNum);
               driver.lefError(1556, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setAntennaSideAreaFactor((yystack_[0].value.doubleVal));
      driver.antennaType = lefiAntennaSAF;
    }
#line 4665 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 559: // $@21: %empty
#line 2164 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (!driver.layerRout && !driver.layerCut && driver.layerMastOver) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.lefError(1557, "ANTENNAMODEL can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
               /*CHKERR();*/
            }
         }
      }
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.5) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "ANTENNAMODEL statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1558, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else if (driver.use5_3) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "ANTENNAMODEL statement is a version 5.4 or earlier syntax.\nYour lef file with version %g, has both old and new ANTENNAMODEL syntax, which is illegal.", driver.versionNum);
               driver.lefError(1559, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      driver.antennaType = lefiAntennaO;
    }
#line 4707 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 561: // layer_option: "ANTENNACUMROUTINGPLUSCUT" ';'
#line 2203 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
         driver.outMsg = (char*)lefMalloc(10000);
         sprintf(driver.outMsg,
           "ANTENNACUMROUTINGPLUSCUT is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
         driver.lefError(1686, driver.outMsg);
         lefFree(driver.outMsg);
         /*CHKERR();*/
      } else {
         if (!driver.layerRout && !driver.layerCut) {
            if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
               if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                  driver.lefError(1560, "ANTENNACUMROUTINGPLUSCUT can only be defined in LAYER with type ROUTING or CUT. Parser stops executions.");
                  /*CHKERR();*/
               }
            }
         }
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setAntennaCumRoutingPlusCut();
      }
    }
#line 4732 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 562: // layer_option: "ANTENNAGATEPLUSDIFF" NUMBER ';'
#line 2224 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
         driver.outMsg = (char*)lefMalloc(10000);
         sprintf(driver.outMsg,
           "ANTENNAGATEPLUSDIFF is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
         driver.lefError(1687, driver.outMsg);
         lefFree(driver.outMsg);
         /*CHKERR();*/
      } else {
         if (!driver.layerRout && !driver.layerCut) {
            if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
               if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                  driver.lefError(1561, "ANTENNAGATEPLUSDIFF can only be defined in LAYER with type ROUTING or CUT. Parser stops executions.");
                  /*CHKERR();*/
               }
            }
         }
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setAntennaGatePlusDiff((yystack_[1].value.doubleVal));
      }
    }
#line 4757 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 563: // layer_option: "ANTENNAAREAMINUSDIFF" NUMBER ';'
#line 2245 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
         driver.outMsg = (char*)lefMalloc(10000);
         sprintf(driver.outMsg,
           "ANTENNAAREAMINUSDIFF is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
         driver.lefError(1688, driver.outMsg);
         lefFree(driver.outMsg);
         /*CHKERR();*/
      } else {
         if (!driver.layerRout && !driver.layerCut) {
            if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
               if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                  driver.lefError(1562, "ANTENNAAREAMINUSDIFF can only be defined in LAYER with type ROUTING or CUT. Parser stops executions.");
                  /*CHKERR();*/
               }
            }
         }
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setAntennaAreaMinusDiff((yystack_[1].value.doubleVal));
      }
    }
#line 4782 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 564: // $@22: %empty
#line 2266 "LefParser.yy"
    {
      if (!driver.layerRout && !driver.layerCut) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.lefError(1563, "ANTENNAAREADIFFREDUCEPWL can only be defined in LAYER with type ROUTING or CUT. Parser stops executions.");
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) { /* require min 2 points, set the 1st 2 */
         driver.lefrAntennaPWLPtr=(lefiAntennaPWL*)lefMalloc(sizeof(lefiAntennaPWL));
         driver.lefrAntennaPWLPtr->lefiAntennaPWL::Init();
         driver.lefrAntennaPWLPtr->lefiAntennaPWL::addAntennaPWL((yystack_[1].value.pt)->x, (yystack_[1].value.pt)->y);
         driver.lefrAntennaPWLPtr->lefiAntennaPWL::addAntennaPWL((yystack_[0].value.pt)->x, (yystack_[0].value.pt)->y);
      }
      delete (yystack_[1].value.pt); delete (yystack_[0].value.pt);
    }
#line 4804 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 565: // $@23: %empty
#line 2284 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
        driver.lefrLayer.lefiLayer::setAntennaPWL(lefiAntennaADR, driver.lefrAntennaPWLPtr);
    }
#line 4813 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 566: // layer_option: "ANTENNAAREADIFFREDUCEPWL" '(' pt pt $@22 layer_diffusion_ratios ')' ';' $@23
#line 2288 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
        driver.outMsg = (char*)lefMalloc(10000);
        sprintf(driver.outMsg,
          "ANTENNAAREADIFFREDUCEPWL is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
        driver.lefError(1689, driver.outMsg);
        lefFree(driver.outMsg);
        /*CHKERR();*/
      }
    }
#line 4828 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 567: // layer_option: "SLOTWIREWIDTH" NUMBER ';'
#line 2299 "LefParser.yy"
    { /* 5.4 syntax */
      if (driver.ignoreVersion) {
         /* do nothing */
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setSlotWireWidth((yystack_[1].value.doubleVal));
      } else if (driver.versionNum >= 5.7) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings)
               driver.lefWarning(2011, "SLOTWIREWIDTH statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "SLOTWIREWIDTH statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1564, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setSlotWireWidth((yystack_[1].value.doubleVal));
    }
#line 4856 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 568: // layer_option: "SLOTWIRELENGTH" NUMBER ';'
#line 2323 "LefParser.yy"
    { /* 5.4 syntax */
      if (driver.ignoreVersion) {
         /* do nothing */
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setSlotWireLength((yystack_[1].value.doubleVal));
      } else if (driver.versionNum >= 5.7) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings)
               driver.lefWarning(2012, "SLOTWIRELENGTH statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "SLOTWIRELENGTH statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1565, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setSlotWireLength((yystack_[1].value.doubleVal));
    }
#line 4884 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 569: // layer_option: "SLOTWIDTH" NUMBER ';'
#line 2347 "LefParser.yy"
    { /* 5.4 syntax */
      if (driver.ignoreVersion) {
         /* do nothing */
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setSlotWidth((yystack_[1].value.doubleVal));
      } else if (driver.versionNum >= 5.7) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings)
               driver.lefWarning(2013, "SLOTWIDTH statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "SLOTWIDTH statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1566, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setSlotWidth((yystack_[1].value.doubleVal));
    }
#line 4912 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 570: // layer_option: "SLOTLENGTH" NUMBER ';'
#line 2371 "LefParser.yy"
    { /* 5.4 syntax */
      if (driver.ignoreVersion) {
         /* do nothing */
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setSlotLength((yystack_[1].value.doubleVal));
      } else if (driver.versionNum >= 5.7) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings)
               driver.lefWarning(2014, "SLOTLENGTH statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "SLOTLENGTH statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1567, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setSlotLength((yystack_[1].value.doubleVal));
    }
#line 4940 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 571: // layer_option: "MAXADJACENTSLOTSPACING" NUMBER ';'
#line 2395 "LefParser.yy"
    { /* 5.4 syntax */
      if (driver.ignoreVersion) {
         /* do nothing */
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setMaxAdjacentSlotSpacing((yystack_[1].value.doubleVal));
      } else if (driver.versionNum >= 5.7) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings)
               driver.lefWarning(2015, "MAXADJACENTSLOTSPACING statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "MAXADJACENTSLOTSPACING statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1568, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setMaxAdjacentSlotSpacing((yystack_[1].value.doubleVal));
    }
#line 4968 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 572: // layer_option: "MAXCOAXIALSLOTSPACING" NUMBER ';'
#line 2419 "LefParser.yy"
    { /* 5.4 syntax */
      if (driver.ignoreVersion) {
         /* do nothing */
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setMaxCoaxialSlotSpacing((yystack_[1].value.doubleVal));
      } else if (driver.versionNum >= 5.7) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings)
                driver.lefWarning(2016, "MAXCOAXIALSLOTSPACING statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "MAXCOAXIALSLOTSPACING statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1569, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setMaxCoaxialSlotSpacing((yystack_[1].value.doubleVal));
    }
#line 4996 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 573: // layer_option: "MAXEDGESLOTSPACING" NUMBER ';'
#line 2443 "LefParser.yy"
    { /* 5.4 syntax */
      if (driver.ignoreVersion) {
         /* do nothing */
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setMaxEdgeSlotSpacing((yystack_[1].value.doubleVal));
      } else if (driver.versionNum >= 5.7) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings)
               driver.lefWarning(2017, "MAXEDGESLOTSPACING statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "MAXEDGESLOTSPACING statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1570, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else
         if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setMaxEdgeSlotSpacing((yystack_[1].value.doubleVal));
    }
#line 5024 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 574: // layer_option: "SPLITWIREWIDTH" NUMBER ';'
#line 2467 "LefParser.yy"
    { /* 5.4 syntax */
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum >= 5.7) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings)
               driver.lefWarning(2018, "SPLITWIREWIDTH statement is obsolete in version 5.7 or later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.7 or later.");
         }
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "SPLITWIREWIDTH statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1571, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setSplitWireWidth((yystack_[1].value.doubleVal));
    }
#line 5051 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 575: // layer_option: "MINIMUMDENSITY" NUMBER ';'
#line 2490 "LefParser.yy"
    { /* 5.4 syntax, pcr 394389 */
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "MINIMUMDENSITY statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1572, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setMinimumDensity((yystack_[1].value.doubleVal));
    }
#line 5073 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 576: // layer_option: "MAXIMUMDENSITY" NUMBER ';'
#line 2508 "LefParser.yy"
    { /* 5.4 syntax, pcr 394389 */
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "MAXIMUMDENSITY statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1573, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setMaximumDensity((yystack_[1].value.doubleVal));
    }
#line 5095 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 577: // layer_option: "DENSITYCHECKWINDOW" NUMBER NUMBER ';'
#line 2526 "LefParser.yy"
    { /* 5.4 syntax, pcr 394389 */
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "DENSITYCHECKWINDOW statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1574, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setDensityCheckWindow((yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
    }
#line 5117 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 578: // layer_option: "DENSITYCHECKSTEP" NUMBER ';'
#line 2544 "LefParser.yy"
    { /* 5.4 syntax, pcr 394389 */
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "DENSITYCHECKSTEP statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1575, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setDensityCheckStep((yystack_[1].value.doubleVal));
    }
#line 5139 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 579: // layer_option: "FILLACTIVESPACING" NUMBER ';'
#line 2562 "LefParser.yy"
    { /* 5.4 syntax, pcr 394389 */
      if (driver.ignoreVersion) {
         /* do nothing */
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "FILLACTIVESPACING statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1576, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setFillActiveSpacing((yystack_[1].value.doubleVal));
    }
#line 5161 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 580: // layer_option: "MAXWIDTH" NUMBER ';'
#line 2580 "LefParser.yy"
    {
      // 5.5 MAXWIDTH, is for routing layer only
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.lefError(1577, "MAXWIDTH statement can only be defined in LAYER with TYPE ROUTING.  Parser stops executions.");
               /*CHKERR();*/
            }
         }
      }
      if (driver.versionNum < 5.5) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "MAXWIDTH statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1578, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setMaxwidth((yystack_[1].value.doubleVal));
    }
#line 5190 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 581: // layer_option: "MINWIDTH" NUMBER ';'
#line 2605 "LefParser.yy"
    {
      // 5.5 MINWIDTH, is for routing layer only
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1579, "MINWIDTH statement can only be defined in LAYER with TYPE ROUTING.  Parser stops executions.");
              /*CHKERR();*/
            }
         }
      }
      if (driver.versionNum < 5.5) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "MINWIDTH statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1580, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setMinwidth((yystack_[1].value.doubleVal));
    }
#line 5219 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 582: // $@24: %empty
#line 2630 "LefParser.yy"
    {
      if (driver.versionNum < 5.5) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "MINENCLOSEDAREA statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1581, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addMinenclosedarea((yystack_[0].value.doubleVal));
    }
#line 5239 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 584: // $@25: %empty
#line 2647 "LefParser.yy"
    { /* pcr 409334 */
      if (/*driver.lefrLayerCbk*/ 1)
         driver.lefrLayer.lefiLayer::addMinimumcut((int)(yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal)); 
      driver.hasLayerMincut = 0;
    }
#line 5249 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 585: // layer_option: "MINIMUMCUT" NUMBER "WIDTH" NUMBER $@25 layer_minimumcut_within layer_minimumcut_from layer_minimumcut_length ';'
#line 2655 "LefParser.yy"
    {
      if (!driver.hasLayerMincut) {   /* FROMABOVE nor FROMBELOW is set */
         if (/*driver.lefrLayerCbk*/ 1)
             driver.lefrLayer.lefiLayer::addMinimumcutConnect((char*)"");
      }
    }
#line 5260 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 586: // $@26: %empty
#line 2662 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addMinstep((yystack_[0].value.doubleVal));
    }
#line 5268 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 587: // layer_option: "MINSTEP" NUMBER $@26 layer_minstep_options ';'
#line 2666 "LefParser.yy"
    {
    }
#line 5275 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 588: // layer_option: "PROTRUSIONWIDTH" NUMBER "LENGTH" NUMBER "WIDTH" NUMBER ';'
#line 2669 "LefParser.yy"
    {
      if (driver.versionNum < 5.5) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "PROTRUSION RULE statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1582, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setProtrusion((yystack_[5].value.doubleVal), (yystack_[3].value.doubleVal), (yystack_[1].value.doubleVal));
    }
#line 5295 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 589: // $@27: %empty
#line 2685 "LefParser.yy"
    {
      if (driver.versionNum < 5.5) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "SPACINGTABLE statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1583, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      // 5.5 either SPACING or SPACINGTABLE in a layer, not both
      if (driver.lefrHasSpacing && driver.layerRout && driver.versionNum < 5.7) {
         if (/*driver.lefrLayerCbk*/ 1)  /* write warning only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefWarning(2010, "It is illegal to have both SPACING rules & SPACINGTABLE rules within a ROUTING layer");
            }
      } 
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addSpacingTable();
      driver.lefrHasSpacingTbl = 1;
    }
#line 5323 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 591: // $@28: %empty
#line 2711 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "ENCLOSURE statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1584, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1)
         driver.lefrLayer.lefiLayer::addEnclosure((*(yystack_[2].value.stringVal)).c_str(), (yystack_[1].value.doubleVal), (yystack_[0].value.doubleVal));
         delete (yystack_[2].value.stringVal);
    }
#line 5345 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 593: // $@29: %empty
#line 2731 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "PREFERENCLOSURE statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1585, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1)
         driver.lefrLayer.lefiLayer::addPreferEnclosure((*(yystack_[2].value.stringVal)).c_str(), (yystack_[1].value.doubleVal), (yystack_[0].value.doubleVal));
         delete (yystack_[2].value.stringVal);
    }
#line 5367 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 595: // layer_option: "RESISTANCE" NUMBER ';'
#line 2750 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "RESISTANCE statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1586, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else {
         if (/*driver.lefrLayerCbk*/ 1)
            driver.lefrLayer.lefiLayer::setResPerCut((yystack_[1].value.doubleVal));
      }
    }
#line 5389 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 596: // layer_option: "DIAGMINEDGELENGTH" NUMBER ';'
#line 2768 "LefParser.yy"
    {
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1587, "DIAGMINEDGELENGTH can only be defined in LAYER with TYPE ROUTING. Parser stops executions.");
              /*CHKERR();*/
            }
         }
      } else if (driver.versionNum < 5.6) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "DIAGMINEDGELENGTH statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1588, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else {
         if (/*driver.lefrLayerCbk*/ 1)
            driver.lefrLayer.lefiLayer::setDiagMinEdgeLength((yystack_[1].value.doubleVal));
      }
    }
#line 5418 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 597: // $@30: %empty
#line 2793 "LefParser.yy"
    {
      // Use the polygon code to retrieve the points for MINSIZE
      driver.lefrGeometriesPtr = (lefiGeometries*)lefMalloc(sizeof(lefiGeometries));
      driver.lefrGeometriesPtr->lefiGeometries::Init();
      driver.lefrDoGeometries = 1;
    }
#line 5429 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 598: // layer_option: "MINSIZE" $@30 firstPt otherPts ';'
#line 2800 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) {
         driver.lefrGeometriesPtr->lefiGeometries::addPolygon();
         driver.lefrLayer.lefiLayer::setMinSize(driver.lefrGeometriesPtr);
      }
     driver.lefrDoGeometries = 0;
      driver.lefrGeometriesPtr->lefiGeometries::Destroy();
                                       // Don't need the object anymore
      lefFree(driver.lefrGeometriesPtr);
    }
#line 5444 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 600: // layer_arraySpacing_long: "LONGARRAY"
#line 2814 "LefParser.yy"
    {
        if (/*driver.lefrLayerCbk*/ 1)
           driver.lefrLayer.lefiLayer::setArraySpacingLongArray();
    }
#line 5453 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 602: // layer_arraySpacing_width: "WIDTH" NUMBER
#line 2822 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
         driver.lefrLayer.lefiLayer::setArraySpacingWidth((yystack_[0].value.doubleVal));
    }
#line 5462 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 605: // layer_arraySpacing_arraycut: "ARRAYCUTS" NUMBER "SPACING" NUMBER
#line 2833 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
         driver.lefrLayer.lefiLayer::addArraySpacingArray((int)(yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
         if (driver.arrayCutsVal > (int)(yystack_[2].value.doubleVal)) {
            /* Mulitiple ARRAYCUTS value needs to me in ascending order */
            if (!driver.arrayCutsWar) {
               if (driver.layerWarnings++ < driver.lefrLayerWarnings)
                  driver.lefWarning(2080, "The number of cut values in multiple ARRAYSPACING ARRAYCUTS are not in increasing order.\nTo be consistent with the documentation, update the cut values to increasing order.");
               driver.arrayCutsWar = 1;
            }
         }
         driver.arrayCutsVal = (int)(yystack_[2].value.doubleVal);
    }
#line 5480 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 606: // $@31: %empty
#line 2849 "LefParser.yy"
    { 
      if (driver.hasInfluence) {  // 5.5 - INFLUENCE table must follow a PARALLEL
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1589, "An INFLUENCE table statement was defined before the PARALLELRUNLENGTH table statement.\nINFLUENCE table statement should be defined following the PARALLELRUNLENGTH.\nChange the LEF file and rerun the parser.");
              /*CHKERR();*/
            }
         }
      }
      if (driver.hasParallel) { // 5.5 - Only one PARALLEL table is allowed per layer
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1590, "There is multiple PARALLELRUNLENGTH table statements are defined within a layer.\nAccording to the LEF Reference Manual, only one PARALLELRUNLENGTH table statement is allowed per layer.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addNumber((yystack_[0].value.doubleVal));
      driver.hasParallel = 1;
    }
#line 5505 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 607: // $@32: %empty
#line 2870 "LefParser.yy"
    {
      driver.spParallelLength = driver.lefrLayer.lefiLayer::getNumber();
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addSpParallelLength();
    }
#line 5514 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 608: // $@33: %empty
#line 2875 "LefParser.yy"
    { 
      if (/*driver.lefrLayerCbk*/ 1) {
         driver.lefrLayer.lefiLayer::addSpParallelWidth((yystack_[0].value.doubleVal));
      }
    }
#line 5524 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 609: // $@34: %empty
#line 2881 "LefParser.yy"
    { 
      if (driver.lefrLayer.lefiLayer::getNumber() != driver.spParallelLength) {
         if (/*driver.lefrLayerCbk*/ 1) {
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1591, "The number of length in PARALLELRUNLENGTH is not the same as the number of spacing in WIDTH.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addSpParallelWidthSpacing();
    }
#line 5540 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 611: // $@35: %empty
#line 2895 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addNumber((yystack_[0].value.doubleVal));
    }
#line 5548 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 612: // $@36: %empty
#line 2899 "LefParser.yy"
    {
      if (driver.hasParallel) { // 5.7 - Either PARALLEL OR TWOWIDTHS per layer
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1592, "A PARALLELRUNLENGTH statement has already defined in the layer.\nOnly either PARALLELRUNLENGTH or TWOWIDTHS is allowed per layer.");
              /*CHKERR();*/
            }
         }
      }
      if (driver.hasTwoWidths) { // 5.7 - only 1 TWOWIDTHS per layer
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1593, "A TWOWIDTHS table statement has already defined in the layer.\nOnly one TWOWIDTHS statement is allowed per layer.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addSpTwoWidths((yystack_[4].value.doubleVal), (yystack_[3].value.doubleVal));
      driver.hasTwoWidths = 1;
    }
#line 5573 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 613: // sp_options: "TWOWIDTHS" "WIDTH" NUMBER layer_sp_TwoWidthsPRL NUMBER $@35 number_list $@36 layer_sp_TwoWidths
#line 2920 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
        driver.outMsg = (char*)lefMalloc(10000);
        sprintf(driver.outMsg,
          "TWOWIDTHS is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
        driver.lefError(1697, driver.outMsg);
        lefFree(driver.outMsg);
        /*CHKERR();*/
      } 
    }
#line 5588 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 614: // $@37: %empty
#line 2931 "LefParser.yy"
    {
      if (driver.hasInfluence) {  // 5.5 - INFLUENCE table must follow a PARALLEL
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1594, "A INFLUENCE table statement has already defined in the layer.\nOnly one INFLUENCE statement is allowed per layer.");
              /*CHKERR();*/
            }
         }
      }
      if (!driver.hasParallel) {  // 5.5 - INFLUENCE must follow a PARALLEL
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1595, "An INFLUENCE table statement has already defined beofre the layer.\nINFLUENCE statement has to be defined after the PARALLELRUNLENGTH table statement in the layer.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) {
         driver.lefrLayer.lefiLayer::setInfluence();
         driver.lefrLayer.lefiLayer::addSpInfluence((yystack_[4].value.doubleVal), (yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
      }
    }
#line 5615 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 618: // layer_spacingtable_opt: "WITHIN" NUMBER "SPACING" NUMBER
#line 2961 "LefParser.yy"
  {
    if (/*driver.lefrLayerCbk*/ 1)
       driver.lefrLayer.lefiLayer::addSpacingTableOrthoWithin((yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
  }
#line 5624 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 619: // layer_enclosure_type_opt: %empty
#line 2967 "LefParser.yy"
    {(yylhs.value.stringVal) = new std::string ("NULL");}
#line 5630 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 620: // layer_enclosure_type_opt: "ABOVE"
#line 2968 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("ABOVE");}
#line 5636 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 621: // layer_enclosure_type_opt: "BELOW"
#line 2969 "LefParser.yy"
             {(yylhs.value.stringVal) = new std::string ("BELOW");}
#line 5642 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 623: // $@38: %empty
#line 2973 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) {
         driver.lefrLayer.lefiLayer::addEnclosureWidth((yystack_[0].value.doubleVal));
      }
    }
#line 5652 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 625: // layer_enclosure_width_opt: "LENGTH" NUMBER
#line 2980 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
         driver.outMsg = (char*)lefMalloc(10000);
         sprintf(driver.outMsg,
           "LENGTH is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
         driver.lefError(1691, driver.outMsg);
         lefFree(driver.outMsg);
         /*CHKERR();*/
      } else {
         if (/*driver.lefrLayerCbk*/ 1) {
            driver.lefrLayer.lefiLayer::addEnclosureLength((yystack_[0].value.doubleVal));
         }
      }
    }
#line 5671 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 627: // layer_enclosure_width_except_opt: "EXCEPTEXTRACUT" NUMBER
#line 2997 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
         driver.outMsg = (char*)lefMalloc(10000);
         sprintf(driver.outMsg,
           "EXCEPTEXTRACUT is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
         driver.lefError(1690, driver.outMsg);
         lefFree(driver.outMsg);
         /*CHKERR();*/
      } else {
         if (/*driver.lefrLayerCbk*/ 1) {
            driver.lefrLayer.lefiLayer::addEnclosureExceptEC((yystack_[0].value.doubleVal));
         }
      }
    }
#line 5690 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 629: // layer_preferenclosure_width_opt: "WIDTH" NUMBER
#line 3014 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) {
         driver.lefrLayer.lefiLayer::addPreferEnclosureWidth((yystack_[0].value.doubleVal));
      }
    }
#line 5700 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 631: // layer_minimumcut_within: "WITHIN" NUMBER
#line 3022 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
        driver.outMsg = (char*)lefMalloc(10000);
        sprintf(driver.outMsg,
          "MINIMUMCUT WITHIN is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
        driver.lefError(1700, driver.outMsg);
        lefFree(driver.outMsg);
        /*CHKERR();*/
      } else {
         if (/*driver.lefrLayerCbk*/ 1) {
            driver.lefrLayer.lefiLayer::addMinimumcutWithin((yystack_[0].value.doubleVal));
         }
      }
    }
#line 5719 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 633: // layer_minimumcut_from: "FROMABOVE"
#line 3039 "LefParser.yy"
    {
      if (driver.versionNum < 5.5) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "FROMABOVE statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1596, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      driver.hasLayerMincut = 1;
      if (/*driver.lefrLayerCbk*/ 1)
         driver.lefrLayer.lefiLayer::addMinimumcutConnect((char*)"FROMABOVE");

    }
#line 5742 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 634: // layer_minimumcut_from: "FROMBELOW"
#line 3058 "LefParser.yy"
    {
      if (driver.versionNum < 5.5) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "FROMBELOW statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1597, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      }
      driver.hasLayerMincut = 1;
      if (/*driver.lefrLayerCbk*/ 1)
         driver.lefrLayer.lefiLayer::addMinimumcutConnect((char*)"FROMBELOW");
    }
#line 5764 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 636: // layer_minimumcut_length: "LENGTH" NUMBER "WITHIN" NUMBER
#line 3078 "LefParser.yy"
    {   
      if (driver.versionNum < 5.5) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "LENGTH WITHIN statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1598, driver.outMsg);
               lefFree(driver.outMsg);
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1)
         driver.lefrLayer.lefiLayer::addMinimumcutLengDis((yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
    }
#line 5785 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 639: // layer_minstep_option: layer_minstep_type
#line 3103 "LefParser.yy"
  {
    if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addMinstepType((*(yystack_[0].value.stringVal)).c_str());
    delete (yystack_[0].value.stringVal);
  }
#line 5794 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 640: // layer_minstep_option: "LENGTHSUM" NUMBER
#line 3108 "LefParser.yy"
  {
    if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addMinstepLengthsum((yystack_[0].value.doubleVal));
  }
#line 5802 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 641: // layer_minstep_option: "MAXEDGES" NUMBER
#line 3112 "LefParser.yy"
  {
    if (driver.versionNum < 5.7) {
      driver.outMsg = (char*)lefMalloc(10000);
      sprintf(driver.outMsg,
        "MAXEDGES is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
      driver.lefError(1685, driver.outMsg);
      lefFree(driver.outMsg);
      /*CHKERR();*/
    } else
       if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addMinstepMaxedges((int)(yystack_[0].value.doubleVal));
  }
#line 5818 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 642: // layer_minstep_type: "INSIDECORNER"
#line 3125 "LefParser.yy"
                 {(yylhs.value.stringVal) = new std::string ("INSIDECORNER");}
#line 5824 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 643: // layer_minstep_type: "OUTSIDECORNER"
#line 3126 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("OUTSIDECORNER");}
#line 5830 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 644: // layer_minstep_type: "STEP"
#line 3127 "LefParser.yy"
           {(yylhs.value.stringVal) = new std::string ("STEP");}
#line 5836 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 645: // layer_antenna_pwl: NUMBER
#line 3131 "LefParser.yy"
      { if (/*driver.lefrLayerCbk*/ 1)
          driver.lefrLayer.lefiLayer::setAntennaValue(driver.antennaType, (yystack_[0].value.doubleVal)); }
#line 5843 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 646: // $@39: %empty
#line 3134 "LefParser.yy"
      { if (/*driver.lefrLayerCbk*/ 1) { /* require min 2 points, set the 1st 2 */
          driver.lefrAntennaPWLPtr = (lefiAntennaPWL*)lefMalloc(sizeof(lefiAntennaPWL));
          driver.lefrAntennaPWLPtr->lefiAntennaPWL::Init();
          driver.lefrAntennaPWLPtr->lefiAntennaPWL::addAntennaPWL((yystack_[1].value.pt)->x, (yystack_[1].value.pt)->y);
          driver.lefrAntennaPWLPtr->lefiAntennaPWL::addAntennaPWL((yystack_[0].value.pt)->x, (yystack_[0].value.pt)->y);
        }
        delete (yystack_[1].value.pt); delete (yystack_[0].value.pt);
      }
#line 5856 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 647: // layer_antenna_pwl: "PWL" '(' pt pt $@39 layer_diffusion_ratios ')'
#line 3143 "LefParser.yy"
      { if (/*driver.lefrLayerCbk*/ 1)
          driver.lefrLayer.lefiLayer::setAntennaPWL(driver.antennaType, driver.lefrAntennaPWLPtr);
      }
#line 5864 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 650: // layer_diffusion_ratio: pt
#line 3156 "LefParser.yy"
  { if (/*driver.lefrLayerCbk*/ 1)
      driver.lefrAntennaPWLPtr->lefiAntennaPWL::addAntennaPWL((yystack_[0].value.pt)->x, (yystack_[0].value.pt)->y);
    delete (yystack_[0].value.pt);
  }
#line 5873 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 652: // layer_antenna_duo: "DIFFUSEONLY"
#line 3163 "LefParser.yy"
      { 
        driver.use5_4 = 1;
        if (driver.ignoreVersion) {
           /* do nothing */
        }
        else if ((driver.antennaType == lefiAntennaAF) && (driver.versionNum <= 5.3)) {
           if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
              if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                 driver.outMsg = (char*)lefMalloc(10000);
                 sprintf (driver.outMsg,
                   "ANTENNAAREAFACTOR with DIFFUSEONLY statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
                 driver.lefError(1599, driver.outMsg);
                 lefFree(driver.outMsg);
                 /*CHKERR();*/
              }
           }
        } else if (driver.use5_3) {
           if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
              if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                 driver.outMsg = (char*)lefMalloc(10000);
                 sprintf (driver.outMsg,
                   "ANTENNAAREAFACTOR with DIFFUSEONLY statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
                 driver.lefError(1599, driver.outMsg);
                 lefFree(driver.outMsg);
                 /*CHKERR();*/
              }
           }
        }
        if (/*driver.lefrLayerCbk*/ 1)
          driver.lefrLayer.lefiLayer::setAntennaDUO(driver.antennaType);
      }
#line 5909 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 653: // layer_table_type: "PEAK"
#line 3196 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("PEAK");}
#line 5915 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 654: // layer_table_type: "AVERAGE"
#line 3197 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("AVERAGE");}
#line 5921 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 655: // layer_table_type: "RMS"
#line 3198 "LefParser.yy"
               {(yylhs.value.stringVal) = new std::string ("RMS");}
#line 5927 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 656: // $@40: %empty
#line 3202 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addNumber((yystack_[0].value.doubleVal)); }
#line 5933 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 657: // $@41: %empty
#line 3204 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addAcFrequency(); }
#line 5939 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 658: // $@42: %empty
#line 3207 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addNumber((yystack_[0].value.doubleVal)); }
#line 5945 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 659: // layer_frequency: "FREQUENCY" NUMBER $@40 number_list ';' $@41 ac_layer_table_opt "TABLEENTRIES" NUMBER $@42 number_list ';'
#line 3209 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addAcTableEntry(); }
#line 5951 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 661: // $@43: %empty
#line 3213 "LefParser.yy"
    {
      if (!driver.layerCut) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1600, "CUTAREA statement can only be defined in LAYER with type CUT.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addNumber((yystack_[0].value.doubleVal));
    }
#line 5967 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 662: // ac_layer_table_opt: "CUTAREA" NUMBER $@43 number_list ';'
#line 3225 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addAcCutarea(); }
#line 5973 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 663: // $@44: %empty
#line 3227 "LefParser.yy"
    {
      if (!driver.layerRout) {
         if (/*driver.lefrLayerCbk*/ 1) { /* write error only if cbk is set */
            if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1601, "WIDTH can only be defined in LAYER with type ROUTING.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addNumber((yystack_[0].value.doubleVal));
    }
#line 5989 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 664: // ac_layer_table_opt: "WIDTH" NUMBER $@44 number_list ';'
#line 3239 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addAcWidth(); }
#line 5995 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 665: // $@45: %empty
#line 3243 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addNumber((yystack_[0].value.doubleVal)); }
#line 6001 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 666: // dc_layer_table: "TABLEENTRIES" NUMBER $@45 number_list ';'
#line 3245 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addDcTableEntry(); }
#line 6007 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 668: // number_list: number_list NUMBER
#line 3249 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addNumber((yystack_[0].value.doubleVal)); }
#line 6013 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 671: // layer_prop: "string" "string"
#line 3258 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) {
        char propTp;
        propTp = driver.lefrLayerProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
        driver.lefrLayer.lefiLayer::addProp((*(yystack_[1].value.stringVal)).c_str(), (*(yystack_[0].value.stringVal)).c_str(), propTp);
      }
      delete (yystack_[1].value.stringVal); delete (yystack_[0].value.stringVal);
    }
#line 6026 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 672: // layer_prop: "string" "qstring"
#line 3267 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) {
        char propTp;
        propTp = driver.lefrLayerProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
        driver.lefrLayer.lefiLayer::addProp((*(yystack_[1].value.stringVal)).c_str(), (*(yystack_[0].value.qstringVal)).c_str(), propTp);
      }
      delete (yystack_[1].value.stringVal); delete (yystack_[0].value.qstringVal);
    }
#line 6039 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 673: // layer_prop: "string" NUMBER
#line 3276 "LefParser.yy"
    {
      char temp[32];
      sprintf(temp, "%.11g", (yystack_[0].value.doubleVal));
      if (/*driver.lefrLayerCbk*/ 1) {
        char propTp;
        propTp = driver.lefrLayerProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
        driver.lefrLayer.lefiLayer::addNumProp((*(yystack_[1].value.stringVal)).c_str(), (yystack_[0].value.doubleVal), temp, propTp);
      }
      delete (yystack_[1].value.stringVal); 
    }
#line 6054 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 674: // current_density_pwl_list: current_density_pwl
#line 3289 "LefParser.yy"
    { }
#line 6060 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 675: // current_density_pwl_list: current_density_pwl_list current_density_pwl
#line 3291 "LefParser.yy"
    { }
#line 6066 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 676: // current_density_pwl: '(' NUMBER NUMBER ')'
#line 3294 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setCurrentPoint((yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal)); }
#line 6072 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 679: // cap_point: '(' NUMBER NUMBER ')'
#line 3302 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setCapacitancePoint((yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal)); }
#line 6078 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 681: // res_points: res_points res_point
#line 3307 "LefParser.yy"
    { }
#line 6084 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 682: // res_point: '(' NUMBER NUMBER ')'
#line 3310 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::setResistancePoint((yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal)); }
#line 6090 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 683: // layer_type: "ROUTING"
#line 3313 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("ROUTING"); driver.layerRout = 1;}
#line 6096 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 684: // layer_type: "CUT"
#line 3314 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("CUT"); driver.layerCut = 1;}
#line 6102 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 685: // layer_type: "OVERLAP"
#line 3315 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("OVERLAP"); driver.layerMastOver = 1;}
#line 6108 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 686: // layer_type: "MASTERSLICE"
#line 3316 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("MASTERSLICE"); driver.layerMastOver = 1;}
#line 6114 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 687: // layer_type: "VIRTUAL"
#line 3317 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("VIRTUAL");}
#line 6120 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 688: // layer_type: "IMPLANT"
#line 3318 "LefParser.yy"
                  {(yylhs.value.stringVal) = new std::string ("IMPLANT");}
#line 6126 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 689: // layer_direction: "HORIZONTAL"
#line 3321 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("HORIZONTAL");}
#line 6132 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 690: // layer_direction: "VERTICAL"
#line 3322 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("VERTICAL");}
#line 6138 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 691: // layer_direction: "DIAG45"
#line 3323 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("DIAG45");}
#line 6144 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 692: // layer_direction: "DIAG135"
#line 3324 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("DIAG135");}
#line 6150 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 694: // layer_minen_width: "WIDTH" NUMBER
#line 3328 "LefParser.yy"
    {
    if (/*driver.lefrLayerCbk*/ 1)
       driver.lefrLayer.lefiLayer::addMinenclosedareaWidth((yystack_[0].value.doubleVal));
    }
#line 6159 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 695: // layer_oxide: "OXIDE1"
#line 3335 "LefParser.yy"
    {
    if (/*driver.lefrLayerCbk*/ 1)
       driver.lefrLayer.lefiLayer::addAntennaModel(1);
    }
#line 6168 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 696: // layer_oxide: "OXIDE2"
#line 3340 "LefParser.yy"
    {
    if (/*driver.lefrLayerCbk*/ 1)
       driver.lefrLayer.lefiLayer::addAntennaModel(2);
    }
#line 6177 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 697: // layer_oxide: "OXIDE3"
#line 3345 "LefParser.yy"
    {
    if (/*driver.lefrLayerCbk*/ 1)
       driver.lefrLayer.lefiLayer::addAntennaModel(3);
    }
#line 6186 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 698: // layer_oxide: "OXIDE4"
#line 3350 "LefParser.yy"
    {
    if (/*driver.lefrLayerCbk*/ 1)
       driver.lefrLayer.lefiLayer::addAntennaModel(4);
    }
#line 6195 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 699: // layer_sp_parallel_widths: %empty
#line 3356 "LefParser.yy"
    { }
#line 6201 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 700: // layer_sp_parallel_widths: layer_sp_parallel_widths layer_sp_parallel_width
#line 3361 "LefParser.yy"
    { }
#line 6207 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 701: // $@46: %empty
#line 3364 "LefParser.yy"
    { 
      if (/*driver.lefrLayerCbk*/ 1) {
         driver.lefrLayer.lefiLayer::addSpParallelWidth((yystack_[0].value.doubleVal));
      }
    }
#line 6217 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 702: // layer_sp_parallel_width: "WIDTH" NUMBER $@46 number_list
#line 3370 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addSpParallelWidthSpacing(); }
#line 6223 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 703: // layer_sp_TwoWidths: %empty
#line 3373 "LefParser.yy"
    { }
#line 6229 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 704: // layer_sp_TwoWidths: layer_sp_TwoWidth layer_sp_TwoWidths
#line 3375 "LefParser.yy"
    { }
#line 6235 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 705: // $@47: %empty
#line 3389 "LefParser.yy"
    {
       if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addNumber((yystack_[0].value.doubleVal));
    }
#line 6243 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 706: // layer_sp_TwoWidth: "WIDTH" NUMBER layer_sp_TwoWidthsPRL NUMBER $@47 number_list
#line 3393 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
         driver.lefrLayer.lefiLayer::addSpTwoWidths((yystack_[4].value.doubleVal), (yystack_[3].value.doubleVal));
    }
#line 6252 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 707: // layer_sp_TwoWidthsPRL: %empty
#line 3399 "LefParser.yy"
    { (yylhs.value.doubleVal) = -1; }
#line 6258 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 708: // layer_sp_TwoWidthsPRL: "PRL" NUMBER
#line 3401 "LefParser.yy"
    { (yylhs.value.doubleVal) = (yystack_[0].value.doubleVal); }
#line 6264 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 709: // layer_sp_influence_widths: %empty
#line 3404 "LefParser.yy"
    { }
#line 6270 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 710: // layer_sp_influence_widths: layer_sp_influence_widths layer_sp_influence_width
#line 3406 "LefParser.yy"
    { }
#line 6276 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 711: // layer_sp_influence_width: "WIDTH" NUMBER "WITHIN" NUMBER "SPACING" NUMBER
#line 3409 "LefParser.yy"
    { if (/*driver.lefrLayerCbk*/ 1) driver.lefrLayer.lefiLayer::addSpInfluence((yystack_[4].value.doubleVal), (yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal)); }
#line 6282 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 712: // maxstack_via: "MAXVIASTACK" NUMBER ';'
#line 3412 "LefParser.yy"
    {
      if (!driver.lefrHasLayer) {  /* 5.5 */
        if (/*lefrMaxStackViaCbk*/ 1) { /* write error only if cbk is set */
           if (driver.maxStackViaWarnings++ < driver.lefrMaxStackViaWarnings) {
             driver.lefError(1602, "MAXVIASTACK statement has to be defined after the LAYER statement.");
             /*CHKERR();*/
           }
        }
      } else if (driver.lefrHasMaxVS) {
        if (/*lefrMaxStackViaCbk*/ 1) { /* write error only if cbk is set */
           if (driver.maxStackViaWarnings++ < driver.lefrMaxStackViaWarnings) {
             driver.lefError(1603, "A MAXVIASTACK has already defined.\nOnly one MAXVIASTACK is allowed per lef file.");
             /*CHKERR();*/
           }
        }
      } else {
        if (/*lefrMaxStackViaCbk*/ 1) {
           driver.lefrMaxStackVia.lefiMaxStackVia::setMaxStackVia((int)(yystack_[1].value.doubleVal));
           driver.lefrMaxStackViaCbk( driver.lefrMaxStackVia);
        }
      }
      if (driver.versionNum < 5.5) {
        if (/*lefrMaxStackViaCbk*/ 1) { /* write error only if cbk is set */
           if (driver.maxStackViaWarnings++ < driver.lefrMaxStackViaWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                "MAXVIASTACK statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1604, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
      driver.lefrHasMaxVS = 1;
    }
#line 6322 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 713: // $@48: %empty
#line 3447 "LefParser.yy"
                                 {driver.lefDumbMode = 2; driver.lefNoNum= 2;}
#line 6328 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 714: // maxstack_via: "MAXVIASTACK" NUMBER "RANGE" $@48 "string" "string" ';'
#line 3449 "LefParser.yy"
    {
      if (!driver.lefrHasLayer) {  /* 5.5 */
        if (/*lefrMaxStackViaCbk*/ 1) { /* write error only if cbk is set */
           if (driver.maxStackViaWarnings++ < driver.lefrMaxStackViaWarnings) {
              driver.lefError(1602, "MAXVIASTACK statement has to be defined after the LAYER statement.");
              /*CHKERR();*/
           }
        }
      } else if (driver.lefrHasMaxVS) {
        if (/*lefrMaxStackViaCbk*/ 1) { /* write error only if cbk is set */
           if (driver.maxStackViaWarnings++ < driver.lefrMaxStackViaWarnings) {
             driver.lefError(1603, "A MAXVIASTACK has already defined.\nOnly one MAXVIASTACK is allowed per lef file.");
             /*CHKERR();*/
           }
        }
      } else {
        if (/*lefrMaxStackViaCbk*/ 1) {
           driver.lefrMaxStackVia.lefiMaxStackVia::setMaxStackVia((int)(yystack_[5].value.doubleVal));
           driver.lefrMaxStackVia.lefiMaxStackVia::setMaxStackViaRange((*(yystack_[2].value.stringVal)).c_str(), (*(yystack_[1].value.stringVal)).c_str());
           driver.lefrMaxStackViaCbk( driver.lefrMaxStackVia);
        }
      }
      driver.lefrHasMaxVS = 1;
      delete (yystack_[2].value.stringVal); delete (yystack_[1].value.stringVal);
    }
#line 6358 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 715: // $@49: %empty
#line 3475 "LefParser.yy"
                { driver.hasViaRule_layer = 0; }
#line 6364 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 716: // via: start_via $@49 via_option end_via
#line 3476 "LefParser.yy"
    { 
      if (/*lefrViaCbk*/ 1)
        driver.lefrViaCbk( driver.lefrVia);
    }
#line 6373 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 717: // via_keyword: "VIA"
#line 3482 "LefParser.yy"
     { driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 6379 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 718: // start_via: via_keyword "string"
#line 3485 "LefParser.yy"
    {
      /* 0 is nodefault */
      if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setName((*(yystack_[0].value.stringVal)).c_str(), 0);
      driver.viaLayer = 0;
      driver.numVia++;
      //strcpy(driver.viaName, $2);
      driver.viaName = *(yystack_[0].value.stringVal);
      delete (yystack_[0].value.stringVal);
    }
#line 6393 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 719: // start_via: via_keyword "string" "DEFAULT"
#line 3495 "LefParser.yy"
    {
      /* 1 is default */
      if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setName((*(yystack_[1].value.stringVal)).c_str(), 1);
      driver.viaLayer = 0;
      //strcpy(driver.viaName, $2);
      driver.viaName = *(yystack_[1].value.stringVal);
      delete (yystack_[1].value.stringVal);
    }
#line 6406 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 720: // start_via: via_keyword "string" "GENERATED"
#line 3504 "LefParser.yy"
    {
      /* 2 is generated */
      if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setName((*(yystack_[1].value.stringVal)).c_str(), 2);
      driver.viaLayer = 0;
      //strcpy(driver.viaName, $2);
      driver.viaName = *(yystack_[1].value.stringVal);
      delete (yystack_[1].value.stringVal);
    }
#line 6419 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 721: // $@50: %empty
#line 3513 "LefParser.yy"
                       {driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 6425 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 722: // $@51: %empty
#line 3515 "LefParser.yy"
           {driver.lefDumbMode = 3; driver.lefNoNum = 1; }
#line 6431 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 723: // $@52: %empty
#line 3518 "LefParser.yy"
    {
       if (driver.versionNum < 5.6) {
         if (/*lefrViaCbk*/ 1) { /* write error only if cbk is set */
            if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                "VIARULE statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1605, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
            }
         }
       }  else
          if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setViaRule((*(yystack_[21].value.stringVal)).c_str(), (yystack_[18].value.doubleVal), (yystack_[17].value.doubleVal), (*(yystack_[13].value.stringVal)).c_str(), (*(yystack_[12].value.stringVal)).c_str(), (*(yystack_[11].value.stringVal)).c_str(),
                          (yystack_[8].value.doubleVal), (yystack_[7].value.doubleVal), (yystack_[4].value.doubleVal), (yystack_[3].value.doubleVal), (yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
       driver.viaLayer++;
       driver.hasViaRule_layer = 1;
       delete (yystack_[21].value.stringVal); delete (yystack_[13].value.stringVal); delete (yystack_[12].value.stringVal); delete (yystack_[11].value.stringVal);
    }
#line 6455 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 727: // via_viarule_option: "ROWCOL" NUMBER NUMBER ';'
#line 3545 "LefParser.yy"
    {
       if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setRowCol((int)(yystack_[2].value.doubleVal), (int)(yystack_[1].value.doubleVal));
    }
#line 6463 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 728: // via_viarule_option: "ORIGIN" NUMBER NUMBER ';'
#line 3549 "LefParser.yy"
    {
       if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setOrigin((yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
    }
#line 6471 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 729: // via_viarule_option: "OFFSET" NUMBER NUMBER NUMBER NUMBER ';'
#line 3553 "LefParser.yy"
    {
       if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setOffset((yystack_[4].value.doubleVal), (yystack_[3].value.doubleVal), (yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
    }
#line 6479 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 730: // $@53: %empty
#line 3556 "LefParser.yy"
              {driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 6485 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 731: // via_viarule_option: "PATTERN" $@53 "string" ';'
#line 3557 "LefParser.yy"
    {
       if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setPattern((*(yystack_[1].value.stringVal)).c_str());
       delete (yystack_[1].value.stringVal);
    }
#line 6494 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 737: // via_other_option: via_foreign
#line 3575 "LefParser.yy"
    { }
#line 6500 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 738: // via_other_option: via_layer_rule
#line 3577 "LefParser.yy"
    { }
#line 6506 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 739: // via_other_option: "RESISTANCE" NUMBER ';'
#line 3579 "LefParser.yy"
    { if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setResistance((yystack_[1].value.doubleVal)); }
#line 6512 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 740: // $@54: %empty
#line 3580 "LefParser.yy"
               { driver.lefDumbMode = 1000000; driver.lefRealNum = 1; driver.lefInProp = 1; }
#line 6518 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 741: // via_other_option: "PROPERTY" $@54 via_prop_list ';'
#line 3581 "LefParser.yy"
    { driver.lefDumbMode = 0;
      driver.lefRealNum = 0;
      driver.lefInProp = 0;
    }
#line 6527 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 742: // via_other_option: "TOPOFSTACKONLY"
#line 3586 "LefParser.yy"
    { 
      if (driver.versionNum < 5.6) {
        if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setTopOfStack();
      } else
        if (/*lefrViaCbk*/ 1)  /* write warning only if cbk is set */
           if (driver.viaWarnings++ < driver.lefrViaWarnings)
              driver.lefWarning(2019, "TOPOFSTACKONLY statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later");
    }
#line 6540 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 745: // via_name_value_pair: "string" NUMBER
#line 3602 "LefParser.yy"
    { 
      char temp[32];
      sprintf(temp, "%.11g", (yystack_[0].value.doubleVal));
      if (/*lefrViaCbk*/ 1) {
         char propTp;
         propTp = driver.lefrViaProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrVia.lefiVia::addNumProp((*(yystack_[1].value.stringVal)).c_str(), (yystack_[0].value.doubleVal), temp, propTp);
      }
     delete (yystack_[1].value.stringVal);
    }
#line 6555 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 746: // via_name_value_pair: "string" "qstring"
#line 3613 "LefParser.yy"
    {
      if (/*lefrViaCbk*/ 1) {
         char propTp;
         propTp = driver.lefrViaProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrVia.lefiVia::addProp((*(yystack_[1].value.stringVal)).c_str(), (*(yystack_[0].value.qstringVal)).c_str(), propTp);
      }
     delete (yystack_[1].value.stringVal); delete (yystack_[0].value.qstringVal);
    }
#line 6568 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 747: // via_name_value_pair: "string" "string"
#line 3622 "LefParser.yy"
    {
      if (/*lefrViaCbk*/ 1) {
         char propTp;
         propTp = driver.lefrViaProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrVia.lefiVia::addProp((*(yystack_[1].value.stringVal)).c_str(), (*(yystack_[0].value.stringVal)).c_str(), propTp);
      }
     delete (yystack_[1].value.stringVal); delete (yystack_[0].value.stringVal);
    }
#line 6581 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 748: // via_foreign: start_foreign ';'
#line 3633 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setForeign((*(yystack_[1].value.stringVal)).c_str(), 0, 0.0, 0.0, -1);
      } else
        if (/*lefrViaCbk*/ 1)  /* write warning only if cbk is set */
           if (driver.viaWarnings++ < driver.lefrViaWarnings)
             driver.lefWarning(2020, "FOREIGN statement in VIA is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
        delete (yystack_[1].value.stringVal);
    }
#line 6595 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 749: // via_foreign: start_foreign pt ';'
#line 3643 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setForeign((*(yystack_[2].value.stringVal)).c_str(), 1, (yystack_[1].value.pt)->x, (yystack_[1].value.pt)->y, -1);
      } else
        if (/*lefrViaCbk*/ 1)  /* write warning only if cbk is set */
           if (driver.viaWarnings++ < driver.lefrViaWarnings)
             driver.lefWarning(2020, "FOREIGN statement in VIA is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
     delete (yystack_[2].value.stringVal); delete (yystack_[1].value.pt);
    }
#line 6609 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 750: // via_foreign: start_foreign pt orientation ';'
#line 3653 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setForeign((*(yystack_[3].value.stringVal)).c_str(), 1, (yystack_[2].value.pt)->x, (yystack_[2].value.pt)->y, (yystack_[1].value.integerVal));
      } else
        if (/*lefrViaCbk*/ 1)  /* write warning only if cbk is set */
           if (driver.viaWarnings++ < driver.lefrViaWarnings)
             driver.lefWarning(2020, "FOREIGN statement in VIA is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
     delete (yystack_[3].value.stringVal); delete (yystack_[2].value.pt);
    }
#line 6623 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 751: // via_foreign: start_foreign orientation ';'
#line 3663 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::setForeign((*(yystack_[2].value.stringVal)).c_str(), 0, 0.0, 0.0, (yystack_[1].value.integerVal));
      } else
        if (/*lefrViaCbk*/ 1)  /* write warning only if cbk is set */
           if (driver.viaWarnings++ < driver.lefrViaWarnings)
             driver.lefWarning(2020, "FOREIGN statement in VIA is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
     delete (yystack_[2].value.stringVal);
    }
#line 6637 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 752: // $@55: %empty
#line 3673 "LefParser.yy"
                          {driver.lefDumbMode = 1; driver.lefNoNum= 1;}
#line 6643 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 753: // start_foreign: "FOREIGN" $@55 "string"
#line 3674 "LefParser.yy"
    { (yylhs.value.stringVal) = (yystack_[0].value.stringVal); }
#line 6649 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 754: // orientation: "N"
#line 3677 "LefParser.yy"
              {(yylhs.value.integerVal) = 0;}
#line 6655 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 755: // orientation: "W"
#line 3678 "LefParser.yy"
              {(yylhs.value.integerVal) = 1;}
#line 6661 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 756: // orientation: "S"
#line 3679 "LefParser.yy"
              {(yylhs.value.integerVal) = 2;}
#line 6667 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 757: // orientation: "E"
#line 3680 "LefParser.yy"
              {(yylhs.value.integerVal) = 3;}
#line 6673 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 758: // orientation: "FN"
#line 3681 "LefParser.yy"
              {(yylhs.value.integerVal) = 4;}
#line 6679 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 759: // orientation: "FW"
#line 3682 "LefParser.yy"
              {(yylhs.value.integerVal) = 5;}
#line 6685 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 760: // orientation: "FS"
#line 3683 "LefParser.yy"
              {(yylhs.value.integerVal) = 6;}
#line 6691 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 761: // orientation: "FE"
#line 3684 "LefParser.yy"
              {(yylhs.value.integerVal) = 7;}
#line 6697 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 762: // orientation: "R0"
#line 3685 "LefParser.yy"
              {(yylhs.value.integerVal) = 0;}
#line 6703 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 763: // orientation: "R90"
#line 3686 "LefParser.yy"
              {(yylhs.value.integerVal) = 1;}
#line 6709 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 764: // orientation: "R180"
#line 3687 "LefParser.yy"
              {(yylhs.value.integerVal) = 2;}
#line 6715 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 765: // orientation: "R270"
#line 3688 "LefParser.yy"
              {(yylhs.value.integerVal) = 3;}
#line 6721 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 766: // orientation: "MY"
#line 3689 "LefParser.yy"
              {(yylhs.value.integerVal) = 4;}
#line 6727 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 767: // orientation: "MYR90"
#line 3690 "LefParser.yy"
              {(yylhs.value.integerVal) = 5;}
#line 6733 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 768: // orientation: "MX"
#line 3691 "LefParser.yy"
              {(yylhs.value.integerVal) = 6;}
#line 6739 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 769: // orientation: "MXR90"
#line 3692 "LefParser.yy"
              {(yylhs.value.integerVal) = 7;}
#line 6745 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 770: // via_layer_rule: via_layer via_geometries
#line 3695 "LefParser.yy"
    { }
#line 6751 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 771: // $@56: %empty
#line 3697 "LefParser.yy"
                   {driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 6757 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 772: // via_layer: "LAYER" $@56 "string" ';'
#line 3698 "LefParser.yy"
    {
      if (/*lefrViaCbk*/ 1) driver.lefrVia.lefiVia::addLayer((*(yystack_[1].value.stringVal)).c_str());
      driver.viaLayer++;
      driver.hasViaRule_layer = 1;
      delete (yystack_[1].value.stringVal);
    }
#line 6768 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 775: // via_geometry: "RECT" pt pt ';'
#line 3712 "LefParser.yy"
    { if (/*lefrViaCbk*/ 1)
        driver.lefrVia.lefiVia::addRectToLayer((yystack_[2].value.pt)->x, (yystack_[2].value.pt)->y, (yystack_[1].value.pt)->x, (yystack_[1].value.pt)->y); 
    delete (yystack_[2].value.pt); delete (yystack_[1].value.pt);
    }
#line 6777 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 776: // $@57: %empty
#line 3717 "LefParser.yy"
    {
      driver.lefrGeometriesPtr = (lefiGeometries*)lefMalloc(sizeof(lefiGeometries));
      driver.lefrGeometriesPtr->lefiGeometries::Init();
      driver.lefrDoGeometries = 1;
    }
#line 6787 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 777: // via_geometry: "POLYGON" $@57 firstPt nextPt nextPt otherPts ';'
#line 3723 "LefParser.yy"
    { 
      if (/*lefrViaCbk*/ 1) {
        driver.lefrGeometriesPtr->lefiGeometries::addPolygon();
        driver.lefrVia.lefiVia::addPolyToLayer(driver.lefrGeometriesPtr);   // 5.6
      }
      driver.lefrGeometriesPtr->lefiGeometries::clearPolyItems(); // free items fields
      lefFree((char*)(driver.lefrGeometriesPtr)); // Don't need anymore, poly data has
      driver.lefrDoGeometries = 0;                // copied
    }
#line 6801 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 778: // $@58: %empty
#line 3733 "LefParser.yy"
               {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 6807 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 779: // end_via: "END" $@58 "string"
#line 3734 "LefParser.yy"
    { 
      // 10/17/2001 - Wanda da Rosa, PCR 404149
      //              Error if no layer in via
      if (!driver.viaLayer) {
         if (/*lefrViaCbk*/ 1) {  /* write error only if cbk is set */
            if (driver.viaWarnings++ < driver.lefrViaWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                "A LAYER statement is missing in the VIA %s.\nAt least one LAYERis required per VIA statement.", (*(yystack_[0].value.stringVal)).c_str());
              driver.lefError(1606, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
            }
         }
      }
      if (driver.viaName != *(yystack_[0].value.stringVal)) {
         if (/*lefrViaCbk*/ 1) { /* write error only if cbk is set */
            if (driver.viaWarnings++ < driver.lefrViaWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                "END VIA name %s is different from the VIA name %s.\nCorrect the LEF file before rerun it through the LEF parser.", (*(yystack_[0].value.stringVal)).c_str(), driver.viaName.c_str());
              driver.lefError(1607, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
            } 
         } 
      } 
      delete (yystack_[0].value.stringVal);
    }
#line 6841 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 780: // $@59: %empty
#line 3764 "LefParser.yy"
                            { driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 6847 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 781: // viarule_keyword: "VIARULE" $@59 "string"
#line 3765 "LefParser.yy"
    { 
      if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::setName((*(yystack_[0].value.stringVal)).c_str());
      driver.viaRuleLayer = 0;
      //strcpy(driver.viaRuleName, $3);
      driver.viaRuleName = (*(yystack_[0].value.stringVal));
      driver.isGenerate = 0;
      delete (yystack_[0].value.stringVal);
    }
#line 6860 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 782: // viarule: viarule_keyword viarule_layer_list via_names opt_viarule_props end_viarule
#line 3776 "LefParser.yy"
    {
      // 10/17/2001 - Wanda da Rosa, PCR 404163
      //              Error if layer number is not equal 2.
      // 11/14/2001 - Wanda da Rosa,
      //              Commented out for pcr 411781
      //if (driver.viaRuleLayer != 2) {
         //driver.lefError("VIARULE requires two layers");
         ///*CHKERR();*/
      //}
      if (driver.viaRuleLayer == 0 || driver.viaRuleLayer > 2) {
         if (/*driver.lefrViaRuleCbk*/ 1) {  /* write error only if cbk is set */
            if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
              driver.lefError(1608, "A VIARULE statement requires two layers.");
              /*CHKERR();*/
            }
         }
      }
      if (/*driver.lefrViaRuleCbk*/ 1)
        driver.lefrViaRuleCbk( driver.lefrViaRule);
      // 2/19/2004 - reset the ENCLOSURE overhang values which may be
      // set by the old syntax OVERHANG -- Not necessary, but just incase
      if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::clearLayerOverhang();
    }
#line 6888 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 783: // $@60: %empty
#line 3802 "LefParser.yy"
    {
      driver.isGenerate = 1;
    }
#line 6896 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 784: // viarule_generate: viarule_keyword "GENERATE" viarule_generate_default $@60 viarule_layer_list opt_viarule_props end_viarule
#line 3806 "LefParser.yy"
    {
      // 10/17/2001 - Wanda da Rosa, PCR 404181
      //              Error if layer number is not equal 3.
      // 11/14/2001 - Wanda da Rosa,
      //              Commented out for pcr 411781
      //if (driver.viaRuleLayer != 3) {
         //driver.lefError("VIARULE requires three layers");
         ///*CHKERR();*/
      //}
      if (driver.viaRuleLayer == 0) {
         if (/*driver.lefrViaRuleCbk*/ 1) {  /* write error only if cbk is set */
            if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
              driver.lefError(1609, "A VIARULE GENERATE requires three layers.");
              /*CHKERR();*/
            }
         }
      } else if ((driver.viaRuleLayer < 3) && (driver.versionNum >= 5.6)) {
         if (/*driver.lefrViaRuleCbk*/ 1)  /* write warning only if cbk is set */
            if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings)
              driver.lefWarning(2021, "turn-via is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
      } else {
         if (/*driver.lefrViaRuleCbk*/ 1) {
            driver.lefrViaRule.lefiViaRule::setGenerate();
            driver.lefrViaRuleCbk( driver.lefrViaRule);
         }
      }
      // 2/19/2004 - reset the ENCLOSURE overhang values which may be
      // set by the old syntax OVERHANG
      if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::clearLayerOverhang();
    }
#line 6931 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 786: // viarule_generate_default: "DEFAULT"
#line 3839 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
         if (/*driver.lefrViaRuleCbk*/ 1) {  /* write error only if cbk is set */
            if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                "DEFAULT statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1605, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
            }
         }
      } else
        if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::setDefault();
    }
#line 6951 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 793: // $@61: %empty
#line 3870 "LefParser.yy"
                         { driver.lefDumbMode = 10000000; driver.lefRealNum = 1; driver.lefInProp = 1; }
#line 6957 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 794: // viarule_prop: "PROPERTY" $@61 viarule_prop_list ';'
#line 3871 "LefParser.yy"
    { driver.lefDumbMode = 0;
      driver.lefRealNum = 0;
      driver.lefInProp = 0;
    }
#line 6966 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 797: // viarule_prop: "string" "string"
#line 3883 "LefParser.yy"
    {
      if (/*driver.lefrViaRuleCbk*/ 1) {
         char propTp;
         propTp = driver.lefrViaRuleProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrViaRule.lefiViaRule::addProp((*(yystack_[1].value.stringVal)).c_str(), (*(yystack_[0].value.stringVal)).c_str(), propTp);
      }
      delete (yystack_[1].value.stringVal); delete (yystack_[0].value.stringVal);
    }
#line 6979 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 798: // viarule_prop: "string" "qstring"
#line 3892 "LefParser.yy"
    {
      if (/*driver.lefrViaRuleCbk*/ 1) {
         char propTp;
         propTp = driver.lefrViaRuleProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrViaRule.lefiViaRule::addProp((*(yystack_[1].value.stringVal)).c_str(), (*(yystack_[0].value.qstringVal)).c_str(), propTp);
      }
      delete (yystack_[1].value.stringVal); delete (yystack_[0].value.qstringVal);
    }
#line 6992 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 799: // viarule_prop: "string" NUMBER
#line 3901 "LefParser.yy"
    {
      char temp[32];
      sprintf(temp, "%.11g", (yystack_[0].value.doubleVal));
      if (/*driver.lefrViaRuleCbk*/ 1) {
         char propTp;
         propTp = driver.lefrViaRuleProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrViaRule.lefiViaRule::addNumProp((*(yystack_[1].value.stringVal)).c_str(), (yystack_[0].value.doubleVal), temp, propTp);
      }
      delete (yystack_[1].value.stringVal);
    }
#line 7007 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 800: // viarule_layer: viarule_layer_name viarule_layer_options
#line 3913 "LefParser.yy"
    {
      // 10/18/2001 - Wanda da Rosa PCR 404181
      //              Make sure the 1st 2 layers in viarule has direction
      // 04/28/2004 - PCR 704072 - DIRECTION in viarule generate is
      //              obsolete in 5.6
      if (driver.versionNum >= 5.6) {
         if (driver.viaRuleLayer < 2 && !driver.viaRuleHasDir && !driver.viaRuleHasEnc &&
             !driver.isGenerate) {
            if (/*driver.lefrViaRuleCbk*/ 1) {  /* write error only if cbk is set */
               if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
                  driver.lefError(1606, "VIARULE statement in a layer, requires a DIRECTION construct statement.");
                  /*CHKERR();*/ 
               }
            }
         }
      } else {
         if (driver.viaRuleLayer < 2 && !driver.viaRuleHasDir && !driver.viaRuleHasEnc &&
             driver.isGenerate) {
            if (/*driver.lefrViaRuleCbk*/ 1) {  /* write error only if cbk is set */
               if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
                  driver.lefError(1606, "VIARULE statement in a layer, requires a DIRECTION construct statement.");
                  /*CHKERR();*/ 
               }
            }
         }
      }
      driver.viaRuleLayer++;
    }
#line 7040 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 803: // via_name: via_keyword "string" ';'
#line 3949 "LefParser.yy"
    { if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::addViaName((*(yystack_[1].value.stringVal)).c_str()); delete (yystack_[1].value.stringVal);}
#line 7046 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 804: // $@62: %empty
#line 3951 "LefParser.yy"
                            {driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 7052 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 805: // viarule_layer_name: "LAYER" $@62 "string" ';'
#line 3952 "LefParser.yy"
    { if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::setLayer((*(yystack_[1].value.stringVal)).c_str());
      driver.viaRuleHasDir = 0;
      driver.viaRuleHasEnc = 0;
      delete (yystack_[1].value.stringVal);
    }
#line 7062 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 808: // viarule_layer_option: "DIRECTION" "HORIZONTAL" ';'
#line 3965 "LefParser.yy"
    {
      if (driver.viaRuleHasEnc) {
        if (/*driver.lefrViaRuleCbk*/ 1) {  /* write error only if cbk is set */
           if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
              driver.lefError(1607, "An ENCLOSRE statement has already defined in the layer.\nOnly either DIRECTION or ENCLOSURE can be specified in a layer.");
              /*CHKERR();*/
           }
        }
      } else {
        if ((driver.versionNum < 5.6) || (!driver.isGenerate)) {
          if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::setHorizontal();
        } else
          if (/*driver.lefrViaRuleCbk*/ 1)  /* write warning only if cbk is set */
             if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings)
               driver.lefWarning(2022, "DIRECTION statement in VIARULE is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
      }
      driver.viaRuleHasDir = 1;
    }
#line 7085 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 809: // viarule_layer_option: "DIRECTION" "VERTICAL" ';'
#line 3984 "LefParser.yy"
    { 
      if (driver.viaRuleHasEnc) {
        if (/*driver.lefrViaRuleCbk*/ 1) { /* write error only if cbk is set */
           if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
              driver.lefError(1607, "An ENCLOSRE statement has already defined in the layer.\nOnly either DIRECTION or ENCLOSURE can be specified in a layer.");
              /*CHKERR();*/
           }
        }
      } else {
        if ((driver.versionNum < 5.6) || (!driver.isGenerate)) {
          if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::setVertical();
        } else
          if (/*driver.lefrViaRuleCbk*/ 1) /* write warning only if cbk is set */
           if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings)
              driver.lefWarning(2022, "DIRECTION statement in VIARULE is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
      }
      driver.viaRuleHasDir = 1;
    }
#line 7108 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 810: // viarule_layer_option: "ENCLOSURE" NUMBER NUMBER ';'
#line 4003 "LefParser.yy"
    {
      if (driver.versionNum < 5.5) {
         if (/*driver.lefrViaRuleCbk*/ 1) { /* write error only if cbk is set */
           if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                "ENCLOSURE statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1608, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
         }
      }
      // 2/19/2004 - Enforced the rule that ENCLOSURE can only be defined
      // in VIARULE GENERATE
      if (!driver.isGenerate) {
         if (/*driver.lefrViaRuleCbk*/ 1) { /* write error only if cbk is set */
           if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
              driver.lefError(1614, "An ENCLOSURE statement is defined in a VIARULE statement only.\nOVERHANG statement can only be defined in VIARULE GENERATE.");
              /*CHKERR();*/
           }
         }
      }
      if (driver.viaRuleHasDir) {
         if (/*driver.lefrViaRuleCbk*/ 1) { /* write error only if cbk is set */
           if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
              driver.lefError(1609, "A DIRECTION statement has already defined in the layer.\nOnly either DIRECTION or ENCLOSURE can be specified in a layer.");
              /*CHKERR();*/
           }
         }
      } else {
         if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::setEnclosure((yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
      }
      driver.viaRuleHasEnc = 1;
    }
#line 7148 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 811: // viarule_layer_option: "WIDTH" NUMBER "TO" NUMBER ';'
#line 4039 "LefParser.yy"
    { if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::setWidth((yystack_[3].value.doubleVal),(yystack_[1].value.doubleVal)); }
#line 7154 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 812: // viarule_layer_option: "RECT" pt pt ';'
#line 4041 "LefParser.yy"
    { if (/*driver.lefrViaRuleCbk*/ 1)
	driver.lefrViaRule.lefiViaRule::setRect((yystack_[2].value.pt)->x, (yystack_[2].value.pt)->y, (yystack_[1].value.pt)->x, (yystack_[1].value.pt)->y); 
    delete (yystack_[2].value.pt); delete (yystack_[1].value.pt);
    }
#line 7163 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 813: // viarule_layer_option: "SPACING" NUMBER "BY" NUMBER ';'
#line 4046 "LefParser.yy"
    { if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::setSpacing((yystack_[3].value.doubleVal),(yystack_[1].value.doubleVal)); }
#line 7169 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 814: // viarule_layer_option: "RESISTANCE" NUMBER ';'
#line 4048 "LefParser.yy"
    { if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::setResistance((yystack_[1].value.doubleVal)); }
#line 7175 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 815: // viarule_layer_option: "OVERHANG" NUMBER ';'
#line 4050 "LefParser.yy"
    {
      if (!driver.viaRuleHasDir) {
         if (/*driver.lefrViaRuleCbk*/ 1) {  /* write error only if cbk is set */
            if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
               driver.lefError(1610, "An OVERHANG statement is defined, but the required DIRECTION statement is not yet defined.\nUpdate the LEF file to define the DIRECTION statement before the OVERHANG.");
               /*CHKERR();*/
            }
         }
      }
      // 2/19/2004 - Enforced the rule that OVERHANG can only be defined
      // in VIARULE GENERATE after 5.3
      if ((driver.versionNum > 5.3) && (!driver.isGenerate)) {
         if (/*driver.lefrViaRuleCbk*/ 1) {  /* write error only if cbk is set */
            if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
               driver.lefError(1611, "An OVERHANG statement is defined in a VIARULE statement only.\nOVERHANG statement can only be defined in VIARULE GENERATE.");
               /*CHKERR();*/
            }
         }
      }
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::setOverhang((yystack_[1].value.doubleVal));
      } else {
        if (/*driver.lefrViaRuleCbk*/ 1)  /* write warning only if cbk is set */
           if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings)
              driver.lefWarning(2023, "OVERHANG statement will be translated into similar ENCLOSURE rule");
        // In 5.6 & later, set it to either ENCLOSURE overhang1 or overhang2
        if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::setOverhangToEnclosure((yystack_[1].value.doubleVal));
      }
    }
#line 7209 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 816: // viarule_layer_option: "METALOVERHANG" NUMBER ';'
#line 4080 "LefParser.yy"
    {
      // 2/19/2004 - Enforced the rule that METALOVERHANG can only be defined
      // in VIARULE GENERATE
      if ((driver.versionNum > 5.3) && (!driver.isGenerate)) {
         if (/*driver.lefrViaRuleCbk*/ 1) {  /* write error only if cbk is set */
            if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
               driver.lefError(1612, "An METALOVERHANG statement is defined in a VIARULE statement only.\nOVERHANG statement can only be defined in VIARULE GENERATE.");
               /*CHKERR();*/
            }
         }
      }
      if (driver.versionNum < 5.6) {
        if (!driver.viaRuleHasDir) {
           if (/*driver.lefrViaRuleCbk*/ 1) {  /* write error only if cbk is set */
             if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
                driver.lefError(1613, "An METALOVERHANG statement is defined, but the required DIRECTION statement is not yet defined.\nUpdate the LEF file to define the DIRECTION statement before the OVERHANG.");
                /*CHKERR();*/
             } 
           }
        }
        if (/*driver.lefrViaRuleCbk*/ 1) driver.lefrViaRule.lefiViaRule::setMetalOverhang((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrViaRuleCbk*/ 1)  /* write warning only if cbk is set */
           if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings)
             driver.lefWarning(2024, "METALOVERHANG statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 7240 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 817: // $@63: %empty
#line 4107 "LefParser.yy"
                   {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 7246 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 818: // end_viarule: "END" $@63 "string"
#line 4108 "LefParser.yy"
    {
      if (driver.viaRuleName != *(yystack_[0].value.stringVal)) {
        if (/*driver.lefrViaRuleCbk*/ 1) {  /* write error only if cbk is set */
           if (driver.viaRuleWarnings++ < driver.lefrViaRuleWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "END VIARULE name %s is different from the VIARULE name %s.\nCorrect the LEF file before rerun it through the LEF parser.", (*(yystack_[0].value.stringVal)).c_str(), driver.viaRuleName.c_str());
              driver.lefError(1615, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           } 
        } 
      } 
      delete (yystack_[0].value.stringVal);
    }
#line 7266 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 819: // spacing_rule: start_spacing spacings end_spacing
#line 4125 "LefParser.yy"
    { }
#line 7272 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 820: // start_spacing: "SPACING"
#line 4128 "LefParser.yy"
    {
      driver.hasSamenet = 0;
      if ((driver.versionNum < 5.6) || (!driver.ndRule)) {
        /* if 5.6 and in nondefaultrule, it should not get in here, */
        /* it should go to the else statement to write out a warning */
        /* if 5.6, not in nondefaultrule, it will get in here */
        /* if 5.5 and earlier in nondefaultrule is ok to get in here */
        if (driver.versionNum >= 5.7) { /* will get to this if statement if */ 
                           /* driver.versionNum is 5.6 and higher but driver.ndRule = 0 */
           if (driver.spacingWarnings == 0) {  /* only print once */
              driver.lefWarning(2077, "A SPACING SAMENET section is defined but it is not legal in a LEF 5.7 version file.\nIt will be ignored which will probably cause real DRC violations to be ignored, and may\ncause false DRC violations to occur.\n\nTo avoid this warning, and correctly handle these DRC rules, you should modify your\nLEF to use the appropriate SAMENET keywords as described in the LEF/DEF 5.7\nmanual under the SPACING statements in the LAYER (Routing) and LAYER (Cut)\nsections listed in the LEF Table of Contents.");
              driver.spacingWarnings++;
           }
        } else if (/*driver.lefrSpacingBeginCbk*/ 1)
	  driver.lefrSpacingBeginCbk( 0);
      } else
        if (/*driver.lefrSpacingBeginCbk*/ 1)  /* write warning only if cbk is set */
           if (driver.spacingWarnings++ < driver.lefrSpacingWarnings)
             driver.lefWarning(2025, "SAMENET statement in NONDEFAULTRULE is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 7297 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 821: // end_spacing: "END" "SPACING"
#line 4150 "LefParser.yy"
    {
      if ((driver.versionNum < 5.6) || (!driver.ndRule)) {
        if ((driver.versionNum <= 5.4) && (!driver.hasSamenet)) {
           driver.lefError(1616, "SAMENET statement is required inside SPACING for any lef file with version 5.4 and earlier, but is not defined in the parsed lef file.");
           /*CHKERR();*/
        } else if (driver.versionNum < 5.7) { /* obsolete in 5.7 and later */
           if (/*driver.lefrSpacingEndCbk*/ 1)
             driver.lefrSpacingEndCbk( 0);
        }
      }
    }
#line 7313 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 824: // spacing: samenet_keyword "string" "string" NUMBER ';'
#line 4168 "LefParser.yy"
    {
      if ((driver.versionNum < 5.6) || (!driver.ndRule)) {
        if (driver.versionNum < 5.7) {
          if (/*driver.lefrSpacingCbk*/ 1) {
            driver.lefrSpacing.lefiSpacing::set((*(yystack_[3].value.stringVal)).c_str(), (*(yystack_[2].value.stringVal)).c_str(), (yystack_[1].value.doubleVal), 0);
            driver.lefrSpacingCbk( driver.lefrSpacing);
          }
        }
      }
      delete (yystack_[3].value.stringVal); delete (yystack_[2].value.stringVal);
    }
#line 7329 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 825: // spacing: samenet_keyword "string" "string" NUMBER "STACK" ';'
#line 4180 "LefParser.yy"
    {
      if ((driver.versionNum < 5.6) || (!driver.ndRule)) {
        if (driver.versionNum < 5.7) {
          if (/*driver.lefrSpacingCbk*/ 1) {
	    driver.lefrSpacing.lefiSpacing::set((*(yystack_[4].value.stringVal)).c_str(), (*(yystack_[3].value.stringVal)).c_str(), (yystack_[2].value.doubleVal), 1);
	    driver.lefrSpacingCbk( driver.lefrSpacing);
          }
        }
      }
      delete (yystack_[4].value.stringVal); delete (yystack_[3].value.stringVal);
    }
#line 7345 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 826: // samenet_keyword: "SAMENET"
#line 4194 "LefParser.yy"
    { driver.lefDumbMode = 2; driver.lefNoNum = 2; driver.hasSamenet = 1; }
#line 7351 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 827: // irdrop: start_irdrop ir_tables end_irdrop
#line 4197 "LefParser.yy"
    { }
#line 7357 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 828: // start_irdrop: "IRDROP"
#line 4200 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrIRDropBeginCbk*/ 1) 
	  driver.lefrIRDropBeginCbk( 0);
      } else
        if (/*driver.lefrIRDropBeginCbk*/ 1) /* write warning only if cbk is set */
          if (driver.iRDropWarnings++ < driver.lefrIRDropWarnings)
            driver.lefWarning(2026, "IRDROP statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 7371 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 829: // end_irdrop: "END" "IRDROP"
#line 4211 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrIRDropEndCbk*/ 1)
	  driver.lefrIRDropEndCbk( 0);
      }
    }
#line 7382 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 832: // ir_table: ir_tablename ir_table_values ';'
#line 4225 "LefParser.yy"
    { 
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrIRDropCbk*/ 1)
          driver.lefrIRDropCbk( driver.lefrIRDrop);
      }
    }
#line 7393 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 835: // ir_table_value: NUMBER NUMBER
#line 4238 "LefParser.yy"
  { if (/*driver.lefrIRDropCbk*/ 1) driver.lefrIRDrop.lefiIRDrop::setValues((yystack_[1].value.doubleVal), (yystack_[0].value.doubleVal)); }
#line 7399 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 836: // ir_tablename: "TABLE" "string"
#line 4241 "LefParser.yy"
  { if (/*driver.lefrIRDropCbk*/ 1) driver.lefrIRDrop.lefiIRDrop::setTableName((*(yystack_[0].value.stringVal)).c_str()); delete (yystack_[0].value.stringVal); }
#line 7405 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 837: // minfeature: "MINFEATURE" NUMBER NUMBER ';'
#line 4244 "LefParser.yy"
  {
    driver.hasMinfeature = 1;
    if (driver.versionNum < 5.4) {
       if (/*driver.lefrMinFeatureCbk*/ 1) {
         driver.lefrMinFeature.lefiMinFeature::set((yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
         driver.lefrMinFeatureCbk( driver.lefrMinFeature);
       }
    } else
       if (/*driver.lefrMinFeatureCbk*/ 1) /* write warning only if cbk is set */
          if (driver.minFeatureWarnings++ < driver.lefrMinFeatureWarnings)
            driver.lefWarning(2027, "MINFEATURE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
  }
#line 7422 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 838: // dielectric: "DIELECTRIC" NUMBER ';'
#line 4258 "LefParser.yy"
  {
    if (driver.versionNum < 5.4) {
       if (/*lefrDielectricCbk*/ 1)
         driver.lefrDielectricCbk( (yystack_[1].value.doubleVal));
    } else
       if (/*lefrDielectricCbk*/ 1) /* write warning only if cbk is set */
         if (driver.dielectricWarnings++ < driver.lefrDielectricWarnings)
           driver.lefWarning(2028, "DIELECTRIC statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
  }
#line 7436 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 839: // $@64: %empty
#line 4268 "LefParser.yy"
                                  {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 7442 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 840: // $@65: %empty
#line 4269 "LefParser.yy"
  {
    driver.lefSetNonDefault((*(yystack_[0].value.stringVal)).c_str());
    if (/*driver.lefrNonDefaultCbk*/ 1) driver.lefrNonDefault.lefiNonDefault::setName((*(yystack_[0].value.stringVal)).c_str());
    driver.ndLayer = 0;
    driver.ndRule = 1;
    driver.numVia = 0;
    //strcpy(driver.nonDefaultRuleName, $3);
    driver.nonDefaultRuleName = *(yystack_[0].value.stringVal);
    delete (yystack_[0].value.stringVal);
  }
#line 7457 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 841: // $@66: %empty
#line 4280 "LefParser.yy"
           {driver.lefNdRule = 1;}
#line 7463 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 842: // nondefault_rule: "NONDEFAULTRULE" $@64 "string" $@65 nd_hardspacing nd_rules $@66 end_nd_rule
#line 4281 "LefParser.yy"
  {
    // 10/18/2001 - Wanda da Rosa, PCR 404189
    //              At least 1 layer is required
    if ((!driver.ndLayer) && (!driver.lefrRelaxMode)) {
       if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
         if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings) {
            driver.lefError(1617, "NONDEFAULTRULE statement requires at least one LAYER statement.");
            /*CHKERR();*/
         }
       }
    }
    if ((!driver.numVia) && (!driver.lefrRelaxMode) && (driver.versionNum < 5.6)) {
       // VIA is no longer a required statement in 5.6
       if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
         if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings) {
            driver.lefError(1618, "NONDEFAULTRULE statement requires at least one VIA statement.");
            /*CHKERR();*/
         }
       }
    }
    if (/*driver.lefrNonDefaultCbk*/ 1) {
      driver.lefrNonDefault.lefiNonDefault::end();
      driver.lefrNonDefaultCbk( driver.lefrNonDefault);
    }
    driver.ndRule = 0;
    driver.lefDumbMode = 0;
    driver.lefUnsetNonDefault();
  }
#line 7496 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 843: // end_nd_rule: "END"
#line 4311 "LefParser.yy"
    {
      if (!driver.nonDefaultRuleName.empty())
        driver.nonDefaultRuleName = "";
    }
#line 7505 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 844: // end_nd_rule: "END" "string"
#line 4316 "LefParser.yy"
    {
      if (driver.nonDefaultRuleName != *(yystack_[0].value.stringVal)) {
        if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
          if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings) {
             driver.outMsg = (char*)lefMalloc(10000);
             sprintf (driver.outMsg,
                "END NONDEFAULTRULE name %s is different from the NONDEFAULTRULE name %s.\nCorrect the LEF file before rerun it through the LEF parser.", (*(yystack_[0].value.stringVal)).c_str(), driver.nonDefaultRuleName.c_str());
             driver.lefError(1619, driver.outMsg);
             lefFree(driver.outMsg);
             /*CHKERR();*/
          } 
        } 
      } 
      delete (yystack_[0].value.stringVal);
    }
#line 7525 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 846: // nd_hardspacing: "HARDSPACING" ';'
#line 4336 "LefParser.yy"
    {
       if (driver.versionNum < 5.6) {
          if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
            if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "HARDSPACING statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1620, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
          }
       } else 
          if (/*driver.lefrNonDefaultCbk*/ 1)
             driver.lefrNonDefault.lefiNonDefault::setHardspacing();
    }
#line 7546 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 856: // usevia: "USEVIA" "string" ';'
#line 4369 "LefParser.yy"
    {
       if (driver.versionNum < 5.6) {
          if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
             driver.outMsg = (char*)lefMalloc(10000);
             sprintf (driver.outMsg,
               "USEVIA statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
             driver.lefError(1621, driver.outMsg);
             lefFree(driver.outMsg);
             /*CHKERR();*/
          }
       } else {
          if (/*driver.lefrNonDefaultCbk*/ 1)
             driver.lefrNonDefault.lefiNonDefault::addUseVia((*(yystack_[1].value.stringVal)).c_str());
       }
       delete (yystack_[1].value.stringVal);
    }
#line 7567 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 857: // useviarule: "USEVIARULE" "string" ';'
#line 4387 "LefParser.yy"
    {
       if (driver.versionNum < 5.6) {
          if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
             if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings) {
                driver.outMsg = (char*)lefMalloc(10000);
                sprintf (driver.outMsg,
                  "USEVIARULE statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
                driver.lefError(1622, driver.outMsg);
                lefFree(driver.outMsg);
                /*CHKERR();*/
             }
          }
       } else {
          if (/*driver.lefrNonDefaultCbk*/ 1)
             driver.lefrNonDefault.lefiNonDefault::addUseViaRule((*(yystack_[1].value.stringVal)).c_str());
       }
       delete (yystack_[1].value.stringVal);
    }
#line 7590 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 858: // mincuts: "MINCUTS" "string" NUMBER ';'
#line 4407 "LefParser.yy"
    {
       if (driver.versionNum < 5.6) {
          if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
             if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings) {
                driver.outMsg = (char*)lefMalloc(10000);
                sprintf (driver.outMsg,
                  "MINCUTS statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
                driver.lefError(1623, driver.outMsg);
                lefFree(driver.outMsg);
                /*CHKERR();*/
             }
          }
       } else {
          if (/*driver.lefrNonDefaultCbk*/ 1)
             driver.lefrNonDefault.lefiNonDefault::addMinCuts((*(yystack_[2].value.stringVal)).c_str(), (int)(yystack_[1].value.doubleVal));
       }
       delete (yystack_[2].value.stringVal);
    }
#line 7613 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 859: // $@67: %empty
#line 4426 "LefParser.yy"
                    { driver.lefDumbMode = 10000000; driver.lefRealNum = 1; driver.lefInProp = 1; }
#line 7619 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 860: // nd_prop: "PROPERTY" $@67 nd_prop_list ';'
#line 4427 "LefParser.yy"
    { driver.lefDumbMode = 0;
      driver.lefRealNum = 0;
      driver.lefInProp = 0;
    }
#line 7628 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 863: // nd_prop: "string" "string"
#line 4439 "LefParser.yy"
    {
      if (/*driver.lefrNonDefaultCbk*/ 1) {
         char propTp;
         propTp = driver.lefrNondefProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrNonDefault.lefiNonDefault::addProp((*(yystack_[1].value.stringVal)).c_str(), (*(yystack_[0].value.stringVal)).c_str(), propTp);
      }
      delete (yystack_[1].value.stringVal); delete (yystack_[0].value.stringVal);
    }
#line 7641 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 864: // nd_prop: "string" "qstring"
#line 4448 "LefParser.yy"
    {
      if (/*driver.lefrNonDefaultCbk*/ 1) {
         char propTp;
         propTp = driver.lefrNondefProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrNonDefault.lefiNonDefault::addProp((*(yystack_[1].value.stringVal)).c_str(), (*(yystack_[0].value.qstringVal)).c_str(), propTp);
      }
      delete (yystack_[1].value.stringVal); delete (yystack_[0].value.qstringVal);
    }
#line 7654 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 865: // nd_prop: "string" NUMBER
#line 4457 "LefParser.yy"
    {
      if (/*driver.lefrNonDefaultCbk*/ 1) {
         char temp[32];
         char propTp;
         sprintf(temp, "%.11g", (yystack_[0].value.doubleVal));
         propTp = driver.lefrNondefProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrNonDefault.lefiNonDefault::addNumProp((*(yystack_[1].value.stringVal)).c_str(), (yystack_[0].value.doubleVal), temp, propTp);
      }
      delete (yystack_[1].value.stringVal);
    }
#line 7669 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 866: // $@68: %empty
#line 4468 "LefParser.yy"
                  {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 7675 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 867: // $@69: %empty
#line 4469 "LefParser.yy"
  {
    if (/*driver.lefrNonDefaultCbk*/ 1) driver.lefrNonDefault.lefiNonDefault::addLayer((*(yystack_[0].value.stringVal)).c_str());
    driver.ndLayer++;
    //strcpy(driver.layerName, $3);
    driver.layerName = *(yystack_[0].value.stringVal);
    driver.ndLayerWidth = 0;
    driver.ndLayerSpace = 0;
    delete (yystack_[0].value.stringVal);
  }
#line 7689 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 868: // $@70: %empty
#line 4479 "LefParser.yy"
  { 
    driver.ndLayerWidth = 1;
    if (/*driver.lefrNonDefaultCbk*/ 1) driver.lefrNonDefault.lefiNonDefault::addWidth((yystack_[1].value.doubleVal));
  }
#line 7698 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 869: // $@71: %empty
#line 4483 "LefParser.yy"
                       {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 7704 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 870: // nd_layer: "LAYER" $@68 "string" $@69 "WIDTH" NUMBER ';' $@70 nd_layer_stmts "END" $@71 "string"
#line 4484 "LefParser.yy"
  {
    if (driver.layerName != *(yystack_[0].value.stringVal)) {
      if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
         if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings) {
            driver.outMsg = (char*)lefMalloc(10000);
            sprintf (driver.outMsg,
               "END LAYER name %s is different from the LAYER name %s.\nCorrect the LEF file before rerun it through the LEF parser.", (*(yystack_[9].value.stringVal)).c_str(), driver.layerName.c_str());
            driver.lefError(1624, driver.outMsg);
            lefFree(driver.outMsg);
            /*CHKERR();*/
         }
      }
    }
    if (!driver.ndLayerWidth) {
      if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
         if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings) {
            driver.lefError(1625, "A WIDTH statement is required in the LAYER statement in NONDEFULTRULE.");
            /*CHKERR();*/
         }
      }
    }
    if (!driver.ndLayerSpace && driver.versionNum < 5.6) {   // 5.6, SPACING is optional
      if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
         if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings) {
            driver.outMsg = (char*)lefMalloc(10000);
            sprintf (driver.outMsg,
               "A SPACING statement is required in the LAYER statement in NONDEFAULTRULE for lef file with version 5.5 and earlier.\nYour lef file is defined with version %g. Update your lef to add a LAYER statement and try again.",
                driver.versionNum);
            driver.lefError(1626, driver.outMsg);
            lefFree(driver.outMsg);
            /*CHKERR();*/
         }
      }
    }
    delete (yystack_[9].value.stringVal); delete (yystack_[0].value.stringVal);
  }
#line 7745 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 873: // nd_layer_stmt: "SPACING" NUMBER ';'
#line 4529 "LefParser.yy"
    {
      driver.ndLayerSpace = 1;
      if (/*driver.lefrNonDefaultCbk*/ 1) driver.lefrNonDefault.lefiNonDefault::addSpacing((yystack_[1].value.doubleVal));
    }
#line 7754 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 874: // nd_layer_stmt: "WIREEXTENSION" NUMBER ';'
#line 4534 "LefParser.yy"
    { if (/*driver.lefrNonDefaultCbk*/ 1)
         driver.lefrNonDefault.lefiNonDefault::addWireExtension((yystack_[1].value.doubleVal)); }
#line 7761 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 875: // nd_layer_stmt: "RESISTANCE" "RPERSQ" NUMBER ';'
#line 4537 "LefParser.yy"
    {
      if (driver.ignoreVersion) {
         if (/*driver.lefrNonDefaultCbk*/ 1)
            driver.lefrNonDefault.lefiNonDefault::addResistance((yystack_[1].value.doubleVal));
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
            if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "RESISTANCE RPERSQ statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1627, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else if (driver.versionNum > 5.5) {  // obsolete in 5.6
         if (/*driver.lefrNonDefaultCbk*/ 1) /* write warning only if cbk is set */
            if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings)
              driver.lefWarning(2029, "RESISTANCE RPERSQ statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
      } else if (/*driver.lefrNonDefaultCbk*/ 1)
         driver.lefrNonDefault.lefiNonDefault::addResistance((yystack_[1].value.doubleVal));
    }
#line 7788 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 876: // nd_layer_stmt: "CAPACITANCE" "CPERSQDIST" NUMBER ';'
#line 4561 "LefParser.yy"
    {
      if (driver.ignoreVersion) {
         if (/*driver.lefrNonDefaultCbk*/ 1)
            driver.lefrNonDefault.lefiNonDefault::addCapacitance((yystack_[1].value.doubleVal));
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
            if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "CAPACITANCE CPERSQDIST statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1628, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
            }
         }
      } else if (driver.versionNum > 5.5) { // obsolete in 5.6
         if (/*driver.lefrNonDefaultCbk*/ 1) /* write warning only if cbk is set */
            if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings)
              driver.lefWarning(2030, "CAPACITANCE CPERSQDIST statement is obsolete in version 5.6. and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
      } else if (/*driver.lefrNonDefaultCbk*/ 1)
         driver.lefrNonDefault.lefiNonDefault::addCapacitance((yystack_[1].value.doubleVal));
    }
#line 7815 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 877: // nd_layer_stmt: "EDGECAPACITANCE" NUMBER ';'
#line 4584 "LefParser.yy"
    {
      if (driver.ignoreVersion) {
         if (/*driver.lefrNonDefaultCbk*/ 1)
            driver.lefrNonDefault.lefiNonDefault::addEdgeCap((yystack_[1].value.doubleVal));
      } else if (driver.versionNum < 5.4) {
         if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
            if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "EDGECAPACITANCE statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1629, driver.outMsg);
               lefFree(driver.outMsg);
              /*CHKERR();*/
            }
         }
      } else if (driver.versionNum > 5.5) {  // obsolete in 5.6
         if (/*driver.lefrNonDefaultCbk*/ 1) /* write warning only if cbk is set */
            if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings)
              driver.lefWarning(2031, "EDGECAPACITANCE statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
      } else if (/*driver.lefrNonDefaultCbk*/ 1)
         driver.lefrNonDefault.lefiNonDefault::addEdgeCap((yystack_[1].value.doubleVal));
    }
#line 7842 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 878: // nd_layer_stmt: "DIAGWIDTH" NUMBER ';'
#line 4607 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {  // 5.6 syntax
         if (/*driver.lefrNonDefaultCbk*/ 1) { /* write error only if cbk is set */
            if (driver.nonDefaultWarnings++ < driver.lefrNonDefaultWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                 "DIAGWIDTH statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
               driver.lefError(1630, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/ 
            }
         }
      } else {
         if (/*driver.lefrNonDefaultCbk*/ 1)
            driver.lefrNonDefault.lefiNonDefault::addDiagWidth((yystack_[1].value.doubleVal));
      }
    }
#line 7864 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 879: // site: start_site site_options end_site
#line 4626 "LefParser.yy"
    { 
      if (/*driver.lefrSiteCbk*/ 1)
        driver.lefrSiteCbk( driver.lefrSite);
    }
#line 7873 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 880: // $@72: %empty
#line 4631 "LefParser.yy"
                   {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 7879 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 881: // start_site: "SITE" $@72 GSTRING
#line 4632 "LefParser.yy"
    { 
      if (/*driver.lefrSiteCbk*/ 1) driver.lefrSite.lefiSite::setName((*(yystack_[0].value.stringVal)).c_str());
      //strcpy(driver.siteName, $3);
      driver.siteName = (*(yystack_[0].value.stringVal));
      driver.hasSiteClass = 0;
      driver.hasSiteSize = 0;
      driver.hasSite = 1;
      delete (yystack_[0].value.stringVal);
    }
#line 7893 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 882: // $@73: %empty
#line 4642 "LefParser.yy"
                {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 7899 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 883: // end_site: "END" $@73 GSTRING
#line 4643 "LefParser.yy"
    {
      if (driver.siteName != *(yystack_[0].value.stringVal)) {
        if (/*driver.lefrSiteCbk*/ 1) { /* write error only if cbk is set */
           if (driver.siteWarnings++ < driver.lefrSiteWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "END SITE name %s is different from the SITE name %s.\nCorrect the LEF file before rerun it through the LEF parser.", (*(yystack_[0].value.stringVal)).c_str(), driver.siteName.c_str());
              driver.lefError(1631, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           } 
        } 
      } else {
        if (/*driver.lefrSiteCbk*/ 1) { /* write error only if cbk is set */
          if (driver.hasSiteClass == 0) {
             driver.lefError(1632, "A CLASS statement is required in the SITE statement.");
             /*CHKERR();*/
          }
          if (driver.hasSiteSize == 0) {
             driver.lefError(1633, "A SIZE  statement is required in the SITE statement.");
             /*CHKERR();*/
          }
        }
      }
      delete (yystack_[0].value.stringVal);
    }
#line 7930 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 886: // site_option: "SIZE" NUMBER "BY" NUMBER ';'
#line 4677 "LefParser.yy"
    {
/* Workaround for pcr 640902
*/
      if (/*driver.lefrSiteCbk*/ 1) driver.lefrSite.lefiSite::setSize((yystack_[3].value.doubleVal),(yystack_[1].value.doubleVal));
      driver.hasSiteSize = 1;
    }
#line 7941 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 887: // site_option: site_symmetry_statement
#line 4684 "LefParser.yy"
    { }
#line 7947 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 888: // site_option: site_class
#line 4686 "LefParser.yy"
    { 
      if (/*driver.lefrSiteCbk*/ 1) driver.lefrSite.lefiSite::setClass((*(yystack_[0].value.stringVal)).c_str());
      driver.hasSiteClass = 1;
      delete (yystack_[0].value.stringVal);
    }
#line 7957 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 889: // site_option: site_rowpattern_statement
#line 4692 "LefParser.yy"
    { }
#line 7963 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 890: // site_class: "CLASS" "PAD" ';'
#line 4695 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("PAD"); }
#line 7969 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 891: // site_class: "CLASS" "CORE" ';'
#line 4696 "LefParser.yy"
                        {(yylhs.value.stringVal) = new std::string ("CORE"); }
#line 7975 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 892: // site_class: "CLASS" "VIRTUAL" ';'
#line 4697 "LefParser.yy"
                           {(yylhs.value.stringVal) = new std::string ("VIRTUAL"); }
#line 7981 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 893: // site_symmetry_statement: "SYMMETRY" site_symmetries ';'
#line 4700 "LefParser.yy"
    { }
#line 7987 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 896: // site_symmetry: "X"
#line 4709 "LefParser.yy"
    { if (/*driver.lefrSiteCbk*/ 1) driver.lefrSite.lefiSite::setXSymmetry(); }
#line 7993 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 897: // site_symmetry: "Y"
#line 4711 "LefParser.yy"
    { if (/*driver.lefrSiteCbk*/ 1) driver.lefrSite.lefiSite::setYSymmetry(); }
#line 7999 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 898: // site_symmetry: "R90"
#line 4713 "LefParser.yy"
    { if (/*driver.lefrSiteCbk*/ 1) driver.lefrSite.lefiSite::set90Symmetry(); }
#line 8005 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 899: // $@74: %empty
#line 4715 "LefParser.yy"
                                        {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 8011 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 900: // site_rowpattern_statement: "ROWPATTERN" $@74 site_rowpatterns ';'
#line 4717 "LefParser.yy"
    { }
#line 8017 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 903: // $@75: %empty
#line 4724 "LefParser.yy"
                                    {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 8023 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 904: // site_rowpattern: "string" orientation $@75
#line 4725 "LefParser.yy"
    { if (/*driver.lefrSiteCbk*/ 1) driver.lefrSite.lefiSite::addRowPattern((*(yystack_[2].value.stringVal)).c_str(), (yystack_[1].value.integerVal)); delete (yystack_[2].value.stringVal);}
#line 8029 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 905: // pt: NUMBER NUMBER
#line 4729 "LefParser.yy"
    { (yylhs.value.pt) = new lefPoint; (yylhs.value.pt)->x = (yystack_[1].value.doubleVal); (yylhs.value.pt)->y = (yystack_[0].value.doubleVal); }
#line 8035 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 906: // pt: '(' NUMBER NUMBER ')'
#line 4731 "LefParser.yy"
    { (yylhs.value.pt) = new lefPoint; (yylhs.value.pt)->x = (yystack_[2].value.doubleVal); (yylhs.value.pt)->y = (yystack_[1].value.doubleVal); }
#line 8041 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 907: // $@76: %empty
#line 4734 "LefParser.yy"
    { 
      if (/*driver.lefrMacroCbk*/ 1)
        driver.lefrMacroCbk( driver.lefrMacro);
      driver.lefrDoSite = 0;
    }
#line 8051 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 909: // $@77: %empty
#line 4741 "LefParser.yy"
                     {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 8057 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 910: // start_macro: "MACRO" $@77 "string"
#line 4742 "LefParser.yy"
    {
      driver.siteDef = 0;
      driver.symDef = 0;
      driver.sizeDef = 0; 
      driver.pinDef = 0; 
      driver.obsDef = 0; 
      driver.origDef = 0;
      driver.lefrMacro.lefiMacro::clear();      
      if (/*driver.lefrMacroBeginCbk || driver.lefrMacroCbk*/ 1) {
        // some reader may not have MacroBeginCB, but has MacroCB set
        driver.lefrMacro.lefiMacro::setName((*(yystack_[0].value.stringVal)).c_str());
        driver.lefrMacroBeginCbk(*(yystack_[0].value.stringVal));
      }
      //strcpy(driver.macroName, $3);
      driver.macroName = (*(yystack_[0].value.stringVal));
      delete (yystack_[0].value.stringVal);
    }
#line 8079 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 911: // $@78: %empty
#line 4760 "LefParser.yy"
                 {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 8085 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 912: // end_macro: "END" $@78 "string"
#line 4761 "LefParser.yy"
    {
      if (driver.macroName != *(yystack_[0].value.stringVal)) {
        if (/*driver.lefrMacroEndCbk*/ 1) { /* write error only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "END MACRO name %s is different from the MACRO name %s.\nCorrect the LEF file before rerun it through the LEF parser.", (*(yystack_[0].value.stringVal)).c_str(), driver.macroName.c_str());
              driver.lefError(1634, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           } 
        } 
      } 
      if (/*driver.lefrMacroEndCbk*/ 1)
        driver.lefrMacroEndCbk(*(yystack_[0].value.stringVal));
        delete (yystack_[0].value.stringVal);
    }
#line 8107 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 920: // macro_option: macro_origin
#line 4794 "LefParser.yy"
      { }
#line 8113 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 921: // macro_option: macro_power
#line 4796 "LefParser.yy"
      { }
#line 8119 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 922: // macro_option: macro_foreign
#line 4798 "LefParser.yy"
      { }
#line 8125 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 925: // macro_option: macro_size
#line 4802 "LefParser.yy"
      { }
#line 8131 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 926: // macro_option: macro_site
#line 4804 "LefParser.yy"
      { }
#line 8137 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 927: // macro_option: macro_pin
#line 4806 "LefParser.yy"
      { }
#line 8143 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 928: // macro_option: "FUNCTION" "BUFFER" ';'
#line 4808 "LefParser.yy"
      { if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setBuffer(); }
#line 8149 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 929: // macro_option: "FUNCTION" "INVERTER" ';'
#line 4810 "LefParser.yy"
      { if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setInverter(); }
#line 8155 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 930: // macro_option: macro_obs
#line 4812 "LefParser.yy"
      { }
#line 8161 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 931: // macro_option: macro_density
#line 4814 "LefParser.yy"
      { }
#line 8167 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 932: // macro_option: macro_clocktype
#line 4816 "LefParser.yy"
      { }
#line 8173 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 933: // macro_option: timing
#line 4818 "LefParser.yy"
      { }
#line 8179 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 934: // $@79: %empty
#line 4819 "LefParser.yy"
               {driver.lefDumbMode = 1000000; driver.lefRealNum = 1; driver.lefInProp = 1; }
#line 8185 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 935: // macro_option: "PROPERTY" $@79 macro_prop_list ';'
#line 4820 "LefParser.yy"
      { driver.lefDumbMode = 0;
        driver.lefRealNum = 0;
        driver.lefInProp = 0;
      }
#line 8194 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 938: // macro_symmetry_statement: "SYMMETRY" macro_symmetries ';'
#line 4831 "LefParser.yy"
    {
      if (driver.siteDef) { /* SITE is defined before SYMMETRY */
          /* pcr 283846 suppress warning */
          if (/*driver.lefrMacroCbk*/ 1) /* write warning only if cbk is set */
             if (driver.macroWarnings++ < driver.lefrMacroWarnings)
               driver.lefWarning(2032, "A SITE statement is defined before SYMMETRY statement.\nTo avoid this warning in the future, define SITE after SYMMETRY");
      }
      driver.symDef = 1;
    }
#line 8208 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 941: // macro_symmetry: "X"
#line 4848 "LefParser.yy"
    { if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setXSymmetry(); }
#line 8214 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 942: // macro_symmetry: "Y"
#line 4850 "LefParser.yy"
    { if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setYSymmetry(); }
#line 8220 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 943: // macro_symmetry: "R90"
#line 4852 "LefParser.yy"
    { if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::set90Symmetry(); }
#line 8226 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 944: // macro_name_value_pair: "string" NUMBER
#line 4856 "LefParser.yy"
    {
      char temp[32];
      sprintf(temp, "%.11g", (yystack_[0].value.doubleVal));
      if (/*driver.lefrMacroCbk*/ 1) {
         char propTp;
         propTp = driver.lefrMacroProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrMacro.lefiMacro::setNumProperty((*(yystack_[1].value.stringVal)).c_str(), (yystack_[0].value.doubleVal), temp,  propTp);
      }
      delete (yystack_[1].value.stringVal);
    }
#line 8241 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 945: // macro_name_value_pair: "string" "qstring"
#line 4867 "LefParser.yy"
    {
      if (/*driver.lefrMacroCbk*/ 1) {
         char propTp;
         propTp = driver.lefrMacroProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrMacro.lefiMacro::setProperty((*(yystack_[1].value.stringVal)).c_str(), (*(yystack_[0].value.qstringVal)).c_str(), propTp);
      }
      delete (yystack_[1].value.stringVal); delete (yystack_[0].value.qstringVal);
    }
#line 8254 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 946: // macro_name_value_pair: "string" "string"
#line 4876 "LefParser.yy"
    {
      if (/*driver.lefrMacroCbk*/ 1) {
         char propTp;
         propTp = driver.lefrMacroProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrMacro.lefiMacro::setProperty((*(yystack_[1].value.stringVal)).c_str(), (*(yystack_[0].value.stringVal)).c_str(), propTp);
      }
      delete (yystack_[1].value.stringVal); delete (yystack_[0].value.stringVal);
    }
#line 8267 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 947: // macro_class: "CLASS" class_type ';'
#line 4886 "LefParser.yy"
    {
       if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setClass((*(yystack_[1].value.stringVal)).c_str());
       if (/*driver.lefrMacroClassTypeCbk*/ 1)
          driver.lefrMacroClassTypeCbk(*(yystack_[1].value.stringVal));
        delete (yystack_[1].value.stringVal);
    }
#line 8278 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 948: // class_type: "COVER"
#line 4894 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("COVER"); }
#line 8284 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 949: // class_type: "COVER" "BUMP"
#line 4896 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("COVER BUMP");
      if (driver.versionNum < 5.5) {
        if (/*driver.lefrMacroCbk*/ 1) { /* write error only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              if (driver.lefrRelaxMode)
                 driver.lefWarning(2033, "The statement COVER BUMP is a LEF verion 5.5 syntax.\nYour LEF file is version 5.4 or earlier which is illegal but will be allowed\nbecause this application does not enforce strict version checking.\nOther tools that enforce strict checking will have a syntax error when reading this file.\nYou can change the VERSION statement in this LEF file to 5.5 or higher to stop this warning.");
              else {
                 driver.outMsg = (char*)lefMalloc(10000);
                 sprintf (driver.outMsg,
                    "COVER BUMP statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
                 driver.lefError(1635, driver.outMsg);
                 lefFree(driver.outMsg);
                 /*CHKERR();*/
              }
           }
        }
      }
    }
#line 8307 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 950: // class_type: "RING"
#line 4914 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("RING"); }
#line 8313 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 951: // class_type: "BLOCK"
#line 4915 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("BLOCK"); }
#line 8319 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 952: // class_type: "BLOCK" "BLACKBOX"
#line 4917 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("BLOCK BLACKBOX");
      if (driver.versionNum < 5.5) {
        if (/*driver.lefrMacroCbk*/ 1) { /* write error only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
             if (driver.lefrRelaxMode)
                driver.lefWarning(2034, "The statement BLOCK BLACKBOX is a LEF verion 5.5 syntax.\nYour LEF file is version 5.4 or earlier which is illegal but will be allowed\nbecause this application does not enforce strict version checking.\nOther tools that enforce strict checking will have a syntax error when reading this file.\nYou can change the VERSION statement in this LEF file to 5.5 or higher to stop this warning.");
              else {
                 driver.outMsg = (char*)lefMalloc(10000);
                 sprintf (driver.outMsg,
                    "BLOCK BLACKBOX statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
                 driver.lefError(1636, driver.outMsg);
                 lefFree(driver.outMsg);
                 /*CHKERR();*/
              }
           }
        }
      }
    }
#line 8342 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 953: // class_type: "BLOCK" "SOFT"
#line 4936 "LefParser.yy"
    {
      if (driver.ignoreVersion) {
        (yylhs.value.stringVal) = new std::string ("BLOCK SOFT");
      } else if (driver.versionNum < 5.6) {
        if (/*driver.lefrMacroCbk*/ 1) { /* write error only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "BLOCK SOFT statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1637, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
      else
        (yylhs.value.stringVal) = new std::string ("BLOCK SOFT");
    }
#line 8365 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 954: // class_type: "NONE"
#line 4954 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("NONE"); }
#line 8371 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 955: // class_type: "BUMP"
#line 4956 "LefParser.yy"
      {
        if (driver.versionNum < 5.7) {
          driver.outMsg = (char*)lefMalloc(10000);
          sprintf(driver.outMsg,
            "BUMP is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
          driver.lefError(1698, driver.outMsg);
          lefFree(driver.outMsg);
          /*CHKERR();*/
      } else
        (yylhs.value.stringVal) = new std::string ("BUMP");
      }
#line 8387 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 956: // class_type: "PAD"
#line 4967 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("PAD"); }
#line 8393 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 957: // class_type: "VIRTUAL"
#line 4968 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("VIRTUAL"); }
#line 8399 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 958: // class_type: "PAD" pad_type
#line 4970 "LefParser.yy"
      {  sprintf(driver.temp_name, "PAD %s", (*(yystack_[0].value.stringVal)).c_str());
        (yylhs.value.stringVal) = new std::string (driver.temp_name); 
        if (driver.versionNum < 5.5) {
           if ("AREAIO" != *(yystack_[0].value.stringVal)) {
             sprintf(driver.temp_name, "PAD %s", (*(yystack_[0].value.stringVal)).c_str());
			 (yylhs.value.stringVal) = new std::string (driver.temp_name); 
           } else if (/*driver.lefrMacroCbk*/ 1) { 
             if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
               if (driver.lefrRelaxMode)
                  driver.lefWarning(2035, "The statement PAD AREAIO is a LEF verion 5.5 syntax.\nYour LEF file is version 5.4 or earlier which is illegal but will be allowed\nbecause this application does not enforce strict version checking.\nOther tools that enforce strict checking will have a syntax error when reading this file.\nYou can change the VERSION statement in this LEF file to 5.5 or higher to stop this warning.");
               else {
                  driver.outMsg = (char*)lefMalloc(10000);
                  sprintf (driver.outMsg,
                     "PAD AREAIO statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
                  driver.lefError(1638, driver.outMsg);
                  lefFree(driver.outMsg);
                  /*CHKERR();*/
               }
            }
          }
        }
        delete (yystack_[0].value.stringVal);
      }
#line 8427 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 959: // class_type: "CORE"
#line 4993 "LefParser.yy"
              {(yylhs.value.stringVal) = new std::string ("CORE"); }
#line 8433 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 960: // class_type: "CORNER"
#line 4995 "LefParser.yy"
      {(yylhs.value.stringVal) = new std::string ("CORNER");
      /* This token is NOT in the spec but has shown up in 
       * some lef files.  This exception came from LEFOUT
       * in 'frameworks'
       */
      }
#line 8444 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 961: // class_type: "CORE" core_type
#line 5002 "LefParser.yy"
      {sprintf(driver.temp_name, "CORE %s", (*(yystack_[0].value.stringVal)).c_str()); delete (yystack_[0].value.stringVal);
      (yylhs.value.stringVal) = new std::string (driver.temp_name);}
#line 8451 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 962: // class_type: "ENDCAP" endcap_type
#line 5005 "LefParser.yy"
      {sprintf(driver.temp_name, "ENDCAP %s", (*(yystack_[0].value.stringVal)).c_str()); delete (yystack_[0].value.stringVal);
      (yylhs.value.stringVal) = new std::string (driver.temp_name);}
#line 8458 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 963: // pad_type: "INPUT"
#line 5010 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("INPUT");}
#line 8464 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 964: // pad_type: "OUTPUT"
#line 5011 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("OUTPUT");}
#line 8470 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 965: // pad_type: "INOUT"
#line 5012 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("INOUT");}
#line 8476 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 966: // pad_type: "POWER"
#line 5013 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("POWER");}
#line 8482 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 967: // pad_type: "SPACER"
#line 5014 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("SPACER");}
#line 8488 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 968: // pad_type: "AREAIO"
#line 5015 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("AREAIO");}
#line 8494 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 969: // core_type: "FEEDTHRU"
#line 5018 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("FEEDTHRU");}
#line 8500 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 970: // core_type: "TIEHIGH"
#line 5019 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("TIEHIGH");}
#line 8506 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 971: // core_type: "TIELOW"
#line 5020 "LefParser.yy"
                    {(yylhs.value.stringVal) = new std::string ("TIELOW");}
#line 8512 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 972: // core_type: "SPACER"
#line 5022 "LefParser.yy"
    { 
      if (driver.ignoreVersion) {
        (yylhs.value.stringVal) = new std::string ("SPACER");
      } else if (driver.versionNum < 5.4) {
        if (/*driver.lefrMacroCbk*/ 1) { /* write error only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "SPACER statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1639, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
      else
        (yylhs.value.stringVal) = new std::string ("SPACER");
    }
#line 8535 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 973: // core_type: "ANTENNACELL"
#line 5041 "LefParser.yy"
    { 
      if (driver.ignoreVersion) {
        (yylhs.value.stringVal) = new std::string ("ANTENNACELL");
      } else if (driver.versionNum < 5.4) {
        if (/*driver.lefrMacroCbk*/ 1) { /* write error only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNACELL statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1640, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
      else
        (yylhs.value.stringVal) = new std::string ("ANTENNACELL");
    }
#line 8558 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 974: // core_type: "WELLTAP"
#line 5060 "LefParser.yy"
    { 
      if (driver.ignoreVersion) {
        (yylhs.value.stringVal) = new std::string ("WELLTAP");
      } else if (driver.versionNum < 5.6) {
        if (/*driver.lefrMacroCbk*/ 1) { /* write error only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "WELLTAP statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1641, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
      else
        (yylhs.value.stringVal) = new std::string ("WELLTAP");
    }
#line 8581 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 975: // endcap_type: "PRE"
#line 5080 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("PRE");}
#line 8587 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 976: // endcap_type: "POST"
#line 5081 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("POST");}
#line 8593 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 977: // endcap_type: "TOPLEFT"
#line 5082 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("TOPLEFT");}
#line 8599 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 978: // endcap_type: "TOPRIGHT"
#line 5083 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("TOPRIGHT");}
#line 8605 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 979: // endcap_type: "BOTTOMLEFT"
#line 5084 "LefParser.yy"
                        {(yylhs.value.stringVal) = new std::string ("BOTTOMLEFT");}
#line 8611 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 980: // endcap_type: "BOTTOMRIGHT"
#line 5085 "LefParser.yy"
                        {(yylhs.value.stringVal) = new std::string ("BOTTOMRIGHT");}
#line 8617 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 981: // macro_generator: "GENERATOR" "string" ';'
#line 5088 "LefParser.yy"
    { if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setGenerator((*(yystack_[1].value.stringVal)).c_str()); delete (yystack_[1].value.stringVal);}
#line 8623 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 982: // macro_generate: "GENERATE" "string" "string" ';'
#line 5091 "LefParser.yy"
    { if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setGenerate((*(yystack_[2].value.stringVal)).c_str(), (*(yystack_[1].value.stringVal)).c_str());  delete (yystack_[2].value.stringVal); delete (yystack_[1].value.stringVal);}
#line 8629 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 983: // macro_source: "SOURCE" "USER" ';'
#line 5095 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setSource("USER");
      } else
        if (/*driver.lefrMacroCbk*/ 1) /* write warning only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings)
             driver.lefWarning(2036, "SOURCE statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 8642 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 984: // macro_source: "SOURCE" "GENERATE" ';'
#line 5104 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setSource("GENERATE");
      } else
        if (/*driver.lefrMacroCbk*/ 1) /* write warning only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings)
             driver.lefWarning(2037, "SOURCE statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 8655 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 985: // macro_source: "SOURCE" "BLOCK" ';'
#line 5113 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setSource("BLOCK");
      } else
        if (/*driver.lefrMacroCbk*/ 1) /* write warning only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings)
             driver.lefWarning(2037, "SOURCE statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    }
#line 8668 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 986: // macro_power: "POWER" NUMBER ';'
#line 5123 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setPower((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrMacroCbk*/ 1) /* write warning only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings)
             driver.lefWarning(2038, "MACRO POWER statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 8681 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 987: // macro_origin: "ORIGIN" pt ';'
#line 5133 "LefParser.yy"
    { 
       if (driver.origDef) { /* Has multiple ORIGIN defined in a macro, stop parsing*/
          if (/*driver.lefrMacroCbk*/ 1) { /* write error only if cbk is set */
             if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
                driver.lefError(1642, "ORIGIN statement has defined more than once in a MACRO statement.\nOnly one ORIGIN statement can be defined in a Macro.\nParser stops executions.");
               /*CHKERR();*/
             }
          }
       }
       driver.origDef = 1;
       if (driver.siteDef) { /* SITE is defined before ORIGIN */
          /* pcr 283846 suppress warning */
          if (/*driver.lefrMacroCbk*/ 1) /* write warning only if cbk is set */
             if (driver.macroWarnings++ < driver.lefrMacroWarnings)
               driver.lefWarning(2039, "A SITE statement is defined before ORIGIN statement.\nTo avoid this warning in the future, define SITE after ORIGIN");
       }
       if (driver.pinDef) { /* PIN is defined before ORIGIN */
          /* pcr 283846 suppress warning */
          if (/*driver.lefrMacroCbk*/ 1) /* write warning only if cbk is set */
             if (driver.macroWarnings++ < driver.lefrMacroWarnings)
               driver.lefWarning(2040, "A PIN statement is defined before ORIGIN statement.\nTo avoid this warning in the future, define PIN after ORIGIN");
       }
       if (driver.obsDef) { /* OBS is defined before ORIGIN */
          /* pcr 283846 suppress warning */
          if (/*driver.lefrMacroCbk*/ 1) /* write warning only if cbk is set */
             if (driver.macroWarnings++ < driver.lefrMacroWarnings)
               driver.lefWarning(2041, "A OBS statement is defined before ORIGIN statement.\nTo avoid this warning in the future, define OBS after ORIGIN");
       }
       /* 11/22/99 - Wanda da Rosa. PCR 283846 
                     can be defined any order.
       if (driver.symDef)  * SYMMETRY is defined before ORIGIN *
          driver.lefWarning("SYMMETRY is defined before ORIGIN.");
       */
       /* Add back it back in per Nora request on PCR 283846 */
       /* 1/14/2000 - Wanda da Rosa, PCR 288770
       if (driver.sizeDef)  * SIZE is defined before ORIGIN *
          driver.lefWarning("SIZE is defined before ORIGIN.");
       */
      
       /* Workaround for pcr 640902 */
       if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setOrigin((yystack_[1].value.pt)->x, (yystack_[1].value.pt)->y);
       if (/*driver.lefrMacroOriginCbk*/ 1) {
          driver.macroNum.x = (yystack_[1].value.pt)->x; 
          driver.macroNum.y = (yystack_[1].value.pt)->y; 
          driver.lefrMacroOriginCbk( driver.macroNum);
       }
       delete (yystack_[1].value.pt);
    }
#line 8734 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 988: // macro_foreign: start_foreign ';'
#line 5184 "LefParser.yy"
    { if (/*driver.lefrMacroCbk*/ 1)
      driver.lefrMacro.lefiMacro::addForeign((*(yystack_[1].value.stringVal)).c_str(), 0, 0.0, 0.0, -1);
      delete (yystack_[1].value.stringVal);
    }
#line 8743 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 989: // macro_foreign: start_foreign pt ';'
#line 5189 "LefParser.yy"
    { if (/*driver.lefrMacroCbk*/ 1)
      driver.lefrMacro.lefiMacro::addForeign((*(yystack_[2].value.stringVal)).c_str(), 1, (yystack_[1].value.pt)->x, (yystack_[1].value.pt)->y, -1);
      delete (yystack_[2].value.stringVal); delete (yystack_[1].value.pt);
    }
#line 8752 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 990: // macro_foreign: start_foreign pt orientation ';'
#line 5194 "LefParser.yy"
    { if (/*driver.lefrMacroCbk*/ 1)
      driver.lefrMacro.lefiMacro::addForeign((*(yystack_[3].value.stringVal)).c_str(), 1, (yystack_[2].value.pt)->x, (yystack_[2].value.pt)->y, (yystack_[1].value.integerVal));
      delete (yystack_[3].value.stringVal); delete (yystack_[2].value.pt);
    }
#line 8761 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 991: // macro_foreign: start_foreign orientation ';'
#line 5199 "LefParser.yy"
    { if (/*driver.lefrMacroCbk*/ 1)
      driver.lefrMacro.lefiMacro::addForeign((*(yystack_[2].value.stringVal)).c_str(), 0, 0.0, 0.0, (yystack_[1].value.integerVal));
      delete (yystack_[2].value.stringVal);
    }
#line 8770 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 992: // $@80: %empty
#line 5204 "LefParser.yy"
                 { driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 8776 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 993: // macro_eeq: "EEQ" $@80 "string" ';'
#line 5205 "LefParser.yy"
    { if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setEEQ((*(yystack_[1].value.stringVal)).c_str()); delete (yystack_[1].value.stringVal);}
#line 8782 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 994: // $@81: %empty
#line 5207 "LefParser.yy"
                 { driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 8788 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 995: // macro_leq: "LEQ" $@81 "string" ';'
#line 5208 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setLEQ((*(yystack_[1].value.stringVal)).c_str());
      } else
        if (/*driver.lefrMacroCbk*/ 1) /* write warning only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings)
             driver.lefWarning(2042, "LEQ statement in MACRO is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
     delete (yystack_[1].value.stringVal);
    }
#line 8802 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 996: // macro_site: macro_site_word GSTRING ';'
#line 5220 "LefParser.yy"
    {
      if (/*driver.lefrMacroCbk*/ 1) {
        driver.lefrMacro.lefiMacro::setSiteName((*(yystack_[1].value.stringVal)).c_str());
/* For later than 5.6 release
        driver.lefrSitePatternPtr = (lefiSitePattern*)lefMalloc(
                              sizeof(lefiSitePattern));
        driver.lefrSitePatternPtr->lefiSitePattern::Init();
        driver.lefrSitePatternPtr->lefiSitePattern::setSiteName($2);
*/
      }
      delete (yystack_[1].value.stringVal);
    }
#line 8819 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 997: // macro_site: macro_site_word sitePattern ';'
#line 5233 "LefParser.yy"
    {
      if (/*driver.lefrMacroCbk*/ 1) {
        /* also set site name in the variable driver.siteName_ in lefiMacro */
        /* this, if user wants to use method driver.siteName will get the name also */
        /* Does not work, it will only set with the last sitename, if multiple
           SITEs are defined.  Therefore, data will not be correct
        driver.lefrMacro.lefiMacro::setSitePatternName(
            driver.lefrSitePatternPtr->lefiSitePattern::name());
        */
	driver.lefrMacro.lefiMacro::setSitePattern(driver.lefrSitePatternPtr);
	driver.lefrSitePatternPtr = 0;
      }
    }
#line 8837 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 998: // macro_site_word: "SITE"
#line 5248 "LefParser.yy"
    { driver.lefDumbMode = 1; driver.lefNoNum = 1; driver.siteDef = 1;
        if (/*driver.lefrMacroCbk*/ 1) driver.lefrDoSite = 1; }
#line 8844 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 999: // site_word: "SITE"
#line 5252 "LefParser.yy"
    { driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 8850 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1000: // macro_size: "SIZE" NUMBER "BY" NUMBER ';'
#line 5255 "LefParser.yy"
    { 
      if (driver.siteDef) { /* SITE is defined before SIZE */
         /* pcr 283846 suppress warning
         if (driver.siteWarnings++ < driver.lefrSiteWarnings)
           driver.lefWarning("SITE is defined before SIZE.");
         return 1; 
         */
      }
      driver.sizeDef = 1;
      if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setSize((yystack_[3].value.doubleVal), (yystack_[1].value.doubleVal));
      if (/*driver.lefrMacroSizeCbk*/ 1) {
         driver.macroNum.x = (yystack_[3].value.doubleVal); 
         driver.macroNum.y = (yystack_[1].value.doubleVal); 
         driver.lefrMacroSizeCbk( driver.macroNum);
      }
    }
#line 8871 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1001: // macro_pin: start_macro_pin macro_pin_options end_macro_pin
#line 5276 "LefParser.yy"
    { 
      if (/*driver.lefrPinCbk*/ 1)
        /*driver.lefrPinCbk( driver.lefrPin);*/
		driver.lefrMacro.addPin(driver.lefrPin);
      driver.lefrPin.lefiPin::clear();
    }
#line 8882 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1002: // $@82: %empty
#line 5283 "LefParser.yy"
                       {driver.lefDumbMode = 1; driver.lefNoNum = 1; driver.pinDef = 1;}
#line 8888 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1003: // start_macro_pin: "PIN" $@82 GSTRING
#line 5284 "LefParser.yy"
    { if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setName((*(yystack_[0].value.stringVal)).c_str());
      //strcpy(driver.pinName, $3);
      driver.pinName = (*(yystack_[0].value.stringVal));
      delete (yystack_[0].value.stringVal);
    }
#line 8898 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1004: // $@83: %empty
#line 5290 "LefParser.yy"
                     {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 8904 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1005: // end_macro_pin: "END" $@83 GSTRING
#line 5291 "LefParser.yy"
    {
      if (driver.pinName != *(yystack_[0].value.stringVal) ) {
        if (/*driver.lefrMacroCbk*/ 1) { /* write error only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "END PIN name %s is different from the PIN name %s.\nCorrect the LEF file before rerun it through the LEF parser.", (*(yystack_[0].value.stringVal)).c_str(), driver.pinName.c_str());
              driver.lefError(1643, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           } 
        } 
      } 
      delete (yystack_[0].value.stringVal);
    }
#line 8924 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1006: // macro_pin_options: %empty
#line 5309 "LefParser.yy"
    { }
#line 8930 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1007: // macro_pin_options: macro_pin_options macro_pin_option
#line 5311 "LefParser.yy"
    { }
#line 8936 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1008: // macro_pin_option: start_foreign ';'
#line 5315 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addForeign((*(yystack_[1].value.stringVal)).c_str(), 0, 0.0, 0.0, -1);
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2043, "FOREIGN statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
     delete (yystack_[1].value.stringVal);
    }
#line 8950 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1009: // macro_pin_option: start_foreign pt ';'
#line 5325 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addForeign((*(yystack_[2].value.stringVal)).c_str(), 1, (yystack_[1].value.pt)->x, (yystack_[1].value.pt)->y, -1);
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2043, "FOREIGN statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    delete (yystack_[2].value.stringVal); delete (yystack_[1].value.pt);
    }
#line 8964 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1010: // macro_pin_option: start_foreign pt orientation ';'
#line 5335 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addForeign((*(yystack_[3].value.stringVal)).c_str(), 1, (yystack_[2].value.pt)->x, (yystack_[2].value.pt)->y, (yystack_[1].value.integerVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2043, "FOREIGN statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    delete (yystack_[3].value.stringVal); delete (yystack_[2].value.pt);
    }
#line 8978 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1011: // macro_pin_option: start_foreign "STRUCTURE" ';'
#line 5345 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addForeign((*(yystack_[2].value.stringVal)).c_str(), 0, 0.0, 0.0, -1);
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2043, "FOREIGN statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    delete (yystack_[2].value.stringVal);
    }
#line 8992 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1012: // macro_pin_option: start_foreign "STRUCTURE" pt ';'
#line 5355 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addForeign((*(yystack_[3].value.stringVal)).c_str(), 1, (yystack_[1].value.pt)->x, (yystack_[1].value.pt)->y, -1);
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2043, "FOREIGN statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    delete (yystack_[3].value.stringVal); delete (yystack_[1].value.pt);
    }
#line 9006 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1013: // macro_pin_option: start_foreign "STRUCTURE" pt orientation ';'
#line 5365 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addForeign((*(yystack_[4].value.stringVal)).c_str(), 1, (yystack_[2].value.pt)->x, (yystack_[2].value.pt)->y, (yystack_[1].value.integerVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2043, "FOREIGN statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    delete (yystack_[4].value.stringVal); delete (yystack_[2].value.pt);
    }
#line 9020 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1014: // $@84: %empty
#line 5374 "LefParser.yy"
          { driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 9026 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1015: // macro_pin_option: "LEQ" $@84 "string" ';'
#line 5375 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setLEQ((*(yystack_[1].value.stringVal)).c_str());
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2044, "LEQ statement in MACRO PIN is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    delete (yystack_[1].value.stringVal);
   }
#line 9040 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1016: // macro_pin_option: "POWER" NUMBER ';'
#line 5385 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setPower((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2045, "MACRO POWER statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9053 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1017: // macro_pin_option: electrical_direction
#line 5394 "LefParser.yy"
    { if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setDirection((*(yystack_[0].value.stringVal)).c_str()); delete (yystack_[0].value.stringVal);}
#line 9059 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1018: // macro_pin_option: "USE" macro_pin_use ';'
#line 5396 "LefParser.yy"
    { if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setUse((*(yystack_[1].value.stringVal)).c_str()); delete (yystack_[1].value.stringVal);}
#line 9065 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1019: // macro_pin_option: "SCANUSE" macro_scan_use ';'
#line 5398 "LefParser.yy"
    { }
#line 9071 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1020: // macro_pin_option: "LEAKAGE" NUMBER ';'
#line 5400 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setLeakage((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2046, "MACRO LEAKAGE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, r emove this statement from the LEF file with version 5.4 or later.");
    }
#line 9084 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1021: // macro_pin_option: "RISETHRESH" NUMBER ';'
#line 5409 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setRiseThresh((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2047, "MACRO RISETHRESH statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9097 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1022: // macro_pin_option: "FALLTHRESH" NUMBER ';'
#line 5418 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setFallThresh((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2048, "MACRO FALLTHRESH statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9110 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1023: // macro_pin_option: "RISESATCUR" NUMBER ';'
#line 5427 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setRiseSatcur((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2049, "MACRO RISESATCUR statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9123 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1024: // macro_pin_option: "FALLSATCUR" NUMBER ';'
#line 5436 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setFallSatcur((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2050, "MACRO FALLSATCUR statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9136 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1025: // macro_pin_option: "VLO" NUMBER ';'
#line 5445 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setVLO((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2051, "MACRO VLO statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9149 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1026: // macro_pin_option: "VHI" NUMBER ';'
#line 5454 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setVHI((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2052, "MACRO VHI statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9162 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1027: // macro_pin_option: "TIEOFFR" NUMBER ';'
#line 5463 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setTieoffr((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2053, "MACRO TIEOFFR statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9175 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1028: // macro_pin_option: "SHAPE" pin_shape ';'
#line 5472 "LefParser.yy"
    { if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setShape((*(yystack_[1].value.stringVal)).c_str()); delete (yystack_[1].value.stringVal);}
#line 9181 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1029: // $@85: %empty
#line 5473 "LefParser.yy"
               {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 9187 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1030: // macro_pin_option: "MUSTJOIN" $@85 "string" ';'
#line 5474 "LefParser.yy"
    { if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setMustjoin((*(yystack_[1].value.stringVal)).c_str()); delete (yystack_[1].value.stringVal);}
#line 9193 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1031: // $@86: %empty
#line 5475 "LefParser.yy"
                        {driver.lefDumbMode = 1;}
#line 9199 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1032: // macro_pin_option: "OUTPUTNOISEMARGIN" $@86 NUMBER NUMBER ';'
#line 5476 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setOutMargin((yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2054, "MACRO OUTPUTNOISEMARGIN statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9212 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1033: // $@87: %empty
#line 5484 "LefParser.yy"
                       {driver.lefDumbMode = 1;}
#line 9218 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1034: // macro_pin_option: "OUTPUTRESISTANCE" $@87 NUMBER NUMBER ';'
#line 5485 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setOutResistance((yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2055, "MACRO OUTPUTRESISTANCE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9231 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1035: // $@88: %empty
#line 5493 "LefParser.yy"
                       {driver.lefDumbMode = 1;}
#line 9237 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1036: // macro_pin_option: "INPUTNOISEMARGIN" $@88 NUMBER NUMBER ';'
#line 5494 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setInMargin((yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2056, "MACRO INPUTNOISEMARGIN statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9250 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1037: // macro_pin_option: "CAPACITANCE" NUMBER ';'
#line 5503 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setCapacitance((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2057, "MACRO CAPACITANCE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9263 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1038: // macro_pin_option: "MAXDELAY" NUMBER ';'
#line 5512 "LefParser.yy"
    { if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setMaxdelay((yystack_[1].value.doubleVal)); }
#line 9269 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1039: // macro_pin_option: "MAXLOAD" NUMBER ';'
#line 5514 "LefParser.yy"
    { if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setMaxload((yystack_[1].value.doubleVal)); }
#line 9275 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1040: // macro_pin_option: "RESISTANCE" NUMBER ';'
#line 5516 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setResistance((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2058, "MACRO RESISTANCE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9288 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1041: // macro_pin_option: "PULLDOWNRES" NUMBER ';'
#line 5525 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setPulldownres((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2059, "MACRO PULLDOWNRES statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9301 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1042: // macro_pin_option: "CURRENTSOURCE" "ACTIVE" ';'
#line 5534 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setCurrentSource("ACTIVE");
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2060, "MACRO CURRENTSOURCE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9314 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1043: // macro_pin_option: "CURRENTSOURCE" "RESISTIVE" ';'
#line 5543 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setCurrentSource("RESISTIVE");
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2061, "MACRO CURRENTSOURCE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9327 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1044: // macro_pin_option: "RISEVOLTAGETHRESHOLD" NUMBER ';'
#line 5552 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setRiseVoltage((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2062, "MACRO RISEVOLTAGETHRESHOLD statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9340 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1045: // macro_pin_option: "FALLVOLTAGETHRESHOLD" NUMBER ';'
#line 5561 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setFallVoltage((yystack_[1].value.doubleVal));
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2063, "MACRO FALLVOLTAGETHRESHOLD statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 9353 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1046: // macro_pin_option: "IV_TABLES" "string" "string" ';'
#line 5570 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setTables((*(yystack_[2].value.stringVal)).c_str(), (*(yystack_[1].value.stringVal)).c_str());
      } else
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2064, "MACRO IV_TABLES statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    delete (yystack_[2].value.stringVal); delete (yystack_[1].value.stringVal);
    }
#line 9367 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1047: // macro_pin_option: "TAPERRULE" "string" ';'
#line 5580 "LefParser.yy"
    { if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setTaperRule((*(yystack_[1].value.stringVal)).c_str()); delete (yystack_[1].value.stringVal);}
#line 9373 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1048: // $@89: %empty
#line 5581 "LefParser.yy"
               {driver.lefDumbMode = 1000000; driver.lefRealNum = 1; driver.lefInProp = 1; }
#line 9379 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1049: // macro_pin_option: "PROPERTY" $@89 pin_prop_list ';'
#line 5582 "LefParser.yy"
    { driver.lefDumbMode = 0;
      driver.lefRealNum = 0;
      driver.lefInProp = 0;
    }
#line 9388 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1050: // macro_pin_option: start_macro_port macro_port_class_option geometries "END"
#line 5587 "LefParser.yy"
    {
      driver.lefDumbMode = 0;
      driver.hasGeoLayer = 0;
      if (/*driver.lefrPinCbk*/ 1) {
	driver.lefrPin.lefiPin::addPort(driver.lefrGeometriesPtr);
	driver.lefrGeometriesPtr = 0;
	driver.lefrDoGeometries = 0;
      }
      if ((driver.needGeometry) && (driver.needGeometry != 2))  // if the last LAYER in PORT
        if (/*driver.lefrPinCbk*/ 1) /* write warning only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings)
             driver.lefWarning(2065, "Either PATH, RECT or POLYGON statement is a required in MACRO/PIN/PORT.");
    }
#line 9406 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1051: // macro_pin_option: start_macro_port "END"
#line 5604 "LefParser.yy"
    {
      // Since in start_macro_port it has call the Init method, here
      // we need to call the Destroy method.
      // Still add a null pointer to set the number of port
      if (/*driver.lefrPinCbk*/ 1) {
        driver.lefrPin.lefiPin::addPort(driver.lefrGeometriesPtr);
        driver.lefrGeometriesPtr = 0;
        driver.lefrDoGeometries = 0;
      }
      driver.hasGeoLayer = 0;
    }
#line 9422 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1052: // macro_pin_option: "ANTENNASIZE" NUMBER opt_layer_name ';'
#line 5616 "LefParser.yy"
    {  /* a pre 5.4 syntax */
      driver.use5_3 = 1;
      if (driver.ignoreVersion) {
        /* do nothing */
      } else if (driver.versionNum >= 5.4) {
        if (driver.use5_4) {
           if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
             if (driver.pinWarnings++ < driver.lefrPinWarnings) {
                driver.outMsg = (char*)lefMalloc(10000);
                sprintf (driver.outMsg,
                   "ANTENNASIZE statement is a version 5.3 and earlier syntax.\nYour lef file is defined with version %g.", driver.versionNum);
                driver.lefError(1644, driver.outMsg);
                lefFree(driver.outMsg);
                /*CHKERR();*/
             }
           }
        }
      }
      if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addAntennaSize((yystack_[2].value.doubleVal), (*(yystack_[1].value.stringVal)).c_str());
      delete (yystack_[1].value.stringVal);
    }
#line 9448 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1053: // macro_pin_option: "ANTENNAMETALAREA" NUMBER opt_layer_name ';'
#line 5638 "LefParser.yy"
    {  /* a pre 5.4 syntax */
      driver.use5_3 = 1;
      if (driver.ignoreVersion) {
        /* do nothing */
      } else if (driver.versionNum >= 5.4) {
        if (driver.use5_4) {
           if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
              if (driver.pinWarnings++ < driver.lefrPinWarnings) {
                 driver.outMsg = (char*)lefMalloc(10000);
                 sprintf (driver.outMsg,
                    "ANTENNAMETALAREA statement is a version 5.3 and earlier syntax.\nYour lef file is defined with version %g.", driver.versionNum);
                 driver.lefError(1645, driver.outMsg);
                 lefFree(driver.outMsg);
                 /*CHKERR();*/
              }
           }
        }
      }
      if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addAntennaMetalArea((yystack_[2].value.doubleVal), (*(yystack_[1].value.stringVal)).c_str());
      delete (yystack_[1].value.stringVal);
    }
#line 9474 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1054: // macro_pin_option: "ANTENNAMETALLENGTH" NUMBER opt_layer_name ';'
#line 5660 "LefParser.yy"
    { /* a pre 5.4 syntax */ 
      driver.use5_3 = 1;
      if (driver.ignoreVersion) {
        /* do nothing */
      } else if (driver.versionNum >= 5.4) {
        if (driver.use5_4) {
           if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
              if (driver.pinWarnings++ < driver.lefrPinWarnings) {
                 driver.outMsg = (char*)lefMalloc(10000);
                 sprintf (driver.outMsg,
                    "ANTENNAMETALLENGTH statement is a version 5.3 and earlier syntax.\nYour lef file is defined with version %g.", driver.versionNum);
                 driver.lefError(1646, driver.outMsg);
                 lefFree(driver.outMsg);
                 /*CHKERR();*/
              }
           }
        }
      }
      if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addAntennaMetalLength((yystack_[2].value.doubleVal), (*(yystack_[1].value.stringVal)).c_str());
      delete (yystack_[1].value.stringVal);
    }
#line 9500 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1055: // macro_pin_option: "RISESLEWLIMIT" NUMBER ';'
#line 5682 "LefParser.yy"
    { if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setRiseSlewLimit((yystack_[1].value.doubleVal)); }
#line 9506 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1056: // macro_pin_option: "FALLSLEWLIMIT" NUMBER ';'
#line 5684 "LefParser.yy"
    { if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setFallSlewLimit((yystack_[1].value.doubleVal)); }
#line 9512 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1057: // macro_pin_option: "ANTENNAPARTIALMETALAREA" NUMBER opt_layer_name ';'
#line 5686 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (driver.ignoreVersion) {
        /* do nothing */
      } else if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAPARTIALMETALAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1647, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      } else if (driver.use5_3) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAPARTIALMETALAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1647, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
      if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addAntennaPartialMetalArea((yystack_[2].value.doubleVal), (*(yystack_[1].value.stringVal)).c_str());
      delete (yystack_[1].value.stringVal);
    }
#line 9547 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1058: // macro_pin_option: "ANTENNAPARTIALMETALSIDEAREA" NUMBER opt_layer_name ';'
#line 5717 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (driver.ignoreVersion) {
        /* do nothing */
      } else if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAPARTIALMETALSIDEAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1648, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      } else if (driver.use5_3) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAPARTIALMETALSIDEAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1648, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
      if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addAntennaPartialMetalSideArea((yystack_[2].value.doubleVal), (*(yystack_[1].value.stringVal)).c_str());
      delete (yystack_[1].value.stringVal);
    }
#line 9582 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1059: // macro_pin_option: "ANTENNAPARTIALCUTAREA" NUMBER opt_layer_name ';'
#line 5748 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (driver.ignoreVersion) {
        /* do nothing */
      } else if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAPARTIALCUTAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1649, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      } else if (driver.use5_3) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAPARTIALCUTAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1649, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
      if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addAntennaPartialCutArea((yystack_[2].value.doubleVal), (*(yystack_[1].value.stringVal)).c_str());
      delete (yystack_[1].value.stringVal);
    }
#line 9617 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1060: // macro_pin_option: "ANTENNADIFFAREA" NUMBER opt_layer_name ';'
#line 5779 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (driver.ignoreVersion) {
        /* do nothing */
      } else if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNADIFFAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1650, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      } else if (driver.use5_3) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNADIFFAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1650, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
      if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addAntennaDiffArea((yystack_[2].value.doubleVal), (*(yystack_[1].value.stringVal)).c_str());
      delete (yystack_[1].value.stringVal);
    }
#line 9652 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1061: // macro_pin_option: "ANTENNAGATEAREA" NUMBER opt_layer_name ';'
#line 5810 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (driver.ignoreVersion) {
        /* do nothing */
      } else if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAGATEAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1651, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      } else if (driver.use5_3) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAGATEAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1651, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
      if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addAntennaGateArea((yystack_[2].value.doubleVal), (*(yystack_[1].value.stringVal)).c_str());
      delete (yystack_[1].value.stringVal);
    }
#line 9687 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1062: // macro_pin_option: "ANTENNAMAXAREACAR" NUMBER req_layer_name ';'
#line 5841 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (driver.ignoreVersion) {
        /* do nothing */
      } else if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAMAXAREACAR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1652, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      } else if (driver.use5_3) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAMAXAREACAR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1652, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
      if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addAntennaMaxAreaCar((yystack_[2].value.doubleVal), (*(yystack_[1].value.stringVal)).c_str());
      delete (yystack_[1].value.stringVal);
    }
#line 9722 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1063: // macro_pin_option: "ANTENNAMAXSIDEAREACAR" NUMBER req_layer_name ';'
#line 5872 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (driver.ignoreVersion) {
        /* do nothing */
      } else if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAMAXSIDEAREACAR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1653, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      } else if (driver.use5_3) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAMAXSIDEAREACAR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1653, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
      if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addAntennaMaxSideAreaCar((yystack_[2].value.doubleVal), (*(yystack_[1].value.stringVal)).c_str());
      delete (yystack_[1].value.stringVal);
    }
#line 9757 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1064: // macro_pin_option: "ANTENNAMAXCUTCAR" NUMBER req_layer_name ';'
#line 5903 "LefParser.yy"
    { /* 5.4 syntax */
      driver.use5_4 = 1;
      if (driver.ignoreVersion) {
        /* do nothing */
      } else if (driver.versionNum < 5.4) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAMAXCUTCAR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1654, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      } else if (driver.use5_3) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAMAXCUTCAR statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1654, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
      if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::addAntennaMaxCutCar((yystack_[2].value.doubleVal), (*(yystack_[1].value.stringVal)).c_str());
      delete (yystack_[1].value.stringVal);
    }
#line 9792 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1065: // $@90: %empty
#line 5934 "LefParser.yy"
    { /* 5.5 syntax */
      driver.use5_4 = 1;
      if (driver.ignoreVersion) {
        /* do nothing */
      } else if (driver.versionNum < 5.5) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAMODEL statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1655, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      } else if (driver.use5_3) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ANTENNAMODEL statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1655, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      }
    }
#line 9825 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1067: // $@91: %empty
#line 5963 "LefParser.yy"
              {driver.lefDumbMode = 2; driver.lefNoNum = 2; }
#line 9831 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1068: // macro_pin_option: "NETEXPR" $@91 "qstring" ';'
#line 5964 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "NETEXPR statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1656, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      } else
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setNetExpr((*(yystack_[1].value.qstringVal)).c_str());
    delete (yystack_[1].value.qstringVal);
    }
#line 9852 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1069: // $@92: %empty
#line 5980 "LefParser.yy"
                        {driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 9858 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1070: // macro_pin_option: "SUPPLYSENSITIVITY" $@92 "string" ';'
#line 5981 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "SUPPLYSENSITIVITY statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1657, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      } else
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setSupplySensitivity((*(yystack_[1].value.stringVal)).c_str());
    delete (yystack_[1].value.stringVal);
    }
#line 9879 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1071: // $@93: %empty
#line 5997 "LefParser.yy"
                        {driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 9885 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1072: // macro_pin_option: "GROUNDSENSITIVITY" $@93 "string" ';'
#line 5998 "LefParser.yy"
    {
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrPinCbk*/ 1) { /* write error only if cbk is set */
           if (driver.pinWarnings++ < driver.lefrPinWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "GROUNDSENSITIVITY statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1658, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      } else
        if (/*driver.lefrPinCbk*/ 1) driver.lefrPin.lefiPin::setGroundSensitivity((*(yystack_[1].value.stringVal)).c_str());
    delete (yystack_[1].value.stringVal);
    }
#line 9906 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1073: // pin_layer_oxide: "OXIDE1"
#line 6017 "LefParser.yy"
    {
    if (/*driver.lefrPinCbk*/ 1)
       driver.lefrPin.lefiPin::addAntennaModel(1);
    }
#line 9915 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1074: // pin_layer_oxide: "OXIDE2"
#line 6022 "LefParser.yy"
    {
    if (/*driver.lefrPinCbk*/ 1)
       driver.lefrPin.lefiPin::addAntennaModel(2);
    }
#line 9924 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1075: // pin_layer_oxide: "OXIDE3"
#line 6027 "LefParser.yy"
    {
    if (/*driver.lefrPinCbk*/ 1)
       driver.lefrPin.lefiPin::addAntennaModel(3);
    }
#line 9933 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1076: // pin_layer_oxide: "OXIDE4"
#line 6032 "LefParser.yy"
    {
    if (/*driver.lefrPinCbk*/ 1)
       driver.lefrPin.lefiPin::addAntennaModel(4);
    }
#line 9942 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1079: // pin_name_value_pair: "string" NUMBER
#line 6044 "LefParser.yy"
    { 
      char temp[32];
      sprintf(temp, "%.11g", (yystack_[0].value.doubleVal));
      if (/*driver.lefrPinCbk*/ 1) {
         char propTp;
         propTp = driver.lefrPinProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrPin.lefiPin::setNumProperty((*(yystack_[1].value.stringVal)).c_str(), (yystack_[0].value.doubleVal), temp, propTp);
      }
      delete (yystack_[1].value.stringVal);
    }
#line 9957 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1080: // pin_name_value_pair: "string" "qstring"
#line 6055 "LefParser.yy"
    {
      if (/*driver.lefrPinCbk*/ 1) {
         char propTp;
         propTp = driver.lefrPinProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrPin.lefiPin::setProperty((*(yystack_[1].value.stringVal)).c_str(), (*(yystack_[0].value.qstringVal)).c_str(), propTp);
      }
      delete (yystack_[1].value.stringVal); delete (yystack_[0].value.qstringVal);
    }
#line 9970 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1081: // pin_name_value_pair: "string" "string"
#line 6064 "LefParser.yy"
    {
      if (/*driver.lefrPinCbk*/ 1) {
         char propTp;
         propTp = driver.lefrPinProp.lefiPropType::propType((*(yystack_[1].value.stringVal)).c_str());
         driver.lefrPin.lefiPin::setProperty((*(yystack_[1].value.stringVal)).c_str(), (*(yystack_[0].value.stringVal)).c_str(), propTp);
      }
      delete (yystack_[1].value.stringVal); delete (yystack_[0].value.stringVal);
    }
#line 9983 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1082: // electrical_direction: "DIRECTION" "INPUT" ';'
#line 6074 "LefParser.yy"
                                  {(yylhs.value.stringVal) = new std::string ("INPUT");}
#line 9989 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1083: // electrical_direction: "DIRECTION" "OUTPUT" ';'
#line 6075 "LefParser.yy"
                                        {(yylhs.value.stringVal) = new std::string ("OUTPUT");}
#line 9995 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1084: // electrical_direction: "DIRECTION" "OUTPUT" "TRISTATE" ';'
#line 6076 "LefParser.yy"
                                        {(yylhs.value.stringVal) = new std::string ("OUTPUT TRISTATE");}
#line 10001 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1085: // electrical_direction: "DIRECTION" "INOUT" ';'
#line 6077 "LefParser.yy"
                                        {(yylhs.value.stringVal) = new std::string ("INOUT");}
#line 10007 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1086: // electrical_direction: "DIRECTION" "FEEDTHRU" ';'
#line 6078 "LefParser.yy"
                                        {(yylhs.value.stringVal) = new std::string ("FEEDTHRU");}
#line 10013 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1087: // start_macro_port: "PORT"
#line 6081 "LefParser.yy"
    {
      if (/*driver.lefrPinCbk*/ 1) {
	driver.lefrDoGeometries = 1;
        driver.hasPRP = 0;
	driver.lefrGeometriesPtr = (lefiGeometries*)lefMalloc( sizeof(lefiGeometries));
	driver.lefrGeometriesPtr->lefiGeometries::Init();
      }
      driver.needGeometry = 0;  // don't need rect/path/poly define yet
      driver.hasGeoLayer = 0;   // make sure LAYER is set before geometry
    }
#line 10028 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1089: // macro_port_class_option: "CLASS" class_type ';'
#line 6094 "LefParser.yy"
    { if (driver.lefrDoGeometries)
        driver.lefrGeometriesPtr->lefiGeometries::addClass((*(yystack_[1].value.stringVal)).c_str()); 
    delete (yystack_[1].value.stringVal);
    }
#line 10037 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1090: // macro_pin_use: "SIGNAL"
#line 6100 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("SIGNAL");}
#line 10043 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1091: // macro_pin_use: "ANALOG"
#line 6101 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("ANALOG");}
#line 10049 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1092: // macro_pin_use: "POWER"
#line 6102 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("POWER");}
#line 10055 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1093: // macro_pin_use: "GROUND"
#line 6103 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("GROUND");}
#line 10061 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1094: // macro_pin_use: "CLOCK"
#line 6104 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("CLOCK");}
#line 10067 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1095: // macro_pin_use: "DATA"
#line 6105 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("DATA");}
#line 10073 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1096: // macro_scan_use: "INPUT"
#line 6108 "LefParser.yy"
          {(yylhs.value.stringVal) = new std::string ("INPUT");}
#line 10079 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1097: // macro_scan_use: "OUTPUT"
#line 6109 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("OUTPUT");}
#line 10085 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1098: // macro_scan_use: "START"
#line 6110 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("START");}
#line 10091 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1099: // macro_scan_use: "STOP"
#line 6111 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("STOP");}
#line 10097 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1100: // pin_shape: %empty
#line 6114 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string (""); }
#line 10103 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1101: // pin_shape: "ABUTMENT"
#line 6115 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("ABUTMENT");}
#line 10109 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1102: // pin_shape: "RING"
#line 6116 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("RING");}
#line 10115 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1103: // pin_shape: "FEEDTHRU"
#line 6117 "LefParser.yy"
                {(yylhs.value.stringVal) = new std::string ("FEEDTHRU");}
#line 10121 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1105: // $@94: %empty
#line 6122 "LefParser.yy"
          {driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 10127 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1106: // $@95: %empty
#line 6123 "LefParser.yy"
    {
      if ((driver.needGeometry) && (driver.needGeometry != 2)) // 1 LAYER follow after another
        if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
          /* geometries is called by MACRO/OBS & MACRO/PIN/PORT */
          if (driver.obsDef)
             driver.lefWarning(2076, "Either PATH, RECT or POLYGON statement is a required in MACRO/OBS.");
          else
             driver.lefWarning(2065, "Either PATH, RECT or POLYGON statement is a required in MACRO/PIN/PORT.");
        }
      if (driver.lefrDoGeometries)
        driver.lefrGeometriesPtr->lefiGeometries::addLayer((*(yystack_[0].value.stringVal)).c_str());
      driver.needGeometry = 1;    // within LAYER it requires either path, rect, poly
      driver.hasGeoLayer = 1;
      delete (yystack_[0].value.stringVal);
    }
#line 10147 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1108: // geometry: "WIDTH" NUMBER ';'
#line 6141 "LefParser.yy"
    { 
      if (driver.lefrDoGeometries) {
        if (driver.hasGeoLayer == 0) {   /* LAYER statement is missing */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              /*CHKERR();*/
           }
        } else
           driver.lefrGeometriesPtr->lefiGeometries::addWidth((yystack_[1].value.doubleVal)); 
      } 
    }
#line 10163 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1109: // geometry: "PATH" firstPt otherPts ';'
#line 6153 "LefParser.yy"
    { if (driver.lefrDoGeometries) {
        if (driver.hasGeoLayer == 0) {   /* LAYER statement is missing */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              /*CHKERR();*/
           }
        } else
           driver.lefrGeometriesPtr->lefiGeometries::addPath();
      }
      driver.hasPRP = 1;
      driver.needGeometry = 2;
    }
#line 10180 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1110: // geometry: "PATH" "ITERATE" firstPt otherPts stepPattern ';'
#line 6166 "LefParser.yy"
    { if (driver.lefrDoGeometries) {
        if (driver.hasGeoLayer == 0) {   /* LAYER statement is missing */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              /*CHKERR();*/
           }
        } else
           driver.lefrGeometriesPtr->lefiGeometries::addPathIter();
      } 
      driver.hasPRP = 1;
      driver.needGeometry = 2;
    }
#line 10197 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1111: // geometry: "RECT" pt pt ';'
#line 6179 "LefParser.yy"
    { if (driver.lefrDoGeometries) {
        if (driver.hasGeoLayer == 0) {   /* LAYER statement is missing */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              /*CHKERR();*/
           }
        } else
           driver.lefrGeometriesPtr->lefiGeometries::addRect((yystack_[2].value.pt)->x, (yystack_[2].value.pt)->y, (yystack_[1].value.pt)->x, (yystack_[1].value.pt)->y);
      }
      driver.needGeometry = 2;
      delete (yystack_[2].value.pt); delete (yystack_[1].value.pt);
    }
#line 10214 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1112: // geometry: "RECT" "ITERATE" pt pt stepPattern ';'
#line 6192 "LefParser.yy"
    { if (driver.lefrDoGeometries) {
        if (driver.hasGeoLayer == 0) {   /* LAYER statement is missing */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              /*CHKERR();*/
           }
        } else
           driver.lefrGeometriesPtr->lefiGeometries::addRectIter((yystack_[3].value.pt)->x, (yystack_[3].value.pt)->y, (yystack_[2].value.pt)->x,
                                                          (yystack_[2].value.pt)->y);
      }
      driver.needGeometry = 2;
      delete (yystack_[3].value.pt); delete (yystack_[2].value.pt);
    }
#line 10232 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1113: // geometry: "POLYGON" firstPt nextPt nextPt otherPts ';'
#line 6206 "LefParser.yy"
    {
      if (driver.lefrDoGeometries) {
        if (driver.hasGeoLayer == 0) {   /* LAYER statement is missing */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              /*CHKERR();*/
           }
        } else
           driver.lefrGeometriesPtr->lefiGeometries::addPolygon();
      }
      driver.hasPRP = 1;
      driver.needGeometry = 2;
    }
#line 10250 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1114: // geometry: "POLYGON" "ITERATE" firstPt nextPt nextPt otherPts stepPattern ';'
#line 6220 "LefParser.yy"
    { if (driver.lefrDoGeometries) {
        if (driver.hasGeoLayer == 0) {   /* LAYER statement is missing */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.lefError(1701, "A LAYER statement is missing in Geometry.\nLAYER is a required statement before any geometry can be defined.");
              /*CHKERR();*/
           }
        } else
           driver.lefrGeometriesPtr->lefiGeometries::addPolygonIter();
      }
      driver.hasPRP = 1;
      driver.needGeometry = 2;
    }
#line 10267 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1115: // geometry: via_placement
#line 6233 "LefParser.yy"
    { }
#line 10273 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1119: // layer_exceptpgnet: "EXCEPTPGNET"
#line 6240 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
        driver.outMsg = (char*)lefMalloc(10000);
        sprintf(driver.outMsg,
          "EXCEPTPGNET is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
        driver.lefError(1699, driver.outMsg);
        lefFree(driver.outMsg);
        /*CHKERR();*/
      } else {
       if (driver.lefrDoGeometries)
        driver.lefrGeometriesPtr->lefiGeometries::addLayerExceptPgNet();
      }
    }
#line 10291 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1121: // layer_spacing: "SPACING" NUMBER
#line 6256 "LefParser.yy"
    { if (driver.lefrDoGeometries) {
        if (driver.zeroOrGt((yystack_[0].value.doubleVal)))
           driver.lefrGeometriesPtr->lefiGeometries::addLayerMinSpacing((yystack_[0].value.doubleVal));
        else {
           driver.outMsg = (char*)lefMalloc(10000);
           sprintf (driver.outMsg,
              "THE SPACING statement has the value %g in MACRO OBS.\nValue has to be 0 or greater.", (yystack_[0].value.doubleVal));
           driver.lefError(1659, driver.outMsg);
           lefFree(driver.outMsg);
           /*CHKERR();*/
        }
      }
    }
#line 10309 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1122: // layer_spacing: "DESIGNRULEWIDTH" NUMBER
#line 6270 "LefParser.yy"
    { if (driver.lefrDoGeometries) {
        if (driver.zeroOrGt((yystack_[0].value.doubleVal)))
           driver.lefrGeometriesPtr->lefiGeometries::addLayerRuleWidth((yystack_[0].value.doubleVal));
        else {
           driver.outMsg = (char*)lefMalloc(10000);
           sprintf (driver.outMsg,
              "THE DESIGNRULEWIDTH statement has the value %g in MACRO OBS.\nValue has to be 0 or greater.", (yystack_[0].value.doubleVal));
           driver.lefError(1660, driver.outMsg);
           lefFree(driver.outMsg);
           /*CHKERR();*/
        }
      }
    }
#line 10327 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1123: // firstPt: pt
#line 6285 "LefParser.yy"
    { if (driver.lefrDoGeometries)
        driver.lefrGeometriesPtr->lefiGeometries::startList((yystack_[0].value.pt)->x, (yystack_[0].value.pt)->y); 
    delete (yystack_[0].value.pt);
    }
#line 10336 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1124: // nextPt: pt
#line 6291 "LefParser.yy"
    { if (driver.lefrDoGeometries)
        driver.lefrGeometriesPtr->lefiGeometries::addToList((yystack_[0].value.pt)->x, (yystack_[0].value.pt)->y); 
    delete (yystack_[0].value.pt);
    }
#line 10345 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1127: // $@96: %empty
#line 6305 "LefParser.yy"
           {driver.lefDumbMode = 1;}
#line 10351 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1128: // via_placement: "VIA" pt $@96 "string" ';'
#line 6306 "LefParser.yy"
    { if (driver.lefrDoGeometries)
        driver.lefrGeometriesPtr->lefiGeometries::addVia((yystack_[3].value.pt)->x, (yystack_[3].value.pt)->y, (*(yystack_[1].value.stringVal)).c_str()); 
    delete (yystack_[3].value.pt); delete (yystack_[1].value.stringVal);
    }
#line 10360 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1129: // $@97: %empty
#line 6310 "LefParser.yy"
                       {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 10366 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1130: // via_placement: "VIA" "ITERATE" pt $@97 "string" stepPattern ';'
#line 6312 "LefParser.yy"
    { if (driver.lefrDoGeometries)
        driver.lefrGeometriesPtr->lefiGeometries::addViaIter((yystack_[4].value.pt)->x, (yystack_[4].value.pt)->y, (*(yystack_[2].value.stringVal)).c_str()); 
    delete (yystack_[4].value.pt); delete (yystack_[2].value.stringVal);
    }
#line 10375 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1131: // stepPattern: "DO" NUMBER "BY" NUMBER "STEP" NUMBER NUMBER
#line 6319 "LefParser.yy"
     { if (driver.lefrDoGeometries)
         driver.lefrGeometriesPtr->lefiGeometries::addStepPattern((yystack_[5].value.doubleVal), (yystack_[3].value.doubleVal), (yystack_[1].value.doubleVal), (yystack_[0].value.doubleVal)); }
#line 10382 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1132: // sitePattern: GSTRING NUMBER NUMBER orientation "DO" NUMBER "BY" NUMBER "STEP" NUMBER NUMBER
#line 6324 "LefParser.yy"
    {
      if (driver.lefrDoSite) {
	driver.lefrSitePatternPtr = (lefiSitePattern*)lefMalloc(
				   sizeof(lefiSitePattern));
	driver.lefrSitePatternPtr->lefiSitePattern::Init();
	driver.lefrSitePatternPtr->lefiSitePattern::set((*(yystack_[10].value.stringVal)).c_str(), (yystack_[9].value.doubleVal), (yystack_[8].value.doubleVal), (yystack_[7].value.integerVal), (yystack_[5].value.doubleVal), (yystack_[3].value.doubleVal),
	  (yystack_[1].value.doubleVal), (yystack_[0].value.doubleVal));
	}
    delete (yystack_[10].value.stringVal);
    }
#line 10397 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1133: // sitePattern: GSTRING NUMBER NUMBER orientation
#line 6335 "LefParser.yy"
    {
      if (driver.lefrDoSite) {
	driver.lefrSitePatternPtr = (lefiSitePattern*)lefMalloc(
				   sizeof(lefiSitePattern));
	driver.lefrSitePatternPtr->lefiSitePattern::Init();
	driver.lefrSitePatternPtr->lefiSitePattern::set((*(yystack_[3].value.stringVal)).c_str(), (yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal), (yystack_[0].value.integerVal), -1, -1,
	  -1, -1);
	}
    delete (yystack_[3].value.stringVal);
    }
#line 10412 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1134: // $@98: %empty
#line 6348 "LefParser.yy"
    { 
      if (driver.lefrDoTrack) {
	driver.lefrTrackPatternPtr = (lefiTrackPattern*)lefMalloc(
				sizeof(lefiTrackPattern));
	driver.lefrTrackPatternPtr->lefiTrackPattern::Init();
	driver.lefrTrackPatternPtr->lefiTrackPattern::set("X", (yystack_[4].value.doubleVal), (int)(yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
      }    
    }
#line 10425 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1135: // $@99: %empty
#line 6356 "LefParser.yy"
            {driver.lefDumbMode = 1000000000;}
#line 10431 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1136: // trackPattern: "X" NUMBER "DO" NUMBER "STEP" NUMBER $@98 "LAYER" $@99 trackLayers
#line 6357 "LefParser.yy"
    { driver.lefDumbMode = 0;}
#line 10437 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1137: // $@100: %empty
#line 6359 "LefParser.yy"
    { 
      if (driver.lefrDoTrack) {
	driver.lefrTrackPatternPtr = (lefiTrackPattern*)lefMalloc(
                                    sizeof(lefiTrackPattern));
	driver.lefrTrackPatternPtr->lefiTrackPattern::Init();
	driver.lefrTrackPatternPtr->lefiTrackPattern::set("Y", (yystack_[4].value.doubleVal), (int)(yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
      }    
    }
#line 10450 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1138: // $@101: %empty
#line 6367 "LefParser.yy"
            {driver.lefDumbMode = 1000000000;}
#line 10456 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1139: // trackPattern: "Y" NUMBER "DO" NUMBER "STEP" NUMBER $@100 "LAYER" $@101 trackLayers
#line 6368 "LefParser.yy"
    { driver.lefDumbMode = 0;}
#line 10462 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1140: // trackPattern: "X" NUMBER "DO" NUMBER "STEP" NUMBER
#line 6370 "LefParser.yy"
    { 
      if (driver.lefrDoTrack) {
	driver.lefrTrackPatternPtr = (lefiTrackPattern*)lefMalloc(
                                    sizeof(lefiTrackPattern));
	driver.lefrTrackPatternPtr->lefiTrackPattern::Init();
	driver.lefrTrackPatternPtr->lefiTrackPattern::set("X", (yystack_[4].value.doubleVal), (int)(yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
      }    
    }
#line 10475 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1141: // trackPattern: "Y" NUMBER "DO" NUMBER "STEP" NUMBER
#line 6379 "LefParser.yy"
    { 
      if (driver.lefrDoTrack) {
	driver.lefrTrackPatternPtr = (lefiTrackPattern*)lefMalloc(
                                    sizeof(lefiTrackPattern));
	driver.lefrTrackPatternPtr->lefiTrackPattern::Init();
	driver.lefrTrackPatternPtr->lefiTrackPattern::set("Y", (yystack_[4].value.doubleVal), (int)(yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
      }    
    }
#line 10488 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1144: // layer_name: "string"
#line 6394 "LefParser.yy"
    { if (driver.lefrDoTrack) driver.lefrTrackPatternPtr->lefiTrackPattern::addLayer((*(yystack_[0].value.stringVal)).c_str()); delete (yystack_[0].value.stringVal);}
#line 10494 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1145: // gcellPattern: "X" NUMBER "DO" NUMBER "STEP" NUMBER
#line 6397 "LefParser.yy"
    {
      if (driver.lefrDoGcell) {
	driver.lefrGcellPatternPtr = (lefiGcellPattern*)lefMalloc(
                                    sizeof(lefiGcellPattern));
	driver.lefrGcellPatternPtr->lefiGcellPattern::Init();
	driver.lefrGcellPatternPtr->lefiGcellPattern::set("X", (yystack_[4].value.doubleVal), (int)(yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
      }    
    }
#line 10507 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1146: // gcellPattern: "Y" NUMBER "DO" NUMBER "STEP" NUMBER
#line 6406 "LefParser.yy"
    {
      if (driver.lefrDoGcell) {
	driver.lefrGcellPatternPtr = (lefiGcellPattern*)lefMalloc(
                                    sizeof(lefiGcellPattern));
	driver.lefrGcellPatternPtr->lefiGcellPattern::Init();
	driver.lefrGcellPatternPtr->lefiGcellPattern::set("Y", (yystack_[4].value.doubleVal), (int)(yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
      }    
    }
#line 10520 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1147: // macro_obs: start_macro_obs geometries "END"
#line 6416 "LefParser.yy"
    { 
      if (/*driver.lefrObstructionCbk*/ 1) {
	driver.lefrObstruction.lefiObstruction::setGeometries(driver.lefrGeometriesPtr);
	driver.lefrGeometriesPtr = 0;
	driver.lefrDoGeometries = 0;
        driver.lefrObstructionCbk( driver.lefrObstruction);
      }
      driver.lefDumbMode = 0;
      driver.hasGeoLayer = 0;       /* reset */
    }
#line 10535 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1148: // macro_obs: start_macro_obs "END"
#line 6434 "LefParser.yy"
    {
       /* The pointer has malloced in start, need to free manually */
       if (driver.lefrGeometriesPtr) {
          driver.lefrGeometriesPtr->lefiGeometries::Destroy();
          lefFree(driver.lefrGeometriesPtr);
          driver.lefrGeometriesPtr = 0;
	  driver.lefrDoGeometries = 0;
       }
       driver.hasGeoLayer = 0;
    }
#line 10550 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1149: // start_macro_obs: "OBS"
#line 6446 "LefParser.yy"
    {
      driver.obsDef = 1;
      if (/*driver.lefrObstructionCbk*/ 1) {
	driver.lefrDoGeometries = 1;
	driver.lefrGeometriesPtr = (lefiGeometries*)lefMalloc(
	    sizeof(lefiGeometries));
	driver.lefrGeometriesPtr->lefiGeometries::Init();
	}
      driver.hasGeoLayer = 0;
    }
#line 10565 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1150: // macro_density: "DENSITY" density_layer density_layers "END"
#line 6458 "LefParser.yy"
    { 
      if (driver.versionNum < 5.6) {
        if (/*driver.lefrDensityCbk*/ 1) { /* write error only if cbk is set */
           if (driver.macroWarnings++ < driver.lefrMacroWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "DENSITY statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1661, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
      } 
      if (/*driver.lefrDensityCbk*/ 1) {
        driver.lefrDensityCbk( driver.lefrDensity);
        driver.lefrDensity.lefiDensity::clear();
      }
      driver.lefDumbMode = 0;
    }
#line 10589 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1153: // $@102: %empty
#line 6482 "LefParser.yy"
                       { driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 10595 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1154: // $@103: %empty
#line 6483 "LefParser.yy"
    {
      if (/*driver.lefrDensityCbk*/ 1)
        driver.lefrDensity.lefiDensity::addLayer((*(yystack_[1].value.stringVal)).c_str());
    delete (yystack_[1].value.stringVal);
    }
#line 10605 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1158: // density_layer_rect: "RECT" pt pt NUMBER ';'
#line 6495 "LefParser.yy"
    {
      if (/*driver.lefrDensityCbk*/ 1)
        driver.lefrDensity.lefiDensity::addRect((yystack_[3].value.pt)->x, (yystack_[3].value.pt)->y, (yystack_[2].value.pt)->x, (yystack_[2].value.pt)->y, (yystack_[1].value.doubleVal)); 
    delete (yystack_[3].value.pt); delete (yystack_[2].value.pt);
    }
#line 10615 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1159: // $@104: %empty
#line 6501 "LefParser.yy"
                             { driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 10621 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1160: // macro_clocktype: "CLOCKTYPE" $@104 "string" ';'
#line 6502 "LefParser.yy"
    { if (/*driver.lefrMacroCbk*/ 1) driver.lefrMacro.lefiMacro::setClockType((*(yystack_[1].value.stringVal)).c_str()); delete (yystack_[1].value.stringVal);}
#line 10627 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1161: // timing: start_timing timing_options end_timing
#line 6505 "LefParser.yy"
    { }
#line 10633 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1162: // start_timing: "TIMING"
#line 6508 "LefParser.yy"
    { /* XXXXX for macros */ }
#line 10639 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1163: // end_timing: "END" "TIMING"
#line 6511 "LefParser.yy"
  {
    if (driver.versionNum < 5.4) {
      if (/*driver.lefrTimingCbk*/ 1 && driver.lefrTiming.lefiTiming::hasData())
        driver.lefrTimingCbk( driver.lefrTiming);
      driver.lefrTiming.lefiTiming::clear();
    } else {
      if (/*driver.lefrTimingCbk*/ 1) /* write warning only if cbk is set */
        if (driver.timingWarnings++ < driver.lefrTimingWarnings)
          driver.lefWarning(2066, "MACRO TIMING statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
      driver.lefrTiming.lefiTiming::clear();
    }
  }
#line 10656 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1166: // $@105: %empty
#line 6531 "LefParser.yy"
    {
    if (driver.versionNum < 5.4) {
      if (/*driver.lefrTimingCbk*/ 1 && driver.lefrTiming.lefiTiming::hasData())
        driver.lefrTimingCbk( driver.lefrTiming);
    }
    driver.lefDumbMode = 1000000000;
    driver.lefrTiming.lefiTiming::clear();
    }
#line 10669 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1167: // timing_option: "FROMPIN" $@105 list_of_from_strings ';'
#line 6540 "LefParser.yy"
    { driver.lefDumbMode = 0;}
#line 10675 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1168: // $@106: %empty
#line 6541 "LefParser.yy"
            {driver.lefDumbMode = 1000000000;}
#line 10681 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1169: // timing_option: "TOPIN" $@106 list_of_to_strings ';'
#line 6542 "LefParser.yy"
    { driver.lefDumbMode = 0;}
#line 10687 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1170: // $@107: %empty
#line 6544 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addRiseFall((*(yystack_[3].value.stringVal)).c_str(),(yystack_[1].value.doubleVal),(yystack_[0].value.doubleVal)); delete (yystack_[3].value.stringVal);}
#line 10693 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1171: // timing_option: risefall "INTRINSIC" NUMBER NUMBER $@107 slew_spec "VARIABLE" NUMBER NUMBER ';'
#line 6546 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addRiseFallVariable((yystack_[2].value.doubleVal),(yystack_[1].value.doubleVal)); }
#line 10699 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1172: // timing_option: risefall delay_or_transition "UNATENESS" unateness "TABLEDIMENSION" NUMBER NUMBER NUMBER ';'
#line 6549 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) {
	if ((*(yystack_[7].value.stringVal))[0] == 'D' || (*(yystack_[7].value.stringVal))[0] == 'd') /* delay */
	  driver.lefrTiming.lefiTiming::addDelay((*(yystack_[8].value.stringVal)).c_str(), (*(yystack_[5].value.stringVal)).c_str(), (yystack_[3].value.doubleVal), (yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
	else
	  driver.lefrTiming.lefiTiming::addTransition((*(yystack_[8].value.stringVal)).c_str(), (*(yystack_[5].value.stringVal)).c_str(), (yystack_[3].value.doubleVal), (yystack_[2].value.doubleVal), (yystack_[1].value.doubleVal));
      }
      delete (yystack_[8].value.stringVal); delete (yystack_[7].value.stringVal); delete (yystack_[5].value.stringVal);
    }
#line 10712 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1173: // timing_option: "TABLEAXIS" list_of_table_axis_numbers ';'
#line 6558 "LefParser.yy"
    { }
#line 10718 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1174: // timing_option: "TABLEENTRIES" list_of_table_entries ';'
#line 6560 "LefParser.yy"
    { }
#line 10724 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1175: // timing_option: "RISERS" NUMBER NUMBER ';'
#line 6562 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::setRiseRS((yystack_[2].value.doubleVal),(yystack_[1].value.doubleVal)); }
#line 10730 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1176: // timing_option: "FALLRS" NUMBER NUMBER ';'
#line 6564 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::setFallRS((yystack_[2].value.doubleVal),(yystack_[1].value.doubleVal)); }
#line 10736 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1177: // timing_option: "RISECS" NUMBER NUMBER ';'
#line 6566 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::setRiseCS((yystack_[2].value.doubleVal),(yystack_[1].value.doubleVal)); }
#line 10742 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1178: // timing_option: "FALLCS" NUMBER NUMBER ';'
#line 6568 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::setFallCS((yystack_[2].value.doubleVal),(yystack_[1].value.doubleVal)); }
#line 10748 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1179: // timing_option: "RISESATT1" NUMBER NUMBER ';'
#line 6570 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::setRiseAtt1((yystack_[2].value.doubleVal),(yystack_[1].value.doubleVal)); }
#line 10754 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1180: // timing_option: "FALLSATT1" NUMBER NUMBER ';'
#line 6572 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::setFallAtt1((yystack_[2].value.doubleVal),(yystack_[1].value.doubleVal)); }
#line 10760 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1181: // timing_option: "RISET0" NUMBER NUMBER ';'
#line 6574 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::setRiseTo((yystack_[2].value.doubleVal),(yystack_[1].value.doubleVal)); }
#line 10766 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1182: // timing_option: "FALLT0" NUMBER NUMBER ';'
#line 6576 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::setFallTo((yystack_[2].value.doubleVal),(yystack_[1].value.doubleVal)); }
#line 10772 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1183: // timing_option: "UNATENESS" unateness ';'
#line 6578 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addUnateness((*(yystack_[1].value.stringVal)).c_str()); delete (yystack_[1].value.stringVal);}
#line 10778 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1184: // timing_option: "STABLE" "SETUP" NUMBER "HOLD" NUMBER risefall ';'
#line 6580 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::setStable((yystack_[4].value.doubleVal),(yystack_[2].value.doubleVal),(*(yystack_[1].value.stringVal)).c_str()); delete (yystack_[1].value.stringVal);}
#line 10784 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1185: // timing_option: two_pin_trigger from_pin_trigger to_pin_trigger "TABLEDIMENSION" NUMBER NUMBER NUMBER ';'
#line 6582 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addSDF2Pins((*(yystack_[7].value.stringVal)).c_str(),(*(yystack_[6].value.stringVal)).c_str(),(*(yystack_[5].value.stringVal)).c_str(),(yystack_[3].value.doubleVal),(yystack_[2].value.doubleVal),(yystack_[1].value.doubleVal)); delete (yystack_[7].value.stringVal); delete (yystack_[6].value.stringVal); delete (yystack_[5].value.stringVal);}
#line 10790 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1186: // timing_option: one_pin_trigger "TABLEDIMENSION" NUMBER NUMBER NUMBER ';'
#line 6584 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addSDF1Pin((*(yystack_[5].value.stringVal)).c_str(),(yystack_[3].value.doubleVal),(yystack_[2].value.doubleVal),(yystack_[2].value.doubleVal)); delete (yystack_[5].value.stringVal);}
#line 10796 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1187: // timing_option: "SDFCONDSTART" "qstring" ';'
#line 6586 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::setSDFcondStart((*(yystack_[1].value.qstringVal)).c_str()); delete (yystack_[1].value.qstringVal);}
#line 10802 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1188: // timing_option: "SDFCONDEND" "qstring" ';'
#line 6588 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::setSDFcondEnd((*(yystack_[1].value.qstringVal)).c_str()); delete (yystack_[1].value.qstringVal);}
#line 10808 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1189: // timing_option: "SDFCOND" "qstring" ';'
#line 6590 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::setSDFcond((*(yystack_[1].value.qstringVal)).c_str()); delete (yystack_[1].value.qstringVal);}
#line 10814 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1190: // timing_option: "EXTENSION" ';'
#line 6592 "LefParser.yy"
    { /* XXXXX */ }
#line 10820 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1191: // one_pin_trigger: "MPWH"
#line 6596 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("MPWH");}
#line 10826 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1192: // one_pin_trigger: "MPWL"
#line 6598 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("MPWL");}
#line 10832 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1193: // one_pin_trigger: "PERIOD"
#line 6600 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("PERIOD");}
#line 10838 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1194: // two_pin_trigger: "SETUP"
#line 6604 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("SETUP");}
#line 10844 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1195: // two_pin_trigger: "HOLD"
#line 6606 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("HOLD");}
#line 10850 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1196: // two_pin_trigger: "RECOVERY"
#line 6608 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("RECOVERY");}
#line 10856 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1197: // two_pin_trigger: "SKEW"
#line 6610 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("SKEW");}
#line 10862 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1198: // from_pin_trigger: "ANYEDGE"
#line 6614 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("ANYEDGE");}
#line 10868 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1199: // from_pin_trigger: "POSEDGE"
#line 6616 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("POSEDGE");}
#line 10874 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1200: // from_pin_trigger: "NEGEDGE"
#line 6618 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("NEGEDGE");}
#line 10880 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1201: // to_pin_trigger: "ANYEDGE"
#line 6622 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("ANYEDGE");}
#line 10886 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1202: // to_pin_trigger: "POSEDGE"
#line 6624 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("POSEDGE");}
#line 10892 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1203: // to_pin_trigger: "NEGEDGE"
#line 6626 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("NEGEDGE");}
#line 10898 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1204: // delay_or_transition: "DELAY"
#line 6630 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("DELAY"); }
#line 10904 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1205: // delay_or_transition: "TRANSITIONTIME"
#line 6632 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("TRANSITION"); }
#line 10910 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1206: // list_of_table_entries: table_entry
#line 6636 "LefParser.yy"
    { }
#line 10916 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1207: // list_of_table_entries: list_of_table_entries table_entry
#line 6638 "LefParser.yy"
    { }
#line 10922 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1208: // table_entry: '(' NUMBER NUMBER NUMBER ')'
#line 6641 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addTableEntry((yystack_[3].value.doubleVal),(yystack_[2].value.doubleVal),(yystack_[1].value.doubleVal)); }
#line 10928 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1209: // list_of_table_axis_numbers: NUMBER
#line 6645 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addTableAxisNumber((yystack_[0].value.doubleVal)); }
#line 10934 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1210: // list_of_table_axis_numbers: list_of_table_axis_numbers NUMBER
#line 6647 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addTableAxisNumber((yystack_[0].value.doubleVal)); }
#line 10940 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1211: // slew_spec: %empty
#line 6651 "LefParser.yy"
    { }
#line 10946 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1212: // slew_spec: NUMBER NUMBER NUMBER NUMBER
#line 6653 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addRiseFallSlew((yystack_[3].value.doubleVal),(yystack_[2].value.doubleVal),(yystack_[1].value.doubleVal),(yystack_[0].value.doubleVal)); }
#line 10952 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1213: // slew_spec: NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER
#line 6655 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addRiseFallSlew((yystack_[6].value.doubleVal),(yystack_[5].value.doubleVal),(yystack_[4].value.doubleVal),(yystack_[3].value.doubleVal));
      if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addRiseFallSlew2((yystack_[2].value.doubleVal),(yystack_[1].value.doubleVal),(yystack_[0].value.doubleVal)); }
#line 10959 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1214: // risefall: "RISE"
#line 6660 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("RISE"); }
#line 10965 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1215: // risefall: "FALL"
#line 6662 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("FALL"); }
#line 10971 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1216: // unateness: "INVERT"
#line 6666 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("INVERT"); }
#line 10977 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1217: // unateness: "NONINVERT"
#line 6668 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("NONINVERT"); }
#line 10983 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1218: // unateness: "NONUNATE"
#line 6670 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string ("NONUNATE"); }
#line 10989 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1219: // list_of_from_strings: "string"
#line 6674 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addFromPin((*(yystack_[0].value.stringVal)).c_str()); delete (yystack_[0].value.stringVal);}
#line 10995 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1220: // list_of_from_strings: list_of_from_strings "string"
#line 6676 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addFromPin((*(yystack_[0].value.stringVal)).c_str()); delete (yystack_[0].value.stringVal);}
#line 11001 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1221: // list_of_to_strings: "string"
#line 6680 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addToPin((*(yystack_[0].value.stringVal)).c_str()); delete (yystack_[0].value.stringVal);}
#line 11007 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1222: // list_of_to_strings: list_of_to_strings "string"
#line 6682 "LefParser.yy"
    { if (/*driver.lefrTimingCbk*/ 1) driver.lefrTiming.lefiTiming::addToPin((*(yystack_[0].value.stringVal)).c_str()); delete (yystack_[0].value.stringVal);}
#line 11013 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1223: // $@108: %empty
#line 6685 "LefParser.yy"
    {
      if (/*driver.lefrArrayCbk*/ 1) driver.lefrArrayCbk( driver.lefrArray);
      driver.lefrArray.lefiArray::clear();
      driver.lefrSitePatternPtr = 0;
      driver.lefrDoSite = 0;
   }
#line 11024 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1225: // $@109: %empty
#line 6693 "LefParser.yy"
                     {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 11030 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1226: // start_array: "ARRAY" $@109 "string"
#line 6694 "LefParser.yy"
    {
      if (/*driver.lefrArrayCbk*/ 1) {
	driver.lefrArray.lefiArray::setName((*(yystack_[0].value.stringVal)).c_str());
	driver.lefrArrayBeginCbk( (*(yystack_[0].value.stringVal)).c_str());
      }
      //strcpy(driver.arrayName, $3);
      driver.arrayName = (*(yystack_[0].value.stringVal));
      delete (yystack_[0].value.stringVal);
    }
#line 11044 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1227: // $@110: %empty
#line 6704 "LefParser.yy"
                 {driver.lefDumbMode = 1; driver.lefNoNum = 1;}
#line 11050 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1228: // end_array: "END" $@110 "string"
#line 6705 "LefParser.yy"
    {
      if (/*driver.lefrArrayCbk && driver.lefrArrayEndCbk*/ 1)
	  driver.lefrArrayEndCbk( (*(yystack_[0].value.stringVal)).c_str());
      if (driver.arrayName != *(yystack_[0].value.stringVal)) {
        if (/*driver.lefrArrayCbk*/ 1) { /* write error only if cbk is set */
           if (driver.arrayWarnings++ < driver.lefrArrayWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "END ARRAY name %s is different from the ARRAY name %s.\nCorrect the LEF file before rerun it through the LEF parser.", (*(yystack_[0].value.stringVal)).c_str(), driver.arrayName.c_str());
              driver.lefError(1662, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           } 
        } 
      } 
      delete (yystack_[0].value.stringVal);
    }
#line 11072 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1229: // array_rules: %empty
#line 6725 "LefParser.yy"
    { }
#line 11078 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1230: // array_rules: array_rules array_rule
#line 6727 "LefParser.yy"
    { }
#line 11084 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1231: // $@111: %empty
#line 6730 "LefParser.yy"
            { if (/*driver.lefrArrayCbk*/ 1) driver.lefrDoSite = 1; driver.lefDumbMode = 1; }
#line 11090 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1232: // array_rule: site_word $@111 sitePattern ';'
#line 6732 "LefParser.yy"
    {
      if (/*driver.lefrArrayCbk*/ 1) {
	driver.lefrArray.lefiArray::addSitePattern(driver.lefrSitePatternPtr);
      }
    }
#line 11100 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1233: // $@112: %empty
#line 6737 "LefParser.yy"
               {driver.lefDumbMode = 1; if (/*driver.lefrArrayCbk*/ 1) driver.lefrDoSite = 1; }
#line 11106 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1234: // array_rule: "CANPLACE" $@112 sitePattern ';'
#line 6739 "LefParser.yy"
    {
      if (/*driver.lefrArrayCbk*/ 1) {
	driver.lefrArray.lefiArray::addCanPlace(driver.lefrSitePatternPtr);
      }
    }
#line 11116 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1235: // $@113: %empty
#line 6744 "LefParser.yy"
                   {driver.lefDumbMode = 1; if (/*driver.lefrArrayCbk*/ 1) driver.lefrDoSite = 1; }
#line 11122 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1236: // array_rule: "CANNOTOCCUPY" $@113 sitePattern ';'
#line 6746 "LefParser.yy"
    {
      if (/*driver.lefrArrayCbk*/ 1) {
	driver.lefrArray.lefiArray::addCannotOccupy(driver.lefrSitePatternPtr);
      }
    }
#line 11132 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1237: // $@114: %empty
#line 6751 "LefParser.yy"
             { if (/*driver.lefrArrayCbk*/ 1) driver.lefrDoTrack = 1; }
#line 11138 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1238: // array_rule: "TRACKS" $@114 trackPattern ';'
#line 6752 "LefParser.yy"
    {
      if (/*driver.lefrArrayCbk*/ 1) {
	driver.lefrArray.lefiArray::addTrack(driver.lefrTrackPatternPtr);
      }
    }
#line 11148 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1239: // array_rule: floorplan_start floorplan_list "END" "string"
#line 6758 "LefParser.yy"
    {
    delete (yystack_[0].value.stringVal);
    }
#line 11156 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1240: // $@115: %empty
#line 6761 "LefParser.yy"
                { if (/*driver.lefrArrayCbk*/ 1) driver.lefrDoGcell = 1; }
#line 11162 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1241: // array_rule: "GCELLGRID" $@115 gcellPattern ';'
#line 6762 "LefParser.yy"
    {
      if (/*driver.lefrArrayCbk*/ 1) {
	driver.lefrArray.lefiArray::addGcell(driver.lefrGcellPatternPtr);
      }
    }
#line 11172 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1242: // array_rule: "DEFAULTCAP" NUMBER cap_list "END" "DEFAULTCAP"
#line 6768 "LefParser.yy"
    {
      if (/*driver.lefrArrayCbk*/ 1) {
	driver.lefrArray.lefiArray::setTableSize((int)(yystack_[3].value.doubleVal));
      }
    }
#line 11182 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1243: // array_rule: def_statement
#line 6774 "LefParser.yy"
    { }
#line 11188 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1244: // floorplan_start: "FLOORPLAN" "string"
#line 6777 "LefParser.yy"
    { if (/*driver.lefrArrayCbk*/ 1) driver.lefrArray.lefiArray::addFloorPlan((*(yystack_[0].value.stringVal)).c_str()); delete (yystack_[0].value.stringVal);}
#line 11194 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1245: // floorplan_list: %empty
#line 6781 "LefParser.yy"
    { }
#line 11200 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1246: // floorplan_list: floorplan_list floorplan_element
#line 6783 "LefParser.yy"
    { }
#line 11206 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1247: // $@116: %empty
#line 6786 "LefParser.yy"
             { driver.lefDumbMode = 1; if (/*driver.lefrArrayCbk*/ 1) driver.lefrDoSite = 1; }
#line 11212 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1248: // floorplan_element: "CANPLACE" $@116 sitePattern ';'
#line 6788 "LefParser.yy"
    {
      if (/*driver.lefrArrayCbk*/ 1)
	driver.lefrArray.lefiArray::addSiteToFloorPlan("CANPLACE",
	driver.lefrSitePatternPtr);
    }
#line 11222 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1249: // $@117: %empty
#line 6793 "LefParser.yy"
                   { if (/*driver.lefrArrayCbk*/ 1) driver.lefrDoSite = 1; driver.lefDumbMode = 1; }
#line 11228 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1250: // floorplan_element: "CANNOTOCCUPY" $@117 sitePattern ';'
#line 6795 "LefParser.yy"
    {
      if (/*driver.lefrArrayCbk*/ 1)
	driver.lefrArray.lefiArray::addSiteToFloorPlan("CANNOTOCCUPY",
	driver.lefrSitePatternPtr);
     }
#line 11238 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1251: // cap_list: %empty
#line 6803 "LefParser.yy"
    { }
#line 11244 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1252: // cap_list: cap_list one_cap
#line 6805 "LefParser.yy"
    { }
#line 11250 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1253: // one_cap: "MINPINS" NUMBER "WIRECAP" NUMBER ';'
#line 6808 "LefParser.yy"
    { if (/*driver.lefrArrayCbk*/ 1) driver.lefrArray.lefiArray::addDefaultCap((int)(yystack_[3].value.doubleVal), (yystack_[1].value.doubleVal)); }
#line 11256 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1254: // $@118: %empty
#line 6811 "LefParser.yy"
            {driver.lefDumbMode=1;driver.lefNlToken=TRUE;}
#line 11262 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1255: // msg_statement: "MESSAGE" $@118 "string" '=' s_expr dtrm
#line 6812 "LefParser.yy"
    { driver.lefAddStringMessage((*(yystack_[3].value.stringVal)).c_str(), (*(yystack_[1].value.stringVal)).c_str()); delete (yystack_[3].value.stringVal); delete (yystack_[1].value.stringVal);}
#line 11268 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1256: // $@119: %empty
#line 6815 "LefParser.yy"
               {driver.lefDumbMode=1;driver.lefNlToken=TRUE;}
#line 11274 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1257: // create_file_statement: "CREATEFILE" $@119 "string" '=' s_expr dtrm
#line 6816 "LefParser.yy"
  {
  delete (yystack_[3].value.stringVal);
  }
#line 11282 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1258: // $@120: %empty
#line 6821 "LefParser.yy"
           {driver.lefDumbMode=1;driver.lefNlToken=TRUE;}
#line 11288 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1259: // def_statement: "DEFINE" $@120 "string" '=' expression dtrm
#line 6822 "LefParser.yy"
    {
      if (driver.versionNum < 5.6)
        driver.lefAddNumDefine((*(yystack_[3].value.stringVal)).c_str(), (yystack_[1].value.doubleVal));
      else
        if (/*driver.lefrArrayCbk*/ 1) /* write warning only if cbk is set */
           if (driver.arrayWarnings++ < driver.lefrArrayWarnings)
             driver.lefWarning(2067, "DEFINE statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    delete (yystack_[3].value.stringVal);
    }
#line 11302 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1260: // $@121: %empty
#line 6831 "LefParser.yy"
               {driver.lefDumbMode=1;driver.lefNlToken=TRUE;}
#line 11308 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1261: // def_statement: "DEFINES" $@121 "string" '=' s_expr dtrm
#line 6832 "LefParser.yy"
    {
      if (driver.versionNum < 5.6)
        driver.lefAddStringDefine((*(yystack_[3].value.stringVal)).c_str(), (*(yystack_[1].value.stringVal)).c_str());
      else
        if (/*driver.lefrArrayCbk*/ 1) /* write warning only if cbk is set */
           if (driver.arrayWarnings++ < driver.lefrArrayWarnings)
             driver.lefWarning(2068, "DEFINES statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    delete (yystack_[3].value.stringVal); delete (yystack_[1].value.stringVal);
    }
#line 11322 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1262: // $@122: %empty
#line 6841 "LefParser.yy"
               {driver.lefDumbMode=1;driver.lefNlToken=TRUE;}
#line 11328 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1263: // def_statement: "DEFINEB" $@122 "string" '=' b_expr dtrm
#line 6842 "LefParser.yy"
    {
      if (driver.versionNum < 5.6)
        driver.lefAddBooleanDefine((*(yystack_[3].value.stringVal)).c_str(), (yystack_[1].value.integerVal));
      else
        if (/*driver.lefrArrayCbk*/ 1) /* write warning only if cbk is set */
           if (driver.arrayWarnings++ < driver.lefrArrayWarnings)
             driver.lefWarning(2069, "DEFINEB statement is obsolete in version 5.6 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.6 or later.");
    delete (yystack_[3].value.stringVal);
    }
#line 11342 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1265: // dtrm: ';'
#line 6854 "LefParser.yy"
         {driver.lefNlToken = FALSE;}
#line 11348 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1266: // dtrm: '\n'
#line 6855 "LefParser.yy"
                {driver.lefNlToken = FALSE;}
#line 11354 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1271: // expression: expression '+' expression
#line 6868 "LefParser.yy"
                              {(yylhs.value.doubleVal) = (yystack_[2].value.doubleVal) + (yystack_[0].value.doubleVal); }
#line 11360 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1272: // expression: expression '-' expression
#line 6869 "LefParser.yy"
                                {(yylhs.value.doubleVal) = (yystack_[2].value.doubleVal) - (yystack_[0].value.doubleVal); }
#line 11366 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1273: // expression: expression '*' expression
#line 6870 "LefParser.yy"
                                {(yylhs.value.doubleVal) = (yystack_[2].value.doubleVal) * (yystack_[0].value.doubleVal); }
#line 11372 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1274: // expression: expression '/' expression
#line 6871 "LefParser.yy"
                                {(yylhs.value.doubleVal) = (yystack_[2].value.doubleVal) / (yystack_[0].value.doubleVal); }
#line 11378 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1275: // expression: '-' expression
#line 6872 "LefParser.yy"
                                {(yylhs.value.doubleVal) = -(yystack_[0].value.doubleVal);}
#line 11384 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1276: // expression: '(' expression ')'
#line 6873 "LefParser.yy"
                                {(yylhs.value.doubleVal) = (yystack_[1].value.doubleVal);}
#line 11390 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1277: // expression: "IF" b_expr then expression else expression
#line 6875 "LefParser.yy"
                {(yylhs.value.doubleVal) = ((yystack_[4].value.integerVal) != 0) ? (yystack_[2].value.doubleVal) : (yystack_[0].value.doubleVal);}
#line 11396 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1278: // expression: NUMBER
#line 6876 "LefParser.yy"
                                {(yylhs.value.doubleVal) = (yystack_[0].value.doubleVal);}
#line 11402 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1279: // b_expr: expression relop expression
#line 6879 "LefParser.yy"
                              {(yylhs.value.integerVal) = driver.comp_num((yystack_[2].value.doubleVal),(yystack_[1].value.integerVal),(yystack_[0].value.doubleVal));}
#line 11408 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1280: // b_expr: expression "AND" expression
#line 6880 "LefParser.yy"
                                {(yylhs.value.integerVal) = (yystack_[2].value.doubleVal) != 0 && (yystack_[0].value.doubleVal) != 0;}
#line 11414 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1281: // b_expr: expression "OR" expression
#line 6881 "LefParser.yy"
                                {(yylhs.value.integerVal) = (yystack_[2].value.doubleVal) != 0 || (yystack_[0].value.doubleVal) != 0;}
#line 11420 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1282: // b_expr: s_expr relop s_expr
#line 6882 "LefParser.yy"
                              {(yylhs.value.integerVal) = driver.comp_str((*(yystack_[2].value.stringVal)).c_str(),(yystack_[1].value.integerVal),(*(yystack_[0].value.stringVal)).c_str()); delete (yystack_[2].value.stringVal); delete (yystack_[0].value.stringVal);}
#line 11426 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1283: // b_expr: s_expr "AND" s_expr
#line 6883 "LefParser.yy"
                              {(yylhs.value.integerVal) = (*(yystack_[2].value.stringVal))[0] != 0 && (*(yystack_[0].value.stringVal))[0] != 0; delete (yystack_[2].value.stringVal); delete (yystack_[0].value.stringVal);}
#line 11432 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1284: // b_expr: s_expr "OR" s_expr
#line 6884 "LefParser.yy"
                              {(yylhs.value.integerVal) = (*(yystack_[2].value.stringVal))[0] != 0 || (*(yystack_[0].value.stringVal))[0] != 0; delete (yystack_[2].value.stringVal); delete (yystack_[0].value.stringVal);}
#line 11438 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1285: // b_expr: b_expr K_EQ b_expr
#line 6885 "LefParser.yy"
                              {(yylhs.value.integerVal) = (yystack_[2].value.integerVal) == (yystack_[0].value.integerVal);}
#line 11444 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1286: // b_expr: b_expr K_NE b_expr
#line 6886 "LefParser.yy"
                              {(yylhs.value.integerVal) = (yystack_[2].value.integerVal) != (yystack_[0].value.integerVal);}
#line 11450 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1287: // b_expr: b_expr "AND" b_expr
#line 6887 "LefParser.yy"
                              {(yylhs.value.integerVal) = (yystack_[2].value.integerVal) && (yystack_[0].value.integerVal);}
#line 11456 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1288: // b_expr: b_expr "OR" b_expr
#line 6888 "LefParser.yy"
                              {(yylhs.value.integerVal) = (yystack_[2].value.integerVal) || (yystack_[0].value.integerVal);}
#line 11462 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1289: // b_expr: K_NOT b_expr
#line 6889 "LefParser.yy"
                                                 {(yylhs.value.integerVal) = !(yylhs.value.integerVal);}
#line 11468 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1290: // b_expr: '(' b_expr ')'
#line 6890 "LefParser.yy"
                              {(yylhs.value.integerVal) = (yystack_[1].value.integerVal);}
#line 11474 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1291: // b_expr: "IF" b_expr then b_expr else b_expr
#line 6892 "LefParser.yy"
          {(yylhs.value.integerVal) = ((yystack_[4].value.integerVal) != 0) ? (yystack_[2].value.integerVal) : (yystack_[0].value.integerVal);}
#line 11480 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1292: // b_expr: "TRUE"
#line 6893 "LefParser.yy"
                              {(yylhs.value.integerVal) = 1;}
#line 11486 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1293: // b_expr: "FALSE"
#line 6894 "LefParser.yy"
                              {(yylhs.value.integerVal) = 0;}
#line 11492 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1294: // s_expr: s_expr '+' s_expr
#line 6898 "LefParser.yy"
    {
	  (yylhs.value.stringVal) = new std::string ((*(yystack_[2].value.stringVal)) + (*(yystack_[0].value.stringVal)));
      delete (yystack_[2].value.stringVal); delete (yystack_[0].value.stringVal);
    }
#line 11501 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1295: // s_expr: '(' s_expr ')'
#line 6903 "LefParser.yy"
    { (yylhs.value.stringVal) = (yystack_[1].value.stringVal); }
#line 11507 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1296: // s_expr: "IF" b_expr then s_expr else s_expr
#line 6905 "LefParser.yy"
    {
      driver.lefDefIf = TRUE;
      if ((yystack_[4].value.integerVal) != 0) {
	(yylhs.value.stringVal) = (yystack_[2].value.stringVal);	
      } else {
	(yylhs.value.stringVal) = (yystack_[0].value.stringVal);
      }
    }
#line 11520 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1297: // s_expr: "qstring"
#line 6914 "LefParser.yy"
    { (yylhs.value.stringVal) = (yystack_[0].value.qstringVal); }
#line 11526 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1298: // relop: K_LE
#line 6917 "LefParser.yy"
       {(yylhs.value.integerVal) = C_LE;}
#line 11532 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1299: // relop: K_LT
#line 6918 "LefParser.yy"
         {(yylhs.value.integerVal) = C_LT;}
#line 11538 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1300: // relop: K_GE
#line 6919 "LefParser.yy"
         {(yylhs.value.integerVal) = C_GE;}
#line 11544 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1301: // relop: K_GT
#line 6920 "LefParser.yy"
         {(yylhs.value.integerVal) = C_GT;}
#line 11550 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1302: // relop: K_EQ
#line 6921 "LefParser.yy"
         {(yylhs.value.integerVal) = C_EQ;}
#line 11556 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1303: // relop: K_NE
#line 6922 "LefParser.yy"
         {(yylhs.value.integerVal) = C_NE;}
#line 11562 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1304: // relop: '='
#line 6923 "LefParser.yy"
         {(yylhs.value.integerVal) = C_EQ;}
#line 11568 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1305: // relop: '<'
#line 6924 "LefParser.yy"
         {(yylhs.value.integerVal) = C_LT;}
#line 11574 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1306: // relop: '>'
#line 6925 "LefParser.yy"
         {(yylhs.value.integerVal) = C_GT;}
#line 11580 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1307: // $@123: %empty
#line 6929 "LefParser.yy"
    { 
      if (/*lefrPropBeginCbk*/ 1)
        driver.lefrPropBeginCbk( 0);
      driver.lefInPropDef = 1;  /* set flag as inside propertydefinitions */
    }
#line 11590 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1308: // prop_def_section: "PROPERTYDEFINITIONS" $@123 prop_stmts "END" "PROPERTYDEFINITIONS"
#line 6935 "LefParser.yy"
    { 
      if (/*lefrPropEndCbk*/ 1)
        driver.lefrPropEndCbk( 0);
      driver.lefRealNum = 0;     /* just want to make sure it is reset */
      driver.lefInPropDef = 0;   /* reset flag */
    }
#line 11601 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1309: // prop_stmts: %empty
#line 6944 "LefParser.yy"
    { }
#line 11607 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1310: // prop_stmts: prop_stmts prop_stmt
#line 6946 "LefParser.yy"
    { }
#line 11613 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1311: // $@124: %empty
#line 6949 "LefParser.yy"
            {driver.lefDumbMode = 1; driver.lefrProp.lefiProp::clear(); }
#line 11619 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1312: // prop_stmt: "LIBRARY" $@124 "string" prop_define ';'
#line 6951 "LefParser.yy"
    { 
      if (/*lefrPropCbk*/ 1) {
        driver.lefrProp.lefiProp::setPropType("library", (*(yystack_[2].value.stringVal)).c_str());
        driver.lefrPropCbk( driver.lefrProp);
      }
      driver.lefrLibProp.lefiPropType::setPropType((*(yystack_[2].value.stringVal)).c_str(), driver.lefPropDefType);
      delete (yystack_[2].value.stringVal);
    }
#line 11632 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1313: // $@125: %empty
#line 6959 "LefParser.yy"
                   {driver.lefDumbMode = 1; driver.lefrProp.lefiProp::clear(); }
#line 11638 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1314: // prop_stmt: "COMPONENTPIN" $@125 "string" prop_define ';'
#line 6961 "LefParser.yy"
    { 
      if (/*lefrPropCbk*/ 1) {
        driver.lefrProp.lefiProp::setPropType("componentpin", (*(yystack_[2].value.stringVal)).c_str());
        driver.lefrPropCbk( driver.lefrProp);
      }
      driver.lefrCompProp.lefiPropType::setPropType((*(yystack_[2].value.stringVal)).c_str(), driver.lefPropDefType);
      delete (yystack_[2].value.stringVal);
    }
#line 11651 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1315: // $@126: %empty
#line 6969 "LefParser.yy"
          {driver.lefDumbMode = 1; driver.lefrProp.lefiProp::clear(); }
#line 11657 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1316: // prop_stmt: "PIN" $@126 "string" prop_define ';'
#line 6971 "LefParser.yy"
    { 
      if (/*lefrPropCbk*/ 1) {
        driver.lefrProp.lefiProp::setPropType("pin", (*(yystack_[2].value.stringVal)).c_str());
        driver.lefrPropCbk( driver.lefrProp);
      }
      driver.lefrPinProp.lefiPropType::setPropType((*(yystack_[2].value.stringVal)).c_str(), driver.lefPropDefType);
      
      delete (yystack_[2].value.stringVal);
    }
#line 11671 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1317: // $@127: %empty
#line 6980 "LefParser.yy"
            {driver.lefDumbMode = 1; driver.lefrProp.lefiProp::clear(); }
#line 11677 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1318: // prop_stmt: "MACRO" $@127 "string" prop_define ';'
#line 6982 "LefParser.yy"
    { 
      if (/*lefrPropCbk*/ 1) {
        driver.lefrProp.lefiProp::setPropType("macro", (*(yystack_[2].value.stringVal)).c_str());
        driver.lefrPropCbk( driver.lefrProp);
      }
      driver.lefrMacroProp.lefiPropType::setPropType((*(yystack_[2].value.stringVal)).c_str(), driver.lefPropDefType);
      delete (yystack_[2].value.stringVal);
    }
#line 11690 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1319: // $@128: %empty
#line 6990 "LefParser.yy"
          {driver.lefDumbMode = 1; driver.lefrProp.lefiProp::clear(); }
#line 11696 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1320: // prop_stmt: "VIA" $@128 "string" prop_define ';'
#line 6992 "LefParser.yy"
    { 
      if (/*lefrPropCbk*/ 1) {
        driver.lefrProp.lefiProp::setPropType("via", (*(yystack_[2].value.stringVal)).c_str());
        driver.lefrPropCbk( driver.lefrProp);
      }
      driver.lefrViaProp.lefiPropType::setPropType((*(yystack_[2].value.stringVal)).c_str(), driver.lefPropDefType);
      delete (yystack_[2].value.stringVal);
    }
#line 11709 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1321: // $@129: %empty
#line 7000 "LefParser.yy"
              {driver.lefDumbMode = 1; driver.lefrProp.lefiProp::clear(); }
#line 11715 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1322: // prop_stmt: "VIARULE" $@129 "string" prop_define ';'
#line 7002 "LefParser.yy"
    { 
      if (/*lefrPropCbk*/ 1) {
        driver.lefrProp.lefiProp::setPropType("viarule", (*(yystack_[2].value.stringVal)).c_str());
        driver.lefrPropCbk( driver.lefrProp);
      }
      driver.lefrViaRuleProp.lefiPropType::setPropType((*(yystack_[2].value.stringVal)).c_str(), driver.lefPropDefType);
      delete (yystack_[2].value.stringVal);
    }
#line 11728 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1323: // $@130: %empty
#line 7010 "LefParser.yy"
            {driver.lefDumbMode = 1; driver.lefrProp.lefiProp::clear(); }
#line 11734 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1324: // prop_stmt: "LAYER" $@130 "string" prop_define ';'
#line 7012 "LefParser.yy"
    { 
      if (/*lefrPropCbk*/ 1) {
        driver.lefrProp.lefiProp::setPropType("layer", (*(yystack_[2].value.stringVal)).c_str());
        driver.lefrPropCbk( driver.lefrProp);
      }
      driver.lefrLayerProp.lefiPropType::setPropType((*(yystack_[2].value.stringVal)).c_str(), driver.lefPropDefType);
      delete (yystack_[2].value.stringVal);
    }
#line 11747 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1325: // $@131: %empty
#line 7020 "LefParser.yy"
                     {driver.lefDumbMode = 1; driver.lefrProp.lefiProp::clear(); }
#line 11753 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1326: // prop_stmt: "NONDEFAULTRULE" $@131 "string" prop_define ';'
#line 7022 "LefParser.yy"
    { 
      if (/*lefrPropCbk*/ 1) {
        driver.lefrProp.lefiProp::setPropType("nondefaultrule", (*(yystack_[2].value.stringVal)).c_str());
        driver.lefrPropCbk( driver.lefrProp);
      }
      driver.lefrNondefProp.lefiPropType::setPropType((*(yystack_[2].value.stringVal)).c_str(), driver.lefPropDefType);
      delete (yystack_[2].value.stringVal);
    }
#line 11766 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1327: // $@132: %empty
#line 7032 "LefParser.yy"
            { driver.lefRealNum = 0; }
#line 11772 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1328: // prop_define: "INTEGER" $@132 opt_def_range opt_def_value
#line 7033 "LefParser.yy"
    { 
      if (/*lefrPropCbk*/ 1) driver.lefrProp.lefiProp::setPropInteger();
      driver.lefPropDefType = 'I';
    }
#line 11781 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1329: // $@133: %empty
#line 7037 "LefParser.yy"
           { driver.lefRealNum = 1; }
#line 11787 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1330: // prop_define: "REAL" $@133 opt_def_range opt_def_value
#line 7038 "LefParser.yy"
    { 
      if (/*lefrPropCbk*/ 1) driver.lefrProp.lefiProp::setPropReal();
      driver.lefPropDefType = 'R';
      driver.lefRealNum = 0;
    }
#line 11797 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1331: // prop_define: "STRING"
#line 7044 "LefParser.yy"
    {
      if (/*lefrPropCbk*/ 1) driver.lefrProp.lefiProp::setPropString();
      driver.lefPropDefType = 'S';
    }
#line 11806 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1332: // prop_define: "STRING" "qstring"
#line 7049 "LefParser.yy"
    {
      if (/*lefrPropCbk*/ 1) driver.lefrProp.lefiProp::setPropQString((*(yystack_[0].value.qstringVal)).c_str());
      driver.lefPropDefType = 'Q';
      delete (yystack_[0].value.qstringVal);
    }
#line 11816 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1333: // prop_define: "NAMEMAPSTRING" "string"
#line 7055 "LefParser.yy"
    {
      if (/*lefrPropCbk*/ 1) driver.lefrProp.lefiProp::setPropNameMapString((*(yystack_[0].value.stringVal)).c_str());
      driver.lefPropDefType = 'S';
      delete (yystack_[0].value.stringVal);
    }
#line 11826 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1334: // opt_range_second: %empty
#line 7063 "LefParser.yy"
    { }
#line 11832 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1335: // opt_range_second: "USELENGTHTHRESHOLD"
#line 7065 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
        driver.lefrLayer.lefiLayer::setSpacingRangeUseLength();
    }
#line 11841 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1336: // opt_range_second: "INFLUENCE" NUMBER
#line 7070 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) {
        driver.lefrLayer.lefiLayer::setSpacingRangeInfluence((yystack_[0].value.doubleVal));
        driver.lefrLayer.lefiLayer::setSpacingRangeInfluenceRange(-1, -1);
      }
    }
#line 11852 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1337: // opt_range_second: "INFLUENCE" NUMBER "RANGE" NUMBER NUMBER
#line 7077 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) {
        driver.lefrLayer.lefiLayer::setSpacingRangeInfluence((yystack_[3].value.doubleVal));
        driver.lefrLayer.lefiLayer::setSpacingRangeInfluenceRange((yystack_[1].value.doubleVal), (yystack_[0].value.doubleVal));
      }
    }
#line 11863 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1338: // opt_range_second: "RANGE" NUMBER NUMBER
#line 7084 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
        driver.lefrLayer.lefiLayer::setSpacingRangeRange((yystack_[1].value.doubleVal), (yystack_[0].value.doubleVal));
    }
#line 11872 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1339: // opt_endofline: %empty
#line 7091 "LefParser.yy"
    { }
#line 11878 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1340: // $@134: %empty
#line 7093 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
        driver.lefrLayer.lefiLayer::setSpacingParSW((yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
    }
#line 11887 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1342: // opt_endofline_twoedges: %empty
#line 7101 "LefParser.yy"
    { }
#line 11893 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1343: // opt_endofline_twoedges: "TWOEDGES"
#line 7103 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
        driver.lefrLayer.lefiLayer::setSpacingParTwoEdges();
    }
#line 11902 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1344: // opt_samenetPGonly: %empty
#line 7110 "LefParser.yy"
    { }
#line 11908 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1345: // opt_samenetPGonly: "PGONLY"
#line 7112 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
        driver.lefrLayer.lefiLayer::setSpacingSamenetPGonly();
    }
#line 11917 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1346: // opt_def_range: %empty
#line 7119 "LefParser.yy"
    { }
#line 11923 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1347: // opt_def_range: "RANGE" NUMBER NUMBER
#line 7121 "LefParser.yy"
    {  if (/*lefrPropCbk*/ 1) driver.lefrProp.lefiProp::setRange((yystack_[1].value.doubleVal), (yystack_[0].value.doubleVal)); }
#line 11929 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1348: // opt_def_value: %empty
#line 7125 "LefParser.yy"
    { }
#line 11935 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1349: // opt_def_value: NUMBER
#line 7127 "LefParser.yy"
    { if (/*lefrPropCbk*/ 1) driver.lefrProp.lefiProp::setNumber((yystack_[0].value.doubleVal)); }
#line 11941 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1352: // layer_spacing_opt: "CENTERTOCENTER"
#line 7134 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) {
         if (driver.hasSpCenter) {
           if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1663, "A CENTERTOCENTER statement has already defined in SPACING\nCENTERTOCENTER can only be defined once per LAYER CUT SPACING.");
              /*CHKERR();*/
           }
        }
        driver.hasSpCenter = 1;
        if (driver.versionNum < 5.6) {
           if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "CENTERTOCENTER statement is a version 5.6 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1664, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
        if (/*driver.lefrLayerCbk*/ 1)
          driver.lefrLayer.lefiLayer::setSpacingCenterToCenter();
      }
    }
#line 11969 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1353: // $@135: %empty
#line 7158 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) {
        if (driver.hasSpSamenet) {
           if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.lefError(1665, "A SAMENET statement has already defined in SPACING\nSAMENET can only be defined once per LAYER CUT SPACING.");
              /*CHKERR();*/
           }
        }
        driver.hasSpSamenet = 1;
        if (/*driver.lefrLayerCbk*/ 1)
          driver.lefrLayer.lefiLayer::setSpacingSamenet();
       }
    }
#line 11987 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1354: // layer_spacing_opt: "SAMENET" $@135 opt_samenetPGonly
#line 7172 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
        driver.outMsg = (char*)lefMalloc(10000);
        sprintf(driver.outMsg,
          "SAMENET is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
        driver.lefError(1684, driver.outMsg);
        lefFree(driver.outMsg);
        /*CHKERR();*/
      }
    }
#line 12002 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1355: // layer_spacing_opt: "PARALLELOVERLAP"
#line 7183 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
        driver.outMsg = (char*)lefMalloc(10000);
        sprintf(driver.outMsg,
          "PARALLELOVERLAP is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
        driver.lefError(1680, driver.outMsg);
        lefFree(driver.outMsg);
        /*CHKERR();*/ 
      } else {
        if (/*driver.lefrLayerCbk*/ 1) {
          if (driver.hasSpParallel) {
             if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                driver.lefError(1666, "A PARALLELOVERLAP statement has already defined in SPACING\nPARALLELOVERLAP can only be defined once per LAYER CUT SPACING.");
                /*CHKERR();*/
             }
          }
          driver.hasSpParallel = 1;
          if (/*driver.lefrLayerCbk*/ 1)
            driver.lefrLayer.lefiLayer::setSpacingParallelOverlap();
        }
      }
    }
#line 12029 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1357: // $@136: %empty
#line 7208 "LefParser.yy"
            {driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 12035 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1358: // $@137: %empty
#line 7209 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
{
        if (driver.versionNum < 5.7) {
           if (driver.hasSpSamenet) {    /* 5.6 and earlier does not allow */
              if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                 driver.lefError(1667, "A SAMENET statement has already defined in SPACING\nEither SAMENET or LAYER can be defined, but not both.");
                 /*CHKERR();*/
              }
           }
        }
        driver.lefrLayer.lefiLayer::setSpacingName((*(yystack_[0].value.stringVal)).c_str());
      }
      delete (yystack_[0].value.stringVal);
    }
#line 12055 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1360: // $@138: %empty
#line 7226 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) {
        if (driver.versionNum < 5.5) {
           if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
              driver.outMsg = (char*)lefMalloc(10000);
              sprintf (driver.outMsg,
                 "ADJACENTCUTS statement is a version 5.5 and later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
              driver.lefError(1668, driver.outMsg);
              lefFree(driver.outMsg);
              /*CHKERR();*/
           }
        }
        if (driver.versionNum < 5.7) {
           if (driver.hasSpSamenet) {    /* 5.6 and earlier does not allow */
              if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                 driver.lefError(1669, "A SAMENET statement has already defined in SPACING\nEither SAMENET or ADJACENTCUTS can be defined, but not both.");
                 /*CHKERR();*/
              }
           }
        }
        driver.lefrLayer.lefiLayer::setSpacingAdjacent((int)(yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
      }
    }
#line 12083 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1362: // layer_spacing_cut_routing: "AREA" NUMBER
#line 7251 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
        driver.outMsg = (char*)lefMalloc(10000);
        sprintf(driver.outMsg,
          "AREA is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
        driver.lefError(1693, driver.outMsg);
        lefFree(driver.outMsg);
        /*CHKERR();*/
      } else {
        if (/*driver.lefrLayerCbk*/ 1) {
          if (driver.versionNum < 5.7) {
             if (driver.hasSpSamenet) {    /* 5.6 and earlier does not allow */
                if (driver.layerWarnings++ < driver.lefrLayerWarnings) {
                   driver.lefError(1670, "A SAMENET statement has already defined in SPACING\nEither SAMENET or AREA can be defined, but not both.");
                   /*CHKERR();*/
                }
             }
          }
          driver.lefrLayer.lefiLayer::setSpacingArea((yystack_[0].value.doubleVal));
        }
      }
    }
#line 12110 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1363: // $@139: %empty
#line 7274 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
        driver.lefrLayer.lefiLayer::setSpacingRange((yystack_[1].value.doubleVal), (yystack_[0].value.doubleVal));
    }
#line 12119 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1365: // layer_spacing_cut_routing: "LENGTHTHRESHOLD" NUMBER
#line 7280 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) {
        driver.lefrLayer.lefiLayer::setSpacingLength((yystack_[0].value.doubleVal));
      }
    }
#line 12129 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1366: // layer_spacing_cut_routing: "LENGTHTHRESHOLD" NUMBER "RANGE" NUMBER NUMBER
#line 7286 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1) {
        driver.lefrLayer.lefiLayer::setSpacingLength((yystack_[3].value.doubleVal));
        driver.lefrLayer.lefiLayer::setSpacingLengthRange((yystack_[1].value.doubleVal), (yystack_[0].value.doubleVal));
      }
    }
#line 12140 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1367: // $@140: %empty
#line 7293 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
        driver.lefrLayer.lefiLayer::setSpacingEol((yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
    }
#line 12149 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1368: // layer_spacing_cut_routing: "ENDOFLINE" NUMBER "WITHIN" NUMBER $@140 opt_endofline
#line 7298 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
        driver.outMsg = (char*)lefMalloc(10000);
        sprintf(driver.outMsg,
          "ENDOFLINE is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
        driver.lefError(1681, driver.outMsg);
        lefFree(driver.outMsg);
        /*CHKERR();*/
      }
    }
#line 12164 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1369: // layer_spacing_cut_routing: "NOTCHLENGTH" NUMBER
#line 7309 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
        driver.outMsg = (char*)lefMalloc(10000);
        sprintf(driver.outMsg,
          "NOTCHLENGTH is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
        driver.lefError(1682, driver.outMsg);
        lefFree(driver.outMsg);
        /*CHKERR();*/
      } else {
        if (/*driver.lefrLayerCbk*/ 1)
          driver.lefrLayer.lefiLayer::setSpacingNotchLength((yystack_[0].value.doubleVal));
      }
    }
#line 12182 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1370: // layer_spacing_cut_routing: "ENDOFNOTCHWIDTH" NUMBER "NOTCHSPACING" NUMBER "NOTCHLENGTH" NUMBER
#line 7323 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
        driver.outMsg = (char*)lefMalloc(10000);
        sprintf(driver.outMsg,
          "ENDOFNOTCHWIDTH is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
        driver.lefError(1696, driver.outMsg);
        lefFree(driver.outMsg);
        /*CHKERR();*/
      } else {
        if (/*driver.lefrLayerCbk*/ 1)
          driver.lefrLayer.lefiLayer::setSpacingEndOfNotchWidth((yystack_[4].value.doubleVal), (yystack_[2].value.doubleVal), (yystack_[0].value.doubleVal));
      }
    }
#line 12200 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1371: // spacing_cut_layer_opt: %empty
#line 7339 "LefParser.yy"
    {}
#line 12206 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1372: // spacing_cut_layer_opt: "STACK"
#line 7341 "LefParser.yy"
    {
      if (/*driver.lefrLayerCbk*/ 1)
        driver.lefrLayer.lefiLayer::setSpacingLayerStack();
    }
#line 12215 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1373: // opt_adjacentcuts_exceptsame: %empty
#line 7348 "LefParser.yy"
    {}
#line 12221 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1374: // opt_adjacentcuts_exceptsame: "EXCEPTSAMEPGNET"
#line 7350 "LefParser.yy"
    {
      if (driver.versionNum < 5.7) {
        driver.outMsg = (char*)lefMalloc(10000);
        sprintf(driver.outMsg,
          "EXCEPTSAMEPGNET is a version 5.7 or later syntax.\nYour lef file is defined with version %g.", driver.versionNum);
        driver.lefError(1683, driver.outMsg);
        lefFree(driver.outMsg);
        /*CHKERR();*/
      } else {
        if (/*driver.lefrLayerCbk*/ 1)
          driver.lefrLayer.lefiLayer::setSpacingAdjacentExcept();
      }
    }
#line 12239 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1375: // opt_layer_name: %empty
#line 7366 "LefParser.yy"
    { (yylhs.value.stringVal) = new std::string (""); }
#line 12245 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1376: // $@141: %empty
#line 7367 "LefParser.yy"
            {driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 12251 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1377: // opt_layer_name: "LAYER" $@141 "string"
#line 7368 "LefParser.yy"
    { (yylhs.value.stringVal) = (yystack_[0].value.stringVal); }
#line 12257 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1378: // $@142: %empty
#line 7372 "LefParser.yy"
           {driver.lefDumbMode = 1; driver.lefNoNum = 1; }
#line 12263 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1379: // req_layer_name: "LAYER" $@142 "string"
#line 7373 "LefParser.yy"
    { (yylhs.value.stringVal) = (yystack_[0].value.stringVal); }
#line 12269 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1380: // universalnoisemargin: "UNIVERSALNOISEMARGIN" NUMBER NUMBER ';'
#line 7377 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*driver.lefrNoiseMarginCbk*/ 1) {
          driver.lefrNoiseMargin.low = (yystack_[2].value.doubleVal);
          driver.lefrNoiseMargin.high = (yystack_[1].value.doubleVal);
          driver.lefrNoiseMarginCbk( driver.lefrNoiseMargin);
        }
      } else
        if (/*driver.lefrNoiseMarginCbk*/ 1) /* write warning only if cbk is set */
          if (driver.noiseMarginWarnings++ < driver.lefrNoiseMarginWarnings)
            driver.lefWarning(2070, "UNIVERSALNOISEMARGIN statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 12286 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1381: // edgeratethreshold1: "EDGERATETHRESHOLD1" NUMBER ';'
#line 7391 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*lefrEdgeRateThreshold1Cbk*/ 1) {
		  driver.lefrEdgeRateThreshold1Cbk((yystack_[1].value.doubleVal));
        }
      } else
        if (/*lefrEdgeRateThreshold1Cbk*/ 1) /* write warning only if cbk is set */
          if (driver.edgeRateThreshold1Warnings++ < driver.lefrEdgeRateThreshold1Warnings)
            driver.lefWarning(2071, "EDGERATETHRESHOLD1 statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 12301 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1382: // edgeratethreshold2: "EDGERATETHRESHOLD2" NUMBER ';'
#line 7403 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*lefrEdgeRateThreshold2Cbk*/ 1) {
		driver.lefrEdgeRateThreshold2Cbk((yystack_[1].value.doubleVal));
        }
      } else
        if (/*lefrEdgeRateThreshold2Cbk*/ 1) /* write warning only if cbk is set */
          if (driver.edgeRateThreshold2Warnings++ < driver.lefrEdgeRateThreshold2Warnings)
            driver.lefWarning(2072, "EDGERATETHRESHOLD2 statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 12316 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1383: // edgeratescalefactor: "EDGERATESCALEFACTOR" NUMBER ';'
#line 7415 "LefParser.yy"
    {
      if (driver.versionNum < 5.4) {
        if (/*lefrEdgeRateScaleFactorCbk*/ 1) {
		driver.lefrEdgeRateScaleFactorCbk((yystack_[1].value.doubleVal));
        }
      } else
        if (/*lefrEdgeRateScaleFactorCbk*/ 1) /* write warning only if cbk is set */
          if (driver.edgeRateScaleFactorWarnings++ < driver.lefrEdgeRateScaleFactorWarnings)
            driver.lefWarning(2073, "EDGERATESCALEFACTOR statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
    }
#line 12331 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1384: // $@143: %empty
#line 7427 "LefParser.yy"
    { if (/*driver.lefrNoiseTableCbk*/ 1) driver.lefrNoiseTable.lefiNoiseTable::setup((int)(yystack_[0].value.doubleVal)); }
#line 12337 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1385: // noisetable: "NOISETABLE" NUMBER $@143 ';' noise_table_list end_noisetable dtrm
#line 7429 "LefParser.yy"
    { }
#line 12343 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1386: // end_noisetable: "END" "NOISETABLE"
#line 7433 "LefParser.yy"
  {
    if (driver.versionNum < 5.4) {
      if (/*driver.lefrNoiseTableCbk*/ 1)
        driver.lefrNoiseTableCbk( driver.lefrNoiseTable);
    } else
      if (/*driver.lefrNoiseTableCbk*/ 1) /* write warning only if cbk is set */
        if (driver.noiseTableWarnings++ < driver.lefrNoiseTableWarnings)
          driver.lefWarning(2074, "NOISETABLE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
  }
#line 12357 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1389: // noise_table_entry: "EDGERATE" NUMBER ';'
#line 7451 "LefParser.yy"
    { if (/*driver.lefrNoiseTableCbk*/ 1)
         {
            driver.lefrNoiseTable.lefiNoiseTable::newEdge();
            driver.lefrNoiseTable.lefiNoiseTable::addEdge((yystack_[1].value.doubleVal));
         }
    }
#line 12368 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1390: // noise_table_entry: output_resistance_entry
#line 7458 "LefParser.yy"
    { }
#line 12374 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1391: // $@144: %empty
#line 7461 "LefParser.yy"
    { if (/*driver.lefrNoiseTableCbk*/ 1) driver.lefrNoiseTable.lefiNoiseTable::addResistance(); }
#line 12380 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1393: // num_list: NUMBER
#line 7467 "LefParser.yy"
    { if (/*driver.lefrNoiseTableCbk*/ 1)
    driver.lefrNoiseTable.lefiNoiseTable::addResistanceNumber((yystack_[0].value.doubleVal)); }
#line 12387 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1394: // num_list: num_list NUMBER
#line 7470 "LefParser.yy"
    { if (/*driver.lefrNoiseTableCbk*/ 1)
    driver.lefrNoiseTable.lefiNoiseTable::addResistanceNumber((yystack_[0].value.doubleVal)); }
#line 12394 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1397: // $@145: %empty
#line 7479 "LefParser.yy"
        { if (/*driver.lefrNoiseTableCbk*/ 1)
	driver.lefrNoiseTable.lefiNoiseTable::addVictimLength((yystack_[1].value.doubleVal)); }
#line 12401 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1398: // victim: "VICTIMLENGTH" NUMBER ';' $@145 "VICTIMNOISE" vnoiselist ';'
#line 7482 "LefParser.yy"
        { }
#line 12407 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1399: // vnoiselist: NUMBER
#line 7486 "LefParser.yy"
    { if (/*driver.lefrNoiseTableCbk*/ 1)
    driver.lefrNoiseTable.lefiNoiseTable::addVictimNoise((yystack_[0].value.doubleVal)); }
#line 12414 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1400: // vnoiselist: vnoiselist NUMBER
#line 7489 "LefParser.yy"
    { if (/*driver.lefrNoiseTableCbk*/ 1)
    driver.lefrNoiseTable.lefiNoiseTable::addVictimNoise((yystack_[0].value.doubleVal)); }
#line 12421 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1401: // $@146: %empty
#line 7493 "LefParser.yy"
    { if (/*driver.lefrCorrectionTableCbk*/ 1)
    driver.lefrCorrectionTable.lefiCorrectionTable::setup((int)(yystack_[1].value.doubleVal)); }
#line 12428 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1402: // correctiontable: "CORRECTIONTABLE" NUMBER ';' $@146 correction_table_list end_correctiontable dtrm
#line 7496 "LefParser.yy"
    { }
#line 12434 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1403: // end_correctiontable: "END" "CORRECTIONTABLE"
#line 7500 "LefParser.yy"
  {
    if (driver.versionNum < 5.4) {
      if (/*driver.lefrCorrectionTableCbk*/ 1)
        driver.lefrCorrectionTableCbk(driver.lefrCorrectionTable);
    } else
      if (/*driver.lefrCorrectionTableCbk*/ 1) /* write warning only if cbk is set */
        if (driver.correctionTableWarnings++ < driver.lefrCorrectionTableWarnings)
          driver.lefWarning(2075, "CORRECTIONTABLE statement is obsolete in version 5.4 and later.\nThe LEF parser will ignore this statement.\nTo avoid this warning in the future, remove this statement from the LEF file with version 5.4 or later.");
  }
#line 12448 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1406: // correction_table_item: "EDGERATE" NUMBER ';'
#line 7517 "LefParser.yy"
    { if (/*driver.lefrCorrectionTableCbk*/ 1)
         {
            driver.lefrCorrectionTable.lefiCorrectionTable::newEdge();
            driver.lefrCorrectionTable.lefiCorrectionTable::addEdge((yystack_[1].value.doubleVal));
         }
    }
#line 12459 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1407: // correction_table_item: output_list
#line 7524 "LefParser.yy"
    { }
#line 12465 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1408: // $@147: %empty
#line 7527 "LefParser.yy"
  { if (/*driver.lefrCorrectionTableCbk*/ 1)
  driver.lefrCorrectionTable.lefiCorrectionTable::addResistance(); }
#line 12472 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1409: // output_list: "OUTPUTRESISTANCE" $@147 numo_list ';' corr_victim_list
#line 7530 "LefParser.yy"
  { }
#line 12478 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1410: // numo_list: NUMBER
#line 7534 "LefParser.yy"
    { if (/*driver.lefrCorrectionTableCbk*/ 1)
    driver.lefrCorrectionTable.lefiCorrectionTable::addResistanceNumber((yystack_[0].value.doubleVal)); }
#line 12485 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1411: // numo_list: numo_list NUMBER
#line 7537 "LefParser.yy"
    { if (/*driver.lefrCorrectionTableCbk*/ 1)
    driver.lefrCorrectionTable.lefiCorrectionTable::addResistanceNumber((yystack_[0].value.doubleVal)); }
#line 12492 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1414: // $@148: %empty
#line 7547 "LefParser.yy"
     { if (/*driver.lefrCorrectionTableCbk*/ 1)
     driver.lefrCorrectionTable.lefiCorrectionTable::addVictimLength((yystack_[1].value.doubleVal)); }
#line 12499 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1415: // corr_victim: "VICTIMLENGTH" NUMBER ';' $@148 "CORRECTIONFACTOR" corr_list ';'
#line 7550 "LefParser.yy"
     { }
#line 12505 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1416: // corr_list: NUMBER
#line 7554 "LefParser.yy"
    { if (/*driver.lefrCorrectionTableCbk*/ 1)
	driver.lefrCorrectionTable.lefiCorrectionTable::addVictimCorrection((yystack_[0].value.doubleVal)); }
#line 12512 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1417: // corr_list: corr_list NUMBER
#line 7557 "LefParser.yy"
    { if (/*driver.lefrCorrectionTableCbk*/ 1)
	driver.lefrCorrectionTable.lefiCorrectionTable::addVictimCorrection((yystack_[0].value.doubleVal)); }
#line 12519 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1418: // input_antenna: "INPUTPINANTENNASIZE" NUMBER ';'
#line 7563 "LefParser.yy"
    { /* 5.3 syntax */
        driver.use5_3 = 1;
        if (driver.ignoreVersion) {
           /* do nothing */
        } else if (driver.versionNum > 5.3) {
           /* A 5.3 syntax in 5.4 */
           if (driver.use5_4) {
              if (/*lefrInputAntennaCbk*/ 1) { /* write warning only if cbk is set */
                if (driver.inputAntennaWarnings++ < driver.lefrInputAntennaWarnings) {
                   driver.outMsg = (char*)lefMalloc(10000);
                   sprintf (driver.outMsg,
                      "INPUTPINANTENNASIZE statement is a version 5.3 or earlier syntax.\nYour lef file with version %g, has both old and new INPUTPINANTENNASIZE syntax, which is illegal.", driver.versionNum);
                   driver.lefError(1671, driver.outMsg);
                   lefFree(driver.outMsg);
                   /*CHKERR();*/
                }
              }
           }
        }
        if (/*lefrInputAntennaCbk*/ 1)
          driver.lefrInputAntennaCbk( (yystack_[1].value.doubleVal));
    }
#line 12546 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1419: // output_antenna: "OUTPUTPINANTENNASIZE" NUMBER ';'
#line 7587 "LefParser.yy"
    { /* 5.3 syntax */
        driver.use5_3 = 1;
        if (driver.ignoreVersion) {
           /* do nothing */
        } else if (driver.versionNum > 5.3) {
           /* A 5.3 syntax in 5.4 */
           if (driver.use5_4) {
              if (/*lefrOutputAntennaCbk*/ 1) { /* write warning only if cbk is set */
                if (driver.outputAntennaWarnings++ < driver.lefrOutputAntennaWarnings) {
                   driver.outMsg = (char*)lefMalloc(10000);
                   sprintf (driver.outMsg,
                      "OUTPUTPINANTENNASIZE statement is a version 5.3 or earlier syntax.\nYour lef file with version %g, has both old and new OUTPUTPINANTENNASIZE syntax, which is illegal.", driver.versionNum);
                   driver.lefError(1672, driver.outMsg);
                   lefFree(driver.outMsg);
                   /*CHKERR();*/
                }
              }
           }
        }
        if (/*lefrOutputAntennaCbk*/ 1)
          driver.lefrOutputAntennaCbk( (yystack_[1].value.doubleVal));
    }
#line 12573 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1420: // inout_antenna: "INOUTPINANTENNASIZE" NUMBER ';'
#line 7611 "LefParser.yy"
    { /* 5.3 syntax */
        driver.use5_3 = 1;
        if (driver.ignoreVersion) {
           /* do nothing */
        } else if (driver.versionNum > 5.3) {
           /* A 5.3 syntax in 5.4 */
           if (driver.use5_4) {
              if (/*lefrInoutAntennaCbk*/ 1) { /* write warning only if cbk is set */
                if (driver.inoutAntennaWarnings++ < driver.lefrInoutAntennaWarnings) {
                   driver.outMsg = (char*)lefMalloc(10000);
                   sprintf (driver.outMsg,
                      "INOUTPINANTENNASIZE statement is a version 5.3 or earlier syntax.\nYour lef file with version %g, has both old and new INOUTPINANTENNASIZE syntax, which is illegal.", driver.versionNum);
                   driver.lefError(1673, driver.outMsg);
                   lefFree(driver.outMsg);
                   /*CHKERR();*/
                }
              }
           }
        }
        if (/*lefrInoutAntennaCbk*/ 1)
          driver.lefrInoutAntennaCbk( (yystack_[1].value.doubleVal));
    }
#line 12600 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1421: // antenna_input: "ANTENNAINPUTGATEAREA" NUMBER ';'
#line 7635 "LefParser.yy"
    { /* 5.4 syntax */
        /* 11/12/2002 - this is obsolete in 5.5, suppose should be ingored */
        /* 12/16/2002 - talked to Dave Noice, leave them in here for debugging*/
        driver.use5_4 = 1;
        if (driver.ignoreVersion) {
           /* do nothing */
        } else if (driver.versionNum < 5.4) {
           if (/*lefrAntennaInputCbk*/ 1) { /* write warning only if cbk is set */
             if (driver.antennaInputWarnings++ < driver.lefrAntennaInputWarnings) {
               driver.outMsg = (char*)lefMalloc(10000);
               sprintf (driver.outMsg,
                  "ANTENNAINPUTGATEAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.\nEither update your VERSION number or use the 5.3 syntax.", driver.versionNum);
               driver.lefError(1674, driver.outMsg);
               lefFree(driver.outMsg);
               /*CHKERR();*/
             }
           }
        } else if (driver.use5_3) {
           if (/*lefrAntennaInputCbk*/ 1) { /* write warning only if cbk is set */
             if (driver.antennaInputWarnings++ < driver.lefrAntennaInputWarnings) {
                driver.outMsg = (char*)lefMalloc(10000);
                sprintf (driver.outMsg,
                   "ANTENNAINPUTGATEAREA statement is a version 5.4 or later syntax.\nYour lef file with version %g, has both old and new ANTENNAINPUTGATEAREA syntax, which is illegal.", driver.versionNum);
                driver.lefError(1675, driver.outMsg);
                lefFree(driver.outMsg);
               /*CHKERR();*/
             }
           }
        }
        if (/*lefrAntennaInputCbk*/ 1)
          driver.lefrAntennaInputCbk( (yystack_[1].value.doubleVal));
    }
#line 12637 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1422: // antenna_inout: "ANTENNAINOUTDIFFAREA" NUMBER ';'
#line 7669 "LefParser.yy"
    { /* 5.4 syntax */
        /* 11/12/2002 - this is obsolete in 5.5, & will be ignored */
        /* 12/16/2002 - talked to Dave Noice, leave them in here for debugging*/
        driver.use5_4 = 1;
        if (driver.ignoreVersion) {
           /* do nothing */
        } else if (driver.versionNum < 5.4) {
           if (/*lefrAntennaInoutCbk*/ 1) { /* write warning only if cbk is set */
              if (driver.antennaInoutWarnings++ < driver.lefrAntennaInoutWarnings) {
                 driver.outMsg = (char*)lefMalloc(10000);
                 sprintf (driver.outMsg,
                    "ANTENNAINOUTDIFFAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.\nEither update your VERSION number or use the 5.3 syntax.", driver.versionNum);
                 driver.lefError(1676, driver.outMsg);
                 lefFree(driver.outMsg);
                 /*CHKERR();*/
              }
           }
        } else if (driver.use5_3) {
           if (/*lefrAntennaInoutCbk*/ 1) { /* write warning only if cbk is set */
              if (driver.antennaInoutWarnings++ < driver.lefrAntennaInoutWarnings) {
                 driver.outMsg = (char*)lefMalloc(10000);
                 sprintf (driver.outMsg,
                    "ANTENNAINOUTDIFFAREA statement is a version 5.4 or later syntax.\nYour lef file with version %g, has both old and new ANTENNAINOUTDIFFAREA syntax, which is illegal.", driver.versionNum);
                 driver.lefError(1677, driver.outMsg);
                 lefFree(driver.outMsg);
                 /*CHKERR();*/
              }
           }
        }
        if (/*lefrAntennaInoutCbk*/ 1)
          driver.lefrAntennaInoutCbk( (yystack_[1].value.doubleVal));
    }
#line 12674 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1423: // antenna_output: "ANTENNAOUTPUTDIFFAREA" NUMBER ';'
#line 7703 "LefParser.yy"
    { /* 5.4 syntax */
        /* 11/12/2002 - this is obsolete in 5.5, & will be ignored */
        /* 12/16/2002 - talked to Dave Noice, leave them in here for debugging*/
        driver.use5_4 = 1;
        if (driver.ignoreVersion) {
           /* do nothing */
        } else if (driver.versionNum < 5.4) {
           if (/*lefrAntennaOutputCbk*/ 1) { /* write warning only if cbk is set */
              if (driver.antennaOutputWarnings++ < driver.lefrAntennaOutputWarnings) {
                 driver.outMsg = (char*)lefMalloc(10000);
                 sprintf (driver.outMsg,
                    "ANTENNAOUTPUTDIFFAREA statement is a version 5.4 and later syntax.\nYour lef file is defined with version %g.\nEither update your VERSION number or use the 5.3 syntax.", driver.versionNum);
                 driver.lefError(1678, driver.outMsg);
                 lefFree(driver.outMsg);
                 /*CHKERR();*/
              }
           }
        } else if (driver.use5_3) {
           if (/*lefrAntennaOutputCbk*/ 1) { /* write warning only if cbk is set */
              if (driver.antennaOutputWarnings++ < driver.lefrAntennaOutputWarnings) {
                 driver.outMsg = (char*)lefMalloc(10000);
                 sprintf (driver.outMsg,
                    "ANTENNAOUTPUTDIFFAREA statement is a version 5.4 or later syntax.\nYour lef file with version %g, has both old and new ANTENNAOUTPUTDIFFAREA syntax, which is illegal.", driver.versionNum);
                 driver.lefError(1679, driver.outMsg);
                 lefFree(driver.outMsg);
                 /*CHKERR();*/
              }
           }
        }
        if (/*lefrAntennaOutputCbk*/ 1)
          driver.lefrAntennaOutputCbk( (yystack_[1].value.doubleVal));
    }
#line 12711 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;

  case 1426: // extension: "BEGINEXT"
#line 7740 "LefParser.yy"
    { 
        if (/*lefrExtensionCbk*/ 1)
          driver.lefrExtensionCbk( driver.Hist_text);
        if (driver.versionNum >= 5.6)
           driver.ge56almostDone = 1;
    }
#line 12722 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"
    break;


#line 12726 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"

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


  const short Parser::yypact_ninf_ = -1821;

  const short Parser::yytable_ninf_ = -1142;

  const short
  Parser::yypact_[] =
  {
    1299, -1821,   147,  2842, -1821, -1821, -1821, -1821, -1821, -1821,
      32,    32,    32, -1821, -1821,  -261,    -7,    32,    32,  -256,
      32,    32,    32,    32,    32, -1821, -1821, -1821,    32,    32,
      32,    68, -1821,    32,   145,    32, -1821, -1821, -1821, -1821,
      32,  -104,  -240, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
    -225, -1821,    20, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821,    84, -1821,  -186,  -182,  -175,  -168,
    -144, -1821, -1821,  -128,  -100,   -98,   -82,   -64, -1821, -1821,
     -44,   -27,    10,   127,   138,   140,   163,   170,   187,  2428,
     -80,   192,  -164,    32,   195,   210,    57, -1821,   229,   245,
     259, -1821,  2428,    32, -1821, -1821,   234,   128,    -5,   141,
     -19,  3081,    49,   -42,   329, -1821,   506, -1821, -1821,   -30,
     -25,   153,  1582,   126,   340, -1821,   152,   267,   273,   276,
     278, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821,   289, -1821, -1821, -1821,   296, -1821, -1821, -1821,
      -9, -1821,   302, -1821, -1821,   319,   332,   335, -1821,   493,
     569,   582,   432,   609,   608,   585,   599,   438, -1821,    39,
      32,   427,    32,    32,    32, -1821, -1821,   418,    32, -1821,
   -1821,    32,    32, -1821,    32,    32,    32, -1821,   774,    32,
      67,   807,    32,    32,    32,    32,    32,    32,   -20,    32,
     142, -1821,    32,    32,    32,    32,    32,    32,    32,    32,
      32,    32, -1821,    32,    32,    32,    32,   142, -1821,    32,
    -193,    32,    32,    32,    32,    32,    32,    32,   612,    32,
       2,    32,    32, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
      32, -1821, -1821, -1821,   748, -1821, -1821, -1821,   650, -1821,
   -1821, -1821, -1821,   464, -1821,  -154,   443,   541, -1821, -1821,
   -1821,   469,   704,   471, -1821, -1821, -1821,   -48, -1821, -1821,
      32, -1821, -1821, -1821, -1821, -1821, -1821,   -18, -1821,   712,
   -1821,    29,   484,   485, -1821, -1821,   132, -1821,    32, -1821,
   -1821,    32,   -33, -1821, -1821,   747,   784, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,  2428,
   -1821, -1821, -1821, -1821,   -45, -1821, -1821, -1821, -1821, -1821,
   -1821,    32,   495, -1821, -1821, -1821, -1821,   786, -1821, -1821,
   -1821, -1821,   104,   135,   -92,   -92,   -92,     9,   509, -1821,
     746,    12, -1821,   633, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821,    32,    32,    32, -1821,
      32,    32,    32,    32,    32, -1821, -1821, -1821,    32, -1821,
     132,   498,   499,   503,   106,   106, -1821,   511,   106,   106,
     513,   514,   672,   523, -1821,   525,   761,   162,   528,   530,
      32,   529,   -66,   534,    32,   535,   536,  -251,   539, -1821,
   -1821, -1821, -1821,   540,   542, -1821, -1821,    32,  2428,   553,
     561,   564,   565,   566,   573,   574, -1821,   604,   581,   132,
   -1821,   587,   164,   287,    32,   594,   823,   175,   589,   590,
     591,   598,   611,   620,   622, -1821,   630,   -55,   623, -1821,
   -1821, -1821, -1821, -1821, -1821,   635,   638,   639,   597,   653,
     654,   643,   661, -1821, -1821,   -50, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821,    32,    32,   655,   759,  -141,   506,   656,
   -1821,   533,   665,   963,  -177, -1821, -1821,   -81,    32,    32,
      32,   132,    32,    32,    32, -1821, -1821,   681, -1821, -1821,
     299,   670,   673,   678,  2428, -1821,  1032,  -192,    31, -1821,
      24, -1821,  1042,    -2, -1821,    34, -1821, -1821,   680,   694,
   -1821, -1821,   695,   684,   685,   698,   688,   701,   690,  2428,
     693,   706,  1056,   697,   699,   700,  -188, -1821,   703,   809,
   -1821, -1821,   308,   705,  3211, -1821, -1821,   -12,   114,   116,
     155,    32,  1003, -1821, -1821,  3305,  2428,  2428, -1821, -1821,
     168,   235,  2428, -1821, -1821,   401,   135,   104,   104, -1821,
     417, -1821,   135, -1821,   135, -1821,   135,    53,    11,   129,
     135,   -92,   326,   326,   326,    32, -1821,    21, -1821, -1821,
     716,   707, -1821,    32, -1821,    30, -1821, -1821,   721, -1821,
     722,   723,   728,   731,   732,   733,   734,   725,   726,   727,
     729,   735,   736,   737,   738,   739,  1000,  1046,   132, -1821,
   -1821, -1821,   730, -1821,   741,   742, -1821,   743,   745, -1821,
   -1821, -1821, -1821, -1821, -1821,   749, -1821,  1046, -1821, -1821,
     773,   751,   750, -1821,   753,    32, -1821,    32,    32,   755,
   -1821,   756, -1821, -1821, -1821,   757, -1821, -1821, -1821,    32,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,   785,    32,
   -1821, -1821, -1821, -1821, -1821, -1821,   762, -1821,   763,    32,
     563,  -259, -1821,    32,   760,   764, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821,    17,    32,   792,    32,   793,   767, -1821,
   -1821, -1821, -1821, -1821,   768,   583,  -212, -1821, -1821,   769,
     744, -1821,    32, -1821, -1821, -1821,   770, -1821,   132, -1821,
     -96, -1821,  -177, -1821, -1821, -1821,   771, -1821, -1821, -1821,
     772,   775,    32,   776,   777,   132,   778,  1133,   836,    32,
   -1821,    32, -1821, -1821, -1821, -1821, -1821,    79,    32, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,   779,
     794,    44,   781, -1821, -1821,   783, -1821,   787, -1821, -1821,
   -1821,   607,    82, -1821,    32, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821,   788,   796, -1821,    32, -1821,
      32,    32,    32,    32,    32,    32,    32, -1821,    32,    32,
      32,    32,    32,    41,    16, -1821,    32,    32,    32,    32,
   -1821, -1821,   797,    32, -1821,    32,    32, -1821, -1821, -1821,
   -1821, -1821,    32, -1821,    32,    32,    32,    32,    32,    32,
     -69,    38, -1821,   798,    32,    27,    32,    32,   157, -1821,
   -1821, -1821,    76,   801,   132, -1821,   132,   132,   132,   132,
     132, -1821,   790, -1821,   -68,   854,   795, -1821,    32,    32,
      32,    32, -1821, -1821, -1821, -1821, -1821, -1821, -1821,    32,
      32,    32,    32,   803,   804,   805, -1821, -1821,   897,    32,
     799, -1821,   160, -1821, -1821,   869,    72,   -28,    32,   802,
     808,    13,    32,    32,   810,    32,    32,   811,   812,   813,
   -1821, -1821,   816, -1821,     8, -1821,   494,   104,   104,   104,
     104, -1821, -1821, -1821,     8, -1821,    25,    45,    92,   104,
     104, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
     104,   135,   135,   135,   135, -1821,   -92,   -92,   -92,   -92,
       8,    87, -1821, -1821, -1821,   814,    32,  1142,  -248, -1821,
     815, -1821,   483,   817,    32,   986,  -248, -1821,   -56,   -56,
     -56,   -56,   -56,   -56,   -56,   -56, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821,    32, -1821, -1821,   819, -1821,
     132, -1821, -1821, -1821, -1821, -1821,   821,    32,  1141,   806,
   -1821,    32,   105, -1821,   827, -1821, -1821, -1821, -1821, -1821,
   -1821,    32,   835, -1821,   180,   -35, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821,   840,   841, -1821, -1821, -1821, -1821,
     663,    17,   902,    32, -1821,    32, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821,  1173, -1821,   842, -1821,   132,   132,   963,
   -1821,  -126, -1821,   849, -1821, -1821,   843, -1821, -1821,   846,
   -1821,    32,    32,  -231, -1821,   847, -1821, -1821,   850, -1821,
     851, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821,   852, -1821, -1821,   847,  1092,  1092,  1096,  1096,  1096,
    1092,  1092,   765,  1092,  1092,  1092,  1092,   855,   862,   863,
     864,   865,   867,  -235,  2428,   873,   879,   880,   881,   894,
      32,   895,   884,   898,   885,   887,   900,   905,    32,    32,
     889,   907,   896,   899,   906,   909,   910,   912, -1821, -1821,
   -1821, -1821,   913, -1821, -1821, -1821,   914,   917,   918,   919,
   -1821, -1821, -1821, -1821, -1821, -1821,   921,   922,   923,   271,
   -1821,  1281,   -18, -1821,   -68, -1821, -1821,   297,   132, -1821,
     132,   132,   929, -1821,   942, -1821, -1821, -1821, -1821,    32,
      32,    32,    32,   944,    32,    32,    32,    32,   933,   934,
     935,    32, -1821,   323,    32,   347, -1821,   948, -1821, -1821,
   -1821,   938,    32, -1821, -1821, -1821,   159, -1821,    32, -1821,
     947, -1821, -1821,  1231,    32, -1821,  1232,  1258, -1821,  1259,
    1260, -1821, -1821, -1821,  2428,  2428, -1821, -1821,  1016,   104,
   -1821,   357,   357, -1821, -1821,   135, -1821, -1821,   621,   621,
     621,  -151,   387, -1821, -1821,   952,   952, -1821,   952,   -92,
   -1821, -1821,   325, -1821, -1821, -1821, -1821,   961, -1821,   962,
     965,   675, -1821, -1821,  1261, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821,   327, -1821, -1821, -1821,   966, -1821,   968,
     958,   959,   960,   969,   970,   971,   972,   973, -1821, -1821,
   -1821,   132, -1821, -1821,    32,    32,   399, -1821,    32,   974,
   -1821,   975, -1821, -1821,   -94, -1821, -1821,   987, -1821, -1821,
   -1821,    32,    32, -1821, -1821, -1821, -1821, -1821,   995,    32,
      32,   409, -1821,  1129,    32,    32,    32,    32, -1821,    32,
      32,    32,   977, -1821,    32,  1005, -1821,  1122,    32, -1821,
     132,   984, -1821, -1821, -1821, -1821, -1821, -1821,   985,   991,
     992, -1821, -1821, -1821, -1821, -1821,  1302, -1821,   993,   997,
   -1821,   998,   999,  1001,  1004,  1006, -1821, -1821, -1821, -1821,
    1008,  1012,  1013,  1014,  1015, -1821, -1821, -1821, -1821, -1821,
   -1821,  1018, -1821, -1821, -1821, -1821, -1821, -1821,  1020,    32,
    1021, -1821,  1022, -1821, -1821,  1026,  1027,    32,    32, -1821,
     679,   183, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821,  1028, -1821, -1821, -1821, -1821, -1821, -1821,  1296, -1821,
    1029,  1030,  1300,  1297,   -34, -1821,   132, -1821,  1318, -1821,
    1017,  1031,  1034,  1035,  1036,  1037, -1821,   203,  1038,  1039,
    1040,  1041, -1821, -1821, -1821,  1277, -1821, -1821,    32, -1821,
   -1821, -1821,   242, -1821,    32, -1821, -1821, -1821,  1104,    32,
     160, -1821,  1067,    32,    32,    32,    32,  1043,  1045, -1821,
     -36,    -3,    37,    81,   -71,  1082, -1821,  1060,    32,  -107,
    1049,  1050, -1821, -1821, -1821,  1063, -1821,  1089, -1821,  1189,
   -1821,  1189, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821,   167, -1821,  1431,    32,  1061, -1821,  1058, -1821,
   -1821,   331,   334,    32,    32,  1062,    32,   702, -1821, -1821,
      32,  1064,  1068,    32,  1069, -1821, -1821, -1821,  1093, -1821,
    1097,  1248,  1084,  1201, -1821,    32, -1821, -1821,    32,    32,
      32,    32,    32,   132, -1821, -1821, -1821, -1821, -1821,  1202,
      32,  1087, -1821, -1821,  1099, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821,  1077, -1821, -1821,
   -1821, -1821,  1085,  1088, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821,  1090, -1821, -1821, -1821, -1821,   -16,    32,  1091, -1821,
     304,  1098,  1318, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821,    32,    32, -1821, -1821,    32,    32,
   -1821,  1161,    32,  1172,  1179,  1180,  1181, -1821, -1821, -1821,
    1410,   104,   135,   -92,    32,  1082, -1821, -1821,  1105, -1821,
    -118, -1821, -1821, -1821,    32,  1089, -1821,    32,    32,    32,
     342,  1106, -1821, -1821, -1821,    32,  1102, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,  1347, -1821,
   -1821, -1821,  1107, -1821,    32,    32,    32, -1821,    32, -1821,
    1136,  1197,  1140, -1821, -1821,  1114, -1821,   132, -1821,  1476,
   -1821, -1821, -1821, -1821, -1821, -1821,    32,    32,  1117,  1478,
   -1821,   -34, -1821, -1821,  1119,   -38,  1120,  1121,    32,    32,
      32,  1123,    32,    32,    32,    32, -1821,   621,    19,   952,
    1124, -1821,  1152, -1821, -1821, -1821,  1126, -1821,    32, -1821,
   -1821, -1821, -1821, -1821,   173,  1214, -1821,  1200,  1200,  1429,
      32,  1131, -1821, -1821, -1821,  1319,  1216,    32,    -6,    32,
    1135,  1136,    32,    32, -1821,  1378,   306,   132,  1202,    32,
   -1821, -1821, -1821,    32,  1137, -1821,  1138, -1821, -1821,    32,
      32,  1183,    32, -1821, -1821, -1821,  1139,  1144, -1821,    32,
   -1821, -1821,   -39, -1821, -1821,    32,    32, -1821, -1821,    32,
   -1821,  1171, -1821,  1448,  1311,    32, -1821, -1821, -1821,    32,
      32, -1821, -1821,  1238, -1821, -1821, -1821, -1821,    32, -1821,
   -1821,    32, -1821,  1233,  1234, -1821, -1821,  1155,    32,    32,
      32,  1399,  1400,  1503,  1159,  1199,    32,    32,  1237, -1821,
   -1821, -1821,    32, -1821, -1821,    32, -1821, -1821,  1298,    32,
      32, -1821, -1821,  1195,  1182,  1169,    32,    32, -1821,    32,
      32,  1175, -1821, -1821,    32, -1821,    32, -1821, -1821,    32,
    1543, -1821, -1821,  1198,    32, -1821, -1821,  1203,    32,    32,
   -1821,  1195,  1185, -1821,    32,    32,    32,  1177, -1821, -1821,
   -1821, -1821,   351, -1821, -1821,   353, -1821, -1821, -1821,    32,
    1178,  1543,   363,    32,    32,    32, -1821, -1821,  1122, -1821,
    1191, -1821, -1821,    32, -1821,  1192,  1192, -1821, -1821,   -49,
   -1821, -1821,   395,   398, -1821,  1267, -1821, -1821, -1821, -1821,
   -1821,  1208,  1210,    32,  1184,    32, -1821, -1821,  1527,    32,
      32, -1821,  1294,    32,    32, -1821, -1821, -1821,   407,    32,
    1240,    32,    32, -1821, -1821,  1524, -1821,    32,  1190,  1204,
    1207,    32,  1209,  1211, -1821, -1821, -1821, -1821,  1286, -1821,
   -1821,    32,  1212, -1821, -1821, -1821,  1215, -1821, -1821,    32,
   -1821,    32,    32, -1821, -1821, -1821,    32,  1217,  1507,    32,
      32,    32,    32,  1218, -1821, -1821,   436,    32,    32, -1821,
      32, -1821,    32,    32,  1219,    32,    32,  1221,  1222,  1223,
      32, -1821, -1821, -1821,  1224, -1821
  };

  const short
  Parser::yydefact_[] =
  {
       0,   425,     0,  1424,     1,  1258,  1262,  1260,  1254,  1256,
       0,     0,     0,  1225,  1426,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   828,   492,   909,     0,     0,
       0,     0,   839,     0,     0,     0,  1307,   880,   820,   480,
       0,     0,   417,   717,   780,   428,   436,   429,   424,   430,
     437,   461,   462,   463,   431,   481,   432,   496,   464,   433,
       0,   715,     0,   434,   435,   439,   822,   442,   830,   441,
     440,   447,   443,   884,   444,   913,   445,  1229,   438,   465,
     446,   448,   449,   450,   452,   451,   453,   454,   455,   456,
     457,   458,   459,   460,   426,  1425,     0,     0,     0,     0,
       0,     2,     3,     0,     0,     0,     0,     0,   474,   473,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1384,     0,     0,
       0,  1309,     0,     0,   475,   476,     0,     0,     0,     0,
       0,     0,   718,     0,   785,   804,   801,   787,   806,     0,
       0,     0,   907,  1223,     0,   416,     0,     0,     0,     0,
       0,  1422,  1421,  1423,  1226,   422,   472,  1401,   838,   421,
    1383,  1381,  1382,  1420,  1418,     5,     6,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    56,    57,    58,    59,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    69,
      70,    71,    72,    73,    74,    75,    76,    77,    78,    79,
      80,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,   108,   109,
     110,   111,   112,   113,   114,   115,   116,   117,   118,   119,
     120,   121,   122,   123,   124,   125,   126,   127,   128,   129,
     130,   131,   132,   133,   134,   135,   136,   137,   138,   139,
     140,   141,   142,   143,   144,   145,   146,   147,   148,   149,
     150,   151,   152,   153,   154,   155,   156,   157,   158,   159,
     160,   161,   162,   163,   164,   165,   166,   167,   168,   169,
     170,   171,   172,   173,   174,   175,   176,   177,   178,   179,
     180,   181,   182,   183,   184,   185,   186,   187,   188,   189,
     190,   191,   192,   193,   194,   195,   196,   197,   198,   199,
     200,   201,   202,   203,   204,   205,   206,   207,   208,   209,
     210,   211,   212,   213,   214,   215,   216,   217,   218,   219,
     220,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,   235,   236,   237,   238,   239,
     240,   241,   242,   243,   244,   245,   246,   247,   248,   249,
     250,   251,   252,   253,   254,   255,   256,   257,   258,   259,
     260,   261,   262,   263,   264,   265,   266,   267,   268,   269,
     270,   271,   272,   273,   274,   275,   276,   277,   278,   279,
     280,   281,   282,   283,   284,   285,   286,   287,   288,   289,
     290,   291,   292,   293,   294,   295,   296,   297,   298,   299,
     300,   301,   302,   303,   304,   305,   306,   307,   308,   309,
     310,   311,   312,   313,   314,   315,   316,   317,   318,   319,
     320,   321,   322,   323,   324,   325,   326,   327,   328,   329,
     330,   331,   332,   333,   334,   335,   336,   337,   338,   339,
     340,   341,   342,   343,   344,   345,   346,   347,   348,   349,
     350,   351,   352,   353,   354,   355,   356,   357,   358,   359,
     360,   361,   362,   363,   364,   365,   366,   367,   368,   369,
     370,   371,   372,   373,   374,   375,   376,   377,   378,   379,
     380,   381,   382,   383,   384,   385,   386,   387,   388,   389,
     390,   391,   392,   393,   394,   395,   396,   397,   398,   399,
     400,   401,   402,   403,   404,   405,   406,   407,   408,   409,
     410,   411,   412,   413,   414,   415,     4,   493,   910,   470,
     713,   712,     0,   467,   466,   840,     0,   469,   468,  1419,
       0,   881,     0,   478,   477,     0,     0,     0,   781,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   482,     0,
       0,     0,     0,     0,     0,   547,   555,     0,     0,   544,
     552,     0,     0,   559,     0,     0,     0,   498,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     619,   494,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   597,     0,     0,     0,     0,   619,   531,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   589,     0,
       0,     0,     0,   491,   497,   719,   720,   752,   771,   740,
       0,   742,   721,   732,     0,   733,   735,   737,     0,   738,
     773,   786,   783,     0,   788,   789,   800,     0,   826,   819,
     823,     0,     0,     0,   827,   831,   833,     0,   882,   899,
       0,   894,   879,   885,   888,   887,   889,     0,  1159,     0,
     992,     0,     0,     0,   994,  1149,     0,  1002,     0,   934,
     998,     0,     0,   939,  1162,     0,     0,   914,   919,   915,
     916,   917,   918,   921,   920,   922,   923,   924,   926,     0,
     925,   927,  1006,   930,     0,   931,   932,   933,  1164,  1235,
    1233,     0,     0,  1240,   999,  1237,  1231,     0,  1230,  1245,
    1243,   427,     0,     0,     0,     0,     0,     0,     0,   837,
     845,     0,  1313,     0,  1323,  1311,  1317,  1325,  1315,  1319,
    1321,  1310,  1380,   471,   418,   420,     0,     0,     0,   479,
       0,     0,     0,     0,     0,   654,   653,   655,   533,   549,
       0,     0,     0,     0,     0,     0,   561,     0,     0,     0,
       0,     0,     0,     0,   557,     0,   599,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   691,
     692,   689,   690,     0,     0,   620,   621,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   582,     0,     0,     0,
     586,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   512,     0,     0,     0,   684,
     688,   686,   685,   683,   687,     0,     0,     0,     0,     0,
       0,     0,     0,   778,   716,   734,   757,   761,   758,   760,
     759,   768,   766,   769,   767,   754,   762,   763,   764,   765,
     756,   755,   748,     0,     0,     0,     0,   770,     0,     0,
     793,     0,     0,     0,   790,   791,   802,     0,     0,     0,
       0,     0,     0,     0,     0,   807,   821,     0,   829,   836,
       0,     0,     0,     0,     0,   901,     0,     0,   951,   955,
     959,   960,   948,     0,   954,   956,   950,   957,     0,     0,
    1153,  1151,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   988,     0,     0,
     911,   908,     0,     0,     0,  1148,  1105,     0,     0,     0,
       0,     0,     0,  1116,  1115,     0,     0,     0,  1251,  1244,
       0,     0,     0,  1227,  1224,     0,     0,     0,     0,  1278,
    1264,  1293,     0,  1292,     0,  1297,     0,     0,  1264,     0,
       0,     0,  1264,  1264,  1264,     0,  1408,     0,  1404,  1407,
       0,     0,   847,     0,  1391,     0,  1387,  1390,     0,  1308,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   651,     0,   563,
     543,   546,     0,   645,     0,     0,   554,     0,     0,   562,
     527,   695,   696,   697,   698,     0,   551,   651,   511,   600,
     601,     0,     0,   525,     0,     0,   528,     0,     0,     0,
     578,     0,   596,   509,   504,     0,   508,   516,   526,     0,
     495,   579,   521,   571,   572,   573,   576,   580,   693,     0,
     575,  1123,  1125,   637,   581,   506,     0,   502,     0,     0,
       0,     0,   669,     0,     0,     0,   595,   524,   570,   569,
     568,   567,   574,  1350,     0,     0,     0,     0,     0,   523,
     501,   510,   522,   753,     0,     0,     0,   743,   739,     0,
       0,   736,     0,   905,   751,   749,     0,   776,     0,   774,
     789,   805,     0,   797,   798,   799,     0,   817,   782,   792,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     832,     0,   834,   891,   890,   892,   883,     0,     0,   898,
     896,   897,   893,   895,   952,   953,   973,   969,   972,   970,
     971,   974,   961,   949,   979,   980,   976,   975,   977,   978,
     962,   968,   965,   963,   964,   966,   967,   958,   947,     0,
       0,     0,     0,   928,   929,     0,   981,     0,   987,  1003,
     986,     0,     0,   936,     0,   985,   984,   983,   943,   941,
     942,   938,   940,   991,   989,     0,     0,   996,     0,   997,
       0,     0,     0,     0,     0,     0,     0,  1065,     0,     0,
       0,     0,     0,     0,     0,  1004,     0,     0,     0,     0,
    1071,  1035,     0,     0,  1014,     0,     0,  1029,  1067,  1031,
    1033,  1087,     0,  1048,     0,     0,     0,     0,     0,     0,
       0,  1100,  1069,     0,     0,     0,     0,     0,     0,  1001,
    1007,  1017,  1088,     0,     0,  1125,     0,     0,     0,     0,
       0,  1127,     0,  1147,  1104,     0,     0,  1215,     0,     0,
       0,     0,  1166,  1195,  1191,  1192,  1193,  1196,  1214,     0,
       0,     0,     0,     0,     0,     0,  1194,  1197,     0,     0,
       0,  1168,     0,  1161,  1165,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1249,  1247,     0,  1246,     0,  1275,     0,     0,     0,     0,
       0,  1265,  1266,  1259,     0,  1289,     0,     0,     0,     0,
       0,  1302,  1303,  1298,  1299,  1300,  1301,  1304,  1305,  1306,
       0,     0,     0,     0,     0,  1263,     0,     0,     0,     0,
       0,     0,  1261,  1255,  1257,     0,     0,     0,  1264,  1405,
       0,   846,   841,     0,     0,     0,  1264,  1388,     0,     0,
       0,     0,     0,     0,     0,     0,   484,   487,   489,   490,
     486,   485,   483,   488,   535,     0,   534,   652,     0,   564,
       0,   548,   556,   545,   553,   560,     0,     0,     0,     0,
     519,     0,     0,   674,     0,   537,   540,   536,   577,   505,
     591,     0,     0,   584,     0,     0,   507,   503,   593,   671,
     672,   673,   532,   670,     0,     0,   517,  1352,  1355,  1353,
    1356,  1350,     0,     0,   606,     0,   590,   772,   747,   746,
     745,   741,   744,     0,   779,     0,   750,     0,     0,     0,
     795,     0,   803,     0,   808,   809,     0,   816,   815,     0,
     814,     0,     0,     0,   835,     0,   900,   902,     0,  1160,
       0,  1150,  1152,   993,   982,   995,   946,   945,   944,   935,
     937,     0,   990,   912,     0,  1375,  1375,     0,     0,     0,
    1375,  1375,     0,  1375,  1375,  1375,  1375,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1096,  1097,
    1098,  1099,     0,  1101,  1103,  1102,     0,     0,     0,     0,
    1091,  1094,  1095,  1093,  1092,  1090,     0,     0,     0,     0,
    1008,     0,     0,  1051,     0,  1106,  1125,     0,     0,  1124,
       0,     0,     0,  1129,     0,  1108,  1117,  1163,  1190,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1209,     0,     0,     0,  1206,     0,  1216,  1217,
    1218,     0,     0,  1198,  1200,  1199,     0,  1204,     0,  1205,
       0,  1236,  1234,     0,     0,  1252,     0,     0,  1241,     0,
       0,  1238,  1232,  1228,     0,     0,  1239,  1267,     0,     0,
    1276,  1272,  1271,  1273,  1274,     0,  1290,  1295,  1280,  1281,
    1279,  1287,  1288,  1285,  1286,  1283,  1284,  1294,  1282,     0,
    1406,  1410,     0,  1403,  1402,   714,   866,     0,   859,     0,
       0,     0,   850,   851,     0,   848,   853,   854,   855,   852,
     849,  1389,  1393,     0,  1386,  1385,  1327,     0,  1329,  1331,
       0,     0,     0,     0,     0,     0,     0,     0,   656,   550,
     648,     0,   558,   602,     0,     0,     0,   677,     0,     0,
     675,     0,   667,   667,   622,   694,   583,   630,   598,  1126,
     642,     0,     0,   643,   644,   587,   638,   639,   628,     0,
       0,     0,   680,  1344,     0,     0,     0,     0,  1357,     0,
       0,     0,     0,  1351,     0,     0,   667,   707,     0,   906,
       0,     0,   784,   794,   796,   818,   810,   812,     0,     0,
       0,   824,   903,   886,  1154,  1000,  1133,  1376,     0,     0,
    1378,     0,     0,     0,     0,     0,  1073,  1074,  1075,  1076,
       0,     0,     0,     0,     0,  1037,  1042,  1043,  1086,  1085,
    1082,     0,  1083,  1005,  1024,  1056,  1022,  1045,     0,     0,
       0,  1020,     0,  1038,  1039,     0,     0,     0,     0,  1016,
       0,     0,  1077,  1041,  1040,  1023,  1055,  1021,  1044,  1019,
    1028,     0,  1047,  1027,  1018,  1026,  1025,  1011,     0,  1009,
       0,     0,     0,  1118,     0,  1109,     0,  1125,     0,  1111,
       0,     0,     0,     0,     0,     0,  1219,     0,     0,     0,
       0,     0,  1189,  1188,  1187,     0,  1173,  1210,     0,  1174,
    1207,  1221,     0,  1183,     0,  1201,  1203,  1202,     0,     0,
       0,  1242,     0,     0,     0,     0,     0,     0,     0,  1268,
       0,     0,     0,     0,     0,     0,  1411,     0,     0,     0,
       0,     0,   863,   864,   865,   843,   842,     0,  1394,  1346,
    1333,  1346,  1332,  1314,  1324,  1312,  1318,  1326,  1316,  1320,
    1322,   667,     0,   646,     0,     0,     0,   678,     0,   529,
     530,     0,     0,     0,     0,     0,     0,   632,   640,   641,
       0,     0,     0,     0,     0,   681,  1345,  1354,     0,  1362,
       0,     0,     0,  1365,  1369,     0,   513,   514,     0,   607,
       0,     0,     0,     0,   775,   813,   811,   825,   904,     0,
       0,     0,  1060,  1061,     0,  1062,  1064,  1063,  1053,  1054,
    1066,  1059,  1057,  1058,  1052,  1084,  1072,     0,  1046,  1015,
    1030,  1068,     0,     0,  1081,  1080,  1079,  1049,  1078,  1070,
    1012,     0,  1010,  1089,  1050,  1119,  1120,     0,     0,  1125,
       0,     0,     0,  1128,  1178,  1176,  1180,  1182,  1220,  1167,
    1177,  1175,  1179,  1181,     0,     0,  1222,  1169,     0,     0,
    1170,     0,     0,     0,     0,     0,     0,  1250,  1248,  1269,
       0,     0,     0,     0,     0,  1409,  1412,   867,     0,   861,
       0,   856,   857,   844,     0,  1392,  1395,     0,  1348,  1348,
       0,     0,   649,   650,   648,     0,     0,   520,   676,   538,
     668,   541,   625,   623,   592,   631,   633,   634,   635,   629,
     594,   588,     0,   518,     0,     0,     0,  1358,     0,  1363,
     616,     0,     0,   708,   611,     0,  1125,     0,  1156,     0,
    1377,  1379,  1036,  1032,  1034,  1013,     0,     0,     0,     0,
    1110,     0,  1113,  1112,     0,     0,     0,     0,     0,  1211,
       0,     0,     0,     0,     0,     0,  1270,  1277,  1291,  1296,
       0,  1413,     0,   858,   860,   862,     0,  1396,     0,  1349,
    1328,  1330,   657,   565,     0,     0,   679,     0,     0,   626,
       0,     0,   682,  1360,  1367,     0,  1371,     0,  1334,     0,
       0,   616,     0,     0,   667,     0,     0,     0,  1155,     0,
    1122,  1121,  1107,     0,     0,  1130,     0,  1208,  1186,     0,
       0,     0,     0,  1253,  1145,  1146,  1134,  1137,  1414,     0,
    1397,  1347,   660,   566,   647,     0,     0,   539,   542,     0,
     624,     0,   585,  1373,  1339,     0,  1372,  1359,  1366,     0,
       0,  1335,  1364,     0,   515,   617,   614,   608,   612,   722,
     777,     0,  1157,     0,     0,  1114,  1184,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   499,
     665,   627,     0,  1374,  1361,     0,  1368,  1370,  1336,     0,
       0,   709,   667,   703,     0,     0,     0,     0,  1185,     0,
       0,     0,  1135,  1138,     0,   868,     0,   661,   663,     0,
     603,   667,   636,     0,     0,  1338,   618,   615,   609,     0,
     613,   703,     0,  1158,     0,     0,  1212,     0,  1172,  1142,
    1142,  1416,     0,   871,  1399,     0,   667,   667,   658,     0,
       0,   603,     0,     0,     0,     0,   710,   699,   707,   704,
       0,  1132,  1131,     0,  1171,  1136,  1139,  1415,  1417,     0,
    1398,  1400,     0,     0,   667,     0,   500,   604,   666,  1340,
    1337,     0,   610,     0,     0,     0,  1144,  1143,     0,     0,
       0,   869,     0,     0,     0,   872,   662,   664,     0,     0,
    1342,     0,     0,   700,   705,     0,  1213,     0,     0,     0,
       0,     0,     0,     0,   659,   605,  1343,  1341,     0,   701,
     667,     0,     0,   878,   877,   870,     0,   873,   874,     0,
     667,   706,     0,   876,   875,   711,   702,     0,     0,     0,
       0,     0,     0,     0,   723,   725,   724,     0,     0,   730,
       0,   726,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   728,   731,   727,     0,   729
  };

  const short
  Parser::yypgoto_[] =
  {
   -1821,   -10,   -84, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,  -695,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,  -553,
   -1821,   980, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821,   -41, -1821,  -450, -1821,   568, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821,  -525, -1821, -1726, -1821,   526,
   -1821,   196, -1821,   -86, -1821,  -110, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821,  -649, -1821, -1821,  -675, -1821, -1821, -1821,
   -1821,   243, -1821,   967, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821,   766, -1821, -1821,   508, -1821,
    -137, -1821,  -739, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821,   740,   496, -1821,
    -918, -1821, -1821,  -143, -1821, -1821, -1821, -1821, -1821, -1821,
     158, -1821,   248, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1820,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821,  -626, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,   431,
   -1821,    62, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,  -174, -1821,
   -1821, -1821, -1821, -1821, -1821,    69,   352, -1821, -1821, -1821,
   -1821, -1821,  -956, -1296, -1310, -1821, -1821, -1821, -1809,  -969,
   -1821, -1821, -1821, -1821, -1821,  -636, -1821, -1821, -1821, -1821,
   -1821, -1821,   454, -1821, -1821, -1821,  -512, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821,    43, -1821, -1821,  -456,  -217, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821,  1541, -1821, -1821, -1821,  -952, -1232,
   -1072,  -773,  -990,  -730, -1025, -1821, -1821, -1821, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821,  -373, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821,  -245,  -392,   207, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,  -662, -1821,
    -690, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821,
     658, -1821, -1821, -1821, -1821,  -386, -1821, -1821, -1821, -1821,
   -1821, -1821,   667, -1821, -1821, -1821, -1821,  -375, -1821, -1821,
   -1821, -1821, -1821, -1821, -1821, -1821, -1821, -1821
  };

  const short
  Parser::yydefgoto_[] =
  {
       0,   924,  1358,     2,    45,   137,   138,    46,    47,     3,
     155,    48,    49,    50,    51,    52,    53,   110,   136,   605,
      54,    55,   140,   618,    56,    57,   119,   683,   858,   141,
     684,   836,  2300,  1143,  2120,   875,  1066,  1752,  2177,  1753,
    2178,   828,   824,  1067,   829,   825,  1087,   832,  1740,  2223,
    1118,  1757,  1123,   887,  1754,  1768,   869,  1090,  1458,  2330,
    2331,  1148,  1786,  2122,  2282,  2337,  2194,  2283,  2281,  2190,
    2191,   857,  1965,  2179,  2230,  1971,  1967,  2108,  2181,  1475,
    1766,  1767,  1074,  2094,  1952,  2092,  1448,   818,  1446,  1951,
    2222,  2354,  2268,  2326,  2327,  2227,  2301,  1961,  1131,  1132,
    1462,  1463,  1746,  1747,  1771,  1772,   895,   853,  1472,  1085,
    2362,  2383,  2410,  2310,  2311,  2400,  1991,  2307,  2336,    58,
     788,    59,   143,    60,    61,   693,   902,  2284,  2425,  2426,
    2431,  2434,   694,   695,   905,   696,   900,  1156,  1157,   697,
     698,   898,   925,   699,   700,   899,   927,  1169,  1507,   904,
    1160,    62,   139,    63,    64,   928,   702,   146,   933,   934,
     935,  1172,  1511,   147,   705,   936,   148,   703,   706,   945,
    1178,  1513,    65,    66,   709,   149,   710,   711,    67,    68,
     714,   150,   715,   950,  1192,   716,    69,    70,    71,   126,
     790,  1714,  1936,  1042,  1422,  1715,  1716,  1717,  1718,  1719,
    1929,  2080,  1720,  1927,  2162,  2323,  2390,  2349,  2375,    72,
      73,   132,   722,   954,   151,   723,   724,   725,   957,  1203,
     726,   955,  1197,  1527,  1998,  1619,    74,   746,    75,   120,
     991,  1256,   152,   747,   981,  1242,   748,   986,  1252,  1243,
     749,   968,  1227,  1212,  1220,   750,   751,   752,   753,   754,
     755,   756,   972,   757,   977,   758,   759,   776,   760,   761,
     762,   979,  1309,  1564,   994,  1310,  1573,  1576,  1578,  1579,
    1570,  1581,  1552,  1577,  1597,  1569,  1820,  1851,  1852,  1311,
    1312,  1614,  1606,  1592,  1596,  1002,  1003,  1313,  1873,  1324,
    2036,  2138,  1122,  1759,  1474,  1004,  1624,  1880,  2038,   993,
    1367,  2261,  2319,  2262,  2320,  2345,  2367,  1364,   763,   764,
     765,  1231,   971,  1230,  1999,  2198,  2128,   766,   969,   767,
     768,  1353,  1005,  1354,  1633,  1647,  2149,  1355,  1356,  1656,
    1908,  1660,  1645,  1646,  1643,  2211,  1357,  1651,  1887,  1902,
      76,   777,    77,   106,  1014,  1369,   153,   778,  1012,  1007,
    1006,  1011,  1010,   779,  1015,  1373,  1675,  1674,  1361,  1665,
      78,    99,    79,   100,    80,    96,    98,    97,  1383,  1679,
    2071,  1027,  1028,  1029,  1400,    81,   131,   600,   801,  1051,
    1048,  1054,  1052,  1055,  1056,  1050,  1053,  1730,  1939,  1941,
    2242,  2276,  2380,  2397,  1977,  2088,  2170,  1490,  1491,  1773,
    1782,  1982,  2186,  2233,  2188,  2234,  2237,  2274,  1808,  2001,
    1811,  2004,    82,    83,    84,    85,    86,   596,  1426,  1045,
    1046,  1047,  1424,  1723,  2085,  2086,  2265,  2325,    87,   787,
    1418,  1037,  1038,  1039,  1416,  1702,  2075,  2076,  2263,  2322,
      88,    89,    90,    91,    92,    93,    94,    95
  };

  const short
  Parser::yytable_[] =
  {
     103,   104,   105,   704,  1409,  1617,   988,   111,   112,  1020,
     114,   115,   116,   117,   118,   745,  1179,  1389,   121,   122,
     123,  1620,  2368,   127,  1097,   130,  1374,  1962,  1401,   983,
     133,  1401,  1384,   951,  1385,   587,  1387,  1359,  1360,  1401,
    1410,  1315,  1317,  1368,   958,  1389,  1600,   959,   601,  1206,
    2069,  2266,   609,  1593,  1032,  1033,  1034,  1401,  1558,  2369,
    1989,  1214,  1215,   960,   961,  1401,  2370,   962,   610,  2041,
    1657,   792,   926,  1389,  2371,   611,  1405,  1030,   995,  2037,
    1412,  1413,  1414,   849,   850,  2069,  2136,  1180,  1221,   145,
    1487,   889,  1204,   707,  1963,  1327,   973,   815,   712,  1402,
     687,  1406,  1199,  1601,   612,  1588,  1248,   963,   687,  2079,
     978,   930,  1406,   592,   793,  1800,  1145,   996,  2069,   989,
    1602,   108,  1726,   602,  1653,  1035,   984,   877,  1043,     5,
       6,     7,   613,   590,   930,   688,  1663,  1035,  1831,  1167,
     996,   134,  1760,   688,  1417,   685,  1043,     4,   851,  1406,
    1612,  1658,  1685,  1425,  1761,   855,   107,  1130,  2069,  1168,
    1560,   113,   930,   101,   102,  2239,  1762,  1531,  1207,  1482,
    1708,  1314,   890,   135,  1727,  -419,   794,  1104,  1699,   144,
    1381,  1708,  1594,   795,  1382,   796,  1589,  1166,  1561,  1603,
    1562,   142,   930,  1832,  1068,   769,   770,  1801,   891,  1613,
     109,   856,  2069,   997,  1155,   145,  1222,   154,  1223,   686,
     974,  1905,   998,  1200,  1201,  1146,  1501,  1249,  1250,   952,
     964,   101,   102,   771,  1763,  2090,   997,   717,   797,   970,
     156,  1664,   999,  2144,   157,   998,  1202,    43,   689,   931,
    1251,   158,  1728,  1121,  1375,  1376,   689,  1390,   159,   965,
    1255,   690,  2372,  1386,  1510,   999,  1403,  1404,  1402,   690,
    2165,  1402,   931,   892,   591,   614,  1338,  1036,   798,  1402,
    1044,  1563,   160,  1016,   772,  1390,   718,  1021,  1590,  1036,
    1591,  1216,   615,  1217,   773,   966,  1488,  1402,  1044,  1224,
     931,  2240,   708,  2373,  1729,  1402,  1388,  1316,  1711,  1318,
     161,  1411,  1793,  1390,  1022,  1654,  1874,  1964,  1181,  1711,
    2164,  1604,   816,  1764,   124,  1185,   893,   125,  1225,   691,
     931,  1147,  1876,  1000,  1877,  1025,  2137,   691,   162,   713,
     163,  1407,  2204,  1001,   164,  1098,   588,  1031,  1320,  1489,
    1648,  1595,  1407,  1559,  1659,   616,  1000,   953,   101,   102,
     692,   985,   817,  1408,  1655,  2374,  1001,  1605,  1616,   839,
    1618,  2070,  2267,  1409,   165,  1208,  1218,  1677,  1219,   852,
    1205,  1121,  1121,  1319,  1321,  1226,  2241,   967,   617,  1407,
     101,   102,   799,   800,   166,  1209,  1210,  1377,  1378,  1379,
    1380,   128,  1906,  1765,   129,   923,  2070,   894,  1072,  1649,
    1650,   167,   101,   102,  1391,  1392,  1393,  1394,  1395,  1396,
     607,  1691,  1692,  1693,  1694,  1403,  1404,   923,  1403,  1404,
    1377,  1378,  1379,  1380,  1211,   701,  1403,  1404,  1397,  2070,
    1398,  1399,  1391,  1392,  1393,  1394,  1395,  1396,   168,  1381,
    1678,  1907,  1449,  1382,  1403,  1404,   101,   102,  1377,  1378,
    1379,  1380,  1403,  1404,  1091,  1680,  1397,   774,  1398,  1399,
    1391,  1392,  1393,  1394,  1395,  1396,  1704,  1134,  2248,  2070,
    1370,  1371,   719,   595,  1725,  1686,  1377,  1378,  1379,  1380,
     603,   101,   102,   604,  1397,   720,  1398,  1399,  1391,  1392,
    1393,  1394,  1395,  1396,  1993,  1525,   840,   775,  1241,  1391,
    1392,  1393,  1394,  1395,  1396,  1408,   721,  1526,  1609,  1023,
    1539,  1408,  1397,  2070,  1398,  1399,  1408,  1687,   101,   102,
     101,   102,  1687,  1397,  1372,  1398,  1399,  1017,   101,   102,
     101,   102,   781,  1018,  1461,  1749,  1391,  1392,  1393,  1394,
    1395,  1396,  1508,   923,   606,   923,   101,   102,  1024,   101,
     102,  1790,  1025,  1408,   937,   169,  2308,   608,  1017,  1519,
    1397,   923,  1398,  1399,  1026,   938,   170,  2040,   171,   101,
     102,   101,   102,  1362,  1363,  2332,   101,   102,   101,   102,
    2039,   101,   102,   782,   923,  1610,   923,   101,   102,   101,
     102,   172,  1125,  1794,   101,   102,   923,  2091,   173,  1850,
    2352,  2353,   923,  2224,  1681,  1682,  1683,  1684,  1758,   923,
     819,  2027,   821,   822,   823,   174,  1688,  1689,   827,  2048,
     589,   830,   831,   593,   833,   834,   835,  1690,  2378,   838,
     841,  2049,   843,   844,   845,   846,   847,   848,   594,   854,
    1365,  1366,   859,   860,   861,   862,   863,   864,   865,   866,
     867,   868,   939,   870,   871,   872,   873,   597,  2056,   876,
     878,   879,   880,   881,   882,   883,   884,   885,  1706,   888,
    2057,   896,   897,   598,  2411,   992,  1695,  1696,  1697,  1698,
     901,  1774,  1611,  2427,  2416,   101,   102,   599,  1121,  2428,
    1121,   145,  1621,  1622,  1623,  1922,  1707,  2126,   783,  1867,
     923,   101,   102,   940,   784,  1917,  1918,   785,  2429,   786,
     956,   101,   102,   101,   102,  1127,  1775,   789,   101,   102,
     101,   102,   101,   102,   791,  1875,   923,  1190,   980,  2141,
     802,   982,  2142,   923,  2250,   923,  1257,   101,   102,   101,
     102,   101,   102,   941,   942,   101,   102,   803,   101,   102,
    1408,  1896,  2430,  1925,  1381,  1937,   101,   102,  1382,  2099,
     804,  1008,  2101,   805,   906,   101,   102,   101,   102,   806,
    2172,  1708,  1019,  1019,  1110,  1899,  1644,   101,   102,  2347,
     807,  2350,  1379,  1380,  1075,   943,  1802,  1077,  1078,  1776,
    1777,  2358,   808,   907,  1403,  1404,  1057,  1058,  1059,   908,
    1060,  1061,  1062,  1063,  1064,  1806,   909,   910,  1065,   101,
     102,   809,   101,   102,  1073,  1073,  2196,   810,  1073,  1073,
     811,   101,   102,  2376,  1741,    38,  2377,  1092,  1745,  1956,
    1095,   813,  1099,   812,  1101,  2394,   814,  1105,  1770,  1974,
    1377,  1378,  1379,  1380,   944,  1381,   826,  1109,  1778,  1382,
    2072,  2073,  2073,  1779,  2106,  2107,   820,  1308,  1812,  1813,
     837,   906,  1126,  1128,  1129,   842,   886,  1135,  1709,  1710,
    1196,   903,  1870,   906,    43,   911,   912,   913,   914,   915,
     929,  1121,  1791,   946,  1809,   947,   948,   949,  1814,  1815,
     907,  1821,  1822,  1823,  1824,  1239,   908,   970,  1409,  1711,
     975,   976,   907,   909,   910,  1780,  1920,   990,   908,  1013,
    1041,  1009,  1921,  1162,  1163,   909,   910,  1377,  1378,  1379,
    1380,  1175,  1049,   906,  1680,  1040,  1069,  1070,  1182,  1183,
    1184,  1071,  1186,  1187,  1188,  1081,  1082,  1083,  1084,  1076,
    1191,  1079,  1080,   916,   917,   918,   919,   101,   102,  1173,
    1174,  1086,   907,  1088,  1089,  1923,  1093,  1096,   908,  1094,
    1781,   906,  1100,  1102,  1103,   909,   910,  1106,  1107,  1924,
    1108,   920,   911,   912,   913,   914,   915,   101,   102,  1479,
    1480,  1111,  1258,  1868,   911,   912,   913,   914,   915,  1112,
     907,  1322,  1113,  1114,  1115,  1878,   908,   101,   102,  1498,
    1499,  1116,  1117,   909,   910,  1119,  1019,  1019,  1019,  1120,
    1130,  1133,  1019,  1153,  1019,  1124,  1019,  1136,  1137,  1138,
    1019,   101,   102,  1536,  1537,  1415,  1139,   704,  1816,  1817,
    1818,  1819,  1144,  1423,   911,   912,   913,   914,   915,  1140,
     916,   917,   918,   919,  1377,  1378,  1379,  1380,  1141,   921,
    1142,  1149,   916,   917,   918,   919,  1731,  1732,  1733,  1734,
    1735,  1736,  1737,  1150,   101,   102,  1151,  1152,   920,  1154,
    1155,  1158,   911,   912,   913,   914,   915,  1159,   922,   923,
     920,  1176,  2158,  1164,  1171,  1464,  1177,  1465,  1466,   101,
     102,  1932,  1933,   101,   102,  2024,  2025,  1189,  1193,  1470,
    1198,  1194,   916,   917,   918,   919,  1195,  1213,  1228,  1473,
    1229,  1232,  1233,  1234,  1235,  1953,  1236,  1237,  1238,  1478,
    1481,  1240,  1241,  1484,  1244,  1245,  1323,  1246,  1247,  2031,
     920,  1253,  1420,  1259,  1492,  1421,  1494,  1428,  1429,  1430,
     916,   917,   918,   919,  1431,  1500,   921,  1432,  1433,  1434,
    1435,  1445,  1505,  1436,  1437,  1438,  1447,  1439,   921,  1450,
    1504,   101,   102,  1440,  1441,  1442,  1443,  1444,   920,  1451,
    1452,  1453,  1516,  1454,  1457,   987,   923,  1455,  1460,  1523,
    1459,  1524,  1461,  1467,  1468,  1469,  1471,  1165,  1528,  1485,
    1476,  1477,  1486,  1493,  1495,  1496,  1497,  1503,  1506,  1512,
    1514,  1521,  1522,  1515,  1517,  1518,  1520,  1529,   921,  1533,
    1530,  1534,  1543,  1571,  1598,  1535,  1542,  1615,  1625,  1627,
    1638,  1639,  1640,  1628,  1641,  1652,  1703,  1724,  1644,  1673,
    1661,  1538,  1676,  1744,  1541,  1745,  1662,  1254,  1668,  1671,
    1672,  1769,  1700,  1705,  1784,  1721,   921,  1739,  1544,  1742,
    1545,  1546,  1547,  1548,  1549,  1550,  1551,  1751,  1553,  1554,
    1555,  1556,  1557,  1756,  1788,  1795,  1565,  1566,  1567,  1568,
    1770,  1796,  1789,  1572,  1797,  1574,  1575,  1807,  1803,  1804,
    1805,  1810,  1580,  1825,  1582,  1583,  1584,  1585,  1586,  1587,
    1826,  1827,  1828,  1829,  1599,  1830,  1607,  1608,  2157,  -423,
       1,  1834,  -423,  -423,  -423,  -423,  -423,  1835,  1836,  1837,
    1838,  1840,  1841,  1843,  1842,  1844,  1845,  1849,  1629,  1630,
    1631,  1632,  1846,  1850,  1853,  1910,  2093,  1854,  1911,  1634,
    1635,  1636,  1637,  1861,  1855,  -423,  -423,  1856,  1857,  1642,
    1858,  1859,  1860,  2159,  -423,  1913,  1862,  1863,  1258,  1864,
    1865,  1866,  1666,  1667,  -423,  1669,  1670,  1879,  1881,  -423,
    1886,  1892,  1893,  1894,  1901,  -423,  1903,  1019,  1019,  1019,
    1019,  1914,  1915,  1916,  -423,  1919,  1408,  1928,  1930,  1019,
    1019,  1931,  1940,  -423,  1935,  1942,  1943,  1944,  1945,  1966,
    1019,  1019,  1019,  1019,  1019,   906,  1970,  1946,  1947,  1948,
    1949,  1950,  1959,  1960,  1976,  1986,  1701,  1988,  -423,  1990,
     906,  -423,  1994,  1995,  1722,  2000,  -423,  -423,  -423,  1996,
    1997,  2002,  -423,  2034,   907,  2003,  2005,  2006,  2035,  2007,
     908,  2037,  2008,  2042,  2009,  1738,  2010,   909,   910,   907,
    2011,  2012,  2013,  2014,  2054,   908,  2015,  1743,  2016,  2018,
    2019,  1748,   909,   910,  2020,  2021,  2029,  2032,  2033,  2043,
    2059,  1755,  2044,  2045,  2046,  2047,  2050,  2051,  2052,  2053,
    2062,  2067,  -423,  2068,  -423,  2074,  2077,  2081,  2082,  2083,
    1833,  -423,  2084,  1785,  -423,  1787,  2087,  2095,  2098,  2097,
    2104,  2116,  2110,  -423,  -423,  2114,  2111,  2113,  2118,  2115,
    2117,  2197,  2127,  2130,  -423,  2132,   911,   912,   913,   914,
     915,  1798,  1799,  2133,  -423,  2131,  2134,  2150,  2135,  2140,
    2152,   911,   912,   913,   914,   915,  2143,  2153,  2154,  2155,
    -423,  2156,  2176,  2163,  2173,  2180,  -423,  2182,  2189,  2192,
    -423,  2193,  2195,  -423,  2199,  2202,  2203,  2205,  2093,  2208,
    2207,  2213,  2218,  2219,  2220,  -423,  2225,  2226,  2229,  2232,
    1839,  2235,  2236,  2244,  2249,  2255,  2256, -1140,  1847,  1848,
    2259,  2251, -1141,  2272,   916,   917,   918,   919,  2273,  2275,
    2280,  2286,  2287,  2288,  2292,  2293,  2294,  2295,  -423,   916,
     917,   918,   919,  2296,  2299,  2304,  2309,  2313,  2312,  2329,
    2333,  2340,   920,  2318,  2335,  2344,  2356,  2364,  2366,  2379,
    2381,  2382,  2385,  2387,  2391,  2396,  2401,   920,  2403,  1882,
    1883,  1884,  1885,  2405,  1888,  1889,  1890,  1891,  2409,  2419,
    -423,  1895,  2404,  1897,  1898,  2438,  2357,  2407,  2245,  2408,
    2413,  -423,  1904,  2414,  2174,  2418,  2424,   874,  1909,  2441,
    2442,  2443,  2445,  2228,  1912,  1456,   727,  1483,  1750,   728,
    1957,  1975,  2339,  2363,  1502,  1712,  1509,  1792,  1170,  1019,
    1713,  1161,   932,  1540,  1871,  1019,  1626,  2028,  -423,  -423,
     921,   729,  -423,  1872,  2346,  1532,  2252,  -423,  1900,  2206,
    -423,  -423,  1926,  2061,   780,   921,  2089,  2171,  1783,  2167,
    2161,  1934,   730,  1427,  1419,     0,     0,     0,     0,  1869,
       0,     0,     0,  1938,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2030,     0,     0,     0,     0,     0,
       0,     0,   687,     0,  1954,  1955,     0,   731,  1958,     0,
       0,   732,     0,   733,     0,     0,     0,     0,     0,     0,
       0,  1968,  1969,     0,     0,     0,     0,     0,     0,  1972,
    1973,     0,     0,     0,  1978,  1979,  1980,  1981,     0,  1983,
    1984,  1985,     0,   734,  1987,     0,     0,     0,  1992,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   735,     0,  2017,
       0,     0,     0,     0,     0,   736,     0,  2022,  2023,     0,
    2026,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   737,
       0,     0,     0,     0,     0,     0,   738,     0,     0,     0,
     739,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2055,     0,
       0,     0,     0,     0,  2058,     0,     0,     0,     0,  2060,
       0,     0,     0,  2063,  2064,  2065,  2066,     0,     0,     0,
       0,     0,     0,   740,   741,     0,     0,     0,  2078,     0,
       0,     0,   742,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   743,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2096,     0,   744,     0,     0,
       0,  2100,  2100,  2102,  2103,     0,  2105,     0,     0,     0,
    2109,     0,     0,  2112,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2119,     0,     0,  2121,  2100,
    2123,  2124,  2125,     0,     0,     0,     0,     0,     0,     0,
    2129,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2139,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2145,  2146,     0,     0,  2147,  2148,
       0,     0,  2151,     0,     0,     0,     0,     0,     0,     0,
       0,  1019,  1019,     0,  2160,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2166,     0,     0,  2168,  2169,  2169,
    2100,     0,     0,     0,     0,  2175,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2183,  2184,  2185,     0,  2187,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2200,  2201,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2209,  2210,
    2212,     0,  2214,  2215,  2216,  2217,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2221,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2231,     0,     0,     0,     0,     0,     0,  2238,     0,  2243,
       0,     0,  2246,  2247,     0,     0,     0,     0,     0,  2253,
       0,     0,     0,  2254,     0,     0,     0,     0,     0,  2257,
    2258,     0,  2260,     0,     0,     0,     0,     0,     0,  2264,
       0,     0,     0,     0,     0,  2269,  2270,     0,     0,  2271,
       0,     0,     0,     0,     0,  2277,     0,     0,     0,  2278,
    2279,     0,     0,     0,     0,     0,     0,     0,  2100,     0,
       0,  2285,     0,     0,     0,     0,     0,     0,  2289,  2290,
    2291,     0,     0,     0,     0,     0,  2297,  2298,     0,     0,
       0,     0,  2302,     0,     0,  2303,     0,     0,     0,  2305,
    2306,     0,     0,     0,     0,     0,  2314,  2315,     0,  2316,
    2317,     0,     0,     0,  2321,     0,  2324,     0,     0,  2328,
       0,     0,     0,     0,  2334,     0,     0,     0,  2100,  2338,
       0,     0,     0,     0,  2341,  2342,  2343,     0,     0,     0,
       0,     0,  2348,     0,     0,  2351,     0,     0,     0,  2355,
       0,     0,  2100,  2359,  2360,  2361,     0,     0,     0,     0,
       0,     0,     0,  2365,     0,     0,     0,     0,     0,     0,
       0,     0,  2100,  2100,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2384,     0,  2386,     0,     0,     0,  2388,
    2389,     0,     0,  2392,  2393,     0,     0,     0,  2100,  2395,
       0,  2398,  2399,     0,     0,     0,     0,  2402,     0,     0,
       0,  2406,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2412,     0,     0,     0,     0,     0,     0,     0,  2415,
       0,  2100,  2417,     0,     0,     0,  2100,     0,     0,  2420,
    2421,  2422,  2423,     0,     0,     0,     0,  2432,  2433,     0,
    2435,     0,  2436,  2437,     0,  2439,  2440,     0,     0,     0,
    2444,   175,   176,   177,   178,   179,   180,   181,   182,   183,
     184,   185,   186,   187,   188,   189,   190,   191,   192,   193,
     194,   195,   196,   197,   198,   199,   200,   201,   202,   203,
     204,   205,   206,   207,   208,   209,   210,   211,   212,   213,
     214,   215,   216,   217,   218,   219,   220,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,   239,   240,   241,   242,   243,
     244,   245,   246,   247,   248,   249,   250,   251,   252,   253,
     254,   255,   256,   257,   258,   259,   260,   261,   262,   263,
     264,   265,   266,   267,   268,   269,   270,   271,   272,   273,
     274,   275,   276,   277,   278,   279,   280,   281,   282,   283,
     284,   285,   286,   287,   288,   289,   290,   291,   292,   293,
     294,   295,   296,   297,   298,   299,   300,   301,   302,   303,
     304,   305,   306,   307,   308,   309,   310,   311,   312,   313,
     314,   315,   316,   317,   318,   319,   320,   321,   322,   323,
     324,   325,   326,   327,   328,   329,   330,   331,   332,   333,
     334,   335,   336,   337,   338,   339,   340,   341,   342,   343,
     344,   345,   346,   347,   348,   349,   350,   351,   352,   353,
     354,   355,   356,   357,   358,   359,   360,   361,   362,   363,
     364,   365,   366,   367,   368,   369,   370,   371,   372,   373,
     374,   375,   376,   377,   378,   379,   380,   381,   382,   383,
     384,   385,   386,   387,   388,   389,   390,   391,   392,   393,
     394,   395,   396,   397,   398,   399,   400,   401,   402,   403,
     404,   405,   406,   407,   408,   409,   410,   411,   412,   413,
     414,   415,   416,   417,   418,   419,   420,   421,   422,   423,
     424,   425,   426,   427,   428,   429,   430,   431,   432,   433,
     434,   435,   436,   437,   438,   439,   440,   441,   442,   443,
     444,   445,   446,   447,   448,   449,   450,   451,   452,   453,
     454,   455,   456,   457,   458,   459,   460,   461,   462,   463,
     464,   465,   466,   467,   468,   469,   470,   471,   472,   473,
     474,   475,   476,   477,   478,   479,   480,   481,   482,   483,
     484,   485,   486,   487,   488,   489,   490,   491,   492,   493,
     494,   495,   496,   497,   498,   499,   500,   501,   502,   503,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,   517,   518,   519,   520,   521,   522,   523,
     524,   525,   526,   527,   528,   529,   530,   531,   532,   533,
     534,   535,   536,   537,   538,   539,   540,   541,   542,   543,
     544,   545,   546,   547,   548,   549,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   575,   576,   577,   578,   579,   580,   581,   582,   583,
     584,   585,     0,     0,   586,     5,     6,     7,     8,     9,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    10,    11,
       0,     0,     0,     0,     0,     0,     0,    12,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    13,     0,     0,
       0,     0,    14,     0,     0,     0,     0,     0,    15,     0,
       0,     0,     0,     0,     0,     0,     0,    16,     0,     0,
       0,     0,     0,     0,     0,     0,    17,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    18,     0,     0,    19,     0,     0,     0,     0,    20,
      21,    22,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    23,     0,    24,     0,     0,
       0,     0,     0,     0,    25,     0,     0,    26,     0,     0,
       0,     0,     0,     0,     0,     0,    27,    28,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    29,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    30,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    31,     0,     0,     0,     0,     0,    32,
       0,     0,     0,    33,     0,     0,    34,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   619,    35,     0,
       0,     0,   620,   621,   622,   623,     0,   624,   625,   626,
     627,   628,     0,   629,   630,     0,   631,     0,     0,   632,
       0,     0,     0,     0,     0,   633,     0,     0,     0,     0,
     634,    36,   635,     0,   636,     0,     0,     0,   637,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   638,   639,     0,     0,     0,     0,     0,     0,
     640,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    37,     0,     0,   641,     0,     0,     0,
       0,   642,   643,     0,    38,     0,   644,   645,   646,   647,
       0,     0,   648,     0,     0,     0,   649,     0,     0,     0,
       0,     0,     0,   650,   651,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    39,    40,     0,     0,    41,   652,     0,     0,     0,
      42,     0,     0,    43,    44,     0,     0,     0,     0,     0,
       0,     0,  1260,     0,     0,  1261,   653,     0,     0,     0,
    1262,  1263,  1264,  1265,  1266,  1267,     0,  1268,  1269,  1270,
       0,  1271,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   654,
     655,     0,  1272,   656,     0,   657,     0,   658,     0,     0,
       0,     0,     0,     0,     0,   659,     0,   660,   661,  1273,
     662,   663,   664,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1274,     0,     0,     0,     0,     0,   665,     0,
       0,     0,     0,     0,  1275,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1276,     0,  1277,
       0,  1278,  1279,     0,     0,     0,     0,     0,     0,   666,
       0,   687,     0,     0,     0,     0,     0,   667,     0,   668,
       0,   669,     0,     0,  1280,     0,     0,     0,     0,     0,
       0,     0,   670,     0,     0,     0,     0,  1281,     0,     0,
       0,     0,     0,     0,     0,  1282,     0,     0,  1283,     0,
       0,     0,  1284,     0,     0,     0,     0,     0,     0,     0,
     671,  1285,     0,     0,  1286,   672,   673,   674,   675,   676,
       0,     0,     0,   677,   678,     0,     0,     0,  1325,     0,
       0,     0,     0,     0,     0,  1287,     0,  1326,  1327,  1328,
    1329,   679,  1330,     0,  1331,  1288,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   680,  1332,
       0,     0,     0,     0,     0,     0,     0,     0,  1289,  1290,
       0,     0,  1333,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   681,     0,     0,   682,     0,     0,     0,     0,
       0,     0,  1291,     0,     0,  1292,     0,     0,     0,  1293,
       0,     0,  1294,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1295,     0,     0,     0,     0,     0,  1296,     0,
    1297,     0,  1298,  1299,     0,     0,     0,  1334,  1335,     0,
       0,     0,     0,     0,  1300,     0,     0,     0,     0,  1301,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1302,     0,     0,     0,     0,     0,  1303,
       0,     0,     0,     0,  1304,     0,     0,     0,     0,  1336,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1305,     0,     0,     0,     0,     0,     0,     0,
       0,  1306,     0,     0,  1337,     0,     0,  1307,     0,  1338,
    1339,  1340,     0,  1341,     0,  1342,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1343,
    1344,  1345,  1346,     0,     0,     0,     0,     0,  1347,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1348,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1349,     0,  1350,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1351,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1352
  };

  const short
  Parser::yycheck_[] =
  {
      10,    11,    12,   146,  1029,  1315,   745,    17,    18,   782,
      20,    21,    22,    23,    24,   152,   934,    20,    28,    29,
      30,  1317,    71,    33,    90,    35,  1016,  1753,    20,    62,
      40,    20,  1022,    81,  1024,   119,  1026,  1006,  1007,    20,
    1030,   997,   998,  1012,    62,    20,    19,    65,   132,    25,
     121,    90,    71,    15,   784,   785,   786,    20,    17,   108,
    1786,    63,    64,    81,    82,    20,   115,    85,    87,  1878,
      98,    80,   698,    20,   123,    94,  1028,   169,   123,   113,
    1032,  1033,  1034,   103,   104,   121,   102,   168,    54,   185,
      73,    89,    61,   123,   188,   133,    67,    58,   123,   250,
     150,    20,   294,    76,   123,   174,   294,   125,   150,  1929,
     736,   288,    20,   123,   123,   346,   171,   185,   121,   745,
      93,   128,   178,   133,    52,   116,   159,   320,   116,     3,
       4,     5,   151,   297,   288,   185,   123,   116,   373,   280,
     185,   245,   177,   185,   123,    96,   116,     0,   168,    20,
      74,   179,  1384,   123,   189,    13,   417,   416,   121,   300,
     144,   417,   288,   414,   415,   171,   201,   123,   144,   428,
     288,   183,   170,   277,   230,   415,   185,   428,  1410,   159,
     428,   288,   144,   192,   432,   194,   255,   926,   172,   162,
     174,   416,   288,   428,   820,    69,    70,   428,   196,   123,
     207,    59,   121,   271,   416,   185,   172,   123,   174,   160,
     181,    52,   280,   405,   406,   270,   428,   405,   406,   267,
     238,   414,   415,    97,   259,  1951,   271,    74,   237,   185,
     416,   218,   300,  2042,   416,   280,   428,   391,   288,   416,
     428,   416,   298,   869,  1017,  1018,   288,   250,   416,   267,
     989,   301,   301,  1026,  1172,   300,   407,   408,   250,   301,
    2080,   250,   416,   261,   428,   284,   304,   258,   277,   250,
     258,   255,   416,   169,   148,   250,   123,   142,   347,   258,
     349,   283,   301,   285,   158,   303,   269,   250,   258,   255,
     416,   297,   322,   342,   350,   250,  1026,   183,   416,   183,
     428,  1031,   428,   250,   169,   233,  1616,   401,   389,   416,
     428,   284,   273,   348,   246,   941,   314,   249,   284,   369,
     416,   376,  1618,   391,  1620,   417,   342,   369,   428,   354,
     428,   250,  2141,   401,   416,   401,   416,   429,   183,   322,
     180,   303,   250,   302,   372,   364,   391,   395,   414,   415,
     392,   384,   313,   424,   282,   404,   401,   330,  1314,   292,
    1316,   432,   401,  1388,   428,   341,   368,   359,   370,   389,
     339,   997,   998,   999,  1000,   341,   382,   395,   397,   250,
     414,   415,   391,   392,   428,   361,   362,   423,   424,   425,
     426,   246,   233,   428,   249,   429,   432,   395,   292,   239,
     240,   428,   414,   415,   407,   408,   409,   410,   411,   412,
     415,  1401,  1402,  1403,  1404,   407,   408,   429,   407,   408,
     423,   424,   425,   426,   400,    96,   407,   408,   431,   432,
     433,   434,   407,   408,   409,   410,   411,   412,   428,   428,
     432,   282,  1068,   432,   407,   408,   414,   415,   423,   424,
     425,   426,   407,   408,   292,   430,   431,   331,   433,   434,
     407,   408,   409,   410,   411,   412,  1418,   292,  2194,   432,
      69,    70,   319,   416,  1426,   430,   423,   424,   425,   426,
     246,   414,   415,   249,   431,   332,   433,   434,   407,   408,
     409,   410,   411,   412,  1790,   416,   429,   371,   416,   407,
     408,   409,   410,   411,   412,   424,   353,   428,   351,   374,
     428,   424,   431,   432,   433,   434,   424,   430,   414,   415,
     414,   415,   430,   431,   123,   433,   434,   423,   414,   415,
     414,   415,   192,   429,   429,   430,   407,   408,   409,   410,
     411,   412,  1168,   429,   416,   429,   414,   415,   413,   414,
     415,  1507,   417,   424,   111,   428,  2282,   416,   423,  1185,
     431,   429,   433,   434,   429,   122,   428,  1877,   428,   414,
     415,   414,   415,   405,   406,  2301,   414,   415,   414,   415,
    1876,   414,   415,   431,   429,   428,   429,   414,   415,   414,
     415,   428,   428,  1511,   414,   415,   429,   430,   428,   416,
    2326,  2327,   429,   430,  1377,  1378,  1379,  1380,   428,   429,
     620,   428,   622,   623,   624,   428,  1389,  1390,   628,   416,
     428,   631,   632,   428,   634,   635,   636,  1400,  2354,   639,
     640,   428,   642,   643,   644,   645,   646,   647,   428,   649,
     405,   406,   652,   653,   654,   655,   656,   657,   658,   659,
     660,   661,   209,   663,   664,   665,   666,   428,   416,   669,
     670,   671,   672,   673,   674,   675,   676,   677,   185,   679,
     428,   681,   682,   428,  2400,   759,  1406,  1407,  1408,  1409,
     690,    18,  1308,   247,  2410,   414,   415,   428,  1314,   253,
    1316,   185,  1318,  1319,  1320,  1685,   213,  1993,   431,   428,
     429,   414,   415,   260,   431,  1674,  1675,   431,   272,   431,
     720,   414,   415,   414,   415,   428,    53,   428,   414,   415,
     414,   415,   414,   415,   428,   428,   429,   428,   738,  2039,
     428,   741,   428,   429,   428,   429,   428,   414,   415,   414,
     415,   414,   415,   300,   301,   414,   415,   428,   414,   415,
     424,   428,   316,   428,   428,   428,   414,   415,   432,   428,
     428,   771,   428,   428,   114,   414,   415,   414,   415,   276,
     428,   288,   782,   783,   858,   428,   429,   414,   415,   428,
     211,   428,   425,   426,   825,   342,  1525,   828,   829,   126,
     127,   428,   210,   143,   407,   408,   806,   807,   808,   149,
     810,   811,   812,   813,   814,  1544,   156,   157,   818,   414,
     415,   379,   414,   415,   824,   825,  2126,   208,   828,   829,
     212,   414,   415,   428,  1450,   342,   428,   837,   429,   430,
     840,   232,   842,   248,   844,   428,   398,   847,   429,   430,
     423,   424,   425,   426,   401,   428,   428,   857,   185,   432,
    1922,  1923,  1924,   190,   152,   153,   429,   994,  1548,  1549,
      86,   114,   872,   873,   874,    58,   254,   877,   385,   386,
     954,   123,  1611,   114,   391,   225,   226,   227,   228,   229,
     416,  1507,  1508,   342,  1546,   416,   182,   416,  1550,  1551,
     143,  1553,  1554,  1555,  1556,   979,   149,   185,  1923,   416,
     416,   416,   143,   156,   157,   242,  1679,   123,   149,   123,
     164,   416,  1685,   923,   924,   156,   157,   423,   424,   425,
     426,   931,   289,   114,   430,   416,   428,   428,   938,   939,
     940,   428,   942,   943,   944,   263,   264,   265,   266,   428,
     950,   428,   428,   293,   294,   295,   296,   414,   415,   416,
     417,   428,   143,   428,   193,  1685,   428,   428,   149,   429,
     297,   114,   428,   428,   428,   156,   157,   428,   428,  1699,
     428,   321,   225,   226,   227,   228,   229,   414,   415,   416,
     417,   428,   992,  1609,   225,   226,   227,   228,   229,   428,
     143,  1001,   428,   428,   428,  1621,   149,   414,   415,   416,
     417,   428,   428,   156,   157,   401,  1016,  1017,  1018,   428,
     416,   188,  1022,   416,  1024,   428,  1026,   428,   428,   428,
    1030,   414,   415,   416,   417,  1035,   428,  1170,   263,   264,
     265,   266,   402,  1043,   225,   226,   227,   228,   229,   428,
     293,   294,   295,   296,   423,   424,   425,   426,   428,   399,
     428,   428,   293,   294,   295,   296,  1429,  1430,  1431,  1432,
    1433,  1434,  1435,   428,   414,   415,   428,   428,   321,   416,
     416,   428,   225,   226,   227,   228,   229,   416,   428,   429,
     321,   416,  2072,   428,   428,  1095,   123,  1097,  1098,   414,
     415,   416,   417,   414,   415,   416,   417,   416,   428,  1109,
      68,   428,   293,   294,   295,   296,   428,    65,   428,  1119,
     416,   416,   428,   428,   416,  1741,   428,   416,   428,  1129,
    1130,   428,   416,  1133,    68,   428,   123,   428,   428,  1868,
     321,   428,   416,   428,  1144,   428,  1146,   416,   416,   416,
     293,   294,   295,   296,   416,  1155,   399,   416,   416,   416,
     416,   151,  1162,   428,   428,   428,   110,   428,   399,   429,
     416,   414,   415,   428,   428,   428,   428,   428,   321,   428,
     428,   428,  1182,   428,   401,   428,   429,   428,   428,  1189,
     429,  1191,   429,   428,   428,   428,   401,   428,  1198,   429,
     428,   428,   428,   401,   401,   428,   428,   428,   428,   428,
     428,    68,   366,   428,   428,   428,   428,   428,   399,   428,
     416,   428,   416,   416,   416,   428,   428,   416,   428,   365,
     417,   417,   417,   428,   327,   356,    84,   241,   429,   416,
     428,  1241,   416,    92,  1244,   429,   428,   428,   428,   428,
     428,   401,   428,   428,   342,   428,   399,   428,  1258,   428,
    1260,  1261,  1262,  1263,  1264,  1265,  1266,   430,  1268,  1269,
    1270,  1271,  1272,   428,    91,   416,  1276,  1277,  1278,  1279,
     429,   428,   430,  1283,   428,  1285,  1286,   185,   428,   428,
     428,   185,  1292,   428,  1294,  1295,  1296,  1297,  1298,  1299,
     428,   428,   428,   428,  1304,   428,  1306,  1307,  2071,     0,
       1,   428,     3,     4,     5,     6,     7,   428,   428,   428,
     416,   416,   428,   428,   416,   428,   416,   428,  1328,  1329,
    1330,  1331,   417,   416,   428,   378,  1952,   428,    97,  1339,
    1340,  1341,  1342,   416,   428,    36,    37,   428,   428,  1349,
     428,   428,   428,  2073,    45,   113,   428,   428,  1358,   428,
     428,   428,  1362,  1363,    55,  1365,  1366,   428,   416,    60,
     416,   428,   428,   428,   416,    66,   428,  1377,  1378,  1379,
    1380,   113,   113,   113,    75,   359,   424,   416,   416,  1389,
    1390,   416,   416,    84,   123,   417,   428,   428,   428,   402,
    1400,  1401,  1402,  1403,  1404,   114,   401,   428,   428,   428,
     428,   428,   428,   428,   275,   428,  1416,   402,   109,   287,
     114,   112,   428,   428,  1424,   113,   117,   118,   119,   428,
     428,   428,   123,   123,   143,   428,   428,   428,   131,   428,
     149,   113,   428,   416,   428,  1445,   428,   156,   157,   143,
     428,   428,   428,   428,   167,   149,   428,  1457,   428,   428,
     428,  1461,   156,   157,   428,   428,   428,   428,   428,   428,
     356,  1471,   428,   428,   428,   428,   428,   428,   428,   428,
     403,   428,   173,   428,   175,   393,   416,   428,   428,   416,
    1564,   182,   393,  1493,   185,  1495,   297,    56,   430,   428,
     428,   243,   428,   194,   195,   402,   428,   428,   297,   402,
     416,  2127,   300,   416,   205,   428,   225,   226,   227,   228,
     229,  1521,  1522,   428,   215,   416,   428,   356,   428,   428,
     348,   225,   226,   227,   228,   229,   428,   348,   348,   348,
     231,   121,   430,   428,   428,   188,   237,   430,   402,   342,
     241,   401,   428,   244,    68,   428,    68,   428,  2174,   428,
     430,   428,   428,   401,   428,   256,   342,   357,   129,   428,
    1570,   242,   346,   428,   186,   428,   428,   428,  1578,  1579,
     387,  2197,   428,   402,   293,   294,   295,   296,   130,   268,
     342,   348,   348,   428,   185,   185,    83,   428,   289,   293,
     294,   295,   296,   394,   357,   297,   401,   428,   416,    56,
     402,   416,   321,   428,   401,   428,   428,   416,   416,   342,
     402,   401,   428,    86,   320,   375,    92,   321,   428,  1629,
    1630,  1631,  1632,   416,  1634,  1635,  1636,  1637,   342,   122,
     331,  1641,   428,  1643,  1644,   416,  2331,   428,  2191,   428,
     428,   342,  1652,   428,  2094,   428,   428,   667,  1658,   428,
     428,   428,   428,  2178,  1664,  1087,    74,  1131,  1462,    77,
    1746,  1771,  2311,  2338,  1156,  1422,  1170,  1509,   928,  1679,
    1422,   905,   705,  1242,  1612,  1685,  1324,  1851,   379,   380,
     399,    99,   383,  1614,  2320,  1231,  2198,   388,  1645,  2145,
     391,   392,  1702,  1910,   153,   399,  1941,  2089,  1491,  2085,
    2075,  1711,   120,  1045,  1037,    -1,    -1,    -1,    -1,   428,
      -1,    -1,    -1,  1723,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   428,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   150,    -1,  1744,  1745,    -1,   155,  1748,    -1,
      -1,   159,    -1,   161,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1761,  1762,    -1,    -1,    -1,    -1,    -1,    -1,  1769,
    1770,    -1,    -1,    -1,  1774,  1775,  1776,  1777,    -1,  1779,
    1780,  1781,    -1,   191,  1784,    -1,    -1,    -1,  1788,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   245,    -1,  1839,
      -1,    -1,    -1,    -1,    -1,   253,    -1,  1847,  1848,    -1,
    1850,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   277,
      -1,    -1,    -1,    -1,    -1,    -1,   284,    -1,    -1,    -1,
     288,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1898,    -1,
      -1,    -1,    -1,    -1,  1904,    -1,    -1,    -1,    -1,  1909,
      -1,    -1,    -1,  1913,  1914,  1915,  1916,    -1,    -1,    -1,
      -1,    -1,    -1,   331,   332,    -1,    -1,    -1,  1928,    -1,
      -1,    -1,   340,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   353,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1955,    -1,   365,    -1,    -1,
      -1,  1961,  1962,  1963,  1964,    -1,  1966,    -1,    -1,    -1,
    1970,    -1,    -1,  1973,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1985,    -1,    -1,  1988,  1989,
    1990,  1991,  1992,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2000,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2037,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2054,  2055,    -1,    -1,  2058,  2059,
      -1,    -1,  2062,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2071,  2072,    -1,  2074,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2084,    -1,    -1,  2087,  2088,  2089,
    2090,    -1,    -1,    -1,    -1,  2095,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2114,  2115,  2116,    -1,  2118,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2136,  2137,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2148,  2149,
    2150,    -1,  2152,  2153,  2154,  2155,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2168,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2180,    -1,    -1,    -1,    -1,    -1,    -1,  2187,    -1,  2189,
      -1,    -1,  2192,  2193,    -1,    -1,    -1,    -1,    -1,  2199,
      -1,    -1,    -1,  2203,    -1,    -1,    -1,    -1,    -1,  2209,
    2210,    -1,  2212,    -1,    -1,    -1,    -1,    -1,    -1,  2219,
      -1,    -1,    -1,    -1,    -1,  2225,  2226,    -1,    -1,  2229,
      -1,    -1,    -1,    -1,    -1,  2235,    -1,    -1,    -1,  2239,
    2240,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2248,    -1,
      -1,  2251,    -1,    -1,    -1,    -1,    -1,    -1,  2258,  2259,
    2260,    -1,    -1,    -1,    -1,    -1,  2266,  2267,    -1,    -1,
      -1,    -1,  2272,    -1,    -1,  2275,    -1,    -1,    -1,  2279,
    2280,    -1,    -1,    -1,    -1,    -1,  2286,  2287,    -1,  2289,
    2290,    -1,    -1,    -1,  2294,    -1,  2296,    -1,    -1,  2299,
      -1,    -1,    -1,    -1,  2304,    -1,    -1,    -1,  2308,  2309,
      -1,    -1,    -1,    -1,  2314,  2315,  2316,    -1,    -1,    -1,
      -1,    -1,  2322,    -1,    -1,  2325,    -1,    -1,    -1,  2329,
      -1,    -1,  2332,  2333,  2334,  2335,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2343,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  2352,  2353,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2363,    -1,  2365,    -1,    -1,    -1,  2369,
    2370,    -1,    -1,  2373,  2374,    -1,    -1,    -1,  2378,  2379,
      -1,  2381,  2382,    -1,    -1,    -1,    -1,  2387,    -1,    -1,
      -1,  2391,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2401,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2409,
      -1,  2411,  2412,    -1,    -1,    -1,  2416,    -1,    -1,  2419,
    2420,  2421,  2422,    -1,    -1,    -1,    -1,  2427,  2428,    -1,
    2430,    -1,  2432,  2433,    -1,  2435,  2436,    -1,    -1,    -1,
    2440,     3,     4,     5,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    51,
      52,    53,    54,    55,    56,    57,    58,    59,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    69,    70,    71,
      72,    73,    74,    75,    76,    77,    78,    79,    80,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    91,
      92,    93,    94,    95,    96,    97,    98,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   109,   110,   111,
     112,   113,   114,   115,   116,   117,   118,   119,   120,   121,
     122,   123,   124,   125,   126,   127,   128,   129,   130,   131,
     132,   133,   134,   135,   136,   137,   138,   139,   140,   141,
     142,   143,   144,   145,   146,   147,   148,   149,   150,   151,
     152,   153,   154,   155,   156,   157,   158,   159,   160,   161,
     162,   163,   164,   165,   166,   167,   168,   169,   170,   171,
     172,   173,   174,   175,   176,   177,   178,   179,   180,   181,
     182,   183,   184,   185,   186,   187,   188,   189,   190,   191,
     192,   193,   194,   195,   196,   197,   198,   199,   200,   201,
     202,   203,   204,   205,   206,   207,   208,   209,   210,   211,
     212,   213,   214,   215,   216,   217,   218,   219,   220,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,   239,   240,   241,
     242,   243,   244,   245,   246,   247,   248,   249,   250,   251,
     252,   253,   254,   255,   256,   257,   258,   259,   260,   261,
     262,   263,   264,   265,   266,   267,   268,   269,   270,   271,
     272,   273,   274,   275,   276,   277,   278,   279,   280,   281,
     282,   283,   284,   285,   286,   287,   288,   289,   290,   291,
     292,   293,   294,   295,   296,   297,   298,   299,   300,   301,
     302,   303,   304,   305,   306,   307,   308,   309,   310,   311,
     312,   313,   314,   315,   316,   317,   318,   319,   320,   321,
     322,   323,   324,   325,   326,   327,   328,   329,   330,   331,
     332,   333,   334,   335,   336,   337,   338,   339,   340,   341,
     342,   343,   344,   345,   346,   347,   348,   349,   350,   351,
     352,   353,   354,   355,   356,   357,   358,   359,   360,   361,
     362,   363,   364,   365,   366,   367,   368,   369,   370,   371,
     372,   373,   374,   375,   376,   377,   378,   379,   380,   381,
     382,   383,   384,   385,   386,   387,   388,   389,   390,   391,
     392,   393,   394,   395,   396,   397,   398,   399,   400,   401,
     402,   403,   404,   405,   406,   407,   408,   409,   410,   411,
     412,   413,    -1,    -1,   416,     3,     4,     5,     6,     7,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    36,    37,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    45,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    55,    -1,    -1,
      -1,    -1,    60,    -1,    -1,    -1,    -1,    -1,    66,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    75,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    84,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   109,    -1,    -1,   112,    -1,    -1,    -1,    -1,   117,
     118,   119,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   173,    -1,   175,    -1,    -1,
      -1,    -1,    -1,    -1,   182,    -1,    -1,   185,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   194,   195,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   205,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   215,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   231,    -1,    -1,    -1,    -1,    -1,   237,
      -1,    -1,    -1,   241,    -1,    -1,   244,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    16,   256,    -1,
      -1,    -1,    21,    22,    23,    24,    -1,    26,    27,    28,
      29,    30,    -1,    32,    33,    -1,    35,    -1,    -1,    38,
      -1,    -1,    -1,    -1,    -1,    44,    -1,    -1,    -1,    -1,
      49,   289,    51,    -1,    53,    -1,    -1,    -1,    57,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    71,    72,    -1,    -1,    -1,    -1,    -1,    -1,
      79,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   331,    -1,    -1,    95,    -1,    -1,    -1,
      -1,   100,   101,    -1,   342,    -1,   105,   106,   107,   108,
      -1,    -1,   111,    -1,    -1,    -1,   115,    -1,    -1,    -1,
      -1,    -1,    -1,   122,   123,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   379,   380,    -1,    -1,   383,   145,    -1,    -1,    -1,
     388,    -1,    -1,   391,   392,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    31,    -1,    -1,    34,   165,    -1,    -1,    -1,
      39,    40,    41,    42,    43,    44,    -1,    46,    47,    48,
      -1,    50,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   198,
     199,    -1,    71,   202,    -1,   204,    -1,   206,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   214,    -1,   216,   217,    88,
     219,   220,   221,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   111,    -1,    -1,    -1,    -1,    -1,   247,    -1,
      -1,    -1,    -1,    -1,   123,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   136,    -1,   138,
      -1,   140,   141,    -1,    -1,    -1,    -1,    -1,    -1,   278,
      -1,   150,    -1,    -1,    -1,    -1,    -1,   286,    -1,   288,
      -1,   290,    -1,    -1,   163,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   301,    -1,    -1,    -1,    -1,   176,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   184,    -1,    -1,   187,    -1,
      -1,    -1,   191,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     329,   200,    -1,    -1,   203,   334,   335,   336,   337,   338,
      -1,    -1,    -1,   342,   343,    -1,    -1,    -1,   123,    -1,
      -1,    -1,    -1,    -1,    -1,   224,    -1,   132,   133,   134,
     135,   360,   137,    -1,   139,   234,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   377,   154,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   257,   258,
      -1,    -1,   167,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   401,    -1,    -1,   404,    -1,    -1,    -1,    -1,
      -1,    -1,   281,    -1,    -1,   284,    -1,    -1,    -1,   288,
      -1,    -1,   291,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   301,    -1,    -1,    -1,    -1,    -1,   307,    -1,
     309,    -1,   311,   312,    -1,    -1,    -1,   222,   223,    -1,
      -1,    -1,    -1,    -1,   323,    -1,    -1,    -1,    -1,   328,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   352,    -1,    -1,    -1,    -1,    -1,   358,
      -1,    -1,    -1,    -1,   363,    -1,    -1,    -1,    -1,   274,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   381,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   390,    -1,    -1,   299,    -1,    -1,   396,    -1,   304,
     305,   306,    -1,   308,    -1,   310,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   324,
     325,   326,   327,    -1,    -1,    -1,    -1,    -1,   333,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     345,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     355,    -1,   357,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   367,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   378
  };

  const short
  Parser::yystos_[] =
  {
       0,     1,   438,   444,     0,     3,     4,     5,     6,     7,
      36,    37,    45,    55,    60,    66,    75,    84,   109,   112,
     117,   118,   119,   173,   175,   182,   185,   194,   195,   205,
     215,   231,   237,   241,   244,   256,   289,   331,   342,   379,
     380,   383,   388,   391,   392,   439,   442,   443,   446,   447,
     448,   449,   450,   451,   455,   456,   459,   460,   554,   556,
     558,   559,   586,   588,   589,   607,   608,   613,   614,   621,
     622,   623,   644,   645,   661,   663,   775,   777,   795,   797,
     799,   810,   847,   848,   849,   850,   851,   863,   875,   876,
     877,   878,   879,   880,   881,   882,   800,   802,   801,   796,
     798,   414,   415,   436,   436,   436,   778,   417,   128,   207,
     452,   436,   436,   417,   436,   436,   436,   436,   436,   461,
     664,   436,   436,   436,   246,   249,   624,   436,   246,   249,
     436,   811,   646,   436,   245,   277,   453,   440,   441,   587,
     457,   464,   416,   557,   159,   185,   592,   598,   601,   610,
     616,   649,   667,   781,   123,   445,   416,   416,   416,   416,
     416,   428,   428,   428,   416,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,     3,     4,     5,     6,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,    60,    61,    62,    63,    64,    65,    66,    67,
      68,    69,    70,    71,    72,    73,    74,    75,    76,    77,
      78,    79,    80,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    91,    92,    93,    94,    95,    96,    97,
      98,    99,   100,   101,   102,   103,   104,   105,   106,   107,
     108,   109,   110,   111,   112,   113,   114,   115,   116,   117,
     118,   119,   120,   121,   122,   123,   124,   125,   126,   127,
     128,   129,   130,   131,   132,   133,   134,   135,   136,   137,
     138,   139,   140,   141,   142,   143,   144,   145,   146,   147,
     148,   149,   150,   151,   152,   153,   154,   155,   156,   157,
     158,   159,   160,   161,   162,   163,   164,   165,   166,   167,
     168,   169,   170,   171,   172,   173,   174,   175,   176,   177,
     178,   179,   180,   181,   182,   183,   184,   185,   186,   187,
     188,   189,   190,   191,   192,   193,   194,   195,   196,   197,
     198,   199,   200,   201,   202,   203,   204,   205,   206,   207,
     208,   209,   210,   211,   212,   213,   214,   215,   216,   217,
     218,   219,   220,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   236,   237,
     238,   239,   240,   241,   242,   243,   244,   245,   246,   247,
     248,   249,   250,   251,   252,   253,   254,   255,   256,   257,
     258,   259,   260,   261,   262,   263,   264,   265,   266,   267,
     268,   269,   270,   271,   272,   273,   274,   275,   276,   277,
     278,   279,   280,   281,   282,   283,   284,   285,   286,   287,
     288,   289,   290,   291,   292,   293,   294,   295,   296,   297,
     298,   299,   300,   301,   302,   303,   304,   305,   306,   307,
     308,   309,   310,   311,   312,   313,   314,   315,   316,   317,
     318,   319,   320,   321,   322,   323,   324,   325,   326,   327,
     328,   329,   330,   331,   332,   333,   334,   335,   336,   337,
     338,   339,   340,   341,   342,   343,   344,   345,   346,   347,
     348,   349,   350,   351,   352,   353,   354,   355,   356,   357,
     358,   359,   360,   361,   362,   363,   364,   365,   366,   367,
     368,   369,   370,   371,   372,   373,   374,   375,   376,   377,
     378,   379,   380,   381,   382,   383,   384,   385,   386,   387,
     388,   389,   390,   391,   392,   393,   394,   395,   396,   397,
     398,   399,   400,   401,   402,   403,   404,   405,   406,   407,
     408,   409,   410,   411,   412,   413,   416,   437,   416,   428,
     297,   428,   436,   428,   428,   416,   852,   428,   428,   428,
     812,   437,   436,   246,   249,   454,   416,   415,   416,    71,
      87,    94,   123,   151,   284,   301,   364,   397,   458,    16,
      21,    22,    23,    24,    26,    27,    28,    29,    30,    32,
      33,    35,    38,    44,    49,    51,    53,    57,    71,    72,
      79,    95,   100,   101,   105,   106,   107,   108,   111,   115,
     122,   123,   145,   165,   198,   199,   202,   204,   206,   214,
     216,   217,   219,   220,   221,   247,   278,   286,   288,   290,
     301,   329,   334,   335,   336,   337,   338,   342,   343,   360,
     377,   401,   404,   462,   465,    96,   160,   150,   185,   288,
     301,   369,   392,   560,   567,   568,   570,   574,   575,   578,
     579,    96,   591,   602,   598,   599,   603,   123,   322,   609,
     611,   612,   123,   354,   615,   617,   620,    74,   123,   319,
     332,   353,   647,   650,   651,   652,   655,    74,    77,    99,
     120,   155,   159,   161,   191,   245,   253,   277,   284,   288,
     331,   332,   340,   353,   365,   575,   662,   668,   671,   675,
     680,   681,   682,   683,   684,   685,   686,   688,   690,   691,
     693,   694,   695,   743,   744,   745,   752,   754,   755,    69,
      70,    97,   148,   158,   331,   371,   692,   776,   782,   788,
     799,   192,   431,   431,   431,   431,   431,   864,   555,   428,
     625,   428,    80,   123,   185,   192,   194,   237,   277,   391,
     392,   813,   428,   428,   428,   428,   276,   211,   210,   379,
     208,   212,   248,   232,   398,    58,   273,   313,   522,   436,
     429,   436,   436,   436,   477,   480,   428,   436,   476,   479,
     436,   436,   482,   436,   436,   436,   466,    86,   436,   292,
     429,   436,    58,   436,   436,   436,   436,   436,   436,   103,
     104,   168,   389,   542,   436,    13,    59,   506,   463,   436,
     436,   436,   436,   436,   436,   436,   436,   436,   436,   491,
     436,   436,   436,   436,   506,   470,   436,   320,   436,   436,
     436,   436,   436,   436,   436,   436,   254,   488,   436,    89,
     170,   196,   261,   314,   395,   541,   436,   436,   576,   580,
     571,   436,   561,   123,   584,   569,   114,   143,   149,   156,
     157,   225,   226,   227,   228,   229,   293,   294,   295,   296,
     321,   399,   428,   429,   436,   577,   660,   581,   590,   416,
     288,   416,   558,   593,   594,   595,   600,   111,   122,   209,
     260,   300,   301,   342,   401,   604,   342,   416,   182,   416,
     618,    81,   267,   395,   648,   656,   436,   653,    62,    65,
      81,    82,    85,   125,   238,   267,   303,   395,   676,   753,
     185,   747,   687,    67,   181,   416,   416,   689,   660,   696,
     436,   669,   436,    62,   159,   384,   672,   428,   577,   660,
     123,   665,   437,   734,   699,   123,   185,   271,   280,   300,
     391,   401,   720,   721,   730,   757,   785,   784,   436,   416,
     787,   786,   783,   123,   779,   789,   169,   423,   429,   436,
     806,   142,   169,   374,   413,   417,   429,   806,   807,   808,
     169,   429,   808,   808,   808,   116,   258,   866,   867,   868,
     416,   164,   628,   116,   258,   854,   855,   856,   815,   289,
     820,   814,   817,   821,   816,   818,   819,   436,   436,   436,
     436,   436,   436,   436,   436,   436,   471,   478,   660,   428,
     428,   428,   292,   436,   517,   517,   428,   517,   517,   428,
     428,   263,   264,   265,   266,   544,   428,   481,   428,   193,
     492,   292,   436,   428,   429,   436,   428,    90,   401,   436,
     428,   436,   428,   428,   428,   436,   428,   428,   428,   436,
     437,   428,   428,   428,   428,   428,   428,   428,   485,   401,
     428,   660,   727,   487,   428,   428,   436,   428,   436,   436,
     416,   533,   534,   188,   292,   436,   428,   428,   428,   428,
     428,   428,   428,   468,   402,   171,   270,   376,   496,   428,
     428,   428,   428,   416,   416,   416,   572,   573,   428,   416,
     585,   570,   436,   436,   428,   428,   577,   280,   300,   582,
     592,   428,   596,   416,   417,   436,   416,   123,   605,   595,
     168,   389,   436,   436,   436,   660,   436,   436,   436,   416,
     428,   436,   619,   428,   428,   428,   437,   657,    68,   294,
     405,   406,   428,   654,    61,   339,    25,   144,   341,   361,
     362,   400,   678,    65,    63,    64,   283,   285,   368,   370,
     679,    54,   172,   174,   255,   284,   341,   677,   428,   416,
     748,   746,   416,   428,   428,   416,   428,   416,   428,   437,
     428,   416,   670,   674,    68,   428,   428,   428,   294,   405,
     406,   428,   673,   428,   428,   577,   666,   428,   436,   428,
      31,    34,    39,    40,    41,    42,    43,    44,    46,    47,
      48,    50,    71,    88,   111,   123,   136,   138,   140,   141,
     163,   176,   184,   187,   191,   200,   203,   224,   234,   257,
     258,   281,   284,   288,   291,   301,   307,   309,   311,   312,
     323,   328,   352,   358,   363,   381,   390,   396,   575,   697,
     700,   714,   715,   722,   183,   727,   183,   727,   183,   660,
     183,   660,   436,   123,   724,   123,   132,   133,   134,   135,
     137,   139,   154,   167,   222,   223,   274,   299,   304,   305,
     306,   308,   310,   324,   325,   326,   327,   333,   345,   355,
     357,   367,   378,   756,   758,   762,   763,   771,   437,   734,
     734,   793,   405,   406,   742,   405,   406,   735,   734,   780,
      69,    70,   123,   790,   807,   806,   806,   423,   424,   425,
     426,   428,   432,   803,   807,   807,   806,   807,   808,    20,
     250,   407,   408,   409,   410,   411,   412,   431,   433,   434,
     809,    20,   250,   407,   408,   803,    20,   250,   424,   809,
     807,   808,   803,   803,   803,   436,   869,   123,   865,   867,
     416,   428,   629,   436,   857,   123,   853,   855,   416,   416,
     416,   416,   416,   416,   416,   416,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   151,   523,   110,   521,   660,
     429,   428,   428,   428,   428,   428,   521,   401,   493,   429,
     428,   429,   535,   536,   436,   436,   436,   428,   428,   428,
     436,   401,   543,   436,   729,   514,   428,   428,   436,   416,
     417,   436,   428,   534,   436,   429,   428,    73,   269,   322,
     832,   833,   436,   401,   436,   401,   428,   428,   416,   417,
     436,   428,   573,   428,   416,   436,   428,   583,   660,   593,
     595,   597,   428,   606,   428,   428,   436,   428,   428,   660,
     428,    68,   366,   436,   436,   416,   428,   658,   436,   428,
     416,   123,   747,   428,   428,   428,   416,   417,   436,   428,
     674,   436,   428,   416,   436,   436,   436,   436,   436,   436,
     436,   436,   707,   436,   436,   436,   436,   436,    17,   302,
     144,   172,   174,   255,   698,   436,   436,   436,   436,   710,
     705,   416,   436,   701,   436,   436,   702,   708,   703,   704,
     436,   706,   436,   436,   436,   436,   436,   436,   174,   255,
     347,   349,   718,    15,   144,   303,   719,   709,   416,   436,
      19,    76,    93,   162,   284,   330,   717,   436,   436,   351,
     428,   660,    74,   123,   716,   416,   727,   729,   727,   660,
     728,   660,   660,   660,   731,   428,   721,   365,   428,   436,
     436,   436,   436,   759,   436,   436,   436,   436,   417,   417,
     417,   327,   436,   769,   429,   767,   768,   760,   180,   239,
     240,   772,   356,    52,   233,   282,   764,    98,   179,   372,
     766,   428,   428,   123,   218,   794,   436,   436,   428,   436,
     436,   428,   428,   416,   792,   791,   416,   359,   432,   804,
     430,   806,   806,   806,   806,   804,   430,   430,   806,   806,
     806,   807,   807,   807,   807,   808,   808,   808,   808,   804,
     428,   436,   870,    84,   803,   428,   185,   213,   288,   385,
     386,   416,   556,   607,   626,   630,   631,   632,   633,   634,
     637,   428,   436,   858,   241,   803,   178,   230,   298,   350,
     822,   822,   822,   822,   822,   822,   822,   822,   436,   428,
     483,   660,   428,   436,    92,   429,   537,   538,   436,   430,
     536,   430,   472,   474,   489,   436,   428,   486,   428,   728,
     177,   189,   201,   259,   348,   428,   515,   516,   490,   401,
     429,   539,   540,   834,    18,    53,   126,   127,   185,   190,
     242,   297,   835,   832,   342,   436,   497,   436,    91,   430,
     727,   660,   605,   428,   595,   416,   428,   428,   436,   436,
     346,   428,   577,   428,   428,   428,   577,   185,   843,   843,
     185,   845,   845,   845,   843,   843,   263,   264,   265,   266,
     711,   843,   843,   843,   843,   428,   428,   428,   428,   428,
     428,   373,   428,   437,   428,   428,   428,   428,   416,   436,
     416,   428,   416,   428,   428,   416,   417,   436,   436,   428,
     416,   712,   713,   428,   428,   428,   428,   428,   428,   428,
     428,   416,   428,   428,   428,   428,   428,   428,   660,   428,
     577,   676,   720,   723,   729,   428,   728,   728,   660,   428,
     732,   416,   436,   436,   436,   436,   416,   773,   436,   436,
     436,   436,   428,   428,   428,   436,   428,   436,   436,   428,
     768,   416,   774,   428,   436,    52,   233,   282,   765,   436,
     378,    97,   436,   113,   113,   113,   113,   734,   734,   359,
     806,   806,   807,   808,   808,   428,   436,   638,   416,   635,
     416,   416,   416,   417,   436,   123,   627,   428,   436,   823,
     416,   824,   417,   428,   428,   428,   428,   428,   428,   428,
     428,   524,   519,   660,   436,   436,   430,   538,   436,   428,
     428,   532,   532,   188,   401,   507,   402,   511,   436,   436,
     401,   510,   436,   436,   430,   540,   275,   829,   436,   436,
     436,   436,   836,   436,   436,   436,   428,   436,   402,   532,
     287,   551,   436,   728,   428,   428,   428,   428,   659,   749,
     113,   844,   428,   428,   846,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   436,   428,   428,
     428,   428,   436,   436,   416,   417,   436,   428,   713,   428,
     428,   577,   428,   428,   123,   131,   725,   113,   733,   728,
     729,   733,   416,   428,   428,   428,   428,   428,   416,   428,
     428,   428,   428,   428,   167,   436,   416,   428,   436,   356,
     436,   772,   403,   436,   436,   436,   436,   428,   428,   121,
     432,   805,   805,   805,   393,   871,   872,   416,   436,   634,
     636,   428,   428,   416,   393,   859,   860,   297,   830,   830,
     532,   430,   520,   660,   518,    56,   436,   428,   430,   428,
     436,   428,   436,   436,   428,   436,   152,   153,   512,   436,
     428,   428,   436,   428,   402,   402,   243,   416,   297,   436,
     469,   436,   498,   436,   436,   436,   728,   300,   751,   436,
     416,   416,   428,   428,   428,   428,   102,   342,   726,   436,
     428,   729,   428,   428,   733,   436,   436,   436,   436,   761,
     356,   436,   348,   348,   348,   348,   121,   806,   807,   808,
     436,   872,   639,   428,   428,   634,   436,   860,   436,   436,
     831,   831,   428,   428,   519,   436,   430,   473,   475,   508,
     188,   513,   430,   436,   436,   436,   837,   436,   839,   402,
     504,   505,   342,   401,   501,   428,   729,   660,   750,    68,
     436,   436,   428,    68,   733,   428,   771,   430,   428,   436,
     436,   770,   436,   428,   436,   436,   436,   436,   428,   401,
     428,   436,   525,   484,   430,   342,   357,   530,   530,   129,
     509,   436,   428,   838,   840,   242,   346,   841,   436,   171,
     297,   382,   825,   436,   428,   504,   436,   436,   532,   186,
     428,   660,   751,   436,   436,   428,   428,   436,   436,   387,
     436,   736,   738,   873,   436,   861,    90,   401,   527,   436,
     436,   436,   402,   130,   842,   268,   826,   436,   436,   436,
     342,   503,   499,   502,   562,   436,   348,   348,   428,   436,
     436,   436,   185,   185,    83,   428,   394,   436,   436,   357,
     467,   531,   436,   436,   297,   436,   436,   552,   532,   401,
     548,   549,   416,   428,   436,   436,   436,   436,   428,   737,
     739,   436,   874,   640,   436,   862,   528,   529,   436,    56,
     494,   495,   532,   402,   436,   401,   553,   500,   436,   548,
     416,   436,   436,   436,   428,   740,   740,   428,   436,   642,
     428,   436,   532,   532,   526,   436,   428,   494,   428,   436,
     436,   436,   545,   551,   416,   436,   416,   741,    71,   108,
     115,   123,   301,   342,   404,   643,   428,   428,   532,   342,
     827,   402,   401,   546,   436,   428,   436,    86,   436,   436,
     641,   320,   436,   436,   428,   436,   375,   828,   436,   436,
     550,    92,   436,   428,   428,   416,   436,   428,   428,   342,
     547,   532,   436,   428,   428,   436,   532,   436,   428,   122,
     436,   436,   436,   436,   428,   563,   564,   247,   253,   272,
     316,   565,   436,   436,   566,   436,   436,   436,   416,   436,
     436,   428,   428,   428,   436,   428
  };

  const short
  Parser::yyr1_[] =
  {
       0,   435,   436,   436,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   438,   440,   439,   441,
     439,   442,   443,   444,   444,   444,   445,   445,   446,   446,
     446,   446,   446,   446,   446,   446,   446,   446,   446,   446,
     446,   446,   446,   446,   446,   446,   446,   446,   446,   446,
     446,   446,   446,   446,   446,   446,   446,   446,   446,   446,
     446,   446,   446,   446,   446,   446,   447,   447,   448,   448,
     449,   450,   451,   452,   452,   453,   453,   454,   454,   455,
     456,   457,   457,   458,   458,   458,   458,   458,   458,   458,
     458,   459,   461,   460,   463,   462,   464,   464,   466,   467,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   468,   465,   469,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   470,   465,   471,   465,   465,   465,   472,   473,   465,
     474,   475,   465,   465,   476,   465,   465,   477,   465,   478,
     465,   465,   479,   465,   465,   480,   465,   481,   465,   482,
     465,   465,   465,   465,   483,   484,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   485,   465,   486,   465,   487,   465,   465,   488,
     465,   489,   465,   490,   465,   465,   465,   491,   465,   492,
     492,   493,   493,   494,   494,   495,   497,   498,   499,   500,
     496,   501,   502,   496,   503,   496,   504,   504,   505,   506,
     506,   506,   507,   508,   507,   507,   509,   509,   510,   510,
     511,   511,   512,   512,   512,   513,   513,   514,   514,   515,
     515,   515,   516,   516,   516,   517,   518,   517,   519,   519,
     520,   521,   521,   522,   522,   522,   524,   525,   526,   523,
     527,   528,   527,   529,   527,   531,   530,   532,   532,   533,
     533,   534,   534,   534,   535,   535,   536,   537,   537,   538,
     539,   539,   540,   541,   541,   541,   541,   541,   541,   542,
     542,   542,   542,   543,   543,   544,   544,   544,   544,   545,
     545,   547,   546,   548,   548,   550,   549,   551,   551,   552,
     552,   553,   554,   555,   554,   557,   556,   558,   559,   559,
     559,   561,   562,   563,   560,   564,   564,   565,   565,   565,
     566,   565,   567,   567,   568,   569,   569,   570,   570,   570,
     571,   570,   570,   572,   572,   573,   573,   573,   574,   574,
     574,   574,   576,   575,   577,   577,   577,   577,   577,   577,
     577,   577,   577,   577,   577,   577,   577,   577,   577,   577,
     578,   580,   579,   581,   581,   582,   583,   582,   585,   584,
     587,   586,   588,   590,   589,   591,   591,   592,   592,   593,
     593,   594,   594,   596,   595,   597,   597,   595,   595,   595,
     598,   599,   599,   600,   602,   601,   603,   603,   604,   604,
     604,   604,   604,   604,   604,   604,   604,   606,   605,   607,
     608,   609,   610,   610,   611,   611,   612,   613,   614,   615,
     616,   616,   617,   618,   618,   619,   620,   621,   622,   624,
     625,   626,   623,   627,   627,   628,   628,   629,   629,   630,
     630,   630,   630,   630,   630,   630,   631,   632,   633,   635,
     634,   636,   636,   634,   634,   634,   638,   639,   640,   641,
     637,   642,   642,   643,   643,   643,   643,   643,   643,   644,
     646,   645,   648,   647,   649,   649,   650,   650,   650,   650,
     651,   651,   651,   652,   653,   653,   654,   654,   654,   656,
     655,   657,   657,   659,   658,   660,   660,   662,   661,   664,
     663,   666,   665,   667,   667,   668,   668,   668,   668,   668,
     668,   668,   668,   668,   668,   668,   668,   668,   668,   668,
     668,   668,   668,   668,   669,   668,   670,   670,   671,   672,
     672,   673,   673,   673,   674,   674,   674,   675,   676,   676,
     676,   676,   676,   676,   676,   676,   676,   676,   676,   676,
     676,   676,   676,   677,   677,   677,   677,   677,   677,   678,
     678,   678,   678,   678,   678,   679,   679,   679,   679,   679,
     679,   680,   681,   682,   682,   682,   683,   684,   685,   685,
     685,   685,   687,   686,   689,   688,   690,   690,   691,   692,
     693,   694,   696,   695,   698,   697,   699,   699,   700,   700,
     700,   700,   700,   700,   701,   700,   700,   700,   700,   700,
     700,   700,   700,   700,   700,   700,   700,   700,   700,   702,
     700,   703,   700,   704,   700,   705,   700,   700,   700,   700,
     700,   700,   700,   700,   700,   700,   700,   700,   706,   700,
     700,   700,   700,   700,   700,   700,   700,   700,   700,   700,
     700,   700,   700,   700,   700,   707,   700,   708,   700,   709,
     700,   710,   700,   711,   711,   711,   711,   712,   712,   713,
     713,   713,   714,   714,   714,   714,   714,   715,   716,   716,
     717,   717,   717,   717,   717,   717,   718,   718,   718,   718,
     719,   719,   719,   719,   720,   722,   723,   721,   721,   721,
     721,   721,   721,   721,   721,   721,   724,   724,   725,   725,
     726,   726,   726,   727,   728,   729,   729,   731,   730,   732,
     730,   733,   734,   734,   736,   737,   735,   738,   739,   735,
     735,   735,   740,   740,   741,   742,   742,   743,   743,   744,
     745,   746,   746,   748,   749,   747,   750,   750,   751,   753,
     752,   754,   755,   756,   757,   757,   759,   758,   760,   758,
     761,   758,   758,   758,   758,   758,   758,   758,   758,   758,
     758,   758,   758,   758,   758,   758,   758,   758,   758,   758,
     758,   762,   762,   762,   763,   763,   763,   763,   764,   764,
     764,   765,   765,   765,   766,   766,   767,   767,   768,   769,
     769,   770,   770,   770,   771,   771,   772,   772,   772,   773,
     773,   774,   774,   776,   775,   778,   777,   780,   779,   781,
     781,   783,   782,   784,   782,   785,   782,   786,   782,   782,
     787,   782,   782,   782,   788,   789,   789,   791,   790,   792,
     790,   793,   793,   794,   796,   795,   798,   797,   800,   799,
     801,   799,   802,   799,   803,   803,   803,   804,   804,   805,
     805,   806,   806,   806,   806,   806,   806,   806,   806,   807,
     807,   807,   807,   807,   807,   807,   807,   807,   807,   807,
     807,   807,   807,   807,   808,   808,   808,   808,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   811,   810,   812,
     812,   814,   813,   815,   813,   816,   813,   817,   813,   818,
     813,   819,   813,   820,   813,   821,   813,   823,   822,   824,
     822,   822,   822,   822,   825,   825,   825,   825,   825,   826,
     827,   826,   828,   828,   829,   829,   830,   830,   831,   831,
     832,   832,   833,   834,   833,   833,   835,   836,   837,   835,
     838,   835,   835,   839,   835,   835,   835,   840,   835,   835,
     835,   841,   841,   842,   842,   843,   844,   843,   846,   845,
     847,   848,   849,   850,   852,   851,   853,   854,   854,   855,
     855,   857,   856,   858,   858,   859,   859,   861,   860,   862,
     862,   864,   863,   865,   866,   866,   867,   867,   869,   868,
     870,   870,   871,   871,   873,   872,   874,   874,   875,   876,
     877,   878,   879,   880,   881,   881,   882
  };

  const signed char
  Parser::yyr2_[] =
  {
       0,     2,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     0,     4,     0,
       4,     3,     3,     0,     2,     1,     0,     2,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     3,     3,     3,
       3,     4,     3,     1,     1,     1,     1,     1,     1,     4,
       1,     0,     2,     4,     4,     4,     4,     4,     4,     4,
       4,     3,     0,     3,     0,     3,     0,     2,     0,     0,
      13,     3,     3,     4,     3,     4,     3,     4,     3,     3,
       3,     3,     0,     6,     0,     9,     3,     4,     7,     4,
       7,     3,     3,     3,     3,     3,     3,     3,     3,     6,
       6,     0,     4,     0,     4,     4,     4,     0,     0,     9,
       0,     0,     9,     3,     0,     4,     3,     0,     4,     0,
       5,     3,     0,     4,     3,     0,     4,     0,     5,     0,
       4,     2,     3,     3,     0,     0,     9,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     4,     3,     3,
       3,     3,     0,     5,     0,     9,     0,     5,     7,     0,
       4,     0,     7,     0,     7,     3,     3,     0,     5,     0,
       1,     0,     2,     0,     2,     4,     0,     0,     0,     0,
      11,     0,     0,     9,     0,     9,     0,     2,     4,     0,
       1,     1,     0,     0,     4,     2,     0,     2,     0,     2,
       0,     2,     0,     1,     1,     0,     4,     0,     2,     1,
       2,     2,     1,     1,     1,     1,     0,     7,     0,     2,
       1,     0,     1,     1,     1,     1,     0,     0,     0,    12,
       0,     0,     5,     0,     5,     0,     5,     0,     2,     1,
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
       2,     0,     4,     0,     2,     4,     0,     7,     0,     3,
       0,     3,     5,     0,     7,     0,     1,     1,     2,     0,
       1,     1,     2,     0,     4,     1,     2,     2,     2,     2,
       2,     0,     2,     3,     0,     4,     0,     2,     3,     3,
       4,     5,     4,     5,     3,     3,     3,     0,     3,     3,
       1,     2,     0,     2,     5,     6,     1,     3,     1,     2,
       0,     2,     3,     0,     2,     2,     2,     4,     3,     0,
       0,     0,     8,     1,     2,     0,     2,     0,     2,     1,
       1,     1,     1,     1,     1,     1,     3,     3,     4,     0,
       4,     1,     2,     2,     2,     2,     0,     0,     0,     0,
      12,     0,     2,     3,     3,     4,     4,     3,     3,     3,
       0,     3,     0,     3,     0,     2,     5,     1,     1,     1,
       3,     3,     3,     3,     0,     2,     1,     1,     1,     0,
       4,     0,     2,     0,     3,     2,     4,     0,     4,     0,
       3,     0,     3,     0,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     3,     3,
       1,     1,     1,     1,     0,     4,     1,     2,     3,     0,
       2,     1,     1,     1,     2,     2,     2,     3,     1,     2,
       1,     1,     2,     2,     1,     1,     1,     1,     2,     1,
       1,     2,     2,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     3,     4,     3,     3,     3,     3,     3,     2,     3,
       4,     3,     0,     4,     0,     4,     3,     3,     1,     1,
       5,     3,     0,     3,     0,     3,     0,     2,     2,     3,
       4,     3,     4,     5,     0,     4,     3,     1,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     0,
       4,     0,     5,     0,     5,     0,     5,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     4,     3,     0,     4,
       4,     2,     4,     4,     4,     3,     3,     4,     4,     4,
       4,     4,     4,     4,     4,     0,     4,     0,     4,     0,
       4,     0,     4,     1,     1,     1,     1,     1,     2,     2,
       2,     2,     3,     3,     4,     3,     3,     1,     0,     3,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       0,     1,     1,     1,     2,     0,     0,     7,     3,     4,
       6,     4,     6,     6,     8,     1,     0,     2,     0,     1,
       0,     2,     2,     1,     1,     0,     2,     0,     5,     0,
       7,     7,    11,     4,     0,     0,    10,     0,     0,    10,
       6,     6,     0,     2,     1,     6,     6,     3,     2,     1,
       4,     0,     2,     0,     0,     7,     0,     2,     5,     0,
       4,     3,     1,     2,     0,     2,     0,     4,     0,     4,
       0,    10,     9,     3,     3,     4,     4,     4,     4,     4,
       4,     4,     4,     3,     7,     8,     6,     3,     3,     3,
       2,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     2,     5,     1,
       2,     0,     4,     7,     1,     1,     1,     1,     1,     1,
       2,     1,     2,     0,     4,     0,     3,     0,     3,     0,
       2,     0,     4,     0,     4,     0,     4,     0,     4,     4,
       0,     4,     5,     1,     2,     0,     2,     0,     4,     0,
       4,     0,     2,     5,     0,     6,     0,     6,     0,     6,
       0,     6,     0,     6,     0,     1,     1,     1,     2,     1,
       2,     3,     3,     3,     3,     2,     3,     6,     1,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     2,
       3,     6,     1,     1,     3,     3,     6,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     0,     5,     0,
       2,     0,     5,     0,     5,     0,     5,     0,     5,     0,
       5,     0,     5,     0,     5,     0,     5,     0,     4,     0,
       4,     1,     2,     2,     0,     1,     2,     5,     3,     0,
       0,     6,     0,     1,     0,     1,     0,     3,     0,     1,
       0,     2,     1,     0,     3,     1,     0,     0,     0,     5,
       0,     6,     2,     0,     5,     2,     5,     0,     6,     2,
       6,     0,     1,     0,     1,     0,     0,     3,     0,     3,
       4,     3,     3,     3,     0,     7,     2,     1,     2,     3,
       1,     0,     5,     1,     2,     1,     2,     0,     7,     1,
       2,     0,     7,     2,     1,     2,     3,     1,     0,     5,
       1,     2,     1,     2,     0,     7,     1,     2,     3,     3,
       3,     3,     3,     3,     0,     1,     1
  };


#if LEFPARSERDEBUG || 1
  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a YYNTOKENS, nonterminals.
  const char*
  const Parser::yytname_[] =
  {
  "\"end of file\"", "error", "\"invalid token\"", "\"DEFINE\"",
  "\"DEFINEB\"", "\"DEFINES\"", "\"MESSAGE\"", "\"CREATEFILE\"",
  "\"OPENFILE\"", "\"CLOSEFILE\"", "\"WARNING\"", "\"ERROR\"",
  "\"FATALERROR\"", "\"ABOVE\"", "\"ABUT\"", "\"ABUTMENT\"",
  "\"ACCURRENTDENSITY\"", "\"ACTIVE\"", "\"ADJACENTCUTS\"", "\"ANALOG\"",
  "\"AND\"", "\"ANTENNAAREAFACTOR\"", "\"ANTENNAAREADIFFREDUCEPWL\"",
  "\"ANTENNAAREAMINUSDIFF\"", "\"ANTENNAAREARATIO\"", "\"ANTENNACELL\"",
  "\"ANTENNACUMAREARATIO\"", "\"ANTENNACUMDIFFAREARATIO\"",
  "\"ANTENNACUMDIFFSIDEAREARATIO\"", "\"ANTENNACUMROUTINGPLUSCUT\"",
  "\"ANTENNACUMSIDEAREARATIO\"", "\"ANTENNADIFFAREA\"",
  "\"ANTENNADIFFAREARATIO\"", "\"ANTENNADIFFSIDEAREARATIO\"",
  "\"ANTENNAGATEAREA\"", "\"ANTENNAGATEPLUSDIFF\"",
  "\"ANTENNAINOUTDIFFAREA\"", "\"ANTENNAINPUTGATEAREA\"",
  "\"ANTENNALENGTHFACTOR\"", "\"ANTENNAMAXAREACAR\"",
  "\"ANTENNAMAXCUTCAR\"", "\"ANTENNAMAXSIDEAREACAR\"",
  "\"ANTENNAMETALAREA\"", "\"ANTENNAMETALLENGTH\"", "\"ANTENNAMODEL\"",
  "\"ANTENNAOUTPUTDIFFAREA\"", "\"ANTENNAPARTIALCUTAREA\"",
  "\"ANTENNAPARTIALMETALAREA\"", "\"ANTENNAPARTIALMETALSIDEAREA\"",
  "\"ANTENNASIDEAREARATIO\"", "\"ANTENNASIZE\"",
  "\"ANTENNASIDEAREAFACTOR\"", "\"ANYEDGE\"", "\"AREA\"", "\"AREAIO\"",
  "\"ARRAY\"", "\"ARRAYCUTS\"", "\"ARRAYSPACING\"", "\"AVERAGE\"",
  "\"BELOW\"", "\"BEGINEXT\"", "\"BLACKBOX\"", "\"BLOCK\"",
  "\"BOTTOMLEFT\"", "\"BOTTOMRIGHT\"", "\"BUMP\"", "\"BUSBITCHARS\"",
  "\"BUFFER\"", "\"BY\"", "\"CANNOTOCCUPY\"", "\"CANPLACE\"",
  "\"CAPACITANCE\"", "\"CAPMULTIPLIER\"", "\"CENTERTOCENTER\"",
  "\"CLASS\"", "\"CLEARANCEMEASURE\"", "\"CLOCK\"", "\"CLOCKTYPE\"",
  "\"COLUMNMAJOR\"", "\"CURRENTDEN\"", "\"COMPONENTPIN\"", "\"CORE\"",
  "\"CORNER\"", "\"CORRECTIONFACTOR\"", "\"CORRECTIONTABLE\"", "\"COVER\"",
  "\"CPERSQDIST\"", "\"CURRENT\"", "\"CURRENTSOURCE\"", "\"CUT\"",
  "\"CUTAREA\"", "\"CUTSIZE\"", "\"CUTSPACING\"", "\"DATA\"",
  "\"DATABASE\"", "\"DCCURRENTDENSITY\"", "\"DEFAULT\"", "\"DEFAULTCAP\"",
  "\"DELAY\"", "\"DENSITY\"", "\"DENSITYCHECKSTEP\"",
  "\"DENSITYCHECKWINDOW\"", "\"DESIGNRULEWIDTH\"", "\"DIAG45\"",
  "\"DIAG135\"", "\"DIAGMINEDGELENGTH\"", "\"DIAGSPACING\"",
  "\"DIAGPITCH\"", "\"DIAGWIDTH\"", "\"DIELECTRIC\"", "\"DIFFUSEONLY\"",
  "\"DIRECTION\"", "\"DIVIDERCHAR\"", "\"DO\"", "\"E\"",
  "\"EDGECAPACITANCE\"", "\"EDGERATE\"", "\"EDGERATESCALEFACTOR\"",
  "\"EDGERATETHRESHOLD1\"", "\"EDGERATETHRESHOLD2\"", "\"EEQ\"",
  "\"ELSE\"", "\"ENCLOSURE\"", "\"END\"", "\"ENDEXT\"", "\"ENDCAP\"",
  "\"ENDOFLINE\"", "\"ENDOFNOTCHWIDTH\"", "\"EUCLIDEAN\"",
  "\"EXCEPTEXTRACUT\"", "\"EXCEPTSAMEPGNET\"", "\"EXCEPTPGNET\"",
  "\"EXTENSION\"", "\"FALL\"", "\"FALLCS\"", "\"FALLRS\"",
  "\"FALLSATCUR\"", "\"FALLSATT1\"", "\"FALLSLEWLIMIT\"", "\"FALLT0\"",
  "\"FALLTHRESH\"", "\"FALLVOLTAGETHRESHOLD\"", "\"FALSE\"", "\"FE\"",
  "\"FEEDTHRU\"", "\"FILLACTIVESPACING\"", "\"FIXED\"", "\"FLIP\"",
  "\"FLOORPLAN\"", "\"FN\"", "\"FOREIGN\"", "\"FREQUENCY\"",
  "\"FROMABOVE\"", "\"FROMBELOW\"", "\"FROMPIN\"", "\"FUNCTION\"",
  "\"FS\"", "\"FW\"", "\"GCELLGRID\"", "\"GENERATE\"", "\"GENERATED\"",
  "\"GENERATOR\"", "\"GROUND\"", "\"GROUNDSENSITIVITY\"",
  "\"HARDSPACING\"", "\"HEIGHT\"", "\"HISTORY\"", "\"HOLD\"",
  "\"HORIZONTAL\"", "\"IF\"", "\"IMPLANT\"", "\"INFLUENCE\"", "\"INOUT\"",
  "\"INOUTPINANTENNASIZE\"", "\"INPUT\"", "\"INPUTPINANTENNASIZE\"",
  "\"INPUTNOISEMARGIN\"", "\"INSIDECORNER\"", "\"INTEGER\"",
  "\"INTRINSIC\"", "\"INVERT\"", "\"INVERTER\"", "\"IRDROP\"",
  "\"ITERATE\"", "\"IV_TABLES\"", "\"LAYER\"", "\"LAYERS\"", "\"LEAKAGE\"",
  "\"LENGTH\"", "\"LENGTHSUM\"", "\"LENGTHTHRESHOLD\"", "\"LEQ\"",
  "\"LIBRARY\"", "\"LONGARRAY\"", "\"MACRO\"", "\"MANUFACTURINGGRID\"",
  "\"MASTERSLICE\"", "\"MATCH\"", "\"MAXADJACENTSLOTSPACING\"",
  "\"MAXCOAXIALSLOTSPACING\"", "\"MAXDELAY\"", "\"MAXEDGES\"",
  "\"MAXEDGESLOTSPACING\"", "\"MAXLOAD\"", "\"MAXIMUMDENSITY\"",
  "\"MAXVIASTACK\"", "\"MAXWIDTH\"", "\"MAXXY\"", "\"MEGAHERTZ\"",
  "\"METALOVERHANG\"", "\"MICRONS\"", "\"MILLIAMPS\"", "\"MILLIWATTS\"",
  "\"MINCUTS\"", "\"MINENCLOSEDAREA\"", "\"MINFEATURE\"", "\"MINIMUMCUT\"",
  "\"MINIMUMDENSITY\"", "\"MINPINS\"", "\"MINSIZE\"", "\"MINSTEP\"",
  "\"MINWIDTH\"", "\"MPWH\"", "\"MPWL\"", "\"MUSTJOIN\"", "\"MX\"",
  "\"MY\"", "\"MXR90\"", "\"MYR90\"", "\"N\"", "\"NAMEMAPSTRING\"",
  "\"NAMESCASESENSITIVE\"", "\"NANOSECONDS\"", "\"NEGEDGE\"",
  "\"NETEXPR\"", "\"NETS\"", "\"NEW\"", "\"NONDEFAULTRULE\"", "\"NONE\"",
  "\"NONINVERT\"", "\"NONUNATE\"", "\"NOISETABLE\"", "\"NOTCHLENGTH\"",
  "\"NOTCHSPACING\"", "\"NOWIREEXTENSIONATPIN\"", "\"OBS\"", "\"OFF\"",
  "\"OFFSET\"", "\"OHMS\"", "\"ON\"", "\"OR\"", "\"ORIENT\"",
  "\"ORIENTATION\"", "\"ORIGIN\"", "\"ORTHOGONAL\"", "\"OUTPUT\"",
  "\"OUTPUTPINANTENNASIZE\"", "\"OUTPUTNOISEMARGIN\"",
  "\"OUTPUTRESISTANCE\"", "\"OUTSIDECORNER\"", "\"OVERHANG\"",
  "\"OVERLAP\"", "\"OVERLAPS\"", "\"OXIDE1\"", "\"OXIDE2\"", "\"OXIDE3\"",
  "\"OXIDE4\"", "\"PAD\"", "\"PARALLELEDGE\"", "\"PARALLELOVERLAP\"",
  "\"PARALLELRUNLENGTH\"", "\"PATH\"", "\"PATTERN\"", "\"PEAK\"",
  "\"PERIOD\"", "\"PGONLY\"", "\"PICOFARADS\"", "\"PIN\"", "\"PITCH\"",
  "\"PLACED\"", "\"POLYGON\"", "\"PORT\"", "\"POSEDGE\"", "\"POST\"",
  "\"POWER\"", "\"PRE\"", "\"PREFERENCLOSURE\"", "\"PRL\"", "\"PROPERTY\"",
  "\"PROPERTYDEFINITIONS\"", "\"PROTRUSIONWIDTH\"", "\"PULLDOWNRES\"",
  "\"PWL\"", "\"R0\"", "\"R90\"", "\"R180\"", "\"R270\"", "\"RANGE\"",
  "\"REAL\"", "\"RECOVERY\"", "\"RECT\"", "\"RESISTANCE\"",
  "\"RESISTIVE\"", "\"RING\"", "\"RISE\"", "\"RISECS\"", "\"RISERS\"",
  "\"RISESATCUR\"", "\"RISESATT1\"", "\"RISESLEWLIMIT\"", "\"RISET0\"",
  "\"RISETHRESH\"", "\"RISEVOLTAGETHRESHOLD\"", "\"RMS\"", "\"ROUTING\"",
  "\"ROWABUTSPACING\"", "\"ROWCOL\"", "\"ROWMAJOR\"", "\"ROWMINSPACING\"",
  "\"ROWPATTERN\"", "\"RPERSQ\"", "\"S\"", "\"SAMENET\"", "\"SCANUSE\"",
  "\"SDFCOND\"", "\"SDFCONDEND\"", "\"SDFCONDSTART\"", "\"SETUP\"",
  "\"SHAPE\"", "\"SHRINKAGE\"", "\"SIGNAL\"", "\"SITE\"", "\"SIZE\"",
  "\"SKEW\"", "\"SLOTLENGTH\"", "\"SLOTWIDTH\"", "\"SLOTWIRELENGTH\"",
  "\"SLOTWIREWIDTH\"", "\"SPLITWIREWIDTH\"", "\"SOFT\"", "\"SOURCE\"",
  "\"SPACER\"", "\"SPACING\"", "\"SPACINGTABLE\"", "\"SPECIALNETS\"",
  "\"STABLE\"", "\"STACK\"", "\"START\"", "\"STEP\"", "\"STOP\"",
  "\"STRING\"", "\"STRUCTURE\"", "\"SUPPLYSENSITIVITY\"", "\"SYMMETRY\"",
  "\"TABLE\"", "\"TABLEAXIS\"", "\"TABLEDIMENSION\"", "\"TABLEENTRIES\"",
  "\"TAPERRULE\"", "\"THEN\"", "\"THICKNESS\"", "\"TIEHIGH\"",
  "\"TIELOW\"", "\"TIEOFFR\"", "\"TIME\"", "\"TIMING\"", "\"TO\"",
  "\"TOPIN\"", "\"TOPLEFT\"", "\"TOPOFSTACKONLY\"", "\"TOPRIGHT\"",
  "\"TRACKS\"", "\"TRANSITIONTIME\"", "\"TRISTATE\"", "\"TRUE\"",
  "\"TWOEDGES\"", "\"TWOWIDTHS\"", "\"TYPE\"", "\"UNATENESS\"",
  "\"UNITS\"", "\"UNIVERSALNOISEMARGIN\"", "\"USE\"",
  "\"USELENGTHTHRESHOLD\"", "\"USEMINSPACING\"", "\"USER\"", "\"USEVIA\"",
  "\"USEVIARULE\"", "\"VARIABLE\"", "\"VERSION\"", "\"VERTICAL\"",
  "\"VHI\"", "\"VIA\"", "\"VIARULE\"", "\"VICTIMLENGTH\"",
  "\"VICTIMNOISE\"", "\"VIRTUAL\"", "\"VLO\"", "\"VOLTAGE\"", "\"VOLTS\"",
  "\"W\"", "\"WELLTAP\"", "\"WIDTH\"", "\"WITHIN\"", "\"WIRECAP\"",
  "\"WIREEXTENSION\"", "\"X\"", "\"Y\"", "K_EQ", "K_NE", "K_LE", "K_LT",
  "K_GE", "K_GT", "K_NOT", "\"integer\"", "\"double\"", "\"string\"",
  "\"qstring\"", "\"binary numbers\"", "\"number\"",
  "\"generalized string\"", "IF", "LNOT", "'-'", "'+'", "'*'", "'/'",
  "UMINUS", "';'", "'('", "')'", "'='", "'\\n'", "'<'", "'>'", "$accept",
  "NUMBER", "GSTRING", "lef_file", "version", "$@1", "$@2", "dividerchar",
  "busbitchars", "rules", "end_library", "rule", "case_sensitivity",
  "wireextension", "manufacturing", "useminspacing", "clearancemeasure",
  "clearance_type", "spacing_type", "spacing_value", "units_section",
  "start_units", "units_rules", "units_rule", "layer_rule", "start_layer",
  "$@3", "end_layer", "$@4", "layer_options", "layer_option", "$@5", "$@6",
  "$@7", "$@8", "$@9", "$@10", "$@11", "$@12", "$@13", "$@14", "$@15",
  "$@16", "$@17", "$@18", "$@19", "$@20", "$@21", "$@22", "$@23", "$@24",
  "$@25", "$@26", "$@27", "$@28", "$@29", "$@30",
  "layer_arraySpacing_long", "layer_arraySpacing_width",
  "layer_arraySpacing_arraycuts", "layer_arraySpacing_arraycut",
  "sp_options", "$@31", "$@32", "$@33", "$@34", "$@35", "$@36", "$@37",
  "layer_spacingtable_opts", "layer_spacingtable_opt",
  "layer_enclosure_type_opt", "layer_enclosure_width_opt", "$@38",
  "layer_enclosure_width_except_opt", "layer_preferenclosure_width_opt",
  "layer_minimumcut_within", "layer_minimumcut_from",
  "layer_minimumcut_length", "layer_minstep_options",
  "layer_minstep_option", "layer_minstep_type", "layer_antenna_pwl",
  "$@39", "layer_diffusion_ratios", "layer_diffusion_ratio",
  "layer_antenna_duo", "layer_table_type", "layer_frequency", "$@40",
  "$@41", "$@42", "ac_layer_table_opt", "$@43", "$@44", "dc_layer_table",
  "$@45", "number_list", "layer_prop_list", "layer_prop",
  "current_density_pwl_list", "current_density_pwl", "cap_points",
  "cap_point", "res_points", "res_point", "layer_type", "layer_direction",
  "layer_minen_width", "layer_oxide", "layer_sp_parallel_widths",
  "layer_sp_parallel_width", "$@46", "layer_sp_TwoWidths",
  "layer_sp_TwoWidth", "$@47", "layer_sp_TwoWidthsPRL",
  "layer_sp_influence_widths", "layer_sp_influence_width", "maxstack_via",
  "$@48", "via", "$@49", "via_keyword", "start_via", "via_viarule", "$@50",
  "$@51", "$@52", "via_viarule_options", "via_viarule_option", "$@53",
  "via_option", "via_other_options", "via_more_options",
  "via_other_option", "$@54", "via_prop_list", "via_name_value_pair",
  "via_foreign", "start_foreign", "$@55", "orientation", "via_layer_rule",
  "via_layer", "$@56", "via_geometries", "via_geometry", "$@57", "end_via",
  "$@58", "viarule_keyword", "$@59", "viarule", "viarule_generate", "$@60",
  "viarule_generate_default", "viarule_layer_list", "opt_viarule_props",
  "viarule_props", "viarule_prop", "$@61", "viarule_prop_list",
  "viarule_layer", "via_names", "via_name", "viarule_layer_name", "$@62",
  "viarule_layer_options", "viarule_layer_option", "end_viarule", "$@63",
  "spacing_rule", "start_spacing", "end_spacing", "spacings", "spacing",
  "samenet_keyword", "irdrop", "start_irdrop", "end_irdrop", "ir_tables",
  "ir_table", "ir_table_values", "ir_table_value", "ir_tablename",
  "minfeature", "dielectric", "nondefault_rule", "$@64", "$@65", "$@66",
  "end_nd_rule", "nd_hardspacing", "nd_rules", "nd_rule", "usevia",
  "useviarule", "mincuts", "nd_prop", "$@67", "nd_prop_list", "nd_layer",
  "$@68", "$@69", "$@70", "$@71", "nd_layer_stmts", "nd_layer_stmt",
  "site", "start_site", "$@72", "end_site", "$@73", "site_options",
  "site_option", "site_class", "site_symmetry_statement",
  "site_symmetries", "site_symmetry", "site_rowpattern_statement", "$@74",
  "site_rowpatterns", "site_rowpattern", "$@75", "pt", "macro", "$@76",
  "start_macro", "$@77", "end_macro", "$@78", "macro_options",
  "macro_option", "$@79", "macro_prop_list", "macro_symmetry_statement",
  "macro_symmetries", "macro_symmetry", "macro_name_value_pair",
  "macro_class", "class_type", "pad_type", "core_type", "endcap_type",
  "macro_generator", "macro_generate", "macro_source", "macro_power",
  "macro_origin", "macro_foreign", "macro_eeq", "$@80", "macro_leq",
  "$@81", "macro_site", "macro_site_word", "site_word", "macro_size",
  "macro_pin", "start_macro_pin", "$@82", "end_macro_pin", "$@83",
  "macro_pin_options", "macro_pin_option", "$@84", "$@85", "$@86", "$@87",
  "$@88", "$@89", "$@90", "$@91", "$@92", "$@93", "pin_layer_oxide",
  "pin_prop_list", "pin_name_value_pair", "electrical_direction",
  "start_macro_port", "macro_port_class_option", "macro_pin_use",
  "macro_scan_use", "pin_shape", "geometries", "geometry", "$@94", "$@95",
  "geometry_options", "layer_exceptpgnet", "layer_spacing", "firstPt",
  "nextPt", "otherPts", "via_placement", "$@96", "$@97", "stepPattern",
  "sitePattern", "trackPattern", "$@98", "$@99", "$@100", "$@101",
  "trackLayers", "layer_name", "gcellPattern", "macro_obs",
  "start_macro_obs", "macro_density", "density_layers", "density_layer",
  "$@102", "$@103", "density_layer_rects", "density_layer_rect",
  "macro_clocktype", "$@104", "timing", "start_timing", "end_timing",
  "timing_options", "timing_option", "$@105", "$@106", "$@107",
  "one_pin_trigger", "two_pin_trigger", "from_pin_trigger",
  "to_pin_trigger", "delay_or_transition", "list_of_table_entries",
  "table_entry", "list_of_table_axis_numbers", "slew_spec", "risefall",
  "unateness", "list_of_from_strings", "list_of_to_strings", "array",
  "$@108", "start_array", "$@109", "end_array", "$@110", "array_rules",
  "array_rule", "$@111", "$@112", "$@113", "$@114", "$@115",
  "floorplan_start", "floorplan_list", "floorplan_element", "$@116",
  "$@117", "cap_list", "one_cap", "msg_statement", "$@118",
  "create_file_statement", "$@119", "def_statement", "$@120", "$@121",
  "$@122", "dtrm", "then", "else", "expression", "b_expr", "s_expr",
  "relop", "prop_def_section", "$@123", "prop_stmts", "prop_stmt", "$@124",
  "$@125", "$@126", "$@127", "$@128", "$@129", "$@130", "$@131",
  "prop_define", "$@132", "$@133", "opt_range_second", "opt_endofline",
  "$@134", "opt_endofline_twoedges", "opt_samenetPGonly", "opt_def_range",
  "opt_def_value", "layer_spacing_opts", "layer_spacing_opt", "$@135",
  "layer_spacing_cut_routing", "$@136", "$@137", "$@138", "$@139", "$@140",
  "spacing_cut_layer_opt", "opt_adjacentcuts_exceptsame", "opt_layer_name",
  "$@141", "req_layer_name", "$@142", "universalnoisemargin",
  "edgeratethreshold1", "edgeratethreshold2", "edgeratescalefactor",
  "noisetable", "$@143", "end_noisetable", "noise_table_list",
  "noise_table_entry", "output_resistance_entry", "$@144", "num_list",
  "victim_list", "victim", "$@145", "vnoiselist", "correctiontable",
  "$@146", "end_correctiontable", "correction_table_list",
  "correction_table_item", "output_list", "$@147", "numo_list",
  "corr_victim_list", "corr_victim", "$@148", "corr_list", "input_antenna",
  "output_antenna", "inout_antenna", "antenna_input", "antenna_inout",
  "antenna_output", "extension_opt", "extension", YY_NULLPTR
  };
#endif


#if LEFPARSERDEBUG
  const short
  Parser::yyrline_[] =
  {
       0,   543,   543,   544,   546,   547,   548,   549,   550,   551,
     552,   553,   554,   555,   556,   557,   558,   559,   560,   561,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,   578,   579,   580,   581,
     582,   583,   584,   585,   586,   587,   588,   589,   590,   591,
     592,   593,   594,   595,   596,   597,   598,   599,   600,   601,
     602,   603,   604,   605,   606,   607,   608,   609,   610,   611,
     612,   613,   614,   615,   616,   617,   618,   619,   620,   621,
     622,   623,   624,   625,   626,   627,   628,   629,   630,   631,
     632,   633,   634,   635,   636,   637,   638,   639,   640,   641,
     642,   643,   644,   645,   646,   647,   648,   649,   650,   651,
     652,   653,   654,   655,   656,   657,   658,   659,   660,   661,
     662,   663,   664,   665,   666,   667,   668,   669,   670,   671,
     672,   673,   674,   675,   676,   677,   678,   679,   680,   681,
     682,   683,   684,   685,   686,   687,   688,   689,   690,   691,
     692,   693,   694,   695,   696,   697,   698,   699,   700,   701,
     702,   703,   704,   705,   706,   707,   708,   709,   710,   711,
     712,   713,   714,   715,   716,   717,   718,   719,   720,   721,
     722,   723,   724,   725,   726,   727,   728,   729,   730,   731,
     732,   733,   734,   735,   736,   737,   738,   739,   740,   741,
     742,   743,   744,   745,   746,   747,   748,   749,   750,   751,
     752,   753,   754,   755,   756,   757,   758,   759,   760,   761,
     762,   763,   764,   765,   766,   767,   768,   769,   770,   771,
     772,   773,   774,   775,   776,   777,   778,   779,   780,   781,
     782,   783,   784,   785,   786,   787,   788,   789,   790,   791,
     792,   793,   794,   795,   796,   797,   798,   799,   800,   801,
     802,   803,   804,   805,   806,   807,   808,   809,   810,   811,
     812,   813,   814,   815,   816,   817,   818,   819,   820,   821,
     822,   823,   824,   825,   826,   827,   828,   829,   830,   831,
     832,   833,   834,   835,   836,   837,   838,   839,   840,   841,
     842,   843,   844,   845,   846,   847,   848,   849,   850,   851,
     852,   853,   854,   855,   856,   857,   858,   859,   860,   861,
     862,   863,   864,   865,   866,   867,   868,   869,   870,   871,
     872,   873,   874,   875,   876,   877,   878,   879,   880,   881,
     882,   883,   884,   885,   886,   887,   888,   889,   890,   891,
     892,   893,   894,   895,   896,   897,   898,   899,   900,   901,
     902,   903,   904,   905,   906,   907,   908,   909,   910,   911,
     912,   913,   914,   915,   916,   917,   918,   919,   920,   921,
     922,   923,   924,   925,   926,   927,   928,   929,   930,   931,
     932,   933,   934,   935,   936,   937,   938,   939,   940,   941,
     942,   943,   944,   945,   946,   947,   948,   949,   950,   951,
     952,   953,   954,   955,   956,   957,   959,   978,   978,  1009,
    1009,  1043,  1058,  1073,  1074,  1075,  1079,  1085,  1113,  1113,
    1113,  1113,  1114,  1114,  1114,  1114,  1114,  1115,  1115,  1116,
    1116,  1116,  1116,  1116,  1116,  1116,  1117,  1117,  1117,  1118,
    1118,  1119,  1119,  1120,  1120,  1120,  1121,  1121,  1122,  1122,
    1122,  1122,  1123,  1123,  1123,  1124,  1127,  1137,  1154,  1165,
    1177,  1183,  1198,  1202,  1203,  1206,  1207,  1210,  1211,  1213,
    1216,  1237,  1238,  1241,  1243,  1245,  1247,  1249,  1251,  1253,
    1259,  1262,  1266,  1266,  1293,  1293,  1334,  1335,  1343,  1349,
    1342,  1365,  1371,  1376,  1381,  1385,  1389,  1393,  1397,  1401,
    1405,  1410,  1415,  1414,  1450,  1449,  1466,  1481,  1493,  1504,
    1516,  1527,  1539,  1551,  1563,  1575,  1587,  1612,  1639,  1658,
    1676,  1695,  1695,  1702,  1701,  1715,  1731,  1747,  1770,  1746,
    1773,  1796,  1772,  1800,  1839,  1838,  1877,  1916,  1915,  1955,
    1954,  1969,  2008,  2007,  2046,  2085,  2084,  2124,  2123,  2164,
    2163,  2202,  2223,  2244,  2266,  2284,  2265,  2298,  2322,  2346,
    2370,  2394,  2418,  2442,  2466,  2489,  2507,  2525,  2543,  2561,
    2579,  2604,  2630,  2629,  2647,  2646,  2662,  2661,  2668,  2685,
    2684,  2711,  2710,  2731,  2730,  2749,  2767,  2793,  2792,  2811,
    2813,  2819,  2821,  2827,  2829,  2832,  2849,  2870,  2875,  2881,
    2848,  2895,  2899,  2894,  2931,  2930,  2955,  2957,  2960,  2967,
    2968,  2969,  2971,  2973,  2972,  2979,  2995,  2996,  3012,  3013,
    3020,  3021,  3037,  3038,  3057,  3076,  3077,  3095,  3096,  3102,
    3107,  3111,  3125,  3126,  3127,  3130,  3134,  3133,  3147,  3148,
    3155,  3161,  3162,  3196,  3197,  3198,  3202,  3204,  3207,  3201,
    3211,  3213,  3212,  3227,  3226,  3243,  3242,  3247,  3248,  3252,
    3253,  3257,  3266,  3275,  3288,  3290,  3293,  3297,  3298,  3301,
    3305,  3306,  3309,  3313,  3314,  3315,  3316,  3317,  3318,  3321,
    3322,  3323,  3324,  3326,  3327,  3334,  3339,  3344,  3349,  3356,
    3357,  3364,  3363,  3373,  3374,  3389,  3388,  3399,  3400,  3404,
    3405,  3408,  3411,  3447,  3447,  3475,  3475,  3481,  3484,  3494,
    3503,  3513,  3515,  3518,  3513,  3540,  3541,  3544,  3548,  3552,
    3556,  3556,  3563,  3564,  3566,  3569,  3570,  3574,  3576,  3578,
    3580,  3580,  3585,  3596,  3597,  3601,  3612,  3621,  3632,  3642,
    3652,  3662,  3673,  3673,  3677,  3678,  3679,  3680,  3681,  3682,
    3683,  3684,  3685,  3686,  3687,  3688,  3689,  3690,  3691,  3692,
    3694,  3697,  3697,  3705,  3707,  3711,  3717,  3716,  3733,  3733,
    3764,  3764,  3775,  3802,  3801,  3837,  3838,  3856,  3857,  3860,
    3862,  3866,  3867,  3870,  3870,  3877,  3878,  3882,  3891,  3900,
    3912,  3943,  3945,  3948,  3951,  3951,  3958,  3960,  3964,  3983,
    4002,  4038,  4040,  4045,  4047,  4049,  4079,  4107,  4107,  4124,
    4127,  4149,  4162,  4164,  4167,  4179,  4192,  4196,  4199,  4210,
    4219,  4221,  4224,  4232,  4234,  4237,  4240,  4243,  4257,  4268,
    4269,  4280,  4268,  4310,  4315,  4333,  4335,  4354,  4355,  4359,
    4360,  4361,  4362,  4363,  4364,  4365,  4368,  4386,  4406,  4426,
    4426,  4433,  4434,  4438,  4447,  4456,  4468,  4469,  4479,  4483,
    4468,  4522,  4524,  4528,  4533,  4536,  4560,  4583,  4606,  4625,
    4631,  4631,  4642,  4642,  4670,  4672,  4676,  4683,  4685,  4691,
    4695,  4696,  4697,  4699,  4702,  4704,  4708,  4710,  4712,  4715,
    4715,  4719,  4721,  4724,  4724,  4728,  4730,  4734,  4733,  4741,
    4741,  4760,  4760,  4779,  4781,  4788,  4789,  4790,  4791,  4792,
    4793,  4795,  4797,  4799,  4800,  4801,  4803,  4805,  4807,  4809,
    4811,  4813,  4815,  4817,  4819,  4819,  4826,  4827,  4830,  4841,
    4843,  4847,  4849,  4851,  4855,  4866,  4875,  4885,  4894,  4895,
    4914,  4915,  4916,  4935,  4954,  4955,  4967,  4968,  4969,  4993,
    4994,  5001,  5004,  5010,  5011,  5012,  5013,  5014,  5015,  5018,
    5019,  5020,  5021,  5040,  5059,  5080,  5081,  5082,  5083,  5084,
    5085,  5087,  5090,  5094,  5103,  5112,  5122,  5132,  5183,  5188,
    5193,  5198,  5204,  5204,  5207,  5207,  5219,  5232,  5247,  5251,
    5254,  5275,  5283,  5283,  5290,  5290,  5309,  5310,  5314,  5324,
    5334,  5344,  5354,  5364,  5374,  5374,  5384,  5393,  5395,  5397,
    5399,  5408,  5417,  5426,  5435,  5444,  5453,  5462,  5471,  5473,
    5473,  5475,  5475,  5484,  5484,  5493,  5493,  5502,  5511,  5513,
    5515,  5524,  5533,  5542,  5551,  5560,  5569,  5579,  5581,  5581,
    5586,  5600,  5615,  5637,  5659,  5681,  5683,  5685,  5716,  5747,
    5778,  5809,  5840,  5871,  5902,  5934,  5933,  5963,  5963,  5980,
    5980,  5997,  5997,  6016,  6021,  6026,  6031,  6038,  6039,  6043,
    6054,  6063,  6074,  6075,  6076,  6077,  6078,  6080,  6092,  6093,
    6100,  6101,  6102,  6103,  6104,  6105,  6108,  6109,  6110,  6111,
    6114,  6115,  6116,  6117,  6119,  6122,  6123,  6122,  6140,  6152,
    6165,  6178,  6191,  6205,  6219,  6232,  6235,  6236,  6238,  6239,
    6254,  6255,  6269,  6284,  6290,  6296,  6301,  6305,  6305,  6310,
    6310,  6318,  6322,  6334,  6348,  6356,  6347,  6359,  6367,  6358,
    6369,  6378,  6388,  6390,  6393,  6396,  6405,  6415,  6426,  6445,
    6457,  6478,  6479,  6482,  6483,  6482,  6490,  6491,  6494,  6501,
    6501,  6504,  6507,  6510,  6524,  6526,  6531,  6530,  6541,  6541,
    6544,  6543,  6547,  6557,  6559,  6561,  6563,  6565,  6567,  6569,
    6571,  6573,  6575,  6577,  6579,  6581,  6583,  6585,  6587,  6589,
    6591,  6595,  6597,  6599,  6603,  6605,  6607,  6609,  6613,  6615,
    6617,  6621,  6623,  6625,  6629,  6631,  6635,  6637,  6640,  6644,
    6646,  6651,  6652,  6654,  6659,  6661,  6665,  6667,  6669,  6673,
    6675,  6679,  6681,  6685,  6684,  6693,  6693,  6704,  6704,  6725,
    6726,  6730,  6730,  6737,  6737,  6744,  6744,  6751,  6751,  6757,
    6761,  6761,  6767,  6773,  6776,  6781,  6782,  6786,  6786,  6793,
    6793,  6803,  6804,  6807,  6811,  6811,  6815,  6815,  6821,  6821,
    6831,  6831,  6841,  6841,  6853,  6854,  6855,  6858,  6859,  6863,
    6864,  6868,  6869,  6870,  6871,  6872,  6873,  6874,  6876,  6879,
    6880,  6881,  6882,  6883,  6884,  6885,  6886,  6887,  6888,  6889,
    6890,  6891,  6893,  6894,  6897,  6902,  6904,  6913,  6917,  6918,
    6919,  6920,  6921,  6922,  6923,  6924,  6925,  6929,  6928,  6944,
    6945,  6949,  6949,  6959,  6959,  6969,  6969,  6980,  6980,  6990,
    6990,  7000,  7000,  7010,  7010,  7020,  7020,  7032,  7032,  7037,
    7037,  7043,  7048,  7054,  7063,  7064,  7069,  7076,  7083,  7091,
    7093,  7092,  7101,  7102,  7110,  7111,  7119,  7120,  7125,  7126,
    7129,  7131,  7133,  7158,  7157,  7182,  7206,  7208,  7209,  7208,
    7226,  7225,  7250,  7274,  7273,  7279,  7285,  7293,  7292,  7308,
    7322,  7339,  7340,  7348,  7349,  7366,  7367,  7367,  7372,  7372,
    7376,  7390,  7402,  7414,  7427,  7426,  7432,  7445,  7446,  7450,
    7457,  7461,  7460,  7466,  7469,  7474,  7475,  7479,  7478,  7485,
    7488,  7493,  7492,  7499,  7511,  7512,  7516,  7523,  7527,  7526,
    7533,  7536,  7541,  7542,  7547,  7546,  7553,  7556,  7562,  7586,
    7610,  7634,  7668,  7702,  7736,  7737,  7739
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
#endif // LEFPARSERDEBUG

  Parser::symbol_kind_type
  Parser::yytranslate_ (int t) YY_NOEXCEPT
  {
    // YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to
    // TOKEN-NUM as returned by yylex.
    static
    const short
    translate_table[] =
    {
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     432,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     429,   430,   425,   424,     2,   423,     2,   426,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,   428,
     433,   431,   434,     2,     2,     2,     2,     2,     2,     2,
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
     415,   416,   417,   418,   419,   420,   421,   422,   427
    };
    // Last valid token kind.
    const int code_max = 678;

    if (t <= 0)
      return symbol_kind::S_YYEOF;
    else if (t <= code_max)
      return static_cast <symbol_kind_type> (translate_table[t]);
    else
      return symbol_kind::S_YYUNDEF;
  }

} // LefParser
#line 15463 "/home/fantastic_na/OpenMPL/OpenMPL/build/thirdparty/Limbo/limbo/parsers/lef/bison/LefParser.cc"

#line 7749 "LefParser.yy"
 /*** Additional Code ***/

void LefParser::Parser::error(const Parser::location_type& l,
			    const std::string& m)
{
    driver.error(l, m);
}

