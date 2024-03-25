#ifndef INFRA_ARGDESCBASE_HH
#define INFRA_ARGDESCBASE_HH

#include "util.hh"

#include <stdlib.h>
#include <string.h>


typedef unsigned int uint;
typedef const char* constcharptr_t;

template<class Targ>
class argdescbase_t {
  public:
    typedef std::vector<argdescbase_t*> argdesc_vt;
  public:
    argdescbase_t(const char* aName, const char* aDesc) : _name(aName), _desc(aDesc) {}
    virtual ~argdescbase_t(){}
  public:
    inline const char* name() const { return _name; }
    inline const char* desc() const { return _desc; }
  public:
    virtual std::ostream& print_usage(std::ostream& os) const = 0;
    virtual int           parse_args(int& currArg, const int noArgs, const char* argv[], Targ&) const = 0;
    virtual int           parse_args2(uint& currArg, const std::vector<char*>& argv, Targ&) const = 0;
  protected:
    const char* _name; // name of the argument (e.g. --option)
    const char* _desc; // verbal description of argument
  private:
    argdescbase_t(const argdescbase_t&);
    argdescbase_t& operator=(const argdescbase_t&);
};

template<class Targ, class Tval>
class argdesc_t : public argdescbase_t<Targ> {
  public:
    typedef Targ arg_t;
    typedef void (arg_t::*setfun_t)(const Tval&);
  public:
    argdesc_t(const char* aName,
              const Tval& aDefault,
              setfun_t    aSetFun, // member function of Targ to set argument/option
              const char* aDesc) : argdescbase_t<Targ>(aName, aDesc), _default(aDefault), _setfun(aSetFun) {}
    virtual ~argdesc_t() {}
  public:
    virtual std::ostream& print_usage(std::ostream& os) const {
              return os;
            }
    // return   0 if argument name does not fit, other argument should be tested
    // returns -1 if argument name fits but error occurred
    // returns  1 if argument name fits and no error occurred
    virtual int parse_args(int& currArg, const int noArgs, const char* argv[], Targ& aArg) const {
              return true;
            }
    virtual int parse_args2(uint& currArg, const std::vector<char*>& argv, Targ&) const {
              return true;
            }
  protected:
    Tval     _default;
    setfun_t _setfun; // member function of Targ to set argument/option
};

template<class Targ>
class argdesc_t<Targ, bool> : public argdescbase_t<Targ> {
  public:
    typedef Targ arg_t;
    typedef bool val_t;
    typedef void (arg_t::*setfun_t)(const val_t&);
    typedef argdescbase_t<Targ> super_t;
  public:
    argdesc_t(const char*  aName,
              const val_t& aDefault,
              setfun_t     aSetFun, // member function of Targ to set argument/option
              const char*  aDesc) : argdescbase_t<Targ>(aName, aDesc), _default(aDefault), _setfun(aSetFun) {}
    virtual ~argdesc_t() {}
  public:
    inline val_t default_value() const { return _default; }
  public:
    virtual std::ostream& print_usage(std::ostream& os) const {
              os << "  "
                 << std::setw(15) << super_t::name() << "  "
                 << super_t::desc() << std::endl
                 << std::string(19, ' ') << "default: " << default_value() 
                 << std::endl;
              return os;
            }
    virtual int parse_args(int& currArg, const int noArgs, const char* argv[], Targ& aArg) const {
                   if(0 == strcmp(argdescbase_t<Targ>::_name, argv[currArg])) {
                     ++currArg;
                     (aArg.*_setfun)(!default_value());
                     return 1;
                   }
                   return 0;
                 }
    virtual int parse_args2(uint& currArg, const std::vector<char*>& argv, Targ& aArg) const {
                   if(0 == strcmp(argdescbase_t<Targ>::_name, argv[currArg])) {
                     ++currArg;
                     (aArg.*_setfun)(!default_value());
                     return 1;
                   }
                   return 0;
                 }
  protected:
    val_t     _default;
    setfun_t  _setfun; // member function of Targ to set argument/option

};

template<class Targ>
class argdesc_t<Targ, char> : public argdescbase_t<Targ> {
  public:
    typedef Targ arg_t;
    typedef char val_t;
    typedef void (arg_t::*setfun_t)(const val_t&);
    typedef argdescbase_t<Targ> super_t;
  public:
    argdesc_t(const char*  aName,
              const val_t& aDefault,
              setfun_t     aSetFun, // member function of Targ to set argument/option
              const char*  aDesc) : argdescbase_t<Targ>(aName, aDesc), _default(aDefault), _setfun(aSetFun) {}
    virtual ~argdesc_t() {}
  public:
    inline val_t default_value() const { return _default; }
  public:
    virtual std::ostream& print_usage(std::ostream& os) const {
              os << "  "
                 << std::setw(15) << super_t::name() << "  "
                 << super_t::desc() << std::endl
                 << std::string(19, ' ') << "default: " << default_value()
                 << std::endl;
              return os;
            }
    virtual int parse_args(int& currArg, const int noArgs, const char* argv[], Targ& aArg) const {
                   if(0 == strcmp(argdescbase_t<Targ>::_name, argv[currArg])) {
                     ++currArg;
                     (aArg.*_setfun)(*argv[currArg]);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }
    virtual int parse_args2(uint& currArg, const std::vector<char*>& argv, Targ& aArg) const {
                   if(0 == strcmp(argdescbase_t<Targ>::_name, argv[currArg])) {
                     ++currArg;
                     (aArg.*_setfun)(*argv[currArg]);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }
  protected:
    val_t     _default;
    setfun_t  _setfun; // member function of Targ to set argument/option
};


template<class Targ>
class argdesc_t<Targ, int> : public argdescbase_t<Targ> {
  public:
    typedef Targ arg_t;
    typedef int  val_t;
    typedef void (arg_t::*setfun_t)(const val_t&);
    typedef argdescbase_t<Targ> super_t;
  public:
    argdesc_t(const char*  aName,
              const val_t& aDefault,
              setfun_t     aSetFun, // member function of Targ to set argument/option
              const char*  aDesc) : argdescbase_t<Targ>(aName, aDesc), _default(aDefault), _setfun(aSetFun) {}
    virtual ~argdesc_t() {}
  public:
    inline val_t default_value() const { return _default; }
  public:
    virtual std::ostream& print_usage(std::ostream& os) const {
              os << "  "
                 << std::setw(15) << super_t::name() << "  "
                 << super_t::desc() << std::endl
                 << std::string(19, ' ') << "default: " << default_value() 
                 << std::endl;
              return os;
            }
    virtual int parse_args(int& currArg, const int noArgs, const char* argv[], Targ& aArg) const {
                   if(0 == strcmp(super_t::_name, argv[currArg])) {
                     ++currArg;
                     if(noArgs <= currArg) {
                       std::cerr << super_t::desc() << " needs integer value" << std::endl;
                       return -1;
                     }
                     char* lEnd = 0;
                     int lVal = strtol(argv[currArg], &lEnd, 10);
                     if(argv[currArg] == lEnd) {
                       std::cerr << super_t::desc() << " needs integer value" << std::endl;
                       return -1;
                     }
                     (aArg.*_setfun)(lVal);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }
    virtual int parse_args2(uint& currArg, const std::vector<char*>& argv, Targ& aArg) const {
                   if(0 == strcmp(super_t::_name, argv[currArg])) {
                     ++currArg;
                     if(argv.size() <= currArg) {
                       std::cerr << super_t::desc() << " needs integer value" << std::endl;
                       return -1;
                     }
                     char* lEnd = 0;
                     int lVal = strtol(argv[currArg], &lEnd, 10);
                     if(argv[currArg] == lEnd) {
                       std::cerr << super_t::desc() << " needs integer value" << std::endl;
                       return -1;
                     }
                     (aArg.*_setfun)(lVal);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }
    
  protected:
    val_t    _default;
    setfun_t _setfun;  // member function of Targ to set argument/option
};

template<class Targ>
class argdesc_t<Targ, unsigned int> : public argdescbase_t<Targ> {
  public:
    typedef Targ arg_t;
    typedef unsigned int  val_t;
    typedef void (arg_t::*setfun_t)(const val_t&);
    typedef argdescbase_t<Targ> super_t;
  public:
    argdesc_t(const char*  aName,
              const val_t& aDefault,
              setfun_t     aSetFun, // member function of Targ to set argument/option
              const char*  aDesc,
              const int    aBase = 10) 
            : argdescbase_t<Targ>(aName, aDesc), _default(aDefault), _setfun(aSetFun), _base(aBase) {}
    virtual ~argdesc_t() {}
  public:
    inline val_t default_value() const { return _default; }
    inline int   base() const { return _base; }
  public:
    virtual std::ostream& print_usage(std::ostream& os) const {
              os << "  "
                 << std::setw(15) << super_t::name() << "  "
                 << super_t::desc() << std::endl
                 << std::string(19, ' ') << "default: " << default_value() 
                 << std::endl;
              return os;
            }
    virtual int parse_args(int& currArg, const int noArgs, const char* argv[], Targ& aArg) const {
                   if(0 == strcmp(super_t::_name, argv[currArg])) {
                     ++currArg;
                     if(noArgs <= currArg) {
                       std::cerr << super_t::desc() << " needs non-negative integer value" << std::endl;
                       return -1;
                     }
                     char* lEnd = 0;
                     long int lVal = strtol(argv[currArg], &lEnd, base());
                     if((argv[currArg] == lEnd) || (0 > lVal)) {
                       std::cerr << super_t::desc() << " needs non-negative integer value" << std::endl;
                       return -1;
                     }
                     (aArg.*_setfun)(lVal);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }
    virtual int parse_args2(uint& currArg, const std::vector<char*>& argv, Targ& aArg) const {
                   if(0 == strcmp(super_t::_name, argv[currArg])) {
                     ++currArg;
                     if(argv.size() <= currArg) {
                       std::cerr << super_t::desc() << " needs non-negative integer value" << std::endl;
                       return -1;
                     }
                     char* lEnd = 0;
                     int lVal = strtol(argv[currArg], &lEnd, 10);
                     if((argv[currArg] == lEnd) || (0 > lVal)) {
                       std::cerr << super_t::desc() << " needs non-negative integer value" << std::endl;
                       return -1;
                     }
                     (aArg.*_setfun)(lVal);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }
  protected:
    val_t    _default;
    setfun_t _setfun;  // member function of Targ to set argument/option
    int      _base;
};


template<class Targ>
class argdesc_t<Targ, float> : public argdescbase_t<Targ> {
  public:
    typedef Targ arg_t;
    typedef float  val_t;
    typedef void (arg_t::*setfun_t)(const val_t&);
    typedef argdescbase_t<Targ> super_t;
  public:
    argdesc_t(const char*  aName,
              const val_t& aDefault,
              setfun_t     aSetFun, // member function of Targ to set argument/option
              const char*  aDesc) : argdescbase_t<Targ>(aName, aDesc), _default(aDefault), _setfun(aSetFun) {}
    virtual ~argdesc_t() {}
  public:
    inline val_t default_value() const { return _default; }
  public:
    virtual std::ostream& print_usage(std::ostream& os) const {
              os << "  "
                 << std::setw(15) << super_t::name() << "  "
                 << super_t::desc() << std::endl
                 << std::string(19, ' ') << "default: " << default_value() 
                 << std::endl;
              return os;
            }
    virtual int parse_args(int& currArg, const int noArgs, const char* argv[], Targ& aArg) const {
                   if(0 == strcmp(super_t::_name, argv[currArg])) {
                     ++currArg;
                     if(noArgs <= currArg) {
                       std::cerr << super_t::desc() << " needs float value" << std::endl;
                       return -1;
                     }
                     char* lEnd = 0;
                     val_t lVal = strtod(argv[currArg], &lEnd);
                     if(argv[currArg] == lEnd) {
                       std::cerr << super_t::desc() << " needs float value" << std::endl;
                       return -1;
                     }
                     (aArg.*_setfun)(lVal);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }
    virtual int parse_args2(uint& currArg, const std::vector<char*>& argv, Targ& aArg) const {
                   if(0 == strcmp(super_t::_name, argv[currArg])) {
                     ++currArg;
                     if(argv.size() <= currArg) {
                       std::cerr << super_t::desc() << " needs float value" << std::endl;
                       return -1;
                     }
                     char* lEnd = 0;
                     val_t lVal = (float) strtod(argv[currArg], &lEnd);
                     if(argv[currArg] == lEnd) {
                       std::cerr << super_t::desc() << " needs float value" << std::endl;
                       return -1;
                     }
                     (aArg.*_setfun)(lVal);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }
  protected:
    val_t    _default;
    setfun_t _setfun;  // member function of Targ to set argument/option
};

template<class Targ>
class argdesc_t<Targ, double> : public argdescbase_t<Targ> {
  public:
    typedef Targ arg_t;
    typedef double  val_t;
    typedef void (arg_t::*setfun_t)(const val_t&);
    typedef argdescbase_t<Targ> super_t;
  public:
    argdesc_t(const char*  aName,
              const val_t& aDefault,
              setfun_t     aSetFun, // member function of Targ to set argument/option
              const char*  aDesc) : argdescbase_t<Targ>(aName, aDesc), _default(aDefault), _setfun(aSetFun) {}
    virtual ~argdesc_t() {}
  public:
    inline val_t default_value() const { return _default; }
  public:
    virtual std::ostream& print_usage(std::ostream& os) const {
              os << "  "
                 << std::setw(15) << super_t::name() << "  "
                 << super_t::desc() << std::endl
                 << std::string(19, ' ') << "default: " << default_value() 
                 << std::endl;
              return os;
            }
    virtual int parse_args(int& currArg, const int noArgs, const char* argv[], Targ& aArg) const {
                   if(0 == strcmp(super_t::_name, argv[currArg])) {
                     ++currArg;
                     if(noArgs <= currArg) {
                       std::cerr << super_t::desc() << " needs float value" << std::endl;
                       return -1;
                     }
                     char* lEnd = 0;
                     val_t lVal = strtod(argv[currArg], &lEnd);
                     if(argv[currArg] == lEnd) {
                       std::cerr << super_t::desc() << " needs float value" << std::endl;
                       return -1;
                     }
                     (aArg.*_setfun)(lVal);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }
    virtual int parse_args2(uint& currArg, const std::vector<char*>& argv, Targ& aArg) const {
                   if(0 == strcmp(super_t::_name, argv[currArg])) {
                     ++currArg;
                     if(argv.size() <= currArg) {
                       std::cerr << super_t::desc() << " needs float value" << std::endl;
                       return -1;
                     }
                     char* lEnd = 0;
                     val_t lVal = strtod(argv[currArg], &lEnd);
                     if(argv[currArg] == lEnd) {
                       std::cerr << super_t::desc() << " needs float value" << std::endl;
                       return -1; }
                     (aArg.*_setfun)(lVal);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }

  protected:
    val_t    _default;
    setfun_t _setfun;  // member function of Targ to set argument/option
};


template<class Targ>
class argdesc_t<Targ, std::string> : public argdescbase_t<Targ> {
  public:
    typedef Targ arg_t;
    typedef std::string  val_t;
    typedef void (arg_t::*setfun_t)(const val_t&);
    typedef argdescbase_t<Targ> super_t;
  public:
    argdesc_t(const char*  aName,
              const val_t& aDefault,
              setfun_t     aSetFun, // member function of Targ to set argument/option
              const char*  aDesc) : argdescbase_t<Targ>(aName, aDesc), _default(aDefault), _setfun(aSetFun) {}
    virtual ~argdesc_t() {}
  public:
    inline const val_t& default_value() const { return _default; }
  public:
    virtual std::ostream& print_usage(std::ostream& os) const {
              os << "  "
                 << std::setw(15) << super_t::name() << "  "
                 << super_t::desc() << std::endl
                 << std::string(19, ' ') << "default: " << default_value() 
                 << std::endl;
              return os;
            }
    virtual int parse_args(int& currArg, const int noArgs, const char* argv[], Targ& aArg) const {
                   if(0 == strcmp(super_t::_name, argv[currArg])) {
                     ++currArg;
                     if(noArgs <= currArg) {
                       std::cerr << super_t::desc() << " needs string value" << std::endl;
                       return -1;
                     }
                     (aArg.*_setfun)(argv[currArg]);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }
    virtual int parse_args2(uint& currArg, const std::vector<char*>& argv, Targ& aArg) const {
                   if(0 == strcmp(super_t::_name, argv[currArg])) {
                     ++currArg;
                     if(argv.size() <= currArg) {
                       std::cerr << super_t::desc() << " needs string value" << std::endl;
                       return -1;
                     }
                     (aArg.*_setfun)(argv[currArg]);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }
  protected:
    val_t    _default;
    setfun_t _setfun;  // member function of Targ to set argument/option
};

template<class Targ>
class argdesc_t<Targ, constcharptr_t> : public argdescbase_t<Targ> {
  public:
    typedef Targ arg_t;
    typedef constcharptr_t val_t;
    typedef void (arg_t::*setfun_t)(const val_t&);
    typedef argdescbase_t<Targ> super_t;
  public:
    argdesc_t(const char*  aName,
                    val_t& aDefault,
              setfun_t     aSetFun, // member function of Targ to set argument/option
              const char*  aDesc) : argdescbase_t<Targ>(aName, aDesc), _default(aDefault), _setfun(aSetFun) {}
    virtual ~argdesc_t() {}
  public:
    inline const val_t& default_value() const { return _default; }
  public:
    virtual std::ostream& print_usage(std::ostream& os) const {
              os << "  "
                 << std::setw(15) << super_t::name() << "  "
                 << super_t::desc() << std::endl
                 << std::string(19, ' ') << "default: " << default_value()
                 << std::endl;
              return os;
            }
    virtual int parse_args(int& currArg, const int noArgs, const char* argv[], Targ& aArg) const {
                   if(0 == strcmp(super_t::_name, argv[currArg])) {
                     ++currArg;
                     if(noArgs <= currArg) {
                       std::cerr << super_t::desc() << " needs string value" << std::endl;
                       return -1;
                     }
                     (aArg.*_setfun)(argv[currArg]);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }
    virtual int parse_args2(uint& currArg, const std::vector<char*>& argv, Targ& aArg) const {
                   if(0 == strcmp(super_t::_name, argv[currArg])) {
                     ++currArg;
                     if(argv.size() <= currArg) {
                       std::cerr << super_t::desc() << " needs string value" << std::endl;
                       return -1;
                     }
                     (aArg.*_setfun)(argv[currArg]);
                     ++currArg;
                     return 1;
                   }
                   return 0;
                 }
  protected:
    val_t    _default;
    setfun_t _setfun;  // member function of Targ to set argument/option
};


template<class Targ>
bool
parse_args(const int   aStart, // starting index
           const int   aNoArgs, 
           const char* argv[], 
           std::vector<argdescbase_t<Targ>*>& aArgDesc,
           Targ& aArgStruct) {
  typedef unsigned int uint;
  int lCurrArg = aStart;


  while(lCurrArg < aNoArgs) {
    // std::cout << "no argdesc: " << aArgDesc.size() << std::endl; 
    // std::cout << "arg: " << lCurrArg << " " << argv[lCurrArg] << std::endl;
    bool lArgCouldBeSuccessfullyHandled = false;
    for(uint i = 0; i < aArgDesc.size(); ++i) {
      // std::cout << "  test: " << aArgDesc[i]->name() << std::endl;

      int lReturnCode = aArgDesc[i]->parse_args(lCurrArg, aNoArgs, argv, aArgStruct);
      if(0 > lReturnCode) { return false; }
      else if(0 == lReturnCode) { continue; }
      else if(0 < lReturnCode) { lArgCouldBeSuccessfullyHandled = true; break; }
    }
    if(!lArgCouldBeSuccessfullyHandled) {
      return false;
    }
  }
  return true;
}


// same just parse an additional string

template<class Targ>
bool
parse_args(const int   aStart, // starting index
           const int   aNoArgs,
           const char* argv[],
           std::vector<argdescbase_t<Targ>*>& aArgDesc,
           Targ&       aArgStruct,
           void (Targ::*aSetStringFun)(const std::string&)) {
  typedef unsigned int uint;
  int lCurrArg = aStart;

  while(lCurrArg < aNoArgs) {
    // std::cout << "no argdesc: " << aArgDesc.size() << std::endl;
    // std::cout << "arg: " << lCurrArg << " " << argv[lCurrArg] << std::endl;
    bool lArgCouldBeSuccessfullyHandled = false;
    for(uint i = 0; i < aArgDesc.size(); ++i) {
      // std::cout << "  test: " << aArgDesc[i]->name() << std::endl;

      int lReturnCode = aArgDesc[i]->parse_args(lCurrArg, aNoArgs, argv, aArgStruct);
      if(0 > lReturnCode) { return false; }
      else if(0 == lReturnCode) { continue; }
      else if(0 < lReturnCode) { lArgCouldBeSuccessfullyHandled = true; break; }
    }
    if(!lArgCouldBeSuccessfullyHandled) {
      if( ((lCurrArg + 1) == aNoArgs) &&
          0 != aSetStringFun) {
        (aArgStruct.*aSetStringFun)(std::string(argv[lCurrArg]));
        ++lCurrArg;
        break;
      }
      return false;
    }
  }
  return true;
}

// the following collects unparsable items into a string container
template<class Targ>
bool
parse_args_collect(const int   aStart, // starting index
                   const int   aNoArgs,
                   const char* argv[],
                   std::vector<argdescbase_t<Targ>*>& aArgDesc,
                   Targ& aArgStruct,
                   std::vector<std::string>& aUnparsed) {
  typedef unsigned int uint;
  int lCurrArg = aStart;


  while(lCurrArg < aNoArgs) {
    // std::cout << "no argdesc: " << aArgDesc.size() << std::endl; 
    // std::cout << "arg: " << lCurrArg << " " << argv[lCurrArg] << std::endl;
    bool lArgCouldBeSuccessfullyHandled = false;
    for(uint i = 0; i < aArgDesc.size(); ++i) {
      // std::cout << "  test: " << aArgDesc[i]->name() << std::endl;

      int lReturnCode = aArgDesc[i]->parse_args(lCurrArg, aNoArgs, argv, aArgStruct);
      if(0 > lReturnCode) { return false; }
      else if(0 == lReturnCode) { continue; }
      else if(0 < lReturnCode) { lArgCouldBeSuccessfullyHandled = true; break; }
    }
    if(!lArgCouldBeSuccessfullyHandled) {
      aUnparsed.push_back(std::string(argv[lCurrArg++]));
    }
  }
  return true;
}



template<class Targ>
bool
parse_args2(const uint   aStart, // starting index
            const std::vector<char*> argv,
            std::vector<argdescbase_t<Targ>*>& aArgDesc,
            Targ& aArgStruct) {
  typedef unsigned int uint;
  uint lCurrArg = aStart;


  while(lCurrArg < argv.size()) {
    // std::cout << "no argdesc: " << aArgDesc.size() << std::endl; 
    // std::cout << "arg: " << lCurrArg << " " << argv[lCurrArg] << std::endl;
    bool lArgCouldBeSuccessfullyHandled = false;
    for(uint i = 0; i < aArgDesc.size(); ++i) {
      // std::cout << "  test: " << aArgDesc[i]->name() << std::endl;

      int lReturnCode = aArgDesc[i]->parse_args2(lCurrArg, argv, aArgStruct);
      if(0 > lReturnCode) { return false; }
      else if(0 == lReturnCode) { continue; }
      else if(0 < lReturnCode) { lArgCouldBeSuccessfullyHandled = true; break; }
    }
    if(!lArgCouldBeSuccessfullyHandled) {
      return false;
    }
  }
  return true;
}



template<class Targ>
void
print_usage(std::ostream& os, const char* aProgrammName, std::vector<argdescbase_t<Targ>*>& aArgDesc) {
  std::cout << "Usage of " << aProgrammName << std::endl;
  for(uint i = 0; i < aArgDesc.size(); ++i) {
    aArgDesc[i]->print_usage(os);
  }
}

template<class Targ>
void
print_usage(std::ostream& os, const char* aProgrammName, std::vector<argdescbase_t<Targ>*>& aArgDesc, const char* x) {
  std::cout << "Usage: " << aProgrammName << ' ' << x << std::endl;
  for(uint i = 0; i < aArgDesc.size(); ++i) {
    aArgDesc[i]->print_usage(os);
  }
}




#endif
