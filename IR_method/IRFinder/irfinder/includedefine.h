#ifndef INCLUDEDEFINE_DEF
#define INCLUDEDEFINE_DEF

#include <cstring>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <sys/types.h>
#include <vector>
#include <map>
#include <algorithm> // std::sort
#include <functional> // std::function

//__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");


using namespace std;


#define DEF_lineLengthMax 10000
#define DEF_adaptLengthMax 500


// infix_iterator.h 
// 
// Lifted from Jerry Coffin's 's prefix_ostream_iterator 
template <class T, 
          class charT=char, 
          class traits=std::char_traits<charT> > 
class infix_ostream_iterator : 
    public std::iterator<std::output_iterator_tag,void,void,void,void> 
{ 
    std::basic_ostream<charT,traits> *os; 
    charT const* delimiter; 
    bool first_elem; 
public: 
    typedef charT char_type; 
    typedef traits traits_type; 
    typedef std::basic_ostream<charT,traits> ostream_type; 
    infix_ostream_iterator(ostream_type& s) 
        : os(&s),delimiter(0), first_elem(true) 
    {} 
    infix_ostream_iterator(ostream_type& s, charT const *d) 
        : os(&s),delimiter(d), first_elem(true) 
    {} 
    infix_ostream_iterator<T,charT,traits>& operator=(T const &item) 
    { 
        // Here's the only real change from ostream_iterator: 
        // Normally, the '*os << item;' would come before the 'if'. 
        if (!first_elem && delimiter != 0) 
            *os << delimiter; 
        *os << item; 
        first_elem = false; 
        return *this; 
    } 
    infix_ostream_iterator<T,charT,traits> &operator*() { 
        return *this; 
    } 
    infix_ostream_iterator<T,charT,traits> &operator++() { 
        return *this; 
    } 
    infix_ostream_iterator<T,charT,traits> &operator++(int) { 
        return *this; 
    } 
};     



#endif