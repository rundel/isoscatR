// Code based on:

//  boost progress.hpp header file 
//  Copyright Beman Dawes 1994-99.  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//  See http://www.boost.org/libs/timer for documentation.


#ifndef RSCAT_PROGRESS_HPP
#define RSCAT_PROGRESS_HPP

#include <boost/timer.hpp>
#include <boost/utility.hpp>  // for noncopyable
#include <iostream>           // for ostream, cout, etc
#include <string>             // for string

class progress_display : private boost::noncopyable
{
public:
    explicit progress_display( unsigned long expected_count,
                               std::ostream & os = std::cout,
                               std::string prefix = "",
                               unsigned int prefix_width = 8,
                               unsigned int width = 50)
      : m_os(os), _prefix(prefix), _prefix_width(prefix_width), _width(width)
    { 
        restart(expected_count); 
    }
private:
    std::ostream &m_os;
    unsigned long _expected_count, _count;
    unsigned int _width, _prefix_width;
    std::string _prefix;
    boost::timer timer;
    
public:

    void restart( unsigned long expected_count )
    {
        _count = 0;
        _expected_count = expected_count;
        if ( !_expected_count ) 
            _expected_count = 1;  
        timer.restart();
        
        display_tic();
        
    }

    unsigned long operator+=( unsigned long increment )
    //  Effects: Display appropriate progress tic if needed.
    //  Postconditions: count()== original count() + increment
    //  Returns: count().
    {
        _count += increment;
        display_tic(); 
        
        return _count;
    }

    unsigned long operator++()
    {
        return operator+=( 1 ); 
    }
    
    unsigned long count() const 
    {
        return _count;
    }
    
    unsigned long expected_count() const 
    {
        return _expected_count;
    }


    void display_tic()
    {   
        double perc = static_cast<double>(_count)/_expected_count;
        unsigned int tics = static_cast<unsigned int>( 50.0*perc );
        
        m_os << "\r";
        if (_prefix!="")
            m_os << std::setw(_prefix_width) << _prefix;
        
        m_os << " |";
        for(unsigned int i=0; i<_width; ++i) {
            char c = (i <= tics) ? '=' : ' ';
            m_os << c;
        }
        
        double time = timer.elapsed();
        double remaining = time / perc;
        
        m_os << "| " << std::setw(3) << static_cast<unsigned int>(perc*100.0) << "%";
        m_os << " [" << (int)(time) << "s of ";
        
        if (perc > 0.03) m_os << (int)(remaining) << "s";
        else             m_os << "--s";
        
        m_os << "]" << std::flush;
        
        if (_count == _expected_count)
            m_os << std::endl;
    }
};

#endif