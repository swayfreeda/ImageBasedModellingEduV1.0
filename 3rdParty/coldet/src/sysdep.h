/*   ColDet - C++ 3D Collision Detection Library
 *   Copyright (C) 2000-2013   Amir Geva
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 * 
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA  02111-1307, USA.
 *
 * Any comments, questions and bug reports send to:
 *   amirgeva@gmail.com
 *
 * Or visit the home page: http://sourceforge.net/projects/coldet/
 */
#ifndef H_SYSDEP
#define H_SYSDEP

#define __CD__BEGIN namespace COLDET {
#define __CD__END }

///////////////////////////////////////////////////
// g++ compiler on most systems
///////////////////////////////////////////////////
#ifdef GCC

unsigned get_tick_count();

///////////////////////////////////////////////////
// Windows compilers
///////////////////////////////////////////////////
#elif defined(WIN32)

  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>

#ifdef _WDLL

  #ifndef EXPORT
    #ifdef COLDET_EXPORTS
      #define EXPORT __declspec(dllexport)
    #else
      #define EXPORT __declspec(dllimport)
    #endif
  #endif

#else

  #ifndef EXPORT
  #define EXPORT
  #endif

#endif

  inline unsigned get_tick_count() { return GetTickCount(); }

///////////////////////////////////////////////////
// MacOS 9.0.4/MacOS X.  CodeWarrior Pro 6
// Thanks to Marco Tenuti for this addition
///////////////////////////////////////////////////
#elif defined(macintosh)
   typedef unsigned long DWORD;
   #include <Events.h>
   #define get_tick_count() ::TickCount()

#elif defined(CUSTOM)
//typedef unsigned DWORD;
double get_tick_count();
#else

#error No system specified (WIN32 GCC macintosh)

#endif

#ifndef EXPORT
  #define EXPORT
#endif

#endif // H_SYSDEP
