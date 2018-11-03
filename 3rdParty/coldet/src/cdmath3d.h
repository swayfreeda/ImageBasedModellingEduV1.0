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
#ifndef H_cdmath3d
#define H_cdmath3d

#include <cmath>
#include <iostream>
#include <fstream>

struct Vector3D;
struct Matrix3;
struct Matrix3D;
struct Plane;

/*
template<class T>
T sqr(const T& a) { return a*a; }
*/

inline float flabs(float f) { return (f>=0.0f?f:-f); }
const float epsilon=1e-8f;
inline bool IsZero(float f) { return flabs(f)<epsilon; }

Vector3D operator*(float scalar, const Vector3D& v);
/** Dot product. */
float    operator*(const Vector3D& v1, const Vector3D& v2);
Vector3D operator+(const Vector3D& v1, const Vector3D& v2);
Vector3D operator-(const Vector3D& v1, const Vector3D& v2);
Vector3D CrossProduct(const Vector3D& v1, const Vector3D& v2);
Matrix3D operator*(const Matrix3D& m1, const Matrix3D& m2);
Matrix3D operator*(float scalar, const Matrix3D& m);

struct Vector3D
{
  float x,y,z;
  static const Vector3D Zero;

  Vector3D() {}
  Vector3D(float X, float Y, float Z) : x(X), y(Y), z(Z) {}
  Vector3D(const Vector3D& v) : x(v.x), y(v.y), z(v.z) {}

  Vector3D& operator+=(const Vector3D& v) { x+=v.x; y+=v.y; z+=v.z; return *this; }
  Vector3D& operator*=(float s) { x*=s; y*=s; z*=s; return *this; }
  Vector3D& operator/=(float s) { return *this *= (1.0f/s); }
  bool      operator==(const Vector3D& v) { return x==v.x && y==v.y && z==v.z; }

  const float* get() const { return &x; }

  Vector3D operator-       () const { return Vector3D(-x,-y,-z); }
  float    SquareMagnitude () const { return x*x+y*y+z*z; }
  float    Magnitude       () const { return (float)sqrt(SquareMagnitude()); }
  Vector3D Normalized      () const { return (1.0f/Magnitude())*(*this); }
  Vector3D Absolute        () const { return Vector3D(fabs(x),fabs(y),fabs(z)); }
  float    operator[] (int i) const { return ((float*)&x)[i]; }
  float&   operator[] (int i)       { return ((float*)&x)[i]; }
};

inline std::ostream& operator<< (std::ostream& os, const Vector3D& v)
{
  return os << v.x << ',' << v.y << ',' << v.z;
}

#define _11 sclr.s11
#define _12 sclr.s12
#define _13 sclr.s13
#define _14 sclr.s14
#define _21 sclr.s21
#define _22 sclr.s22
#define _23 sclr.s23
#define _24 sclr.s24
#define _31 sclr.s31
#define _32 sclr.s32
#define _33 sclr.s33
#define _34 sclr.s34
#define _41 sclr.s41
#define _42 sclr.s42
#define _43 sclr.s43
#define _44 sclr.s44

/** 3x3 matrix */
struct Matrix3
{
  union {
    struct { float s11,s12,s13,
                   s21,s22,s23,
                   s31,s32,s33; } sclr;
    float m[3][3];
  };
  static const Matrix3 Identity;

  Vector3D& baseRow(int i) { return *((Vector3D*)m[i]); }
  float  operator() (int i, int j) const { return m[i][j]; }
  float& operator() (int i, int j)       { return m[i][j]; }
};

/** 4x4 matrix, used for transformations. */
struct Matrix3D
{
  union {
    struct { float s11,s12,s13,s14,
                   s21,s22,s23,s24,
                   s31,s32,s33,s34,
                   s41,s42,s43,s44; } sclr;
    float m[4][4];
  };
  static const Matrix3D Identity;

  Matrix3D() {}

  Matrix3D(float f11, float f12, float f13, float f14,
           float f21, float f22, float f23, float f24,
           float f31, float f32, float f33, float f34,
           float f41, float f42, float f43, float f44)
  {
    _11=f11; _12=f12; _13=f13; _14=f14;
    _21=f21; _22=f22; _23=f23; _24=f24;
    _31=f31; _32=f32; _33=f33; _34=f34;
    _41=f41; _42=f42; _43=f43; _44=f44;
  }

  Matrix3D& operator*= (const Matrix3D& m)
  {
    return *this = *this * m;
  }

  friend Matrix3D PitchMatrix3D(const float theta);
  friend Matrix3D YawMatrix3D(const float theta);
  friend Matrix3D RollMatrix3D(const float theta);
  void rotate(const Vector3D& v);

  Matrix3D Inverse() const;
  Matrix3D Adjoint() const;
  float Determinant() const;

  float  operator() (int i, int j) const { return m[i][j]; }
  float& operator() (int i, int j)       { return m[i][j]; }
};

/** 3D Plane.  Used in conjunction with triangles. */
struct Plane
{
  Vector3D normal;
  float    d;

  Plane(const Vector3D& a, const Vector3D& b, const Vector3D& c)
  {
    normal = CrossProduct(b - a, c - a).Normalized();
    d = -normal * a;
  }

  float Classify(const Vector3D& v)
  {
    return v * normal + d;
  }
};

inline Vector3D operator* (float scalar, const Vector3D& v)
{
  return Vector3D(scalar*v.x,scalar*v.y,scalar*v.z);
}

inline Vector3D operator+ (const Vector3D& v1, const Vector3D& v2)
{
  return Vector3D(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z);
}

inline Vector3D operator- (const Vector3D& v1, const Vector3D& v2)
{
  return Vector3D(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z);
}

inline float operator* (const Vector3D& v1, const Vector3D& v2)
{
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

inline Vector3D CrossProduct(const Vector3D& v1, const Vector3D& v2)
{
  return Vector3D(v1.y*v2.z-v2.y*v1.z,
                  v1.z*v2.x-v2.z*v1.x,
                  v1.x*v2.y-v2.x*v1.y);
}

inline Vector3D Transform(const Vector3D& v, const Matrix3D& m)
{
  return Vector3D(v.x*m._11 + v.y*m._21 + v.z*m._31 + m._41,
                  v.x*m._12 + v.y*m._22 + v.z*m._32 + m._42,
                  v.x*m._13 + v.y*m._23 + v.z*m._33 + m._43);
}

inline Vector3D rotateVector(const Vector3D& v, const Matrix3D& m)
{
  return Vector3D(v.x*m._11 + v.y*m._21 + v.z*m._31,
                  v.x*m._12 + v.y*m._22 + v.z*m._32,
                  v.x*m._13 + v.y*m._23 + v.z*m._33);
}

inline Matrix3D operator*(float scalar, const Matrix3D& m)
{
  return Matrix3D(scalar*m(0,0),scalar*m(0,1),scalar*m(0,2),scalar*m(0,3),
                  scalar*m(1,0),scalar*m(1,1),scalar*m(1,2),scalar*m(1,3),
                  scalar*m(2,0),scalar*m(2,1),scalar*m(2,2),scalar*m(2,3),
                  scalar*m(3,0),scalar*m(3,1),scalar*m(3,2),scalar*m(3,3));
}

inline Matrix3D operator*(const Matrix3D& m1, const Matrix3D& m2)
{
  return Matrix3D(
    m1._11*m2._11 + m1._12*m2._21 + m1._13*m2._31 + m1._14*m2._41,
    m1._11*m2._12 + m1._12*m2._22 + m1._13*m2._32 + m1._14*m2._42,
    m1._11*m2._13 + m1._12*m2._23 + m1._13*m2._33 + m1._14*m2._43,
    m1._11*m2._14 + m1._12*m2._24 + m1._13*m2._34 + m1._14*m2._44,
    m1._21*m2._11 + m1._22*m2._21 + m1._23*m2._31 + m1._24*m2._41,
    m1._21*m2._12 + m1._22*m2._22 + m1._23*m2._32 + m1._24*m2._42,
    m1._21*m2._13 + m1._22*m2._23 + m1._23*m2._33 + m1._24*m2._43,
    m1._21*m2._14 + m1._22*m2._24 + m1._23*m2._34 + m1._24*m2._44,
    m1._31*m2._11 + m1._32*m2._21 + m1._33*m2._31 + m1._34*m2._41,
    m1._31*m2._12 + m1._32*m2._22 + m1._33*m2._32 + m1._34*m2._42,
    m1._31*m2._13 + m1._32*m2._23 + m1._33*m2._33 + m1._34*m2._43,
    m1._31*m2._14 + m1._32*m2._24 + m1._33*m2._34 + m1._34*m2._44,
    m1._41*m2._11 + m1._42*m2._21 + m1._43*m2._31 + m1._44*m2._41,
    m1._41*m2._12 + m1._42*m2._22 + m1._43*m2._32 + m1._44*m2._42,
    m1._41*m2._13 + m1._42*m2._23 + m1._43*m2._33 + m1._44*m2._43,
    m1._41*m2._14 + m1._42*m2._24 + m1._43*m2._34 + m1._44*m2._44);
}

inline Matrix3D
TranslateMatrix3D(const Vector3D& v)
{
  return Matrix3D(1.0f,0.0f,0.0f,0.0f,
                  0.0f,1.0f,0.0f,0.0f,
                  0.0f,0.0f,1.0f,0.0f,
                   v.x, v.y, v.z,1.0f);
}


inline Matrix3D
ScaleMatrix3D(const Vector3D& v)
{
   return Matrix3D( v.x,0.0f,0.0f,0.0f,
                   0.0f, v.y,0.0f,0.0f,
                   0.0f,0.0f, v.z,0.0f,
                   0.0f,0.0f,0.0f,1.0f);
}


inline Matrix3D
ScaleMatrix3D(const float s)
{
   return ScaleMatrix3D(Vector3D(s,s,s));
}


inline Matrix3D
PitchMatrix3D(const float c, const float s)
{
   return Matrix3D(1.0f, 0.0f, 0.0f, 0.0f,
                   0.0f,    c,   -s, 0.0f,
                   0.0f,    s,    c, 0.0f,
                   0.0f, 0.0f, 0.0f, 1.0f);
}


inline Matrix3D
PitchMatrix3D(const float theta)
{
   return PitchMatrix3D((float) cos(theta), (float) sin(theta));
}


inline Matrix3D
YawMatrix3D(const float c, const float s)
{
   return Matrix3D(   c, 0.0f,    s, 0.0f,
                   0.0f, 1.0f, 0.0f, 0.0f,
                     -s, 0.0f,    c, 0.0f,
                   0.0f, 0.0f, 0.0f, 1.0f);
}


inline Matrix3D
YawMatrix3D(const float theta)
{
   return YawMatrix3D((float) cos(theta), (float) sin(theta));
}


inline Matrix3D
RollMatrix3D(const float c, const float s)
{
   return Matrix3D(c,   -s,    0.0f, 0.0f,
                   s,    c,    0.0f, 0.0f,
                   0.0f, 0.0f, 1.0f, 0.0f,
                   0.0f, 0.0f, 0.0f, 1.0f);
}


inline Matrix3D
RollMatrix3D(const float theta)
{
   return RollMatrix3D((float) cos(theta), (float) sin(theta));
}

inline void
Matrix3D::rotate(const Vector3D& v)
{
   if (v.x!=0.0f) *this = PitchMatrix3D(v.x) * (*this);
   if (v.y!=0.0f) *this = YawMatrix3D  (v.y) * (*this);
   if (v.z!=0.0f) *this = RollMatrix3D (v.z) * (*this);
}



template<class T>
inline T Max(T a, T b)
{
  return (a>b ? a : b);
}

template<class T>
inline T Min(T a, T b)
{
  return (a<b ? a : b);
}

#undef _11
#undef _12
#undef _13
#undef _14
#undef _21
#undef _22
#undef _23
#undef _24
#undef _31
#undef _32
#undef _33
#undef _34
#undef _41
#undef _42
#undef _43
#undef _44

#endif // H_cdmath3d
