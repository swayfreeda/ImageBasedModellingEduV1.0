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
#ifndef H_COLDET_IMPL
#define H_COLDET_IMPL

#include "coldet.h"
#include "box.h"
#include "cdmath3d.h"
#include <vector>

__CD__BEGIN

class CollisionModel3DImpl : public CollisionModel3D
{
  float m_Radius;
public:
  CollisionModel3DImpl(bool Static);
  void setTriangleNumber(int num) { if (!m_Final) m_Triangles.reserve(num); }

  void addTriangle(float x1, float y1, float z1,
                   float x2, float y2, float z2,
                   float x3, float y3, float z3)
  {
    addTriangle(Vector3D(x1,y1,z1),
                Vector3D(x2,y2,z2),
                Vector3D(x3,y3,z3));
  }
  void addTriangle(const float v1[3], const float v2[3], const float v3[3])
  {
    addTriangle(Vector3D(v1[0],v1[1],v1[2]),
                Vector3D(v2[0],v2[1],v2[2]),
                Vector3D(v3[0],v3[1],v3[2]));
  }
  void addTriangle(const Vector3D& v1, const Vector3D& v2, const Vector3D& v3);
  void finalize();

  float getRadius() { return m_Radius; }

  void setTransform(const float m[16]) { setTransform(*(Matrix3D*)m); }
  void setTransform(const Matrix3D& m);

  bool collision(CollisionModel3D* other, 
                 int AccuracyDepth, 
                 int MaxProcessingTime,
                 float* other_transform);

  bool rayCollision(const float origin[3], const float direction[3], bool closest,
                    float segmin, float segmax);
  bool sphereCollision(const float origin[3], float radius);

  bool getCollidingTriangles(float t1[9], float t2[9], bool ModelSpace);
  bool getCollidingTriangles(int& t1, int& t2);
  bool getCollisionPoint(float p[3], bool ModelSpace);


  int getTriangleIndex(BoxedTriangle* bt)
  {
    return int(bt-&(*m_Triangles.begin()));
  }

  /** Stores all the actual triangles.  Other objects will use
      pointers into this array.
  */
  std::vector<BoxedTriangle> m_Triangles;
  /** Root of the hierarchy tree */
  BoxTreeInnerNode           m_Root;
  /** The current transform and its inverse */
  Matrix3D                   m_Transform,m_InvTransform;
  /** The triangles that last collided */
  Triangle                   m_ColTri1,m_ColTri2;
  /** The indices of the triangles that last collided */
  int                        m_iColTri1,m_iColTri2;
  /** The collision point of the last test */
  Vector3D                   m_ColPoint;

  /** Type of the last collision test */
  enum { Models, Ray, Sphere }       
                             m_ColType;
  /** Flag for indicating the model is finalized. */
  bool                       m_Final;
  /** Static models will maintain the same transform for a while
      so the inverse transform is calculated each set instead
      of in the collision test. */
  bool                       m_Static;
};

__CD__END

#endif // H_COLDET_IMPL
