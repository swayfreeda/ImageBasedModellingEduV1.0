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
#ifndef H_BOX
#define H_BOX

#include <vector>
#include "cdmath3d.h"
#include "sysdep.h"

__CD__BEGIN

/** Stores rotation vectors used in the intersection tests, 
    to avoid recalculating them each time. */
class RotationState
{
public:
  RotationState(const Matrix3D& transform);
  Vector3D N[3];
  Matrix3D t;
};

/** AABB class, with support for testing against OBBs. */
class Box
{
public:
  /** Default constructor */
  Box() {}
  /** Construct from scalar corner position and size */
  Box(float x, float y, float z, float sx, float sy, float sz) 
    : m_Pos(x,y,z), m_Size(sx,sy,sz), 
      m_Center(x+0.5f*sx,y+0.5f*sy,z+0.5f*sz) {}
  /** Construct from corner position and size */
  Box(const Vector3D& pos, const Vector3D& size) 
    : m_Pos(pos), m_Size(size), m_Center(pos+0.5f*size) {}
  /** Copy constructor */
  Box(const Box& b) : m_Pos(b.m_Pos), m_Size(b.m_Size), m_Center(b.m_Center) {}
  virtual ~Box() {}
  /** Returns the box's position */
  const Vector3D& getPosition() const { return m_Pos; }
  /** Returns the sizes of the box's edges */
  const Vector3D& getSize() const { return m_Size; }
  /** Returns the center position of the box */
  const Vector3D& getCenter() const { return m_Center; }
  /** Returns the volume of the box */
  float getVolume() const { return m_Size.x*m_Size.y*m_Size.z; }
  /** Ray intersection */
  bool intersect(const Vector3D& O, const Vector3D& D);
  /** Line segment intersection */
  bool intersect(const Vector3D& O, const Vector3D& D, float segmax);
  /** Sphere intersection */
  bool intersect(const Vector3D& O, float radius);
  /** Point in box */
  bool intersect(const Vector3D& p) const;
  /** Aligned box intersection */
  bool intersect(const Box& b);
  /** Oriented box intersection. */
  bool intersect(const Box& b, RotationState& rs);

  /** Position of box corner */
  Vector3D m_Pos;
  /** Size of box box edges */
  Vector3D m_Size;
  /** Position of box center.  m_Pos+0.5f*m_Size;  */
  Vector3D m_Center;
};

/** A single triangle in the model */
class Triangle
{
public:
  /** Default constructor */
  Triangle() {}
  /** Constructor to build a triangle from 3 points */
  Triangle(const Vector3D& _1, const Vector3D& _2, const Vector3D& _3);
  /** Tests for intersection with another triangle. */
  bool intersect(const Triangle& t) const;
  /** Tests for intersection with a ray (O origin, D direction)
      Returns true if collision occured.
      Outputs collision point in cp
      Outputs the distance from the origin to the collision point in tparm
      This distance is relative to the magnitude of D
      Allows testing against a finite segment, by specifying 
      the maximum length of the ray in segmax
      This length is also relative to the magnitude of D
  */
  bool intersect(const Vector3D& O, const Vector3D& D, Vector3D& cp, 
                 float& tparm, float segmax);
  /** Test for intersection with a sphere (O origin) 
      Returns true if collision occured.
      Outputs collision point in cp
  */
  bool intersect(const Vector3D& O, float radius, Vector3D& cp);

  Vector3D v1,v2,v3;
  Vector3D center;
};

class BoxedTriangle;

/** Base class for hierarchy tree nodes. */
class BoxTreeNode : public Box
{
public:
  /** Default constructor */
  BoxTreeNode() : Box() {}
  /** Constructor for a box from position and size */
  BoxTreeNode(const Vector3D& pos, const Vector3D& size) 
    : Box(pos,size) {}
  /** Returns true if the node is a leaf node. */
  virtual bool isLeaf() const  = 0;
  /** Returns the number of sons this node has */
  virtual int            getSonsNumber() = 0;
  /** Returns a son node, by index */ 
  virtual BoxTreeNode*   getSon(int which) = 0;
  /** Returns the number of triangles in this node.
      Only non-zero for leaf nodes. */
  virtual int            getTrianglesNumber() = 0;
  /** Returns the boxed triangle contained in this node
      by its index 
  */
  virtual BoxedTriangle* getTriangle(int which) = 0;
};

/** Inner node, containing other nodes. */
class BoxTreeInnerNode : public BoxTreeNode
{
public:
  BoxTreeInnerNode(const Vector3D& pos, const Vector3D& size, int logdepth) 
    : BoxTreeNode(pos,size), m_First(NULL), m_Second(NULL), 
      m_logdepth(logdepth), m_OwnFirst(true), m_OwnSecond(true) {}
  ~BoxTreeInnerNode()
  {
    if (m_OwnFirst)  delete m_First;
    if (m_OwnSecond) delete m_Second;
  }

  virtual bool isLeaf() const { return false; }
  /** Create the sons that will divide this box */
  int  createSons(const Vector3D& center);
  /** Recalculate the bounds of this box to fully contain
      all of its triangles
  */
  void recalcBounds(Vector3D& center);
  /** Recursively divide this box */
  int  divide(int p_depth);

  int getSonsNumber()
  {
    int n=0;
    if (m_First!=NULL) n++;
    if (m_Second!=NULL) n++;
    return n;
  }

  int getTrianglesNumber();
  BoxedTriangle* getTriangle(int which);

  BoxTreeNode* getSon(int which)
  {
    if (which==0) return m_First;
    if (which==1) return m_Second;
    return NULL;
  }

  BoxTreeNode*                m_First;
  bool                        m_OwnFirst;
  BoxTreeNode*                m_Second;
  bool                        m_OwnSecond;
  int                         m_logdepth;
  std::vector<BoxedTriangle*> m_Boxes;
};

/** Leaf node, containing 1 triangle. */
class BoxedTriangle : public BoxTreeNode, public Triangle
{
public:
  BoxedTriangle(const Vector3D& _1, const Vector3D& _2, const Vector3D& _3);
  virtual bool isLeaf() const { return true; }
  int getSonsNumber() { return 0; }
  BoxTreeNode* getSon(int which) { return NULL; }
  int getTrianglesNumber() { return 1; }
  BoxedTriangle* getTriangle(int which)
  {
    if (which==0) return this;
    return NULL;
  }

};

__CD__END

#endif // H_BOX
