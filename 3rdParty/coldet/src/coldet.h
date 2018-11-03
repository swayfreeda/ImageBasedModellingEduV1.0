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
/** \file coldet.h
    3D Collision Detection

    Interface for the library.  
    Isolated from any implementation details.
*/
#ifndef H_COLDET
#define H_COLDET

#ifndef EXPORT
#define EXPORT
#endif

/** Collision Model.  Will represent the mesh to be tested for
    collisions.  It has to be notified of all triangles, via
    addTriangle()
    After all triangles are added, a call to finalize() will
    process the information and prepare for collision tests.
    Call collision() to check for a collision

    Note: Transformations must not contain scaling.
*/
class CollisionModel3D
{
public:
  virtual ~CollisionModel3D() {}

  /** Optional: Optimization for construction speed.
      If you know the number of triangles. */
  virtual void setTriangleNumber(int num) = 0;

  /** Use any of the forms of this functions to enter the coordinates
      of the model's triangles. */
  virtual void addTriangle(float x1, float y1, float z1,
                           float x2, float y2, float z2,
                           float x3, float y3, float z3) = 0;
  virtual void addTriangle(const float v1[3], const float v2[3], const float v3[3]) = 0;

  /** All triangles have been added, process model. */
  virtual void finalize() = 0;

  /** Returns the bounding sphere radius
      Note that this is not the optimal bounding sphere, but centered
      in the origin of the coordinate system of the triangles */
  virtual float getRadius() = 0;

  /** The the current affine matrix for the model.
      See transform.txt for format information */
  virtual void setTransform(const float m[16]) = 0;

  /** Check for collision with another model.
      Do not mix model types here.

      MaxProcessingTime determines the maximum time in milliseconds
      to check for collision.  If a rejection is not found by that
      time, the function will return true.

      AccuracyDepth is not yet supported.

      other_transform allows overriding the other model's 
      transform, by supplying an alternative one.
      This can be useful when testing a model against itself
      with different orientations.
  */
  virtual bool collision(CollisionModel3D* other,
                         int AccuracyDepth=-1,
                         int MaxProcessingTime=0,
                         float* other_transform=0) = 0;

  /** Returns true if the ray given in world space coordinates
      intersects with the object.  
      getCollidingTriangles() and getCollisionPoint() can be
      used to retrieve information about a collision.
      If closest if false, the first triangle that collides with
      the ray is used.  Otherwise the closest one will be used.
      Closest triangle searching will slow the test considerably.
      The default ray is a standard infinite ray.  However, using
      segmin and segmax you can define a line segment along the
      ray.
  */
  virtual bool rayCollision(const float origin[3],
                            const float direction[3],
                            bool closest=false,
                            float segmin=0.0f,
                            float segmax=3.4e+38F) = 0;

  /** Returns true if the given sphere collides with the model.
      getCollidingTriangles() and getCollisionPoint() can be
      used to retrieve information about a collision.
  */
  virtual bool sphereCollision(const float origin[3],
                               float radius) = 0;

  /** Retrieve the pair of triangles that collided.
      Only valid after a call to collision() that returned true.
      t1 is this model's triangle and t2 is the other one.
      In case of ray or sphere collision, only t1 will be valid.
      The coordinates will be in _this_ model's coordinate space,
      unless ModelSpace is false, in which case, coordinates will
      be transformed by the model's current transform to world space.
  */
  virtual bool getCollidingTriangles(float t1[9], float t2[9], bool ModelSpace=true) = 0;

  /** Retrieve the pair of triangles indices that collided.
      Only valid after a call to collision() that returned true.
      t1 belongs to _this_ model, while t2 is in the other one. 
  */
  virtual bool getCollidingTriangles(int& t1, int& t2) = 0;

  /** Retrieve the detected collision point. 
      Only valid after a call to collision()
      that returned true.
      The coordinates will be in _this_ model's coordinate space,
      unless ModelSpace is false, in which case, coordinates will
      be transformed by the model's current transform to world space.
  */
  virtual bool getCollisionPoint(float p[3], bool ModelSpace=true) = 0;
};

/** Timeout exception class.  Exception will be thrown if
    the detection algorithm could not complete within
    the given time limit. */
class TimeoutExpired {};

/** Inconsistency exception. Exception will be thrown if
    the model is inconsistent.  
    Examples: 
      Checking for collisions before calling finalize()
      Trying to add triangles after calling finalize()  */
class Inconsistency {};

/** Create a new collision model object.
    Use delete when finished with it. 

    Setting Static to true indicates that the model does not
    move a lot, and certain calculations can be done every time
    its transform changes instead of every collision test. 
*/
EXPORT CollisionModel3D* newCollisionModel3D(bool Static=false);



//////////////////////////////////////////////
// Utility Functions
//////////////////////////////////////////////

/** Checks for intersection between a ray and a sphere.
    center, radius define the sphere
    origin, direction define the ray
    point will contain point of intersection, if one is found.
*/
EXPORT bool SphereRayCollision(const float* center, float radius,
                               const float* origin, const float* direction,
                               float* point);

/** Checks for intersection between 2 spheres. */
EXPORT bool SphereSphereCollision(float c1[3], float r1,
                           float c2[3], float r2, float point[3]);



#endif // H_COLDET
