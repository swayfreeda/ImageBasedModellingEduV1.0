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
#ifndef H_COLDET_MULTI_OBJECT
#define H_COLDET_MULTI_OBJECT

#include "coldet.h"
#include "sysdep.h"

class TransformUpdater
{
public:
  virtual ~TransformUpdater() {}
  virtual const float* update() const = 0;
};

__CD__BEGIN

struct CollisionDetails
{
  int   id1,id2;
  float point[3];
  float t1[9],t2[9];
};

/** Generic class to handle multiple objects scenes */
class MultiObjectSystem
{
public:
  virtual ~MultiObjectSystem() {}

  /** Adds a model to the system at a specific position.  model is not owned by the system
      Returns a unique ID of the model, for later move / remove */
  virtual int add_object(CollisionModel3D* model, const float* position) = 0;

  /** Adds a model to the system at a specific position.  model is not owned by the system
      Returns a unique ID of the model, for later move / remove */
  virtual int add_object(CollisionModel3D* model, const TransformUpdater* updater) = 0;

  /** Adds an object that is not composed of polygons, but can be
      represented by a sphere.
      Returns a unique ID of the model, for later move / remove */
  virtual int add_sphere_object(float radius, const float* position) = 0;

  /** Remove model from the system.
      'id' is the value returned when the model was added */
  virtual void remove_object(int id) = 0;

  /** Move the model in space.  Will incrementally update the system.
      'id' is the value returned when the model was added
      'position' is the new position (3D vector)   */
  virtual void move_object(int id, const float* new_position) = 0;

  /** Retrieve collision model from id */
  virtual CollisionModel3D* get_collision_model(int id) = 0;

  /** Check for collisions between all pairs.
      'exact' indicates activation of triangle accurate intersection test
      without which, collisions are only estimated by proximity
      Returns the number of collisions found */
  virtual int check_collisions(bool exact) = 0;

  /** Retrieve collision details.
      'index' is between 0..(N-1)  N was the value returned from check_collisions
      id1,id2 are the two colliding models.  If exact collision was used,
      the collision details can be retrieved by supplying the details struct */
  virtual bool get_collision(int index, CollisionDetails& details) = 0;

};

MultiObjectSystem* newSpheresSystem();
MultiObjectSystem* newSweepPruneSystem(int max_objects);

__CD__END

#endif // H_COLDET_MULTI_OBJECT

