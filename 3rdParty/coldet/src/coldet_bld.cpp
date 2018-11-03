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
#include "sysdep.h"
#include "coldetimpl.h"

EXPORT CollisionModel3D* newCollisionModel3D(bool Static)
{
  return new COLDET::CollisionModel3DImpl(Static);
}

__CD__BEGIN

CollisionModel3DImpl::CollisionModel3DImpl(bool Static)
: m_Root(Vector3D::Zero, Vector3D::Zero,0),
  m_Transform(Matrix3D::Identity),
  m_InvTransform(Matrix3D::Identity),
  m_ColTri1(Vector3D::Zero,Vector3D::Zero,Vector3D::Zero),
  m_ColTri2(Vector3D::Zero,Vector3D::Zero,Vector3D::Zero),
  m_iColTri1(0),
  m_iColTri2(0),
  m_Final(false),
  m_Static(Static),
  m_Radius(0)
{}

void CollisionModel3DImpl::addTriangle(const Vector3D& v1, const Vector3D& v2, const Vector3D& v3)
{
  if (m_Final) throw Inconsistency();
  m_Triangles.push_back(BoxedTriangle(v1,v2,v3));
  m_Radius=Max(m_Radius,v1.SquareMagnitude());
  m_Radius=Max(m_Radius,v2.SquareMagnitude());
  m_Radius=Max(m_Radius,v3.SquareMagnitude());
}

void CollisionModel3DImpl::setTransform(const Matrix3D& m)
{
  m_Transform=m;
  if (m_Static) m_InvTransform=m_Transform.Inverse();
}

void CollisionModel3DImpl::finalize()
{
  if (m_Final) throw Inconsistency();
  // Prepare initial triangle list
  m_Final=true;
  m_Radius=sqrt(m_Radius);
  for(unsigned i=0;i<m_Triangles.size();i++)
  {
    BoxedTriangle& bt=m_Triangles[i];
    m_Root.m_Boxes.push_back(&bt);
  }
  int logdepth=0;
  for(int num=m_Triangles.size();num>0;num>>=1,logdepth++);
  m_Root.m_logdepth=int(logdepth*1.5f);
  m_Root.divide(0);
}

__CD__END
