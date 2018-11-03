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
#include "box.h"

__CD__BEGIN

// point in box test
bool Box::intersect(const Vector3D& p) const
{
  const Vector3D& pos=getPosition();
  const Vector3D& s=getSize();
  if (p.x<pos.x || p.x>(pos.x+s.x)) return false;
  if (p.y<pos.y || p.y>(pos.y+s.y)) return false;
  if (p.z<pos.z || p.z>(pos.z+s.z)) return false;
  return true;
}

// Non-oriented intersection test
bool Box::intersect(const Box& b)
{
  const Vector3D& t1=getPosition();
  Vector3D t2=getPosition()+getSize();
  const Vector3D& p1=b.getPosition();
  Vector3D p2=b.getPosition()+b.getSize();
  return (Max(p1.x,t1.x) <= Min(p2.x,t2.x)  &&
          Max(p1.y,t1.y) <= Min(p2.y,t2.y)  &&
          Max(p1.z,t1.z) <= Min(p2.z,t2.z));
}

BoxedTriangle::BoxedTriangle(const Vector3D& _1, 
                             const Vector3D& _2, 
                             const Vector3D& _3)
                             : BoxTreeNode(), Triangle(_1,_2,_3)
{
  m_Pos.x=Min(Min(_1.x,_2.x),_3.x);
  m_Pos.y=Min(Min(_1.y,_2.y),_3.y);
  m_Pos.z=Min(Min(_1.z,_2.z),_3.z);
  Vector3D mx;
  mx.x=Max(Max(_1.x,_2.x),_3.x);
  mx.y=Max(Max(_1.y,_2.y),_3.y);
  mx.z=Max(Max(_1.z,_2.z),_3.z);
  m_Size=mx-getPosition();
  m_Center=getPosition()+0.5f*getSize();
}

int BoxTreeInnerNode::createSons(const Vector3D& center)
{
  int longest=0;
  Vector3D p=getPosition();
  Vector3D s=getSize();

  Vector3D dist(Vector3D::Zero);
  for(unsigned i=0;i<m_Boxes.size();i++)
  {
    BoxedTriangle* bt=m_Boxes[i];
    dist.x+=flabs(bt->center.x - center.x);
    dist.y+=flabs(bt->center.y - center.y);
    dist.z+=flabs(bt->center.z - center.z);
  }
  if (dist.y>dist.x && dist.y>dist.z) longest=1;
  else
  if (dist.z>dist.x && dist.z>dist.y) longest=2;

  float s1=center[longest]-p[longest];
  float s2=s[longest]-s1;
  s[longest]=s1;
  m_First=new BoxTreeInnerNode(p,s,m_logdepth);
  p[longest]+=s1;
  s[longest]=s2;
  m_Second=new BoxTreeInnerNode(p,s,m_logdepth);
  return longest;
}

void BoxTreeInnerNode::recalcBounds(Vector3D& center)
{
  if (m_Boxes.empty()) return;
  center=Vector3D::Zero;
  Vector3D mn(9e9f,9e9f,9e9f),mx(-9e9f,-9e9f,-9e9f);
  for(unsigned i=0;i<m_Boxes.size();i++)
  {
    BoxedTriangle* bt=m_Boxes[i];
    center+=bt->center;
    mn.x=Min(Min(Min(bt->v1.x,bt->v2.x),bt->v3.x),mn.x);
    mn.y=Min(Min(Min(bt->v1.y,bt->v2.y),bt->v3.y),mn.y);
    mn.z=Min(Min(Min(bt->v1.z,bt->v2.z),bt->v3.z),mn.z);
    mx.x=Max(Max(Max(bt->v1.x,bt->v2.x),bt->v3.x),mx.x);
    mx.y=Max(Max(Max(bt->v1.y,bt->v2.y),bt->v3.y),mx.y);
    mx.z=Max(Max(Max(bt->v1.z,bt->v2.z),bt->v3.z),mx.z);
  }
  center/=float(m_Boxes.size());
  m_Pos=mn;
  m_Size=mx-mn;
  if (m_Size.x==0.0f) { m_Size.x=0.002f; m_Pos.x-=0.001f; }
  if (m_Size.y==0.0f) { m_Size.y=0.002f; m_Pos.y-=0.001f; }
  if (m_Size.z==0.0f) { m_Size.z=0.002f; m_Pos.z-=0.001f; }
  m_Center=getPosition()+0.5f*getSize();
}

int BoxTreeInnerNode::divide(int p_depth)
{
  if (m_Boxes.empty()) return 0;
  Vector3D center;
  recalcBounds(center);
  int longest=createSons(center);
  BoxTreeInnerNode* f=static_cast<BoxTreeInnerNode*>(m_First);
  BoxTreeInnerNode* s=static_cast<BoxTreeInnerNode*>(m_Second);
  int depth=1;
  int bnum=m_Boxes.size();
#ifdef _DEBUG
  int fnum=0;
#endif
  for(int i=0;i<bnum;i++)
  {
    BoxedTriangle* bt=m_Boxes[i];
    if (bt->center[longest]<center[longest])
    {
      f->m_Boxes.push_back(bt);
      #ifdef _DEBUG
        fnum++;
      #endif
    }
    else
    {
      s->m_Boxes.push_back(bt);
    }
  }
  
  int b1num=f->m_Boxes.size();
  int b2num=s->m_Boxes.size();
  if ((b1num==bnum  ||  b2num==bnum))// && p_depth>m_logdepth)
  {
    delete m_First;  m_First=NULL;
    delete m_Second; m_Second=NULL;
    return depth+1;
  }
  
  m_Boxes.clear();
  if (f->m_Boxes.empty()) { delete m_First; m_First=NULL; }
  else
  if (f->m_Boxes.size()==1)
  {
    BoxedTriangle* bt=f->m_Boxes.back();
    delete m_First;
    m_OwnFirst=false;
    m_First=bt;
  } else depth=f->divide(p_depth+1);
  if (s->m_Boxes.empty()) { delete m_Second; m_Second=NULL; }
  else
  if (s->m_Boxes.size()==1)
  {
    BoxedTriangle* bt=s->m_Boxes.back();
    delete m_Second;
    m_OwnSecond=false;
    m_Second=bt;
  } else depth=Max(depth,s->divide(p_depth+1));
  return depth+1;
}

__CD__END
