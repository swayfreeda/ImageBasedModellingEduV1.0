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
#ifndef H_COLDET_SWEEP_PRUNE
#define H_COLDET_SWEEP_PRUNE

#include <iostream>
#include <vector>
#include <list>
#include "bitmatrix.h"
#include <set>
#include "cdmath3d.h"
#include "sysdep.h"

__CD__BEGIN

namespace SNP
{

struct Intersection
{
  int id1,id2;

  Intersection(int i1, int i2)
	: id1(i1), id2(i2) 
  {
    if (id2<id1) std::swap(id1,id2);
  }

  bool operator< (const Intersection& rhs) const
  {
    if (id1==rhs.id1) return id2<rhs.id2;
    return id1<rhs.id1;
  }
};

inline bool operator== (const Intersection& a, const Intersection& b)
{
  return a.id1==b.id1 && a.id2==b.id2;
}
inline bool operator!= (const Intersection& a, const Intersection& b)
{
  return !(a==b);
}

typedef std::set<Intersection> intr_set;


typedef unsigned short ushort;

struct EndPoint
{
  float      pos;
  ushort     id;
  bool       start;
};

inline bool operator== (const EndPoint& p1, const EndPoint& p2)
{
  return p1.pos==p2.pos && p1.id==p2.id && p1.start==p2.start;
}
inline bool operator!= (const EndPoint& p1, const EndPoint& p2)
{
  return !(p1==p2);
}


struct Index
{
  Index() : start(0), stop(0) {}
  void update(int index, bool strt)
  {
    if (strt) start=index;
    else      stop=index;
  }
  ushort start,stop;
};

inline bool operator== (const Index& i1, const Index& i2)
{
  return i1.start==i2.start && i1.stop==i2.stop;
}
inline bool operator!= (const Index& i1, const Index& i2)
{
  return !(i1==i2);
}

typedef std::vector<EndPoint> ep_vec;
typedef std::vector<Index> indices;

class SweepPruneAxe
{
  typedef  BitMatrixRow kbits;
  typedef  BitMatrix matrix;
  ep_vec   m_EndPoints;
  indices  m_Indices;
  matrix   m_Matrix;
  int      m_N;

  void set_cell(int i, int j, bool state=true)
  {
    m_Matrix.set(i,j,state);
    m_Matrix.set(j,i,state);
  }

  void reset_cell(int i, int j)
  {
    m_Matrix.reset(j,i);
    m_Matrix.reset(i,j);
  }

  void move_end_point_forward(int which, float d)
  {
    int n=m_EndPoints.size();
    EndPoint& ep=m_EndPoints[which];
    ep.pos+=d;
    int to_index=n-1;
    for(int i=which+1;i<n;++i)
    {
      if (m_EndPoints[i].pos>ep.pos)
      {
        if (i==(which+1)) return; // point does not move. no update needed
        to_index=i-1;
        break;
      }
    }
    EndPoint copy=ep;
    for(int j=which;j<to_index;++j)
    {
      EndPoint& cur=m_EndPoints[j];
      cur=m_EndPoints[j+1];
      m_Indices[cur.id].update(j,cur.start);
      if (copy.start)
      {
        if (!cur.start) 
          reset_cell(copy.id,cur.id);
      }
      else
      {
        set_cell(copy.id,cur.id);
      }
    }
    m_EndPoints[to_index]=copy;
    m_Indices[copy.id].update(to_index,copy.start);
  }

  void move_end_point_backward(int which, float d)
  {
    EndPoint& ep=m_EndPoints[which];
    ep.pos+=d;
    int to_index=0;
    for(int i=which-1;i>=0;--i)
    {
      if (m_EndPoints[i].pos<=ep.pos)
      {
        if (i==(which-1)) return; // point does not move. no update needed
        to_index=i+1;
        break;
      }
    }
    EndPoint copy=ep;
    for(int j=which;j>to_index;--j)
    {
      EndPoint& cur=m_EndPoints[j];
      cur=m_EndPoints[j-1];
      m_Indices[cur.id].update(j,cur.start);
      if (!copy.start)
      {
        if (cur.start) 
          reset_cell(copy.id,cur.id);
      }
      else
      {
        set_cell(copy.id,cur.id);
      }
    }
    m_EndPoints[to_index]=copy;
    m_Indices[copy.id].update(to_index,copy.start);
  }

  friend class SweepPrune;
public:
  SweepPruneAxe()
  {}
    
  void initialize(unsigned N)
  {
    m_Matrix.resize(N,N);
    m_Indices.resize(N);
    m_N=N;
    m_EndPoints.reserve(N*2);
  }

  const ep_vec&  get_end_points() const       { return m_EndPoints;  }
  const indices& get_indices()    const       { return m_Indices;    }
  const matrix&  get_matrix()     const       { return m_Matrix;     }
  const kbits&   get_matrix_row(int id) const { return m_Matrix[id]; }

  bool move_interval(int id, float d)
  {
    if (id<0 || id>=m_N) return false;
    Index& ind=m_Indices[id];
    if (ind.start==ind.stop) return false;
    if (d>0)
    {
      move_end_point_forward(ind.stop,d);
      move_end_point_forward(ind.start,d);
    }
    else
    {
      move_end_point_backward(ind.start,d);
      move_end_point_backward(ind.stop,d);
    }
    return true;
  }

  void remove_interval(int id)
  {
    if (id<0 || id>=m_N) return;
    Index& ind=m_Indices[id];
    if (ind.start==ind.stop) return; // No such interval
    for(int i=0;i<m_N;++i)
      m_Matrix[i].reset(unsigned(id));
    unsigned maxi=m_EndPoints.size()-2;
    unsigned step=1;
    for(unsigned i=ind.start;i<maxi;++i)
    {
      if ((i+1)==ind.stop) step=2;
      EndPoint& cur=m_EndPoints[i];
      cur=m_EndPoints[i+step];
      Index& cur_index=m_Indices[cur.id];
      if (cur.start) cur_index.start=i; else cur_index.stop=i;
    }
    m_EndPoints.pop_back();
    m_EndPoints.pop_back();
    ind.start=ind.stop=0;
  }

  bool add_interval(int id, float from, float to)
  {
    if (id<0 || id>=m_N) return false;
    Index& ind=m_Indices[id];
    if (ind.start!=ind.stop) return false;
    kbits& row=m_Matrix[id];
    row.reset();
    EndPoint np[] = { {float(from),ushort(id),true}, {float(to),ushort(id),false} };
    int n=m_EndPoints.size();
    EndPoint* cur_ins=np;
    int phase=0;
    for(int i=0;i<n;++i)
    {
      EndPoint& cur=m_EndPoints[i];
      kbits& sec_row=m_Matrix[cur.id];
      if (phase<2 && cur.pos>cur_ins->pos)
      {
        m_EndPoints.insert(m_EndPoints.begin()+i,*cur_ins);
        if (phase==0) ind.start=i;
        else          ind.stop=i;
        ++n;
        ++cur_ins;
        ++phase;
      }
      else
      {
        if (phase>0)
        {
          if (cur.start)
            m_Indices[cur.id].start+=phase;
          else
            m_Indices[cur.id].stop+=phase;
        }
        if (phase==1 || (phase==0 && cur.start))
        { 
          row.set(cur.id); sec_row.set(id); 
        }
        else
        if (phase==0)
        {
          row.reset(cur.id); sec_row.reset(id);
        }
        else
        if (!cur.start) // closing endpoint in phase 2
        {
          if (m_Indices[cur.id].start<ind.stop)
          {
            row.set(cur.id); sec_row.set(id);
          }
        }
      }
    }
    for(;phase<2;++phase)
    {
      if (phase==0) ind.start=m_EndPoints.size();
      else          ind.stop=m_EndPoints.size();
      m_EndPoints.push_back(*cur_ins);
      ++cur_ins;
    }
    return true;
  }
};

class SweepPrune
{
  typedef SweepPruneAxe::kbits kbits;
  SweepPruneAxe    m_Axes[3];
  intr_set         m_Intersections;
  int              m_N;

  void toggle_intersections(int id, const kbits& row)
  {
    if (row.any())
    {
      for(int i=0;i<m_N;++i)
      {
        if (row.test(i))
        {
          Intersection intr(id,i);
          std::pair<intr_set::iterator,bool> res=m_Intersections.insert(intr);
          if (!res.second) m_Intersections.erase(res.first);
        }
      }
    }
  }

  void get_collisions(int id, kbits& row)
  {
    row=m_Axes[0].get_matrix_row(id);
    row&=m_Axes[1].get_matrix_row(id);
    row&=m_Axes[2].get_matrix_row(id);
  }

public:
  SweepPrune(unsigned N)
    : m_N(N)
  {
    for(int i=0;i<3;++i)
      m_Axes[i].initialize(N);
  }

  void clear_intersections()
  {
    m_Intersections.clear();
  }

  const intr_set& get_intersections() const
  {
    return m_Intersections;
  }

  void add_box(int id, const float* pos, float radius)
  {
    kbits pre_row(m_N),post_row(m_N);
    get_collisions(id,pre_row);
    for(int i=0;i<3;++i)
      m_Axes[i].add_interval(id,pos[i]-radius,pos[i]+radius);
    get_collisions(id,post_row);
    pre_row^=post_row;
    toggle_intersections(id,pre_row);
  }

  void remove_box(int id)
  {
    for(int i=0;i<3;++i)
      m_Axes[i].remove_interval(id);
  }

  void move_box(int id, const float* delta)
  {
    kbits pre_row(m_N),post_row(m_N);
    get_collisions(id,pre_row);
    for(int i=0;i<3;++i)
      m_Axes[i].move_interval(id,delta[i]);
    get_collisions(id,post_row);
    pre_row^=post_row;
    toggle_intersections(id,pre_row);
  }
};

} // namespace SNP

__CD__END

#endif // H_COLDET_SWEEP_PRUNE
