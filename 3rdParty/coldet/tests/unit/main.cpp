#include <iostream>
#include <coldet.h>
#include <cdmath3d.h>
#include <multiobject.h>
#include <cstdlib>
#include <map>
#include <vector>

using namespace COLDET;

inline float rand_float()
{
  const float norm=1.0f/RAND_MAX;
  return ((rand()*norm)-0.5f);
}

/// Generate a random triangle soup
CollisionModel3D* create_test_model()
{
  static int seed=1;
  srand(seed++);
  CollisionModel3D* model=newCollisionModel3D(false);
  float v[9];
  for(int i=0;i<32;++i)
  {
    for(int j=0;j<9;++j)
      v[j]=rand_float();
    model->addTriangle(&v[0],&v[3],&v[6]);
  }
  model->finalize();
  return model;
}

/// Represents a game model.  A collision model and a transform
class PositionedModel : public TransformUpdater
{
  CollisionModel3D* model;
  Matrix3D          matrix;
  int               ID;
public:
  PositionedModel() : model(0){}
  ~PositionedModel() { delete model; }

  void create_model()
  {
    model=create_test_model();
    matrix=Matrix3D::Identity;
  }

  CollisionModel3D* get_cm()
  {
    return model;
  }

  const float* update() const
  {
    return reinterpret_cast<const float*>(&matrix);
  }

  void set_id(int id)
  {
    ID=id;
  }

  Vector3D get_pos() const
  {
    return Vector3D(matrix(3,0),matrix(3,1),matrix(3,2));
  }

  void move(const float* pos)
  {
    matrix(3,0)=pos[0];
    matrix(3,1)=pos[1];
    matrix(3,2)=pos[2];
  }
};

typedef std::vector<PositionedModel> pm_vec;

int main(int argc, char* argv[])
{
  // Test the prune & sweep against a brute force n^2 sphere system
  MultiObjectSystem* mos1=newSpheresSystem();
  MultiObjectSystem* mos2=newSweepPruneSystem(1024);
  const float POS_RANGE=100;
  const int MODELS = 1000;
  pm_vec models(MODELS);
  const int ITERATIONS = 10000;
  for(int i=0;i<MODELS;++i)
  {
    models[i].create_model();
    float pos[3];
    for(int j=0;j<3;++j)
      pos[j]=POS_RANGE*rand_float();
    models[i].move(pos);
    mos1->add_object(models[i].get_cm(),&models[i]);
    mos2->add_object(models[i].get_cm(),&models[i]);
  }
  for(int i=0;i<ITERATIONS;++i)
  {
    std::cout << i << "\r";
    std::cout.flush();
    int n1=mos1->check_collisions(true);
    int n2=mos2->check_collisions(true);
    if (n1!=n2) std::cout << "Mismatch in number of collisions: " << n1 << " : " << n2 << "\n";
    else
    {
      std::multimap<int,int> cols1,cols2;
      for(int i=0;i<n1;++i)
      {
        CollisionDetails details;
        mos1->get_collision(i,details);
        cols1.insert(std::make_pair(details.id1,details.id2));
        mos2->get_collision(i,details);
        cols2.insert(std::make_pair(details.id1,details.id2));
      }
      if (cols1!=cols2)
      {
        std::cout << "Collisions mismatch\n";
      }
      /// Move a random model 
      int id=(rand()%MODELS);
      float pos[3];
      for(int j=0;j<3;++j)
        pos[j]=POS_RANGE*rand_float();
      models[id].move(pos);
    }
  }

  delete mos1;
  delete mos2;
  std::cout << "Done\n";
  return 0;
}