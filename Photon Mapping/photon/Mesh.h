#ifndef __Mesh_h__
#define __Mesh_h__


#include <vector>
#include "algebra3.h"

using namespace std;

class SimpleMesh{

    public :

    bool isCurrentTraceMesh = false;
    bool isCurrentRayIn = false;

    vec3 staticAmbient; // Calculate every time the transform of the mesh changed
    vec3 staticDiffuse; // Calculate every time the transform of the mesh changed
    vec3 tempIntersectionPoint;

    vec3 surfColor;
    vec3 emissionColor;
    float Ks , Kd , Ka;
    float specularity;
    float reflect;
    float refract;
    float Nr;

    SimpleMesh() : staticAmbient(0,0,0),staticDiffuse(0,0,0),
    surfColor(1.0f,1.0f,1.0f),Ks(1.0f),Kd(1.0f),Ka(1.0f),
    specularity(6.0f),reflect(0.9f),refract(0),Nr(1.0f),
    emissionColor(0,0,0)
    {


    }

    virtual void PrintElements(){

    };

    virtual vec3 GetNormal(vec3 point){

    };

    virtual vector<vec3> GetPositions(){

    };

    virtual bool IntersectionTest(vec3 rayDir, vec3 rayOri){
        return false;
    };

    /*
    virtual void PreComputeAmbientDiffuse(){
    }
    */
    virtual vec3 GetIntersectionPoint(vec3 rayDir, vec3 rayOri){

    }

};


class Sphere : public SimpleMesh{
    public :
    vec3 center;
    float radius;

    Sphere() : SimpleMesh(){

    }

    Sphere(vec3 center,float radius) : SimpleMesh(){
        this->center = center;
        this->radius = radius;
    }

    vector<vec3> GetPositions();
    void PrintElements();
    //void PreComputeAmbientDiffuse();

    vec3 GetNormal(vec3 point);

    bool IntersectionTest(vec3 rayDir , vec3 rayOri);
    vec3 GetIntersectionPoint(vec3 rayDir, vec3 rayOri);
};

class Triangle : public SimpleMesh{
    public:
    vec3 vertices[3];
    vec3 normal[3];

    Triangle(vec3 v1,vec3 v2,vec3 v3):SimpleMesh(){
        vertices[0] = v1;
        vertices[1] = v2;
        vertices[2] = v3;
    }

    Triangle() : SimpleMesh(){

    }


    vec3 GetNormal(vec3 point);

    vector<vec3> GetPositions();
    void PrintElements();
    //void PreComputeAmbientDiffuse();
    bool IntersectionTest(vec3 rayDir, vec3 rayOri);
    vec3 GetIntersectionPoint(vec3 rayDir, vec3 rayOri);
};

#endif // __Mesh_h__
