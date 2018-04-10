#include "algebra3.h"








class SimpleMesh{

    public :

    SimpleMesh(){

    }

    virtual void PrintElements(){

    };

    virtual bool IntersectionTest(vec3 rayDir, vec3 rayOri){
        return false;
    };

};


class Sphere : public SimpleMesh{
    public :
    vec3 center;
    float radius;

    Sphere(){

    }

    Sphere(vec3 center,float radius){
        this->center = center;
        this->radius = radius;
    }

    void PrintElements();

    bool IntersectionTest(vec3 rayDir , vec3 rayOri);
};

class Triangle : public SimpleMesh{
    public:
    vec3 vertices[3];
    Triangle(vec3 v1,vec3 v2,vec3 v3){
        vertices[0] = v1;
        vertices[1] = v2;
        vertices[2] = v3;
    }
    Triangle(){

    }

    void PrintElements();
    bool IntersectionTest(vec3 rayDir, vec3 rayOri);
};

class Camera{

public :

    vec3 center;
    vec3 up;
    vec3 viewDirection;
    float fovX , fovY ;
    float planeDistance;
    int resX;
    int resY;

    vec3 planeCenter;
    vec3 planeNormal;
    vec3 planeXAxis;
    vec3 planeYAxis;
    vec3 planeCorners[4]; // topLeft ,topRight, bottomLeft , bottomRight


    Camera(){
        planeDistance = 1.0f; // planeDistance好像沒在input格式中定義...QAQ
    }

    Camera(vec3 center , vec3 up, vec3 viewDirection,float fovX,float fovY,int resX,int resY){
        this->center = center;
        this->up = up;
        this->viewDirection = viewDirection;
        this->fovX = fovX;
        this->fovY = fovY;
        this->resX = resX;
        this->resY = resY;
    }

    vec3 GetPixelCenter(vec2 offset);
    void PrintElements();
    void CalcPlaneParameters();
};

bool RayCastInterSectionTest(Camera  cam , vec2 xy );
void GenerateBinaryIntersectionTestImage(Camera cam);

