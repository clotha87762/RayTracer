#include "algebra3.h"
#include <vector>

using namespace std;


class SimpleLight{
    public :
    float intensity;
    vec3 pos;
    vec3 color;

    SimpleLight() : intensity(1.0f) , pos(0,0,0) , color(1.0f,1.0f,1.0f){

    }



};

class PointLight : public SimpleLight{
    public :
    PointLight() : SimpleLight(){
    }
};

class SpotLight : public SimpleLight{
     public :
    vec3 dir;
    float angle;
    SpotLight() : SimpleLight(),dir(0,0,1.0f),angle(60.0f){
    }
};

class DirectionalLight : public SimpleLight{
     public :
    vec3 dir;
    DirectionalLight() : SimpleLight(),dir(0,0,1.0f){
    }
};

class SimpleMesh{

    public :

    bool isCurrentTraceMesh = false;
    bool isCurrentRayIn = false;

    vec3 staticAmbient; // Calculate every time the transform of the mesh changed
    vec3 staticDiffuse; // Calculate every time the transform of the mesh changed
    vec3 tempIntersectionPoint;

    vec3 surfColor;
    float Ks , Kd , Ka;
    float specularity;
    float reflect;
    float refract;
    float Nr;

    SimpleMesh() : staticAmbient(0,0,0),staticDiffuse(0,0,0),
    surfColor(1.0f,1.0f,1.0f),Ks(1.0f),Kd(1.0f),Ka(1.0f),
    specularity(6.0f),reflect(0.9f),refract(0),Nr(1.0f)
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

class Camera{

public :

    int ssaaRatio;

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
        ssaaRatio = 1;
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

class KDNode{

    public :
    KDNode* leftChild;
    KDNode* rightChild;
    vector<SimpleMesh*> kdMeshes;
    int splitAxis;
    float splitValue;
    vec3 center , extend;
    bool splitLeft,splitRight;
    bool nonZero;
    static int index;
    int myIndex;
    int layer;
    bool isLeaf;
    float modelScale;


    KDNode():splitAxis(0),splitLeft(false),splitRight(false),nonZero(false),isLeaf(false),modelScale(2.0f){
        leftChild = NULL;
        rightChild = NULL;
        myIndex = index++;
    };
    void Clean(){
        if(leftChild!=NULL){
            leftChild->Clean();
            delete leftChild;
        }
        if(rightChild!=NULL){
            rightChild->Clean();
            delete rightChild;
        }
    }
    void CalcBound();
    void BuildKD(int depth);
    bool HitTest(vec3 ori,vec3 ray);

};

class SimpleObject{
    public:
    //vector<SimpleMesh*> meshes;
    KDNode* kdRoot;
    //void BuildKDTree();
};


void MeasureRenderTime(Camera cam, char* inputName);

bool RayCastInterSectionTest(Camera  cam , vec2 xy );
void GenerateBinaryIntersectionTestImage(Camera cam);



bool IsInShadow(vec3 point , vec3 lightPos , SimpleMesh* currentMesh , vec3 normal);
vec3 SimpleMeshRayTrace(vec3 ori , vec3 dir , int depth ,Camera cam , float currentNr);
vec3 RayCastSynthesisColorImage(Camera  cam , vec2 xy);
void GenerateColorImage(Camera cam ,char* inputName);


