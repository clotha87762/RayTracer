#ifndef __final_h__
#define __final_h__

#include <vector>

#include "algebra3.h"
#include "photon.h"
#include "Mesh.h"
#include "Camera.h"

using namespace std;


const float PLANE_DISTANCE = 1.0f;
const float ASPECT_RATIO = 1.0f; // 想要不同的fov就可以設成非0的值  x:y

int kdMaxDepth = 10;
int kdThreshold = 10;
int tracingDepth = 10;

vec3 backGroundColor(0,0,0);
float ambientIntensity = 1.0f;
float shadowBias = 0.0001f;
float refractBias = 0.01f;
float boundingBoxMinimunSize = 0.01f;

float airReractivity = 1.0f;
int SSAA_RATIO = 1;


// ======================================================



vec3 upVec(0,1.0f,0);

//Camera camera;




bool BARYCENTRIC_INTERPOLATION_NORMAL =false;
bool rayHitAnything = false;











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


bool TestVisualizeLight();

bool IsInShadow(vec3 point , vec3 lightPos , SimpleMesh* currentMesh , vec3 normal);
vec3 SimpleMeshRayTrace( vec3 ori , vec3 dir , int depth ,Camera cam , float currentNr);
vec3 RayCastSynthesisColorImage(Camera  cam , vec2 xy);
void GenerateColorImage(Camera cam ,char* inputName);
void PhotonMappingImageSynthesis(Camera cam );

#endif
